//! Production-native optimized weighted graph construction.
//!
//! This keeps the legacy pair-scoring arithmetic and discovery orientation while
//! replacing repeated record decoding, interval scans, and the global checked-pair
//! hash set with immutable per-record views and an offline interval sweep.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::sync::Arc;
use std::time::Instant;

use log::info;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Record;

use crate::haplotype_determination::{
    extract_error_vector, extract_hap_vector, extract_query_seq, get_read_id,
};
use crate::structs::{AlleleDepthMap, HaplotypeConfig, PhasingGraphResult, ReadPair, ReadPairMap};

const GAP_MARKER: i16 = -10;
const OFFSET_UNIT: i16 = 10;
const ROLE_INDEX: u8 = 1;
const ROLE_QUERY: u8 = 2;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairStatus {
    Compatible,
    Blocked,
    Indeterminate,
}

struct SegmentView {
    start: i64,
    end: i64,
    stream: u32,
    symbols: Vec<u8>,
    nongap_index: Vec<u32>,
    nongap_pos: Vec<i64>,
    axis_token: Vec<u32>,
    state: Arc<Vec<i16>>,
    uncertainty: Arc<Vec<f32>>,
}

impl SegmentView {
    fn build(
        record: &Record,
        state: Arc<Vec<i16>>,
        uncertainty: Arc<Vec<f32>>,
    ) -> Result<Self, String> {
        let (symbols, positions) = extract_query_seq(record).map_err(|error| error.to_string())?;
        let mut nongap_index = Vec::with_capacity(positions.len());
        let mut nongap_pos = Vec::with_capacity(positions.len());
        for (index, &position) in positions.iter().enumerate() {
            if position != -1 {
                nongap_index.push(index as u32);
                nongap_pos.push(position);
            }
        }

        let start = record.reference_start();
        let end = record.reference_end();
        let span = usize::try_from(end - start)
            .map_err(|_| format!("invalid record interval {start}..{end}"))?;
        let cigar = record.cigar();
        let axis_token = (0..span)
            .map(|offset| {
                cigar
                    .read_pos((start + offset as i64) as u32, true, true)
                    .ok()
                    .flatten()
                    .unwrap_or(u32::MAX)
            })
            .collect();

        let stream = u32::try_from(record.tid())
            .map_err(|_| "record has no reference stream".to_string())?;
        Ok(Self {
            start,
            end,
            stream,
            symbols,
            nongap_index,
            nongap_pos,
            axis_token,
            state,
            uncertainty,
        })
    }

    fn slice_symbols(&self, start: i64, end: i64) -> &[u8] {
        let lo = self
            .nongap_pos
            .partition_point(|&position| position < start);
        let hi = self.nongap_pos.partition_point(|&position| position < end);
        if lo >= hi {
            return &[];
        }
        let first = self.nongap_index[lo] as usize;
        let last = self.nongap_index[hi - 1] as usize;
        &self.symbols[first..=last]
    }

    fn slice_state(&self, start: i64, end: i64) -> &[i16] {
        let start_offset = (start - self.start) as usize;
        let end_offset = (end - self.start) as usize;
        if start_offset >= self.state.len() {
            return &[];
        }
        &self.state[start_offset..end_offset.min(self.state.len())]
    }

    fn symbol_at(&self, target: i64) -> Option<u8> {
        let slot = self
            .nongap_pos
            .partition_point(|&position| position < target);
        (self.nongap_pos.get(slot) == Some(&target))
            .then(|| self.symbols[self.nongap_index[slot] as usize])
    }
}

struct EntityView {
    segments: [SegmentView; 2],
}

struct BuiltEntity {
    view: EntityView,
    first_id: String,
    second_id: String,
}

#[inline]
fn is_snv_value(value: i16) -> bool {
    value == -4 || (value > 1 && value % OFFSET_UNIT == 6)
}

#[inline]
fn is_indel_value(value: i16) -> bool {
    value == GAP_MARKER || value > 1
}

fn symbols_match(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(&left, &right)| left == right || left == 4 || right == 4)
}

fn is_low_confidence_difference(
    view: &SegmentView,
    position: i64,
    chrom: &str,
    support: &AlleleDepthMap,
) -> bool {
    let Ok(offset) = usize::try_from(position - view.start) else {
        return false;
    };
    let Some(&token_position) = view.axis_token.get(offset) else {
        return false;
    };
    if token_position == u32::MAX {
        return false;
    }
    let token_position = token_position as usize;
    let Some(&source_symbol) = view.symbols.get(token_position) else {
        return false;
    };
    let Some(&uncertainty) = view.uncertainty.get(offset) else {
        return false;
    };
    let approximate_confidence = if uncertainty > 0.0 {
        (-10.0 * uncertainty.log10()) as u8
    } else {
        40
    };
    if approximate_confidence >= 20 {
        return false;
    }
    let Some(values) = support.get(chrom, position as u32) else {
        return false;
    };
    let item_support = values.get(source_symbol as usize).copied().unwrap_or(0);
    let total_support = values[5];
    if total_support == 0 {
        return false;
    }
    let fraction = item_support as f32 / total_support as f32;
    (fraction <= 0.02 || (item_support == 1 && total_support >= 10)) && approximate_confidence < 13
}

// Keep the parity-critical scoring inputs explicit in this hot path.
#[allow(clippy::too_many_arguments)]
fn compute_score(
    left_state: &[i16],
    right_state: &[i16],
    start: i64,
    left: &SegmentView,
    right: &SegmentView,
    chrom: &str,
    intrinsic_support: &AlleleDepthMap,
    config: &HaplotypeConfig,
) -> f32 {
    let count = left_state.len().min(right_state.len());

    let mut equal_span = 0usize;
    let mut indel_blocks = 0usize;
    let mut in_indel_run = false;
    for index in 0..count {
        if left_state[index] == right_state[index] {
            equal_span += 1;
            if is_indel_value(left_state[index]) {
                if !in_indel_run {
                    indel_blocks += 1;
                    in_indel_run = true;
                }
            } else {
                in_indel_run = false;
            }
        }
    }

    let mut shared = 0usize;
    let mut supported_shared = 0usize;
    for index in 0..count {
        if !is_snv_value(left_state[index]) || !is_snv_value(right_state[index]) {
            continue;
        }
        shared += 1;
        let position = start + index as i64;
        if position < 0 {
            continue;
        }
        let Some(left_symbol) = left.symbol_at(position) else {
            continue;
        };
        let Some(right_symbol) = right.symbol_at(position) else {
            continue;
        };
        if left_symbol != right_symbol || left_symbol == 4 {
            continue;
        }
        if intrinsic_support
            .get(chrom, position as u32)
            .is_some_and(|values| values[left_symbol as usize] > 0)
        {
            supported_shared += 1;
        }
    }

    let mut score = equal_span as f32;
    for value in config
        .score_array
        .iter()
        .take(shared.min(config.score_array.len()))
    {
        score += *value;
    }
    score += config.mean_read_length * (shared as f32 - supported_shared as f32) * 0.75;
    score += config.mean_read_length * 3.0 * indel_blocks as f32;
    score
}

// The comparison kernel shares the same explicit scoring context.
#[allow(clippy::too_many_arguments)]
fn compare_views(
    left: &SegmentView,
    right: &SegmentView,
    start: i64,
    end: i64,
    chrom: &str,
    support: &AlleleDepthMap,
    intrinsic_support: &AlleleDepthMap,
    config: &HaplotypeConfig,
) -> (PairStatus, Option<f32>) {
    let left_interval = left.slice_symbols(start, end);
    let right_interval = right.slice_symbols(start, end);
    if left_interval.is_empty() || right_interval.is_empty() {
        return (PairStatus::Indeterminate, None);
    }

    let left_state = left.slice_state(start, end);
    let right_state = right.slice_state(start, end);
    if symbols_match(left_interval, right_interval) {
        let raw_score = compute_score(
            left_state,
            right_state,
            start,
            left,
            right,
            chrom,
            intrinsic_support,
            config,
        );
        return (
            PairStatus::Compatible,
            Some(raw_score / (config.mean_read_length * 10.0)),
        );
    }

    if left_interval.len() != right_interval.len() {
        return (PairStatus::Blocked, None);
    }

    let length = (end - start) as usize;
    let count = left_state.len().min(right_state.len()).min(length);
    let mut positions = [0i64; 2];
    let mut difference_count = 0usize;
    for index in 0..count {
        if left_state[index] == right_state[index] {
            continue;
        }
        if is_indel_value(left_state[index]) || is_indel_value(right_state[index]) {
            return (PairStatus::Blocked, None);
        }
        if difference_count < 2 {
            positions[difference_count] = start + index as i64;
        }
        difference_count += 1;
        if difference_count >= 3 {
            return (PairStatus::Blocked, None);
        }
    }

    let shared_count = left_state.len().min(right_state.len());
    let mut distinct_count = 0usize;
    for index in 0..shared_count {
        if !is_snv_value(left_state[index]) || !is_snv_value(right_state[index]) {
            continue;
        }
        let position = start + index as i64;
        let Some(left_symbol) = left.symbol_at(position) else {
            continue;
        };
        let Some(right_symbol) = right.symbol_at(position) else {
            continue;
        };
        if left_symbol != 4 && right_symbol != 4 && left_symbol != right_symbol {
            if difference_count + distinct_count < 2 {
                positions[difference_count + distinct_count] = position;
            }
            distinct_count += 1;
            if difference_count + distinct_count >= 3 {
                return (PairStatus::Blocked, None);
            }
        }
    }

    if difference_count == 0 && distinct_count == 0 {
        return (PairStatus::Blocked, None);
    }
    let total = difference_count + distinct_count;
    debug_assert!(total <= 2);

    let mut accepted_count = 0usize;
    for &position in &positions[..total] {
        let index = (position - start) as usize;
        let left_value = left_state.get(index).copied().unwrap_or(1);
        let right_value = right_state.get(index).copied().unwrap_or(1);
        if is_indel_value(left_value) || is_indel_value(right_value) {
            return (PairStatus::Blocked, None);
        }
        let left_accepted = is_low_confidence_difference(left, position, chrom, support);
        let right_accepted = is_low_confidence_difference(right, position, chrom, support);
        if (left_accepted && is_snv_value(left_value))
            || (right_accepted && is_snv_value(right_value))
        {
            accepted_count += 1;
        } else {
            return (PairStatus::Blocked, None);
        }
    }

    let raw_score = compute_score(
        left_state,
        right_state,
        start,
        left,
        right,
        chrom,
        intrinsic_support,
        config,
    );
    let penalized = (raw_score - accepted_count as f32 * 20.0).max(0.0);
    (
        PairStatus::Compatible,
        Some(penalized / (config.mean_read_length * 10.0)),
    )
}

struct SweepItem {
    start: i64,
    end: i64,
    entity: u32,
    roles: u8,
}

fn enumerate_candidates(
    read_pair_map: &ReadPairMap,
    read_pairs: &[&ReadPair],
    stream_count: usize,
) -> Result<Vec<(u32, u32)>, String> {
    if read_pairs.len() >= (1usize << 30) {
        return Err("candidate keys support fewer than 2^30 entities".to_string());
    }

    let mut per_stream: Vec<Vec<SweepItem>> = (0..stream_count).map(|_| Vec::new()).collect();
    for read_pair in read_pairs {
        let second = read_pair
            .read2
            .as_ref()
            .ok_or_else(|| format!("incomplete read pair {}", read_pair.qname_idx))?;
        let first_stream = usize::try_from(read_pair.read1.tid())
            .map_err(|_| format!("unmapped first record for pair {}", read_pair.qname_idx))?;
        let second_stream = usize::try_from(second.tid())
            .map_err(|_| format!("unmapped second record for pair {}", read_pair.qname_idx))?;
        if first_stream >= stream_count || second_stream >= stream_count {
            return Err(format!(
                "reference stream is out of range for pair {}",
                read_pair.qname_idx
            ));
        }
        let roles = if first_stream == second_stream {
            ROLE_INDEX | ROLE_QUERY
        } else {
            ROLE_QUERY
        };
        for record in [&read_pair.read1, second] {
            per_stream[first_stream].push(SweepItem {
                start: record.reference_start(),
                end: record.reference_end(),
                entity: read_pair.qname_idx as u32,
                roles,
            });
        }
    }

    let mut events: Vec<u64> = per_stream
        .into_par_iter()
        .flat_map(|mut items| {
            items.sort_unstable_by_key(|item| (item.start, item.end, item.entity));
            let mut active: BinaryHeap<Reverse<(i64, usize)>> = BinaryHeap::new();
            let mut stream_events = Vec::new();
            for (slot, item) in items.iter().enumerate() {
                while let Some(&Reverse((end, _))) = active.peek() {
                    if end <= item.start {
                        active.pop();
                    } else {
                        break;
                    }
                }
                for &Reverse((_, other_slot)) in &active {
                    let other = &items[other_slot];
                    if other.entity == item.entity || other.start >= item.end {
                        continue;
                    }
                    let (low, high, low_item, high_item) = if item.entity < other.entity {
                        (item.entity, other.entity, item, other)
                    } else {
                        (other.entity, item.entity, other, item)
                    };
                    let mut flags = 0u64;
                    if low_item.roles & ROLE_QUERY != 0 && high_item.roles & ROLE_INDEX != 0 {
                        flags |= 1;
                    }
                    if high_item.roles & ROLE_QUERY != 0 && low_item.roles & ROLE_INDEX != 0 {
                        flags |= 2;
                    }
                    if flags != 0 {
                        stream_events.push((u64::from(low) << 34) | (u64::from(high) << 2) | flags);
                    }
                }
                active.push(Reverse((item.end, slot)));
            }
            stream_events
        })
        .collect();
    events.par_sort_unstable();

    let mut rank = vec![0u32; read_pairs.len()];
    for (order, (&id, _)) in read_pair_map.readpair_dict.iter().enumerate() {
        rank[id] = order as u32;
    }

    let mut candidates = Vec::new();
    let mut cursor = 0usize;
    while cursor < events.len() {
        let key = events[cursor] >> 2;
        let mut flags = events[cursor] & 3;
        cursor += 1;
        while cursor < events.len() && events[cursor] >> 2 == key {
            flags |= events[cursor] & 3;
            cursor += 1;
        }
        let low = (key >> 32) as u32;
        let high = (key & 0xffff_ffff) as u32;
        let oriented = match flags {
            1 => (low, high),
            2 => (high, low),
            _ if rank[low as usize] < rank[high as usize] => (low, high),
            _ => (high, low),
        };
        candidates.push(oriented);
    }
    Ok(candidates)
}

struct ChunkOutput {
    assignments: Vec<(u32, u32, f32)>,
    edges: Vec<(u32, u32, f32)>,
    blocked: usize,
}

// Pair classification is a tight loop over immutable run-level context.
#[allow(clippy::too_many_arguments)]
fn classify_pair(
    entity_views: &[EntityView],
    entity_id: u32,
    other_id: u32,
    stream_names: &[String],
    support: &AlleleDepthMap,
    intrinsic_support: &AlleleDepthMap,
    config: &HaplotypeConfig,
    output: &mut ChunkOutput,
) {
    let entity = &entity_views[entity_id as usize];
    let other = &entity_views[other_id as usize];
    let chrom = &stream_names[entity.segments[0].stream as usize];

    let mut overlaps: [(i64, i64, u8, u8); 4] = [(0, 0, 0, 0); 4];
    let mut overlap_count = 0usize;
    for left_slot in 0..2u8 {
        for right_slot in 0..2u8 {
            let left = &entity.segments[left_slot as usize];
            let right = &other.segments[right_slot as usize];
            let start = left.start.max(right.start);
            let end = left.end.min(right.end);
            if start < end {
                overlaps[overlap_count] = (start, end, left_slot, right_slot);
                overlap_count += 1;
            }
        }
    }
    if overlap_count == 0 {
        return;
    }

    let mut has_blocked = false;
    let mut pair_score: Option<f32> = None;
    let mut inspected: [(i64, i64); 4] = [(0, 0); 4];
    for (inspected_count, &(overlap_start, overlap_end, left_slot, right_slot)) in
        overlaps[..overlap_count].iter().enumerate()
    {
        let mut covering: [(i64, i64); 4] = [(0, 0); 4];
        let mut covering_count = 0usize;
        for &(start, end) in &inspected[..inspected_count] {
            if start <= overlap_end && end >= overlap_start {
                covering[covering_count] = (start, end);
                covering_count += 1;
            }
        }
        covering[..covering_count].sort_unstable_by_key(|&(start, _)| start);

        let mut uncovered: [(i64, i64); 5] = [(0, 0); 5];
        let mut uncovered_count = 0usize;
        let mut current = overlap_start;
        for &(start, end) in &covering[..covering_count] {
            if start > current {
                uncovered[uncovered_count] = (current, start);
                uncovered_count += 1;
            }
            current = current.max(end);
        }
        if current < overlap_end {
            uncovered[uncovered_count] = (current, overlap_end);
            uncovered_count += 1;
        }

        for &(start, end) in &uncovered[..uncovered_count] {
            let (status, score) = compare_views(
                &entity.segments[left_slot as usize],
                &other.segments[right_slot as usize],
                start,
                end,
                chrom,
                support,
                intrinsic_support,
                config,
            );
            if let Some(score) = score {
                let score = score.max(0.0);
                pair_score = Some(pair_score.map_or(score, |current| current + score));
            }
            if status == PairStatus::Blocked {
                has_blocked = true;
                break;
            }
        }
        inspected[inspected_count] = (overlap_start, overlap_end);
    }

    if has_blocked {
        output.assignments.push((entity_id, other_id, -1.0));
        output.assignments.push((other_id, entity_id, -1.0));
        output.blocked += 1;
    } else {
        let score = pair_score.unwrap_or(0.0).max(1.0e-4);
        output.assignments.push((entity_id, other_id, score));
        output.assignments.push((other_id, entity_id, score));
        output.edges.push((entity_id, other_id, score));
    }
}

fn build_entity(read_pair: &ReadPair) -> Result<BuiltEntity, String> {
    let second = read_pair
        .read2
        .as_ref()
        .ok_or_else(|| format!("incomplete read pair {}", read_pair.qname_idx))?;
    let first_id = get_read_id(&read_pair.read1).map_err(|error| error.to_string())?;
    let second_id = get_read_id(second).map_err(|error| error.to_string())?;

    let first_state = Arc::new(extract_hap_vector(&read_pair.read1));
    let first_uncertainty = Arc::new(extract_error_vector(&read_pair.read1));
    let (second_state, second_uncertainty) = if first_id == second_id {
        (Arc::clone(&first_state), Arc::clone(&first_uncertainty))
    } else {
        (
            Arc::new(extract_hap_vector(second)),
            Arc::new(extract_error_vector(second)),
        )
    };
    let first = SegmentView::build(&read_pair.read1, first_state, first_uncertainty)?;
    let second = SegmentView::build(second, second_state, second_uncertainty)?;
    Ok(BuiltEntity {
        view: EntityView {
            segments: [first, second],
        },
        first_id,
        second_id,
    })
}

pub(crate) fn build_phasing_graph(
    read_pair_map: &ReadPairMap,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_ad_map: &AlleleDepthMap,
    header: &rust_htslib::bam::HeaderView,
    config: &HaplotypeConfig,
) -> Result<PhasingGraphResult, Box<dyn std::error::Error>> {
    if allele_depth_map.is_empty() {
        return Ok(PhasingGraphResult::new());
    }

    let item_count = read_pair_map.readpair_dict.len();
    if u32::try_from(item_count).is_err() {
        return Err("graph has more vertices than compact indices can represent".into());
    }
    let mut result = PhasingGraphResult::new();
    let mut read_pairs = Vec::with_capacity(item_count);
    for entity_id in 0..item_count {
        let node = result.graph.add_node(());
        if node.index() != entity_id {
            return Err("graph node identifiers are not contiguous".into());
        }
        let read_pair = read_pair_map
            .readpair_dict
            .get(&entity_id)
            .ok_or_else(|| format!("missing read pair {entity_id}"))?;
        read_pairs.push(read_pair);
    }

    let precompute_start = Instant::now();
    let built: Vec<Result<BuiltEntity, String>> = read_pairs
        .par_iter()
        .map(|pair| build_entity(pair))
        .collect();
    let mut entity_views = Vec::with_capacity(item_count);
    for built in built {
        let built = built.map_err(|error| -> Box<dyn std::error::Error> { error.into() })?;
        result
            .read_hap_vectors
            .entry(built.first_id.clone())
            .or_insert_with(|| built.view.segments[0].state.as_ref().clone());
        result
            .read_error_vectors
            .entry(built.first_id.clone())
            .or_insert_with(|| built.view.segments[0].uncertainty.as_ref().clone());
        result
            .read_hap_vectors
            .entry(built.second_id.clone())
            .or_insert_with(|| built.view.segments[1].state.as_ref().clone());
        result
            .read_error_vectors
            .entry(built.second_id.clone())
            .or_insert_with(|| built.view.segments[1].uncertainty.as_ref().clone());
        result
            .node_read_ids
            .push((built.first_id, Some(built.second_id)));
        entity_views.push(built.view);
    }
    let precompute_time = precompute_start.elapsed();

    let stream_count = header.target_count() as usize;
    let stream_names: Vec<String> = (0..header.target_count())
        .map(|stream| String::from_utf8_lossy(header.tid2name(stream)).into_owned())
        .collect();
    let candidate_start = Instant::now();
    let candidates = enumerate_candidates(read_pair_map, &read_pairs, stream_count)
        .map_err(|error| -> Box<dyn std::error::Error> { error.into() })?;
    let candidate_time = candidate_start.elapsed();

    result.initialize_weight_store(item_count);
    let classify_start = Instant::now();
    let mut blocked_pairs = 0usize;
    const CHUNK: usize = 4096;
    const BATCH: usize = 4096 * CHUNK;
    for batch in candidates.chunks(BATCH) {
        let outputs: Vec<ChunkOutput> = batch
            .par_chunks(CHUNK)
            .map(|chunk| {
                let mut output = ChunkOutput {
                    assignments: Vec::with_capacity(chunk.len() * 2),
                    edges: Vec::new(),
                    blocked: 0,
                };
                for &(entity_id, other_id) in chunk {
                    classify_pair(
                        &entity_views,
                        entity_id,
                        other_id,
                        &stream_names,
                        allele_depth_map,
                        intrinsic_ad_map,
                        config,
                        &mut output,
                    );
                }
                output
            })
            .collect();
        for output in outputs {
            result.extend_weight_entries(output.assignments);
            for (left, right, weight) in output.edges {
                result.graph.add_edge(
                    NodeIndex::new(left as usize),
                    NodeIndex::new(right as usize),
                    weight,
                );
            }
            blocked_pairs += output.blocked;
        }
    }
    let classify_time = classify_start.elapsed();

    result.lowqual_qnames = read_pair_map.noisy_qnames.keys().cloned().collect();
    info!(
        "[build_phasing_graph] fast engine vertices={} candidates={} blocked={} edges={} assignments={} precompute_s={:.3} candidate_s={:.3} classify_s={:.3}",
        result.vertex_count(),
        candidates.len(),
        blocked_pairs,
        result.edge_count(),
        result.sparse_weight_entry_count(),
        precompute_time.as_secs_f64(),
        candidate_time.as_secs_f64(),
        classify_time.as_secs_f64(),
    );
    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use petgraph::visit::EdgeRef;
    use rust_htslib::bam::ext::BamRecordExtensions;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView, Record};

    use super::build_phasing_graph;
    use crate::graph_builder::build_phasing_graph_legacy;
    use crate::structs::{
        AlleleDepthMap, HaplotypeConfig, PhasingGraphResult, ReadPair, ReadPairMap,
        SortedVecIntervals,
    };

    struct Lcg(u64);

    impl Lcg {
        fn new(seed: u64) -> Self {
            Self(
                seed.wrapping_mul(0x9e37_79b9_7f4a_7c15)
                    .wrapping_add(0x1234_5678),
            )
        }

        fn next_u64(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0
        }

        fn below(&mut self, bound: u64) -> u64 {
            self.next_u64() % bound
        }

        fn chance(&mut self, percent: u64) -> bool {
            self.below(100) < percent
        }
    }

    fn make_header(stream_count: usize) -> HeaderView {
        let mut header = Header::new();
        for stream in 0..stream_count {
            let mut record = HeaderRecord::new(b"SQ");
            record.push_tag(b"SN", format!("stream_{stream}"));
            record.push_tag(b"LN", 10_000);
            header.push_record(&record);
        }
        HeaderView::from_header(&header)
    }

    fn query_len(cigar: &[Cigar]) -> usize {
        cigar
            .iter()
            .map(|operation| match operation {
                Cigar::Match(length)
                | Cigar::Ins(length)
                | Cigar::SoftClip(length)
                | Cigar::Equal(length)
                | Cigar::Diff(length) => *length as usize,
                _ => 0,
            })
            .sum()
    }

    fn random_record(rng: &mut Lcg, qname: &[u8], flags: u16, stream: i32, start: i64) -> Record {
        let mut cigar = Vec::new();
        if rng.chance(30) {
            cigar.push(Cigar::SoftClip(1 + rng.below(3) as u32));
        }
        cigar.push(Cigar::Equal(5 + rng.below(8) as u32));
        if rng.chance(35) {
            cigar.push(Cigar::Ins(1 + rng.below(2) as u32));
        }
        cigar.push(Cigar::Equal(4 + rng.below(7) as u32));
        if rng.chance(45) {
            cigar.push(Cigar::Diff(1 + rng.below(2) as u32));
        }
        if rng.chance(25) {
            cigar.push(Cigar::Del(1 + rng.below(2) as u32));
        }
        if rng.chance(15) {
            cigar.push(Cigar::RefSkip(1 + rng.below(2) as u32));
        }
        cigar.push(Cigar::Equal(6 + rng.below(8) as u32));
        if rng.chance(25) {
            cigar.push(Cigar::SoftClip(1 + rng.below(2) as u32));
        }

        let length = query_len(&cigar);
        let alphabet = [b'A', b'T', b'C', b'G', b'N'];
        let sequence: Vec<u8> = (0..length)
            .map(|_| alphabet[rng.below(alphabet.len() as u64) as usize])
            .collect();
        let quality: Vec<u8> = (0..length)
            .map(|_| if rng.chance(20) { 8 } else { 30 })
            .collect();
        let mut record = Record::new();
        record.set(qname, Some(&CigarString(cigar)), &sequence, &quality);
        record.set_flags(flags);
        record.set_tid(stream);
        record.set_pos(start);
        record.set_mapq(60);
        record
    }

    fn fixture(seed: u64, entity_count: usize, stream_count: usize) -> (ReadPairMap, HeaderView) {
        let mut rng = Lcg::new(seed);
        let header = make_header(stream_count);
        let mut map = ReadPairMap::new();
        for stream in 0..stream_count {
            map.interval_trees
                .insert(format!("stream_{stream}"), SortedVecIntervals::new());
        }

        for entity_id in 0..entity_count {
            let qname = format!("entity_{entity_id:06}");
            let first_stream = rng.below(stream_count as u64) as i32;
            let second_stream = if stream_count > 1 && rng.chance(20) {
                (first_stream + 1) % stream_count as i32
            } else {
                first_stream
            };
            let first_start = 100 + rng.below(35) as i64;
            let second_start = 100 + rng.below(35) as i64;
            let first = random_record(&mut rng, qname.as_bytes(), 65, first_stream, first_start);
            let second_flags = if rng.chance(20) { 65 } else { 129 };
            let second = random_record(
                &mut rng,
                qname.as_bytes(),
                second_flags,
                second_stream,
                second_start,
            );

            for record in [&first, &second] {
                let chrom = format!("stream_{}", record.tid());
                map.interval_trees
                    .get_mut(&chrom)
                    .unwrap()
                    .add_interval(record.reference_start(), record.reference_end(), entity_id)
                    .unwrap();
            }
            map.qname_to_idx.insert(qname.clone(), entity_id);
            map.idx_to_qname.insert(entity_id, qname.clone());
            map.readpair_dict.insert(
                entity_id,
                ReadPair::new_complete(first, second, qname, entity_id),
            );
        }
        for intervals in map.interval_trees.values_mut() {
            intervals.finalize().unwrap();
        }
        (map, header)
    }

    fn support(stream_count: usize) -> (AlleleDepthMap, AlleleDepthMap) {
        let mut primary = AlleleDepthMap::new();
        let mut secondary = AlleleDepthMap::new();
        for stream in 0..stream_count {
            let chrom = format!("stream_{stream}");
            for position in 0..400u32 {
                primary.insert(&chrom, position, [1, 2, 3, 4, 0, 20]);
                if position % 3 == 0 {
                    secondary.insert(&chrom, position, [2, 0, 1, 0, 0, 20]);
                }
            }
        }
        (primary, secondary)
    }

    fn sorted_assignments(result: &PhasingGraphResult) -> Vec<(u32, u32, u32)> {
        let mut assignments: Vec<_> = result
            .sparse_weight_entries()
            .iter()
            .map(|&(left, right, weight)| (left, right, weight.to_bits()))
            .collect();
        assignments.sort_unstable();
        assignments
    }

    fn sorted_edges(result: &PhasingGraphResult) -> Vec<(usize, usize, u32)> {
        let mut edges: Vec<_> = result
            .graph
            .edge_references()
            .map(|edge| {
                let left = edge.source().index();
                let right = edge.target().index();
                (left.min(right), left.max(right), edge.weight().to_bits())
            })
            .collect();
        edges.sort_unstable();
        edges
    }

    fn normalized_map<T: Clone>(map: &ahash::AHashMap<String, Vec<T>>) -> HashMap<String, Vec<T>> {
        map.iter()
            .map(|(key, values)| (key.clone(), values.clone()))
            .collect()
    }

    fn assert_equivalent(seed: u64, entity_count: usize, stream_count: usize) {
        let (map, header) = fixture(seed, entity_count, stream_count);
        let (primary, secondary) = support(stream_count);
        let config = HaplotypeConfig::new(100.0);
        let legacy =
            build_phasing_graph_legacy(&map, &primary, &secondary, &header, &config).unwrap();
        let fast = build_phasing_graph(&map, &primary, &secondary, &header, &config).unwrap();

        assert_eq!(fast.vertex_count(), legacy.vertex_count(), "seed={seed}");
        assert_eq!(fast.edge_count(), legacy.edge_count(), "seed={seed}");
        assert_eq!(
            fast.content_digests(),
            legacy.content_digests(),
            "seed={seed}"
        );
        assert_eq!(
            sorted_assignments(&fast),
            sorted_assignments(&legacy),
            "seed={seed}"
        );
        assert_eq!(sorted_edges(&fast), sorted_edges(&legacy), "seed={seed}");
        assert_eq!(fast.node_read_ids, legacy.node_read_ids, "seed={seed}");
        assert_eq!(fast.lowqual_qnames, legacy.lowqual_qnames, "seed={seed}");
        assert_eq!(
            normalized_map(&fast.read_hap_vectors),
            normalized_map(&legacy.read_hap_vectors),
            "seed={seed}"
        );
        assert_eq!(
            normalized_map(&fast.read_error_vectors),
            normalized_map(&legacy.read_error_vectors),
            "seed={seed}"
        );
    }

    #[test]
    fn randomized_records_match_legacy_exactly() {
        for seed in 0..16 {
            assert_equivalent(seed, 32, 3);
        }
    }

    #[test]
    fn dense_single_stream_records_match_legacy_exactly() {
        for seed in 100..104 {
            assert_equivalent(seed, 48, 1);
        }
    }
}

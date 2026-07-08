//! Faithful Rust ports of the numba/scipy kernels that drive Greedy-Clique-Expansion.
//!
//! These mirror `fp_control/numba_operators.py` and the `@njit` helpers at the top of
//! `fp_control/gce_algorithm.py`. Each kernel is intentionally a 1:1 translation (same
//! tie-breaking, same `-1`/`0` semantics) so the partition matches Python exactly.
//!
//! Weight-matrix value convention (from `build_phasing_graph`):
//!   * `-1.0` : the two read-pairs carry conflicting variants -> may NOT share a haplotype.
//!   * `0.0`  : no informative overlap (dropped from the sparse representation).
//!   * `0..1` : edge weight (shared-variant evidence).

use ndarray::Array2;

/// Compressed Sparse Row view of a square weight matrix.
///
/// Built to match `scipy.sparse.csr_matrix(dense)`: every structurally **non-zero**
/// entry is stored (so `-1.0` incompatibility markers are kept; exact `0.0` is dropped),
/// with each row's column indices in ascending order.
#[derive(Clone, Debug)]
pub struct Csr {
    pub data: Vec<f32>,
    pub indices: Vec<i32>,
    pub indptr: Vec<i32>,
    pub size: usize, // square: n_rows == n_cols
}

impl Csr {
    /// Build CSR from a dense square matrix, scipy-style (drop exact zeros).
    pub fn from_dense(m: &Array2<f32>) -> Self {
        let size = m.nrows();
        assert_eq!(size, m.ncols(), "weight matrix must be square");
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = Vec::with_capacity(size + 1);
        indptr.push(0i32);
        for i in 0..size {
            for j in 0..size {
                let v = m[[i, j]];
                if v != 0.0 {
                    data.push(v);
                    indices.push(j as i32);
                }
            }
            indptr.push(data.len() as i32);
        }
        Csr {
            data,
            indices,
            indptr,
            size,
        }
    }

    /// Build CSR directly from matrix entries, without materializing the dense matrix.
    ///
    /// Entries with exact value `0.0` are dropped, matching [`Csr::from_dense`]. Duplicate
    /// `(row, col)` entries keep the last value in input order, matching repeated assignment
    /// into a dense matrix before conversion. Stored columns are sorted ascending within rows.
    pub fn from_entries<I>(size: usize, entries: I) -> Self
    where
        I: IntoIterator<Item = (usize, usize, f32)>,
    {
        let mut entries: Vec<(usize, usize, usize, f32)> = entries
            .into_iter()
            .enumerate()
            .map(|(order, (row, col, value))| {
                assert!(
                    row < size,
                    "CSR row index {row} out of bounds for size {size}"
                );
                assert!(
                    col < size,
                    "CSR col index {col} out of bounds for size {size}"
                );
                (row, col, order, value)
            })
            .collect();

        entries.sort_unstable_by_key(|&(row, col, order, _)| (row, col, order));

        let mut data = Vec::with_capacity(entries.len());
        let mut indices = Vec::with_capacity(entries.len());
        let mut indptr = Vec::with_capacity(size + 1);
        indptr.push(0i32);

        let mut pos = 0usize;
        for row in 0..size {
            let row_start = pos;
            while pos < entries.len() && entries[pos].0 == row {
                let col = entries[pos].1;
                let mut value = entries[pos].3;
                pos += 1;
                while pos < entries.len() && entries[pos].0 == row && entries[pos].1 == col {
                    value = entries[pos].3;
                    pos += 1;
                }
                if value != 0.0 {
                    data.push(value);
                    indices.push(col as i32);
                }
            }
            debug_assert!(
                pos >= row_start,
                "entry cursor should not move backwards while building CSR"
            );
            indptr.push(data.len() as i32);
        }

        Csr {
            data,
            indices,
            indptr,
            size,
        }
    }

    #[inline]
    pub fn row(&self, i: usize) -> (&[f32], &[i32]) {
        let s = self.indptr[i] as usize;
        let e = self.indptr[i + 1] as usize;
        (&self.data[s..e], &self.indices[s..e])
    }

    /// Dense-matrix-compatible lookup. Missing structural entries are exact `0.0`.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f32 {
        assert!(
            row < self.size,
            "CSR row index {row} out of bounds for size {}",
            self.size
        );
        assert!(
            col < self.size,
            "CSR col index {col} out of bounds for size {}",
            self.size
        );
        let (rd, rc) = self.row(row);
        match rc.binary_search(&(col as i32)) {
            Ok(k) => rd[k],
            Err(_) => 0.0,
        }
    }

    /// `np.all(weight_matrix <= threshold, axis=1)` without materializing dense rows.
    ///
    /// Since absent entries are zero, they are always `<= threshold` for the positive
    /// thresholds used by phasing. Stored `-1.0` incompatibilities also count as small.
    pub fn small_row_mask(&self, threshold: f32) -> Vec<bool> {
        (0..self.size)
            .map(|row| {
                let (rd, _) = self.row(row);
                rd.iter().all(|&v| v <= threshold)
            })
            .collect()
    }

    /// Reindex to the sub-matrix over rows/cols where `mask` is true, compacting indices.
    /// Mirrors `sparse_matrix[mask, :][:, mask]`.
    #[allow(clippy::needless_range_loop)] // `i` both filters the mask and selects the CSR row
    pub fn select(&self, mask: &[bool]) -> Csr {
        debug_assert_eq!(mask.len(), self.size);
        // new column index for each old column that survives (others = -1)
        let mut remap = vec![-1i32; self.size];
        let mut new_size = 0usize;
        for (old, &keep) in mask.iter().enumerate() {
            if keep {
                remap[old] = new_size as i32;
                new_size += 1;
            }
        }
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = Vec::with_capacity(new_size + 1);
        indptr.push(0i32);
        for i in 0..self.size {
            if !mask[i] {
                continue;
            }
            let (rd, rc) = self.row(i);
            for (k, &col) in rc.iter().enumerate() {
                let nc = remap[col as usize];
                if nc >= 0 {
                    data.push(rd[k]);
                    indices.push(nc);
                }
            }
            indptr.push(data.len() as i32);
        }
        Csr {
            data,
            indices,
            indptr,
            size: new_size,
        }
    }
}

/// Extract the dense sub-matrix over rows/cols where `mask` is true (`np.ix_` semantics).
pub fn dense_submatrix(m: &Array2<f32>, mask: &[bool]) -> Array2<f32> {
    let idx: Vec<usize> = mask
        .iter()
        .enumerate()
        .filter(|(_, &b)| b)
        .map(|(i, _)| i)
        .collect();
    let n = idx.len();
    let mut out = Array2::<f32>::zeros((n, n));
    for (ri, &oi) in idx.iter().enumerate() {
        for (ci, &oj) in idx.iter().enumerate() {
            out[[ri, ci]] = m[[oi, oj]];
        }
    }
    out
}

/// `numba_max_idx_mem`: argmax with `>=` tie-breaking (ties resolve to the LAST index).
/// Returns `(index, value)`; `index = -1` if no element is considered.
pub fn max_idx_mem(data: &[f32], index_mask: Option<&[bool]>) -> (i32, f32) {
    let mut max_val = f32::NEG_INFINITY;
    let mut max_idx = -1i32;
    match index_mask {
        None => {
            for (i, &d) in data.iter().enumerate() {
                if d >= max_val {
                    max_idx = i as i32;
                    max_val = d;
                }
            }
        }
        Some(mask) => {
            for (i, (&d, &keep)) in data.iter().zip(mask).enumerate() {
                if keep && d >= max_val {
                    max_idx = i as i32;
                    max_val = d;
                }
            }
        }
    }
    (max_idx, max_val)
}

/// `efficient_mask`: for one CSR row, mark columns whose stored value equals `mask_value`
/// (default `-1.0`) as `false`; everything else stays `true`.
pub fn efficient_mask(
    row_data: &[f32],
    row_cols: &[i32],
    size: usize,
    mask_value: f32,
) -> Vec<bool> {
    let mut mask = vec![true; size];
    for (k, &v) in row_data.iter().enumerate() {
        if v == mask_value {
            mask[row_cols[k] as usize] = false;
        }
    }
    mask
}

/// `efficient_row_max`: best edge from one row to an in-mask, non-`-1` column.
/// Returns `(column_index, value)`, or `(-1, -2.0)` when no valid edge exists.
pub fn efficient_row_max(row_data: &[f32], row_cols: &[i32], index_mask: &[bool]) -> (i32, f32) {
    // valid_mask = (data != -1) AND index_mask[col]
    let mut valid = vec![false; row_data.len()];
    let mut any = false;
    for k in 0..row_data.len() {
        if row_data[k] != -1.0 && index_mask[row_cols[k] as usize] {
            valid[k] = true;
            any = true;
        }
    }
    if any {
        let (short, val) = max_idx_mem(row_data, Some(&valid));
        let col = row_cols[short as usize];
        (col, val)
    } else {
        (-1, -2.0)
    }
}

/// `row_wise_max_with_mask_sparse`: per-row max over valid (non-`-1`, in-mask, non-self)
/// entries. Rows outside `index_mask` and rows with no valid entry yield `0.0`.
#[allow(clippy::needless_range_loop)] // `i` indexes the mask, selects the CSR row, and excludes self
pub fn row_wise_max_with_mask(csr: &Csr, index_mask: &[bool], mask_value: f32) -> Vec<f32> {
    let mut out = vec![0.0f32; csr.size];
    for i in 0..csr.size {
        if !index_mask[i] {
            continue;
        }
        let (rd, rc) = csr.row(i);
        let mut best = f32::NEG_INFINITY;
        let mut found = false;
        for k in 0..rd.len() {
            let col = rc[k] as usize;
            // exclude self (col == i) and -1 markers and out-of-mask columns
            if col != i && rd[k] != mask_value && index_mask[col] {
                if rd[k] > best {
                    best = rd[k];
                }
                found = true;
            }
        }
        if found {
            out[i] = best;
        }
    }
    out
}

/// `reverse_boolean_mask`: all-true, then set the given indices to false.
pub fn reverse_boolean_mask(total: usize, indices: &[i32]) -> Vec<bool> {
    let mut mask = vec![true; total];
    for &idx in indices {
        mask[idx as usize] = false;
    }
    mask
}

/// `boolean_mask`: all-false, then set the given indices to true.
pub fn boolean_mask_from_indices(total: usize, indices: &[i32]) -> Vec<bool> {
    let mut mask = vec![false; total];
    for &idx in indices {
        mask[idx as usize] = true;
    }
    mask
}

/// Elementwise AND of two boolean slices.
pub fn and_masks(a: &[bool], b: &[bool]) -> Vec<bool> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x && y).collect()
}

/// Count of true entries (`numba_sum` over a bool array).
#[inline]
pub fn count_true(mask: &[bool]) -> usize {
    mask.iter().filter(|&&b| b).count()
}

/// `apply_index_mask`: the indices (0..size) where `mask` is true, as i32.
pub fn apply_index_mask(mask: &[bool]) -> Vec<i32> {
    mask.iter()
        .enumerate()
        .filter(|(_, &b)| b)
        .map(|(i, _)| i as i32)
        .collect()
}

/// `numba_isin` over a 0..size arange: true where the index is in `set`.
pub fn isin_arange(size: usize, set: &std::collections::HashSet<i32>) -> Vec<bool> {
    (0..size).map(|i| set.contains(&(i as i32))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn sample_matrix() -> Array2<f32> {
        // 4x4 symmetric. -1 = incompatible, 0 = no overlap, positive = weight.
        array![
            [0.0, 0.8, -1.0, 0.0],
            [0.8, 0.0, 0.3, 0.5],
            [-1.0, 0.3, 0.0, 0.9],
            [0.0, 0.5, 0.9, 0.0]
        ]
    }

    #[test]
    fn csr_from_dense_drops_zeros_keeps_neg1() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        // row 0: cols 1 (0.8) and 2 (-1.0) are stored; col 0,3 are zero -> dropped
        let (rd, rc) = csr.row(0);
        assert_eq!(rc, &[1, 2]);
        assert_eq!(rd, &[0.8, -1.0]);
        // round-trip a couple of entries
        let (rd3, rc3) = csr.row(3);
        assert_eq!(rc3, &[1, 2]);
        assert_eq!(rd3, &[0.5, 0.9]);
    }

    #[test]
    fn csr_from_entries_matches_dense_and_keeps_last_duplicate() {
        let mut m = sample_matrix();
        let mut entries = Vec::new();
        for i in 0..m.nrows() {
            for j in 0..m.ncols() {
                entries.push((i, j, m[[i, j]]));
            }
        }
        entries.push((1, 3, 0.2));
        entries.push((1, 3, 0.5)); // final dense assignment wins
        entries.push((2, 3, 0.0)); // final zero assignment deletes a stored non-zero
        m[[2, 3]] = 0.0;

        let csr = Csr::from_entries(m.nrows(), entries);
        let dense_csr = Csr::from_dense(&m);
        assert_eq!(csr.data, dense_csr.data);
        assert_eq!(csr.indices, dense_csr.indices);
        assert_eq!(csr.indptr, dense_csr.indptr);
        assert_eq!(csr.size, dense_csr.size);
    }

    #[test]
    fn csr_get_matches_dense_semantics() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        assert_eq!(csr.get(0, 1), 0.8);
        assert_eq!(csr.get(0, 2), -1.0);
        assert_eq!(csr.get(0, 3), 0.0);
    }

    #[test]
    fn csr_small_row_mask_matches_dense_rule() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        assert_eq!(csr.small_row_mask(0.1), vec![false, false, false, false]);

        let only_weak = Csr::from_entries(3, [(0, 1, 0.05), (1, 0, 0.05), (1, 2, -1.0)]);
        assert_eq!(only_weak.small_row_mask(0.1), vec![true, true, true]);
    }

    #[test]
    fn max_idx_mem_ties_pick_last() {
        let data = [0.5f32, 0.9, 0.9, 0.2];
        let (idx, val) = max_idx_mem(&data, None);
        assert_eq!(idx, 2); // last of the tied 0.9s
        assert_eq!(val, 0.9);
        // with a mask excluding index 2, the earlier 0.9 wins
        let mask = [true, true, false, true];
        let (idx2, val2) = max_idx_mem(&data, Some(&mask));
        assert_eq!(idx2, 1);
        assert_eq!(val2, 0.9);
    }

    #[test]
    fn efficient_mask_marks_neg1_columns_false() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        let (rd, rc) = csr.row(0);
        let mask = efficient_mask(rd, rc, csr.size, -1.0);
        // col 2 is -1 in row 0 -> false; others true
        assert_eq!(mask, vec![true, true, false, true]);
    }

    #[test]
    fn efficient_row_max_skips_neg1_and_out_of_mask() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        let (rd, rc) = csr.row(2); // cols: 1 (0.3), 3 (0.9); col0 is -1
        let full = vec![true; 4];
        let (col, val) = efficient_row_max(rd, rc, &full);
        assert_eq!(col, 3);
        assert_eq!(val, 0.9);
        // mask out col 3 -> falls back to 0.3 at col 1
        let mut mask = vec![true; 4];
        mask[3] = false;
        let (col2, val2) = efficient_row_max(rd, rc, &mask);
        assert_eq!(col2, 1);
        assert_eq!(val2, 0.3);
        // mask out everything valid -> sentinel
        let none = vec![false; 4];
        assert_eq!(efficient_row_max(rd, rc, &none), (-1, -2.0));
    }

    #[test]
    fn row_wise_max_excludes_self_and_neg1() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        let mask = vec![true; 4];
        let rm = row_wise_max_with_mask(&csr, &mask, -1.0);
        // row0: only valid edge is 0.8 (col2 is -1). row3: max(0.5,0.9)=0.9
        assert_eq!(rm[0], 0.8);
        assert_eq!(rm[1], 0.8);
        assert_eq!(rm[2], 0.9);
        assert_eq!(rm[3], 0.9);
    }

    #[test]
    fn select_submatrix_reindexes() {
        let m = sample_matrix();
        let csr = Csr::from_dense(&m);
        // keep vertices {1,2,3} -> new indices 0,1,2
        let mask = vec![false, true, true, true];
        let sub = csr.select(&mask);
        assert_eq!(sub.size, 3);
        // old row1 (cols 0:0.8 dropped, 2:0.3, 3:0.5) -> new row0 cols (1:0.3, 2:0.5)
        let (rd, rc) = sub.row(0);
        assert_eq!(rc, &[1, 2]);
        assert_eq!(rd, &[0.3, 0.5]);
    }

    #[test]
    fn dense_submatrix_matches_ix() {
        let m = sample_matrix();
        let mask = vec![false, true, true, true];
        let sub = dense_submatrix(&m, &mask);
        assert_eq!(sub.shape(), &[3, 3]);
        assert_eq!(sub[[0, 1]], 0.3); // old [1,2]
        assert_eq!(sub[[1, 2]], 0.9); // old [2,3]
    }

    #[test]
    fn mask_helpers() {
        assert_eq!(
            reverse_boolean_mask(4, &[1, 3]),
            vec![true, false, true, false]
        );
        assert_eq!(
            boolean_mask_from_indices(4, &[1, 3]),
            vec![false, true, false, true]
        );
        assert_eq!(apply_index_mask(&[false, true, true, false]), vec![1, 2]);
        assert_eq!(count_true(&[true, false, true, true]), 3);
        assert_eq!(
            and_masks(&[true, true, false], &[true, false, false]),
            vec![true, false, false]
        );
    }
}

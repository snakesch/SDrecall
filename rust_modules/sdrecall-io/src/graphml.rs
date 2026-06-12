//! GraphML read/write — the ONE graphml unit.
//!
//! [`write_graphml`] wraps `petgraph_graphml::GraphMl`, projecting each node /
//! edge weight to a list of `(attr_name, value)` string pairs (graph-tool's
//! GraphML attributes are all string-typed, so this matches the shape the Python
//! pipeline writes via `graph.save(...)`). [`read_graphml`] parses the file back
//! with `quick-xml` event reading — `petgraph-graphml` has no reader — into a
//! `petgraph::Graph` whose node / edge weights are `AHashMap<String, String>`
//! attribute maps.
//!
//! ## Format produced by `petgraph-graphml` (what the reader must accept)
//!
//! ```xml
//! <graphml xmlns="...">
//!   <key id="ATTR" for="node" attr.name="ATTR" attr.type="string"/>
//!   <graph edgedefault="undirected|directed">
//!     <node id="n0"><data key="ATTR">value</data></node>
//!     <edge id="e0" source="n0" target="n1"><data key="W">0.5</data></edge>
//!   </graph>
//! </graphml>
//! ```
//!
//! Node ids are `n<petgraph-index>`; the reader maps each `n<k>` to a freshly
//! added petgraph node, preserving index order (nodes are emitted in index
//! order), so a Rust write → read round-trip is index-stable.

use ahash::AHashMap;
use petgraph::graph::Graph;
use petgraph::visit::{GraphProp, IntoEdgeReferences, IntoNodeReferences, NodeIndexable};
use petgraph::{Directed, EdgeType, Undirected};
use quick_xml::events::Event;
use quick_xml::{Reader, XmlVersion};
use sdrecall_utils::{Result, SdError};
use std::borrow::Cow;
use std::io::BufWriter;
use std::path::Path;

/// A read-back GraphML graph: node and edge weights are string→string attribute
/// maps (`<data key="...">value</data>` pairs). `Ty` is the edge direction
/// (`Directed` from [`read_graphml`], `Undirected` from
/// [`read_graphml_undirected`]). A named alias keeps the public signatures
/// readable (and satisfies `clippy::type_complexity`).
pub type AttrGraph<Ty> = Graph<AHashMap<String, String>, AHashMap<String, String>, Ty>;

/// Write a petgraph graph to GraphML at `path`.
///
/// `node_attrs` / `edge_attrs` project a weight to its `(attr_name, value)`
/// string pairs (empty list = no attributes for that element). The graph's
/// directedness is encoded in `edgedefault`. The `EdgeType` bound keeps this
/// generic over both `Graph<_,_, Directed>` and `Graph<_,_, Undirected>`.
pub fn write_graphml<N, E, Ty>(
    path: &Path,
    g: &Graph<N, E, Ty>,
    node_attrs: impl Fn(&N) -> Vec<(String, String)> + 'static,
    edge_attrs: impl Fn(&E) -> Vec<(String, String)> + 'static,
) -> Result<()>
where
    Ty: EdgeType,
    for<'a> &'a Graph<N, E, Ty>:
        GraphProp + IntoNodeReferences<NodeWeight = N> + IntoEdgeReferences<EdgeWeight = E> + NodeIndexable,
{
    // petgraph-graphml's closures must return `Vec<(Cow<'static,str>, Cow<'a,str>)>`.
    // The attr NAME is owned (static-ified via Cow::Owned); the VALUE borrows the
    // produced String for the call's lifetime.
    let to_cow = |pairs: Vec<(String, String)>| -> Vec<(Cow<'static, str>, Cow<'_, str>)> {
        pairs
            .into_iter()
            .map(|(k, v)| (Cow::Owned(k), Cow::Owned(v)))
            .collect()
    };

    let ml = petgraph_graphml::GraphMl::new(g)
        .pretty_print(true)
        .export_node_weights(Box::new(move |w| to_cow(node_attrs(w))))
        .export_edge_weights(Box::new(move |w| to_cow(edge_attrs(w))));

    let file = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    ml.to_writer(BufWriter::new(file))
        .map_err(|e| SdError::GraphMl(format!("write {}: {e}", path.display())))?;
    Ok(())
}

/// Parse a GraphML file into a petgraph graph with string attribute maps on each
/// node and edge.
///
/// Returns a `Graph<.., .., Directed>` regardless of `edgedefault` (an undirected
/// GraphML simply yields edges in the listed order; callers that need the
/// undirected type read via the dedicated helper below if added). Nodes are added
/// in document order and `n<id>` → node-index is recorded so edge `source`/
/// `target` references resolve.
pub fn read_graphml(path: &Path) -> Result<AttrGraph<Directed>> {
    read_graphml_impl::<Directed>(path)
}

/// Same as [`read_graphml`] but yields an explicitly undirected graph (for
/// callers that round-trip an `edgedefault="undirected"` graph and want the type
/// to reflect it). The parser logic is shared — this is the only place the type
/// parameter differs.
pub fn read_graphml_undirected(path: &Path) -> Result<AttrGraph<Undirected>> {
    read_graphml_impl::<Undirected>(path)
}

fn read_graphml_impl<Ty: EdgeType>(path: &Path) -> Result<AttrGraph<Ty>> {
    let xml = std::fs::read_to_string(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut reader = Reader::from_str(&xml);
    reader.config_mut().trim_text(true);

    let mut graph: AttrGraph<Ty> = Graph::with_capacity(0, 0);
    // GraphML node id (e.g. "n0") → petgraph node index.
    let mut id_to_node: AHashMap<String, petgraph::graph::NodeIndex> = AHashMap::new();

    // Parser state for the element currently being built.
    enum Elem {
        None,
        Node {
            id: String,
            attrs: AHashMap<String, String>,
        },
        Edge {
            source: String,
            target: String,
            attrs: AHashMap<String, String>,
        },
    }
    let mut cur = Elem::None;
    // The <data key="..."> currently open, awaiting its text.
    let mut cur_data_key: Option<String> = None;

    let attr_str = |e: &quick_xml::events::BytesStart, name: &[u8]| -> Result<Option<String>> {
        for a in e.attributes().with_checks(false) {
            let a = a.map_err(|err| SdError::GraphMl(format!("attribute parse: {err}")))?;
            if a.key.as_ref() == name {
                let v = a
                    .normalized_value(XmlVersion::Implicit1_0)
                    .map_err(|err| SdError::GraphMl(format!("attribute unescape: {err}")))?;
                return Ok(Some(v.into_owned()));
            }
        }
        Ok(None)
    };

    // Resolve (or create) the node index for a GraphML node id.
    let node_for = |id: &str,
                    graph: &mut AttrGraph<Ty>,
                    map: &mut AHashMap<String, petgraph::graph::NodeIndex>|
     -> petgraph::graph::NodeIndex {
        if let Some(&n) = map.get(id) {
            n
        } else {
            let n = graph.add_node(AHashMap::new());
            map.insert(id.to_string(), n);
            n
        }
    };

    loop {
        match reader
            .read_event()
            .map_err(|e| SdError::GraphMl(format!("read event: {e}")))?
        {
            Event::Eof => break,
            // A normal opening tag: open a node/edge element or arm `cur_data_key`.
            Event::Start(e) => match e.name().as_ref() {
                b"node" => {
                    let id = attr_str(&e, b"id")?
                        .ok_or_else(|| SdError::GraphMl("node without id".into()))?;
                    cur = Elem::Node {
                        id,
                        attrs: AHashMap::new(),
                    };
                }
                b"edge" => {
                    let source = attr_str(&e, b"source")?
                        .ok_or_else(|| SdError::GraphMl("edge without source".into()))?;
                    let target = attr_str(&e, b"target")?
                        .ok_or_else(|| SdError::GraphMl("edge without target".into()))?;
                    cur = Elem::Edge {
                        source,
                        target,
                        attrs: AHashMap::new(),
                    };
                }
                b"data" => {
                    cur_data_key = attr_str(&e, b"key")?;
                }
                _ => {}
            },
            // A self-closing tag (`<node id="n0"/>`, `<edge .../>` with no data,
            // or an empty `<data .../>`). No matching `End` follows, so finalize
            // here. graph-tool / petgraph-graphml emit attr-less nodes this way.
            Event::Empty(e) => match e.name().as_ref() {
                b"node" => {
                    let id = attr_str(&e, b"id")?
                        .ok_or_else(|| SdError::GraphMl("node without id".into()))?;
                    node_for(&id, &mut graph, &mut id_to_node);
                }
                b"edge" => {
                    let source = attr_str(&e, b"source")?
                        .ok_or_else(|| SdError::GraphMl("edge without source".into()))?;
                    let target = attr_str(&e, b"target")?
                        .ok_or_else(|| SdError::GraphMl("edge without target".into()))?;
                    let s = node_for(&source, &mut graph, &mut id_to_node);
                    let t = node_for(&target, &mut graph, &mut id_to_node);
                    graph.add_edge(s, t, AHashMap::new());
                }
                // `<data key="K"/>` — an empty-valued attribute on the open element.
                b"data" => {
                    if let Some(key) = attr_str(&e, b"key")? {
                        match &mut cur {
                            Elem::Node { attrs, .. } | Elem::Edge { attrs, .. } => {
                                attrs.insert(key, String::new());
                            }
                            Elem::None => {}
                        }
                    }
                }
                _ => {}
            },
            Event::Text(t) => {
                if let Some(key) = cur_data_key.take() {
                    let decoded = t
                        .decode()
                        .map_err(|e| SdError::GraphMl(format!("text decode: {e}")))?;
                    let val = quick_xml::escape::unescape(&decoded)
                        .map_err(|e| SdError::GraphMl(format!("text unescape: {e}")))?
                        .into_owned();
                    match &mut cur {
                        Elem::Node { attrs, .. } | Elem::Edge { attrs, .. } => {
                            attrs.insert(key, val);
                        }
                        Elem::None => {}
                    }
                }
            }
            Event::End(e) => match e.name().as_ref() {
                b"data" => {
                    // A <data></data> with no text leaves cur_data_key set; an
                    // empty value is recorded as "".
                    if let Some(key) = cur_data_key.take() {
                        match &mut cur {
                            Elem::Node { attrs, .. } | Elem::Edge { attrs, .. } => {
                                attrs.insert(key, String::new());
                            }
                            Elem::None => {}
                        }
                    }
                }
                b"node" => {
                    if let Elem::Node { id, attrs } = std::mem::replace(&mut cur, Elem::None) {
                        let n = node_for(&id, &mut graph, &mut id_to_node);
                        graph[n] = attrs;
                    }
                }
                b"edge" => {
                    if let Elem::Edge {
                        source,
                        target,
                        attrs,
                    } = std::mem::replace(&mut cur, Elem::None)
                    {
                        let s = node_for(&source, &mut graph, &mut id_to_node);
                        let t = node_for(&target, &mut graph, &mut id_to_node);
                        graph.add_edge(s, t, attrs);
                    }
                }
                _ => {}
            },
            _ => {}
        }
    }

    Ok(graph)
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::Graph;
    use petgraph::Undirected;

    /// A tiny undirected graph: 3 nodes with a "qname" attr, 2 weighted edges.
    fn sample() -> Graph<String, f64, Undirected> {
        let mut g: Graph<String, f64, Undirected> = Graph::new_undirected();
        let a = g.add_node("qa".to_string());
        let b = g.add_node("qb".to_string());
        let c = g.add_node("qc".to_string());
        g.add_edge(a, b, 0.5);
        g.add_edge(b, c, 1.25);
        g
    }

    #[test]
    fn write_then_read_round_trip_nodes_edges_attrs() {
        let g = sample();
        let tmp = tempfile::Builder::new().suffix(".graphml").tempfile().unwrap();
        write_graphml(
            tmp.path(),
            &g,
            |q| vec![("qname".to_string(), q.clone())],
            |w| vec![("weight".to_string(), format!("{w}"))],
        )
        .unwrap();

        // structural read-back (Directed container is fine; we count + check attrs)
        let back = read_graphml(tmp.path()).unwrap();
        assert_eq!(back.node_count(), 3, "3 nodes");
        assert_eq!(back.edge_count(), 2, "2 edges");

        // node attrs: every node carries a "qname" in {qa,qb,qc}
        let mut qnames: Vec<String> = back
            .node_weights()
            .map(|m| m.get("qname").cloned().unwrap_or_default())
            .collect();
        qnames.sort();
        assert_eq!(qnames, vec!["qa", "qb", "qc"]);

        // edge attrs: weights round-trip as strings
        let mut weights: Vec<String> = back
            .edge_weights()
            .map(|m| m.get("weight").cloned().unwrap_or_default())
            .collect();
        weights.sort();
        assert_eq!(weights, vec!["0.5", "1.25"]);
    }

    #[test]
    fn read_undirected_round_trip_preserves_topology() {
        let g = sample();
        let tmp = tempfile::Builder::new().suffix(".graphml").tempfile().unwrap();
        write_graphml(
            tmp.path(),
            &g,
            |q| vec![("qname".to_string(), q.clone())],
            |w| vec![("weight".to_string(), format!("{w}"))],
        )
        .unwrap();
        let back = read_graphml_undirected(tmp.path()).unwrap();
        assert_eq!(back.node_count(), 3);
        assert_eq!(back.edge_count(), 2);
        // The middle node (qb) is incident to both edges.
        let qb = back
            .node_indices()
            .find(|&n| back[n].get("qname").map(|s| s == "qb").unwrap_or(false))
            .unwrap();
        assert_eq!(back.edges(qb).count(), 2);
    }

    #[test]
    fn node_without_attrs_round_trips_empty_map() {
        let mut g: Graph<String, f64, Undirected> = Graph::new_undirected();
        g.add_node("x".to_string());
        let tmp = tempfile::Builder::new().suffix(".graphml").tempfile().unwrap();
        // emit no attrs at all
        write_graphml(tmp.path(), &g, |_| Vec::new(), |_| Vec::new()).unwrap();
        let back = read_graphml(tmp.path()).unwrap();
        assert_eq!(back.node_count(), 1);
        assert!(back.node_weights().next().unwrap().is_empty());
    }

    #[test]
    fn read_missing_file_is_io_error() {
        let err = read_graphml(Path::new("/no/such/file.graphml")).unwrap_err();
        assert!(matches!(err, SdError::Io { .. }), "got {err:?}");
    }
}

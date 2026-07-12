//! FrankenNetworkx interop (bead gtql): converters between the canonical
//! [`Csr`] adjacency and fnx's `GraphSnapshot`.
//!
//! fnx-classes' `Graph` is string-keyed and OWNED, so — unlike the borrowed
//! [`Csr::row`] neighbor view — these converters necessarily COPY (fnx's
//! ownership model offers no zero-copy path in that direction; the honest
//! deviation from wsbf's "wrap, do not copy" phrasing). Mapping:
//! node `i` ⇄ its decimal string key `"i"`; a stored entry `(r, c, v)` ⇄ a
//! directed [`EdgeSnapshot`] `r→c` carrying weight `v` under [`WEIGHT_KEY`].
//!
//! Gated behind the `fnx-interop` feature so the default L1 crate stays
//! dependency-lean. Round-trip `Csr → GraphSnapshot → Csr` is the identity
//! (bitwise on values), tested in this module's `tests` (run with
//! `cargo test -p fs-sparse --features fnx-interop`).

use crate::Csr;
use fnx_classes::{AttrMap, EdgeSnapshot, GraphSnapshot};
use fnx_runtime::{CgseValue, CompatibilityMode};
use std::collections::BTreeMap;

/// The edge-attribute key under which scalar weights are carried.
pub const WEIGHT_KEY: &str = "weight";

/// Interop failure — a total, honest error surface.
#[derive(Debug, Clone, PartialEq)]
pub enum InteropError {
    /// A CSR graph adjacency must be square (the node set is rows = cols).
    NotSquare {
        /// Row count of the offending matrix.
        nrows: usize,
        /// Column count of the offending matrix.
        ncols: usize,
    },
    /// An edge referenced a node key absent from `snap.nodes`.
    UnknownNode(String),
}

impl std::fmt::Display for InteropError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InteropError::NotSquare { nrows, ncols } => {
                write!(f, "CSR graph adjacency must be square, got {nrows}x{ncols}")
            }
            InteropError::UnknownNode(k) => {
                write!(f, "edge references unknown node key {k:?}")
            }
        }
    }
}

impl std::error::Error for InteropError {}

/// Read a weight attribute as `f64`: `Float` verbatim, `Int` coerced, anything
/// else (or absent) defaults to `1.0` — the unweighted-graph convention.
fn weight_of(v: &CgseValue) -> f64 {
    match v {
        CgseValue::Float(x) => *x,
        // Interop coercion for integer weights; graph weights are far below the
        // 2^53 exact-integer boundary in every practical case.
        #[allow(clippy::cast_precision_loss)]
        CgseValue::Int(i) => *i as f64,
        _ => 1.0,
    }
}

/// `Csr → GraphSnapshot`. The matrix MUST be square (a graph adjacency); node
/// `i` becomes key `"i"`, and each stored `(r, c, v)` becomes a directed
/// [`EdgeSnapshot`] `r→c` carrying weight `v` under [`WEIGHT_KEY`]. Edge order
/// is row-major, per-row ascending-column (the CSR canonical order), so the
/// output is deterministic.
///
/// # Errors
/// [`InteropError::NotSquare`] when `nrows != ncols` — a non-square CSR is not
/// a graph adjacency and is refused rather than silently reinterpreted.
pub fn csr_to_graph_snapshot(
    a: &Csr,
    mode: CompatibilityMode,
) -> Result<GraphSnapshot, InteropError> {
    if a.nrows() != a.ncols() {
        return Err(InteropError::NotSquare {
            nrows: a.nrows(),
            ncols: a.ncols(),
        });
    }
    let n = a.nrows();
    let nodes: Vec<String> = (0..n).map(|i| i.to_string()).collect();
    let mut edges = Vec::with_capacity(a.nnz());
    for r in 0..n {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            let mut attrs = AttrMap::new();
            attrs.insert(WEIGHT_KEY.to_string(), CgseValue::Float(v));
            edges.push(EdgeSnapshot {
                left: r.to_string(),
                right: c.to_string(),
                attrs,
            });
        }
    }
    Ok(GraphSnapshot {
        mode,
        nodes,
        node_attrs: BTreeMap::new(),
        edges,
    })
}

/// `GraphSnapshot → Csr`. Nodes map to indices BY THEIR ORDER in `snap.nodes`
/// (0-based), so node keys may be arbitrary strings. Each edge's weight is read
/// from [`WEIGHT_KEY`] (`Float`/`Int` coerced; absent or non-numeric ⇒ `1.0`).
/// Parallel edges on the same `(r, c)` SUM their weights (multigraph → simple
/// accumulation, deterministic). The result is square (`N` = node count) and
/// satisfies the canonical CSR invariant.
///
/// # Errors
/// [`InteropError::UnknownNode`] when an edge references a key not in
/// `snap.nodes`.
pub fn graph_snapshot_to_csr(snap: &GraphSnapshot) -> Result<Csr, InteropError> {
    let n = snap.nodes.len();
    let index: BTreeMap<&str, usize> = snap
        .nodes
        .iter()
        .enumerate()
        .map(|(i, k)| (k.as_str(), i))
        .collect();
    // Accumulate into an ordered (r, c) → weight map: BTreeMap iteration is
    // ascending (r, c), which is exactly row-major + per-row ascending-column,
    // so the CSR arrays fall out already canonical (unique, sorted) with
    // parallel edges summed.
    let mut acc: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for e in &snap.edges {
        let r = *index
            .get(e.left.as_str())
            .ok_or_else(|| InteropError::UnknownNode(e.left.clone()))?;
        let c = *index
            .get(e.right.as_str())
            .ok_or_else(|| InteropError::UnknownNode(e.right.clone()))?;
        let w = e.attrs.get(WEIGHT_KEY).map_or(1.0, weight_of);
        // Assign the FIRST contribution verbatim (`or_insert(w)`), accumulate
        // only on a genuine duplicate. Seeding at 0.0 and adding lost the sign
        // of a stored `-0.0` (`0.0 + -0.0 == +0.0`) and canonicalized a NaN
        // payload, breaking the CONTRACT's "bitwise round-trip" guarantee. For
        // finite values `0.0 + w == w` exactly, so summed duplicates are
        // unchanged.
        acc.entry((r, c)).and_modify(|a| *a += w).or_insert(w);
    }
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::with_capacity(acc.len());
    let mut vals = Vec::with_capacity(acc.len());
    for (&(r, c), &v) in &acc {
        row_ptr[r + 1] += 1;
        col_idx.push(c);
        vals.push(v);
    }
    for r in 0..n {
        row_ptr[r + 1] += row_ptr[r];
    }
    Ok(Csr::from_parts(n, n, row_ptr, col_idx, vals))
}

#[cfg(test)]
mod tests {
    use super::{
        InteropError, WEIGHT_KEY, csr_to_graph_snapshot, graph_snapshot_to_csr, weight_of,
    };
    use crate::Csr;
    use fnx_classes::{AttrMap, EdgeSnapshot, GraphSnapshot};
    use fnx_runtime::{CgseValue, CompatibilityMode};
    use std::collections::BTreeMap;

    /// A 4×4 fixture with an empty row (2) and a self-loop (1→1).
    fn sample() -> Csr {
        Csr::from_parts(
            4,
            4,
            vec![0, 2, 3, 3, 5],
            vec![1, 3, 1, 0, 2],
            vec![2.0, -1.5, 3.0, 0.5, 4.0],
        )
    }

    fn edge(left: &str, right: &str, w: Option<f64>) -> EdgeSnapshot {
        let mut attrs = AttrMap::new();
        if let Some(w) = w {
            attrs.insert(WEIGHT_KEY.to_string(), CgseValue::Float(w));
        }
        EdgeSnapshot {
            left: left.to_string(),
            right: right.to_string(),
            attrs,
        }
    }

    fn snapshot(nodes: &[&str], edges: Vec<EdgeSnapshot>) -> GraphSnapshot {
        GraphSnapshot {
            mode: CompatibilityMode::Strict,
            nodes: nodes.iter().map(|s| (*s).to_string()).collect(),
            node_attrs: BTreeMap::new(),
            edges,
        }
    }

    /// Value-bitwise equality of two CSR matrices (guards against ±0/precision).
    fn bits_equal(a: &Csr, b: &Csr) -> bool {
        a.nrows() == b.nrows()
            && a.ncols() == b.ncols()
            && a.nnz() == b.nnz()
            && (0..a.nrows()).all(|r| {
                let (ca, va) = a.row(r);
                let (cb, vb) = b.row(r);
                ca == cb && va.iter().zip(vb).all(|(x, y)| x.to_bits() == y.to_bits())
            })
    }

    #[test]
    fn csr_to_snapshot_shape_and_weights() {
        let a = sample();
        let snap = csr_to_graph_snapshot(&a, CompatibilityMode::Strict).expect("square");
        assert!(
            snap.nodes
                .iter()
                .map(String::as_str)
                .eq(["0", "1", "2", "3"])
        );
        assert_eq!(snap.edges.len(), a.nnz());
        // Row-major, ascending-column: first edge is 0→1 weight 2.0.
        let e0 = &snap.edges[0];
        assert_eq!((e0.left.as_str(), e0.right.as_str()), ("0", "1"));
        assert_eq!(e0.attrs.get(WEIGHT_KEY), Some(&CgseValue::Float(2.0)));
    }

    #[test]
    fn round_trip_is_bitwise_identity() {
        let fixtures = [
            Csr::identity(5),
            sample(),
            Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![7.5]),
            // Regression: a stored -0.0 must round-trip BITWISE (the old
            // `0.0 + w` accumulator canonicalized it to +0.0, silently breaking
            // the CONTRACT's unconditional "bitwise on values" round-trip).
            Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![-0.0]),
        ];
        for a in fixtures {
            let snap = csr_to_graph_snapshot(&a, CompatibilityMode::Strict).expect("square");
            let back = graph_snapshot_to_csr(&snap).expect("valid keys");
            assert!(bits_equal(&a, &back), "round trip changed the matrix");
        }
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact small-integer weights
    fn parallel_edges_sum_weights() {
        let snap = snapshot(
            &["0", "1"],
            vec![edge("0", "1", Some(2.0)), edge("0", "1", Some(3.0))],
        );
        let csr = graph_snapshot_to_csr(&snap).expect("valid");
        assert_eq!(csr.get(0, 1), 5.0, "parallel edges must sum");
        assert_eq!(csr.nnz(), 1, "summed into a single stored entry");
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact default weight
    fn missing_weight_defaults_to_one() {
        let snap = snapshot(&["0", "1"], vec![edge("0", "1", None)]);
        let csr = graph_snapshot_to_csr(&snap).expect("valid");
        assert_eq!(csr.get(0, 1), 1.0);
    }

    #[test]
    fn unknown_node_key_errors() {
        let snap = snapshot(&["0"], vec![edge("0", "ghost", Some(1.0))]);
        assert_eq!(
            graph_snapshot_to_csr(&snap),
            Err(InteropError::UnknownNode("ghost".to_string()))
        );
    }

    #[test]
    fn non_square_csr_refused() {
        let a = Csr::from_parts(2, 3, vec![0, 1, 2], vec![0, 2], vec![1.0, 1.0]);
        assert_eq!(
            csr_to_graph_snapshot(&a, CompatibilityMode::Strict),
            Err(InteropError::NotSquare { nrows: 2, ncols: 3 })
        );
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact weights
    fn arbitrary_node_keys_map_by_order() {
        let snap = snapshot(
            &["alpha", "beta", "gamma"],
            vec![
                edge("alpha", "gamma", Some(9.0)),
                edge("beta", "alpha", Some(1.0)),
            ],
        );
        let csr = graph_snapshot_to_csr(&snap).expect("valid");
        assert_eq!(csr.get(0, 2), 9.0); // alpha→gamma
        assert_eq!(csr.get(1, 0), 1.0); // beta→alpha
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact coercions
    fn weight_coercions() {
        assert_eq!(weight_of(&CgseValue::Float(2.5)), 2.5);
        assert_eq!(weight_of(&CgseValue::Int(3)), 3.0);
        assert_eq!(weight_of(&CgseValue::Bool(true)), 1.0); // non-numeric → 1.0
    }
}

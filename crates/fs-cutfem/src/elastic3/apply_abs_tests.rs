//! `CutElasticity3::apply_abs` against an INDEPENDENT dense recomputation of the absolute-term action and of the rounding chain, written from the private cell and ghost data
//! with different loops from the ones under test, so that a later edit of `apply` that forgets `apply_abs` (or the reverse) fails here.

use super::{CutElasticity3, ElasticityOptions3};
use crate::quad3::{QuadratureControl3, QuadratureOptions3};
use crate::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
use std::ops::ControlFlow;

pub(in crate::elastic3) struct Slab(pub(in crate::elastic3) f64);

impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(self.0, self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}

/// A 2 x 2 x 2 Cartesian operator clamped at x = 0, cut by the plane z = `cut` (a cut above 1 leaves every cell whole and no ghost face), with varied cell scales.
pub(in crate::elastic3) fn cartesian(cut: f64) -> CutElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    let mut op = CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Slab(cut),
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut control,
    )
    .unwrap();
    let varied: Vec<f64> = (0..op.cells.len()).map(|i| 0.4 + 0.07 * (i % 5) as f64).collect();
    op.set_scales(&varied).unwrap();
    op
}

/// A vector of mixed signs and no zero entry.
pub(in crate::elastic3) fn mixed(n: usize) -> Vec<f64> {
    (0..n).map(|i| ((i * 7 + 3) % 11) as f64 / 11.0 - 0.45).collect()
}

/// T[row * n + column] = the sum, over every term that `apply` adds into `row` through `column`, of the absolute coefficient; a fixed row is its own copy.
pub(in crate::elastic3) fn dense_abs_terms(op: &CutElasticity3) -> Vec<f64> {
    let n = op.n();
    let mut t = vec![0.0; n * n];
    for (id, cell) in op.cells.iter().enumerate() {
        for i in 0..24 {
            for j in 0..24 {
                let (row_node, column_node) = (cell.nodes[i / 3], cell.nodes[j / 3]);
                if op.fixed[row_node] || op.fixed[column_node] {
                    continue;
                }
                t[(3 * row_node + i % 3) * n + 3 * column_node + j % 3] += op.scales[id].abs() * cell.stiffness[i][j].abs();
            }
        }
    }
    for face in &op.ghosts {
        let weight = face.scaled_weight(&op.scales).abs();
        for c in 0..3 {
            for (&a, &ja) in face.nodes.iter().zip(&face.jump) {
                for (&b, &jb) in face.nodes.iter().zip(&face.jump) {
                    if !op.fixed[a] && !op.fixed[b] {
                        t[(3 * a + c) * n + 3 * b + c] += weight * ja.abs() * jb.abs();
                    }
                }
            }
        }
    }
    for (node, &fixed) in op.fixed.iter().enumerate() {
        if fixed {
            for c in 0..3 {
                t[(3 * node + c) * n + 3 * node + c] = 1.0;
            }
        }
    }
    t
}

/// The rounding chain recounted from the cell and face lists: 25 for a cell term, `n_f + 2` for a face term, plus the other terms of the row.
pub(in crate::elastic3) fn recount(op: &CutElasticity3) -> usize {
    let mut longest = 0;
    for node in 0..op.nodes.len() {
        if op.fixed[node] {
            continue;
        }
        let cells = op.cells.iter().filter(|cell| cell.nodes.contains(&node)).count();
        let faces: Vec<usize> = op.ghosts.iter().filter(|face| face.nodes.contains(&node)).map(|face| face.nodes.len() + 2).collect();
        let terms = cells + faces.len();
        if terms == 0 {
            continue;
        }
        let chain = faces.iter().copied().chain((cells > 0).then_some(25)).max().unwrap();
        longest = longest.max(chain + terms - 1);
    }
    longest
}

#[test]
fn the_absolute_term_action_equals_an_independent_dense_recomputation_and_the_recounted_chain() {
    for cut in [0.73, 1.1] {
        let op = cartesian(cut);
        if cut < 1.0 {
            assert!(!op.ghosts.is_empty(), "the cut fixture has ghost faces");
        }
        let n = op.n();
        let x = mixed(n);
        let mut y = vec![0.0; n];
        let m = op.apply_abs(&x, &mut y).expect("the cut operator states its chain");
        let t = dense_abs_terms(&op);
        for i in 0..n {
            let reference: f64 = (0..n).map(|j| t[i * n + j] * x[j].abs()).sum();
            assert!((y[i] - reference).abs() <= 1.0e-12 * reference, "cut {cut} row {i}: {} against {reference}", y[i]);
        }
        assert_eq!(m, recount(&op), "cut {cut}");
        assert!(m >= 25, "a cell term alone is a 24-term dot times a scale: {m}");
        let mut ax = vec![0.0; n];
        op.apply(&x, &mut ax);
        let slack = 1.0 + 4.0 * (m as f64 + 1.0) * (f64::EPSILON / 2.0);
        for i in 0..n {
            assert!(ax[i].abs() <= y[i] * slack, "cut {cut} row {i}: |apply| {:e} above apply_abs {:e}", ax[i].abs(), y[i]);
        }
    }
}

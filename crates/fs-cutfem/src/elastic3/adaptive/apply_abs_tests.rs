//! `AdaptiveElasticity3::apply_abs` (hanging nodes, ghost faces) against an INDEPENDENT dense recomputation: |E|^T T |E| |x| from the stored transfer rows and the raw operator's dense absolute terms,
//! and the three stage chains added.

use super::AdaptiveElasticity3;
use crate::elastic3::apply_abs_tests::{Slab, dense_abs_terms, mixed, recount};
use crate::elastic3::ElasticityOptions3;
use crate::octree3::Octree3;
use crate::quad3::{QuadratureControl3, QuadratureOptions3};
use crate::HexCell;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
use std::ops::ControlFlow;

/// A level-1 octree with one leaf refined, so some nodes hang, cut by z = 0.73, clamped at x = 0.
fn hanging() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let mark = *tree.leaves().iter().find(|c| c.index() == [0, 0, 1]).unwrap();
    let tree = tree.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Slab(0.73),
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut control,
    )
    .unwrap()
}

#[test]
fn the_composed_absolute_term_action_equals_a_dense_recomputation_and_adds_the_three_chains() {
    let op = hanging();
    assert!(op.rows.iter().any(|row| row.len() > 1), "the refined fixture has hanging-node rows");
    assert!(!op.raw.ghosts.is_empty(), "the cut fixture has ghost faces");
    let n = op.n();
    let nr = op.raw.n();
    let x = mixed(n);
    // |E|: row 3 node + c, column 3 m + c, the weights of `physical_coefficient` in absolute value, a fixed master contributing nothing.
    let mut e = vec![0.0; nr * n];
    for (node, row) in op.rows.iter().enumerate() {
        for &(m, w) in row {
            if !op.fixed[m] {
                for c in 0..3 {
                    e[(3 * node + c) * n + 3 * m + c] += w.abs();
                }
            }
        }
    }
    let t = dense_abs_terms(&op.raw);
    let ex: Vec<f64> = (0..nr).map(|i| (0..n).map(|j| e[i * n + j] * x[j].abs()).sum()).collect();
    let tex: Vec<f64> = (0..nr).map(|i| (0..nr).map(|j| t[i * nr + j] * ex[j]).sum()).collect();
    let mut reference: Vec<f64> = (0..n).map(|j| (0..nr).map(|i| e[i * n + j] * tex[i]).sum()).collect();
    for (m, &fixed) in op.fixed.iter().enumerate() {
        if fixed {
            for c in 0..3 {
                reference[3 * m + c] = x[3 * m + c].abs();
            }
        }
    }
    let mut y = vec![0.0; n];
    let chain = op.apply_abs(&x, &mut y).expect("the adaptive operator states its chain");
    for i in 0..n {
        assert!((y[i] - reference[i]).abs() <= 1.0e-12 * reference[i], "row {i}: {} against {}", y[i], reference[i]);
    }
    let expand = op.rows.iter().map(Vec::len).max().unwrap();
    let mut references = vec![0_usize; op.fixed.len()];
    for row in &op.rows {
        for &(m, _) in row {
            if !op.fixed[m] {
                references[m] += 1;
            }
        }
    }
    let restrict = references.into_iter().max().unwrap();
    assert_eq!(chain, expand + recount(&op.raw) + restrict);
    let mut ax = vec![0.0; n];
    op.apply(&x, &mut ax);
    let slack = 1.0 + 4.0 * (chain as f64 + 1.0) * (f64::EPSILON / 2.0);
    for i in 0..n {
        assert!(ax[i].abs() <= y[i] * slack, "row {i}: |apply| {:e} above apply_abs {:e}", ax[i].abs(), y[i]);
    }
}

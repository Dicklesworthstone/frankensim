//! `LinearOp::apply_abs` on the matrix-free elasticity operators, against a dense probe of each (n unit-vector applies): |A_probed||x| <= apply_abs(x) and |apply(x)| <= apply_abs(x), entrywise,
//! within gamma of the stated count. Equality is not asserted: wherever two cells add into one entry with opposite signs (mixed-component entries do) the assembled entry is smaller than the sum
//! of its terms. The in-crate tests compare against an independent dense recomputation of the terms themselves.

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;
use std::ops::ControlFlow;

const U: f64 = f64::EPSILON / 2.0;

struct Slab(f64);

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

fn cartesian() -> CutElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Slab(0.73),
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut control,
    )
    .unwrap()
}

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

/// Probe the operator with n unit vectors and compare, returning the stated chain.
fn probed(op: &impl LinearOp, label: &str) -> usize {
    let n = op.n();
    let columns: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let mut column = vec![0.0; n];
            op.apply(&e, &mut column);
            column
        })
        .collect();
    let x: Vec<f64> = (0..n).map(|i| ((i * 7 + 3) % 11) as f64 / 11.0 - 0.45).collect();
    let mut y = vec![0.0; n];
    let m = op.apply_abs(&x, &mut y).unwrap_or_else(|| panic!("{label} states no chain"));
    let mut ax = vec![0.0; n];
    op.apply(&x, &mut ax);
    assert!(m >= 25, "{label}: a cell term alone is a 24-term dot times a scale: {m}");
    let slack = 1.0 + 4.0 * (m as f64 + 1.0) * U;
    for i in 0..n {
        let dense: f64 = (0..n).map(|j| columns[j][i].abs() * x[j].abs()).sum();
        assert!(dense <= y[i] * slack, "{label} row {i}: |A_probed||x| {dense:e} above apply_abs {:e}", y[i]);
        assert!(ax[i].abs() <= y[i] * slack, "{label} row {i}: |apply| {:e} above apply_abs {:e}", ax[i].abs(), y[i]);
    }
    m
}

#[test]
fn cartesian_cut_operator_apply_abs_bounds_the_probed_matrix_and_the_action() {
    probed(&cartesian(), "cartesian");
}

#[test]
fn hanging_node_operator_apply_abs_bounds_the_probed_matrix_and_the_action() {
    probed(&hanging(), "hanging");
}

#[test]
fn the_solve_space_states_what_its_operator_states() {
    let space = AdaptiveSolveSpace3::jacobi(hanging(), 1_000_000);
    let n = space.n();
    let x: Vec<f64> = (0..n).map(|i| ((i * 5 + 1) % 13) as f64 / 13.0 - 0.4).collect();
    let (mut through, mut direct) = (vec![0.0; n], vec![0.0; n]);
    let a = space.apply_abs(&x, &mut through);
    let b = space.elasticity().apply_abs(&x, &mut direct);
    assert!(a.is_some() && a == b);
    assert_eq!(through.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), direct.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
}

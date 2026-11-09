//! Strong support force/moment from the original retained three-dimensional FEM.
use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::equilibrium::EquilibriumResponse3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticityError3, ElasticityFormulation3, ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;

struct Slab(usize, f64);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[self.0] - self.1
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[self.0], hi[self.0]) - Interval::new(self.1, self.1)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let index = match axis {
            HeightAxis::X => 0,
            HeightAxis::Y => 1,
            HeightAxis::Z => 2,
        };
        let value = if index == self.0 { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a:e} != {b:e}");
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn balanced(response: &EquilibriumResponse3, applied: [f64; 3], centroid: [f64; 3]) {
    let moment = cross(
        std::array::from_fn(|i| centroid[i] - response.moment_origin[i]),
        applied,
    );
    for i in 0..3 {
        near(response.support_force[i], -applied[i]);
        near(response.support_moment[i], -moment[i]);
    }
    assert!(
        response.free_residual_norm < 1e-9,
        "{}",
        response.free_residual_norm
    );
}
fn cartesian(cut: f64, nu: f64) -> CutElasticity3 {
    let mut p = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut p,
    )
    .unwrap();
    CutElasticity3::build(
        domain(),
        [2; 3],
        &Slab(2, cut),
        &IsotropicElastic::new(10.0, nu, 0.1).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut q,
    )
    .unwrap()
}

#[test]
fn affine_extension_reaction_and_supported_external_force_are_distinct() {
    let op = cartesian(1.1, 0.0);
    let mut force = vec![0.0; op.n()];
    for (node, p) in op.nodes().iter().enumerate() {
        if p[0] == 1.0 {
            let wy = if p[1] == 0.0 || p[1] == 1.0 { 0.5 } else { 1.0 };
            let wz = if p[2] == 0.0 || p[2] == 1.0 { 0.5 } else { 1.0 };
            force[3 * node] = 0.01 * wy * wz / 4.0;
        }
    }
    let solve = |rhs: &[f64]| {
        op.solve_physical_controlled(rhs, 1e-11, 20_000, 32, |_| poll())
            .unwrap()
    };
    let solution = solve(&force);
    for (p, u) in op
        .nodes()
        .iter()
        .zip(solution.coefficients().chunks_exact(3))
    {
        near(u[0], 0.001 * p[0]);
        near(u[1], 0.0);
        near(u[2], 0.0);
    }
    let origin = [0.2, -0.1, 0.3];
    let report = op
        .equilibrium_response(solution.coefficients(), &force, origin, poll)
        .unwrap();
    balanced(&report, [0.01, 0.0, 0.0], [1.0, 0.5, 0.5]);
    let node = op
        .nodes()
        .iter()
        .position(|p| *p == [0.0, 0.5, 0.5])
        .unwrap();
    let supported_force = [0.004, -0.003, 0.002];
    for i in 0..3 {
        force[3 * node + i] += supported_force[i];
    }
    let replay = solve(&force);
    assert_eq!(replay.coefficients(), solution.coefficients());
    let changed = op
        .equilibrium_response(replay.coefficients(), &force, origin, poll)
        .unwrap();
    let moment = cross(
        std::array::from_fn(|i| op.nodes()[node][i] - origin[i]),
        supported_force,
    );
    for i in 0..3 {
        near(
            changed.support_force[i],
            report.support_force[i] - supported_force[i],
        );
        near(
            changed.support_moment[i],
            report.support_moment[i] - moment[i],
        );
        near(
            changed.residual[3 * node + i],
            report.residual[3 * node + i] - supported_force[i],
        );
    }
}

#[test]
fn full_cut_body_load_and_ghost_forces_close_force_and_moment_balance() {
    let mut op = cartesian(0.73, 0.3);
    op.set_scales(
        &(0..op.cells())
            .map(|i| 0.4 + 0.03 * i as f64)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let density = [0.001, -0.002, -0.003];
    let full = op
        .equilibrium_load(ReferenceLoad3::body(&|_| density), None, poll)
        .unwrap();
    let masked = op.body_load(&|_| density, poll).unwrap();
    for (node, &fixed) in op.fixed().iter().enumerate() {
        for i in 0..3 {
            if fixed {
                assert_eq!(masked[3 * node + i], 0.0);
            } else {
                assert_eq!(full[3 * node + i].to_bits(), masked[3 * node + i].to_bits());
            }
        }
    }
    let total: [f64; 3] = std::array::from_fn(|i| full.iter().skip(i).step_by(3).sum());
    for i in 0..3 {
        near(total[i], 0.73 * density[i]);
    }
    let solution = op
        .solve_physical_controlled(&full, 1e-11, 20_000, 32, |_| poll())
        .unwrap();
    let report = op
        .equilibrium_response(solution.coefficients(), &full, [0.2, 0.1, -0.1], poll)
        .unwrap();
    balanced(&report, density.map(|v| 0.73 * v), [0.5, 0.5, 0.365]);
    let zero = op
        .equilibrium_response(&vec![0.0; op.n()], &full, [0.0; 3], poll)
        .unwrap();
    assert!(
        zero.free_residual_norm > 1e-5,
        "an unsolved field must retain its free imbalance"
    );
}

#[test]
fn adaptive_hanging_transpose_preserves_bbar_support_balance() {
    for graded in [false, true] {
        let tree = Octree3::uniform(1, 3, 128).unwrap();
        let tree = if graded {
            let mark = *tree
                .leaves()
                .iter()
                .find(|c| c.index() == [0, 0, 1])
                .unwrap();
            tree.refined(&[mark], poll).unwrap()
        } else {
            tree
        };
        let mut p = |_| poll();
        let mut q = QuadratureControl3::new(
            QuadratureOptions3 {
                depth: 1,
                ..Default::default()
            },
            &mut p,
        )
        .unwrap();
        let mut op = AdaptiveElasticity3::build(
            domain(),
            &tree,
            &Slab(2, 0.73),
            &IsotropicElastic::new(10.0, 0.49, 0.1).unwrap(),
            &|p| p[0] == 0.0,
            ElasticityOptions3 {
                formulation: ElasticityFormulation3::BBar,
                ..Default::default()
            },
            &mut q,
        )
        .unwrap();
        if graded {
            assert!(op.physical_nodes().len() > op.nodes().len());
        }
        op.set_scales(
            &(0..op.cells())
                .map(|i| 0.4 + 0.02 * i as f64)
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let density = [0.001, -0.002, -0.003];
        let full = op
            .equilibrium_load(ReferenceLoad3::body(&|_| density), None, poll)
            .unwrap();
        let masked = op.body_load(&|_| density, poll).unwrap();
        for (i, (&a, &b)) in full.iter().zip(&masked).enumerate() {
            if !op.fixed()[i / 3] {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
        let solution = op
            .solve_physical_controlled(&full, 1e-11, 40_000, 32, |_| poll())
            .unwrap();
        let report = op
            .equilibrium_response(solution.coefficients(), &full, [0.2, 0.1, -0.1], poll)
            .unwrap();
        balanced(&report, density.map(|v| 0.73 * v), [0.5, 0.5, 0.365]);
    }
}

#[test]
fn embedded_motion_lifting_is_retained_at_strong_supports() {
    let mut p = |_| poll();
    let mut q = QuadratureControl3::new(Default::default(), &mut p).unwrap();
    // One x cell makes the embedded lifting act directly on fixed x=0 rows.
    let mut op = CutElasticity3::build_with_embedded_dirichlet(
        domain(),
        [1, 2, 2],
        &Slab(0, 0.83),
        &IsotropicElastic::new(10.0, 0.0, 0.1).unwrap(),
        &|p| p[0] == 0.0,
        &|_, _| true,
        ElasticityOptions3::default(),
        Default::default(),
        Default::default(),
        &mut q,
    )
    .unwrap();
    op.set_scales(&vec![0.6; op.cells()]).unwrap();
    let motion = |p: [f64; 3], _: [f64; 3]| [0.001 * p[0], 0.0, 0.0];
    let full = op
        .equilibrium_load(ReferenceLoad3::default(), Some(&motion), poll)
        .unwrap();
    let masked = op.prescribed_displacement_load(&motion, poll).unwrap();
    assert!(
        full.iter()
            .enumerate()
            .any(|(i, v)| op.fixed()[i / 3] && v.abs() > 1e-6)
    );
    for (i, (&a, &b)) in full.iter().zip(&masked).enumerate() {
        if !op.fixed()[i / 3] {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    let solution = op
        .solve_physical_controlled(&full, 1e-11, 20_000, 32, |_| poll())
        .unwrap();
    let report = op
        .equilibrium_response(solution.coefficients(), &full, [0.0; 3], poll)
        .unwrap();
    near(report.support_force[0], -0.006);
    let embedded = op
        .embedded_reaction(
            solution.coefficients(),
            Some(&motion),
            &|_, _| [1.0, 0.0, 0.0],
            poll,
        )
        .unwrap();
    near(report.support_force[0] + embedded.value, 0.0);
    assert!(report.free_residual_norm < 1e-9);
}

#[test]
fn pressure_at_fully_fixed_nodes_still_has_a_physical_reaction() {
    let mut p = |_| poll();
    let mut q = QuadratureControl3::new(Default::default(), &mut p).unwrap();
    let op = CutElasticity3::build_with_surface(
        domain(),
        [1; 3],
        &Slab(2, 0.73),
        &IsotropicElastic::new(10.0, 0.3, 0.1).unwrap(),
        &|_| true,
        ElasticityOptions3::default(),
        Default::default(),
        &mut q,
    )
    .unwrap();
    assert!(op.fixed().iter().all(|v| *v));
    let full = op
        .equilibrium_load(ReferenceLoad3::pressure(&|_| 0.003), None, poll)
        .unwrap();
    assert!(
        op.pressure_load(&|_| 0.003, poll)
            .unwrap()
            .rhs
            .iter()
            .all(|v| *v == 0.0)
    );
    let report = op
        .equilibrium_response(&vec![0.0; op.n()], &full, [0.2, 0.1, -0.1], poll)
        .unwrap();
    balanced(&report, [0.0, 0.0, -0.003], [0.5, 0.5, 0.73]);
}

#[test]
fn invalid_or_cancelled_response_publishes_nothing_and_preserves_scales() {
    let op = cartesian(0.73, 0.3);
    let u = vec![0.0; op.n()];
    let rhs = op
        .equilibrium_load(ReferenceLoad3::body(&|_| [0.0, 0.0, -0.001]), None, poll)
        .unwrap();
    let scales = op.scales().to_vec();
    let mut bad = rhs.clone();
    bad[0] = f64::NAN;
    assert!(op.equilibrium_response(&u, &bad, [0.0; 3], poll).is_err());
    assert!(
        op.equilibrium_response(&u, &rhs[..1], [0.0; 3], poll)
            .is_err()
    );
    assert!(
        op.equilibrium_response(&u, &rhs, [f64::INFINITY; 3], poll)
            .is_err()
    );
    let outside: Vec<_> = op.nodes().iter().flat_map(|p| [p[0], 0.0, 0.0]).collect();
    assert!(
        matches!(op.equilibrium_response(&outside, &rhs, [0.0; 3], poll),
        Err(ElasticityError3::Invalid(message)) if message.contains("strain"))
    );
    let mut calls = 0;
    let reference = op
        .equilibrium_response(&u, &rhs, [0.0; 3], || {
            calls += 1;
            poll()
        })
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        op.equilibrium_response(&u, &rhs, [0.0; 3], || {
            current += 1;
            if current == calls {
                ControlFlow::Break(())
            } else {
                poll()
            }
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(matches!(
        op.equilibrium_load(ReferenceLoad3::body(&|_| [0.0; 3]), None, || {
            ControlFlow::Break(())
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(op.scales(), scales);
    assert_eq!(
        op.equilibrium_response(&u, &rhs, [0.0; 3], poll).unwrap(),
        reference
    );
}

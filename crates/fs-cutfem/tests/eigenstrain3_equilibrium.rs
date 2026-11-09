//! Complete thermal force balance with strong and embedded displacement supports.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::equilibrium::EquilibriumResponse3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::op::LinearOp;

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
fn control<T>(f: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    f(&mut q)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a:e} != {b:e}");
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn balanced(report: &EquilibriumResponse3, applied: [f64; 3], centroid: [f64; 3]) {
    let moment = cross(
        std::array::from_fn(|i| centroid[i] - report.moment_origin[i]),
        applied,
    );
    for i in 0..3 {
        near(report.support_force[i], -applied[i]);
        near(report.support_moment[i], -moment[i]);
    }
    assert!(
        report.free_residual_norm < 1e-9,
        "{}",
        report.free_residual_norm
    );
}
fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation,
        ..Default::default()
    }
}
fn cartesian(formulation: ElasticityFormulation3, limit: f64) -> CutElasticity3 {
    control(|q| {
        CutElasticity3::build(
            domain(),
            [2; 3],
            &Slab(2, 0.73),
            &IsotropicElastic::new(10.0, 0.49, limit).unwrap(),
            &|p| p[0] == 0.0,
            options(formulation),
            q,
        )
        .unwrap()
    })
}
fn adaptive(formulation: ElasticityFormulation3, limit: f64) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 3, 256).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    control(|q| {
        AdaptiveElasticity3::build(
            domain(),
            &tree,
            &Slab(2, 0.73),
            &IsotropicElastic::new(10.0, 0.49, limit).unwrap(),
            &|p| p[0] == 0.0,
            options(formulation),
            q,
        )
        .unwrap()
    })
}

#[test]
fn compatible_expansion_has_no_strong_reaction_beyond_narrow_elastic_calibration() {
    let expansion = 0.003;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = cartesian(formulation, 1e-5);
        op.set_scales(
            &(0..op.cells())
                .map(|i| 0.4 + 0.02 * i as f64)
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let samples = Cell::new(0);
        let prepared = op
            .prepare_eigenstrain(
                &|_, _| {
                    samples.set(samples.get() + 1);
                    Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0])
                },
                100_000,
                poll,
            )
            .unwrap();
        let sample_count = samples.get();
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [expansion * p[0], 0.0, 0.0])
            .collect();
        let full = op
            .eigenstrain_equilibrium_load(&prepared, ReferenceLoad3::default(), None, poll)
            .unwrap();
        let masked = op.eigenstrain_load(&prepared, poll).unwrap();
        assert!(
            full.iter()
                .enumerate()
                .any(|(i, f)| op.fixed()[i / 3] && f.abs() > 1e-4)
        );
        for (i, (&a, &b)) in full.iter().zip(&masked).enumerate() {
            if op.fixed()[i / 3] {
                assert_eq!(b, 0.0);
            } else {
                near(a, b);
            }
        }
        assert!(op.equilibrium_response(&u, &full, [0.0; 3], poll).is_err());
        let mut calls = 0;
        let report = op
            .eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                100_000,
                [0.2, -0.1, 0.3],
                || {
                    calls += 1;
                    poll()
                },
            )
            .unwrap();
        for &force in &report.residual {
            near(force, 0.0);
        }
        balanced(&report, [0.0; 3], [0.0; 3]);
        assert_eq!(samples.get(), sample_count);
        assert!(
            op.eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                prepared.bulk_points() - 1,
                [0.0; 3],
                poll,
            )
            .is_err()
        );
        assert!(
            op.eigenstrain_equilibrium_response(
                &prepared, &u, &full, 0.0001, 100_000, [0.0; 3], poll,
            )
            .is_err()
        );
        assert!(
            op.eigenstrain_equilibrium_response(
                &prepared,
                &vec![0.0; op.n()],
                &full,
                0.02,
                100_000,
                [0.0; 3],
                poll,
            )
            .is_err()
        );
        let other = cartesian(formulation, 1e-5);
        let external_calls = Cell::new(0);
        assert!(
            other
                .eigenstrain_equilibrium_load(
                    &prepared,
                    ReferenceLoad3::body(&|_| {
                        external_calls.set(external_calls.get() + 1);
                        [0.0; 3]
                    }),
                    None,
                    poll,
                )
                .is_err()
        );
        assert_eq!(external_calls.get(), 0);
        let mut current = 0;
        assert!(matches!(
            op.eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                100_000,
                [0.2, -0.1, 0.3],
                || {
                    current += 1;
                    if current == calls {
                        ControlFlow::Break(())
                    } else {
                        poll()
                    }
                },
            ),
            Err(ElasticityError3::Cancelled)
        ));
        assert_eq!(
            op.eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                100_000,
                [0.2, -0.1, 0.3],
                poll,
            )
            .unwrap(),
            report
        );

        let mut op = adaptive(formulation, 1e-5);
        assert!(op.physical_nodes().len() > op.nodes().len());
        op.set_scales(
            &(0..op.cells())
                .map(|i| 0.4 + 0.02 * i as f64)
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let prepared = op
            .prepare_eigenstrain(
                &|_, _| Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0]),
                100_000,
                poll,
            )
            .unwrap();
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [expansion * p[0], 0.0, 0.0])
            .collect();
        let full = op
            .eigenstrain_equilibrium_load(&prepared, ReferenceLoad3::default(), None, poll)
            .unwrap();
        assert!(op.equilibrium_response(&u, &full, [0.0; 3], poll).is_err());
        let report = op
            .eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                100_000,
                [0.2, -0.1, 0.3],
                poll,
            )
            .unwrap();
        for &force in &report.residual {
            near(force, 0.0);
        }
        balanced(&report, [0.0; 3], [0.0; 3]);
    }
}

#[test]
fn restrained_source_material_preserves_supported_body_forces_and_thermal_face_forces() {
    let young = 12.0;
    let density = 2700.0;
    let identity = ContentHash([91; 32]);
    let stiffness =
        std::array::from_fn(|i| std::array::from_fn(|j| if i == j { young } else { 0.0 }));
    let material =
        ElasticAcousticState::try_new_mandel(stiffness, density, identity, 0.05).unwrap();
    let mut op = control(|q| {
        CutElasticity3::build_with_material_field(
            domain(),
            [1; 3],
            &Slab(2, 0.73),
            ElasticMaterialField3::Uniform(&material),
            &|_| true,
            options(ElasticityFormulation3::Displacement),
            q,
        )
        .unwrap()
    });
    let expansion = 0.002;
    let prepared = op
        .prepare_eigenstrain(
            &|_, _| Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0]),
            100_000,
            poll,
        )
        .unwrap();
    let acceleration: [f64; 3] = [2e-6, -1e-6, 3e-6];
    let body = acceleration.map(|a| density * a);
    let mass = op.reference_cell_masses_kg().unwrap()[0];
    assert!((mass - density * 0.73).abs() < 1e-6);
    let u = vec![0.0; op.n()];
    for scale in [0.4, 0.75] {
        op.set_scales(&[scale]).unwrap();
        let full = op
            .eigenstrain_equilibrium_load(&prepared, ReferenceLoad3::body(&|_| body), None, poll)
            .unwrap();
        assert!(
            op.eigenstrain_load(&prepared, poll)
                .unwrap()
                .iter()
                .all(|f| *f == 0.0)
        );
        let report = op
            .eigenstrain_equilibrium_response(
                &prepared,
                &u,
                &full,
                0.02,
                100_000,
                [0.2, -0.1, 0.3],
                poll,
            )
            .unwrap();
        balanced(&report, body.map(|f| 0.73 * f), [0.5, 0.5, 0.365]);
        assert_eq!(report.free_residual_norm, 0.0);
        for (coordinate, sign) in [(0.0, 1.0), (1.0, -1.0)] {
            let force: f64 = op
                .nodes()
                .iter()
                .enumerate()
                .filter(|(_, p)| p[0] == coordinate)
                .map(|(node, _)| report.residual[3 * node])
                .sum();
            near(
                force,
                sign * scale * young * expansion * 0.73 - 0.5 * body[0] * 0.73,
            );
        }
        assert_eq!(op.reference_cell_masses_kg().unwrap()[0], mass);
        assert_eq!(
            op.cell_material_state(0).unwrap().material_state_identity(),
            identity
        );
    }
}

#[test]
fn graded_bbar_thermal_equilibrium_balances_external_force_and_origin_dependent_moment() {
    let mut op = adaptive(ElasticityFormulation3::BBar, 0.1);
    op.set_scales(
        &(0..op.cells())
            .map(|i| 0.4 + 0.02 * i as f64)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let prepared = op
        .prepare_eigenstrain(
            &|_, p| {
                Ok([
                    0.0004 * (1.0 + p[0]),
                    -0.0002 * p[1],
                    0.0003 * p[2],
                    0.00005 * p[1],
                    -0.00003 * p[0],
                    0.00004 * p[2],
                ])
            },
            100_000,
            poll,
        )
        .unwrap();
    let body = [0.0001, -0.0002, -0.0003];
    let full = op
        .eigenstrain_equilibrium_load(&prepared, ReferenceLoad3::body(&|_| body), None, poll)
        .unwrap();
    let masked_thermal = op.eigenstrain_load(&prepared, poll).unwrap();
    let masked_external = op.body_load(&|_| body, poll).unwrap();
    for (i, &force) in full.iter().enumerate() {
        if !op.fixed()[i / 3] {
            near(force, masked_thermal[i] + masked_external[i]);
        }
    }
    let solution = op
        .solve_controlled(&full, 1e-11, 40_000, 32, |_| poll())
        .unwrap();
    let u = solution.coefficients();
    for origin in [[0.2, 0.1, -0.1], [-0.3, 0.4, 0.2]] {
        let report = op
            .eigenstrain_equilibrium_response(&prepared, u, &full, 0.05, 100_000, origin, poll)
            .unwrap();
        balanced(&report, body.map(|f| 0.73 * f), [0.5, 0.5, 0.365]);
        let node = op.fixed().iter().position(|fixed| *fixed).unwrap();
        let extra = [0.0007, -0.0005, 0.0002];
        let mut with_nodal_load = full.clone();
        for i in 0..3 {
            with_nodal_load[3 * node + i] += extra[i];
        }
        let changed = op
            .eigenstrain_equilibrium_response(
                &prepared,
                u,
                &with_nodal_load,
                0.05,
                100_000,
                origin,
                poll,
            )
            .unwrap();
        let extra_moment = cross(
            std::array::from_fn(|i| op.nodes()[node][i] - origin[i]),
            extra,
        );
        for i in 0..3 {
            near(changed.support_force[i], report.support_force[i] - extra[i]);
            near(
                changed.support_moment[i],
                report.support_moment[i] - extra_moment[i],
            );
        }
        near(changed.free_residual_norm, report.free_residual_norm);
    }
}

#[test]
fn strong_and_embedded_thermal_reactions_balance_without_counting_either_twice() {
    let expansion = 0.002;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = control(|q| {
            CutElasticity3::build_with_embedded_dirichlet(
                domain(),
                [1, 2, 2],
                &Slab(0, 0.83),
                &IsotropicElastic::new(10.0, 0.0, 0.1).unwrap(),
                &|p| p[0] == 0.0,
                &|_, _| true,
                options(formulation),
                Default::default(),
                Default::default(),
                q,
            )
            .unwrap()
        });
        op.set_scales(&vec![0.6; op.cells()]).unwrap();
        let prepared = op
            .prepare_eigenstrain(
                &|_, _| Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0]),
                100_000,
                poll,
            )
            .unwrap();
        let origin = [0.2, -0.1, 0.3];
        for imposed in [0.0, 0.5 * expansion, expansion] {
            let motion = |p: [f64; 3], _: [f64; 3]| [imposed * p[0], 0.0, 0.0];
            let u: Vec<_> = op
                .nodes()
                .iter()
                .flat_map(|p| [imposed * p[0], 0.0, 0.0])
                .collect();
            let full = op
                .eigenstrain_equilibrium_load(
                    &prepared,
                    ReferenceLoad3::default(),
                    Some(&motion),
                    poll,
                )
                .unwrap();
            let report = op
                .eigenstrain_equilibrium_response(&prepared, &u, &full, 0.02, 100_000, origin, poll)
                .unwrap();
            let force = 0.6 * 10.0 * (expansion - imposed);
            near(report.support_force[0], force);
            near(report.support_force[1], 0.0);
            near(report.support_force[2], 0.0);
            assert!(report.free_residual_norm < 1e-9);
            let expected_moment = cross(
                [0.0 - origin[0], 0.5 - origin[1], 0.5 - origin[2]],
                [force, 0.0, 0.0],
            );
            for i in 0..3 {
                near(report.support_moment[i], expected_moment[i]);
            }
            let embedded = op
                .eigenstrain_embedded_reaction(
                    &prepared,
                    &u,
                    Some(&motion),
                    &|_, _| [1.0, 0.0, 0.0],
                    poll,
                )
                .unwrap();
            near(report.support_force[0] + embedded.value, 0.0);
            let embedded_moment = op
                .eigenstrain_embedded_reaction(
                    &prepared,
                    &u,
                    Some(&motion),
                    &|p, _| [-(p[1] - origin[1]), p[0] - origin[0], 0.0],
                    poll,
                )
                .unwrap();
            near(report.support_moment[2] + embedded_moment.value, 0.0);
        }
    }
}

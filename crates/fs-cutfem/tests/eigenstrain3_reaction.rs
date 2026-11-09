//! G1/G3/G4: thermal support force, moment and exact design derivatives.

use std::cell::Cell;
use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticityError3, ElasticityFormulation3, ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::LinearOp;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        (p[0] - 0.17) * (p[0] - 0.83)
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X {
            Interval::new(2.0 * lo[0] - 1.0, 2.0 * hi[0] - 1.0)
        } else {
            Interval::new(0.0, 0.0)
        }
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}
fn right(_: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    if normal[0] > 0.0 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0; 3]
    }
}
fn right_moment_z(p: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    if normal[0] > 0.0 {
        [-p[1], p[0], 0.0]
    } else {
        [0.0; 3]
    }
}
fn cartesian(formulation: ElasticityFormulation3, poisson: f64) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    CutElasticity3::build_with_embedded_dirichlet(
        domain(),
        [2; 3],
        &Slab,
        &IsotropicElastic::new(30.0, poisson, 0.1).unwrap(),
        &|_| false,
        &|_, _| true,
        ElasticityOptions3 {
            formulation,
            ..Default::default()
        },
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap()
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let mut callback = |_| poll();
    let mut control =
        QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(
        domain(),
        &tree,
        &Slab,
        &IsotropicElastic::new(30.0, 0.2, 0.1).unwrap(),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap()
}

#[test]
fn affine_free_expansion_has_no_thermal_support_force_or_moment() {
    let amplitude = 0.001;
    let motion = |p: [f64; 3], _: [f64; 3]| p.map(|x| amplitude * x);
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = cartesian(formulation, 0.25);
        op.set_scales(&vec![0.6; op.cells()]).unwrap();
        let calls = Cell::new(0);
        let prepared = op
            .prepare_eigenstrain(
                &|_, _| {
                    calls.set(calls.get() + 1);
                    Ok([amplitude, amplitude, amplitude, 0.0, 0.0, 0.0])
                },
                100_000,
                poll,
            )
            .unwrap();
        let samples = calls.get();
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| p.map(|x| amplitude * x))
            .collect();
        let uncorrected = op
            .embedded_reaction(&u, Some(&motion), &right, poll)
            .unwrap();
        assert!(uncorrected.value.abs() > 0.01);
        let reaction = op
            .eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right, poll)
            .unwrap();
        close(reaction.value, 0.0, 0.0, 2e-12);
        close(
            op.eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right_moment_z, poll)
                .unwrap()
                .value,
            0.0,
            0.0,
            2e-12,
        );
        assert_eq!(
            reaction.displacement_gradient,
            uncorrected.displacement_gradient
        );
        assert_eq!(
            calls.get(),
            samples,
            "reaction must consume retained free-strain stresses"
        );
    }
    let op = adaptive();
    assert!(op.physical_nodes().len() > op.nodes().len());
    let prepared = op
        .prepare_eigenstrain(
            &|_, _| Ok([amplitude, amplitude, amplitude, 0.0, 0.0, 0.0]),
            100_000,
            poll,
        )
        .unwrap();
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|p| p.map(|x| amplitude * x))
        .collect();
    close(
        op.eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right, poll)
            .unwrap()
            .value,
        0.0,
        0.0,
        3e-12,
    );
    close(
        op.eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right_moment_z, poll)
            .unwrap()
            .value,
        0.0,
        0.0,
        3e-12,
    );
}

#[test]
fn restrained_axial_expansion_has_analytic_compressive_force_and_moment() {
    let expansion = 0.002;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = cartesian(formulation, 0.0);
        op.set_scales(&vec![0.6; op.cells()]).unwrap();
        let prepared = op
            .prepare_eigenstrain(
                &|_, _| Ok([expansion, 0.0, 0.0, 0.0, 0.0, 0.0]),
                100_000,
                poll,
            )
            .unwrap();
        // At nu=0 this confined axial strain has no lateral traction, so u=0
        // satisfies equilibrium with the two original embedded planes held.
        assert!(
            op.eigenstrain_load(&prepared, poll)
                .unwrap()
                .iter()
                .all(|f| f.abs() < 2e-12)
        );
        let u = vec![0.0; op.n()];
        let force = -0.6 * 30.0 * expansion; // unit face area
        close(
            op.eigenstrain_embedded_reaction(&prepared, &u, None, &right, poll)
                .unwrap()
                .value,
            force,
            2e-12,
            2e-13,
        );
        let left = |_: [f64; 3], normal: [f64; 3]| {
            if normal[0] < 0.0 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0; 3]
            }
        };
        close(
            op.eigenstrain_embedded_reaction(&prepared, &u, None, &left, poll)
                .unwrap()
                .value,
            -force,
            2e-12,
            2e-13,
        );
        close(
            op.eigenstrain_embedded_reaction(&prepared, &u, None, &right_moment_z, poll)
                .unwrap()
                .value,
            -0.5 * force,
            2e-12,
            2e-13,
        );
        let about_offset = |p: [f64; 3], normal: [f64; 3]| {
            if normal[0] > 0.0 {
                [-(p[1] - 0.2), p[0], 0.0]
            } else {
                [0.0; 3]
            }
        };
        close(
            op.eigenstrain_embedded_reaction(&prepared, &u, None, &about_offset, poll)
                .unwrap()
                .value,
            -0.3 * force,
            2e-12,
            2e-13,
        );
    }
}

#[test]
fn adaptive_reaction_partial_and_equilibrium_gradients_include_original_thermal_load() {
    let mut op = adaptive();
    let prepared = op
        .prepare_eigenstrain(
            &|_, p| {
                Ok([
                    0.001 * (1.0 + p[0]),
                    -0.0002 * p[2],
                    0.0003 * p[1],
                    0.0001 * p[1],
                    0.0002 * p[0],
                    -0.0001 * p[2],
                ])
            },
            100_000,
            poll,
        )
        .unwrap();
    let base: Vec<_> = (0..op.cells()).map(|i| 0.5 + 0.01 * i as f64).collect();
    op.set_scales(&base).unwrap();
    let motion = |p: [f64; 3], _: [f64; 3]| [0.0003 * p[1], 0.0002 * p[2], 0.0001];
    let external = op.body_load(&|_| [0.0001, 0.0002, 0.0], poll).unwrap();
    let solve = |op: &AdaptiveElasticity3, rhs: &[f64]| {
        op.solve_controlled(rhs, 1e-11, 30_000, 32, |_| poll())
            .unwrap()
            .coefficients()
            .to_vec()
    };
    let forward = |op: &AdaptiveElasticity3| {
        let mut rhs = op.eigenstrain_load(&prepared, poll).unwrap();
        let lifting = op.prescribed_displacement_load(&motion, poll).unwrap();
        for ((rhs, external), lifting) in rhs.iter_mut().zip(&external).zip(lifting) {
            *rhs += external + lifting;
        }
        solve(op, &rhs)
    };
    let u = forward(&op);
    let reaction = op
        .eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right_moment_z, poll)
        .unwrap();
    let original = op
        .embedded_reaction(&u, Some(&motion), &right_moment_z, poll)
        .unwrap();
    assert_eq!(
        reaction.displacement_gradient,
        original.displacement_gradient
    );
    let direction: Vec<_> = (0..op.n())
        .map(|i| 0.0001 * ((i % 5) as f64 - 2.0))
        .collect();
    let h = 1e-4;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    let observed = |op: &AdaptiveElasticity3, field: &[f64]| {
        op.eigenstrain_embedded_reaction(&prepared, field, Some(&motion), &right_moment_z, poll)
            .unwrap()
            .value
    };
    close(
        (observed(&op, &plus) - observed(&op, &minus)) / (2.0 * h),
        dot(&reaction.displacement_gradient, &direction),
        2e-7,
        1e-10,
    );

    let z = solve(&op, &reaction.displacement_gradient);
    let stiffness = op.scale_bilinear_forms(&z, &u, poll).unwrap();
    let thermal = op
        .eigenstrain_load_scale_pullback(&prepared, &z, poll)
        .unwrap();
    let lifting = op
        .prescribed_displacement_scale_work(&motion, &z, poll)
        .unwrap();
    let total: Vec<_> = (0..op.cells())
        .map(|i| reaction.scale_gradient[i] + thermal[i] + lifting[i] - stiffness[i])
        .collect();
    let ds: Vec<_> = (0..op.cells()).map(|i| 0.1 * (1 + i % 3) as f64).collect();
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s + h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let partial_plus = observed(&op, &u);
    let total_plus = observed(&op, &forward(&op));
    op.set_scales(
        &base
            .iter()
            .zip(&ds)
            .map(|(s, d)| s - h * d)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let partial_minus = observed(&op, &u);
    let total_minus = observed(&op, &forward(&op));
    op.set_scales(&base).unwrap();
    close(
        (partial_plus - partial_minus) / (2.0 * h),
        dot(&reaction.scale_gradient, &ds),
        2e-7,
        1e-10,
    );
    close(
        (total_plus - total_minus) / (2.0 * h),
        dot(&total, &ds),
        5e-5,
        1e-9,
    );

    let other = adaptive();
    assert!(
        other
            .eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right, poll)
            .is_err()
    );
    assert!(
        op.eigenstrain_embedded_reaction(
            &prepared,
            &u,
            Some(&motion),
            &|_, _| [f64::NAN, 0.0, 0.0],
            poll,
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        op.eigenstrain_embedded_reaction(&prepared, &u, Some(&motion), &right_moment_z, || {
            polls += 1;
            if polls > 20 {
                ControlFlow::Break(())
            } else {
                poll()
            }
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert_eq!(polls, 21);
    assert_eq!(op.scales(), base.as_slice());
    assert_eq!(observed(&op, &u).to_bits(), reaction.value.to_bits());
}

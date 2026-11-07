//! Physical-volume B-bar energy, locking relief, and retained material/design data.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
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

struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(self.0, self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let derivative = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(derivative, derivative)
    }
}

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn options(ghost_gamma: f64) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation: ElasticityFormulation3::BBar,
        ghost_gamma,
        ..Default::default()
    }
}
fn with_quadrature<T>(build: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    build(&mut control)
}
fn legacy(
    counts: [usize; 3],
    cut: f64,
    law: &IsotropicElastic,
    options: ElasticityOptions3,
) -> CutElasticity3 {
    with_quadrature(|control| {
        CutElasticity3::build(
            domain(),
            counts,
            &Slab(cut),
            law,
            &|p| p[0] == 0.0,
            options,
            control,
        )
        .unwrap()
    })
}
fn native(
    counts: [usize; 3],
    cut: f64,
    state: &ElasticAcousticState,
    options: ElasticityOptions3,
) -> Result<CutElasticity3, ElasticityError3> {
    with_quadrature(|control| {
        CutElasticity3::build_with_material_field(
            domain(),
            counts,
            &Slab(cut),
            ElasticMaterialField3::Uniform(state),
            &|p| p[0] == 0.0,
            options,
            control,
        )
    })
}
fn adaptive(
    tree: &Octree3,
    law: &IsotropicElastic,
    options: ElasticityOptions3,
) -> AdaptiveElasticity3 {
    with_quadrature(|control| {
        AdaptiveElasticity3::build(
            domain(),
            tree,
            &Slab(0.73),
            law,
            &|p| p[0] == 0.0,
            options,
            control,
        )
        .unwrap()
    })
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()),
        "{a:e} != {b:e} (relative tolerance {relative:e})"
    );
}
fn quadratic(op: &impl LinearOp, u: &[f64]) -> f64 {
    let mut ku = vec![0.0; op.n()];
    op.apply(u, &mut ku);
    u.iter().zip(ku).map(|(u, ku)| u * ku).sum()
}
fn probe(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.01 * x * x + 0.002 * x * y,
                -0.007 * x * z,
                0.004 * x * y * z,
            ]
        })
        .collect()
}
fn same_action(a: &impl LinearOp, b: &impl LinearOp, u: &[f64]) {
    assert_eq!(a.n(), b.n());
    let mut au = vec![0.0; a.n()];
    let mut bu = vec![0.0; b.n()];
    a.apply(u, &mut au);
    b.apply(u, &mut bu);
    let scale = au.iter().chain(&bu).map(|v| v.abs()).fold(0.0, f64::max);
    let error = au
        .iter()
        .zip(&bu)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= 2e-11 * scale,
        "action error {error:e}, scale {scale:e}"
    );
}

#[test]
fn projected_dilatation_uses_the_physical_cut_volume() {
    let amplitude: f64 = 0.02;
    for poisson in [0.3, 0.4999] {
        let law = IsotropicElastic::new(1200.0, poisson, 0.1).unwrap();
        let mu = law.youngs / (2.0 * (1.0 + poisson));
        let bulk = law.youngs / (3.0 * (1.0 - 2.0 * poisson));
        for cut in [0.2_f64, 0.73] {
            let op = legacy([1; 3], cut, &law, options(0.0));
            let u: Vec<_> = op
                .nodes()
                .iter()
                .flat_map(|p| [amplitude * p[0] * p[2], 0.0, 0.0])
                .collect();
            // u_x = a*x*z: div(u)/a = z, whose PHYSICAL mean is cut/2.
            // Integrate fully resolved deviatoric/shear work plus K*V*mean².
            let expected = amplitude.powi(2)
                * (4.0 * mu * cut.powi(3) / 9.0 + mu * cut / 3.0 + bulk * cut.powi(3) / 4.0);
            close(quadratic(&op, &u), expected, 2e-11);
            close(op.volumes().iter().sum(), cut, 2e-12);
        }
    }
}

#[test]
fn nearly_incompressible_pure_bending_energy_converges_without_volumetric_locking() {
    let law = IsotropicElastic::new(1200.0, 0.4999, 0.1).unwrap();
    let state = ElasticAcousticState::from_isotropic(&law, 1.0, ContentHash([51; 32])).unwrap();
    let curvature: f64 = 0.02;
    // Pure bending has sigma_xx = -E*kappa*z and all other stresses zero.
    // For this unit cube, u^T K u = E*kappa²*integral(z²) = E*kappa²/12.
    let exact = law.youngs * curvature.powi(2) / 12.0;
    let mut errors = [[0.0; 3]; 2];
    for (level, n) in [2, 4, 8].into_iter().enumerate() {
        for (kind, formulation) in [
            ElasticityFormulation3::Displacement,
            ElasticityFormulation3::BBar,
        ]
        .into_iter()
        .enumerate()
        {
            let op = with_quadrature(|control| {
                CutElasticity3::build_with_material_field(
                    HexCell::try_new([0.0, -0.5, -0.5], [1.0, 0.5, 0.5]).unwrap(),
                    [n; 3],
                    &Slab(1.0),
                    ElasticMaterialField3::Uniform(&state),
                    &|p| p == [0.0; 3],
                    ElasticityOptions3 {
                        formulation,
                        ..options(0.0)
                    },
                    control,
                )
                .unwrap()
            });
            // Exact compatible bending displacement, interpolated by Q1.
            // The sole clamped node is zero in this field; no solve is claimed.
            let u: Vec<_> = op
                .nodes()
                .iter()
                .flat_map(|&[x, y, z]| {
                    [
                        -curvature * x * z,
                        law.poisson * curvature * y * z,
                        0.5 * curvature * (x * x + law.poisson * (z * z - y * y)),
                    ]
                })
                .collect();
            errors[kind][level] = (quadratic(&op, &u) / exact - 1.0).abs();
        }
    }
    assert!(errors[1][0] < 0.06, "B-bar coarse errors {errors:?}");
    assert!(errors[1][2] < 0.004, "B-bar fine errors {errors:?}");
    for pair in errors[1].windows(2) {
        assert!(
            pair[1] < 0.26 * pair[0],
            "second-order bending convergence: {errors:?}"
        );
    }
    assert!(
        errors[0][2] > 5.0,
        "fixture must expose full-integration locking: {errors:?}"
    );
}

#[test]
fn native_isotropy_matches_legacy_bbar_and_retains_density_and_identity() {
    assert_eq!(
        ElasticityOptions3::default().formulation,
        ElasticityFormulation3::Displacement
    );
    let law = IsotropicElastic::new(1200.0, 0.4999, 0.1).unwrap();
    let state = ElasticAcousticState::from_isotropic(&law, 1700.0, ContentHash([52; 32])).unwrap();
    let mut original = legacy([2; 3], 0.73, &law, options(0.2));
    let mut sourced = native([2; 3], 0.73, &state, options(0.2)).unwrap();
    assert_eq!(original.nodes(), sourced.nodes());
    let masses = sourced.reference_cell_masses_kg().unwrap().to_vec();
    let scales: Vec<_> = (0..sourced.cells())
        .map(|cell| 0.45 + 0.02 * cell as f64)
        .collect();
    original.set_scales(&scales).unwrap();
    sourced.set_scales(&scales).unwrap();
    same_action(&original, &sourced, &probe(original.nodes()));
    close(masses.iter().sum(), 1700.0 * 0.73, 2e-12);
    assert_eq!(sourced.reference_cell_masses_kg().unwrap(), masses);
    assert_eq!(sourced.formulation(), ElasticityFormulation3::BBar);
    for cell in 0..sourced.cells() {
        assert_eq!(sourced.cell_material_state(cell), Some(&state));
    }
}

#[test]
fn bbar_requires_explicit_isotropy_instead_of_fitting_an_isotropic_proxy() {
    let law = IsotropicElastic::new(1200.0, 0.3, 0.1).unwrap();
    let isotropic =
        ElasticAcousticState::from_isotropic(&law, 1700.0, ContentHash([53; 32])).unwrap();
    for perturb in [false, true] {
        let mut tensor = *isotropic.stiffness_mandel_pa();
        if perturb {
            tensor[0][0] *= 1.1;
        }
        // Even identical isotropic-looking entries retain the generic tensor's
        // admission semantics; numerical resemblance cannot authorize B-bar.
        let generic =
            ElasticAcousticState::try_new_mandel(tensor, 1700.0, ContentHash([54; 32]), 0.1)
                .unwrap();
        assert!(!generic.is_explicitly_isotropic());
        assert!(
            native(
                [1; 3],
                0.73,
                &generic,
                ElasticityOptions3 {
                    ghost_gamma: 0.0,
                    ..Default::default()
                }
            )
            .is_ok()
        );
        assert!(matches!(
            native([1; 3], 0.73, &generic, options(0.0)),
            Err(ElasticityError3::Invalid(what)) if what.contains("explicitly isotropic")
        ));
    }
}

#[test]
fn bbar_scale_derivatives_include_nonzero_cut_ghost_work() {
    let law = IsotropicElastic::new(1200.0, 0.49, 0.1).unwrap();
    let mut stabilized = legacy([2; 3], 0.73, &law, options(0.25));
    let mut bulk_only = legacy([2; 3], 0.73, &law, options(0.0));
    let base: Vec<_> = (0..stabilized.cells())
        .map(|cell| 0.4 + 0.03 * cell as f64)
        .collect();
    stabilized.set_scales(&base).unwrap();
    bulk_only.set_scales(&base).unwrap();
    let u = probe(stabilized.nodes());
    let energy = quadratic(&stabilized, &u);
    assert!(
        energy - quadratic(&bulk_only, &u) > 1e-8,
        "fixture must excite the ghost penalty"
    );
    let derivatives = stabilized.scale_quadratic_forms(&u).unwrap();
    close(
        derivatives.iter().zip(&base).map(|(d, s)| d * s).sum(),
        energy,
        2e-11,
    );
    for (cell, derivative) in derivatives.iter().enumerate() {
        let h = 1e-5;
        let mut scales = base.clone();
        scales[cell] += h;
        stabilized.set_scales(&scales).unwrap();
        let plus = quadratic(&stabilized, &u);
        scales[cell] -= 2.0 * h;
        stabilized.set_scales(&scales).unwrap();
        let minus = quadratic(&stabilized, &u);
        close((plus - minus) / (2.0 * h), *derivative, 2e-7);
    }
    stabilized.set_scales(&base).unwrap();
    assert_eq!(stabilized.formulation(), ElasticityFormulation3::BBar);
}

#[test]
fn adaptive_bbar_matches_cartesian_and_refuses_a_changed_transfer_formulation() {
    let law = IsotropicElastic::new(1200.0, 0.3, 0.1).unwrap();
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|leaf| leaf.index() == [1, 0, 0])
        .unwrap();
    let refined = tree.refined(&[mark], poll).unwrap();
    let coarse = adaptive(&tree, &law, options(0.15));
    let fine = adaptive(&refined, &law, options(0.15));
    let cartesian = legacy([2; 3], 0.73, &law, options(0.15));
    assert_eq!(coarse.nodes(), cartesian.nodes());
    same_action(&coarse, &cartesian, &probe(coarse.nodes()));
    assert_eq!(fine.formulation(), ElasticityFormulation3::BBar);
    assert!(fine.physical_nodes().len() > fine.nodes().len());
    let transfer = AdaptiveTransfer3::new(&coarse, &fine, 10_000, poll).unwrap();
    let affine = |p: &[f64; 3]| [0.01 * p[0], -0.004 * p[0], 0.003 * p[0]];
    let u: Vec<_> = coarse.nodes().iter().flat_map(affine).collect();
    let interpolated = transfer.prolongate(&u, poll).unwrap();
    let physical = fine.physical_displacements(&interpolated).unwrap();
    for (p, u) in fine.physical_nodes().iter().zip(physical.chunks_exact(3)) {
        for (actual, expected) in u.iter().zip(affine(p)) {
            close(*actual, expected, 2e-13);
        }
    }
    let changed = adaptive(
        &refined,
        &law,
        ElasticityOptions3 {
            formulation: ElasticityFormulation3::Displacement,
            ..options(0.15)
        },
    );
    assert!(matches!(
        AdaptiveTransfer3::new(&coarse, &changed, 10_000, poll),
        Err(ElasticityError3::Invalid(_))
    ));
}

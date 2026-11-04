//! G1/G3: reference-mass modes, signed thermal prestress and the complete adjoint.

#[path = "modal/simp.rs"]
mod simp;

use super::*;
use fs_spectral::generalized::GeneralizedError;
use fs_topopt::elasticity::second_order::thermal::modal::{
    ThermalFrequencyGoal, ThermalFrequencyOptions,
};
use fs_topopt::modal::MatrixFreeEigenOptions;

fn modal_options(count: usize) -> ThermalFrequencyOptions {
    ThermalFrequencyOptions {
        solver: MatrixFreeEigenOptions {
            count,
            oversampling: 12,
            max_iterations: 80,
            relative_tolerance: 1e-9,
            linear_tolerance: 1e-11,
            seed: 197,
        },
        aggregation_beta: 0.02,
    }
}

fn frequency<'g, 'a>(
    goal: &'g ThermalComplianceGoal<'a>,
    case: usize,
    mass: &[f64],
    count: usize,
) -> ThermalFrequencyGoal<'g, 'a> {
    let mut callback = |_| ControlFlow::Continue(());
    goal.prestressed_frequency_goal(
        case,
        mass,
        modal_options(count),
        &mut SolveControl::new(SolveBudget::default(), &mut callback),
    )
    .unwrap()
}

// Two identical tetrahedral cantilevers meet at their only free vertex. Their
// support planes are z=+/-1, so the free gradients are exactly +/-e_z.
fn two_sided() -> (DensityElasticity, TetElasticMaterial) {
    let mesh = TetComplex::from_tets(7, vec![[0, 1, 3, 2], [0, 4, 5, 6]]);
    let points = [
        [0., 0., 0.],
        [0., 0., -1.],
        [1., 0., -1.],
        [0., 1., -1.],
        [0., 0., 1.],
        [1., 0., 1.],
        [0., 1., 1.],
    ];
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.0, [71; 32]).unwrap();
    let elasticity = DensityElasticity::with_material_field(
        &mesh,
        &points,
        TetMaterialField::Uniform(&material),
        &|p| p[2] != 0.0,
    )
    .unwrap();
    (elasticity, material)
}

fn two_equal_low_modes(low: f64, high: f64, beta: f64) -> (f64, f64, f64) {
    let tail = (-beta * (high - low)).exp();
    let denominator = 2.0 + tail;
    (
        low - denominator.ln() / beta,
        1.0 / denominator,
        tail / denominator,
    )
}

#[test]
fn traction_prestress_uses_original_mass_and_cancels_false_modulus_dependence() {
    let (elasticity, material) = one_tet(false);
    let thermal = state(&material, [0., 0., 0.01, 0., 0., 0.], 0.1);
    let m = 0.37;
    let r = 0.65;
    for applied in [-0.08, 0.0, 0.08] {
        let mut force = vec![0.0; elasticity.n()];
        force[11] = applied;
        // Zero compliance weight is deliberately retained in physical vibration.
        let loads = [ThermalLoadCase {
            force: &force,
            weight: 0.0,
            strains: TetThermalStrainField::Uniform(&thermal),
        }];
        let goal = solve(&elasticity, &[m], &loads);
        let result = frequency(&goal, 0, &[r], 3);
        assert!(std::ptr::eq(result.source(), &goal));
        assert_eq!(result.case_index(), 0);
        assert_eq!(result.mass_fractions(), &[r]);
        // K_transverse=10m, K_axial=20m, G=F on each free component;
        // M=reference_density*V/10*r=0.1r. Free axial expansion adds no stress.
        let low = (100.0 * m + 10.0 * applied) / r;
        let high = (200.0 * m + 10.0 * applied) / r;
        let (value, w, z) = two_equal_low_modes(low, high, modal_options(3).aggregation_beta);
        for (&actual, expected) in result.modal().values.iter().zip([low, low, high]) {
            close(actual, expected, 1e-9);
        }
        close(result.value(), value, 1e-9);
        // dG/dm=0 for the actual traction state: direct stress and equilibrium
        // derivatives cancel. Freezing equilibrium instead would add F/m.
        close(
            result.modulus_gradient()[0],
            (200.0 * w + 200.0 * z) / r,
            2e-8,
        );
        close(
            result.mass_fraction_gradient()[0],
            -(2.0 * w * low + z * high) / r,
            2e-8,
        );
        assert!(result.adjoint_relative_residual() < options().tolerance);
        for mode in &result.modal().modes {
            for dof in 0..9 {
                assert_eq!(mode[dof], 0.0);
            }
            close(0.1 * r * mode.iter().map(|v| v * v).sum::<f64>(), 1.0, 2e-9);
        }
        assert!(result.modal().relative_residuals.iter().all(|r| *r <= 1e-9));
        // The full three-dimensional trial block is already complete; only
        // the one shared equilibrium adjoint is needed for all three modes.
        assert_eq!(result.modal().work.linear_solves, 1);
    }
}

#[test]
fn restrained_thermal_expansion_softens_modes_without_changing_reference_mass() {
    let (elasticity, material) = two_sided();
    let force = vec![0.0; elasticity.n()];
    let m = 0.6;
    let r = 0.7;
    for alpha in [-0.02, 0.0, 0.02] {
        let thermal = state(&material, [0., 0., alpha, 0., 0., 0.], 0.1);
        let loads = [ThermalLoadCase {
            force: &force,
            weight: 1.0,
            strains: TetThermalStrainField::Uniform(&thermal),
        }];
        let goal = solve(&elasticity, &[m, m], &loads);
        assert!(goal.displacements()[0].iter().all(|u| u.abs() < 1e-13));
        let result = frequency(&goal, 0, &[r, r], 3);
        // Each cell has sigma_zz=-120m alpha. G=-40m alpha, whereas
        // K=(20m,20m,40m) and M=0.2r at the shared vertex.
        let low = (100.0 - 200.0 * alpha) * m / r;
        let high = (200.0 - 200.0 * alpha) * m / r;
        let (value, w, z) = two_equal_low_modes(low, high, modal_options(3).aggregation_beta);
        for (&actual, expected) in result.modal().values.iter().zip([low, low, high]) {
            close(actual, expected, 1e-9);
        }
        close(result.value(), value, 1e-9);
        let uniform_m = (2.0 * w * low + z * high) / m;
        let uniform_r = -(2.0 * w * low + z * high) / r;
        close(result.modulus_gradient().iter().sum(), uniform_m, 2e-8);
        close(
            result.mass_fraction_gradient().iter().sum(),
            uniform_r,
            2e-8,
        );
        for mode in &result.modal().modes {
            close(
                0.2 * r * mode[..3].iter().map(|v| v * v).sum::<f64>(),
                1.0,
                2e-9,
            );
            assert!(mode[3..].iter().all(|u| *u == 0.0));
        }
    }
}

#[test]
fn heterogeneous_full_tensor_frequency_adjoint_matches_both_independent_design_fields() {
    let (_, _, elasticity, thermal, mut a, b) = cube();
    for f in &mut a {
        *f *= 8.0;
    }
    let loads = [
        ThermalLoadCase {
            force: &a,
            weight: 0.0,
            strains: TetThermalStrainField::PerElement(&thermal),
        },
        ThermalLoadCase {
            force: &b,
            weight: 1.0,
            strains: TetThermalStrainField::PerElement(&thermal),
        },
    ];
    let m: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.48 + 0.04 * i as f64)
        .collect();
    let r: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.72 - 0.03 * i as f64)
        .collect();
    let goal = solve(&elasticity, &m, &loads);
    let count = elasticity.free().iter().filter(|f| **f).count();
    for case in 0..2 {
        let actual = frequency(&goal, case, &r, count);
        for cell in 0..m.len() {
            let h = 2e-5;
            let mut plus = m.clone();
            plus[cell] += h;
            let mut minus = m.clone();
            minus[cell] -= h;
            let plus_goal = solve(&elasticity, &plus, &loads);
            let minus_goal = solve(&elasticity, &minus, &loads);
            let fd = (frequency(&plus_goal, case, &r, count).value()
                - frequency(&minus_goal, case, &r, count).value())
                / (2.0 * h);
            close(actual.modulus_gradient()[cell], fd, 3e-6);
            let mut plus = r.clone();
            plus[cell] += h;
            let mut minus = r.clone();
            minus[cell] -= h;
            let fd = (frequency(&goal, case, &plus, count).value()
                - frequency(&goal, case, &minus, count).value())
                / (2.0 * h);
            close(actual.mass_fraction_gradient()[cell], fd, 3e-6);
        }
        // Complete spectrum allows a source-independent check of mass
        // orthogonality using the original unit-density element mass blocks.
        let density: Vec<_> = r
            .iter()
            .enumerate()
            .map(|(e, f)| f * elasticity.material_density_kg_m3(e).unwrap())
            .collect();
        for (i, mode) in actual.modal().modes.iter().enumerate() {
            let mut image = vec![0.0; elasticity.n()];
            elasticity.apply_mass(&density, mode, &mut image);
            for (j, other) in actual.modal().modes.iter().enumerate() {
                let product = other.iter().zip(&image).map(|(a, b)| a * b).sum::<f64>();
                close(product, if i == j { 1.0 } else { 0.0 }, 2e-5);
            }
        }
    }
}

#[test]
fn nonpositive_tangent_invalid_mass_and_stopped_work_do_not_publish_frequency_results() {
    let (elasticity, material) = two_sided();
    let force = vec![0.0; elasticity.n()];
    let thermal = state(&material, [0., 0., 0.01, 0., 0., 0.], 1.0);
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 1.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    let goal = solve(&elasticity, &[0.6, 0.7], &loads);
    let mut keep = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut keep);
    for mass in [
        &[0.5][..],
        &[-0.1, 0.5],
        &[0.0, 0.0],
        &[1.1, 0.5],
        &[f64::NAN, 0.5],
    ] {
        assert!(
            goal.prestressed_frequency_goal(0, mass, modal_options(3), &mut control)
                .is_err()
        );
    }
    assert!(
        goal.prestressed_frequency_goal(1, &[0.6, 0.6], modal_options(3), &mut control)
            .is_err()
    );
    let mut stop = |_| ControlFlow::Break(());
    let mut cancelled = SolveControl::new(SolveBudget::default(), &mut stop);
    assert!(matches!(
        goal.prestressed_frequency_goal(0, &[0.6, 0.6], modal_options(3), &mut cancelled),
        Err(GeneralizedError::Operator(ComplianceError::Stopped(
            EvaluationStop::Cancelled
        )))
    ));
    assert_eq!(cancelled.work().linear_solves, 0);

    // A complete trial space still needs the equilibrium adjoint for gradients.
    // Its budget refusal must not mutate the retained thermal point.
    let mut limited_keep = |_| ControlFlow::Continue(());
    let mut limited = SolveControl::new(
        SolveBudget {
            per_solve_iterations: 0,
            total_iterations: 0,
        },
        &mut limited_keep,
    );
    assert!(
        goal.prestressed_frequency_goal(0, &[0.6, 0.6], modal_options(3), &mut limited)
            .is_err()
    );
    let before = frequency(&goal, 0, &[0.6, 0.6], 3).value();
    close(frequency(&goal, 0, &[0.6, 0.6], 3).value(), before, 1e-12);

    // Deliberately admitted large small-strain fixture exposes negative signed
    // geometric stiffness; the modal owner must refuse instead of taking abs.
    let hot = state(&material, [0., 0., 0.6, 0., 0., 0.], 1.0);
    let hot_loads = [ThermalLoadCase {
        force: &force,
        weight: 0.0,
        strains: TetThermalStrainField::Uniform(&hot),
    }];
    let hot_goal = solve(&elasticity, &[0.6, 0.6], &hot_loads);
    assert!(matches!(
        hot_goal.prestressed_frequency_goal(0, &[0.6, 0.6], modal_options(3), &mut control),
        Err(GeneralizedError::Operator(ComplianceError::InvalidInput(_)))
    ));
}

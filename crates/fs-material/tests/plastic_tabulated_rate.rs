//! Independent consistency, tensor derivatives, dissipation and held-history
//! relaxation for the measured-curve Perzyna model. Synthetic coefficients
//! test the numerical law, not calibration of a particular alloy.

use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{LinearOverstressJ2, PowerLawOverstressJ2};
use fs_material::plastic::tabulated::rate::TabulatedOverstressJ2;
use fs_material::plastic::tabulated::{
    HardeningPoint, TabulatedError, TabulatedHardening, TabulatedJ2,
};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const MU: f64 = 500.0;
const HK: f64 = 160.0;
const ETA: f64 = 40.0;
const S0: f64 = 2.0;
const DT: f64 = 0.2;
const KNOTS: [(f64, f64); 5] = [
    (0.0, 3.0),
    (0.002, 4.0),
    (0.006, 4.0),
    (0.012, 7.0),
    (0.04, 8.4),
];
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn native(hardening: f64, kinematic: f64) -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        3.0,
        hardening,
        kinematic,
    )
    .unwrap()
}
fn table(points: &[(f64, f64)]) -> TabulatedHardening {
    TabulatedHardening::new(
        points
            .iter()
            .map(|&(plastic_strain, yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .collect(),
    )
    .unwrap()
}
fn backbone(kinematic: f64) -> TabulatedJ2 {
    TabulatedJ2::new(native(0.0, kinematic), table(&KNOTS)).unwrap()
}
fn law(exponent: f64, duration: f64) -> TabulatedOverstressJ2 {
    TabulatedOverstressJ2::new_power_law(backbone(HK), ETA, S0, exponent, duration).unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}, difference {:.8e}",
        (actual - expected).abs()
    );
}
fn deviatoric(axial: f64) -> Voigt {
    [axial, -0.5 * axial, -0.5 * axial, 0.0, 0.0, 0.0]
}
fn radius(alpha: f64) -> f64 {
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= b {
            return ya + (yb - ya) * ((alpha - a) / (b - a));
        }
    }
    panic!("independent curve oracle outside support");
}
fn storage(alpha: f64) -> f64 {
    let mut result = 0.0;
    for pair in KNOTS.windows(2) {
        let [(a, ya), (b, yb)] = [pair[0], pair[1]];
        if alpha <= a {
            break;
        }
        let width = alpha.min(b) - a;
        result += (ya - 3.0) * width + 0.5 * (yb - ya) / (b - a) * width * width;
    }
    result
}
fn viscous(delta: f64, exponent: f64, duration: f64) -> f64 {
    S0 * (ETA * delta / (duration * S0)).powf(1.0 / exponent)
}

#[test]
fn scalar_returns_cross_knots_and_plateaus_with_independent_rate_equations() {
    for exponent in [1.0, 1.7, 4.0] {
        let material = law(exponent, DT);
        for old in [0.0, 0.001, 0.002, 0.005, 0.011, 0.012] {
            for target in [0.0015, 0.002, 0.003, 0.006, 0.009, 0.012, 0.019, 0.04] {
                if target <= old {
                    continue;
                }
                let increment = target - old;
                for correction in [E, 3.0 * MU] {
                    let trial = radius(target)
                        + (correction + HK) * increment
                        + viscous(increment, exponent, DT);
                    let actual = material.plastic_increment(trial, old, correction).unwrap();
                    near(actual, increment, 3e-16, 2e-13);
                    near(
                        material.endpoint_overstress_pa(actual).unwrap(),
                        viscous(increment, exponent, DT),
                        2e-13,
                        2e-13,
                    );
                }
            }
        }
        assert_eq!(material.plastic_increment(2.0, 0.0, 3.0 * MU).unwrap(), 0.0);
        assert_eq!(material.endpoint_overstress_pa(0.0).unwrap(), 0.0);
        assert_eq!(material.duration_s(), DT);
        assert_eq!(material.stress_scale_pa(), Some(S0));
    }
}

#[test]
fn native_response_separates_table_storage_yield_viscosity_and_integration_loss() {
    for exponent in [1.0, 3.2] {
        let material = law(exponent, DT);
        for alpha in [0.001, 0.004, 0.008, 0.014, 0.022] {
            let overstress = viscous(alpha, exponent, DT);
            let total = ((3.0 * MU + HK) * alpha + radius(alpha) + overstress) / (3.0 * MU);
            let strain = deviatoric(total);
            let state = material.initial_state();
            let response = material.response(&strain, &state).unwrap();
            let increment = material
                .advance_with_energy(&[0.0; 6], &state, &strain, TOL)
                .unwrap();
            near(response.state.alpha, alpha, 3e-16, 2e-13);
            near(response.state.plastic_strain[0], alpha, 3e-16, 2e-13);
            near(trace(&response.state.plastic_strain), 0.0, 3e-17, 0.0);
            near(
                response.stress_pa[0],
                (2.0 / 3.0) * (radius(alpha) + HK * alpha + overstress),
                3e-13,
                3e-13,
            );
            near(
                increment.energy.final_state.hardening_j_m3,
                storage(alpha),
                3e-16,
                3e-13,
            );
            near(
                increment.energy.final_state.kinematic_hardening_j_m3,
                0.5 * HK * alpha * alpha,
                3e-16,
                3e-13,
            );
            near(
                increment.energy.physical_dissipation_j_m3,
                (3.0 + overstress) * alpha,
                3e-16,
                3e-13,
            );
            let elastic = deviatoric(total - alpha);
            let elastic_stress = native(0.0, HK).elastic.stress(&elastic, &());
            let expected_numerical = 0.5 * contract(&elastic, &elastic_stress)
                + 0.5 * HK * alpha * alpha
                + (radius(alpha) - 3.0) * alpha
                - storage(alpha);
            near(
                increment.energy.numerical_dissipation_j_m3,
                expected_numerical,
                3e-15,
                3e-12,
            );
            near(increment.energy.residual_j_m3, 0.0, 5e-14, 0.0);
            let (stress, potential) = material.stress_and_potential(&strain, &state).unwrap();
            assert_eq!(stress, response.stress_pa);
            near(
                potential,
                increment.energy.final_state.total_j_m3()
                    + (3.0 + exponent / (exponent + 1.0) * overstress) * alpha,
                3e-15,
                3e-12,
            );
            // A viscous endpoint is admissible even though its equilibrium
            // yield residual is positive. Its history must remain usable.
            material
                .validate_history(&strain, &increment.state)
                .unwrap();
            assert!(
                material
                    .backbone()
                    .validate_committed_history(&strain, &increment.state, 1e-9)
                    .is_err()
            );
        }
    }
}

fn derivatives(material: &TabulatedOverstressJ2, strain: &Voigt, old: &J2State) {
    let response = material.response(strain, old).unwrap();
    let h = 1e-7;
    for column in 0..6 {
        let mut plus = *strain;
        let mut minus = *strain;
        plus[column] += h;
        minus[column] -= h;
        let (sp, pp) = material.stress_and_potential(&plus, old).unwrap();
        let (sm, pm) = material.stress_and_potential(&minus, old).unwrap();
        let weight = if column < 3 { 1.0 } else { 2.0 };
        near(
            (pp - pm) / (2.0 * h),
            weight * response.stress_pa[column],
            3e-8,
            2e-7,
        );
        for row in 0..6 {
            near(
                (sp[row] - sm[row]) / (2.0 * h),
                response.tangent_pa[row][column],
                3e-5,
                2e-7,
            );
        }
    }
}

#[test]
fn noncoaxial_prager_history_retains_full_tangent_potential_and_energy_on_reversal() {
    for exponent in [1.0, 3.2] {
        let mut state = law(exponent, DT).initial_state();
        let mut old_strain = [0.0; 6];
        for (duration, strain) in [
            (0.2, [0.0002, -0.0001, -0.0001, 0.0001, 0.0, 0.0]),
            (0.4, [0.009, -0.003, -0.006, 0.002, 0.001, 0.0]),
            (0.1, [0.0085, -0.0028, -0.0057, 0.0019, 0.0009, 0.0]),
            (0.25, [-0.004, 0.001, 0.003, -0.002, 0.002, 0.001]),
            (0.35, [0.010, -0.006, -0.004, 0.005, -0.002, 0.001]),
        ] {
            let material = law(exponent, duration);
            derivatives(&material, &strain, &state);
            let retained = state.clone();
            let step = material
                .advance_with_energy(&old_strain, &state, &strain, TOL)
                .unwrap();
            assert_eq!(state, retained);
            assert!(step.state.alpha >= state.alpha);
            assert!(step.energy.physical_dissipation_j_m3 >= 0.0);
            assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
            let delta = step.state.alpha - state.alpha;
            near(
                step.energy.physical_dissipation_j_m3,
                (3.0 + viscous(delta, exponent, duration)) * delta,
                3e-15,
                3e-12,
            );
            near(step.energy.residual_j_m3, 0.0, 5e-14, 0.0);
            near(
                material
                    .stored_energy(&strain, &step.state)
                    .unwrap()
                    .total_j_m3(),
                step.energy.final_state.total_j_m3(),
                3e-15,
                3e-12,
            );
            state = step.state;
            old_strain = strain;
        }
    }
}

#[test]
fn affine_and_perfect_tables_reduce_to_the_original_native_rate_models() {
    for hardening in [0.0, 180.0] {
        for exponent in [1.0, 3.2] {
            let backbone = TabulatedJ2::new(
                native(0.0, HK),
                table(&[(0.0, 3.0), (0.3, 3.0 + 0.3 * hardening)]),
            )
            .unwrap();
            let material =
                TabulatedOverstressJ2::new_power_law(backbone.clone(), ETA, S0, exponent, DT)
                    .unwrap();
            let original =
                PowerLawOverstressJ2::new(native(hardening, HK), ETA, S0, exponent, DT).unwrap();
            let mut state = material.initial_state();
            let mut old = [0.0; 6];
            for strain in [
                deviatoric(0.009),
                deviatoric(0.012),
                deviatoric(0.011),
                deviatoric(-0.004),
            ] {
                let a = material
                    .advance_with_energy(&old, &state, &strain, TOL)
                    .unwrap();
                let b = original
                    .advance_with_energy(&old, &state, &strain, TOL)
                    .unwrap();
                let ta = material.tangent(&strain, &state);
                let tb = original.tangent(&strain, &state);
                near(a.state.alpha, b.state.alpha, 3e-16, 2e-12);
                near(
                    a.energy.final_state.total_j_m3(),
                    b.energy.final_state.total_j_m3(),
                    3e-15,
                    2e-12,
                );
                near(
                    a.energy.physical_dissipation_j_m3,
                    b.energy.physical_dissipation_j_m3,
                    3e-15,
                    2e-12,
                );
                near(
                    a.energy.numerical_dissipation_j_m3,
                    b.energy.numerical_dissipation_j_m3,
                    3e-15,
                    2e-12,
                );
                for i in 0..6 {
                    near(a.stress_pa[i], b.stress_pa[i], 3e-13, 2e-12);
                    near(
                        a.state.plastic_strain[i],
                        b.state.plastic_strain[i],
                        3e-16,
                        2e-12,
                    );
                    for j in 0..6 {
                        near(ta[i][j], tb[i][j], 3e-10, 2e-12);
                    }
                }
                old = strain;
                state = a.state;
            }
            if exponent == 1.0 {
                let linear = TabulatedOverstressJ2::new_linear(backbone, ETA, DT).unwrap();
                let original = LinearOverstressJ2::new(native(hardening, HK), ETA, DT).unwrap();
                let strain = deviatoric(0.009);
                for (a, b) in linear
                    .stress(&strain, &linear.initial_state())
                    .iter()
                    .zip(original.stress(&strain, &original.initial_state()))
                {
                    near(*a, b, 3e-13, 2e-12);
                }
                assert_eq!(linear.stress_scale_pa(), None);
            }
        }
    }
}

#[test]
fn held_strain_relaxes_on_a_measured_plateau_without_erasing_permanent_deformation() {
    let material = TabulatedOverstressJ2::new_linear(backbone(HK), ETA, DT).unwrap();
    let strain = deviatoric(0.009);
    let mut step = material
        .advance_with_energy(&[0.0; 6], &material.initial_state(), &strain, TOL)
        .unwrap();
    let correction = 3.0 * MU + HK;
    let rate_modulus = ETA / DT;
    let mut alpha = (3.0 * MU * 0.009 - 4.0) / (correction + rate_modulus);
    near(step.state.alpha, alpha, 3e-16, 2e-13);
    let mut excess = rate_modulus * alpha;
    for _ in 0..5 {
        let old = step.state.clone();
        let delta = excess / (correction + rate_modulus);
        alpha += delta;
        excess = rate_modulus * delta;
        step = material
            .advance_with_energy(&strain, &old, &strain, TOL)
            .unwrap();
        near(step.state.alpha, alpha, 3e-16, 2e-13);
        assert!(step.state.alpha > old.alpha && step.state.alpha < 0.006);
        near(
            material
                .backbone()
                .yield_function_with_state(&step.stress_pa, &step.state),
            excess,
            3e-13,
            2e-12,
        );
        near(step.energy.mechanical_work_j_m3, 0.0, 0.0, 0.0);
        near(
            -step.energy.stored_change_j_m3(),
            step.energy.physical_dissipation_j_m3 + step.energy.numerical_dissipation_j_m3,
            3e-14,
            2e-12,
        );
    }
    // Fully release total stress. The remaining Prager translation lies well
    // inside the supported surface; no sub-yield creep is introduced.
    let released_strain = step.state.plastic_strain;
    let unloaded = material
        .advance_with_energy(&strain, &step.state, &released_strain, TOL)
        .unwrap();
    assert_eq!(unloaded.state, step.state);
    near(unloaded.energy.physical_dissipation_j_m3, 0.0, 0.0, 0.0);
    let held = material
        .advance_with_energy(&released_strain, &unloaded.state, &released_strain, TOL)
        .unwrap();
    assert_eq!(held.state, unloaded.state);
    assert_eq!(held.energy.physical_dissipation_j_m3, 0.0);
    assert!(held.state.alpha > 0.005);
}

#[test]
fn viscosity_can_keep_a_return_inside_support_and_true_domain_exits_are_atomic() {
    let hardening = table(&[(0.0, 3.0), (0.002, 4.0), (0.005, 4.0)]);
    let base = TabulatedJ2::new(native(0.0, HK), hardening).unwrap();
    for exponent in [1.0, 3.2] {
        let duration = 0.02;
        let material =
            TabulatedOverstressJ2::new_power_law(base.clone(), ETA, S0, exponent, duration)
                .unwrap();
        let alpha = 0.00475;
        let trial = 4.0 + (3.0 * MU + HK) * alpha + viscous(alpha, exponent, duration);
        let strain = deviatoric(trial / (3.0 * MU));
        let old = material.initial_state();
        assert!(base.response(&strain, &old).is_err());
        let actual = material.response(&strain, &old).unwrap();
        near(actual.state.alpha, alpha, 3e-16, 2e-13);
        let unsupported = 0.0052;
        let excess_trial =
            4.0 + (3.0 * MU + HK) * unsupported + viscous(unsupported, exponent, duration);
        assert!(matches!(
            material.plastic_increment(excess_trial, 0.0, 3.0 * MU),
            Err(TabulatedError::Domain { .. })
        ));
        assert!(
            material
                .response(&deviatoric(excess_trial / (3.0 * MU)), &old)
                .is_err()
        );
        assert_eq!(old, material.initial_state());
        let mut invalid = old;
        invalid.alpha = 0.0051;
        assert!(material.validate_history(&[0.0; 6], &invalid).is_err());
    }
    assert!(TabulatedOverstressJ2::new_linear(base.clone(), ETA, 0.0).is_err());
    assert!(TabulatedOverstressJ2::new_power_law(base.clone(), ETA, 0.0, 2.0, DT).is_err());
    assert!(TabulatedOverstressJ2::new_power_law(base.clone(), ETA, S0, 0.99, DT).is_err());
    assert!(TabulatedOverstressJ2::new_power_law(base, -1.0, S0, 2.0, DT).is_err());
}

#[test]
fn final_measured_rate_endpoint_remains_inside_support_in_every_native_return() {
    let old = J2State {
        plastic_strain: [0.0; 6],
        alpha: 0.012,
    };
    let delta = 0.04 - old.alpha;
    for exponent in [1.0, 3.2] {
        let material = law(exponent, DT);
        let trial = radius(0.04) + (3.0 * MU + HK) * delta + viscous(delta, exponent, DT);
        let axial = (trial / (3.0 * MU)).next_up();
        let strain = deviatoric(axial);
        let response = material.response(&strain, &old).unwrap();
        let step = material
            .advance_with_energy(&[0.0; 6], &old, &strain, TOL)
            .unwrap();
        assert!(response.state.alpha <= 0.04 && step.state.alpha <= 0.04);
        near(response.state.alpha, 0.04, 3e-16, 2e-13);
        near(step.state.alpha, response.state.alpha, 3e-16, 2e-13);
        near(step.energy.residual_j_m3, 0.0, 1e-13, 0.0);
        assert_eq!(material.update_state(&strain, &old), response.state);
        assert!(material.response(&deviatoric(axial + 1e-6), &old).is_err());
        assert_eq!(old.alpha, 0.012);
    }
}

#[test]
fn a_kinetic_upper_bound_resolves_very_slow_flow_with_a_fixed_iteration_budget() {
    let viscosity = 1e120;
    let target_overstress = 0.4;
    for exponent in [1.1_f64, 3.2, 40.0] {
        let material =
            TabulatedOverstressJ2::new_power_law(backbone(HK), viscosity, S0, exponent, DT)
                .unwrap();
        let expected = (DT * S0 / viscosity) * (target_overstress / S0).powf(exponent);
        // At this scale the elastic/hardening corrections are far below the
        // driving stress's resolution; the independent inverse kinetic law
        // gives the representable positive flow instead of an artificial zero.
        let strain = deviatoric((3.0 + target_overstress) / (3.0 * MU));
        let response = material
            .response(&strain, &material.initial_state())
            .unwrap();
        near(response.state.alpha, expected, 0.0, 5e-12);
        near(
            material
                .endpoint_overstress_pa(response.state.alpha)
                .unwrap(),
            target_overstress,
            3e-13,
            3e-13,
        );
        assert!(response.state.alpha > 0.0);
    }
}

//! Free-sided primary resistance recovery against independent scalar references.
//! Manufactured coefficients establish implementation behavior, not alloy data.

use fs_material::creep::recovery::uniaxial::UniaxialRecoveringNorton;
use fs_material::creep::recovery::{CreepRecoveryParameters, RecoveringNortonCreep};
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::creep::{
    CreepEnergy, CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::{IsotropicElastic, SmallStrainLaw, Uniaxial};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const STRESS: f64 = 3e6;
const RATE: f64 = 2e-3;
const A: f64 = 8e-4;

fn primary(n: f64, dt: f64) -> NortonCreep {
    primary_at_poisson(n, dt, NU)
}

fn primary_at_poisson(n: f64, dt: f64, poisson: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(E, poisson, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: A,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap()
}

fn law(n: f64, dt: f64, recovery_rate: f64) -> UniaxialRecoveringNorton {
    UniaxialRecoveringNorton::new(
        RecoveringNortonCreep::new(
            primary(n, dt),
            CreepRecoveryParameters {
                reference_rate_s: recovery_rate,
                activation_energy_j_mol: 0.0,
                reference_temperature_k: 500.0,
            },
        )
        .unwrap(),
    )
}

fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 3e-6,
        energy_j_m3: 2e-7,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn close_energy(actual: CreepEnergy, expected: CreepEnergy) {
    for (actual, expected) in [
        (actual.initial_stored_j_m3, expected.initial_stored_j_m3),
        (actual.final_stored_j_m3, expected.final_stored_j_m3),
        (actual.mechanical_work_j_m3, expected.mechanical_work_j_m3),
        (actual.physical_dissipation_j_m3, expected.physical_dissipation_j_m3),
        (actual.creep_dissipation_j_m3, expected.creep_dissipation_j_m3),
        (actual.plastic_dissipation_j_m3, expected.plastic_dissipation_j_m3),
        (actual.numerical_dissipation_j_m3, expected.numerical_dissipation_j_m3),
        (actual.thermoelastic_storage_change_j_m3, expected.thermoelastic_storage_change_j_m3),
        (actual.reversible_heat_j_m3, expected.reversible_heat_j_m3),
        (actual.residual_j_m3, expected.residual_j_m3),
    ] {
        close(actual, expected, tolerances().energy_j_m3);
    }
}

#[test]
fn free_sided_recovery_matches_quadratic_flow_tangent_and_native_tensor_heat() {
    let dt = 0.7;
    let rate = 1.3;
    let law = law(1.0, dt, rate);
    let old_strain = 0.002;
    let old = law
        .response(old_strain, &law.initial_state())
        .unwrap()
        .state;
    let frozen = old.clone();
    for strain in [-0.003, 0.006] {
        let trial = E * (strain - old.creep_strain());
        let divisor = 1.0 + rate * dt;
        let effective_a = A * divisor;
        let coefficient = dt * RATE * effective_a / STRESS;
        let b = effective_a + old.resistance_strain() + E * coefficient;
        let c = coefficient * trial.abs();
        let delta = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
        let stress = trial.signum() * (trial.abs() - E * delta);
        let tangent =
            E / (1.0 + E * coefficient / (effective_a + old.resistance_strain() + 2.0 * delta));
        let response = law.response(strain, &old).unwrap();
        let checked = law
            .advance_with_energy(old_strain, &old, strain, tolerances())
            .unwrap();
        close(response.stress_pa, stress, 3e-8);
        close(response.tangent_pa, tangent, 3e-5);
        close(
            response.state.creep_strain(),
            old.creep_strain() + trial.signum() * delta,
            3e-17,
        );
        close(
            response.state.equivalent_creep_strain(),
            old.equivalent_creep_strain() + delta,
            3e-17,
        );
        close(
            response.state.resistance_strain(),
            (old.resistance_strain() + delta) / divisor,
            3e-17,
        );
        close(
            response.strain[1],
            -NU * stress / E - 0.5 * response.state.creep_strain(),
            3e-17,
        );
        assert_eq!(response.strain[1], response.strain[2]);
        assert_eq!(checked.state, response.state);
        close(
            checked.energy.creep_dissipation_j_m3,
            stress.abs() * delta,
            1e-8,
        );
        assert_eq!(
            checked.energy.creep_dissipation_j_m3,
            checked.energy.physical_dissipation_j_m3
        );
        assert_eq!(checked.energy.plastic_dissipation_j_m3, 0.0);
        close(checked.energy.residual_j_m3, 0.0, tolerances().energy_j_m3);
        let native = law
            .backbone()
            .advance_with_energy(
                &law.tensor_strain(old_strain, &old).unwrap(),
                &old.tensor_history(),
                &response.strain,
                tolerances(),
            )
            .unwrap();
        // The axial and tensor returns solve different scalar equations for
        // the same endpoint. Their rounding need not be bit-identical; compare
        // stress and the complete work/heat account within numerical budgets.
        close_energy(checked.energy, native.energy);
        for (actual, expected) in checked.stress_pa.into_iter().zip(native.stress_pa) {
            close(actual, expected, tolerances().stress_pa);
        }
        for transverse in &checked.stress_pa[1..] {
            close(*transverse, 0.0, 3e-8);
        }
    }
    assert_eq!(old, frozen);
    assert_eq!(
        law.kinetics().strain_hardening().unwrap().reference_strain,
        A
    );
    assert_eq!(law.kinetics(), law.backbone().norton());
}

#[test]
fn nonlinear_recovery_tangent_is_the_native_free_sided_schur_complement() {
    for exponent in [1.0, 3.2, 8.0] {
        let native = primary(exponent, 0.4)
            .with_strain_hardening(StrainHardeningCreep {
                reference_strain: A,
                exponent: 0.8,
                residual_rate_fraction: 0.06,
            })
            .unwrap();
        let law = UniaxialRecoveringNorton::new(
            RecoveringNortonCreep::new(
                native,
                CreepRecoveryParameters {
                    reference_rate_s: 1.7,
                    activation_energy_j_mol: 0.0,
                    reference_temperature_k: 500.0,
                },
            )
            .unwrap(),
        );
        let old = law.response(0.004, &law.initial_state()).unwrap().state;
        let frozen = old.clone();
        for strain in [-0.003, 0.0, 0.006] {
            let response = law.response(strain, &old).unwrap();
            let h = 1e-8;
            let difference =
                (law.stress(strain + h, &old) - law.stress(strain - h, &old)) / (2.0 * h);
            close(response.tangent_pa, difference, 2e-5 * E);
            let native = law
                .backbone()
                .tangent(&response.strain, &old.tensor_history());
            let determinant = native[1][1] * native[2][2] - native[1][2] * native[2][1];
            let dy = (native[2][2] * native[1][0] - native[1][2] * native[2][0]) / determinant;
            let dz = (native[1][1] * native[2][0] - native[2][1] * native[1][0]) / determinant;
            let condensed = native[0][0] - native[0][1] * dy - native[0][2] * dz;
            close(response.tangent_pa, condensed, 3e-5);
            assert_eq!(response, law.response(strain, &old).unwrap());
        }
        assert_eq!(old, frozen);
    }
}

#[test]
fn checked_recovery_trajectory_preserves_the_independent_tensor_return() {
    let recovery_rate = 1.3;
    for poisson in [-0.8, 0.2, 0.49] {
        for exponent in [1.0, 3.0] {
            for duration in [1e-3, 0.2, 5.0] {
                let native = primary_at_poisson(exponent, duration, poisson)
                    .with_strain_hardening(StrainHardeningCreep {
                        reference_strain: A,
                        exponent: 0.8,
                        residual_rate_fraction: 0.06,
                    })
                    .unwrap();
                let tensor_law = RecoveringNortonCreep::new(
                    native,
                    CreepRecoveryParameters {
                        reference_rate_s: recovery_rate,
                        activation_energy_j_mol: 0.0,
                        reference_temperature_k: 500.0,
                    },
                )
                .unwrap();
                let axial = UniaxialRecoveringNorton::new(tensor_law.clone());
                let mut axial_state = axial.initial_state();
                let mut tensor_state = tensor_law.initial_state();
                let mut old_axial_strain = 0.0;
                let mut old_tensor_strain = [0.0; 6];
                // Numerical equivalence, not identical roots or a constitutive
                // accuracy claim: scale history by E, tangent by elasticity.
                let strain_tolerance = tolerances().stress_pa / E;
                let tangent_tolerance =
                    1e-11 * E / (1.0 - 2.0 * poisson).min(1.0 + poisson);
                for target in [Some(0.004), Some(0.004), Some(-0.003), None, None, Some(0.005)] {
                    let strain = target.unwrap_or(axial_state.creep_strain());
                    let response = axial.response(strain, &axial_state).unwrap();
                    let checked = axial
                        .advance_with_energy(old_axial_strain, &axial_state, strain, tolerances())
                        .unwrap();
                    let tensor = tensor_law
                        .advance_with_energy(
                            &old_tensor_strain,
                            &tensor_state,
                            &response.strain,
                            tolerances(),
                        )
                        .unwrap();
                    let history = checked.state.tensor_history();
                    for i in 0..6 {
                        close(checked.stress_pa[i], tensor.stress_pa[i], tolerances().stress_pa);
                        close(
                            history.native.creep_strain[i],
                            tensor.state.native.creep_strain[i],
                            strain_tolerance,
                        );
                        for j in 0..6 {
                            close(
                                checked.tangent_pa[i][j],
                                tensor.tangent_pa[i][j],
                                tangent_tolerance,
                            );
                        }
                    }
                    close(
                        history.native.equivalent_creep_strain,
                        tensor.state.native.equivalent_creep_strain,
                        strain_tolerance,
                    );
                    close(
                        history.resistance_strain,
                        tensor.state.resistance_strain,
                        strain_tolerance,
                    );
                    close_energy(checked.energy, tensor.energy);
                    if target.is_none() {
                        assert_eq!(response.stress_pa, 0.0);
                        assert_eq!(checked.state.native(), axial_state.native());
                        close(
                            checked.state.resistance_strain(),
                            axial_state.resistance_strain() / (1.0 + recovery_rate * duration),
                            strain_tolerance,
                        );
                        assert_eq!(checked.energy.physical_dissipation_j_m3, 0.0);
                    }
                    old_axial_strain = strain;
                    old_tensor_strain = tensor.strain;
                    axial_state = checked.state;
                    tensor_state = tensor.state;
                }
            }
        }
    }
}

#[test]
fn zero_recovery_is_exact_original_axial_primary_creep_through_reversals() {
    let recovering = law(3.2, 0.4, 0.0);
    let native = UniaxialNorton::new(recovering.kinetics().clone());
    let mut state = recovering.initial_state();
    let mut old_strain = 0.0;
    for strain in [0.004, 0.004, -0.003, 0.006] {
        let expected = native.response(strain, state.native()).unwrap();
        let response = recovering.response(strain, &state).unwrap();
        assert_eq!(response.state.native(), &expected.state);
        assert_eq!(
            response.state.resistance_strain(),
            expected.state.equivalent_creep_strain()
        );
        assert_eq!(response.stress_pa, expected.stress_pa);
        assert_eq!(response.tangent_pa, expected.tangent_pa);
        assert_eq!(response.strain, expected.strain);
        let expected = native
            .advance_with_energy(old_strain, state.native(), strain, tolerances())
            .unwrap();
        let checked = recovering
            .advance_with_energy(old_strain, &state, strain, tolerances())
            .unwrap();
        assert_eq!(checked.energy, expected.energy);
        assert_eq!(checked.stress_pa, expected.stress_pa);
        assert_eq!(checked.tangent_pa, expected.tangent_pa);
        state = response.state;
        old_strain = strain;
    }
}

#[test]
fn unloaded_dwell_recovers_only_resistance_and_renews_flow_on_reload() {
    let law = law(3.0, 0.2, 1.0);
    let loaded = law.response(0.004, &law.initial_state()).unwrap().state;
    let free = loaded.creep_strain();
    let mut rested = loaded.clone();
    for step in 1..=8 {
        let checked = law
            .advance_with_energy(free, &rested, free, tolerances())
            .unwrap();
        assert_eq!(checked.state.native(), loaded.native());
        close(
            checked.state.resistance_strain(),
            loaded.resistance_strain() / 1.2_f64.powi(step),
            2e-18,
        );
        assert_eq!(checked.stress_pa, [0.0; 6]);
        assert_eq!(checked.energy.creep_dissipation_j_m3, 0.0);
        assert_eq!(checked.energy.mechanical_work_j_m3, 0.0);
        assert_eq!(checked.energy.numerical_dissipation_j_m3, 0.0);
        rested = checked.state;
    }
    let immediate = law.response(free + 0.003, &loaded).unwrap();
    let after_dwell = law.response(free + 0.003, &rested).unwrap();
    assert!(
        after_dwell.state.equivalent_creep_strain() > immediate.state.equivalent_creep_strain()
    );
    assert!(after_dwell.stress_pa < immediate.stress_pa);
    assert!(after_dwell.state.equivalent_creep_strain() > loaded.equivalent_creep_strain());
    assert!(after_dwell.state.resistance_strain() < after_dwell.state.equivalent_creep_strain());
}

#[test]
fn exact_flow_increment_survives_roundoff_in_accumulated_and_signed_history() {
    let loading = law(1.0, 0.2, 1.0);
    let loaded = loading
        .response(0.004, &loading.initial_state())
        .unwrap()
        .state;
    let free = loaded.creep_strain();
    // A strong stress-free dwell leaves permanent history intact but removes
    // enough resistance to expose a later increment below alpha's spacing.
    let dwell = law(1.0, 1e3, 1e15);
    let rested = dwell.response(free, &loaded).unwrap().state;
    assert_eq!(rested.native(), loaded.native());
    let short = law(1.0, 1e-17, 1e17);
    let strain = free + 1e6 / E;
    let checked = short
        .advance_with_energy(strain, &rested, strain, tolerances())
        .unwrap();
    assert_eq!(
        checked.state.equivalent_creep_strain(),
        rested.equivalent_creep_strain()
    );
    assert_eq!(checked.state.creep_strain(), rested.creep_strain());
    assert!(checked.state.resistance_strain() > rested.resistance_strain());
    assert!(checked.energy.creep_dissipation_j_m3 > 0.0);
    let delta = 2.0 * checked.state.resistance_strain() - rested.resistance_strain();
    let sigma = checked.stress_pa[0].abs();
    let rate = short
        .kinetics()
        .equivalent_creep_rate_at_strain_s(sigma, checked.state.resistance_strain())
        .unwrap();
    close(delta, 1e-17 * rate, 2e-32);
    close(checked.energy.creep_dissipation_j_m3, sigma * delta, 2e-26);
    for strain in [f64::NAN, f64::INFINITY, 1.0] {
        assert!(short.response(strain, &rested).is_err());
        assert!(short.stress(strain, &rested).is_nan());
    }
}

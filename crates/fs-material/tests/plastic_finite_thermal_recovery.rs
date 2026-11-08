//! Independent finite thermal/recovery references on authored synthetic laws.
//! These tests verify mechanics and thermodynamics, not empirical calibration.

use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::creep::recovery::FiniteRecoveringCreepState;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, AffineThermalRecoveringCreep, FiniteCreepThermalRecovery,
    FiniteJ2TemperatureSlopes, FiniteJ2ThermalError, FiniteJ2ThermalIncrement,
    FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters, FiniteNortonThermalParameters,
    FiniteThermalInelasticMaterial,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const A: f64 = 0.025;
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.1,
    maximum_elastic_log_strain: 0.8,
    maximum_plastic_log_strain: 3.0,
    maximum_accumulated_plastic_strain: 3.0,
};

fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 600.0],
        reference: FiniteJ2Parameters {
            shear_modulus_pa: 80_000.0,
            bulk_modulus_pa: 160_000.0,
            yield_stress_pa: 1000.0,
            hardening_modulus_pa: 12_000.0,
        },
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: -40.0,
            bulk_modulus_pa_k: -60.0,
            yield_stress_pa_k: -0.5,
            hardening_modulus_pa_k: 3.0,
        },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}

fn kinetics(n: f64) -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 400_000.0,
        stress_scale_pa: 2000.0,
        exponent: n,
        activation_temperature_k: 1800.0,
    }
}

fn base_without_primary(combined: bool, n: f64) -> AffineThermalFiniteJ2 {
    let p = parameters();
    if combined {
        AffineThermalFiniteJ2::new(p, LIMITS)
            .unwrap()
            .with_norton_creep(kinetics(n), 3.0)
            .unwrap()
    } else {
        AffineThermalFiniteJ2::new_norton(
            FiniteNortonThermalParameters {
                reference_temperature_k: p.reference_temperature_k,
                temperature_range_k: p.temperature_range_k,
                shear_modulus_pa: p.reference.shear_modulus_pa,
                bulk_modulus_pa: p.reference.bulk_modulus_pa,
                shear_modulus_slope_pa_k: p.slopes.shear_modulus_pa_k,
                bulk_modulus_slope_pa_k: p.slopes.bulk_modulus_pa_k,
                logarithmic_expansion_per_k: p.logarithmic_expansion_per_k,
                logarithmic_expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
            },
            LIMITS,
            kinetics(n),
        )
        .unwrap()
    }
}

fn base(combined: bool, n: f64) -> AffineThermalFiniteJ2 {
    base_without_primary(combined, n)
        .with_norton_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap()
}

fn recovery() -> FiniteCreepThermalRecovery {
    FiniteCreepThermalRecovery {
        reference_recovery_rate_s: 2.0,
        activation_temperature_k: 900.0,
    }
}

fn family(combined: bool, n: f64) -> AffineThermalRecoveringCreep {
    AffineThermalRecoveringCreep::new(base(combined, n), recovery()).unwrap()
}

fn theta(t: f64) -> f64 {
    let p = parameters();
    let delta = t - p.reference_temperature_k;
    (p.logarithmic_expansion_per_k * delta
        + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta)
        .exp()
}

fn coaxial(axial: f64) -> [f64; 9] {
    [
        axial.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * axial).exp(),
    ]
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}",
    );
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

#[test]
fn activated_recovery_uses_the_independent_quadratic_and_original_thermal_potential() {
    let family = family(false, 1.0);
    let (t, dt, q, net) = (420.0, 0.4, 2000.0, 0.025_f64);
    let old = FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [(-2.0 * net).exp(), net.exp(), net.exp(), 0.0, 0.0, 0.0],
            alpha: 0.3,
            creep_alpha: 0.0,
        },
        resistance_strain: 0.04,
    };
    let snapshot = old.clone();
    let mu = 80_000.0 - 40.0 * (t - 300.0);
    let eta = 400_000.0 * (1800.0_f64 * (1.0 / t - 1.0 / 300.0)).exp();
    let k = 2.0 * (900.0_f64 * (1.0 / 300.0 - 1.0 / t)).exp();
    let g = 1.0 + k * dt;
    let b = A * g + old.resistance_strain;
    let forcing = dt * q / eta * A * g;
    let dc = 2.0 * forcing / (b + (b * b + 4.0 * forcing).sqrt());
    let elastic_axial = q / (3.0 * mu);
    let f = coaxial(net + dc + elastic_axial).map(|x| theta(t) * x);
    let frozen = family.at_temperature_for_step(t, dt).unwrap();
    let response = frozen.response(&f, &old).unwrap();
    near(family.recovery_rate_s(t).unwrap(), k, 1e-13, 2e-13);
    near(response.equivalent_kirchhoff_pa, q, 1e-7, 2e-10);
    near(response.plastic_increment, dc, 1e-13, 2e-10);
    near(
        response.state.native.alpha,
        old.native.alpha + dc,
        1e-13,
        2e-10,
    );
    near(
        response.state.resistance_strain,
        (old.resistance_strain + dc) / g,
        1e-13,
        2e-10,
    );
    near(response.inelastic_dissipation_j_m3(), q * dc, 1e-8, 2e-10);
    let flow_potential = eta / dt
        * (0.5 * (1.0 + old.resistance_strain / (A * g)) * dc * dc + dc * dc * dc / (3.0 * A * g));
    near(response.flow_potential_j_m3, flow_potential, 1e-8, 2e-10);
    near(response.stored_energy_j_m3, q * q / (6.0 * mu), 1e-8, 2e-10);

    let thermodynamics = frozen.thermodynamics(&f, &response.state).unwrap();
    let entropy = 40.0 * 1.5 * elastic_axial * elastic_axial;
    near(thermodynamics.entropy_j_m3_k, entropy, 1e-10, 2e-9);
    let increment = family
        .checked_increment_for_step(&coaxial(net), &old, 300.0, &f, t, 1e-6, dt)
        .unwrap();
    near(increment.balance_error_j_m3, 0.0, 1e-8, 0.0);
    near(increment.coefficient_storage_change_j_m3, 0.0, 1e-8, 0.0);
    near(
        increment.thermoelastic_storage_change_j_m3,
        t * entropy,
        1e-8,
        2e-9,
    );
    assert_eq!(increment.mechanical.state, response.state);
    assert_eq!(old, snapshot);
    assert_eq!(frozen.duration_s(), Some(dt));
}

fn same_increment(
    actual: &FiniteJ2ThermalIncrement<FiniteRecoveringCreepState>,
    expected: &FiniteJ2ThermalIncrement,
) {
    assert_eq!(actual.mechanical.state.native, expected.mechanical.state);
    assert_eq!(
        actual.mechanical.first_piola_pa,
        expected.mechanical.first_piola_pa
    );
    assert_eq!(
        actual.mechanical.kirchhoff_pa,
        expected.mechanical.kirchhoff_pa
    );
    assert_eq!(actual.mechanical.tangent_pa, expected.mechanical.tangent_pa);
    assert_eq!(actual.old, expected.old);
    assert_eq!(actual.new, expected.new);
    let fields = |value: &FiniteJ2ThermalIncrement| {
        [
            value.mechanical.incremental_potential_j_m3,
            value.mechanical.flow_potential_j_m3,
            value.mechanical.stored_energy_j_m3,
            value.mechanical.plastic_dissipation_j_m3,
            value.mechanical.plastic_increment,
            value.mechanical.creep_dissipation_j_m3,
            value.mechanical.creep_increment,
            value.mechanical.equivalent_kirchhoff_pa,
            value.mechanical.jacobian,
            value.coefficient_storage_change_j_m3,
            value.thermoelastic_storage_change_j_m3,
            value.mechanical_work_j_m3,
            value.mechanical_coordinate_work_j_m3,
            value.thermal_deformation_work_j_m3,
            value.endpoint_work_defect_j_m3,
            value.reversible_heat_j_m3,
            value.balance_error_j_m3,
        ]
    };
    assert_eq!(
        [
            actual.mechanical.incremental_potential_j_m3,
            actual.mechanical.flow_potential_j_m3,
            actual.mechanical.stored_energy_j_m3,
            actual.mechanical.plastic_dissipation_j_m3,
            actual.mechanical.plastic_increment,
            actual.mechanical.creep_dissipation_j_m3,
            actual.mechanical.creep_increment,
            actual.mechanical.equivalent_kirchhoff_pa,
            actual.mechanical.jacobian,
            actual.coefficient_storage_change_j_m3,
            actual.thermoelastic_storage_change_j_m3,
            actual.mechanical_work_j_m3,
            actual.mechanical_coordinate_work_j_m3,
            actual.thermal_deformation_work_j_m3,
            actual.endpoint_work_defect_j_m3,
            actual.reversible_heat_j_m3,
            actual.balance_error_j_m3,
        ],
        fields(expected)
    );
}

#[test]
fn zero_recovery_preserves_the_complete_original_pure_combined_and_voce_accounts() {
    for base in [
        base(false, 1.7),
        base(true, 1.7),
        base(true, 1.7)
            .with_voce_hardening(800.0, 20.0, 0.3)
            .unwrap(),
    ] {
        let family = AffineThermalRecoveringCreep::new(
            base.clone(),
            FiniteCreepThermalRecovery {
                reference_recovery_rate_s: 0.0,
                activation_temperature_k: 900.0,
            },
        )
        .unwrap();
        let mut old_f = ID;
        let mut old_t = 300.0;
        let mut old = family.initial_state();
        for (t, axial, dt) in [(380.0, 0.06, 0.3), (430.0, 0.10, 0.2)] {
            let f = coaxial(axial).map(|x| theta(t) * x);
            let expected = base
                .checked_increment_for_step(&old_f, &old.native, old_t, &f, t, 1e-6, dt)
                .unwrap();
            let actual = FiniteThermalInelasticMaterial::checked_increment_for_step(
                &family, &old_f, &old, old_t, &f, t, 1e-6, dt,
            )
            .unwrap();
            same_increment(&actual, &expected);
            let accumulation = if base.is_norton_creep() {
                actual.mechanical.state.native.alpha
            } else {
                actual.mechanical.state.native.creep_alpha
            };
            assert_eq!(actual.mechanical.state.resistance_strain, accumulation);
            old = actual.mechanical.state;
            old_f = f;
            old_t = t;
        }
        assert_eq!(
            family.minimum_mechanical_heat_capacity_j_m3_k(),
            base.minimum_mechanical_heat_capacity_j_m3_k()
        );
    }
}

#[test]
fn heated_stress_free_dwell_keeps_the_permanent_metric_and_accumulated_path() {
    let family = family(false, 1.0);
    let permanent_f = [2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];
    let mut old = FiniteRecoveringCreepState {
        native: FiniteJ2State {
            inverse_plastic_cauchy_green: [0.25, 4.0, 1.0, 0.0, 0.0, 0.0],
            alpha: 1.2,
            creep_alpha: 0.0,
        },
        resistance_strain: 0.3,
    };
    let permanent = old.native.clone();
    let mut old_f = permanent_f;
    let mut old_t = 300.0;
    let mut expected_resistance = 0.3;
    for (t, dt) in [(340.0, 0.5), (510.0, 0.2), (300.0, 0.4)] {
        let inspection = family.at_temperature_for_inspection(t).unwrap();
        assert_eq!(inspection.duration_s(), None);
        let f = permanent_f.map(|x| inspection.thermal_stretch() * x);
        assert!(matches!(
            inspection.response(&f, &old),
            Err(FiniteJ2ThermalError::DurationRequired)
        ));
        assert!(matches!(
            inspection.stress_and_potential(&f, &old),
            Err(FiniteJ2ThermalError::DurationRequired)
        ));
        let before = old.clone();
        let next = family
            .checked_increment_for_step(&old_f, &old, old_t, &f, t, 1e-6, dt)
            .unwrap();
        let k = 2.0 * (900.0_f64 * (1.0 / 300.0 - 1.0 / t)).exp();
        expected_resistance /= 1.0 + k * dt;
        near(
            next.mechanical.state.resistance_strain,
            expected_resistance,
            1e-13,
            2e-10,
        );
        near(
            next.mechanical.state.native.alpha,
            permanent.alpha,
            1e-13,
            0.0,
        );
        for (actual, expected) in next
            .mechanical
            .state
            .native
            .inverse_plastic_cauchy_green
            .iter()
            .zip(permanent.inverse_plastic_cauchy_green)
        {
            near(*actual, expected, 1e-12, 2e-10);
        }
        near(next.mechanical.inelastic_dissipation_j_m3(), 0.0, 1e-9, 0.0);
        near(next.caloric_heat_j_m3(), 0.0, 1e-8, 0.0);
        near(next.balance_error_j_m3, 0.0, 1e-8, 0.0);
        assert_eq!(old, before);
        old = next.mechanical.state;
        old_t = t;
        old_f = f;
    }
    assert!(old.resistance_strain < 0.3 && old.native.alpha == 1.2);
}

#[test]
fn total_f_tangent_potential_rotation_and_fixed_history_entropy_remain_consistent() {
    for combined in [false, true] {
        let family = family(combined, 2.5);
        let old = family
            .at_temperature_for_step(340.0, 0.1)
            .unwrap()
            .response(
                &coaxial(0.025).map(|x| theta(340.0) * x),
                &family.initial_state(),
            )
            .unwrap()
            .state;
        let old_snapshot = old.clone();
        let frozen = family.at_temperature_for_step(430.0, 0.2).unwrap();
        let f = [1.11, 0.04, -0.01, 0.015, 0.97, 0.02, 0.0, -0.01, 0.96];
        let response = frozen.response(&f, &old).unwrap();
        let (stress, potential) = frozen.stress_and_potential(&f, &old).unwrap();
        assert_eq!(stress, response.first_piola_pa);
        assert_eq!(potential, response.incremental_potential_j_m3);
        let h = 2e-6;
        for column in 0..9 {
            let mut plus = f;
            let mut minus = f;
            plus[column] += h;
            minus[column] -= h;
            let (sp, wp) = frozen.stress_and_potential(&plus, &old).unwrap();
            let (sm, wm) = frozen.stress_and_potential(&minus, &old).unwrap();
            near((wp - wm) / (2.0 * h), stress[column], 2e-4, 2e-6);
            for row in 0..9 {
                near(
                    (sp[row] - sm[row]) / (2.0 * h),
                    response.tangent_pa[row][column],
                    0.03,
                    3e-5,
                );
                near(
                    response.tangent_pa[row][column],
                    response.tangent_pa[column][row],
                    1e-7,
                    2e-10,
                );
            }
        }
        let angle = 0.7_f64;
        let r = [
            angle.cos(),
            -angle.sin(),
            0.0,
            angle.sin(),
            angle.cos(),
            0.0,
            0.0,
            0.0,
            1.0,
        ];
        let rotated = frozen.response(&multiply(&r, &f), &old).unwrap();
        let expected_piola = multiply(&r, &response.first_piola_pa);
        for (actual, expected) in rotated.first_piola_pa.iter().zip(expected_piola) {
            near(*actual, expected, 1e-7, 2e-9);
        }
        near(rotated.incremental_potential_j_m3, potential, 1e-8, 2e-10);
        near(
            rotated.state.resistance_strain,
            response.state.resistance_strain,
            1e-12,
            2e-9,
        );

        // Different resistance at fixed permanent history changes kinetics only.
        let mut recovered = response.state.clone();
        recovered.resistance_strain = 0.0;
        let thermo = frozen.thermodynamics(&f, &response.state).unwrap();
        assert_eq!(thermo, frozen.thermodynamics(&f, &recovered).unwrap());
        let h_t = 0.002;
        let upper = family
            .thermodynamics(&f, &response.state, 430.0 + h_t)
            .unwrap();
        let lower = family
            .thermodynamics(&f, &response.state, 430.0 - h_t)
            .unwrap();
        near(
            -(upper.free_energy_j_m3 - lower.free_energy_j_m3) / (2.0 * h_t),
            thermo.entropy_j_m3_k,
            1e-5,
            2e-7,
        );
        near(
            thermo.internal_energy_j_m3,
            thermo.free_energy_j_m3 + 430.0 * thermo.entropy_j_m3_k,
            1e-10,
            2e-13,
        );
        assert_eq!(old, old_snapshot);
    }
}

#[test]
fn invalid_selection_activation_duration_and_complete_history_refuse_without_reset() {
    for combined in [false, true] {
        assert!(
            AffineThermalRecoveringCreep::new(base_without_primary(combined, 1.0), recovery())
                .is_err()
        );
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            AffineThermalRecoveringCreep::new(
                base(false, 1.0),
                FiniteCreepThermalRecovery {
                    reference_recovery_rate_s: bad,
                    ..recovery()
                }
            )
            .is_err()
        );
        assert!(
            AffineThermalRecoveringCreep::new(
                base(false, 1.0),
                FiniteCreepThermalRecovery {
                    activation_temperature_k: bad,
                    ..recovery()
                }
            )
            .is_err()
        );
    }
    assert!(
        AffineThermalRecoveringCreep::new(
            base(false, 1.0),
            FiniteCreepThermalRecovery {
                activation_temperature_k: 1e7,
                ..recovery()
            }
        )
        .is_err()
    );
    let disabled = AffineThermalRecoveringCreep::new(
        base(false, 1.0),
        FiniteCreepThermalRecovery {
            reference_recovery_rate_s: 0.0,
            activation_temperature_k: 1e7,
        },
    )
    .unwrap();
    assert_eq!(disabled.recovery_rate_s(600.0).unwrap(), 0.0);

    let family = family(true, 1.0);
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(family.at_temperature_for_step(300.0, dt).is_err());
    }
    for t in [0.0, 279.0, 601.0, f64::NAN, f64::INFINITY] {
        assert!(family.at_temperature_for_inspection(t).is_err());
        assert!(family.recovery_rate_s(t).is_err());
    }
    let mut bad = family.initial_state();
    bad.native.alpha = 0.2;
    bad.resistance_strain = 0.01;
    let before = bad.clone();
    // Plastic accumulation cannot license resistance with no permanent creep.
    assert!(
        family
            .at_temperature_for_inspection(300.0)
            .unwrap()
            .thermodynamics(&ID, &bad)
            .is_err()
    );
    assert!(
        family
            .checked_increment_for_step(&ID, &bad, 300.0, &ID, 350.0, 1e-6, 0.1)
            .is_err()
    );
    assert_eq!(bad, before);
}

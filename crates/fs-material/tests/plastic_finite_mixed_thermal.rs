//! Independent total-volume, thermal and kinetic references for mixed finite flow.
//! Authored synthetic coefficients verify the equations, not material calibration.

use fs_material::creep::StrainHardeningCreep;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::plastic::finite::creep::recovery::FiniteRecoveringCreepState;
use fs_material::plastic::finite::thermal::mixed::{
    MixedThermalInelasticPoint, checked_mixed_increment_for_step,
};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, AffineThermalRecoveringCreep, FiniteCreepThermalRecovery,
    FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters,
    FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const F: [f64; 9] = [1.18, 0.15, 0.03, -0.04, 0.94, 0.07, 0.02, -0.03, 1.0];
const DT: f64 = 0.25;
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.1,
    maximum_elastic_log_strain: 0.8,
    maximum_plastic_log_strain: 3.0,
    maximum_accumulated_plastic_strain: 3.0,
};

fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
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

fn kinetics() -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 60_000.0,
        stress_scale_pa: 2000.0,
        exponent: 1.4,
        activation_temperature_k: 1800.0,
    }
}

fn norton(kinetics: FiniteJ2ThermalOverstress) -> AffineThermalFiniteJ2 {
    let p = parameters();
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
        kinetics,
    )
    .unwrap()
}

fn families() -> [AffineThermalFiniteJ2; 4] {
    let base = AffineThermalFiniteJ2::new(parameters(), LIMITS).unwrap();
    [
        base.clone(),
        base.clone().with_overstress(kinetics()).unwrap(),
        norton(kinetics()),
        base.with_norton_creep(kinetics(), 3.0).unwrap(),
    ]
}

fn theta(t: f64) -> f64 {
    let p = parameters();
    let delta = t - p.reference_temperature_k;
    (p.logarithmic_expansion_per_k * delta
        + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta)
        .exp()
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn pressure(f: &[f64; 9], t: f64) -> f64 {
    let bulk = 160_000.0 - 60.0 * (t - 300.0);
    bulk * (determinant(f) / theta(t).powi(3)).ln() / determinant(f)
}

fn coaxial(axial: f64, mean: f64) -> [f64; 9] {
    [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ]
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}, difference {:.8e}",
        (actual - expected).abs()
    );
}

#[test]
fn frozen_history_total_f_potential_and_all_pressure_tangent_blocks_have_exact_derivatives() {
    for family in families() {
        let old = family
            .at_temperature_for_step(300.0, DT)
            .unwrap()
            .response(&coaxial(0.04, 0.0), &family.initial_state())
            .unwrap()
            .state;
        let snapshot = old.clone();
        let point = family.at_temperature_for_step(410.0, DT).unwrap();
        let frozen = point.freeze_mixed(&old).unwrap();
        let f = F.map(|x| theta(410.0) * x);
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let p = 3000.0;
            let center = frozen.mixed_response(&f, p, treatment).unwrap();
            for col in 0..9 {
                let h = 1e-6;
                let (mut plus, mut minus) = (f, f);
                plus[col] += h;
                minus[col] -= h;
                let plus = frozen.mixed_state(&plus, p, treatment).unwrap();
                let minus = frozen.mixed_state(&minus, p, treatment).unwrap();
                near(
                    (plus.potential_pa - minus.potential_pa) / (2.0 * h),
                    center.state.first_piola_pa[col],
                    5e-3,
                    3e-7,
                );
                near(
                    (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
                    center.pressure_coupling[col],
                    3e-9,
                    3e-9,
                );
                for row in 0..9 {
                    near(
                        (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                        center.tangent_pa[row][col],
                        0.03,
                        3e-6,
                    );
                    near(
                        center.tangent_pa[row][col],
                        center.tangent_pa[col][row],
                        2e-7,
                        3e-11,
                    );
                }
            }
            let h = 1.0;
            let plus = frozen.mixed_state(&f, p + h, treatment).unwrap();
            let minus = frozen.mixed_state(&f, p - h, treatment).unwrap();
            near(
                (plus.potential_pa - minus.potential_pa) / (2.0 * h),
                center.state.pressure_residual,
                3e-9,
                3e-8,
            );
            near(
                (plus.pressure_residual - minus.pressure_residual) / (2.0 * h),
                center.pressure_tangent_per_pa,
                3e-14,
                3e-8,
            );
            for row in 0..9 {
                near(
                    (plus.first_piola_pa[row] - minus.first_piola_pa[row]) / (2.0 * h),
                    center.pressure_coupling[row],
                    2e-10,
                    2e-10,
                );
            }
        }
        assert_eq!(old, snapshot);
    }
}

#[test]
fn local_pressure_condensation_recovers_native_thermal_endpoints_and_schur_tangent() {
    let treatment = VolumetricTreatment::OriginalMaterial;
    let old_f = coaxial(0.06, 0.015);
    let (old_t, new_t) = (300.0, 410.0);
    let new_f = F.map(|x| theta(new_t) * x);
    for family in families() {
        let old = family
            .at_temperature_for_step(old_t, DT)
            .unwrap()
            .response(&old_f, &family.initial_state())
            .unwrap()
            .state;
        let native = family
            .checked_increment_for_step(&old_f, &old, old_t, &new_f, new_t, 1e-6, DT)
            .unwrap();
        let point = family.at_temperature_for_step(new_t, DT).unwrap();
        let frozen = point.freeze_mixed(&old).unwrap();
        let mixed = frozen
            .mixed_response(&new_f, pressure(&new_f, new_t), treatment)
            .unwrap();
        let step = checked_mixed_increment_for_step(
            &family,
            &old_f,
            &old,
            old_t,
            pressure(&old_f, old_t),
            &new_f,
            new_t,
            pressure(&new_f, new_t),
            treatment,
            1e-6,
            DT,
        )
        .unwrap();
        assert_eq!(step.mechanical.state, native.mechanical.state);
        assert_eq!(
            step.mechanical.flow_potential_j_m3,
            native.mechanical.flow_potential_j_m3
        );
        assert_eq!(
            step.mechanical.inelastic_dissipation_j_m3(),
            native.mechanical.inelastic_dissipation_j_m3()
        );
        near(mixed.state.pressure_residual, 0.0, 3e-14, 0.0);
        near(
            mixed.state.potential_pa,
            native.mechanical.incremental_potential_j_m3,
            3e-8,
            3e-11,
        );
        near(
            step.mechanical.stored_energy_j_m3 + mixed.state.bulk_stored_pa,
            native.mechanical.stored_energy_j_m3,
            3e-8,
            3e-11,
        );
        for (a, b) in [
            (step.old.free_energy_j_m3, native.old.free_energy_j_m3),
            (step.new.free_energy_j_m3, native.new.free_energy_j_m3),
            (step.old.entropy_j_m3_k, native.old.entropy_j_m3_k),
            (step.new.entropy_j_m3_k, native.new.entropy_j_m3_k),
            (
                step.old.internal_energy_j_m3,
                native.old.internal_energy_j_m3,
            ),
            (
                step.new.internal_energy_j_m3,
                native.new.internal_energy_j_m3,
            ),
        ] {
            near(a, b, 3e-8, 3e-11);
        }
        for row in 0..9 {
            near(
                mixed.state.first_piola_pa[row],
                native.mechanical.first_piola_pa[row],
                3e-8,
                3e-11,
            );
            near(
                step.new.first_piola_pa[row],
                native.new.first_piola_pa[row],
                3e-8,
                3e-11,
            );
            for col in 0..9 {
                let condensed = mixed.tangent_pa[row][col]
                    - mixed.pressure_coupling[row] * mixed.pressure_coupling[col]
                        / mixed.pressure_tangent_per_pa;
                near(
                    condensed,
                    native.mechanical.tangent_pa[row][col],
                    3e-5,
                    3e-10,
                );
            }
        }
        let endpoint_work: f64 = (0..9)
            .map(|i| native.new.first_piola_pa[i] * (new_f[i] - old_f[i]))
            .sum();
        near(step.mechanical_work_j_m3, endpoint_work, 3e-8, 3e-11);
        // Mixed bulk work uses reconstructed scalar volume; its finite-step
        // thermal/work-defect split need not equal native matrix-Fm quadrature.
        near(
            endpoint_work,
            step.internal_energy_change_j_m3()
                + step.caloric_heat_j_m3()
                + step.endpoint_work_defect_j_m3,
            3e-8,
            3e-11,
        );
        near(step.balance_error_j_m3, 0.0, 3e-8, 0.0);
    }
}

#[test]
fn exact_thermal_volume_has_pressure_entropy_and_scalar_work_without_penalty_storage() {
    let treatment = VolumetricTreatment::ExactIncompressible;
    let t = 410.0;
    let stretch = theta(t);
    let volume = stretch.powi(3);
    let expansion = 2e-4 + 1e-7 * (t - 300.0);
    let f = ID.map(|x| stretch * x);
    // The exact path must never form unused K*log(Jm)^2 from rounded total F.
    for bulk in [160_000.0, 1e300] {
        let mut parameters = parameters();
        parameters.reference.bulk_modulus_pa = bulk;
        parameters.slopes.bulk_modulus_pa_k = 0.0;
        let family = AffineThermalFiniteJ2::new(parameters, LIMITS).unwrap();
        let old = family.initial_state();
        let point = family.at_temperature_for_step(t, DT).unwrap();
        let frozen = point.freeze_mixed(&old).unwrap();
        for p in [-4000.0, 5000.0] {
            let response = frozen.mixed_response(&f, p, treatment).unwrap();
            let thermodynamics = point.mixed_thermodynamics(&f, &old, p, treatment).unwrap();
            let entropy = 3.0 * expansion * p * volume;
            near(response.state.pressure_residual, 0.0, 3e-14, 0.0);
            near(response.state.bulk_volume_ratio, volume, 3e-14, 0.0);
            near(
                response.state.bulk_complementary_pa,
                p * (volume - 1.0),
                3e-10,
                3e-12,
            );
            assert_eq!(response.state.bulk_stored_pa, 0.0);
            assert_eq!(response.pressure_tangent_per_pa, 0.0);
            assert_eq!(thermodynamics.bulk_free_energy_j_m3, 0.0);
            near(thermodynamics.free_energy_j_m3, 0.0, 1e-18, 0.0);
            near(thermodynamics.bulk_entropy_j_m3_k, entropy, 3e-11, 3e-12);
            near(thermodynamics.entropy_j_m3_k, entropy, 3e-11, 3e-12);
            near(
                thermodynamics.internal_energy_j_m3,
                t * entropy,
                3e-8,
                3e-12,
            );
            for i in 0..9 {
                near(
                    response.state.first_piola_pa[i],
                    p * stretch * stretch * ID[i],
                    3e-8,
                    3e-12,
                );
            }

            // At fixed total F and pressure, entropy differentiates the
            // constraint term even though its physical bulk storage is zero.
            let h = 0.02;
            let plus_point = family.at_temperature_for_step(t + h, DT).unwrap();
            let minus_point = family.at_temperature_for_step(t - h, DT).unwrap();
            let plus = plus_point
                .freeze_mixed(&old)
                .unwrap()
                .mixed_state(&f, p, treatment)
                .unwrap();
            let minus = minus_point
                .freeze_mixed(&old)
                .unwrap()
                .mixed_state(&f, p, treatment)
                .unwrap();
            near(
                -(plus.potential_pa - minus.potential_pa) / (2.0 * h),
                entropy,
                2e-7,
                2e-8,
            );

            let step = checked_mixed_increment_for_step(
                &family, &ID, &old, 300.0, p, &f, t, p, treatment, 1e-6, DT,
            )
            .unwrap();
            let work = 3.0 * p * stretch * stretch * (stretch - 1.0);
            let thermal_work = p * (volume - 1.0);
            let internal_change = t * entropy - 300.0 * (3.0 * 2e-4 * p);
            near(step.mechanical_work_j_m3, work, 3e-8, 3e-12);
            near(
                step.thermal_deformation_work_j_m3,
                thermal_work,
                3e-8,
                3e-12,
            );
            near(
                step.mechanical_coordinate_work_j_m3,
                work - thermal_work,
                3e-8,
                3e-12,
            );
            near(
                step.endpoint_work_defect_j_m3,
                work - thermal_work,
                3e-8,
                3e-12,
            );
            near(step.coefficient_storage_change_j_m3, 0.0, 1e-18, 0.0);
            near(
                step.internal_energy_change_j_m3(),
                internal_change,
                3e-8,
                3e-12,
            );
            near(
                step.reversible_heat_j_m3,
                thermal_work - internal_change,
                3e-8,
                3e-12,
            );
            near(
                step.caloric_heat_j_m3(),
                thermal_work - internal_change,
                3e-8,
                3e-12,
            );
            near(step.balance_error_j_m3, 0.0, 3e-8, 0.0);
            assert_eq!(step.mechanical.inelastic_dissipation_j_m3(), 0.0);
            assert_eq!(step.mechanical.state, old);
        }
    }
}

#[test]
fn thermal_softening_checks_old_history_at_old_temperature_and_preserves_first_law() {
    let mut parameters = parameters();
    parameters.slopes = FiniteJ2TemperatureSlopes {
        shear_modulus_pa_k: 0.0,
        bulk_modulus_pa_k: 0.0,
        yield_stress_pa_k: -4.0,
        hardening_modulus_pa_k: -10.0,
    };
    parameters.logarithmic_expansion_per_k = 0.0;
    parameters.logarithmic_expansion_slope_per_k2 = 0.0;
    let family = AffineThermalFiniteJ2::new(parameters, LIMITS).unwrap();
    let treatment = VolumetricTreatment::ExactIncompressible;
    let axial = 0.15;
    let f = coaxial(axial, 0.0);
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&f, &family.initial_state())
        .unwrap()
        .state;
    let snapshot = old.clone();
    let hot = family.at_temperature_for_step(420.0, DT).unwrap();
    assert!(
        hot.validate_mixed_committed(&f, &old, 0.0, treatment, 1e-6)
            .is_err()
    );
    let step = checked_mixed_increment_for_step(
        &family, &f, &old, 300.0, 0.0, &f, 420.0, 0.0, treatment, 1e-6, DT,
    )
    .unwrap();
    // Independent radial return at fixed total deformation after Y and H soften.
    let alpha0 = (3.0 * 80_000.0 * axial - 1000.0) / (3.0 * 80_000.0 + 12_000.0);
    let (yield_hot, hardening_hot) = (520.0, 10_800.0);
    let increment = (3.0 * 80_000.0 * (axial - alpha0) - yield_hot - hardening_hot * alpha0)
        / (3.0 * 80_000.0 + hardening_hot);
    let alpha1 = alpha0 + increment;
    let psi0 = 1.5 * 80_000.0 * (axial - alpha0).powi(2) + 0.5 * 12_000.0 * alpha0.powi(2);
    let psi1 = 1.5 * 80_000.0 * (axial - alpha1).powi(2) + 0.5 * hardening_hot * alpha1.powi(2);
    let (s0, s1) = (5.0 * alpha0.powi(2), 5.0 * alpha1.powi(2));
    let coefficient_change = 0.5 * (hardening_hot - 12_000.0) * alpha0.powi(2);
    let entropy_storage_change = 420.0 * s1 - 300.0 * s0;
    let physical_heat = yield_hot * increment;
    let defect = 0.5 * (3.0 * 80_000.0 + hardening_hot) * increment.powi(2);
    assert!(increment > 0.0);
    for (actual, expected) in [
        (old.alpha, alpha0),
        (step.mechanical.plastic_increment, increment),
        (step.mechanical.state.alpha, alpha1),
        (step.old.free_energy_j_m3, psi0),
        (step.new.free_energy_j_m3, psi1),
        (step.old.entropy_j_m3_k, s0),
        (step.new.entropy_j_m3_k, s1),
        (step.coefficient_storage_change_j_m3, coefficient_change),
        (
            step.thermoelastic_storage_change_j_m3,
            entropy_storage_change,
        ),
        (step.mechanical.inelastic_dissipation_j_m3(), physical_heat),
        (step.endpoint_work_defect_j_m3, defect),
        (
            step.reversible_heat_j_m3,
            -coefficient_change - entropy_storage_change,
        ),
        (
            step.internal_energy_change_j_m3(),
            psi1 - psi0 + entropy_storage_change,
        ),
    ] {
        near(actual, expected, 2e-8, 3e-10);
    }
    assert_eq!(step.mechanical_work_j_m3, 0.0);
    near(
        step.caloric_heat_j_m3(),
        physical_heat - coefficient_change - entropy_storage_change,
        2e-8,
        3e-10,
    );
    near(
        step.internal_energy_change_j_m3() + step.caloric_heat_j_m3() + defect,
        0.0,
        2e-8,
        0.0,
    );
    assert!(
        checked_mixed_increment_for_step(
            &family, &f, &old, 420.0, 0.0, &f, 420.0, 0.0, treatment, 1e-6, DT,
        )
        .is_err()
    );

    // A coefficient change at stress-free retained plastic deformation stores
    // -Delta(H alpha^2/2) in T*s. It creates neither flow nor caloric heat.
    let relaxed = coaxial(alpha0, 0.0);
    let held = checked_mixed_increment_for_step(
        &family, &relaxed, &old, 300.0, 0.0, &relaxed, 420.0, 0.0, treatment, 1e-6, DT,
    )
    .unwrap();
    assert_eq!(held.mechanical.state, old);
    near(
        held.coefficient_storage_change_j_m3,
        coefficient_change,
        2e-8,
        3e-10,
    );
    near(
        held.thermoelastic_storage_change_j_m3,
        -coefficient_change,
        2e-8,
        3e-10,
    );
    near(held.caloric_heat_j_m3(), 0.0, 2e-8, 0.0);
    near(held.internal_energy_change_j_m3(), 0.0, 2e-8, 0.0);
    assert_eq!(old, snapshot);
}

#[test]
fn activated_primary_recovery_keeps_its_clock_complete_memory_and_physical_heat() {
    let primary_strain = 0.025;
    let rate = FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 400_000.0,
        exponent: 1.0,
        ..kinetics()
    };
    let recovery = FiniteCreepThermalRecovery {
        reference_recovery_rate_s: 2.0,
        activation_temperature_k: 900.0,
    };
    let (t, dt, stress, net) = (420.0, 0.4, 2000.0, 0.025_f64);
    let mu = 80_000.0 - 40.0 * (t - 300.0);
    let hardening = 12_000.0 + 3.0 * (t - 300.0);
    let yield_stress = 1000.0 - 0.5 * (t - 300.0);
    let viscosity = 400_000.0 * (1800.0_f64 * (1.0 / t - 1.0 / 300.0)).exp();
    let recovery_rate = 2.0 * (900.0_f64 * (1.0 / 300.0 - 1.0 / t)).exp();
    let relaxation = 1.0 + recovery_rate * dt;
    let treatment = VolumetricTreatment::ExactIncompressible;
    for combined in [false, true] {
        let base = if combined {
            AffineThermalFiniteJ2::new(parameters(), LIMITS)
                .unwrap()
                .with_norton_creep(rate, 3.0)
                .unwrap()
        } else {
            norton(rate)
        };
        let base = base
            .with_norton_strain_hardening(StrainHardeningCreep {
                reference_strain: primary_strain,
                exponent: 1.0,
                residual_rate_fraction: 0.0,
            })
            .unwrap();
        let family = AffineThermalRecoveringCreep::new(base, recovery).unwrap();
        let old = FiniteRecoveringCreepState {
            native: FiniteJ2State {
                inverse_plastic_cauchy_green: [
                    (-2.0 * net).exp(),
                    net.exp(),
                    net.exp(),
                    0.0,
                    0.0,
                    0.0,
                ],
                alpha: if combined { 0.03 } else { 0.3 },
                creep_alpha: if combined { 0.3 } else { 0.0 },
            },
            resistance_strain: 0.04,
        };
        let snapshot = old.clone();
        let b = primary_strain * relaxation + old.resistance_strain;
        let forcing = dt * stress / viscosity * primary_strain * relaxation;
        let dc = 2.0 * forcing / (b + (b * b + 4.0 * forcing).sqrt());
        let dp = if combined {
            (stress - yield_stress - hardening * old.native.alpha) / hardening
        } else {
            0.0
        };
        let elastic_axial = stress / (3.0 * mu);
        let old_f = coaxial(net, 0.0);
        let new_f = coaxial(net + dc + dp + elastic_axial, 0.0).map(|x| theta(t) * x);
        let p = 5000.0;
        let point = family.at_temperature_for_step(t, dt).unwrap();
        assert_eq!(point.duration_s(), Some(dt));
        let step = checked_mixed_increment_for_step(
            &family, &old_f, &old, 300.0, 0.0, &new_f, t, p, treatment, 1e-6, dt,
        )
        .unwrap();
        let new_state = &step.mechanical.state;
        near(step.mechanical.equivalent_kirchhoff_pa, stress, 2e-7, 3e-10);
        near(
            new_state.resistance_strain,
            (old.resistance_strain + dc) / relaxation,
            3e-13,
            3e-10,
        );
        near(
            new_state.native.alpha,
            old.native.alpha + if combined { dp } else { dc },
            3e-13,
            3e-10,
        );
        near(
            new_state.native.creep_alpha,
            old.native.creep_alpha + if combined { dc } else { 0.0 },
            3e-13,
            3e-10,
        );
        near(
            new_state.native.inverse_plastic_cauchy_green[0],
            (-2.0 * (net + dp + dc)).exp(),
            3e-13,
            3e-10,
        );
        let physical_heat = stress * dc + yield_stress * dp;
        let flow = yield_stress * dp
            + viscosity / dt
                * (0.5 * (1.0 + old.resistance_strain / (primary_strain * relaxation)) * dc * dc
                    + dc * dc * dc / (3.0 * primary_strain * relaxation));
        near(
            step.mechanical.inelastic_dissipation_j_m3(),
            physical_heat,
            2e-7,
            3e-10,
        );
        near(step.mechanical.flow_potential_j_m3, flow, 2e-7, 3e-10);
        assert!(physical_heat > flow);
        let alpha0 = if combined { old.native.alpha } else { 0.0 };
        let alpha1 = alpha0 + dp;
        let psi0 = 0.5 * 12_000.0 * alpha0 * alpha0;
        let psi1 = stress * stress / (6.0 * mu) + 0.5 * hardening * alpha1 * alpha1;
        let s0 = -1.5 * alpha0 * alpha0;
        let s1 = 60.0 * elastic_axial * elastic_axial - 1.5 * alpha1 * alpha1
            + 3.0 * (2e-4 + 1e-7 * (t - 300.0)) * p * theta(t).powi(3);
        near(step.old.free_energy_j_m3, psi0, 2e-7, 3e-10);
        near(step.new.free_energy_j_m3, psi1, 2e-7, 3e-10);
        near(step.old.entropy_j_m3_k, s0, 2e-9, 3e-10);
        near(step.new.entropy_j_m3_k, s1, 2e-9, 3e-10);
        let coefficient_change = 0.5 * (hardening - 12_000.0) * alpha0 * alpha0;
        let entropy_storage_change = t * s1 - 300.0 * s0;
        let thermal_ratio = theta(t) - 1.0;
        let remainder_work: f64 = [
            (0, 2.0 * stress / 3.0),
            (4, -stress / 3.0),
            (8, -stress / 3.0),
        ]
        .into_iter()
        .map(|(i, tau)| tau / new_f[i] * (new_f[i] - old_f[i]))
        .sum();
        let remainder_thermal_work: f64 = [
            (0, 2.0 * stress / 3.0),
            (4, -stress / 3.0),
            (8, -stress / 3.0),
        ]
        .into_iter()
        .map(|(i, tau)| tau / new_f[i] * old_f[i] * thermal_ratio)
        .sum();
        let bulk_work: f64 = [0, 4, 8]
            .into_iter()
            .map(|i| p * theta(t).powi(3) / new_f[i] * (new_f[i] - old_f[i]))
            .sum();
        let thermal_work = remainder_thermal_work + p * (theta(t).powi(3) - 1.0);
        let total_work = remainder_work + bulk_work;
        let defect = total_work - thermal_work - (psi1 - psi0 - coefficient_change) - physical_heat;
        near(step.mechanical_work_j_m3, total_work, 2e-7, 3e-10);
        near(
            step.thermal_deformation_work_j_m3,
            thermal_work,
            2e-7,
            3e-10,
        );
        near(
            step.coefficient_storage_change_j_m3,
            coefficient_change,
            2e-7,
            3e-10,
        );
        near(step.endpoint_work_defect_j_m3, defect, 2e-7, 3e-10);
        near(
            step.reversible_heat_j_m3,
            thermal_work - coefficient_change - entropy_storage_change,
            2e-7,
            3e-10,
        );
        near(
            step.internal_energy_change_j_m3(),
            psi1 - psi0 + entropy_storage_change,
            2e-7,
            3e-10,
        );
        near(
            total_work,
            step.internal_energy_change_j_m3() + step.caloric_heat_j_m3() + defect,
            2e-7,
            3e-10,
        );
        near(step.balance_error_j_m3, 0.0, 2e-7, 0.0);

        let inspection = family.at_temperature_for_inspection(t).unwrap();
        let missing_interval_refused = match inspection.freeze_mixed(&old) {
            Ok(frozen) => frozen.mixed_response(&new_f, p, treatment).is_err(),
            Err(_) => true,
        };
        assert!(missing_interval_refused);
        for duration in [0.0, -dt, f64::NAN] {
            assert!(
                checked_mixed_increment_for_step(
                    &family, &old_f, &old, 300.0, 0.0, &new_f, t, p, treatment, 1e-6, duration,
                )
                .is_err()
            );
        }
        assert_eq!(old, snapshot);
        let retry = checked_mixed_increment_for_step(
            &family, &old_f, &old, 300.0, 0.0, &new_f, t, p, treatment, 1e-6, dt,
        )
        .unwrap();
        assert_eq!(retry.mechanical.state, step.mechanical.state);
        assert_eq!(
            retry.mechanical.flow_potential_j_m3,
            step.mechanical.flow_potential_j_m3
        );
    }
}

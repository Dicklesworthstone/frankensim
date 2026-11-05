//! G1/G3 finite Voce thermal energy, endpoint freezing and kinetic composition.
//! Synthetic coefficients check the constitutive model, not material calibration.

use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters, FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const Q: f64 = 12e6;
const B: f64 = 9.0;
const DQ: f64 = -20_000.0;

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        reference: FiniteJ2Parameters {
            shear_modulus_pa: 80e6,
            bulk_modulus_pa: 160e6,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 8e6,
        },
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: -20_000.0,
            bulk_modulus_pa_k: -30_000.0,
            yield_stress_pa_k: -4_000.0,
            hardening_modulus_pa_k: -10_000.0,
        },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}

fn constant() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: 0.0,
            bulk_modulus_pa_k: 0.0,
            yield_stress_pa_k: 0.0,
            hardening_modulus_pa_k: 0.0,
        },
        logarithmic_expansion_per_k: 0.0,
        logarithmic_expansion_slope_per_k2: 0.0,
        ..parameters()
    }
}

fn coaxial(e: f64) -> [f64; 9] {
    [
        e.exp(), 0.0, 0.0, 0.0, (-0.5 * e).exp(), 0.0, 0.0, 0.0, (-0.5 * e).exp(),
    ]
}

fn shape(alpha: f64) -> f64 {
    alpha + (-B * alpha).exp_m1() / B
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}

fn kinetics() -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 40e6,
        stress_scale_pa: 3e6,
        exponent: 2.0,
        activation_temperature_k: 600.0,
    }
}

#[test]
fn constant_saturation_recovers_native_response_and_zero_saturation_keeps_linear_arithmetic() {
    let base = AffineThermalFiniteJ2::new(constant(), limits()).unwrap();
    let family = base.clone().with_voce_hardening(Q, B, 0.0).unwrap();
    let native = FiniteJ2::new(constant().reference, limits())
        .unwrap()
        .with_voce_hardening(Q, B)
        .unwrap();
    let f = [1.2, 0.12, 0.01, 0.0, 0.94, 0.04, 0.02, 0.0, 0.95];
    let old = family.initial_state();
    let expected = native.response(&f, &old).unwrap();
    assert!(expected.plastic_increment > 0.0);
    let step = family
        .checked_increment(&ID, &old, 300.0, &f, 410.0, 1e-3)
        .unwrap();
    assert_eq!(step.mechanical, expected);
    assert_eq!(step.coefficient_storage_change_j_m3, 0.0);
    assert_eq!(step.thermoelastic_storage_change_j_m3, 0.0);
    assert_eq!(step.reversible_heat_j_m3, 0.0);
    assert_eq!(step.caloric_heat_j_m3(), expected.plastic_dissipation_j_m3);
    near(step.balance_error_j_m3, 0.0, 2e-7, 0.0);

    let zero = base.clone().with_voce_hardening(0.0, B, 0.0).unwrap();
    assert_eq!(
        zero.checked_increment(&ID, &old, 300.0, &f, 410.0, 1e-3).unwrap(),
        base.checked_increment(&ID, &old, 300.0, &f, 410.0, 1e-3).unwrap()
    );
}

#[test]
fn affine_saturation_entropy_and_internal_energy_are_fixed_history_derivatives() {
    let base = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    let family = base.clone().with_voce_hardening(Q, B, DQ).unwrap();
    let history = family
        .at_temperature(300.0).unwrap()
        .response(&coaxial(0.15), &family.initial_state()).unwrap().state;
    let f = [1.12, 0.1, 0.02, -0.03, 0.95, 0.05, 0.01, 0.0, 1.01];
    let t = 370.0;
    let baseline = base.thermodynamics(&f, &history, t).unwrap();
    let center = family.thermodynamics(&f, &history, t).unwrap();
    let primitive = shape(history.alpha);
    let q = Q + DQ * (t - 300.0);
    near(center.free_energy_j_m3 - baseline.free_energy_j_m3, q * primitive, 1e-8, 2e-13);
    near(center.entropy_j_m3_k - baseline.entropy_j_m3_k, -DQ * primitive, 1e-10, 2e-13);
    near(center.internal_energy_j_m3 - baseline.internal_energy_j_m3, (q - t * DQ) * primitive, 2e-8, 2e-13);
    assert_eq!(center.heat_capacity_j_m3_k, baseline.heat_capacity_j_m3_k);
    assert_eq!(center.first_piola_pa, baseline.first_piola_pa);
    assert_eq!(
        center.first_piola_temperature_derivative_fixed_history_pa_k,
        baseline.first_piola_temperature_derivative_fixed_history_pa_k
    );
    assert_eq!(family.minimum_mechanical_heat_capacity_j_m3_k(), base.minimum_mechanical_heat_capacity_j_m3_k());
    let h = 0.05;
    let plus = family.thermodynamics(&f, &history, t + h).unwrap();
    let minus = family.thermodynamics(&f, &history, t - h).unwrap();
    near(-(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * h), center.entropy_j_m3_k, 0.01, 2e-7);
    near(-t * (plus.free_energy_j_m3 - 2.0 * center.free_energy_j_m3 + minus.free_energy_j_m3) / (h * h), center.heat_capacity_j_m3_k, 0.05, 2e-5);
}

#[test]
fn zero_strength_endpoint_retains_entropy_and_closed_frozen_history_thermal_cycle() {
    let base = AffineThermalFiniteJ2::new(constant(), limits()).unwrap();
    let q0 = 6e6;
    let slope = -30_000.0;
    let family = base.clone().with_voce_hardening(q0, B, slope).unwrap();
    let history = family
        .at_temperature(300.0).unwrap()
        .response(&coaxial(0.18), &family.initial_state()).unwrap().state;
    // Removing the returned elastic strain leaves the retained plastic metric
    // stress free, so changing coefficients cannot trigger another return.
    let f = coaxial(history.alpha);
    let at_zero_q = family.thermodynamics(&f, &history, 500.0).unwrap();
    let linear = base.thermodynamics(&f, &history, 500.0).unwrap();
    assert_eq!(at_zero_q.free_energy_j_m3, linear.free_energy_j_m3);
    near(at_zero_q.entropy_j_m3_k, -slope * shape(history.alpha), 1e-10, 2e-13);
    assert!(at_zero_q.entropy_j_m3_k > 0.0);
    assert_eq!(
        family.at_temperature(500.0).unwrap().response(&coaxial(0.25), &history).unwrap(),
        base.at_temperature(500.0).unwrap().response(&coaxial(0.25), &history).unwrap()
    );
    for (old_t, new_t) in [(300.0, 500.0), (500.0, 300.0)] {
        let step = family.checked_increment(&f, &history, old_t, &f, new_t, 1e-3).unwrap();
        assert_eq!(step.mechanical.state, history);
        assert_eq!(step.mechanical.plastic_increment, 0.0);
        near(step.coefficient_storage_change_j_m3, slope * (new_t - old_t) * shape(history.alpha), 1e-8, 2e-13);
        near(step.thermoelastic_storage_change_j_m3, -step.coefficient_storage_change_j_m3, 1e-8, 2e-13);
        near(step.internal_energy_change_j_m3(), 0.0, 2e-8, 0.0);
        near(step.caloric_heat_j_m3(), 0.0, 2e-8, 0.0);
        near(step.balance_error_j_m3, 0.0, 2e-8, 0.0);
    }
}

#[test]
fn thermal_saturation_remains_native_for_perzyna_and_simultaneous_norton_j2() {
    let parameters = constant();
    let base = AffineThermalFiniteJ2::new(parameters, limits()).unwrap();
    let temperature = 410.0;
    let duration = 0.12;
    let f = coaxial(0.2);
    for combined in [false, true] {
        let with_hardening = base.clone().with_voce_hardening(Q, B, DQ).unwrap();
        let (family, opposite_order) = if combined {
            (
                with_hardening.with_norton_creep(kinetics(), 3.0).unwrap(),
                base.clone().with_norton_creep(kinetics(), 3.0).unwrap().with_voce_hardening(Q, B, DQ).unwrap(),
            )
        } else {
            (
                with_hardening.with_overstress(kinetics()).unwrap(),
                base.clone().with_overstress(kinetics()).unwrap().with_voce_hardening(Q, B, DQ).unwrap(),
            )
        };
        assert_eq!(family, opposite_order);
        assert!(family.at_temperature(temperature).is_err());
        let old = family.initial_state();
        let native = FiniteJ2::new(parameters.reference, limits()).unwrap()
            .with_voce_hardening(DQ.mul_add(temperature - 300.0, Q), B).unwrap();
        let rate = FiniteOverstressParameters {
            viscosity_pa_s: family.overstress_viscosity_pa_s(temperature).unwrap().unwrap(),
            stress_scale_pa: kinetics().stress_scale_pa,
            exponent: kinetics().exponent,
        };
        let expected = if combined {
            FiniteNortonJ2::new(native, rate, duration, 3.0).unwrap().response(&f, &old).unwrap()
        } else {
            FiniteOverstressJ2::new(native, rate, duration).unwrap().response(&f, &old).unwrap()
        };
        let step = family.checked_increment_for_step(&ID, &old, 300.0, &f, temperature, 1e-3, duration).unwrap();
        assert_eq!(step.mechanical, expected);
        assert!(step.mechanical.plastic_increment > 0.0);
        if combined {
            assert!(step.mechanical.state.creep_alpha > 0.0);
        }
        near(step.reversible_heat_j_m3, DQ * temperature * shape(step.mechanical.state.alpha), 1e-7, 2e-13);
        near(step.caloric_heat_j_m3(), expected.inelastic_dissipation_j_m3() + step.reversible_heat_j_m3, 1e-8, 2e-13);
        near(step.balance_error_j_m3, 0.0, 2e-7, 0.0);
    }
}

#[test]
fn saturation_admission_covers_the_whole_thermal_span_and_refuses_pure_creep() {
    let family = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    for (q, b, slope) in [
        (-1.0, B, 0.0), (Q, 0.0, 0.0), (Q, f64::NAN, 0.0),
        (Q, B, f64::INFINITY), (Q, B, -70_000.0), (Q, B, 700_000.0),
    ] {
        assert!(family.clone().with_voce_hardening(q, b, slope).is_err());
    }
    let selected = family.with_voce_hardening(Q, B, DQ).unwrap();
    assert_eq!(selected.clone().with_voce_hardening(Q, B, DQ).unwrap(), selected);
    assert!(selected.with_voce_hardening(2.0 * Q, B, DQ).is_err());
    let p = parameters();
    let creep = AffineThermalFiniteJ2::new_norton(
        FiniteNortonThermalParameters {
            reference_temperature_k: p.reference_temperature_k,
            temperature_range_k: p.temperature_range_k,
            shear_modulus_pa: p.reference.shear_modulus_pa,
            bulk_modulus_pa: p.reference.bulk_modulus_pa,
            shear_modulus_slope_pa_k: p.slopes.shear_modulus_pa_k,
            bulk_modulus_slope_pa_k: p.slopes.bulk_modulus_pa_k,
            logarithmic_expansion_per_k: p.logarithmic_expansion_per_k,
            logarithmic_expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        }, limits(), kinetics(),
    ).unwrap();
    assert!(creep.with_voce_hardening(Q, B, DQ).is_err());
}

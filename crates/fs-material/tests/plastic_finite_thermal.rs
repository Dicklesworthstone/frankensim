//! Native thermal mappings, fixed-history thermodynamics and physical heat.
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
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
        e.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}, difference {:.8e}",
        (actual - expected).abs()
    );
}

#[test]
fn zero_thermal_coefficients_recover_native_response_and_never_heat_signed_work_defect() {
    let family = AffineThermalFiniteJ2::new(constant(), limits()).unwrap();
    let native = FiniteJ2::new(parameters().reference, limits()).unwrap();
    let f = [1.2, 0.12, 0.01, 0.0, 0.94, 0.04, 0.02, 0.0, 0.95];
    let old = native.initial_state();
    let expected = native.response(&f, &old).unwrap();
    assert_eq!(
        family
            .at_temperature(410.0)
            .unwrap()
            .response(&f, &old)
            .unwrap(),
        expected
    );
    let increment = family
        .checked_increment(&ID, &old, 300.0, &f, 410.0, 1e-4)
        .unwrap();
    assert_eq!(increment.coefficient_storage_change_j_m3, 0.0);
    assert_eq!(increment.thermoelastic_storage_change_j_m3, 0.0);
    assert_eq!(increment.reversible_heat_j_m3, 0.0);
    assert_eq!(
        increment.caloric_heat_j_m3(),
        expected.plastic_dissipation_j_m3
    );
    near(increment.balance_error_j_m3, 0.0, 2e-7, 0.0);

    // A deliberately wide synthetic Hencky envelope exposes nonconvexity in F.
    // Its negative endpoint quadrature defect must not become heat or be clipped.
    let wide = AffineThermalFiniteJ2::new(
        constant(),
        FiniteJ2Limits {
            maximum_elastic_log_strain: 3.0,
            ..limits()
        },
    )
    .unwrap();
    let f0 = ID.map(|x| 3.0 * x);
    let f1 = ID.map(|x| 4.0 * x);
    let step = wide
        .checked_increment(&f0, &old, 300.0, &f1, 300.0, 1e-3)
        .unwrap();
    assert!(step.endpoint_work_defect_j_m3 < -1e6);
    assert_eq!(step.mechanical.plastic_increment, 0.0);
    assert_eq!(step.caloric_heat_j_m3(), 0.0);
    near(step.balance_error_j_m3, 0.0, 2e-6, 0.0);
}

#[test]
fn free_expansion_and_clamped_heating_have_independent_hencky_reference_states() {
    let p = parameters();
    let family = AffineThermalFiniteJ2::new(p, limits()).unwrap();
    let t = 380.0;
    let law = family.at_temperature(t).unwrap();
    let delta = t - 300.0;
    let log_theta = p.logarithmic_expansion_per_k * delta
        + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta;
    let theta = log_theta.exp();
    near(law.thermal_stretch(), theta, 2e-16, 2e-16);
    let (s, c) = 0.73_f64.sin_cos();
    let q = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let free = q.map(|x| theta * x);
    let initial = family.initial_state();
    let response = law.response(&free, &initial).unwrap();
    assert_eq!(response.state, initial);
    assert!(response.first_piola_pa.iter().all(|x| x.abs() < 1e-6));
    let state = law.thermodynamics(&free, &initial).unwrap();
    near(state.free_energy_j_m3, 0.0, 1e-18, 0.0);
    near(state.entropy_j_m3_k, 0.0, 2e-9, 0.0);
    near(state.jacobian, theta.powi(3), 1e-15, 1e-15);
    near(state.mechanical_jacobian, 1.0, 1e-15, 0.0);
    let bulk = p.reference.bulk_modulus_pa + p.slopes.bulk_modulus_pa_k * delta;
    let a = p.logarithmic_expansion_per_k + p.logarithmic_expansion_slope_per_k2 * delta;
    near(
        state.heat_capacity_j_m3_k,
        -9.0 * t * bulk * a * a,
        1e-8,
        2e-13,
    );

    let clamp = law.thermodynamics(&ID, &initial).unwrap();
    let g = -3.0 * log_theta;
    for i in [0, 4, 8] {
        near(clamp.first_piola_pa[i], bulk * g, 1e-7, 3e-13);
    }
    near(clamp.free_energy_j_m3, 0.5 * bulk * g * g, 1e-8, 3e-13);
    near(
        clamp.entropy_j_m3_k,
        3.0 * bulk * a * g - 0.5 * p.slopes.bulk_modulus_pa_k * g * g,
        1e-8,
        3e-13,
    );
    assert!(clamp.heat_capacity_j_m3_k >= family.minimum_mechanical_heat_capacity_j_m3_k());
}

#[test]
fn frozen_thermal_piola_tangent_and_potential_use_the_exact_total_f_chain() {
    let family = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    let initial = family.initial_state();
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&coaxial(0.15), &initial)
        .unwrap()
        .state;
    let law = family.at_temperature(410.0).unwrap();
    let f = [1.18, 0.15, 0.03, -0.04, 0.94, 0.07, 0.02, -0.03, 1.0];
    let original = old.clone();
    let r = law.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0);
    assert_eq!(
        law.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    for col in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        let mut minus = f;
        plus[col] += h;
        minus[col] -= h;
        let a = law.stress_and_potential(&plus, &old).unwrap();
        let b = law.stress_and_potential(&minus, &old).unwrap();
        near((a.1 - b.1) / (2.0 * h), r.first_piola_pa[col], 0.1, 4e-7);
        for row in 0..9 {
            near(
                (a.0[row] - b.0[row]) / (2.0 * h),
                r.tangent_pa[row][col],
                8.0,
                3e-6,
            );
            near(r.tangent_pa[row][col], r.tangent_pa[col][row], 2e-6, 3e-12);
        }
    }
    assert_eq!(old, original);
}

#[test]
fn mechanical_entropy_internal_energy_capacity_and_maxwell_derivatives_hold_at_fixed_history() {
    let family = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    let f = [1.12, 0.1, 0.02, -0.03, 0.95, 0.05, 0.01, 0.0, 1.01];
    let history = family
        .at_temperature(300.0)
        .unwrap()
        .response(&f, &family.initial_state())
        .unwrap()
        .state;
    let t = 370.0;
    let h = 0.05;
    let center = family
        .at_temperature(t)
        .unwrap()
        .thermodynamics(&f, &history)
        .unwrap();
    let plus = family
        .at_temperature(t + h)
        .unwrap()
        .thermodynamics(&f, &history)
        .unwrap();
    let minus = family
        .at_temperature(t - h)
        .unwrap()
        .thermodynamics(&f, &history)
        .unwrap();
    near(
        -(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * h),
        center.entropy_j_m3_k,
        0.01,
        2e-7,
    );
    near(
        -t * (plus.free_energy_j_m3 - 2.0 * center.free_energy_j_m3 + minus.free_energy_j_m3)
            / (h * h),
        center.heat_capacity_j_m3_k,
        0.03,
        2e-5,
    );
    near(
        center.internal_energy_j_m3,
        center.free_energy_j_m3 + t * center.entropy_j_m3_k,
        2e-9,
        2e-15,
    );
    for i in 0..9 {
        near(
            (plus.first_piola_pa[i] - minus.first_piola_pa[i]) / (2.0 * h),
            center.first_piola_temperature_derivative_fixed_history_pa_k[i],
            0.01,
            2e-7,
        );
        let mut fp = f;
        let mut fm = f;
        fp[i] += 1e-6;
        fm[i] -= 1e-6;
        let fixed = family.at_temperature(t).unwrap();
        let ep = fixed.thermodynamics(&fp, &history).unwrap().entropy_j_m3_k;
        let em = fixed.thermodynamics(&fm, &history).unwrap().entropy_j_m3_k;
        near(
            -(ep - em) / 2e-6,
            center.first_piola_temperature_derivative_fixed_history_pa_k[i],
            0.03,
            4e-7,
        );
    }
    assert!(center.heat_capacity_j_m3_k >= family.minimum_mechanical_heat_capacity_j_m3_k());
}

#[test]
fn heating_can_reactivate_yield_and_affine_hardening_storage_is_not_deposited_as_heat() {
    let family = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    let cold = family.at_temperature(300.0).unwrap();
    let f = coaxial(0.15);
    let old = cold.response(&f, &family.initial_state()).unwrap().state;
    let saved = old.clone();
    assert!(
        family
            .at_temperature(450.0)
            .unwrap()
            .validate_committed(&f, &old, 1e-4)
            .is_err()
    );
    let step = family
        .checked_increment(&f, &old, 300.0, &f, 450.0, 1e-4)
        .unwrap();
    assert!(step.mechanical.plastic_increment > 0.0);
    let yield_at_hot =
        parameters().reference.yield_stress_pa + 150.0 * parameters().slopes.yield_stress_pa_k;
    near(
        step.mechanical.plastic_dissipation_j_m3,
        yield_at_hot * step.mechanical.plastic_increment,
        1e-8,
        2e-14,
    );
    assert_eq!(step.mechanical_work_j_m3, 0.0);
    near(
        step.mechanical_work_j_m3,
        step.internal_energy_change_j_m3()
            + step.caloric_heat_j_m3()
            + step.endpoint_work_defect_j_m3,
        2e-7,
        2e-13,
    );
    near(step.balance_error_j_m3, 0.0, 2e-7, 0.0);
    assert_eq!(old, saved);

    let mut only_h = constant();
    only_h.slopes.hardening_modulus_pa_k = -10_000.0;
    let hardening = AffineThermalFiniteJ2::new(only_h, limits()).unwrap();
    let retained = hardening
        .at_temperature(300.0)
        .unwrap()
        .response(&f, &hardening.initial_state())
        .unwrap()
        .state;
    let relaxed = coaxial(retained.alpha);
    let change = hardening
        .checked_increment(&relaxed, &retained, 300.0, &relaxed, 400.0, 1e-4)
        .unwrap();
    assert_eq!(change.mechanical.state, retained);
    assert!(change.coefficient_storage_change_j_m3 < 0.0);
    near(
        change.coefficient_storage_change_j_m3,
        -0.5 * 1e6 * retained.alpha.powi(2),
        1e-8,
        1e-14,
    );
    near(
        change.thermoelastic_storage_change_j_m3,
        -change.coefficient_storage_change_j_m3,
        1e-8,
        1e-14,
    );
    near(change.caloric_heat_j_m3(), 0.0, 1e-8, 0.0);
    near(change.internal_energy_change_j_m3(), 0.0, 1e-8, 0.0);
}

#[test]
fn temperature_and_mechanical_jacobian_domains_are_explicit_and_history_is_unchanged_on_refusal() {
    let family = AffineThermalFiniteJ2::new(parameters(), limits()).unwrap();
    for t in [0.0, 279.0, 501.0, f64::NAN] {
        assert!(family.at_temperature(t).is_err());
    }
    let mut bad = parameters();
    bad.slopes.shear_modulus_pa_k = -1e6;
    assert!(AffineThermalFiniteJ2::new(bad, limits()).is_err());
    bad = parameters();
    bad.logarithmic_expansion_per_k = 20.0;
    assert!(AffineThermalFiniteJ2::new(bad, limits()).is_err());
    let mut expanded = constant();
    expanded.logarithmic_expansion_per_k = 0.002;
    let family = AffineThermalFiniteJ2::new(
        expanded,
        FiniteJ2Limits {
            maximum_elastic_log_strain: 2.0,
            ..limits()
        },
    )
    .unwrap();
    let state = family.initial_state();
    let stretch = 0.095_f64.powf(1.0 / 3.0);
    let total = ID.map(|x| stretch * x);
    let cool = family.at_temperature(280.0).unwrap();
    let r = cool.response(&total, &state).unwrap();
    assert!(r.jacobian < limits().minimum_jacobian);
    assert!(
        cool.thermodynamics(&total, &state)
            .unwrap()
            .mechanical_jacobian
            > limits().minimum_jacobian
    );
    let snapshot = state.clone();
    assert!(
        family
            .at_temperature(400.0)
            .unwrap()
            .response(&total, &state)
            .is_err()
    );
    assert_eq!(state, snapshot);
}

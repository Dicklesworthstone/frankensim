//! G0/G1/G3: authored porous thermomechanics, not experimental calibration.
use fs_material::plastic::finite::porous::thermal::{
    AffineThermalFinitePorous, FinitePorousTemperatureSlopes, FinitePorousThermalParameters,
};
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const Q: f64 = 12e6;
const B: f64 = 3.0;
const Q_T: f64 = -30_000.0;

fn limits() -> FinitePorousLimits {
    FinitePorousLimits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}
fn parameters() -> FinitePorousThermalParameters {
    FinitePorousThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        reference: FinitePorousParameters {
            shear_modulus_pa: 80e6,
            bulk_modulus_pa: 160e6,
            yield_stress_pa: 2e6,
            hardening_modulus_pa: 8e6,
            pressure_sensitivity: 1.2,
        },
        slopes: FinitePorousTemperatureSlopes {
            shear_modulus_pa_k: -20_000.0,
            bulk_modulus_pa_k: -30_000.0,
            yield_stress_pa_k: -4_000.0,
            hardening_modulus_pa_k: -10_000.0,
        },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}
fn constant() -> FinitePorousThermalParameters {
    FinitePorousThermalParameters {
        slopes: FinitePorousTemperatureSlopes {
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
fn family() -> AffineThermalFinitePorous {
    AffineThermalFinitePorous::new(parameters(), limits())
        .unwrap()
        .with_voce_hardening(Q, B, Q_T)
        .unwrap()
}
fn diagonal(stretch: f64) -> [f64; 9] {
    ID.map(|x| stretch * x)
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}; absolute {absolute}, relative {relative}"
    );
}

#[test]
fn heated_hydrostatic_crush_and_unload_keep_separate_thermal_elastic_and_plastic_volumes() {
    let p = parameters();
    let family = family();
    let old = family.initial_state();
    let stretch = 0.8_f64;
    let temperature = 400.0;
    let delta = temperature - p.reference_temperature_k;
    let log_theta = delta
        * (p.logarithmic_expansion_per_k + 0.5 * delta * p.logarithmic_expansion_slope_per_k2);
    let theta = log_theta.exp();
    let bulk = p.reference.bulk_modulus_pa + delta * p.slopes.bulk_modulus_pa_k;
    let yield_stress = p.reference.yield_stress_pa + delta * p.slopes.yield_stress_pa_k;
    let hardening = p.reference.hardening_modulus_pa + delta * p.slopes.hardening_modulus_pa_k;
    let q = Q + delta * Q_T;
    let beta = p.reference.pressure_sensitivity
        / (1.0 + p.reference.pressure_sensitivity.powi(2) / 9.0).sqrt();
    let trial_g = 3.0 * (stretch.ln() - log_theta);
    let radius = |alpha: f64| yield_stress + hardening * alpha + q * -(-B * alpha).exp_m1();
    // Independent scalar hydrostatic stress/volume equation, with both the
    // thermal stretch and accumulated plastic volume explicit.
    let (mut lo, mut hi) = (0.0, -trial_g / beta);
    for _ in 0..100 {
        let alpha = 0.5 * (lo + hi);
        if bulk * beta * (-trial_g - beta * alpha) > radius(alpha) {
            lo = alpha;
        } else {
            hi = alpha;
        }
    }
    let alpha = 0.5 * (lo + hi);
    let tau = -radius(alpha) / beta;
    let g = tau / bulk;
    let jp = (-beta * alpha).exp();
    let shape = alpha + (-B * alpha).exp_m1() / B;
    let a = p.logarithmic_expansion_per_k + delta * p.logarithmic_expansion_slope_per_k2;
    let expected_entropy = 3.0 * bulk * a * g
        - 0.5 * p.slopes.bulk_modulus_pa_k * g * g
        - 0.5 * p.slopes.hardening_modulus_pa_k * alpha * alpha
        - Q_T * shape;
    let step = family
        .checked_increment(&ID, &old, 300.0, &diagonal(stretch), temperature, 1e-4)
        .unwrap();
    near(step.mechanical.state.alpha, alpha, 2e-12, 0.0);
    near(step.mechanical.plastic_jacobian, jp, 3e-12, 0.0);
    near(
        step.mechanical.first_piola_pa[0],
        tau / stretch,
        2e-5,
        2e-12,
    );
    near(step.new.jacobian, stretch.powi(3), 2e-14, 0.0);
    near(
        step.new.mechanical_jacobian,
        (stretch / theta).powi(3),
        2e-14,
        0.0,
    );
    near(step.new.elastic_jacobian, g.exp(), 2e-12, 0.0);
    near(step.new.plastic_jacobian, jp, 3e-12, 0.0);
    near(
        step.new.free_energy_j_m3,
        0.5 * tau * tau / bulk + 0.5 * hardening * alpha * alpha + q * shape,
        2e-5,
        3e-12,
    );
    near(step.new.entropy_j_m3_k, expected_entropy, 1e-7, 3e-12);
    near(
        step.new.internal_energy_j_m3,
        step.new.free_energy_j_m3 + temperature * expected_entropy,
        3e-5,
        3e-12,
    );
    near(
        step.mechanical.plastic_dissipation_j_m3,
        yield_stress * alpha,
        1e-6,
        3e-12,
    );
    near(step.balance_error_j_m3, 0.0, 1e-6, 0.0);
    assert_eq!(old, family.initial_state());

    // Free unloading at the SAME hot temperature preserves compaction, while
    // total current volume still includes reversible thermal expansion.
    let released_f = diagonal(theta * jp.cbrt());
    let released = family
        .checked_increment(
            &diagonal(stretch),
            &step.mechanical.state,
            temperature,
            &released_f,
            temperature,
            1e-4,
        )
        .unwrap();
    assert_eq!(released.mechanical.state, step.mechanical.state);
    assert_eq!(released.mechanical.plastic_dissipation_j_m3, 0.0);
    near(released.new.jacobian, theta.powi(3) * jp, 3e-12, 0.0);
    near(released.new.elastic_jacobian, 1.0, 3e-12, 0.0);
    for stress in released.new.first_piola_pa {
        near(stress, 0.0, 2e-5, 0.0);
    }
    near(
        released.new.free_energy_j_m3,
        0.5 * hardening * alpha * alpha + q * shape,
        3e-5,
        3e-12,
    );
}

#[test]
fn thermal_softening_at_fixed_total_deformation_returns_old_history_once() {
    let p = FinitePorousThermalParameters {
        slopes: FinitePorousTemperatureSlopes {
            yield_stress_pa_k: -4_000.0,
            ..constant().slopes
        },
        ..constant()
    };
    let family = AffineThermalFinitePorous::new(p, limits()).unwrap();
    let f = diagonal(0.9);
    let cold = family
        .at_temperature(300.0)
        .unwrap()
        .response(&f, &family.initial_state())
        .unwrap();
    let original = cold.state.clone();
    assert!(
        family
            .validate_committed(&f, &cold.state, 450.0, 1e-4)
            .is_err()
    );
    let hot = family
        .checked_increment(&f, &cold.state, 300.0, &f, 450.0, 1e-4)
        .unwrap();
    assert!(hot.mechanical.plastic_increment > 0.0);
    assert!(hot.mechanical.plastic_jacobian < cold.plastic_jacobian);
    assert_eq!(hot.mechanical_work_j_m3, 0.0);
    assert_eq!(hot.coefficient_storage_change_j_m3, 0.0);
    assert_eq!(hot.reversible_heat_j_m3, 0.0);
    let yield_hot = p.reference.yield_stress_pa + 150.0 * p.slopes.yield_stress_pa_k;
    near(
        hot.caloric_heat_j_m3(),
        yield_hot * hot.mechanical.plastic_increment,
        1e-7,
        2e-13,
    );
    assert!(hot.endpoint_work_defect_j_m3.abs() > 1.0);
    near(
        hot.mechanical_work_j_m3
            - hot.internal_energy_change_j_m3()
            - hot.caloric_heat_j_m3()
            - hot.endpoint_work_defect_j_m3,
        0.0,
        1e-6,
        0.0,
    );
    assert_eq!(cold.state, original);
}

#[test]
fn noncoaxial_hot_return_and_fixed_history_thermodynamics_have_consistent_derivatives() {
    let family = family();
    let crushed = family
        .at_temperature(300.0)
        .unwrap()
        .response(&diagonal(0.8), &family.initial_state())
        .unwrap()
        .state;
    let old_f = [0.89, 0.11, -0.02, 0.03, 0.93, 0.06, -0.01, 0.04, 0.84];
    let old = family
        .at_temperature(300.0)
        .unwrap()
        .response(&old_f, &crushed)
        .unwrap()
        .state;
    let f = [1.03, -0.06, 0.04, 0.09, 0.91, 0.03, 0.02, -0.01, 0.87];
    let law = family.at_temperature(410.0).unwrap();
    let response = law.response(&f, &old).unwrap();
    assert!(response.plastic_increment > 0.0);
    assert!(response.plastic_jacobian < 0.9);
    for column in 0..9 {
        let h = 2e-6;
        let (mut plus, mut minus) = (f, f);
        plus[column] += h;
        minus[column] -= h;
        let (pa, wa) = law.stress_and_potential(&plus, &old).unwrap();
        let (pb, wb) = law.stress_and_potential(&minus, &old).unwrap();
        near(
            (wa - wb) / (2.0 * h),
            response.first_piola_pa[column],
            1.0,
            2e-6,
        );
        for row in 0..9 {
            near(
                (pa[row] - pb[row]) / (2.0 * h),
                response.tangent_pa[row][column],
                8.0,
                3e-6,
            );
            near(
                response.tangent_pa[row][column],
                response.tangent_pa[column][row],
                2e-5,
                5e-12,
            );
        }
    }
    let history = &response.state;
    let center = family.thermodynamics(&f, history, 410.0).unwrap();
    let h = 0.025;
    let plus = family.thermodynamics(&f, history, 410.0 + h).unwrap();
    let minus = family.thermodynamics(&f, history, 410.0 - h).unwrap();
    near(
        -(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * h),
        center.entropy_j_m3_k,
        1e-3,
        2e-7,
    );
    near(
        (plus.internal_energy_j_m3 - minus.internal_energy_j_m3) / (2.0 * h),
        center.heat_capacity_j_m3_k,
        1e-2,
        3e-7,
    );
    assert!(center.heat_capacity_j_m3_k >= family.minimum_mechanical_heat_capacity_j_m3_k());
    for column in 0..9 {
        let derivative = center.first_piola_temperature_derivative_fixed_history_pa_k[column];
        near(
            (plus.first_piola_pa[column] - minus.first_piola_pa[column]) / (2.0 * h),
            derivative,
            2e-3,
            3e-7,
        );
        let (mut a, mut b) = (f, f);
        a[column] += 1e-6;
        b[column] -= 1e-6;
        let entropy_a = family
            .thermodynamics(&a, history, 410.0)
            .unwrap()
            .entropy_j_m3_k;
        let entropy_b = family
            .thermodynamics(&b, history, 410.0)
            .unwrap()
            .entropy_j_m3_k;
        near(-(entropy_a - entropy_b) / 2e-6, derivative, 0.1, 3e-6);
    }
}

#[test]
fn isothermal_limit_is_exact_and_zero_saturation_endpoint_keeps_its_entropy_derivative() {
    let family = AffineThermalFinitePorous::new(constant(), limits())
        .unwrap()
        .with_voce_hardening(Q, B, 0.0)
        .unwrap();
    let native = FinitePorous::new(constant().reference, limits())
        .unwrap()
        .with_voce_hardening(Q, B)
        .unwrap();
    let f = [0.91, 0.08, -0.02, 0.01, 0.94, 0.04, 0.0, 0.03, 0.87];
    let initial = native.initial_state();
    let expected = native.response(&f, &initial).unwrap();
    assert_eq!(
        family
            .at_temperature(410.0)
            .unwrap()
            .response(&f, &initial)
            .unwrap(),
        expected
    );
    let step = family
        .checked_increment(&ID, &initial, 300.0, &f, 410.0, 1e-4)
        .unwrap();
    assert_eq!(step.caloric_heat_j_m3(), expected.plastic_dissipation_j_m3);
    assert_eq!(step.reversible_heat_j_m3, 0.0);
    assert_eq!(family.minimum_mechanical_heat_capacity_j_m3_k(), 0.0);

    let plain = AffineThermalFinitePorous::new(constant(), limits()).unwrap();
    let vanishing = plain.clone().with_voce_hardening(6e6, B, Q_T).unwrap();
    let state = vanishing
        .at_temperature(300.0)
        .unwrap()
        .response(&diagonal(0.85), &vanishing.initial_state())
        .unwrap()
        .state;
    let a = vanishing.thermodynamics(&f, &state, 500.0).unwrap();
    let b = plain.thermodynamics(&f, &state, 500.0).unwrap();
    assert_eq!(a.free_energy_j_m3, b.free_energy_j_m3);
    let shape = state.alpha + (-B * state.alpha).exp_m1() / B;
    near(
        a.entropy_j_m3_k - b.entropy_j_m3_k,
        -Q_T * shape,
        1e-8,
        2e-13,
    );
    near(
        a.internal_energy_j_m3 - b.internal_energy_j_m3,
        -500.0 * Q_T * shape,
        1e-5,
        2e-13,
    );
    assert_eq!(a.heat_capacity_j_m3_k, b.heat_capacity_j_m3_k);
}

#[test]
fn whole_temperature_domain_coefficients_and_state_transitions_are_admitted_explicitly() {
    let mut invalid = parameters();
    invalid.slopes.bulk_modulus_pa_k = -1e6;
    assert!(AffineThermalFinitePorous::new(invalid, limits()).is_err());
    invalid = parameters();
    invalid.logarithmic_expansion_per_k = f64::NAN;
    assert!(AffineThermalFinitePorous::new(invalid, limits()).is_err());
    let family = family();
    let saved = family.clone();
    for t in [0.0, 279.0, 501.0, f64::NAN, f64::INFINITY] {
        assert!(family.at_temperature(t).is_err());
    }
    assert!(family.clone().with_voce_hardening(Q + 1.0, B, Q_T).is_err());
    assert!(
        AffineThermalFinitePorous::new(parameters(), limits())
            .unwrap()
            .with_voce_hardening(1e6, B, -10_000.0)
            .is_err()
    );
    assert!(
        AffineThermalFinitePorous::new(parameters(), limits())
            .unwrap()
            .with_voce_hardening(Q, 0.0, Q_T)
            .is_err()
    );
    let old = family.initial_state();
    let mut bad_history = old.clone();
    bad_history.inverse_plastic_cauchy_green[0] = 2.0;
    assert!(family.thermodynamics(&ID, &bad_history, 300.0).is_err());
    assert!(
        family
            .checked_increment(&ID, &old, 300.0, &diagonal(0.8), 501.0, 1e-4)
            .is_err()
    );
    assert_eq!(family, saved);
    assert_eq!(old, family.initial_state());
}

mod source {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
    };
    use fs_material::plastic::finite::porous::thermal::FiniteThermalPorousMaterial;
    use fs_material::plastic::finite::porous::thermal::source::{
        ResolvedAffineThermalFinitePorous, THERMAL_FINITE_POROUS_LAW_ID,
        THERMAL_FINITE_POROUS_VOCE_LAW_ID,
    };
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Dims, QuantitySpec};

    const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
    const PA_K: Dims = Dims([-1, 1, -2, -1, 0, 0]);
    const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
    fn quantity(kind: QuantityKind) -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
    }
    fn context() -> QueryPoint {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 1e5)
            .unwrap()
            .with("lot", 23.0)
            .unwrap()
    }
    fn card() -> ConstitutiveModelCard {
        let p = parameters();
        let b = limits();
        let source =
            fs_blake3::hash_bytes(b"authored finite porous thermal test; no empirical calibration");
        ConstitutiveModelCard {
            law: LawId(THERMAL_FINITE_POROUS_VOCE_LAW_ID.into()),
            law_version: 1,
            parameters: [
                ("reference_temperature_k", p.reference_temperature_k, K),
                ("minimum_temperature_k", p.temperature_range_k[0], K),
                ("maximum_temperature_k", p.temperature_range_k[1], K),
                ("shear_modulus_pa", p.reference.shear_modulus_pa, PA),
                ("bulk_modulus_pa", p.reference.bulk_modulus_pa, PA),
                ("yield_stress_pa", p.reference.yield_stress_pa, PA),
                ("hardening_modulus_pa", p.reference.hardening_modulus_pa, PA),
                (
                    "pressure_sensitivity",
                    p.reference.pressure_sensitivity,
                    Dims::NONE,
                ),
                (
                    "shear_modulus_slope_pa_k",
                    p.slopes.shear_modulus_pa_k,
                    PA_K,
                ),
                ("bulk_modulus_slope_pa_k", p.slopes.bulk_modulus_pa_k, PA_K),
                ("yield_stress_slope_pa_k", p.slopes.yield_stress_pa_k, PA_K),
                (
                    "hardening_modulus_slope_pa_k",
                    p.slopes.hardening_modulus_pa_k,
                    PA_K,
                ),
                (
                    "logarithmic_expansion_per_k",
                    p.logarithmic_expansion_per_k,
                    Dims([0, 0, 0, -1, 0, 0]),
                ),
                (
                    "logarithmic_expansion_slope_per_k2",
                    p.logarithmic_expansion_slope_per_k2,
                    Dims([0, 0, 0, -2, 0, 0]),
                ),
                (
                    "minimum_mechanical_jacobian",
                    b.minimum_jacobian,
                    Dims::NONE,
                ),
                (
                    "maximum_elastic_log_strain",
                    b.maximum_elastic_log_strain,
                    Dims::NONE,
                ),
                (
                    "maximum_plastic_log_strain",
                    b.maximum_plastic_log_strain,
                    Dims::NONE,
                ),
                (
                    "maximum_accumulated_plastic_strain",
                    b.maximum_accumulated_plastic_strain,
                    Dims::NONE,
                ),
                ("voce_saturation_pa", Q, PA),
                ("voce_rate_per_strain", B, Dims::NONE),
                ("voce_saturation_slope_pa_k", Q_T, PA_K),
            ]
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
            state_schema_version: 1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    280.0,
                    500.0,
                )
                .with_quantity("pressure", quantity(QuantityKind::Pressure), 5e4, 2e5),
            sources: vec![source],
            provenance: Provenance {
                source: "synthetic porous thermal law".into(),
                license: "test fixture".into(),
                artifact: Some(source),
            },
        }
    }

    #[test]
    fn source_family_retains_original_context_and_native_hot_compaction_response() {
        let card = card();
        let context = context();
        let resolved = ResolvedAffineThermalFinitePorous::try_from_model_card_at(
            &card,
            &context,
            "temperature",
        )
        .unwrap();
        assert_eq!(resolved.model_card(), &card);
        assert_eq!(resolved.context(), &context);
        assert_eq!(resolved.temperature_axis(), "temperature");
        let original = family();
        assert_eq!(resolved.law(), &original);
        let state = original.initial_state();
        let expected = original
            .checked_increment(&ID, &state, 300.0, &diagonal(0.8), 400.0, 1e-4)
            .unwrap();
        let sourced = resolved
            .thermal_finite_porous()
            .checked_increment(&ID, &state, 300.0, &diagonal(0.8), 400.0, 1e-4)
            .unwrap();
        assert_eq!(sourced, expected);
        assert!(sourced.new.plastic_jacobian < 0.9);
        assert_eq!(resolved.model_card(), &card);
        assert_eq!(resolved.context(), &context);
    }

    #[test]
    fn source_requires_exact_coefficients_units_selector_and_full_typed_temperature_domain() {
        let original = card();
        let context = context();
        let refuses = |card: &ConstitutiveModelCard, point: &QueryPoint| {
            assert!(
                ResolvedAffineThermalFinitePorous::try_from_model_card_at(
                    card,
                    point,
                    "temperature"
                )
                .is_err()
            );
        };
        let mut missing = original.clone();
        missing.parameters.remove("pressure_sensitivity");
        refuses(&missing, &context);
        let mut mistyped = original.clone();
        mistyped
            .parameters
            .get_mut("logarithmic_expansion_per_k")
            .unwrap()
            .dims = Dims::NONE;
        refuses(&mistyped, &context);
        let mut linear = original.clone();
        linear.law = LawId(THERMAL_FINITE_POROUS_LAW_ID.into());
        refuses(&linear, &context);
        for name in [
            "voce_saturation_pa",
            "voce_rate_per_strain",
            "voce_saturation_slope_pa_k",
        ] {
            linear.parameters.remove(name);
        }
        assert!(
            ResolvedAffineThermalFinitePorous::try_from_model_card_at(
                &linear,
                &context,
                "temperature"
            )
            .is_ok()
        );
        let mut limited = original.clone();
        limited.validity = limited.validity.with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            280.0,
            499.0,
        );
        refuses(&limited, &context);
        refuses(
            &original,
            &context.clone().with("temperature", 350.0).unwrap(),
        );
        let untyped = QueryPoint::new()
            .with("temperature", 300.0)
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 1e5)
            .unwrap();
        refuses(&original, &untyped);
        let mut second = original.clone();
        second.validity = second.validity.with_quantity(
            "second_temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            280.0,
            500.0,
        );
        let second_context = context
            .clone()
            .with_quantity(
                "second_temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
            )
            .unwrap();
        refuses(&second, &second_context);
        let mut version = original.clone();
        version.state_schema_version = 2;
        refuses(&version, &context);
        let mut isothermal = original.clone();
        isothermal.law =
            LawId(fs_material::plastic::finite::porous::source::FINITE_POROUS_VOCE_LAW_ID.into());
        refuses(&isothermal, &context);
    }
}

//! Synthetic secondary-creep references, not empirical material calibrations.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::FiniteJ2Limits;
use fs_material::plastic::finite::thermal::source::{
    ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_J2_LAW_ID, THERMAL_FINITE_NORTON_LAW_ID,
    THERMAL_PERZYNA_FINITE_J2_LAW_ID,
};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteJ2ThermalOverstress,
    FiniteNortonThermalParameters,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.1,
    maximum_elastic_log_strain: 0.8,
    maximum_plastic_log_strain: 3.0,
    maximum_accumulated_plastic_strain: 3.0,
};
fn parameters() -> FiniteNortonThermalParameters {
    FiniteNortonThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        shear_modulus_pa: 80e6,
        bulk_modulus_pa: 160e6,
        shear_modulus_slope_pa_k: -20_000.0,
        bulk_modulus_slope_pa_k: -30_000.0,
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}
fn constant() -> FiniteNortonThermalParameters {
    FiniteNortonThermalParameters {
        shear_modulus_slope_pa_k: 0.0,
        bulk_modulus_slope_pa_k: 0.0,
        logarithmic_expansion_per_k: 0.0,
        logarithmic_expansion_slope_per_k2: 0.0,
        ..parameters()
    }
}
fn kinetics(exponent: f64) -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 10e6,
        stress_scale_pa: 2e6,
        exponent,
        activation_temperature_k: 1800.0,
    }
}
fn coaxial(mean: f64, axial: f64) -> [f64; 9] {
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
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn activated_finite_creep_flows_below_plastic_yield_with_independent_thermal_volume_oracle() {
    let p = parameters();
    let k = kinetics(1.0);
    let family = AffineThermalFiniteJ2::new_norton(p, LIMITS, k).unwrap();
    let mut threshold_parameters = *family.parameters();
    threshold_parameters.reference.yield_stress_pa = 2e6;
    let plastic = AffineThermalFiniteJ2::new(threshold_parameters, LIMITS)
        .unwrap()
        .with_overstress(k)
        .unwrap();
    let (mean, axial, dt) = (0.03, 0.002, 0.7);
    let mut previous_gamma = 0.0;
    for temperature in [300.0, 420.0] {
        let delta = temperature - p.reference_temperature_k;
        let mu = p.shear_modulus_pa + p.shear_modulus_slope_pa_k * delta;
        let bulk = p.bulk_modulus_pa + p.bulk_modulus_slope_pa_k * delta;
        let theta = (p.logarithmic_expansion_per_k * delta
            + 0.5 * p.logarithmic_expansion_slope_per_k2 * delta * delta)
            .exp();
        let eta = k.reference_viscosity_pa_s
            * (k.activation_temperature_k * (1.0 / temperature - 1.0 / 300.0)).exp();
        let f = coaxial(mean, axial).map(|x| theta * x);
        let gamma = 3.0 * mu * axial / (3.0 * mu + eta / dt);
        let elastic_axial = axial - gamma;
        let q = eta * gamma / dt;
        let stored = 1.5 * mu * elastic_axial * elastic_axial + 0.5 * bulk * (3.0 * mean).powi(2);
        let loss = q * gamma;
        let initial = family.initial_state();
        let response = family
            .at_temperature_for_step(temperature, dt)
            .unwrap()
            .response(&f, &initial)
            .unwrap();
        assert!(q < threshold_parameters.reference.yield_stress_pa);
        assert!(response.plastic_increment > previous_gamma);
        previous_gamma = response.plastic_increment;
        near(response.plastic_increment, gamma, 2e-14, 2e-11);
        near(response.state.alpha, gamma, 2e-14, 2e-11);
        near(
            response.state.inverse_plastic_cauchy_green[0],
            (-2.0 * gamma).exp(),
            2e-14,
            2e-11,
        );
        for i in [1, 2] {
            near(
                response.state.inverse_plastic_cauchy_green[i],
                gamma.exp(),
                2e-14,
                2e-11,
            );
        }
        near(response.equivalent_kirchhoff_pa, q, 2e-6, 2e-11);
        near(
            response.jacobian,
            theta.powi(3) * (3.0 * mean).exp(),
            2e-14,
            2e-12,
        );
        near(response.stored_energy_j_m3, stored, 2e-6, 2e-12);
        near(response.plastic_dissipation_j_m3, loss, 2e-7, 2e-11);
        near(
            response.incremental_potential_j_m3,
            stored + 0.5 * loss,
            2e-6,
            2e-12,
        );
        for (i, tau) in [
            (0, 3.0 * bulk * mean + 2.0 * mu * elastic_axial),
            (4, 3.0 * bulk * mean - mu * elastic_axial),
            (8, 3.0 * bulk * mean - mu * elastic_axial),
        ] {
            near(response.first_piola_pa[i], tau / f[i], 2e-6, 2e-12);
        }
        let no_plastic_flow = plastic
            .at_temperature_for_step(temperature, dt)
            .unwrap()
            .response(&f, &plastic.initial_state())
            .unwrap();
        assert_eq!(no_plastic_flow.plastic_increment, 0.0);
        assert_eq!(no_plastic_flow.plastic_dissipation_j_m3, 0.0);
        assert_eq!(initial, family.initial_state());
    }
}

#[test]
fn sustained_small_kirchhoff_stress_accumulates_large_creep_without_stored_hardening() {
    let p = constant();
    let k = FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: 4e6,
        activation_temperature_k: 0.0,
        ..kinetics(1.0)
    };
    let family = AffineThermalFiniteJ2::new_norton(p, LIMITS, k).unwrap();
    let (q, dt) = (0.4e6, 0.3);
    let elastic_axial = q / (3.0 * p.shear_modulus_pa);
    let gamma = q * dt / k.reference_viscosity_pa_s;
    let mut state = family.initial_state();
    let mut old_f = coaxial(0.0, elastic_axial);
    let mut alpha = 0.0;
    for _ in 0..16 {
        alpha += gamma;
        let f = coaxial(0.0, alpha + elastic_axial);
        let saved = state.clone();
        let step = family
            .checked_increment_for_step(&old_f, &state, 300.0, &f, 300.0, 0.0, dt)
            .unwrap();
        near(step.mechanical.state.alpha, alpha, 3e-13, 3e-12);
        near(step.mechanical.equivalent_kirchhoff_pa, q, 2e-5, 3e-11);
        near(step.caloric_heat_j_m3(), q * gamma, 2e-6, 3e-11);
        near(
            step.new.free_energy_j_m3,
            q * q / (6.0 * p.shear_modulus_pa),
            2e-6,
            3e-11,
        );
        near(step.internal_energy_change_j_m3(), 0.0, 2e-6, 0.0);
        near(step.balance_error_j_m3, 0.0, 2e-6, 0.0);
        near(
            step.mechanical.incremental_potential_j_m3,
            step.new.free_energy_j_m3 + 0.5 * step.caloric_heat_j_m3(),
            2e-6,
            3e-11,
        );
        assert_eq!(state, saved);
        state = step.mechanical.state;
        old_f = f;
    }
    assert!(old_f[0] > 1.5 && state.alpha > 0.45);
    // Complete creep history remains in a stress-free permanent shape; alpha
    // contributes no hidden hardening or thermal storage to this mechanism.
    let free_shape = coaxial(0.0, alpha);
    let stored = family.thermodynamics(&free_shape, &state, 300.0).unwrap();
    near(stored.free_energy_j_m3, 0.0, 1e-12, 0.0);
    assert!(stored.first_piola_pa.iter().all(|x| x.abs() < 1e-6));
    let mut longer_history = state.clone();
    longer_history.alpha += 0.2;
    assert_eq!(
        family
            .thermodynamics(&old_f, &longer_history, 300.0)
            .unwrap(),
        family.thermodynamics(&old_f, &state, 300.0).unwrap()
    );
}

#[test]
fn activated_power_creep_keeps_exact_thermal_tangent_and_physical_two_temperature_work() {
    let family = AffineThermalFiniteJ2::new_norton(parameters(), LIMITS, kinetics(2.4)).unwrap();
    let old_f = coaxial(0.01, 0.06);
    let old = family
        .at_temperature_for_step(300.0, 0.1)
        .unwrap()
        .response(&old_f, &family.initial_state())
        .unwrap()
        .state;
    let f = [1.12, 0.09, 0.02, -0.02, 0.96, 0.04, 0.01, -0.01, 1.03];
    let (temperature, duration) = (410.0, 0.4);
    let law = family
        .at_temperature_for_step(temperature, duration)
        .unwrap();
    let response = law.response(&f, &old).unwrap();
    let gamma = response.plastic_increment;
    let eta =
        kinetics(2.4).reference_viscosity_pa_s * (1800.0 * (1.0 / temperature - 1.0 / 300.0)).exp();
    let q = 2e6 * (eta * gamma / (duration * 2e6)).powf(1.0 / 2.4);
    near(response.equivalent_kirchhoff_pa, q, 3e-6, 3e-12);
    near(response.plastic_dissipation_j_m3, q * gamma, 2e-6, 3e-12);
    near(
        response.incremental_potential_j_m3,
        response.stored_energy_j_m3 + (2.4 / 3.4) * q * gamma,
        2e-6,
        3e-12,
    );
    for col in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        plus[col] += h;
        let mut minus = f;
        minus[col] -= h;
        let a = law.stress_and_potential(&plus, &old).unwrap();
        let b = law.stress_and_potential(&minus, &old).unwrap();
        near(
            (a.1 - b.1) / (2.0 * h),
            response.first_piola_pa[col],
            0.1,
            4e-7,
        );
        for row in 0..9 {
            near(
                (a.0[row] - b.0[row]) / (2.0 * h),
                response.tangent_pa[row][col],
                8.0,
                3e-6,
            );
            near(
                response.tangent_pa[row][col],
                response.tangent_pa[col][row],
                3e-6,
                3e-12,
            );
        }
    }
    let step = family
        .checked_increment_for_step(&old_f, &old, 300.0, &f, temperature, 0.0, duration)
        .unwrap();
    assert_eq!(step.mechanical, response);
    assert_ne!(step.reversible_heat_j_m3, 0.0);
    assert_ne!(step.thermoelastic_storage_change_j_m3, 0.0);
    assert_eq!(
        step.caloric_heat_j_m3(),
        step.mechanical.plastic_dissipation_j_m3 + step.reversible_heat_j_m3
    );
    near(
        step.mechanical_work_j_m3,
        step.internal_energy_change_j_m3()
            + step.caloric_heat_j_m3()
            + step.endpoint_work_defect_j_m3,
        2e-6,
        3e-12,
    );
    let other_kinetics = family
        .clone()
        .with_overstress(FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: 1e3,
            activation_temperature_k: 0.0,
            ..kinetics(1.0)
        })
        .unwrap();
    let fixed = family.thermodynamics(&f, &old, temperature).unwrap();
    assert_eq!(
        other_kinetics
            .thermodynamics(&f, &old, temperature)
            .unwrap(),
        fixed
    );
    let h = 0.01;
    let plus = family.thermodynamics(&f, &old, temperature + h).unwrap();
    let minus = family.thermodynamics(&f, &old, temperature - h).unwrap();
    near(
        -(plus.free_energy_j_m3 - minus.free_energy_j_m3) / (2.0 * h),
        fixed.entropy_j_m3_k,
        1e-3,
        1e-7,
    );
    near(
        (plus.internal_energy_j_m3 - minus.internal_energy_j_m3) / (2.0 * h),
        fixed.heat_capacity_j_m3_k,
        1e-3,
        1e-7,
    );
    for i in 0..9 {
        near(
            (plus.first_piola_pa[i] - minus.first_piola_pa[i]) / (2.0 * h),
            fixed.first_piola_temperature_derivative_fixed_history_pa_k[i],
            2e-3,
            2e-7,
        );
    }
}

#[test]
fn norton_is_explicit_requires_time_and_preserves_old_plastic_admission_and_creep_bounds() {
    let family = AffineThermalFiniteJ2::new_norton(constant(), LIMITS, kinetics(1.0)).unwrap();
    assert!(family.is_norton_creep());
    assert!(AffineThermalFiniteJ2::new(*family.parameters(), LIMITS).is_err());
    let mut plastic_parameters = *family.parameters();
    plastic_parameters.reference.yield_stress_pa = 2e6;
    assert!(
        !AffineThermalFiniteJ2::new(plastic_parameters, LIMITS)
            .unwrap()
            .is_norton_creep()
    );
    assert_eq!(family.parameters().reference.hardening_modulus_pa, 0.0);
    assert!(matches!(
        family.at_temperature(300.0),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    let old = family.initial_state();
    assert!(matches!(
        family.checked_increment(&ID, &old, 300.0, &ID, 300.0, 0.0),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    let inspector = family.at_temperature_for_inspection(300.0).unwrap();
    inspector
        .validate_committed(&coaxial(0.0, 0.001), &old, 0.0)
        .unwrap();
    assert!(matches!(
        inspector.response(&ID, &old),
        Err(FiniteJ2ThermalError::DurationRequired)
    ));
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(family.at_temperature_for_step(300.0, dt).is_err());
    }
    let limited = AffineThermalFiniteJ2::new_norton(
        constant(),
        FiniteJ2Limits {
            maximum_accumulated_plastic_strain: 1e-4,
            ..LIMITS
        },
        kinetics(1.0),
    )
    .unwrap();
    assert!(
        limited
            .checked_increment_for_step(&ID, &old, 300.0, &coaxial(0.0, 0.01), 300.0, 0.0, 1.0)
            .is_err()
    );
    assert_eq!(old, family.initial_state());
    for k in [
        FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: 0.0,
            ..kinetics(1.0)
        },
        FiniteJ2ThermalOverstress {
            exponent: 0.9,
            ..kinetics(1.0)
        },
        FiniteJ2ThermalOverstress {
            activation_temperature_k: 1e9,
            ..kinetics(1.0)
        },
    ] {
        assert!(AffineThermalFiniteJ2::new_norton(parameters(), LIMITS, k).is_err());
    }
}

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature_quantity(), 300.0)
        .unwrap()
}
fn card() -> ConstitutiveModelCard {
    let p = parameters();
    let k = kinetics(2.4);
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let pa_k = Dims([-1, 1, -2, -1, 0, 0]);
    let kelvin = Dims([0, 0, 0, 1, 0, 0]);
    let hash =
        fs_blake3::hash_bytes(b"authored finite Norton thermal fixture, no measured calibration");
    ConstitutiveModelCard {
        law: LawId(THERMAL_FINITE_NORTON_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("reference_temperature_k", p.reference_temperature_k, kelvin),
            ("minimum_temperature_k", p.temperature_range_k[0], kelvin),
            ("maximum_temperature_k", p.temperature_range_k[1], kelvin),
            ("shear_modulus_pa", p.shear_modulus_pa, pa),
            ("bulk_modulus_pa", p.bulk_modulus_pa, pa),
            ("shear_modulus_slope_pa_k", p.shear_modulus_slope_pa_k, pa_k),
            ("bulk_modulus_slope_pa_k", p.bulk_modulus_slope_pa_k, pa_k),
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
                LIMITS.minimum_jacobian,
                Dims::NONE,
            ),
            (
                "maximum_elastic_log_strain",
                LIMITS.maximum_elastic_log_strain,
                Dims::NONE,
            ),
            (
                "maximum_plastic_log_strain",
                LIMITS.maximum_plastic_log_strain,
                Dims::NONE,
            ),
            (
                "maximum_accumulated_plastic_strain",
                LIMITS.maximum_accumulated_plastic_strain,
                Dims::NONE,
            ),
            (
                "creep_reference_viscosity_pa_s",
                k.reference_viscosity_pa_s,
                Dims([-1, 1, -1, 0, 0, 0]),
            ),
            ("creep_stress_scale_pa", k.stress_scale_pa, pa),
            ("creep_exponent", k.exponent, Dims::NONE),
            (
                "creep_activation_temperature_k",
                k.activation_temperature_k,
                kelvin,
            ),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            temperature_quantity(),
            p.temperature_range_k[0],
            p.temperature_range_k[1],
        ),
        sources: vec![hash],
        provenance: Provenance {
            source: "authored finite Norton thermal fixture".into(),
            license: "test fixture".into(),
            artifact: Some(hash),
        },
    }
}

#[test]
fn explicit_seventeen_parameter_norton_source_preserves_identity_and_executes_the_shared_law() {
    let card = card();
    assert_eq!(card.parameters.len(), 17);
    let context = context();
    let source =
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context, "temperature")
            .unwrap();
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.model_card().content_hash(), card.content_hash());
    assert_eq!(source.context(), &context);
    assert_eq!(source.temperature_axis(), "temperature");
    assert!(source.law().is_norton_creep());
    let native = AffineThermalFiniteJ2::new_norton(parameters(), LIMITS, kinetics(2.4)).unwrap();
    assert_eq!(source.law(), &native);
    let f = coaxial(0.02, 0.003);
    let result = source
        .law()
        .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 0.0, 0.4)
        .unwrap();
    assert!(result.mechanical.plastic_increment > 0.0);
    assert!(result.mechanical.plastic_dissipation_j_m3 > 0.0);
    assert_eq!(
        result,
        native
            .checked_increment_for_step(&ID, &native.initial_state(), 300.0, &f, 410.0, 0.0, 0.4,)
            .unwrap()
    );
}

#[test]
fn norton_source_rejects_missing_coefficients_units_plastic_aliases_and_uncovered_temperatures() {
    let original = card();
    let resolve = |card: &ConstitutiveModelCard| {
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(card, &context(), "temperature")
    };
    for name in original.parameters.keys() {
        let mut missing = original.clone();
        missing.parameters.remove(name);
        assert!(resolve(&missing).is_err(), "missing {name} accepted");
        let mut wrong = original.clone();
        wrong.parameters.get_mut(name).unwrap().dims = Dims([1, 0, 0, 0, 0, 0]);
        assert!(
            resolve(&wrong).is_err(),
            "wrong dimensions accepted for {name}"
        );
    }
    for name in [
        "yield_stress_pa",
        "hardening_modulus_pa",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
        "overstress_exponent",
    ] {
        let mut extra = original.clone();
        extra.parameters.insert(
            name.into(),
            LawParameter {
                value: 0.0,
                dims: Dims::NONE,
            },
        );
        assert!(resolve(&extra).is_err(), "unexpected {name} accepted");
    }
    for name in [THERMAL_FINITE_J2_LAW_ID, THERMAL_PERZYNA_FINITE_J2_LAW_ID] {
        let mut legacy = original.clone();
        legacy.law = LawId(name.into());
        assert!(resolve(&legacy).is_err());
    }
    for (name, value) in [
        ("creep_reference_viscosity_pa_s", 0.0),
        ("creep_stress_scale_pa", -1.0),
        ("creep_exponent", 0.99),
        ("creep_activation_temperature_k", -1.0),
        ("creep_activation_temperature_k", 1e9),
        ("shear_modulus_slope_pa_k", -1e6),
    ] {
        let mut invalid = original.clone();
        invalid.parameters.get_mut(name).unwrap().value = value;
        assert!(
            resolve(&invalid).is_err(),
            "invalid {name}={value} accepted"
        );
    }
    let mut narrow = original.clone();
    narrow.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        temperature_quantity(),
        300.0,
        450.0,
    );
    assert!(resolve(&narrow).is_err());
    let shifted = context().with("temperature", 310.0).unwrap();
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&original, &shifted, "temperature")
            .is_err()
    );
    let mut version = original.clone();
    version.law_version = 2;
    assert!(resolve(&version).is_err());
    version.law_version = 1;
    version.state_schema_version = 2;
    assert!(resolve(&version).is_err());
}

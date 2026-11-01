//! Synthetic finite thermoplastic source transport, not experimental calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::thermal::source::{
    ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_J2_LAW_ID,
};
use fs_material::plastic::finite::thermal::{
    FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters, FiniteThermalJ2Material,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_qty::semantic::{
    FrequencyConvention, QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm,
};
use fs_qty::{Dims, QuantitySpec};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PA_K: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const INV_K: Dims = Dims([0, 0, 0, -1, 0, 0]);
const INV_K2: Dims = Dims([0, 0, 0, -2, 0, 0]);
const T_REF: f64 = 300.0;
const T_SPAN: [f64; 2] = [280.0, 420.0];
const REFERENCE: FiniteJ2Parameters = FiniteJ2Parameters {
    shear_modulus_pa: 80.0e6,
    bulk_modulus_pa: 160.0e6,
    yield_stress_pa: 2.0e6,
    hardening_modulus_pa: 8.0e6,
};
const SLOPES: FiniteJ2TemperatureSlopes = FiniteJ2TemperatureSlopes {
    shear_modulus_pa_k: -40.0e3,
    bulk_modulus_pa_k: -60.0e3,
    yield_stress_pa_k: -2.0e3,
    hardening_modulus_pa_k: -5.0e3,
};
const EXPANSION: f64 = 1.0e-4;
const EXPANSION_SLOPE: f64 = 1.0e-7;
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.25,
    maximum_elastic_log_strain: 0.5,
    maximum_plastic_log_strain: 2.0,
    maximum_accumulated_plastic_strain: 3.0,
};

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            T_REF,
        )
        .unwrap()
}

fn fixture() -> ConstitutiveModelCard {
    let source = fs_blake3::hash_bytes(
        b"synthetic finite thermal J2 original-volume Kirchhoff fixture; no measured calibration",
    );
    ConstitutiveModelCard {
        law: LawId(THERMAL_FINITE_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("reference_temperature_k", T_REF, K),
            ("minimum_temperature_k", T_SPAN[0], K),
            ("maximum_temperature_k", T_SPAN[1], K),
            ("shear_modulus_pa", REFERENCE.shear_modulus_pa, PA),
            ("bulk_modulus_pa", REFERENCE.bulk_modulus_pa, PA),
            ("yield_stress_pa", REFERENCE.yield_stress_pa, PA),
            ("hardening_modulus_pa", REFERENCE.hardening_modulus_pa, PA),
            ("shear_modulus_slope_pa_k", SLOPES.shear_modulus_pa_k, PA_K),
            ("bulk_modulus_slope_pa_k", SLOPES.bulk_modulus_pa_k, PA_K),
            ("yield_stress_slope_pa_k", SLOPES.yield_stress_pa_k, PA_K),
            (
                "hardening_modulus_slope_pa_k",
                SLOPES.hardening_modulus_pa_k,
                PA_K,
            ),
            ("logarithmic_expansion_per_k", EXPANSION, INV_K),
            (
                "logarithmic_expansion_slope_per_k2",
                EXPANSION_SLOPE,
                INV_K2,
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
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            T_SPAN[0],
            T_SPAN[1],
        ),
        sources: vec![source],
        provenance: Provenance {
            source: "authored synthetic finite thermoplastic source fixture".into(),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    }
}

fn resolve(card: &ConstitutiveModelCard) -> ResolvedAffineThermalFiniteJ2 {
    ResolvedAffineThermalFiniteJ2::try_from_model_card_at(card, &context(), "temperature").unwrap()
}

fn coefficients(temperature: f64) -> FiniteJ2Parameters {
    let delta = temperature - T_REF;
    FiniteJ2Parameters {
        shear_modulus_pa: REFERENCE.shear_modulus_pa + delta * SLOPES.shear_modulus_pa_k,
        bulk_modulus_pa: REFERENCE.bulk_modulus_pa + delta * SLOPES.bulk_modulus_pa_k,
        yield_stress_pa: REFERENCE.yield_stress_pa + delta * SLOPES.yield_stress_pa_k,
        hardening_modulus_pa: REFERENCE.hardening_modulus_pa
            + delta * SLOPES.hardening_modulus_pa_k,
    }
}

fn log_thermal_stretch(temperature: f64) -> f64 {
    let delta = temperature - T_REF;
    EXPANSION * delta + 0.5 * EXPANSION_SLOPE * delta * delta
}

fn deformation(temperature: f64, volume_log: f64, axial_deviator: f64) -> [f64; 9] {
    let log_theta = log_thermal_stretch(temperature);
    [
        (log_theta + volume_log + axial_deviator).exp(),
        0.0,
        0.0,
        0.0,
        (log_theta + volume_log - 0.5 * axial_deviator).exp(),
        0.0,
        0.0,
        0.0,
        (log_theta + volume_log - 0.5 * axial_deviator).exp(),
    ]
}

fn close_at(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn close(actual: f64, expected: f64) {
    close_at(actual, expected, 2.0e-10);
}

#[test]
fn exact_source_family_drives_large_plastic_return_at_expanded_nonunit_volume() {
    let card = fixture();
    let source = resolve(&card);
    assert_eq!(source.model_card(), &card);
    assert_eq!(source.model_card().content_hash(), card.content_hash());
    assert_eq!(source.context(), &context());
    assert_eq!(source.temperature_axis(), "temperature");
    assert_eq!(
        source.law().parameters(),
        &FiniteJ2ThermalParameters {
            reference_temperature_k: T_REF,
            temperature_range_k: T_SPAN,
            reference: REFERENCE,
            slopes: SLOPES,
            logarithmic_expansion_per_k: EXPANSION,
            logarithmic_expansion_slope_per_k2: EXPANSION_SLOPE,
        }
    );
    assert_eq!(source.law().limits(), &LIMITS);
    assert!(std::ptr::eq(source.thermal_finite_j2(), source.law()));

    let (temperature, v, a) = (350.0, 0.04, 0.22);
    let p = coefficients(temperature);
    let frozen = source.law().at_temperature(temperature).unwrap();
    close(frozen.temperature_k(), temperature);
    close(
        frozen.log_thermal_stretch(),
        log_thermal_stretch(temperature),
    );
    close(
        frozen.thermal_stretch(),
        log_thermal_stretch(temperature).exp(),
    );
    let old = source.law().initial_state();
    let f = deformation(temperature, v, a);
    let response = frozen.response(&f, &old).unwrap();
    let increment = (3.0 * p.shear_modulus_pa * a - p.yield_stress_pa)
        / (3.0 * p.shear_modulus_pa + p.hardening_modulus_pa);
    let returned = a - increment;
    let pressure = 3.0 * p.bulk_modulus_pa * v;
    close(
        response.jacobian,
        (3.0 * (v + log_thermal_stretch(temperature))).exp(),
    );
    assert!(response.jacobian > 1.1);
    assert!(increment > 0.2);
    close(response.plastic_increment, increment);
    close(response.state.alpha, increment);
    close(
        response.state.inverse_plastic_cauchy_green[0],
        (-2.0 * increment).exp(),
    );
    close(
        response.state.inverse_plastic_cauchy_green[1],
        increment.exp(),
    );
    close(
        response.state.inverse_plastic_cauchy_green[2],
        increment.exp(),
    );
    close(
        response.equivalent_kirchhoff_pa,
        p.yield_stress_pa + p.hardening_modulus_pa * increment,
    );
    for (index, tau) in [
        (0, pressure + 2.0 * p.shear_modulus_pa * returned),
        (4, pressure - p.shear_modulus_pa * returned),
        (8, pressure - p.shear_modulus_pa * returned),
    ] {
        close(response.kirchhoff_pa[index], tau);
        close(response.first_piola_pa[index], tau / f[index]);
    }
    // All terms use the original volume. Multiplying by theta^3 here would
    // alter both thermal storage and the explicitly selected Kirchhoff law.
    let stored = 1.5 * p.shear_modulus_pa * returned.powi(2)
        + 0.5 * p.bulk_modulus_pa * (3.0 * v).powi(2)
        + 0.5 * p.hardening_modulus_pa * increment.powi(2);
    close(response.stored_energy_j_m3, stored);
    close(
        response.plastic_dissipation_j_m3,
        p.yield_stress_pa * increment,
    );
    close(
        response.incremental_potential_j_m3,
        stored + p.yield_stress_pa * increment,
    );
    assert_eq!(old, source.law().initial_state());
    frozen
        .validate_committed(&f, &response.state, 1.0e-5)
        .unwrap();
}

#[test]
fn fixed_history_entropy_internal_energy_and_stress_derivatives_follow_original_source_slopes() {
    let source = resolve(&fixture());
    let (temperature, v, a) = (350.0, 0.04, 0.22);
    let p = coefficients(temperature);
    let frozen = source.law().at_temperature(temperature).unwrap();
    let f = deformation(temperature, v, a);
    let accepted = frozen
        .response(&f, &source.law().initial_state())
        .unwrap()
        .state;
    let saved = accepted.clone();
    let thermo = frozen.thermodynamics(&f, &accepted).unwrap();
    let g = 3.0 * v;
    let shape = 1.5 * (a - accepted.alpha).powi(2);
    let expansion = EXPANSION + EXPANSION_SLOPE * (temperature - T_REF);
    let free = p.shear_modulus_pa * shape
        + 0.5 * p.bulk_modulus_pa * g * g
        + 0.5 * p.hardening_modulus_pa * accepted.alpha.powi(2);
    let entropy = 3.0 * p.bulk_modulus_pa * expansion * g
        - SLOPES.shear_modulus_pa_k * shape
        - 0.5 * SLOPES.bulk_modulus_pa_k * g * g
        - 0.5 * SLOPES.hardening_modulus_pa_k * accepted.alpha.powi(2);
    let capacity = temperature
        * ((6.0 * SLOPES.bulk_modulus_pa_k * expansion
            + 3.0 * p.bulk_modulus_pa * EXPANSION_SLOPE)
            * g
            - 9.0 * p.bulk_modulus_pa * expansion * expansion);
    close(thermo.free_energy_j_m3, free);
    close(thermo.entropy_j_m3_k, entropy);
    close(thermo.internal_energy_j_m3, free + temperature * entropy);
    close(thermo.heat_capacity_j_m3_k, capacity);
    close(thermo.mechanical_jacobian, g.exp());
    close(
        thermo.jacobian,
        (g + 3.0 * log_thermal_stretch(temperature)).exp(),
    );
    assert!(capacity < 0.0);
    assert!(source.law().minimum_mechanical_heat_capacity_j_m3_k() <= capacity);

    // Differentiate at fixed TOTAL F and complete retained history. Performing
    // a new return at either perturbed T would not compute the entropy derivative.
    let h = 0.01;
    let inspect = |t| {
        source
            .law()
            .at_temperature(t)
            .unwrap()
            .thermodynamics(&f, &accepted)
            .unwrap()
    };
    let low = inspect(temperature - h);
    let high = inspect(temperature + h);
    close_at(
        -(high.free_energy_j_m3 - low.free_energy_j_m3) / (2.0 * h),
        entropy,
        2.0e-7,
    );
    close_at(
        (high.internal_energy_j_m3 - low.internal_energy_j_m3) / (2.0 * h),
        capacity,
        2.0e-7,
    );
    close_at(
        temperature * (high.entropy_j_m3_k - low.entropy_j_m3_k) / (2.0 * h),
        capacity,
        2.0e-7,
    );
    for index in 0..9 {
        close_at(
            (high.first_piola_pa[index] - low.first_piola_pa[index]) / (2.0 * h),
            thermo.first_piola_temperature_derivative_fixed_history_pa_k[index],
            2.0e-7,
        );
    }
    assert_eq!(accepted, saved);

    // The source's expansion convention gives an exactly stress-free expanded
    // virgin state, including at the family endpoints.
    for t in [T_SPAN[0], temperature, T_SPAN[1]] {
        let frozen = source.law().at_temperature(t).unwrap();
        let free = frozen
            .thermodynamics(&deformation(t, 0.0, 0.0), &source.law().initial_state())
            .unwrap();
        close(free.free_energy_j_m3, 0.0);
        close(free.entropy_j_m3_k, 0.0);
        for stress in free.first_piola_pa {
            close_at(stress, 0.0, 1.0e-6);
        }
    }
}

#[test]
fn source_domain_covers_the_entire_native_temperature_span_and_exact_fixed_context() {
    let pressure = semantic(QuantityKind::Pressure);
    let mut card = fixture();
    card.validity = card
        .validity
        .with_quantity("environment_pressure", pressure, 1.0e5, 2.0e5);
    let point = context()
        .with_quantity("environment_pressure", pressure, 1.5e5)
        .unwrap()
        .with("unconstrained_note_coordinate", 17.0)
        .unwrap();
    let source =
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &point, "temperature")
            .unwrap();
    assert_eq!(source.context(), &point);
    assert_eq!(source.model_card(), &card);
    let mut other = card.clone();
    other.provenance.source = "independently authored synthetic calibration identity".into();
    let other =
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&other, &point, "temperature")
            .unwrap();
    assert_eq!(source.law(), other.law());
    assert_ne!(
        source.model_card().content_hash(),
        other.model_card().content_hash()
    );
    for bad_point in [
        QueryPoint::new(),
        point.clone().with("temperature", T_REF + 1.0).unwrap(),
        point.clone().with("environment_pressure", 0.9e5).unwrap(),
        point
            .clone()
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                T_REF,
            )
            .unwrap(),
        point
            .clone()
            .with_quantity(
                "environment_pressure",
                semantic(QuantityKind::Stress),
                1.5e5,
            )
            .unwrap(),
    ] {
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &bad_point, "temperature")
                .is_err()
        );
    }
    for (lo, hi) in [(281.0, 420.0), (280.0, 419.0), (299.0, 301.0)] {
        let mut restricted = fixture();
        restricted.validity = ValidityDomain::unconstrained().with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            lo,
            hi,
        );
        // The reference point is covered, but that cannot grant a longer path.
        assert!(
            restricted
                .validity
                .contains_typed(context().axes(), context().axis_quantities())
        );
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
                &restricted,
                &context(),
                "temperature"
            )
            .is_err()
        );
    }
    let mut untyped = fixture();
    untyped.validity = ValidityDomain::unconstrained().with("temperature", T_SPAN[0], T_SPAN[1]);
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&untyped, &context(), "temperature")
            .is_err()
    );
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&fixture(), &context(), "T").is_err()
    );
}

#[test]
fn evolving_mechanical_rate_and_secondary_temperature_source_axes_are_refused() {
    for (axis, quantity, value) in [
        ("stress", semantic(QuantityKind::Stress), 1.0),
        (
            "strain",
            semantic(QuantityKind::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Normal,
            }),
            0.1,
        ),
        ("stretch", QuantitySpec::dimensional(Dims::NONE), 1.0),
        (
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            1.0,
        ),
        (
            "temperature_rate",
            QuantitySpec::dimensional(Dims([0, 0, -1, 1, 0, 0])),
            1.0,
        ),
        (
            "other_temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            T_REF,
        ),
        (
            "pressure_sample",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Instantaneous,
            )),
            1.0,
        ),
    ] {
        let mut card = fixture();
        card.validity = card.validity.with_quantity(axis, quantity, 0.0, 500.0);
        let point = context().with_quantity(axis, quantity, value).unwrap();
        assert!(
            card.validity
                .contains_typed(point.axes(), point.axis_quantities())
        );
        let error =
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &point, "temperature")
                .unwrap_err();
        assert!(
            error.to_string().contains("fixed external condition"),
            "{axis}: {error}"
        );
    }
    let mut card = fixture();
    card.validity = card.validity.with("T", 290.0, 310.0);
    let point = context().with("T", T_REF).unwrap();
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &point, "temperature")
            .is_err()
    );
}

#[test]
fn all_seventeen_parameters_are_explicit_dimensioned_finite_and_native_admissible() {
    let original = fixture();
    assert_eq!(original.parameters.len(), 17);
    for name in original.parameters.keys() {
        let mut card = original.clone();
        card.parameters.remove(name);
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "missing {name}"
        );
        let mut card = original.clone();
        let parameter = card.parameters.get_mut(name).unwrap();
        parameter.dims = if parameter.dims == PA { Dims::NONE } else { PA };
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "units {name}"
        );
        let mut card = original.clone();
        card.parameters.get_mut(name).unwrap().value = f64::NAN;
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "finite {name}"
        );
    }
    for (name, value) in [
        ("reference_temperature_k", 430.0),
        ("minimum_temperature_k", 0.0),
        ("maximum_temperature_k", 280.0),
        ("shear_modulus_slope_pa_k", -1.0e6),
        ("bulk_modulus_slope_pa_k", -2.0e6),
        ("yield_stress_slope_pa_k", -1.0e5),
        ("hardening_modulus_slope_pa_k", -1.0e5),
        ("minimum_mechanical_jacobian", 1.0),
        ("maximum_elastic_log_strain", 0.0),
        ("maximum_plastic_log_strain", 0.0),
        ("maximum_accumulated_plastic_strain", 0.0),
    ] {
        let mut card = original.clone();
        card.parameters.get_mut(name).unwrap().value = value;
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "native domain {name}"
        );
    }

    let mut card = fixture();
    card.parameters
        .get_mut("minimum_mechanical_jacobian")
        .unwrap()
        .value = 0.995;
    let source = resolve(&card);
    let hot = source.law().at_temperature(350.0).unwrap();
    let old = source.law().initial_state();
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    // Total J=1 passes a total-volume bound, but thermal expansion leaves
    // mechanical J=theta^-3 below this explicitly mechanical source limit.
    assert!(hot.response(&identity, &old).is_err());
    assert!(hot.response(&deformation(350.0, 0.0, 0.0), &old).is_ok());
    assert!(source.law().at_temperature(279.0).is_err());
    assert!(source.law().at_temperature(421.0).is_err());
    assert_eq!(old, source.law().initial_state());
}

#[test]
fn small_strain_static_aliases_and_unowned_history_physics_cannot_select_this_family() {
    for id in [
        "j2-plasticity",
        "finite-j2-hencky-kirchhoff",
        "finite-j2-hencky-cauchy-thermal-affine",
        "finite-j2-hencky-kirchhoff-thermal-voce",
    ] {
        let mut card = fixture();
        card.law = LawId(id.into());
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "{id}"
        );
    }
    let mut version = fixture();
    version.law_version = 2;
    let mut schema = fixture();
    schema.state_schema_version = 2;
    let mut history = fixture();
    history.initial_state = InitialStatePolicy::RequiresDeclaredState;
    for card in [version, schema, history] {
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err()
        );
    }
    for name in [
        "young_modulus_pa",
        "current_volume_shear_modulus_pa",
        "thermal_expansion_per_k",
        "engineering_expansion_per_k",
        "kinematic_hardening_modulus_pa",
        "voce_saturation_pa",
        "taylor_quinney_fraction",
        "plastic_entropy",
    ] {
        let mut card = fixture();
        card.parameters.insert(
            name.into(),
            LawParameter {
                value: 1.0,
                dims: Dims::NONE,
            },
        );
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature")
                .is_err(),
            "{name}"
        );
    }
}

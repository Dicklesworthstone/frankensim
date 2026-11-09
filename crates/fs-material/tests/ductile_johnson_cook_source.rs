//! G0/G1/G3: one original fit drives physical ductility, with no inferred data.
//! All coefficients and source identities below are synthetic software fixtures.
#![deny(unsafe_code)]

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MaterialCard,
    MaterialStateId, Provenance, QueryPoint,
};
use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileJ2, DuctileJ2Conditions, DuctileJ2Settings,
};
use fs_material::ductile::source::{
    DuctileSourceError, JOHNSON_COOK_INITIATION_LAW_ID, PLASTIC_STRAIN_RATE_AXIS,
    STRESS_TRIAXIALITY_AXIS, resolve_johnson_cook_initiation,
};
use fs_material::ductile::{DuctileError, DuctileInitiationLaw, DuctileState};
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_qty::semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const TEMPERATURE: &str = "specimen-temperature";
const PRESSURE: &str = "ambient-pressure";
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn fixed_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, absolute_temperature(), temperature)
        .unwrap()
        .with_quantity(PRESSURE, pressure_quantity(), 101325.0)
        .unwrap()
}

fn point() -> QueryPoint {
    fixed_point(300.0)
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            1.0 / 3.0,
        )
        .unwrap()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            0.001,
        )
        .unwrap()
}

fn model() -> ConstitutiveModelCard {
    let source = hash_bytes(b"synthetic Johnson-Cook fit; not an alloy calibration");
    ConstitutiveModelCard {
        law: LawId(JOHNSON_COOK_INITIATION_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("d1", 0.05, Dims::NONE),
            ("d2", 0.25, Dims::NONE),
            ("d3", 1.5, Dims::NONE),
            ("d4", 0.03, Dims::NONE),
            ("d5", 0.4, Dims::NONE),
            ("reference_strain_rate_s_inv", 0.001, RATE),
            ("transition_temperature_k", 300.0, Temperature::DIMS),
            ("melting_temperature_k", 1500.0, Temperature::DIMS),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                STRESS_TRIAXIALITY_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -1.0,
                2.0,
            )
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(RATE),
                1e-5,
                1.0,
            )
            .with_quantity(TEMPERATURE, absolute_temperature(), 290.0, 1200.0)
            .with_quantity(PRESSURE, pressure_quantity(), 1e4, 1e6),
        sources: vec![source],
        provenance: Provenance {
            source: "synthetic complete minus-triaxiality Johnson-Cook fit".into(),
            license: "test-only".into(),
            artifact: Some(source),
        },
    }
}

fn card(models: Vec<ConstitutiveModelCard>) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "synthetic JC initiation".into(),
            revision: 0,
        },
        ClaimSet::new(),
        models,
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64) {
    assert!((a - b).abs() <= relative * b.abs().max(1e-12), "{a} != {b}");
}

#[test]
fn complete_original_fit_drives_actual_rate_and_temperature_without_source_queries() {
    let original = model();
    let identity = original.content_hash();
    let source = card(vec![original.clone()]);
    let resolved =
        resolve_johnson_cook_initiation(&source, identity, &point(), TEMPERATURE).unwrap();
    assert_eq!(resolved.source().card_identity(), source.content_hash());
    assert_eq!(resolved.source().material_state(), source.id());
    assert_eq!(resolved.source().model_card(), &original);
    assert_eq!(resolved.source().model_identity(), identity);
    assert_eq!(resolved.source().query(), &point());
    assert_eq!(resolved.source().temperature_axis(), TEMPERATURE);
    assert_eq!(resolved.domain().temperature_k, [290.0, 1200.0]);
    assert_eq!(resolved.domain().strain_rate_s_inv, [1e-5, 1.0]);
    for eta in [-1.0, 1.0 / 3.0, 2.0] {
        for rate in [1e-5_f64, 0.001, 0.1, 1.0] {
            for temperature in [290.0_f64, 300.0, 900.0, 1200.0] {
                let expected = (0.05 + 0.25 * (-1.5_f64 * eta).exp())
                    * (1.0 + 0.03 * (rate / 0.001).ln())
                    * (1.0 + 0.4 * ((temperature - 300.0) / 1200.0).clamp(0.0, 1.0));
                close(
                    resolved.fracture_strain(eta, rate, temperature).unwrap(),
                    expected,
                    3e-14,
                );
            }
        }
    }
    let initial = DuctileState::virgin();
    let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let cold_slow = initial
        .trial(&resolved, 0.02, 20.0, 300.0, &stress)
        .unwrap();
    let hot_fast = initial.trial(&resolved, 0.02, 0.2, 900.0, &stress).unwrap();
    assert!(hot_fast.state.utilization() < cold_slow.state.utilization());
    assert_eq!(initial, DuctileState::virgin());
    assert_eq!(resolved.source().model_card(), &original);
}

#[test]
fn source_model_derivatives_match_finite_differences_and_progressive_native_law() {
    let mut original = model();
    original.parameters.get_mut("d1").unwrap().value *= 0.01;
    original.parameters.get_mut("d2").unwrap().value *= 0.01;
    let source = card(vec![original.clone()]);
    let resolved =
        resolve_johnson_cook_initiation(&source, original.content_hash(), &point(), TEMPERATURE)
            .unwrap();
    let (eta, rate, temperature) = (0.7, 0.02, 850.0);
    let gradient = resolved
        .fracture_strain_gradient(eta, rate, temperature)
        .unwrap();
    let e = 1e-6;
    let f = |eta, rate| resolved.fracture_strain(eta, rate, temperature).unwrap();
    close(
        gradient.triaxiality,
        (f(eta + e, rate) - f(eta - e, rate)) / (2.0 * e),
        2e-9,
    );
    close(
        gradient.strain_rate_s_inv,
        (f(eta, rate + e) - f(eta, rate - e)) / (2.0 * e),
        2e-8,
    );
    assert_eq!(gradient.lode_parameter, 0.0);

    let native =
        J2Plasticity::new(IsotropicElastic::new(200e9, 0.3, 0.05).unwrap(), 250e6, 1e9).unwrap();
    let settings = DuctileJ2Settings {
        characteristic_length_m: 0.001,
        post_initiation_plastic_work_j_m2: 1e6,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-4,
            yield_pa: 1.0,
        },
    };
    let numeric = DuctileJ2::new(native.clone(), resolved.law().clone(), settings).unwrap();
    let sourced = DuctileJ2::new(native, resolved, settings).unwrap();
    let strain = [0.01, -0.004, -0.003, 0.002, 0.0, 0.0];
    let conditions = DuctileJ2Conditions {
        duration_s: 1.0,
        temperature_k: 900.0,
    };
    let a = numeric
        .trial(&strain, &numeric.initial_state(), conditions)
        .unwrap();
    let b = sourced
        .trial(&strain, &sourced.initial_state(), conditions)
        .unwrap();
    assert_eq!(a.stress_pa, b.stress_pa);
    assert_eq!(a.branch_tangent_pa, b.branch_tangent_pa);
    assert_eq!(a.energy, b.energy);
    assert_eq!(a.state.initiation(), b.state.initiation());
    assert!(b.state.damage() > 0.0);
}

#[test]
fn material_identity_and_held_environment_remain_bound_while_runtime_temperature_changes() {
    let original = model();
    let source = card(vec![original.clone()]);
    let resolved =
        resolve_johnson_cook_initiation(&source, original.content_hash(), &point(), TEMPERATURE)
            .unwrap();
    let matches = |p: &QueryPoint| {
        resolved
            .source()
            .matches_material_state(source.content_hash(), source.id(), p)
    };
    assert!(matches(&fixed_point(900.0)));
    assert!(matches(&point().with(TEMPERATURE, 900.0).unwrap()));
    assert!(!matches(&fixed_point(1200.1)));
    assert!(!matches(
        &fixed_point(300.0).with(PRESSURE, 101326.0).unwrap()
    ));
    assert!(!matches(&QueryPoint::new()));
    assert!(!matches(
        &fixed_point(300.0).with("extra-coordinate", 1.0).unwrap()
    ));
    assert!(!matches(
        &point().with(PLASTIC_STRAIN_RATE_AXIS, 2.0).unwrap()
    ));
    assert!(!matches(
        &point().with(STRESS_TRIAXIALITY_AXIS, 3.0).unwrap()
    ));
    let wrong_temperature = fixed_point(300.0)
        .with_quantity(
            TEMPERATURE,
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
        )
        .unwrap();
    assert!(!matches(&wrong_temperature));
    let mut other = source.id().clone();
    other.revision = 1;
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        &other,
        &fixed_point(300.0)
    ));
    assert!(!resolved.source().matches_material_state(
        hash_bytes(b"another card"),
        source.id(),
        &fixed_point(300.0)
    ));
}

#[test]
fn explicit_pin_selects_one_coherent_fit_and_never_falls_back() {
    let first = model();
    let mut second = first.clone();
    second.parameters.get_mut("d2").unwrap().value = 0.5;
    second.provenance.source = "second synthetic whole fit".into();
    let source = card(vec![first.clone(), second.clone()]);
    for original in [first, second] {
        let resolved = resolve_johnson_cook_initiation(
            &source,
            original.content_hash(),
            &point(),
            TEMPERATURE,
        )
        .unwrap();
        assert_eq!(resolved.source().model_card(), &original);
        assert_eq!(
            resolved.law().parameters().d2,
            original.parameters["d2"].value
        );
    }
    assert!(
        resolve_johnson_cook_initiation(&source, hash_bytes(b"absent fit"), &point(), TEMPERATURE,)
            .is_err()
    );
}

fn refuses(model: ConstitutiveModelCard, at: &QueryPoint) {
    let pin = model.content_hash();
    assert!(resolve_johnson_cook_initiation(&card(vec![model]), pin, at, TEMPERATURE).is_err());
}

#[test]
fn missing_parameters_units_versions_conventions_and_calibration_sources_refuse() {
    let original = model();
    for key in original.parameters.keys() {
        let mut absent = original.clone();
        absent.parameters.remove(key);
        refuses(absent, &point());
        let mut wrong = original.clone();
        wrong.parameters.get_mut(key).unwrap().dims = Pressure::DIMS;
        refuses(wrong, &point());
    }
    let mut extra = original.clone();
    extra.parameters.insert(
        "d3_original_plus_sign".into(),
        LawParameter {
            value: 1.5,
            dims: Dims::NONE,
        },
    );
    refuses(extra, &point());
    let mut convention = original.clone();
    convention.law = LawId("ductile-johnson-cook-original-plus-triaxiality".into());
    refuses(convention, &point());
    let mut version = original.clone();
    version.law_version += 1;
    refuses(version, &point());
    let mut schema = original.clone();
    schema.state_schema_version += 1;
    refuses(schema, &point());
    let mut history = original.clone();
    history.initial_state = InitialStatePolicy::RequiresDeclaredState;
    refuses(history, &point());
    let mut sources = original;
    sources.sources.clear();
    refuses(sources, &point());
}

#[test]
fn original_complete_calibration_support_and_exact_runtime_quantity_kinds_are_required() {
    let original = model();
    let mut no_bounds = original.clone();
    no_bounds.validity = ValidityDomain::unconstrained();
    refuses(no_bounds, &point());
    let mut invalid_corner = original.clone();
    // The selected reference rate is harmless, but the lowest calibrated rate
    // makes 1+d4*ln(rate/reference) negative. Checking one sample is insufficient.
    invalid_corner.parameters.get_mut("d4").unwrap().value = 1.0;
    refuses(invalid_corner, &point());
    let mut second_temperature = original.clone();
    second_temperature.validity = second_temperature.validity.with_quantity(
        "other-temperature",
        absolute_temperature(),
        300.0,
        300.0,
    );
    refuses(
        second_temperature,
        &point()
            .with_quantity("other-temperature", absolute_temperature(), 300.0)
            .unwrap(),
    );
    let mut evolving = original.clone();
    evolving.validity =
        evolving
            .validity
            .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    refuses(
        evolving,
        &point()
            .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.0)
            .unwrap(),
    );
    let wrong_rate = point()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Frequency(FrequencyConvention::Cyclic),
                ValueForm::Static,
            )),
            0.001,
        )
        .unwrap();
    refuses(original.clone(), &wrong_rate);
    let wrong_temperature = point()
        .with_quantity(
            TEMPERATURE,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap();
    refuses(original.clone(), &wrong_temperature);
    let source = card(vec![original.clone()]);
    let resolved =
        resolve_johnson_cook_initiation(&source, original.content_hash(), &point(), TEMPERATURE)
            .unwrap();
    for (eta, rate, temperature) in [
        (2.01, 0.001, 300.0),
        (0.0, 1.01, 300.0),
        (0.0, 0.001, 1200.1),
    ] {
        assert!(matches!(
            resolved.fracture_strain(eta, rate, temperature),
            Err(DuctileError::Domain { .. })
        ));
    }
    assert!(matches!(
        resolve_johnson_cook_initiation(&source, ContentHash([0; 32]), &point(), TEMPERATURE,),
        Err(DuctileSourceError::Invalid(_))
    ));
}

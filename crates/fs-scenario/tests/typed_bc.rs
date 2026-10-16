//! G0/G3 coverage for typed multiphysics boundary payloads, their admission
//! matrix, semantic-validation accounting, and canonical scenario-IR envelope.

use fs_qty::chemistry::SpeciesId;
use fs_qty::{Dims, QtyAny};
use fs_scenario::bc::{Expectation, expectation};
use fs_scenario::ir::{check_round_trip, parse_ir, write_ir};
use fs_scenario::payload::{
    CharacteristicComponent, CharacteristicDirection, CharacteristicState, OrientationParity,
    OutsideDomainPolicy, Payload, PayloadId, PayloadKind, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, ScalarPayload, SpeciesBundle, SpeciesValue,
    TableInterpolation, VectorPayload,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Compat, Environment, FrameId, IrSourceSpan, LoadCase,
    Physics, Scenario, ScenarioError, ValidationBudget, ValidationError, Violation,
};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const MASS_FLOW: Dims = Dims([0, 1, -1, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const MAGNETIC_VECTOR_POTENTIAL: Dims = Dims([1, 1, -2, 0, -1, 0]);
const MAGNETIC_FLUX_DENSITY: Dims = Dims([0, 1, -2, 0, -1, 0]);
const ELECTRIC_POTENTIAL: Dims = Dims([2, 1, -3, 0, -1, 0]);
const CURRENT_DENSITY: Dims = Dims([-2, 0, 0, 0, 1, 0]);
const SPECIES_AMOUNT_FLUX: Dims = Dims([-2, 0, -1, 0, 0, 1]);
const SPECIES_MASS_FLUX: Dims = Dims([-2, 1, -1, 0, 0, 0]);

fn id(value: &str) -> PayloadId {
    PayloadId::new(value).expect("canonical fixture identifier")
}

fn dimensional_meta(dims: Dims, frame: u32) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Dimensions(dims),
        id("basis/world-cartesian"),
        FrameId(frame),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .expect("valid dimensional metadata")
}

fn heterogeneous_meta(frame: u32) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Heterogeneous,
        id("basis/gas-characteristics"),
        FrameId(frame),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .expect("valid heterogeneous metadata")
}

fn scalar_payload(dims: Dims, frame: u32) -> Payload {
    Payload::Scalar(
        ScalarPayload::new(
            dimensional_meta(dims, frame),
            SampleSource::fixed(QtyAny::new(1.0, dims)),
        )
        .expect("valid scalar payload"),
    )
}

fn vector_payload(dims: Dims, frame: u32) -> Payload {
    Payload::Vector(
        VectorPayload::new(
            dimensional_meta(dims, frame),
            SampleSource::fixed(vec![
                QtyAny::new(1.0, dims),
                QtyAny::new(2.0, dims),
                QtyAny::new(3.0, dims),
            ]),
        )
        .expect("valid vector payload"),
    )
}

fn vector_table_payload(dims: Dims, frame: u32) -> Payload {
    Payload::Vector(
        VectorPayload::new(
            dimensional_meta(dims, frame),
            SampleSource::table(
                vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
                vec![
                    vec![
                        QtyAny::new(1.0, dims),
                        QtyAny::new(2.0, dims),
                        QtyAny::new(3.0, dims),
                    ],
                    vec![
                        QtyAny::new(4.0, dims),
                        QtyAny::new(5.0, dims),
                        QtyAny::new(6.0, dims),
                    ],
                ],
                TableInterpolation::Linear,
                OutsideDomainPolicy::Refuse,
            )
            .expect("valid vector table"),
        )
        .expect("valid vector-table payload"),
    )
}

fn species_payload(dims: Dims, frame: u32) -> Payload {
    let values = vec![
        SpeciesValue::new(
            SpeciesId::new("CO2").expect("species id"),
            QtyAny::new(0.4, dims),
        ),
        SpeciesValue::new(
            SpeciesId::new("H2O").expect("species id"),
            QtyAny::new(0.6, dims),
        ),
    ];
    Payload::SpeciesBundle(
        SpeciesBundle::new(dimensional_meta(dims, frame), SampleSource::fixed(values))
            .expect("valid species payload"),
    )
}

fn characteristic_payload(frame: u32) -> Payload {
    Payload::CharacteristicState(
        CharacteristicState::new(
            heterogeneous_meta(frame),
            vec![
                CharacteristicComponent::new(
                    id("pressure"),
                    CharacteristicDirection::Incoming,
                    QuantityContract::Dimensions(PRESSURE),
                )
                .expect("pressure characteristic"),
                CharacteristicComponent::new(
                    id("temperature"),
                    CharacteristicDirection::Outgoing,
                    QuantityContract::Dimensions(TEMPERATURE),
                )
                .expect("temperature characteristic"),
            ],
            SampleSource::fixed(vec![
                QtyAny::new(101_325.0, PRESSURE),
                QtyAny::new(293.15, TEMPERATURE),
            ]),
        )
        .expect("valid characteristic payload"),
    )
}

fn typed_bc(physics: Physics, kind: BcKind, payload: Payload) -> BoundaryCondition {
    BoundaryCondition {
        region: format!("{physics:?}-{kind:?}"),
        physics,
        kind,
        value: Some(BcValue::Typed(payload)),
        compatibility: None,
        frame: 0,
    }
}

fn codes(findings: &[Violation]) -> Vec<&'static str> {
    findings.iter().map(|finding| finding.code).collect()
}

#[test]
fn g0_typed_expectation_matrix_is_closed_and_exact() {
    let rows = [
        (
            Physics::Magnetics,
            BcKind::MagneticVectorPotential,
            Expectation::Typed {
                kind: PayloadKind::Vector,
                dims: Some(MAGNETIC_VECTOR_POTENTIAL),
            },
        ),
        (
            Physics::Magnetics,
            BcKind::NormalMagneticFluxDensity,
            Expectation::Typed {
                kind: PayloadKind::Scalar,
                dims: Some(MAGNETIC_FLUX_DENSITY),
            },
        ),
        (
            Physics::Electrics,
            BcKind::ElectricPotential,
            Expectation::Typed {
                kind: PayloadKind::Scalar,
                dims: Some(ELECTRIC_POTENTIAL),
            },
        ),
        (
            Physics::Electrics,
            BcKind::NormalCurrentDensity,
            Expectation::Typed {
                kind: PayloadKind::Scalar,
                dims: Some(CURRENT_DENSITY),
            },
        ),
        (
            Physics::GasExchange,
            BcKind::SpeciesAmountFlux,
            Expectation::Typed {
                kind: PayloadKind::SpeciesBundle,
                dims: Some(SPECIES_AMOUNT_FLUX),
            },
        ),
        (
            Physics::GasExchange,
            BcKind::SpeciesMassFlux,
            Expectation::Typed {
                kind: PayloadKind::SpeciesBundle,
                dims: Some(SPECIES_MASS_FLUX),
            },
        ),
        (
            Physics::GasExchange,
            BcKind::GasCharacteristicInlet,
            Expectation::Typed {
                kind: PayloadKind::CharacteristicState,
                dims: None,
            },
        ),
        (
            Physics::GasExchange,
            BcKind::GasCharacteristicOutlet,
            Expectation::Typed {
                kind: PayloadKind::CharacteristicState,
                dims: None,
            },
        ),
    ];

    for (physics, kind, expected) in rows {
        assert_eq!(expectation(physics, kind), expected, "{physics:?}/{kind:?}");
    }
    assert_eq!(
        expectation(Physics::GasExchange, BcKind::MassFlowInlet),
        Expectation::Unsupported
    );
    assert_eq!(
        expectation(Physics::IncompressibleFlow, BcKind::GasCharacteristicOutlet),
        Expectation::Unsupported
    );
    assert_eq!(
        expectation(Physics::Magnetics, BcKind::ElectricPotential),
        Expectation::Unsupported
    );
}

#[test]
fn g0_matching_payloads_satisfy_every_new_expectation_row() {
    let rows = [
        (
            Physics::Magnetics,
            BcKind::MagneticVectorPotential,
            vector_payload(MAGNETIC_VECTOR_POTENTIAL, 0),
        ),
        (
            Physics::Magnetics,
            BcKind::NormalMagneticFluxDensity,
            scalar_payload(MAGNETIC_FLUX_DENSITY, 0),
        ),
        (
            Physics::Electrics,
            BcKind::ElectricPotential,
            scalar_payload(ELECTRIC_POTENTIAL, 0),
        ),
        (
            Physics::Electrics,
            BcKind::NormalCurrentDensity,
            scalar_payload(CURRENT_DENSITY, 0),
        ),
        (
            Physics::GasExchange,
            BcKind::SpeciesAmountFlux,
            species_payload(SPECIES_AMOUNT_FLUX, 0),
        ),
        (
            Physics::GasExchange,
            BcKind::SpeciesMassFlux,
            species_payload(SPECIES_MASS_FLUX, 0),
        ),
        (
            Physics::GasExchange,
            BcKind::GasCharacteristicInlet,
            characteristic_payload(0),
        ),
        (
            Physics::GasExchange,
            BcKind::GasCharacteristicOutlet,
            characteristic_payload(0),
        ),
    ];

    for (physics, kind, payload) in rows {
        let condition = typed_bc(physics, kind, payload);
        let mut findings = Vec::new();
        condition.check(&mut findings);
        assert!(findings.is_empty(), "{physics:?}/{kind:?}: {findings:#?}");
    }
}

#[test]
fn g0_typed_rows_report_kind_dimension_and_frame_disagreements() {
    let condition = typed_bc(
        Physics::Magnetics,
        BcKind::MagneticVectorPotential,
        scalar_payload(TEMPERATURE, 7),
    );
    let mut findings = Vec::new();
    condition.check(&mut findings);
    let observed = codes(&findings);
    assert!(observed.contains(&"bc-payload-kind"), "{findings:#?}");
    assert!(observed.contains(&"bc-payload-dims"), "{findings:#?}");
    assert!(observed.contains(&"bc-payload-frame"), "{findings:#?}");

    let dims = findings
        .iter()
        .find(|finding| finding.code == "bc-payload-dims")
        .expect("dimension finding");
    assert!(dims.what.contains(&format!("{:?}", TEMPERATURE.0)));
    assert!(
        dims.what
            .contains(&format!("{:?}", MAGNETIC_VECTOR_POTENTIAL.0))
    );
}

#[test]
fn g0_heterogeneous_characteristic_states_are_admitted_without_dimension_coercion() {
    let payload = characteristic_payload(0);
    assert_eq!(payload.homogeneous_dims(), None);
    for kind in [
        BcKind::GasCharacteristicInlet,
        BcKind::GasCharacteristicOutlet,
    ] {
        let condition = typed_bc(Physics::GasExchange, kind, payload.clone());
        let mut findings = Vec::new();
        condition.check(&mut findings);
        assert!(findings.is_empty(), "{kind:?}: {findings:#?}");
    }
}

#[test]
fn g0_legacy_and_typed_carriers_are_never_implicitly_coerced() {
    let mut legacy_on_typed = typed_bc(
        Physics::Electrics,
        BcKind::ElectricPotential,
        scalar_payload(ELECTRIC_POTENTIAL, 0),
    );
    legacy_on_typed.value = Some(BcValue::Uniform(QtyAny::new(12.0, ELECTRIC_POTENTIAL)));
    let mut findings = Vec::new();
    legacy_on_typed.check(&mut findings);
    assert_eq!(codes(&findings), ["bc-typed-payload-required"]);

    let typed_on_legacy = BoundaryCondition {
        region: "heated-wall".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Typed(scalar_payload(TEMPERATURE, 0))),
        compatibility: None,
        frame: 0,
    };
    findings.clear();
    typed_on_legacy.check(&mut findings);
    assert_eq!(codes(&findings), ["bc-legacy-value-required"]);

    let typed_mass_flow = BoundaryCondition {
        region: "inlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Typed(scalar_payload(MASS_FLOW, 0))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    findings.clear();
    typed_mass_flow.check(&mut findings);
    assert_eq!(codes(&findings), ["bc-mass-flow-typed"]);

    let mut missing = typed_bc(
        Physics::Electrics,
        BcKind::ElectricPotential,
        scalar_payload(ELECTRIC_POTENTIAL, 0),
    );
    missing.value = None;
    findings.clear();
    missing.check(&mut findings);
    assert_eq!(codes(&findings), ["bc-value-missing"]);
}

#[test]
fn g0_typed_mass_flow_is_refused_by_direct_and_whole_scenario_evaluation() {
    let typed_mass_flow = BoundaryCondition {
        region: "inlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Typed(scalar_payload(MASS_FLOW, 0))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    assert!(matches!(
        typed_mass_flow.mass_flow_at(0.0),
        Err(ScenarioError::Evaluate { what })
            if what.contains("typed Scalar payload data") && what.contains("kg/s")
    ));

    let pressure_outlet = BoundaryCondition {
        region: "outlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        value: Some(BcValue::Uniform(QtyAny::new(101_325.0, PRESSURE))),
        compatibility: None,
        frame: 0,
    };
    let mut scenario = Scenario::new("typed-mass-flow", 11, Environment::earth_lab());
    scenario.base_bcs.extend([typed_mass_flow, pressure_outlet]);
    let findings = scenario.validate();
    let observed = codes(&findings);
    assert!(observed.contains(&"bc-mass-flow-typed"), "{findings:#?}");
    assert!(observed.contains(&"flux-evaluation"), "{findings:#?}");
}

#[test]
fn g3_gas_boundary_names_never_change_incompressible_compatibility() {
    let mut gas_outlet = typed_bc(
        Physics::GasExchange,
        BcKind::GasCharacteristicOutlet,
        characteristic_payload(0),
    );
    let mut findings = Vec::new();
    gas_outlet.check(&mut findings);
    assert!(!codes(&findings).contains(&"bc-compat-missing"));

    gas_outlet.compatibility = Some(Compat::Incompressible);
    findings.clear();
    gas_outlet.check(&mut findings);
    assert!(codes(&findings).contains(&"bc-compat-forbidden"));
    gas_outlet.compatibility = None;

    let mass_flow = BoundaryCondition {
        region: "legacy-inlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Uniform(QtyAny::new(1.0, MASS_FLOW))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    let mut scenario = Scenario::new("gas-not-pressure", 12, Environment::earth_lab());
    scenario.base_bcs.extend([mass_flow, gas_outlet]);
    let findings = scenario.validate();
    assert!(
        codes(&findings).contains(&"flux-imbalance"),
        "{findings:#?}"
    );

    let mut species_only = Scenario::new("species-not-total", 13, Environment::earth_lab());
    species_only.base_bcs.push(typed_bc(
        Physics::GasExchange,
        BcKind::SpeciesMassFlux,
        species_payload(SPECIES_MASS_FLUX, 0),
    ));
    assert_eq!(species_only.validate(), Vec::new());
}

#[test]
fn g0_legacy_expectation_rows_remain_unchanged() {
    assert_eq!(
        expectation(Physics::IncompressibleFlow, BcKind::Dirichlet),
        Expectation::Value(VELOCITY)
    );
    assert_eq!(
        expectation(Physics::IncompressibleFlow, BcKind::MassFlowInlet),
        Expectation::Value(MASS_FLOW)
    );
    assert_eq!(
        expectation(Physics::IncompressibleFlow, BcKind::PressureOutlet),
        Expectation::Value(PRESSURE)
    );
    assert_eq!(
        expectation(Physics::IncompressibleFlow, BcKind::WallNoSlip),
        Expectation::NoValue
    );
    assert_eq!(
        expectation(Physics::Thermal, BcKind::Dirichlet),
        Expectation::Value(TEMPERATURE)
    );
    assert_eq!(
        expectation(Physics::Elasticity, BcKind::Traction),
        Expectation::Value(PRESSURE)
    );
}

#[test]
fn g0_typed_payload_scenario_ir_round_trips_and_versions_fail_closed() {
    let mut scenario = Scenario::new("typed-ir", 14, Environment::earth_lab());
    scenario.base_bcs.push(typed_bc(
        Physics::Electrics,
        BcKind::ElectricPotential,
        scalar_payload(ELECTRIC_POTENTIAL, 0),
    ));
    let canonical = write_ir(&scenario);
    const MARKER: &str = "(typed :version 1 \"";
    let hex_start = canonical.find(MARKER).expect("typed wire marker") + MARKER.len();
    let hex_end = hex_start
        + canonical[hex_start..]
            .find('"')
            .expect("typed wire terminator");
    assert_eq!(
        parse_ir(&canonical).expect("canonical typed IR").scenario(),
        &scenario
    );

    let wrong_inner_version = canonical.replacen("(typed :version 1 ", "(typed :version 2 ", 1);
    assert!(
        parse_ir(&wrong_inner_version)
            .expect_err("unknown payload version must refuse")
            .to_string()
            .contains("unsupported typed payload version 2")
    );

    let mut uppercase = canonical.clone().into_bytes();
    let lowercase_offset = uppercase[hex_start..hex_end]
        .iter()
        .position(u8::is_ascii_lowercase)
        .expect("payload magic contains a lowercase hex digit");
    uppercase[hex_start + lowercase_offset].make_ascii_uppercase();
    let uppercase = String::from_utf8(uppercase).expect("ASCII mutation remains UTF-8");
    let error = parse_ir(&uppercase).expect_err("uppercase hex must refuse");
    let ScenarioError::Parse { span, path, what } = error else {
        panic!("uppercase payload hex must be a located parse refusal");
    };
    assert_eq!(
        span,
        IrSourceSpan {
            start: hex_start + lowercase_offset,
            end: hex_start + lowercase_offset + 1,
        }
    );
    assert_eq!(path, "$.base_bcs[0].value.typed.bytes");
    assert!(what.contains("canonical lowercase digits"));

    let mut wrong_tag = canonical.clone();
    let typed_start = wrong_tag
        .find("(typed :version ")
        .expect("typed payload header is present");
    let tag_start = typed_start + "(typed ".len();
    wrong_tag.replace_range(tag_start..tag_start + ":version".len(), ":versoin");
    let error = parse_ir(&wrong_tag).expect_err("misspelled typed version tag must refuse");
    let ScenarioError::Parse { span, path, what } = error else {
        panic!("typed version tag must be a located parse refusal");
    };
    assert_eq!(
        span,
        IrSourceSpan {
            start: tag_start,
            end: tag_start + ":versoin".len(),
        }
    );
    assert_eq!(path, "$.base_bcs[0].value.typed.version_tag");
    assert!(what.contains("must be :version"));

    let mut trailing = canonical.clone();
    trailing.insert_str(hex_end, "00");
    assert!(
        parse_ir(&trailing)
            .expect_err("payload trailing byte must refuse")
            .to_string()
            .contains("trailing")
    );

    let form_start = canonical.find("(typed ").expect("typed form");
    let form_end = form_start
        + canonical[form_start..]
            .find(')')
            .expect("typed form terminator")
        + 1;
    let typed_form = &canonical[form_start..form_end];
    let explicit_legacy = format!(
        "(scenario :version 1 \"legacy-typed\" 7 (environment \
         (qty 0 1 0 -2 0 0) (qty 0 1 0 -2 0 0) \
         (qty -9.80665 1 0 -2 0 0) (qty 293.15 0 0 0 1 0) \
         (qty 101325 -1 1 -2 0 0)) (frames) \
         (bcs (bc \"boundary\" thermal dirichlet 0 {typed_form} none)) \
         (cases) (combos) (ensembles) (contacts))"
    );
    let unversioned_legacy = explicit_legacy.replacen("(scenario :version 1 ", "(scenario ", 1);
    for legacy in [&explicit_legacy, &unversioned_legacy] {
        assert!(
            parse_ir(legacy)
                .expect_err("historical IR cannot embed typed payloads")
                .to_string()
                .contains("typed payloads require scenario IR version 2"),
            "legacy form did not reach the typed-payload version gate: {legacy}"
        );
    }
}

#[test]
fn g0_default_parser_accepts_a_canonical_typed_atom_above_the_old_one_mib_ceiling() {
    let values = (0..40_000)
        .map(|index| QtyAny::new(f64::from(index), MAGNETIC_VECTOR_POTENTIAL))
        .collect();
    let payload = Payload::Vector(
        VectorPayload::new(
            dimensional_meta(MAGNETIC_VECTOR_POTENTIAL, 0),
            SampleSource::fixed(values),
        )
        .expect("large admitted vector payload"),
    );
    let mut scenario = Scenario::new("large-typed-atom", 16, Environment::earth_lab());
    scenario.base_bcs.push(typed_bc(
        Physics::Magnetics,
        BcKind::MagneticVectorPotential,
        payload,
    ));

    let canonical = write_ir(&scenario);
    const MARKER: &str = "(typed :version 1 \"";
    let hex_start = canonical.find(MARKER).expect("typed marker") + MARKER.len();
    let hex_len = canonical[hex_start..].find('"').expect("typed terminator");
    assert!(
        hex_len > 1024 * 1024,
        "fixture must cross the former per-atom ceiling: {hex_len}"
    );
    let decoded = parse_ir(&canonical).expect("default parser admits canonical writer output");
    assert_eq!(decoded.scenario(), &scenario);
    let mut findings = Vec::new();
    check_round_trip(&scenario, &mut findings);
    assert!(findings.is_empty(), "{findings:#?}");
}

#[test]
fn g0_semantic_plan_charges_typed_payload_scalars_and_identities() {
    let payload = vector_table_payload(MAGNETIC_VECTOR_POTENTIAL, 0);
    let dynamic_scalars = payload
        .bounded_dynamic_scalar_count()
        .expect("bounded fixture scalar count");
    let (payload_identity_bytes, payload_component_bytes) =
        payload.identity_stats().expect("bounded identity count");
    assert_eq!(dynamic_scalars, 8, "two times plus two three-vectors");

    let mut baseline = Scenario::new("typed-budget", 15, Environment::earth_lab());
    baseline.cases.push(LoadCase {
        name: "magnetic-case".to_string(),
        bcs: Vec::new(),
    });
    let baseline_plan = baseline
        .validation_plan(ValidationBudget::default())
        .expect("baseline plan");

    let base_region = "base";
    let case_region = "case";
    let mut scenario = baseline.clone();
    let mut control = baseline;
    let mut base_condition = typed_bc(
        Physics::Magnetics,
        BcKind::MagneticVectorPotential,
        payload.clone(),
    );
    base_condition.region = base_region.to_string();
    let mut control_base = base_condition.clone();
    control_base.value = Some(BcValue::Uniform(QtyAny::new(
        0.0,
        MAGNETIC_VECTOR_POTENTIAL,
    )));
    scenario.base_bcs.push(base_condition);
    control.base_bcs.push(control_base);
    let mut case_condition = typed_bc(Physics::Magnetics, BcKind::MagneticVectorPotential, payload);
    case_condition.region = case_region.to_string();
    let mut control_case = case_condition.clone();
    control_case.value = Some(BcValue::Uniform(QtyAny::new(
        0.0,
        MAGNETIC_VECTOR_POTENTIAL,
    )));
    scenario.cases[0].bcs.push(case_condition);
    control.cases[0].bcs.push(control_case);
    let plan = scenario
        .validation_plan(ValidationBudget::default())
        .expect("typed plan");
    let control_plan = control
        .validation_plan(ValidationBudget::default())
        .expect("same-shape legacy-carrier control plan");

    assert_eq!(
        plan.signal_scalars - baseline_plan.signal_scalars,
        2 * dynamic_scalars
    );
    assert_eq!(
        plan.identity_bytes - baseline_plan.identity_bytes,
        base_region.len() + case_region.len() + 2 * payload_identity_bytes
    );
    assert_eq!(
        plan.identity_component_bytes,
        baseline_plan
            .identity_component_bytes
            .max(base_region.len())
            .max(case_region.len())
            .max(payload_component_bytes)
    );
    assert_eq!(plan.identity_component_bytes, payload_component_bytes);
    assert_eq!(
        plan.planned_work - control_plan.planned_work,
        2 * (u128::try_from(dynamic_scalars).expect("fixture count fits u128")
            + u128::try_from(payload_identity_bytes).expect("fixture bytes fit u128")
            + 1),
        "each typed payload adds its scalar slots, identity bytes, and one retained basis-id visit"
    );

    let mut short_scalars = ValidationBudget::default();
    short_scalars.max_signal_scalars = plan.signal_scalars - 1;
    assert!(matches!(
        scenario.validation_plan(short_scalars),
        Err(ValidationError::LimitExceeded {
            resource: "signal scalars",
            requested,
            limit,
        }) if requested == plan.signal_scalars && limit + 1 == requested
    ));

    let mut short_identities = ValidationBudget::default();
    short_identities.max_identity_bytes = plan.identity_bytes - 1;
    assert!(matches!(
        scenario.validation_plan(short_identities),
        Err(ValidationError::LimitExceeded {
            resource: "identity bytes",
            requested,
            limit,
        }) if requested == plan.identity_bytes && limit + 1 == requested
    ));

    let mut short_component = ValidationBudget::default();
    short_component.max_identity_component_bytes = plan.identity_component_bytes - 1;
    assert!(matches!(
        scenario.validation_plan(short_component),
        Err(ValidationError::LimitExceeded {
            resource: "identity component bytes",
            requested,
            limit,
        }) if requested == plan.identity_component_bytes && limit + 1 == requested
    ));

    let mut exact_work = ValidationBudget::default();
    exact_work.max_work = plan.planned_work;
    assert!(scenario.validation_plan(exact_work).is_ok());
    exact_work.max_work -= 1;
    assert!(matches!(
        scenario.validation_plan(exact_work),
        Err(ValidationError::WorkExceeded { requested, limit })
            if requested == plan.planned_work && limit + 1 == requested
    ));
}

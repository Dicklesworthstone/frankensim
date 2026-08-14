//! G0/G3 thermal-QoI and eight-term-budget integration battery.

use std::collections::BTreeMap;

use fs_airflow::qoi::{
    DiscretizationReceipt, FanPowerSpec, JunctionRegion, QoiError, SafetyFactorAuthority,
    SurfaceRegion, ThermalOutputAuditError, ThermalQoiCardUse, ThermalQoiDeclarations,
    ThermalQoiKind, ThermalRequirement, extract_thermal_qois,
};
use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_conduction::ResidualClaim;
use fs_conduction::fixtures::unit_cube;
use fs_conduction::{
    ConductionMesh, ConductionReport, ConductionSolution, EnergyBalance, LinearSolveEvidence,
    ProvenanceClass, StopReason,
};
use fs_convection::{CorrelationId, correlation_catalog};
use fs_evidence::uncertainty::{
    BudgetTotal, ENGINEERING_UNCERTAINTY_TERM_COUNT, EngineeringUncertaintyKind, TermValue,
};
use fs_evidence::{ModelCard, NumericalKind, ValidityDomain};
use fs_qty::{Pressure, Temperature, VolumetricFlowRate};
use fs_regime::{
    EnvelopeCoverage, OperatingPoint as RegimeOperatingPoint, OverrideAcknowledgement,
    RegimeAuditCard,
};

fn source(id: &str) -> SourceProvenance {
    SourceProvenance::new("retained synthetic G0 source", id)
}

fn fan_curve() -> FanCurve {
    FanCurve::new(
        "qoi-fixture-fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(160.0)),
            FanPoint::new(VolumetricFlowRate::new(0.04), Pressure::new(130.0)),
            FanPoint::new(VolumetricFlowRate::new(0.08), Pressure::new(70.0)),
            FanPoint::new(VolumetricFlowRate::new(0.12), Pressure::new(0.0)),
        ],
        source("qoi-fan-v1"),
        0.08,
        ToleranceBasis::EngineeringAllowance,
        VolumetricFlowRate::new(0.01),
        (0.7, 1.3),
    )
    .expect("valid fan fixture")
}

fn loss(name: &str, resistance: f64, uncertainty: f64) -> LossElement {
    LossElement::new(
        name,
        LossResistance::new(resistance),
        uncertainty,
        source(&format!("qoi-loss-{name}")),
        ToleranceBasis::EngineeringAllowance,
    )
    .expect("valid loss fixture")
}

fn network() -> EnclosureNetwork {
    let primary = LossNetwork::series(vec![
        LossNetwork::Element(loss("inlet", 40_000.0, 0.10)),
        LossNetwork::Element(loss("heatsink", 30_000.0, 0.12)),
        LossNetwork::Element(loss("outlet", 12_000.0, 0.08)),
    ])
    .expect("series network");
    EnclosureNetwork::new(
        primary,
        LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
    )
}

fn network_with_validated_heatsink_loss() -> EnclosureNetwork {
    let heatsink = loss("heatsink", 30_000.0, 0.12)
        .with_regime_validity(ValidityDomain::unconstrained().with(
            "loss_reynolds",
            2_000.0,
            80_000.0,
        ))
        .expect("explicit finite loss validity");
    let primary = LossNetwork::series(vec![
        LossNetwork::Element(loss("inlet", 40_000.0, 0.10)),
        LossNetwork::Element(heatsink),
        LossNetwork::Element(loss("outlet", 12_000.0, 0.08)),
    ])
    .expect("series network");
    EnclosureNetwork::new(
        primary,
        LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
    )
}

fn solve_fixture_network(network: &EnclosureNetwork) -> fs_airflow::OperatingPoint {
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, 1.0).expect("fan bank");
    solve_operating_point(&fan, network).expect("operating point")
}

fn operating_point() -> fs_airflow::OperatingPoint {
    solve_fixture_network(&network())
}

fn mesh_and_solution() -> (ConductionMesh, ConductionSolution) {
    let (complex, positions) = unit_cube(1);
    let mesh = ConductionMesh::new(complex, positions).expect("unit cube mesh");
    let temperature = vec![300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 360.0];
    let solution = ConductionSolution {
        temperature,
        report: ConductionReport {
            iterations: 2,
            residual_history: vec![1.0, 1.0e-10],
            final_residual: 1.0e-12,
            residual_threshold: 1.0e-10,
            stop_reason: StopReason::ResidualTolerance,
            linear: Vec::new(),
            energy: EnergyBalance {
                source_w: 10.0,
                neumann_out_w: 0.0,
                robin_out_w: 9.999_999_999_999,
                dirichlet_in_w: 0.0,
                closure_w: 1.0e-12,
                scale_w: 10.0,
            },
            material_provenance: ProvenanceClass::MatdbReceipts,
            material_receipts: 3,
            interface_fluxes: Vec::new(),
            robin_fluxes: Vec::new(),
            free_dofs: 8,
            elements: mesh.element_count(),
        },
    };
    (mesh, solution)
}

fn declarations(mesh: &ConductionMesh) -> (JunctionRegion, SurfaceRegion, FanPowerSpec) {
    let junction = JunctionRegion::try_new("package", vec![7, 0, 6]).expect("junction region");
    let surface =
        SurfaceRegion::try_new("case", (0..mesh.boundary().len()).rev().collect::<Vec<_>>())
            .expect("surface region");
    let power = FanPowerSpec::try_new(0.72, 0.04, source("efficiency-v1")).expect("fan efficiency");
    (junction, surface, power)
}

fn extract_fixture_run() -> (fs_airflow::qoi::ThermalQoiSet, fs_airflow::OperatingPoint) {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let qois = extract_fixture_qois(&mesh, &solution, &operating);
    (qois, operating)
}

/// A declared derating policy. The factor is retained authority only: the
/// effective limit handed to [`ThermalRequirement`] is already post-factor.
fn safety_factor(factor: f64) -> SafetyFactorAuthority {
    SafetyFactorAuthority::try_new(factor, source("derating-policy-v1")).expect("safety factor")
}

fn requirement_at(effective_limit_k: f64, factor: f64) -> ThermalRequirement {
    ThermalRequirement::try_new(
        Temperature::new(effective_limit_k),
        safety_factor(factor),
        source("component-datasheet-limit-v1"),
    )
    .expect("requirement")
}

fn extract_fixture_qois(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating: &fs_airflow::OperatingPoint,
) -> fs_airflow::qoi::ThermalQoiSet {
    extract_fixture_qois_with(
        mesh,
        solution,
        operating,
        &requirement_at(380.0, 1.25),
        None,
    )
}

fn extract_fixture_qois_with(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating: &fs_airflow::OperatingPoint,
    requirement: &ThermalRequirement,
    discretization: Option<&DiscretizationReceipt>,
) -> fs_airflow::qoi::ThermalQoiSet {
    let (junction, surface, power) = declarations(mesh);
    extract_thermal_qois(
        mesh,
        solution,
        operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(requirement),
            discretization,
        },
    )
    .expect("QoI extraction")
}

/// Replace every nodal temperature with `a + b*z` from the mesh's own
/// coordinates, so the expected surface mean and spread follow from calculus
/// rather than from a mirror of the production loop.
fn linear_in_z(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    a: f64,
    b: f64,
) -> ConductionSolution {
    let mut linear = solution.clone();
    linear.temperature = mesh
        .positions()
        .iter()
        .map(|position| a + b * position[2])
        .collect();
    linear
}

fn fan_regime_card() -> ModelCard {
    fan_curve().model_card()
}

fn convection_regime_card() -> ModelCard {
    correlation_catalog()
        .into_iter()
        .find(|card| card.id == CorrelationId::DittusBoelter)
        .expect("catalog retains Dittus-Boelter")
        .model
}

fn thermal_regime_point(id: &str, flow_m3_s: f64, reynolds: f64) -> RegimeOperatingPoint {
    RegimeOperatingPoint {
        id: id.to_string(),
        groups: BTreeMap::from([
            ("L_over_Dh".to_string(), 100.0),
            ("Pr".to_string(), 7.0),
            ("Re".to_string(), reynolds),
            ("flow_m3_s".to_string(), flow_m3_s),
            ("speed_ratio".to_string(), 1.0),
        ]),
    }
}

fn thermal_regime_point_with_loss_reynolds(
    id: &str,
    flow_m3_s: f64,
    reynolds: f64,
    loss_reynolds: f64,
) -> RegimeOperatingPoint {
    let mut point = thermal_regime_point(id, flow_m3_s, reynolds);
    point
        .groups
        .insert("loss_reynolds".to_string(), loss_reynolds);
    point
}

fn card_uses(
    qois: &fs_airflow::qoi::ThermalQoiSet,
    model_cards: &[ModelCard],
) -> Vec<ThermalQoiCardUse> {
    qois.budgets()
        .into_iter()
        .map(|budget| ThermalQoiCardUse {
            qoi: budget.qoi().to_string(),
            model_cards: model_cards.iter().map(|card| card.name.clone()).collect(),
            override_acknowledgement: None,
        })
        .collect()
}

fn regime_card_uses(
    qois: &fs_airflow::qoi::ThermalQoiSet,
    model_cards: &[RegimeAuditCard],
) -> Vec<ThermalQoiCardUse> {
    qois.budgets()
        .into_iter()
        .map(|budget| ThermalQoiCardUse {
            qoi: budget.qoi().to_string(),
            model_cards: model_cards.iter().map(|card| card.name.clone()).collect(),
            override_acknowledgement: None,
        })
        .collect()
}

#[test]
fn every_reference_qoi_emits_an_eight_term_budget_without_laundering_unknowns() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (junction, surface, power) = declarations(&mesh);
    let requirement = requirement_at(380.0, 1.25);

    let qois = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("QoI extraction");

    assert_eq!(qois.junction_maximum.vertex, 6, "lowest-index tie wins");
    assert_eq!(qois.junction_maximum.qoi.evidence.value.value(), 360.0);
    assert_eq!(qois.thermal_margin.evidence.value.value(), 20.0);
    assert_eq!(
        qois.junction_maximum.qoi.evidence.numerical.kind,
        NumericalKind::NoClaim,
        "a raw nodal maximum has no DWR enclosure"
    );
    assert!(qois.fan_power.evidence.value.value() > 0.0);
    assert!(
        qois.uniformity
            .mean_temperature
            .evidence
            .value
            .value()
            .is_finite()
    );
    assert!(qois.uniformity.spread.evidence.value.value() > 0.0);

    for budget in qois.budgets() {
        assert_eq!(budget.terms().len(), ENGINEERING_UNCERTAINTY_TERM_COUNT);
        assert!(matches!(
            budget.term(EngineeringUncertaintyKind::ModelForm).value(),
            TermValue::Unknown { .. }
        ));
        assert!(matches!(budget.total(), BudgetTotal::Unknown { .. }));
        let report = budget.render_report();
        assert!(report.contains("model-form"));
        assert!(report.contains("provenance="));
    }
    assert!(qois.all_totals_are_honestly_unknown());
    assert_eq!(qois.junction_maximum.qoi.uncertainty.unit(), "kelvin");
    assert_eq!(qois.pressure_drop.uncertainty.unit(), "pascal");
    assert_eq!(qois.fan_power.uncertainty.unit(), "watt");

    assert!(matches!(
        qois.pressure_drop
            .uncertainty
            .term(EngineeringUncertaintyKind::BoundaryConditions)
            .value(),
        TermValue::IntervalBound { .. }
    ));
    assert!(matches!(
        qois.fan_power
            .uncertainty
            .term(EngineeringUncertaintyKind::Parameters)
            .value(),
        TermValue::IntervalBound { .. }
    ));
}

#[test]
fn every_emitted_qoi_declares_what_its_scalar_means_beyond_its_unit() {
    let (qois, _) = extract_fixture_run();

    assert_eq!(
        qois.junction_maximum.qoi.kind,
        ThermalQoiKind::AbsoluteTemperature
    );
    assert_eq!(
        qois.uniformity.mean_temperature.kind,
        ThermalQoiKind::AbsoluteTemperature
    );
    assert_eq!(
        qois.uniformity.spread.kind,
        ThermalQoiKind::TemperatureDifference
    );
    assert_eq!(
        qois.uniformity.face_mean_standard_deviation.kind,
        ThermalQoiKind::TemperatureDifference
    );
    assert_eq!(
        qois.thermal_margin.kind,
        ThermalQoiKind::TemperatureDifference
    );
    assert_eq!(qois.pressure_drop.kind, ThermalQoiKind::Pressure);
    assert_eq!(qois.fan_power.kind, ThermalQoiKind::Power);

    // This is the ambiguity the discriminant exists to remove: all five
    // temperature QoIs are `Temperature` carrying the unit label `kelvin`, so
    // unit and Rust type together still cannot separate the two absolute ones
    // from the three intervals.
    for budget in [
        &qois.junction_maximum.qoi.uncertainty,
        &qois.uniformity.mean_temperature.uncertainty,
        &qois.uniformity.spread.uncertainty,
        &qois.uniformity.face_mean_standard_deviation.uncertainty,
        &qois.thermal_margin.uncertainty,
    ] {
        assert_eq!(budget.unit(), "kelvin");
    }
}

#[test]
fn a_negative_thermal_margin_is_an_admissible_interval_never_an_absolute_temperature() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    // The fixture's junction maximum is 360 K, so a 350 K effective limit is
    // the design MISSING its requirement -- the case the product exists to
    // report, and the one where absolute/interval confusion is most costly.
    let qois = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(350.0, 1.0),
        None,
    );

    let margin = qois.thermal_margin.evidence.value.value();
    assert_eq!(margin, -10.0, "the fixture must actually violate its limit");
    assert_eq!(
        qois.thermal_margin.kind,
        ThermalQoiKind::TemperatureDifference
    );
    assert!(qois.thermal_margin.kind.admits_negative());
    // The same scalar under the other kind is inadmissible, which is what
    // makes the discriminant load-bearing rather than decorative.
    assert!(!ThermalQoiKind::AbsoluteTemperature.admits_negative());
    assert_eq!(
        ThermalQoiKind::TemperatureDifference.as_str(),
        "temperature-difference"
    );
}

#[test]
fn an_absolute_temperature_below_zero_kelvin_refuses_instead_of_travelling_downstream() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (junction, surface, power) = declarations(&mesh);
    let requirement = requirement_at(380.0, 1.25);

    let below_zero = linear_in_z(&mesh, &solution, -50.0, 0.0);
    let error = extract_thermal_qois(
        &mesh,
        &below_zero,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect_err("a negative absolute temperature is not a physical solve result");
    match error {
        QoiError::InvalidInput { field, detail } => {
            assert_eq!(field, "temperature quantity kind");
            assert!(
                detail.contains("absolute-temperature"),
                "the refusal names the kind it enforced: {detail}"
            );
        }
        other => panic!("expected a quantity-kind refusal, got {other:?}"),
    }

    // Positive control: the same CONSTANT field moved above zero admits, so
    // the refusal keys on the sign and not on the field being uniform.
    let above_zero = linear_in_z(&mesh, &solution, 300.0, 0.0);
    let admitted = extract_fixture_qois(&mesh, &above_zero, &operating);
    assert_eq!(
        admitted.junction_maximum.qoi.evidence.value.value(),
        300.0,
        "the control must reach the same code path and pass it"
    );
    assert_eq!(admitted.uniformity.spread.evidence.value.value(), 0.0);
}

#[test]
fn final_audit_demotes_every_e05_10_qoi_and_rebinds_each_model_budget() {
    let (qois, operating) = extract_fixture_run();
    let nominal_flow = operating.flow.value.value();
    let card = fan_regime_card();
    let mut uses = card_uses(&qois, core::slice::from_ref(&card));
    uses[0].override_acknowledgement = Some(OverrideAcknowledgement {
        actor: "thermal-reviewer".to_string(),
        reason: "retain estimate for redesign triage".to_string(),
    });
    let audited = qois
        .audit_operating_envelope(
            &[card],
            &[
                thermal_regime_point("nominal", nominal_flow, 50_000.0),
                thermal_regime_point("high-flow", 0.13, 50_000.0),
            ],
            &uses,
        )
        .expect("complete card declarations admit the final audit");

    assert_eq!(audited.audit.receipts.len(), 7);
    assert!(audited.audit.receipts.iter().all(|receipt| {
        receipt.coverage == EnvelopeCoverage::Partial
            && receipt.in_domain_points == ["nominal"]
            && receipt.out_of_domain_points == ["high-flow"]
            && receipt.model_cards.len() == 1
            && receipt.model_cards[0].name == "airflow.fan.qoi-fixture-fan"
            && receipt.model_cards[0].version == "1"
            && receipt.violations.len() == 1
            && receipt.violations[0].point == "high-flow"
            && receipt.violations[0].card == "airflow.fan.qoi-fixture-fan"
            && receipt.violations[0].axis == "flow_m3_s"
            && receipt.violations[0].observed == Some(0.13)
            && receipt.violations[0].hi == 0.12
            && receipt.violations[0].distance > 0.0
            && matches!(
                receipt.effective_color,
                fs_evidence::Color::Estimated { dispersion, .. }
                    if dispersion.is_infinite()
            )
    }));
    assert!(audited.audit.receipts.iter().any(|receipt| {
        receipt
            .override_acknowledgement
            .as_ref()
            .is_some_and(|ack| {
                ack.actor == "thermal-reviewer"
                    && matches!(
                        receipt.effective_color,
                        fs_evidence::Color::Estimated { dispersion, .. }
                            if dispersion.is_infinite()
                    )
            })
    }));
    for budget in audited.qois.budgets() {
        let model = budget.term(EngineeringUncertaintyKind::ModelForm);
        assert!(matches!(model.value(), TermValue::Unknown { .. }));
        assert_eq!(model.provenance().role(), "regime-output-audit");
        let receipt = audited
            .audit
            .receipts
            .iter()
            .find(|receipt| receipt.qoi == budget.qoi())
            .expect("matching final receipt");
        assert_eq!(model.provenance().digest(), receipt.content_id());
    }
}

#[test]
fn actual_convection_card_alone_demotes_the_complete_qoi_set() {
    let (qois, operating) = extract_fixture_run();
    let nominal_flow = operating.flow.value.value();
    let fan = fan_regime_card();
    let convection = convection_regime_card();
    let uses = card_uses(&qois, &[fan.clone(), convection.clone()]);
    let audited = qois
        .audit_operating_envelope(
            &[fan, convection],
            &[
                thermal_regime_point("nominal", nominal_flow, 50_000.0),
                thermal_regime_point("low-reynolds", nominal_flow, 1_000.0),
            ],
            &uses,
        )
        .expect("actual fan and convection cards admit a complete audit");

    assert_eq!(audited.audit.receipts.len(), 7);
    for receipt in &audited.audit.receipts {
        assert_eq!(receipt.coverage, EnvelopeCoverage::Partial);
        assert_eq!(receipt.in_domain_points, ["nominal"]);
        assert_eq!(receipt.out_of_domain_points, ["low-reynolds"]);
        assert_eq!(
            receipt
                .model_cards
                .iter()
                .map(|card| card.name.as_str())
                .collect::<Vec<_>>(),
            [
                "airflow.fan.qoi-fixture-fan",
                CorrelationId::DittusBoelter.name(),
            ]
        );
        assert!(matches!(
            receipt.effective_color,
            fs_evidence::Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ));
        assert_eq!(receipt.violations.len(), 1);
        let violation = &receipt.violations[0];
        assert_eq!(violation.point, "low-reynolds");
        assert_eq!(violation.card, CorrelationId::DittusBoelter.name());
        assert_eq!(violation.axis, "Re");
        assert_eq!(violation.observed, Some(1_000.0));
        assert_eq!(violation.lo, 10_000.0);
        assert!(violation.distance > 0.0);

        let budget = audited
            .qois
            .budgets()
            .into_iter()
            .find(|budget| budget.qoi() == receipt.qoi)
            .expect("matching QoI budget");
        let model = budget.term(EngineeringUncertaintyKind::ModelForm);
        assert!(matches!(model.value(), TermValue::Unknown { .. }));
        assert_eq!(model.provenance().digest(), receipt.content_id());
    }
}

#[test]
fn actual_loss_card_alone_demotes_the_complete_qoi_set() {
    assert!(
        network().regime_audit_cards().is_empty(),
        "legacy loss coefficients have no validated-domain authority"
    );
    assert!(matches!(
        loss("unvalidated", 1_000.0, 0.10)
            .with_regime_validity(ValidityDomain::unconstrained()),
        Err(AirflowError::EmptyLossRegimeDomain { element })
            if element == "unvalidated"
    ));

    let network = network_with_validated_heatsink_loss();
    let cards = network.regime_audit_cards();
    assert_eq!(cards.len(), 1);
    assert_eq!(cards[0].name, "airflow.loss.heatsink");
    assert_eq!(cards[0].version, "qoi-loss-heatsink");
    assert_eq!(
        cards[0].validity.bound("loss_reynolds"),
        Some((2_000.0, 80_000.0))
    );

    let operating = solve_fixture_network(&network);
    let nominal_flow = operating.flow.value.value();
    let (mesh, solution) = mesh_and_solution();
    let qois = extract_fixture_qois(&mesh, &solution, &operating);
    let uses = regime_card_uses(&qois, &cards);
    let audited = qois
        .audit_operating_envelope_with_cards(
            &cards,
            &[
                thermal_regime_point_with_loss_reynolds(
                    "nominal",
                    nominal_flow,
                    50_000.0,
                    50_000.0,
                ),
                thermal_regime_point_with_loss_reynolds(
                    "high-loss-reynolds",
                    nominal_flow,
                    50_000.0,
                    90_000.0,
                ),
            ],
            &uses,
        )
        .expect("actual validated loss card admits the complete audit");

    assert_eq!(audited.audit.receipts.len(), 7);
    for receipt in &audited.audit.receipts {
        assert_eq!(receipt.coverage, EnvelopeCoverage::Partial);
        assert_eq!(receipt.in_domain_points, ["nominal"]);
        assert_eq!(receipt.out_of_domain_points, ["high-loss-reynolds"]);
        assert_eq!(receipt.model_cards.len(), 1);
        assert_eq!(receipt.model_cards[0].name, "airflow.loss.heatsink");
        assert_eq!(receipt.model_cards[0].version, "qoi-loss-heatsink");
        assert_eq!(receipt.violations.len(), 1);
        let violation = &receipt.violations[0];
        assert_eq!(violation.point, "high-loss-reynolds");
        assert_eq!(violation.card, "airflow.loss.heatsink");
        assert_eq!(violation.axis, "loss_reynolds");
        assert_eq!(violation.observed, Some(90_000.0));
        assert_eq!(violation.hi, 80_000.0);
        assert!(violation.distance > 0.0);
        assert!(matches!(
            receipt.effective_color,
            fs_evidence::Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ));

        let budget = audited
            .qois
            .budgets()
            .into_iter()
            .find(|budget| budget.qoi() == receipt.qoi)
            .expect("matching QoI budget");
        let model = budget.term(EngineeringUncertaintyKind::ModelForm);
        assert!(matches!(model.value(), TermValue::Unknown { .. }));
        assert_eq!(model.provenance().digest(), receipt.content_id());
    }
}

#[test]
fn owner_neutral_card_path_is_exactly_the_evidence_card_path() {
    let (qois, operating) = extract_fixture_run();
    let nominal_flow = operating.flow.value.value();
    let cards = vec![fan_regime_card(), convection_regime_card()];
    let audit_cards = cards.iter().map(RegimeAuditCard::from).collect::<Vec<_>>();
    let uses = card_uses(&qois, &cards);
    let points = [
        thermal_regime_point("nominal", nominal_flow, 50_000.0),
        thermal_regime_point("low-reynolds", nominal_flow, 1_000.0),
    ];

    let evidence_audit = qois
        .clone()
        .audit_operating_envelope(&cards, &points, &uses)
        .expect("evidence-card audit");
    let owner_neutral_audit = qois
        .audit_operating_envelope_with_cards(&audit_cards, &points, &uses)
        .expect("owner-neutral audit");

    assert_eq!(owner_neutral_audit, evidence_audit);
}

#[test]
fn final_audit_is_exact_in_domain_and_refuses_incomplete_card_use_maps() {
    let (qois, operating) = extract_fixture_run();
    let nominal_flow = operating.flow.value.value();
    let card = fan_regime_card();
    let uses = card_uses(&qois, core::slice::from_ref(&card));
    let baseline = qois.clone();
    let admitted = qois
        .audit_operating_envelope(
            &[card.clone()],
            &[thermal_regime_point("nominal", nominal_flow, 50_000.0)],
            &uses,
        )
        .expect("in-domain final audit");
    assert_eq!(admitted.qois, baseline);
    assert!(
        admitted
            .audit
            .receipts
            .iter()
            .all(|receipt| !receipt.demoted())
    );

    let mut missing = uses.clone();
    missing.pop();
    assert!(matches!(
        baseline.clone().audit_operating_envelope(
            &[card.clone()],
            &[thermal_regime_point("nominal", nominal_flow, 50_000.0)],
            &missing,
        ),
        Err(ThermalOutputAuditError::MissingCardUse { .. })
    ));

    let mut duplicate = uses.clone();
    duplicate.push(uses[0].clone());
    assert!(matches!(
        baseline.clone().audit_operating_envelope(
            &[card.clone()],
            &[thermal_regime_point("nominal", nominal_flow, 50_000.0)],
            &duplicate,
        ),
        Err(ThermalOutputAuditError::DuplicateCardUse { .. })
    ));

    let mut foreign = uses;
    foreign.push(ThermalQoiCardUse {
        qoi: "foreign-qoi".to_string(),
        model_cards: vec![card.name.clone()],
        override_acknowledgement: None,
    });
    assert!(matches!(
        baseline.audit_operating_envelope(
            &[card],
            &[thermal_regime_point("nominal", nominal_flow, 50_000.0)],
            &foreign,
        ),
        Err(ThermalOutputAuditError::UnknownQoi { .. })
    ));
}

#[test]
fn region_order_is_canonical_and_maximum_tie_break_is_stable() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let requirement = ThermalRequirement::try_new(
        Temperature::new(380.0),
        safety_factor(1.25),
        source("limit-v1"),
    )
    .expect("requirement");
    let power = FanPowerSpec::try_new(0.72, 0.04, source("efficiency-v1")).expect("efficiency");
    let ascending =
        SurfaceRegion::try_new("case", (0..mesh.boundary().len()).collect()).expect("ascending");
    let descending = SurfaceRegion::try_new("case", (0..mesh.boundary().len()).rev().collect())
        .expect("descending");
    let first = JunctionRegion::try_new("package", vec![7, 6, 0]).expect("first");
    let second = JunctionRegion::try_new("package", vec![0, 6, 7]).expect("second");

    let a = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &first,
            surface_region: &ascending,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("first extraction");
    let b = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &second,
            surface_region: &descending,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("second extraction");

    assert_eq!(a, b);
    assert_eq!(a.junction_maximum.vertex, 6);
}

#[test]
fn missing_requirement_and_malformed_regions_refuse() {
    let duplicate =
        JunctionRegion::try_new("package", vec![1, 1]).expect_err("duplicate vertices must refuse");
    assert!(matches!(duplicate, QoiError::InvalidInput { .. }));
    assert!(SurfaceRegion::try_new("", vec![0]).is_err());

    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (junction, surface, power) = declarations(&mesh);
    let missing = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: None,
            discretization: None,
        },
    )
    .expect_err("margin cannot invent a requirement");
    assert_eq!(missing, QoiError::MissingRequirement);
}

#[test]
fn widening_an_upstream_operating_envelope_cannot_shrink_qoi_terms() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let mut wider = operating.clone();
    wider.pressure.numerical.lo *= 0.9;
    wider.pressure.numerical.hi *= 1.1;
    wider.flow.numerical.lo *= 0.9;
    wider.flow.numerical.hi *= 1.1;
    let (junction, surface, power) = declarations(&mesh);
    let requirement = ThermalRequirement::try_new(
        Temperature::new(380.0),
        safety_factor(1.25),
        source("limit-v1"),
    )
    .expect("requirement");

    let base = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("base");
    let enlarged = extract_thermal_qois(
        &mesh,
        &solution,
        &wider,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("wider");

    let upper = |value: &TermValue| match value {
        TermValue::IntervalBound { upper, .. } => *upper,
        other => panic!("expected interval term, got {other:?}"),
    };
    assert!(
        upper(
            &enlarged
                .pressure_drop
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value()
        ) >= upper(
            &base
                .pressure_drop
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value()
        )
    );
    assert!(
        upper(
            &enlarged
                .fan_power
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value()
        ) >= upper(
            &base
                .fan_power
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value()
        )
    );
}

#[test]
fn source_changes_rebind_fan_power_and_margin_identities() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (junction, surface, power_a) = declarations(&mesh);
    let power_b = FanPowerSpec::try_new(0.72, 0.04, source("efficiency-v2"))
        .expect("alternate efficiency source");
    let requirement_a = ThermalRequirement::try_new(
        Temperature::new(380.0),
        safety_factor(1.25),
        source("limit-v1"),
    )
    .expect("first requirement");
    let requirement_b = ThermalRequirement::try_new(
        Temperature::new(380.0),
        safety_factor(1.25),
        source("limit-v2"),
    )
    .expect("second requirement");

    let a = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power_a,
            requirement: Some(&requirement_a),
            discretization: None,
        },
    )
    .expect("first");
    let b = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power_b,
            requirement: Some(&requirement_b),
            discretization: None,
        },
    )
    .expect("second");

    assert_eq!(a.fan_power.evidence.value, b.fan_power.evidence.value);
    assert_ne!(
        a.fan_power.uncertainty.content_id(),
        b.fan_power.uncertainty.content_id()
    );
    assert_eq!(
        a.thermal_margin.evidence.value,
        b.thermal_margin.evidence.value
    );
    assert_ne!(
        a.thermal_margin.uncertainty.content_id(),
        b.thermal_margin.uncertainty.content_id()
    );
}

#[test]
fn geometry_changes_rebind_temperature_qoi_identities() {
    let (mesh, solution) = mesh_and_solution();
    let (complex, mut positions) = unit_cube(1);
    for position in &mut positions {
        for coordinate in position {
            *coordinate *= 2.0;
        }
    }
    let scaled_mesh = ConductionMesh::new(complex, positions).expect("scaled unit cube mesh");
    let operating = operating_point();
    let (junction, surface, power) = declarations(&mesh);
    let (scaled_junction, scaled_surface, scaled_power) = declarations(&scaled_mesh);
    let requirement = ThermalRequirement::try_new(
        Temperature::new(380.0),
        safety_factor(1.25),
        source("limit-v1"),
    )
    .expect("requirement");

    let base = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("base geometry");
    let scaled = extract_thermal_qois(
        &scaled_mesh,
        &solution,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &scaled_junction,
            surface_region: &scaled_surface,
            fan_power: &scaled_power,
            requirement: Some(&requirement),
            discretization: None,
        },
    )
    .expect("scaled geometry");

    assert_eq!(
        base.uniformity.mean_temperature.evidence.value,
        scaled.uniformity.mean_temperature.evidence.value,
        "uniform scaling preserves the area-weighted temperature mean"
    );
    assert_ne!(
        base.uniformity.mean_temperature.uncertainty.content_id(),
        scaled.uniformity.mean_temperature.uncertainty.content_id(),
        "the semantic identity must still bind the physical mesh"
    );
}

/// Analytic oracle. For `T = a + b*z` the P1 face mean is exact, so the
/// area-weighted surface mean over the closed unit cube is `a + b/2` by
/// calculus (bottom `a`, top `a+b`, four sides `a+b/2`, equal areas), the
/// surface-vertex spread is exactly `b`, and the junction maximum is
/// `a + b*max(z)` over the declared region. None of these expectations is
/// computed by re-running the production loop.
#[test]
fn a_linear_field_pins_the_surface_mean_spread_and_junction_maximum_analytically() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (a, b) = (300.0_f64, 60.0_f64);
    let linear = linear_in_z(&mesh, &solution, a, b);
    let qois = extract_fixture_qois(&mesh, &linear, &operating);

    let mean = qois.uniformity.mean_temperature.evidence.value.value();
    assert!(
        (mean - (a + 0.5 * b)).abs() < 1.0e-9,
        "closed-cube area-weighted mean of a+b*z must be a+b/2: got {mean}, want {}",
        a + 0.5 * b
    );

    let spread = qois.uniformity.spread.evidence.value.value();
    assert!(
        (spread - b).abs() < 1.0e-9,
        "surface-vertex spread of a+b*z over the unit cube must be b: got {spread}"
    );

    // The declared junction region is {0, 6, 7}; its maximum is a + b*max(z).
    let junction_max_z = [0usize, 6, 7]
        .into_iter()
        .map(|vertex| mesh.positions()[vertex][2])
        .fold(f64::NEG_INFINITY, f64::max);
    let maximum = qois.junction_maximum.qoi.evidence.value.value();
    assert!(
        (maximum - (a + b * junction_max_z)).abs() < 1.0e-9,
        "junction maximum must be a+b*max(z) over the declared region: got {maximum}"
    );

    // A dispersion cannot exceed half the range it is drawn from.
    let deviation = qois
        .uniformity
        .face_mean_standard_deviation
        .evidence
        .value
        .value();
    assert!(
        deviation > 0.0 && deviation <= 0.5 * spread + 1.0e-12,
        "face-mean dispersion {deviation} must be positive and within half the spread {spread}"
    );
}

/// Degenerate analytic case: a constant field has zero spread and zero
/// dispersion exactly, and the margin is the whole effective limit gap.
#[test]
fn a_constant_field_has_exactly_zero_spread_and_zero_dispersion() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let constant = linear_in_z(&mesh, &solution, 330.0, 0.0);
    let qois = extract_fixture_qois(&mesh, &constant, &operating);

    // Bit equality, not an epsilon: exactness is the claim being made, and
    // `to_bits` also distinguishes the -0.0 a sign-losing mutant would produce.
    assert_eq!(
        qois.uniformity
            .mean_temperature
            .evidence
            .value
            .value()
            .to_bits(),
        330.0_f64.to_bits(),
        "the area-weighted mean of a constant field is that constant"
    );
    assert_eq!(
        qois.uniformity.spread.evidence.value.value().to_bits(),
        0.0_f64.to_bits(),
        "a constant field has exactly zero spread"
    );
    assert_eq!(
        qois.uniformity
            .face_mean_standard_deviation
            .evidence
            .value
            .value()
            .to_bits(),
        0.0_f64.to_bits(),
        "a constant field has exactly zero face-mean dispersion"
    );
    assert_eq!(
        qois.thermal_margin.evidence.value.value().to_bits(),
        (380.0_f64 - 330.0_f64).to_bits(),
        "margin is the effective limit minus the junction maximum"
    );
}

/// The margin must be a functional of the junction maximum, not an
/// independently fabricated budget. A declared refinement-ladder receipt is the
/// only route by which the maximum acquires a known Discretization term; the
/// margin must inherit exactly that term value.
///
/// The `None` arm is the non-vacuity guard: it proves the assertion below is
/// sensitive, because without a receipt both budgets really are Unknown and an
/// equality check alone would pass even for a fabricated budget.
#[test]
fn the_margin_budget_inherits_every_junction_maximum_term() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();

    let without = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(380.0, 1.25),
        None,
    );
    assert!(
        matches!(
            without
                .junction_maximum
                .qoi
                .uncertainty
                .term(EngineeringUncertaintyKind::Discretization)
                .value(),
            TermValue::Unknown { .. }
        ),
        "absent a receipt the maximum's discretization term must stay a named Unknown"
    );

    let receipt = DiscretizationReceipt::try_new(2.5, source("refinement-ladder-v1"))
        .expect("declared ladder half-width");
    let with = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(380.0, 1.25),
        Some(&receipt),
    );

    let maximum_term = with
        .junction_maximum
        .qoi
        .uncertainty
        .term(EngineeringUncertaintyKind::Discretization)
        .value()
        .clone();
    assert_eq!(
        maximum_term,
        TermValue::interval(0.0, 2.5).expect("declared interval"),
        "a declared receipt must populate the maximum's discretization term"
    );

    // Every term transfers 1:1 because margin = limit - maximum in kelvin.
    for kind in EngineeringUncertaintyKind::ALL {
        assert_eq!(
            with.thermal_margin.uncertainty.term(kind).value(),
            with.junction_maximum.qoi.uncertainty.term(kind).value(),
            "the margin must inherit the junction maximum's {} term",
            kind.name()
        );
    }

    // ...and admitting a ladder observation must not upgrade the certificate.
    assert_eq!(
        with.junction_maximum.qoi.evidence.numerical.kind,
        NumericalKind::NoClaim,
        "a declared ladder observation is not a numerical certificate"
    );
}

/// "Safety-factor single application": the factor was applied once, upstream,
/// by the declaring authority. Holding the effective limit fixed and changing
/// only the factor must rebind identity while leaving the margin value
/// bit-identical. A mutant that re-applies the factor here dies on the value
/// assertion; a mutant that drops the factor from identity dies on the
/// identity assertion.
#[test]
fn a_safety_factor_rebinds_identity_without_moving_the_margin() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();

    let lenient = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(380.0, 1.0),
        None,
    );
    let strict = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(380.0, 2.0),
        None,
    );

    assert_eq!(
        lenient.thermal_margin.evidence.value.value().to_bits(),
        strict.thermal_margin.evidence.value.value().to_bits(),
        "the effective limit is already post-factor; the factor must never be re-applied"
    );
    assert_ne!(
        lenient.thermal_margin.uncertainty.content_id(),
        strict.thermal_margin.uncertainty.content_id(),
        "the declared factor and its policy must still bind requirement identity"
    );

    // A factor below one, or non-finite, is not an admissible derating.
    assert!(SafetyFactorAuthority::try_new(0.99, source("bad-policy")).is_err());
    assert!(SafetyFactorAuthority::try_new(f64::NAN, source("bad-policy")).is_err());
    assert!(SafetyFactorAuthority::try_new(f64::INFINITY, source("bad-policy")).is_err());
}

/// DONE-WHEN clause: the validity intersection must demonstrably NARROW when an
/// upstream card narrows. The existing battery moves the operating POINT out of
/// a fixed card; this holds the point fixed and shrinks the CARD around it,
/// which is the direction the clause actually names.
#[test]
fn narrowing_one_upstream_card_narrows_the_intersection_and_demotes() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let qois = extract_fixture_qois(&mesh, &solution, &operating);

    // One fixed point, comfortably inside the published fan-card flow bound.
    let point = thermal_regime_point("nominal", 0.10, 20_000.0);

    let published = fan_regime_card();
    let uses = card_uses(&qois, std::slice::from_ref(&published));
    let wide = qois
        .clone()
        .audit_operating_envelope(
            std::slice::from_ref(&published),
            std::slice::from_ref(&point),
            &uses,
        )
        .expect("published card admits the point");
    assert_eq!(
        wide.audit.receipts[0].coverage,
        EnvelopeCoverage::FullyInDomain,
        "non-vacuity: the point must start inside the published domain"
    );
    assert!(
        wide.audit.receipts.iter().all(|receipt| !receipt.demoted()),
        "nothing may demote while the point is in domain"
    );

    // Narrow ONLY the flow bound, around the same unchanged point.
    let mut narrowed = published.clone();
    narrowed.validity = narrowed.validity.with("flow_m3_s", 0.0, 0.05);
    let tight = qois
        .clone()
        .audit_operating_envelope(std::slice::from_ref(&narrowed), &[point], &uses)
        .expect("narrowed card still audits");

    assert_eq!(
        tight.audit.receipts[0].coverage,
        EnvelopeCoverage::FullyOutOfDomain,
        "narrowing the card must push the unchanged point out of the intersection"
    );
    assert!(
        tight
            .audit
            .receipts
            .iter()
            .all(fs_regime::OutputClaimReceipt::demoted),
        "every QoI must demote once its consumed card no longer covers the point"
    );
    for budget in tight.qois.budgets() {
        assert!(
            matches!(
                budget.term(EngineeringUncertaintyKind::ModelForm).value(),
                TermValue::Unknown { .. }
            ),
            "{} must lose its model-form authority when the domain narrows past the point",
            budget.qoi()
        );
    }
}

/// Metamorphic law: shifting the requirement and every nodal temperature by the
/// same offset leaves the margin invariant, because margin is a difference of
/// two temperatures. This kills any mutant that scales rather than subtracts,
/// or that applies a factor to one side only.
#[test]
fn a_common_temperature_offset_leaves_the_margin_invariant() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();

    let base = extract_fixture_qois_with(
        &mesh,
        &solution,
        &operating,
        &requirement_at(380.0, 1.25),
        None,
    );

    let offset = 17.5_f64;
    let mut shifted_solution = solution.clone();
    for value in &mut shifted_solution.temperature {
        *value += offset;
    }
    let shifted = extract_fixture_qois_with(
        &mesh,
        &shifted_solution,
        &operating,
        &requirement_at(380.0 + offset, 1.25),
        None,
    );

    assert!(
        (base.thermal_margin.evidence.value.value()
            - shifted.thermal_margin.evidence.value.value())
        .abs()
            < 1.0e-9,
        "a common offset on limit and field must leave the margin invariant"
    );
    assert!(
        (base.uniformity.spread.evidence.value.value()
            - shifted.uniformity.spread.evidence.value.value())
        .abs()
            < 1.0e-9,
        "a common offset must leave the surface spread invariant"
    );
    assert!(
        (shifted.uniformity.mean_temperature.evidence.value.value()
            - base.uniformity.mean_temperature.evidence.value.value()
            - offset)
            .abs()
            < 1.0e-9,
        "a common offset must shift the mean by exactly that offset"
    );
}

/// Upstream-evidence mutation sweep. Each case degrades exactly one upstream
/// evidence slice and asserts the producer refuses or demotes with a stable
/// diagnostic, rather than emitting a QoI that reads like the undegraded one.
#[test]
fn degrading_any_upstream_evidence_slice_refuses_or_demotes() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();

    // Baseline is admissible, so every refusal below is caused by the single
    // mutation and not by the fixture.
    let base = extract_fixture_qois(&mesh, &solution, &operating);

    // (1) Material authority: a caller-declared conductivity has NO material
    // provenance. It must not read like a receipt-backed solve.
    let mut declared = solution.clone();
    declared.report.material_provenance = ProvenanceClass::Declared;
    declared.report.material_receipts = 0;
    let declared_qois = extract_fixture_qois(&mesh, &declared, &operating);

    assert!(
        base.junction_maximum
            .qoi
            .evidence
            .model
            .cards
            .iter()
            .any(|card| card == "fs-conduction:material-matdb-receipts"),
        "a receipt-backed solve must name its material authority"
    );
    assert!(
        declared_qois
            .junction_maximum
            .qoi
            .evidence
            .model
            .cards
            .iter()
            .any(|card| card == "fs-conduction:material-declared"),
        "a declared-conductivity solve must name that it has no material receipt"
    );
    assert_ne!(
        base.junction_maximum.qoi.evidence.model.cards,
        declared_qois.junction_maximum.qoi.evidence.model.cards,
        "material authority must be visible in the evidence, not only in the digest"
    );

    // The Parameters gap must name WHICH gap it is, in every temperature QoI.
    for (label, base_qoi, declared_qoi) in [
        (
            "junction maximum",
            &base.junction_maximum.qoi,
            &declared_qois.junction_maximum.qoi,
        ),
        (
            "surface mean",
            &base.uniformity.mean_temperature,
            &declared_qois.uniformity.mean_temperature,
        ),
        (
            "thermal margin",
            &base.thermal_margin,
            &declared_qois.thermal_margin,
        ),
    ] {
        let reason = |qoi: &fs_airflow::qoi::ThermalQoi<Temperature>| match qoi
            .uncertainty
            .term(EngineeringUncertaintyKind::Parameters)
            .value()
        {
            TermValue::Unknown { reason } => reason.clone(),
            other => panic!("{label} parameters term must stay Unknown, got {other:?}"),
        };
        assert_ne!(
            reason(base_qoi),
            reason(declared_qoi),
            "{label} must distinguish declared material authority from receipt-backed"
        );
        assert!(
            reason(declared_qoi).contains("caller-declared"),
            "{label} must name the declared-authority gap, got: {}",
            reason(declared_qoi)
        );
    }
    // The margin reaches the correct reason only by inheriting it, which is
    // the propagation seam doing real work rather than being decorative.
}

/// The second half of the upstream-evidence sweep: a retained report that is
/// self-contradictory or algebraically unconverged must refuse, not produce a
/// QoI that reads like a supported one.
#[test]
fn a_contradictory_or_unconverged_conduction_report_refuses() {
    let (mesh, solution) = mesh_and_solution();
    let operating = operating_point();
    let (junction, surface, power) = declarations(&mesh);

    // (2) A claim of retained receipts with zero receipts is self-contradictory.
    let mut contradictory = solution.clone();
    contradictory.report.material_receipts = 0;
    let refusal = extract_thermal_qois(
        &mesh,
        &contradictory,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement_at(380.0, 1.25)),
            discretization: None,
        },
    )
    .expect_err("receipts claimed with none retained must refuse");
    match &refusal {
        QoiError::InvalidInput { field, detail } => {
            assert_eq!(*field, "conduction report");
            assert!(
                detail.contains("zero"),
                "diagnostic must name the contradiction, got: {detail}"
            );
        }
        other => panic!("expected a typed refusal, got {other:?}"),
    }

    // (3) An algebraically unconverged linear solve is unsupported, not weaker.
    let mut unconverged = solution.clone();
    unconverged.report.linear.push(LinearSolveEvidence {
        nonlinear_iteration: 1,
        method: "pcg",
        iterations: 500,
        reported: ResidualClaim::RecursiveEstimate(1.0e-11),
        true_relative_residual: 3.2e-2,
        converged_true: false,
        stall: None,
    });
    let refusal = extract_thermal_qois(
        &mesh,
        &unconverged,
        &operating,
        &ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &power,
            requirement: Some(&requirement_at(380.0, 1.25)),
            discretization: None,
        },
    )
    .expect_err("a non-converged linear solve must refuse");
    match &refusal {
        QoiError::InvalidInput { field, detail } => {
            assert_eq!(*field, "conduction report");
            assert!(
                detail.contains("without converging") && detail.contains("pcg"),
                "diagnostic must name the failing solve, got: {detail}"
            );
        }
        other => panic!("expected a typed refusal, got {other:?}"),
    }

    // A converged record with the same shape must still be admissible, which
    // proves the guard keys on convergence and not merely on the vec length.
    let mut converged = solution.clone();
    converged.report.linear.push(LinearSolveEvidence {
        nonlinear_iteration: 1,
        method: "pcg",
        iterations: 12,
        reported: ResidualClaim::RecursiveEstimate(1.0e-11),
        true_relative_residual: 1.0e-11,
        converged_true: true,
        stall: None,
    });
    extract_fixture_qois(&mesh, &converged, &operating);
}

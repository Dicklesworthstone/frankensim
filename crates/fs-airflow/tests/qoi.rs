//! G0/G3 thermal-QoI and eight-term-budget integration battery.

use std::collections::BTreeMap;

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, QoiError, SurfaceRegion, ThermalOutputAuditError,
    ThermalQoiCardUse, ThermalRequirement, extract_thermal_qois,
};
use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::{
    ConductionMesh, ConductionReport, ConductionSolution, EnergyBalance, ProvenanceClass,
    StopReason,
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

fn extract_fixture_qois(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating: &fs_airflow::OperatingPoint,
) -> fs_airflow::qoi::ThermalQoiSet {
    let (junction, surface, power) = declarations(mesh);
    let requirement = ThermalRequirement::try_new(
        Temperature::new(380.0),
        source("component-datasheet-limit-v1"),
    )
    .expect("requirement");
    extract_thermal_qois(
        mesh,
        solution,
        operating,
        &junction,
        &surface,
        &power,
        Some(&requirement),
    )
    .expect("QoI extraction")
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
    let requirement = ThermalRequirement::try_new(
        Temperature::new(380.0),
        source("component-datasheet-limit-v1"),
    )
    .expect("requirement");

    let qois = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &junction,
        &surface,
        &power,
        Some(&requirement),
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
    let requirement = ThermalRequirement::try_new(Temperature::new(380.0), source("limit-v1"))
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
        &first,
        &ascending,
        &power,
        Some(&requirement),
    )
    .expect("first extraction");
    let b = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &second,
        &descending,
        &power,
        Some(&requirement),
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
        &mesh, &solution, &operating, &junction, &surface, &power, None,
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
    let requirement = ThermalRequirement::try_new(Temperature::new(380.0), source("limit-v1"))
        .expect("requirement");

    let base = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &junction,
        &surface,
        &power,
        Some(&requirement),
    )
    .expect("base");
    let enlarged = extract_thermal_qois(
        &mesh,
        &solution,
        &wider,
        &junction,
        &surface,
        &power,
        Some(&requirement),
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
    let requirement_a = ThermalRequirement::try_new(Temperature::new(380.0), source("limit-v1"))
        .expect("first requirement");
    let requirement_b = ThermalRequirement::try_new(Temperature::new(380.0), source("limit-v2"))
        .expect("second requirement");

    let a = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &junction,
        &surface,
        &power_a,
        Some(&requirement_a),
    )
    .expect("first");
    let b = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &junction,
        &surface,
        &power_b,
        Some(&requirement_b),
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
    let requirement = ThermalRequirement::try_new(Temperature::new(380.0), source("limit-v1"))
        .expect("requirement");

    let base = extract_thermal_qois(
        &mesh,
        &solution,
        &operating,
        &junction,
        &surface,
        &power,
        Some(&requirement),
    )
    .expect("base geometry");
    let scaled = extract_thermal_qois(
        &scaled_mesh,
        &solution,
        &operating,
        &scaled_junction,
        &scaled_surface,
        &scaled_power,
        Some(&requirement),
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

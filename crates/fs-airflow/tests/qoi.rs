//! G0/G3 thermal-QoI and eight-term-budget integration battery.

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, QoiError, SurfaceRegion, ThermalRequirement, extract_thermal_qois,
};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, SourceProvenance, ToleranceBasis, solve_operating_point,
};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::{
    ConductionMesh, ConductionReport, ConductionSolution, EnergyBalance, ProvenanceClass,
    StopReason,
};
use fs_evidence::NumericalKind;
use fs_evidence::uncertainty::{
    BudgetTotal, ENGINEERING_UNCERTAINTY_TERM_COUNT, EngineeringUncertaintyKind, TermValue,
};
use fs_qty::{Pressure, Temperature, VolumetricFlowRate};

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

fn operating_point() -> fs_airflow::OperatingPoint {
    let primary = LossNetwork::series(vec![
        LossNetwork::Element(loss("inlet", 40_000.0, 0.10)),
        LossNetwork::Element(loss("heatsink", 30_000.0, 0.12)),
        LossNetwork::Element(loss("outlet", 12_000.0, 0.08)),
    ])
    .expect("series network");
    let network = EnclosureNetwork::new(
        primary,
        LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
    );
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, 1.0).expect("fan bank");
    solve_operating_point(&fan, &network).expect("operating point")
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
            budget.term(EngineeringUncertaintyKind::ModelForm).value,
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
            .value,
        TermValue::IntervalBound { .. }
    ));
    assert!(matches!(
        qois.fan_power
            .uncertainty
            .term(EngineeringUncertaintyKind::Parameters)
            .value,
        TermValue::IntervalBound { .. }
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
                .value
        ) >= upper(
            &base
                .pressure_drop
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value
        )
    );
    assert!(
        upper(
            &enlarged
                .fan_power
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value
        ) >= upper(
            &base
                .fan_power
                .uncertainty
                .term(EngineeringUncertaintyKind::BoundaryConditions)
                .value
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

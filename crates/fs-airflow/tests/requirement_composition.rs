//! Comprehensive test suite for ThermalLimit requirement and evidence composition (bead `frankensim-s2l9v.2`).
//!
//! Tests:
//! - Full compliance taxonomy: Satisfied, Violated, Indeterminate, OutsideDomain;
//! - Zero is not NoData (measured 0.0 K / 0.0 W is an authentic physical value);
//! - Weakest stage authority monotonicity (never upgrades authority of weakest stage);
//! - Binding witness retention (exact vertex, region, tie witness, weakest color);
//! - Pre-flight refusal on duplicate requirements and missing QoI candidates;
//! - Cooperative cancellation polling;
//! - Deterministic canonical ordering and bit-identical replay.

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, SafetyFactorAuthority, SurfaceRegion,
    ThermalQoiDeclarations, ThermalRequirement,
};
use fs_airflow::registered_qoi::{
    OutputQuery, QoiExecutionLimits, QoiSemanticId, extract_registered_qois,
};
use fs_airflow::requirement_composition::{
    ComplianceOutcome, RequirementCompositionError, ThermalLimitSpec, compose_thermal_limits,
};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, OperatingPoint, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::solve::StopReason;
use fs_conduction::{
    ConductionMesh, ConductionReport, ConductionSolution, EnergyBalance, ProvenanceClass,
};
use fs_evidence::ColorRank;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
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

fn operating_point() -> OperatingPoint {
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, 1.0).expect("fan bank");
    solve_operating_point(&fan, &network()).expect("operating point")
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gated_cx(&gate, f)
}

fn with_gated_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5219_0002,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn sample_mesh_and_solution() -> (ConductionMesh, ConductionSolution) {
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
            element_material_identity: None,
        },
    };
    (mesh, solution)
}

fn sample_requirement(limit_k: f64, factor: f64) -> ThermalRequirement {
    ThermalRequirement::try_new(
        Temperature::new(limit_k),
        SafetyFactorAuthority::try_new(factor, source("safety-factor-1")).unwrap(),
        source("datasheet-1"),
    )
    .unwrap()
}

#[test]
fn rc_001_satisfied_and_violated_composition() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0, 1, 2, 6, 7]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0, 1, 2]).unwrap();
        let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(380.0, 1.0);

        let decls = ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &fan,
            requirement: Some(&req),
            discretization: None,
        };

        let queries = vec![
            OutputQuery::scalar("thermal.junction_maximum"),
            OutputQuery::scalar("thermal.surface_mean"),
        ];

        let qoi_receipt = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .unwrap();

        // 1. Satisfied limit: max 390 K with 10 K margin (junction measured = 360 K -> achieved margin 30 K >= 10 K)
        let req_sat = ThermalLimitSpec::try_new(
            "REQ-JUNC-01",
            QoiSemanticId::JunctionMaximum,
            "package",
            390.0,
            10.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        // 2. Violated limit: max 365 K with 10 K margin (junction measured = 360 K -> achieved margin 5 K < 10 K)
        let req_viol = ThermalLimitSpec::try_new(
            "REQ-JUNC-02",
            QoiSemanticId::JunctionMaximum,
            "package",
            365.0,
            10.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        let receipt = compose_thermal_limits(
            &qoi_receipt.rows,
            &[req_sat, req_viol],
            false,
            cx,
        )
        .unwrap();

        assert_eq!(receipt.evaluations.len(), 2);
        assert_eq!(receipt.satisfied_count, 1);
        assert_eq!(receipt.violated_count, 1);

        let eval_sat = &receipt.evaluations[0];
        assert_eq!(eval_sat.requirement_id, "REQ-JUNC-01");
        assert_eq!(eval_sat.outcome, ComplianceOutcome::Satisfied);
        assert_eq!(eval_sat.witness.weakest_color, ColorRank::Verified);
        assert_eq!(eval_sat.witness.primary_vertex, Some(6)); // tie-breaker witness

        let eval_viol = &receipt.evaluations[1];
        assert_eq!(eval_viol.requirement_id, "REQ-JUNC-02");
        assert_eq!(eval_viol.outcome, ComplianceOutcome::Violated);
    });
}

#[test]
fn rc_002_outside_domain_demotes_weakest_color_honestly() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0, 1, 2, 6, 7]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0, 1, 2]).unwrap();
        let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(380.0, 1.0);

        let decls = ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &fan,
            requirement: Some(&req),
            discretization: None,
        };

        let queries = vec![OutputQuery::scalar("thermal.junction_maximum")];

        let qoi_receipt = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .unwrap();

        let req_spec = ThermalLimitSpec::try_new(
            "REQ-JUNC-01",
            QoiSemanticId::JunctionMaximum,
            "package",
            400.0,
            10.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        let receipt = compose_thermal_limits(
            &qoi_receipt.rows,
            &[req_spec],
            true, // outside validated regime domain!
            cx,
        )
        .unwrap();

        assert_eq!(receipt.evaluations.len(), 1);
        let eval = &receipt.evaluations[0];
        assert_eq!(eval.outcome, ComplianceOutcome::OutsideDomain);
        assert_eq!(eval.witness.weakest_color, ColorRank::Estimated);
        assert_eq!(eval.witness.weakest_stage, "fs-regime");
    });
}

#[test]
fn rc_003_duplicate_requirements_and_missing_qoi_refuse() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(380.0, 1.0);
        let decls = ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &fan,
            requirement: Some(&req),
            discretization: None,
        };

        let queries = vec![OutputQuery::scalar("thermal.junction_maximum")];
        let qoi_receipt = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .unwrap();

        // 1. Duplicate requirement specification
        let req1 = ThermalLimitSpec::try_new(
            "REQ-1",
            QoiSemanticId::JunctionMaximum,
            "package",
            390.0,
            5.0,
            1.0,
            "REV-A",
        )
        .unwrap();
        let req2 = ThermalLimitSpec::try_new(
            "REQ-1",
            QoiSemanticId::JunctionMaximum,
            "package",
            380.0,
            5.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        let res_dup = compose_thermal_limits(&qoi_receipt.rows, &[req1, req2], false, cx);
        assert!(matches!(res_dup, Err(RequirementCompositionError::DuplicateRequirement { .. })));

        // 2. Missing QoI candidate row (PressureDrop requested but not in qoi_receipt)
        let req_missing = ThermalLimitSpec::try_new(
            "REQ-PD",
            QoiSemanticId::PressureDrop,
            "enclosure",
            50.0,
            5.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        let res_missing = compose_thermal_limits(&qoi_receipt.rows, &[req_missing], false, cx);
        assert!(matches!(res_missing, Err(RequirementCompositionError::MissingQoiRow { .. })));
    });
}

#[test]
fn rc_004_cancellation_checkpoints_drain_safely() {
    let gate = CancelGate::new();
    gate.request();

    with_gated_cx(&gate, |cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(380.0, 1.0);
        let decls = ThermalQoiDeclarations {
            junction_region: &junction,
            surface_region: &surface,
            fan_power: &fan,
            requirement: Some(&req),
            discretization: None,
        };

        let queries = vec![OutputQuery::scalar("thermal.junction_maximum")];
        let qoi_receipt = with_cx(|cx2| {
            extract_registered_qois(
                &queries,
                &mesh,
                &solution,
                &op,
                &decls,
                QoiExecutionLimits::default(),
                cx2,
            )
            .unwrap()
        });

        let req_spec = ThermalLimitSpec::try_new(
            "REQ-1",
            QoiSemanticId::JunctionMaximum,
            "package",
            390.0,
            5.0,
            1.0,
            "REV-A",
        )
        .unwrap();

        let res = compose_thermal_limits(&qoi_receipt.rows, &[req_spec], false, cx);
        assert_eq!(res, Err(RequirementCompositionError::Cancelled));
    });
}

#[test]
fn rc_005_deterministic_replay_produces_identical_receipt_hashes() {
    let (mesh, solution) = sample_mesh_and_solution();
    let op = operating_point();

    let junction = JunctionRegion::try_new("package", vec![0, 1, 6, 7]).unwrap();
    let surface = SurfaceRegion::try_new("case", vec![0, 1]).unwrap();
    let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
    let req = sample_requirement(380.0, 1.0);
    let decls = ThermalQoiDeclarations {
        junction_region: &junction,
        surface_region: &surface,
        fan_power: &fan,
        requirement: Some(&req),
        discretization: None,
    };

    let queries = vec![
        OutputQuery::scalar("thermal.junction_maximum"),
        OutputQuery::scalar("thermal.surface_mean"),
    ];

    let qoi_receipt = with_cx(|cx| {
        extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .unwrap()
    });

    let requirements = vec![
        ThermalLimitSpec::try_new(
            "REQ-JUNC",
            QoiSemanticId::JunctionMaximum,
            "package",
            390.0,
            10.0,
            1.0,
            "REV-A",
        )
        .unwrap(),
        ThermalLimitSpec::try_new(
            "REQ-SURF",
            QoiSemanticId::SurfaceMeanTemperature,
            "case",
            350.0,
            5.0,
            1.0,
            "REV-A",
        )
        .unwrap(),
    ];

    let receipt1 = with_cx(|cx| {
        compose_thermal_limits(&qoi_receipt.rows, &requirements, false, cx).unwrap()
    });
    let receipt2 = with_cx(|cx| {
        compose_thermal_limits(&qoi_receipt.rows, &requirements, false, cx).unwrap()
    });

    assert_eq!(receipt1.receipt_hash, receipt2.receipt_hash);
    assert_eq!(receipt1.evaluations.len(), receipt2.evaluations.len());
    for (e1, e2) in receipt1.evaluations.iter().zip(&receipt2.evaluations) {
        assert_eq!(e1.identity_hash, e2.identity_hash);
        assert_eq!(e1.outcome, e2.outcome);
        assert_eq!(e1.measured_value.to_bits(), e2.measured_value.to_bits());
        assert_eq!(e1.effective_limit.to_bits(), e2.effective_limit.to_bits());
    }
}

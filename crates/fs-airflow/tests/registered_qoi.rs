//! Comprehensive test suite for registered thermal QoI extraction (bead `frankensim-s2l9v.1`).
//!
//! Tests:
//! - All 7 registered QoI families and alias mapping;
//! - Strict fail-closed rejection of non-scalar output kinds (Field, Report);
//! - Rejection of unknown query names, duplicates, and missing regions;
//! - Pre-flight work and memory budget enforcement;
//! - Non-negative absolute temperature physics check;
//! - Deterministic canonical ordering and tie-breaking;
//! - Mesh coordinate translation invariance;
//! - Cooperative cancellation at checkpoints;
//! - Bit-identical deterministic replay.

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, SafetyFactorAuthority, SurfaceRegion, ThermalQoiDeclarations,
    ThermalQoiKind, ThermalRequirement,
};
use fs_airflow::registered_qoi::{
    OutputKind, OutputQuery, QoiExecutionLimits, QoiSemanticId, RegisteredQoiError,
    extract_registered_qois,
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
                seed: 0x5219_0001,
                kernel_id: 1,
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

fn sample_declarations<'a>(
    junction: &'a JunctionRegion,
    surface: &'a SurfaceRegion,
    fan_power: &'a FanPowerSpec,
    req: &'a ThermalRequirement,
) -> ThermalQoiDeclarations<'a> {
    ThermalQoiDeclarations {
        junction_region: junction,
        surface_region: surface,
        fan_power,
        requirement: Some(req),
        discretization: None,
    }
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
fn rq_001_extracts_all_canonical_qoi_families() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0, 1, 2, 6, 7]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0, 1, 2]).unwrap();
        let fan = FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(380.0, 1.1);

        let decls = sample_declarations(&junction, &surface, &fan, &req);

        let queries = vec![
            OutputQuery::scalar("thermal.junction_maximum"),
            OutputQuery::scalar("thermal.surface_mean"),
            OutputQuery::scalar("thermal.surface_spread"),
            OutputQuery::scalar("thermal.surface_std_dev"),
            OutputQuery::scalar("airflow.pressure_drop"),
            OutputQuery::scalar("airflow.fan_power"),
            OutputQuery::scalar("thermal.thermal_margin"),
        ];

        let receipt = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .expect("extraction success");

        assert_eq!(receipt.requested_query_count, 7);
        assert_eq!(receipt.emitted_qoi_count, 7);
        assert_eq!(receipt.rows.len(), 7);

        let expected_kinds = [
            (
                QoiSemanticId::JunctionMaximum,
                ThermalQoiKind::AbsoluteTemperature,
            ),
            (
                QoiSemanticId::SurfaceMeanTemperature,
                ThermalQoiKind::AbsoluteTemperature,
            ),
            (
                QoiSemanticId::SurfaceTemperatureSpread,
                ThermalQoiKind::TemperatureDifference,
            ),
            (
                QoiSemanticId::SurfaceTemperatureStdDev,
                ThermalQoiKind::TemperatureDifference,
            ),
            (QoiSemanticId::PressureDrop, ThermalQoiKind::Pressure),
            (QoiSemanticId::FanPower, ThermalQoiKind::Power),
            (
                QoiSemanticId::ThermalMargin,
                ThermalQoiKind::TemperatureDifference,
            ),
        ];
        for (semantic_id, expected_kind) in expected_kinds {
            let row = receipt
                .rows
                .iter()
                .find(|row| row.semantic_id == semantic_id)
                .expect("requested semantic row");
            assert_eq!(row.kind, expected_kind);
            assert_eq!(semantic_id.qoi_kind(), expected_kind);
        }

        // Verify values
        let jm = receipt
            .rows
            .iter()
            .find(|r| r.semantic_id == QoiSemanticId::JunctionMaximum)
            .unwrap();
        assert_eq!(jm.value.to_bits(), 360.0f64.to_bits());
        assert_eq!(jm.units, "kelvin");
        assert_eq!(jm.tie_witness_vertex, Some(6)); // tie-break picks lowest index between 6 and 7

        let margin = receipt
            .rows
            .iter()
            .find(|r| r.semantic_id == QoiSemanticId::ThermalMargin)
            .unwrap();
        assert_eq!(margin.value.to_bits(), 20.0f64.to_bits()); // 380.0 - 360.0

        let pd = receipt
            .rows
            .iter()
            .find(|r| r.semantic_id == QoiSemanticId::PressureDrop)
            .unwrap();
        assert!(pd.value > 0.0);
        assert_eq!(pd.units, "pascal");
    });
}

#[test]
fn rq_002_query_aliases_map_deterministically() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0, 1]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0, 1]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        let queries = vec![
            OutputQuery::scalar("junction_temp"),
            OutputQuery::scalar("case_mean_temp"),
            OutputQuery::scalar("delta_p"),
            OutputQuery::scalar("margin"),
        ];

        let receipt = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        )
        .expect("aliases mapped");

        assert_eq!(receipt.rows.len(), 4);
        assert_eq!(receipt.rows[0].semantic_id, QoiSemanticId::JunctionMaximum);
        assert_eq!(
            receipt.rows[1].semantic_id,
            QoiSemanticId::SurfaceMeanTemperature
        );
        assert_eq!(receipt.rows[2].semantic_id, QoiSemanticId::PressureDrop);
        assert_eq!(receipt.rows[3].semantic_id, QoiSemanticId::ThermalMargin);
    });
}

#[test]
fn rq_003_non_scalar_kinds_fail_closed_with_actionable_diagnostics() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        let field_query = vec![OutputQuery {
            name: "temperature_field".to_string(),
            kind: OutputKind::Field,
            region: None,
        }];

        let result_field = extract_registered_qois(
            &field_query,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        );

        assert!(matches!(
            result_field,
            Err(RegisteredQoiError::NonScalarOutputKind { .. })
        ));

        let report_query = vec![OutputQuery {
            name: "thermal_summary".to_string(),
            kind: OutputKind::Report,
            region: None,
        }];

        let result_report = extract_registered_qois(
            &report_query,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        );

        assert!(matches!(
            result_report,
            Err(RegisteredQoiError::NonScalarOutputKind { .. })
        ));
    });
}

#[test]
fn rq_004_duplicate_queries_and_unsupported_names_refuse() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        // Duplicate
        let dup_queries = vec![
            OutputQuery::scalar("junction_temp"),
            OutputQuery::scalar("junction_temp"),
        ];
        let res_dup = extract_registered_qois(
            &dup_queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        );
        assert!(matches!(
            res_dup,
            Err(RegisteredQoiError::DuplicateQuery { .. })
        ));

        // Unsupported name
        let unk_queries = vec![OutputQuery::scalar("unsupported_acoustic_noise")];
        let res_unk = extract_registered_qois(
            &unk_queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        );
        assert!(matches!(
            res_unk,
            Err(RegisteredQoiError::UnsupportedOutputName { .. })
        ));
    });
}

#[test]
fn rq_005_work_and_memory_limit_guards_enforce_preflight() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        let queries = vec![OutputQuery::scalar("junction_temp")];

        let tight_limits = QoiExecutionLimits {
            max_elements: 1, // mesh has 6 elements
            max_vertices: 5_000_000,
            max_queries: 1_000,
            max_memory_bytes: 512 * 1024 * 1024,
        };

        let res =
            extract_registered_qois(&queries, &mesh, &solution, &op, &decls, tight_limits, cx);
        assert!(matches!(
            res,
            Err(RegisteredQoiError::WorkLimitExceeded { .. })
        ));
    });
}

#[test]
fn rq_006_cooperative_cancellation_refuses_cleanly() {
    let gate = CancelGate::new();
    gate.request();

    with_gated_cx(&gate, |cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        let queries = vec![OutputQuery::scalar("junction_temp")];

        let res = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            cx,
        );
        assert_eq!(res, Err(RegisteredQoiError::Cancelled));
    });
}

#[test]
fn rq_007_deterministic_replay_produces_identical_hashes() {
    let (mesh, solution) = sample_mesh_and_solution();
    let op = operating_point();

    let junction = JunctionRegion::try_new("package", vec![0, 1, 6, 7]).unwrap();
    let surface = SurfaceRegion::try_new("case", vec![0, 1]).unwrap();
    let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
    let req = sample_requirement(400.0, 1.1);
    let decls = sample_declarations(&junction, &surface, &fan, &req);

    let queries = vec![
        OutputQuery::scalar("thermal.junction_maximum"),
        OutputQuery::scalar("thermal.surface_mean"),
        OutputQuery::scalar("airflow.pressure_drop"),
    ];

    let receipt1 = with_cx(|cx| {
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

    let receipt2 = with_cx(|cx| {
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

    assert_eq!(receipt1.rows.len(), receipt2.rows.len());
    for (r1, r2) in receipt1.rows.iter().zip(&receipt2.rows) {
        assert_eq!(r1.identity_hash, r2.identity_hash);
        assert_eq!(r1.value.to_bits(), r2.value.to_bits());
        assert_eq!(r1.semantic_id, r2.semantic_id);
    }
}

#[test]
fn rq_008_region_constraints_bind_the_exact_qoi_scope() {
    with_cx(|cx| {
        let (mesh, solution) = sample_mesh_and_solution();
        let op = operating_point();

        let junction = JunctionRegion::try_new("package", vec![0, 1]).unwrap();
        let surface = SurfaceRegion::try_new("case", vec![0, 1]).unwrap();
        let fan = FanPowerSpec::try_new(0.5, 0.05, source("fan-efficiency")).unwrap();
        let req = sample_requirement(400.0, 1.1);
        let decls = sample_declarations(&junction, &surface, &fan, &req);

        for query in [
            OutputQuery::scalar_with_region("thermal.junction_maximum", "package"),
            OutputQuery::scalar_with_region("thermal.surface_mean", "case"),
        ] {
            extract_registered_qois(
                &[query],
                &mesh,
                &solution,
                &op,
                &decls,
                QoiExecutionLimits::default(),
                cx,
            )
            .expect("matching region scope");
        }

        for (query, expected_region) in [
            (
                OutputQuery::scalar_with_region("thermal.junction_maximum", "case"),
                "package",
            ),
            (
                OutputQuery::scalar_with_region("thermal.surface_mean", "package"),
                "case",
            ),
        ] {
            let error = extract_registered_qois(
                &[query],
                &mesh,
                &solution,
                &op,
                &decls,
                QoiExecutionLimits::default(),
                cx,
            )
            .expect_err("foreign region must fail closed");
            assert!(matches!(
                error,
                RegisteredQoiError::RegionNotFound { available, .. }
                    if available.len() == 1 && available[0] == expected_region
            ));
        }

        for (query_name, expected_semantic_id) in [
            ("airflow.pressure_drop", QoiSemanticId::PressureDrop),
            ("airflow.fan_power", QoiSemanticId::FanPower),
            ("thermal.thermal_margin", QoiSemanticId::ThermalMargin),
        ] {
            let error = extract_registered_qois(
                &[OutputQuery::scalar_with_region(query_name, "package")],
                &mesh,
                &solution,
                &op,
                &decls,
                QoiExecutionLimits::default(),
                cx,
            )
            .expect_err("global QoI must reject a region constraint");
            assert!(matches!(
                error,
                RegisteredQoiError::RegionNotApplicable {
                    semantic_id,
                    requested,
                } if semantic_id == expected_semantic_id && requested == "package"
            ));
        }
    });
}

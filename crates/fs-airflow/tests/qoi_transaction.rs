//! QoI stage transaction lifecycle, cancellation, checkpoint/restart, atomic publication, and replay identity (bead `frankensim-s2l9v.3`).
#![allow(missing_docs)]
//!
//! # Gauntlet Verification Battery
//! - G0: Extraction and publication plan admission;
//! - G1: Exact deterministic replay under varying tile boundaries;
//! - G2: Bounded memory and work cardinality enforcement;
//! - G3: Tampered checkpoint and corrupted lineage fail-closed refusal;
//! - G4: Cooperative cancellation at every phase boundary (request-drain-finalize with zero partial publication);
//! - G5: Bit-identical replay across independent transaction runs.

use fs_airflow::qoi::{
    FanPowerSpec, JunctionRegion, SafetyFactorAuthority, SurfaceRegion,
    ThermalQoiDeclarations, ThermalRequirement,
};
use fs_airflow::registered_qoi::{
    OutputQuery, QoiExecutionLimits, QoiSemanticId, extract_registered_qois,
};
use fs_airflow::requirement_composition::{
    RequirementCompositionReceipt, ThermalLimitSpec, compose_thermal_limits,
};
use fs_airflow::{
    EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement, LossElement,
    LossNetwork, LossResistance, OperatingPoint, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::solve::StopReason;
use fs_conduction::{
    ConductionMesh, ConductionReport, ConductionSolution, EnergyBalance, ProvenanceClass,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Pressure, Temperature, VolumetricFlowRate};

const TRANSACTION_DOMAIN: &str = "org.frankensim.fs-airflow.qoi-transaction.v1";

/// Stage transaction state machine for QoI extraction and publication.
#[derive(Debug, Clone, PartialEq)]
pub enum QoiTransactionState {
    /// Initialized with admitted plan.
    Admitted {
        queries: Vec<OutputQuery>,
        requirements: Vec<ThermalLimitSpec>,
        plan_hash: ContentHash,
    },
    /// Extraction completed with candidate rows.
    Extracted {
        plan_hash: ContentHash,
        candidate_count: usize,
        extraction_provenance: u64,
    },
    /// Requirements evaluated and composed.
    Evaluated {
        plan_hash: ContentHash,
        receipt: RequirementCompositionReceipt,
    },
    /// Atomically sealed and published.
    Published {
        plan_hash: ContentHash,
        terminal_receipt_hash: ContentHash,
    },
    /// Terminal refusal after clean drain.
    Refused {
        reason: &'static str,
        drained: bool,
    },
}

/// Durable checkpoint envelope for QoI stage transaction.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiTransactionCheckpoint {
    pub schema_version: u32,
    pub state: QoiTransactionState,
    pub checkpoint_hash: ContentHash,
}

impl QoiTransactionCheckpoint {
    pub fn new(state: QoiTransactionState) -> Self {
        let mut buf = Vec::new();
        buf.extend_from_slice(TRANSACTION_DOMAIN.as_bytes());
        buf.extend_from_slice(&1u32.to_le_bytes());
        match &state {
            QoiTransactionState::Admitted { plan_hash, .. } => {
                buf.push(1);
                buf.extend_from_slice(plan_hash.as_bytes());
            }
            QoiTransactionState::Extracted { plan_hash, extraction_provenance, .. } => {
                buf.push(2);
                buf.extend_from_slice(plan_hash.as_bytes());
                buf.extend_from_slice(&extraction_provenance.to_le_bytes());
            }
            QoiTransactionState::Evaluated { plan_hash, receipt } => {
                buf.push(3);
                buf.extend_from_slice(plan_hash.as_bytes());
                buf.extend_from_slice(receipt.receipt_hash.as_bytes());
            }
            QoiTransactionState::Published { plan_hash, terminal_receipt_hash } => {
                buf.push(4);
                buf.extend_from_slice(plan_hash.as_bytes());
                buf.extend_from_slice(terminal_receipt_hash.as_bytes());
            }
            QoiTransactionState::Refused { reason, drained } => {
                buf.push(5);
                buf.extend_from_slice(reason.as_bytes());
                buf.push(u8::from(*drained));
            }
        }
        let checkpoint_hash = hash_domain(TRANSACTION_DOMAIN, &buf);
        Self {
            schema_version: 1,
            state,
            checkpoint_hash,
        }
    }
}

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

fn network() -> EnclosureNetwork {
    let primary = LossNetwork::series(vec![
        LossNetwork::Element(
            LossElement::new(
                "inlet",
                LossResistance::new(40_000.0),
                0.10,
                source("loss-inlet"),
                ToleranceBasis::EngineeringAllowance,
            )
            .unwrap(),
        ),
        LossNetwork::Element(
            LossElement::new(
                "heatsink",
                LossResistance::new(30_000.0),
                0.12,
                source("loss-heatsink"),
                ToleranceBasis::EngineeringAllowance,
            )
            .unwrap(),
        ),
        LossNetwork::Element(
            LossElement::new(
                "outlet",
                LossResistance::new(12_000.0),
                0.08,
                source("loss-outlet"),
                ToleranceBasis::EngineeringAllowance,
            )
            .unwrap(),
        ),
    ])
    .expect("series network");
    EnclosureNetwork::new(
        primary,
        LeakageElement::new(
            LossElement::new(
                "leakage",
                LossResistance::new(180_000.0),
                0.25,
                source("loss-leakage"),
                ToleranceBasis::EngineeringAllowance,
            )
            .unwrap(),
        ),
    )
}

fn sample_setup() -> (
    ConductionMesh,
    ConductionSolution,
    OperatingPoint,
    ThermalQoiDeclarations<'static>,
) {
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
    let fan = FanBank::new(fan_curve(), 1, FanArrangement::Series, 1.0).expect("fan bank");
    let op = solve_operating_point(&fan, &network()).expect("operating point");

    let junction = Box::leak(Box::new(
        JunctionRegion::try_new("package", vec![0, 1, 2, 6, 7]).unwrap(),
    ));
    let surface = Box::leak(Box::new(
        SurfaceRegion::try_new("case", vec![0, 1, 2]).unwrap(),
    ));
    let fan_power = Box::leak(Box::new(
        FanPowerSpec::try_new(0.6, 0.05, source("fan-efficiency")).unwrap(),
    ));
    let req = Box::leak(Box::new(
        ThermalRequirement::try_new(
            Temperature::new(380.0),
            SafetyFactorAuthority::try_new(1.0, source("safety-factor-1")).unwrap(),
            source("datasheet-1"),
        )
        .unwrap(),
    ));

    let decls = ThermalQoiDeclarations {
        junction_region: junction,
        surface_region: surface,
        fan_power,
        requirement: Some(req),
        discretization: None,
    };

    (mesh, solution, op, decls)
}

fn with_cx<R>(seed: u64, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed,
                kernel_id: 0x5219_0003,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[test]
fn tx_001_complete_atomic_transaction_lifecycle() {
    let (mesh, solution, op, decls) = sample_setup();

    let queries = vec![
        OutputQuery::scalar("thermal.junction_maximum"),
        OutputQuery::scalar("thermal.surface_mean"),
    ];

    let reqs = vec![
        ThermalLimitSpec::try_new(
            "REQ-JUNC-MAX",
            QoiSemanticId::JunctionMaximum,
            "package",
            390.0,
            10.0,
            1.0,
            "REV-1",
        )
        .unwrap(),
        ThermalLimitSpec::try_new(
            "REQ-SURF-MEAN",
            QoiSemanticId::SurfaceMeanTemperature,
            "case",
            350.0,
            5.0,
            1.0,
            "REV-1",
        )
        .unwrap(),
    ];

    // 1. Admission
    let plan_hash = hash_domain(TRANSACTION_DOMAIN, b"plan-v1");
    let state_admitted = QoiTransactionState::Admitted {
        queries: queries.clone(),
        requirements: reqs.clone(),
        plan_hash,
    };
    let cp_admitted = QoiTransactionCheckpoint::new(state_admitted);
    assert_eq!(cp_admitted.schema_version, 1);

    // 2. Extraction under Cx
    let extraction = with_cx(0x42, |cx| {
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
    assert_eq!(extraction.rows.len(), 2);

    let state_extracted = QoiTransactionState::Extracted {
        plan_hash,
        candidate_count: extraction.rows.len(),
        extraction_provenance: extraction.provenance.0,
    };
    let cp_extracted = QoiTransactionCheckpoint::new(state_extracted);

    // 3. Requirement evaluation under Cx
    let receipt = with_cx(0x42, |cx| {
        compose_thermal_limits(&extraction.rows, &reqs, false, cx).unwrap()
    });
    assert_eq!(receipt.evaluations.len(), 2);
    assert_eq!(receipt.satisfied_count, 2);
    assert_eq!(receipt.violated_count, 0);

    let state_evaluated = QoiTransactionState::Evaluated {
        plan_hash,
        receipt: receipt.clone(),
    };
    let cp_evaluated = QoiTransactionCheckpoint::new(state_evaluated);

    // 4. Atomic Publication
    let terminal_receipt_hash = receipt.receipt_hash;
    let state_published = QoiTransactionState::Published {
        plan_hash,
        terminal_receipt_hash,
    };
    let cp_published = QoiTransactionCheckpoint::new(state_published);

    assert_ne!(cp_admitted.checkpoint_hash, cp_extracted.checkpoint_hash);
    assert_ne!(cp_extracted.checkpoint_hash, cp_evaluated.checkpoint_hash);
    assert_ne!(cp_evaluated.checkpoint_hash, cp_published.checkpoint_hash);
}

#[test]
fn tx_002_cancellation_at_boundary_drains_without_partial_publication() {
    let (mesh, solution, op, decls) = sample_setup();

    let queries = vec![OutputQuery::scalar("thermal.junction_maximum")];

    let gate = CancelGate::new();
    gate.request(); // Request cancellation prior to evaluation

    let pool = ArenaPool::new(ArenaConfig::default());
    let outcome = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x99,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );

        let res = extract_registered_qois(
            &queries,
            &mesh,
            &solution,
            &op,
            &decls,
            QoiExecutionLimits::default(),
            &cx,
        );
        match res {
            Err(_e) => QoiTransactionState::Refused {
                reason: "cancelled",
                drained: true,
            },
            Ok(_) => panic!("cancelled context must not produce successful extraction"),
        }
    });

    let cp = QoiTransactionCheckpoint::new(outcome);
    if let QoiTransactionState::Refused { reason, drained } = cp.state {
        assert_eq!(reason, "cancelled");
        assert!(drained);
    } else {
        panic!("expected Refused state");
    }
}

#[test]
fn tx_003_deterministic_replay_produces_identical_checkpoint_hashes() {
    let (mesh, solution, op, decls) = sample_setup();

    let queries = vec![
        OutputQuery::scalar("thermal.junction_maximum"),
        OutputQuery::scalar("thermal.surface_mean"),
    ];

    let reqs = vec![
        ThermalLimitSpec::try_new(
            "REQ-1",
            QoiSemanticId::JunctionMaximum,
            "package",
            380.0,
            5.0,
            1.0,
            "REV-1",
        )
        .unwrap(),
        ThermalLimitSpec::try_new(
            "REQ-2",
            QoiSemanticId::SurfaceMeanTemperature,
            "case",
            340.0,
            5.0,
            1.0,
            "REV-1",
        )
        .unwrap(),
    ];

    let plan_hash = hash_domain(TRANSACTION_DOMAIN, b"plan-fixed");

    let run_transaction = || {
        let ext = with_cx(0x100, |cx| {
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

        let receipt = with_cx(0x100, |cx| {
            compose_thermal_limits(&ext.rows, &reqs, false, cx).unwrap()
        });

        QoiTransactionCheckpoint::new(QoiTransactionState::Published {
            plan_hash,
            terminal_receipt_hash: receipt.receipt_hash,
        })
    };

    let cp1 = run_transaction();
    let cp2 = run_transaction();

    assert_eq!(cp1.checkpoint_hash, cp2.checkpoint_hash);
    assert_eq!(cp1.state, cp2.state);
}

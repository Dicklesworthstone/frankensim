//! G0/G3 battery for floating-to-target deployment refinement contracts.

use fs_contract::Interval;
use fs_contract::deployment::{
    ArtifactIdentity, AxisEvidence, BoundedSet, DeployedTargetIdentity, DeploymentRefinementError,
    DeploymentRefinementProblem, DeploymentRefinementSpec, EvidenceBasis, EvidenceVerdict,
    FaultContract, InterfaceRelation, NumericContract, OfflineProofBudget, PermittedError,
    ProofAxis, RefinementRelation, TimingContract, TransitionSystemIdentity,
    discharge_universal_claim,
};
use std::collections::{BTreeMap, BTreeSet};

fn strings(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

fn artifact(name: &str, version: &str, byte: u8) -> ArtifactIdentity {
    ArtifactIdentity {
        name: name.to_string(),
        version: version.to_string(),
        content_hash: [byte; 32],
    }
}

fn spec() -> DeploymentRefinementSpec {
    let source = TransitionSystemIdentity {
        artifact: artifact("floating-controller", "1.3.0", 1),
        state_schema: "controller-state/f64.v2".to_string(),
        observation_schema: "shaft-observation/f64.v1".to_string(),
    };
    let target_system = TransitionSystemIdentity {
        artifact: artifact("deployed-controller", "1.3.0-q31", 2),
        state_schema: "controller-state/q31.v2".to_string(),
        observation_schema: "shaft-observation/adc16.v1".to_string(),
    };
    DeploymentRefinementSpec {
        interface: InterfaceRelation {
            source_state_schema: source.state_schema.clone(),
            target_state_schema: target_system.state_schema.clone(),
            state_map_id: "state-map/f64-to-q31.v1".to_string(),
            source_observation_schema: source.observation_schema.clone(),
            target_observation_schema: target_system.observation_schema.clone(),
            observation_map_id: "observation-map/adc16.v1".to_string(),
            source_unit: "rad/s".to_string(),
            target_unit: "rad/s".to_string(),
            source_frame: "motor-shaft".to_string(),
            target_frame: "motor-shaft".to_string(),
        },
        source,
        target: DeployedTargetIdentity {
            system: target_system,
            target_triple: "thumbv7em-none-eabihf".to_string(),
            device_revision: "drive-board-r4/mcu-b2".to_string(),
            toolchain: artifact("rust-codegen-bundle", "nightly-2026-07-01", 3),
            capabilities: strings(&["fixed-point-q31", "monotonic-clock"]),
        },
        timing: TimingContract {
            source_period_ns: 1_000_000,
            target_period_ns: 100_000,
            max_latency_ns: 80_000,
            max_jitter_ns: 10_000,
        },
        numeric: NumericContract {
            quantization_step: 1.0 / 2_147_483_648.0,
            saturation: Interval::new(-1.0, 1.0).unwrap(),
        },
        plant: artifact("motor-plant-abstraction", "4.2.1", 4),
        environment: BoundedSet::new("drive-environment.v1")
            .with("temperature_k", Interval::new(250.0, 350.0).unwrap())
            .with("bus_voltage_v", Interval::new(42.0, 58.0).unwrap()),
        disturbances: BoundedSet::new("shaft-disturbances.v2")
            .with("load_torque_nm", Interval::new(-2.0, 2.0).unwrap()),
        faults: FaultContract {
            model_id: "drive-faults.v3".to_string(),
            required_faults: strings(&["sensor-dropout", "deadline-overrun"]),
            modeled_faults: strings(&["sensor-dropout", "deadline-overrun", "watchdog-reset"]),
            source_safe_state: "zero-torque-latched".to_string(),
            target_safe_state: "zero-torque-latched".to_string(),
        },
        horizon_steps: 10_000,
        invariant: "no-overcurrent-and-bounded-speed.v2".to_string(),
        control_objective: "speed-tracking-with-current-limit.v3".to_string(),
        relation: RefinementRelation::ApproximateSimulation,
        permitted_error: PermittedError {
            numeric_abs: 1.0e-5,
            temporal_ns: 100_000,
            functional_relation: "bounded-output-trace-distance.v2".to_string(),
            safety_relation: "robust-invariant-preservation.v1".to_string(),
        },
        proof_budget: OfflineProofBudget {
            max_work_units: 10_000_000,
            max_memory_bytes: 512 * 1024 * 1024,
            max_wall_time_ns: 3_600_000_000_000,
            cancellation_poll_stride: 1_024,
            required_capabilities: strings(&["fixed-point-q31", "monotonic-clock"]),
        },
        assumptions: BTreeMap::from([
            (
                "scheduler".to_string(),
                "single-core-static-priority.v2".to_string(),
            ),
            ("rounding".to_string(), "nearest-ties-even".to_string()),
        ]),
    }
}

fn problem() -> DeploymentRefinementProblem {
    DeploymentRefinementProblem::admit(spec()).expect("valid deployment refinement problem")
}

fn evidence(root: u64) -> Vec<AxisEvidence> {
    [
        ProofAxis::Numeric,
        ProofAxis::Temporal,
        ProofAxis::Functional,
        ProofAxis::Safety,
    ]
    .into_iter()
    .enumerate()
    .map(|(index, axis)| AxisEvidence {
        axis,
        manifest_root: root,
        checker: artifact("independent-refinement-checker", "2.1.0", 10 + index as u8),
        evidence_hash: [20 + index as u8; 32],
        basis: if axis == ProofAxis::Functional {
            EvidenceBasis::ExhaustiveModelCheck
        } else {
            EvidenceBasis::StaticProof
        },
        verdict: EvidenceVerdict::Established,
    })
    .collect()
}

#[test]
fn admitted_problem_exposes_four_distinct_proof_axes() {
    let problem = problem();
    let obligations = problem.proof_obligations();
    assert_eq!(
        obligations.each_ref().map(|obligation| obligation.axis),
        [
            ProofAxis::Numeric,
            ProofAxis::Temporal,
            ProofAxis::Functional,
            ProofAxis::Safety,
        ]
    );
    assert!(
        obligations
            .iter()
            .all(|obligation| obligation.relation == RefinementRelation::ApproximateSimulation)
    );
    assert_eq!(
        obligations
            .iter()
            .map(|obligation| obligation.requirement)
            .collect::<BTreeSet<_>>()
            .len(),
        4,
        "no proof axis is an alias for another"
    );
}

#[test]
fn missing_maps_and_incompatible_schemas_clocks_units_and_frames_refuse() {
    let mut candidate = spec();
    candidate.interface.state_map_id.clear();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::MissingStateMap
    );

    let mut candidate = spec();
    candidate.interface.observation_map_id.clear();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::MissingObservationMap
    );

    let mut candidate = spec();
    candidate.interface.target_observation_schema = "wrong-observation.v9".to_string();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::RelationSchemaMismatch {
            relation: "observation"
        }
    );

    let mut candidate = spec();
    candidate.timing.target_period_ns = 333_333;
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::IncompatibleClocks
    );

    let mut candidate = spec();
    candidate.interface.target_unit = "rpm".to_string();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::UnitMismatch
    );

    let mut candidate = spec();
    candidate.interface.target_frame = "gear-output".to_string();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::FrameMismatch
    );
}

#[test]
fn fault_omission_safe_state_mismatch_and_missing_capability_refuse() {
    let mut candidate = spec();
    candidate.faults.modeled_faults.remove("sensor-dropout");
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::FaultOmission {
            fault: "sensor-dropout".to_string()
        }
    );

    let mut candidate = spec();
    candidate.faults.target_safe_state = "coast-unlatched".to_string();
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::SafeStateMismatch
    );

    let mut candidate = spec();
    candidate.target.capabilities.remove("fixed-point-q31");
    assert_eq!(
        DeploymentRefinementProblem::admit(candidate).unwrap_err(),
        DeploymentRefinementError::MissingCapability {
            capability: "fixed-point-q31".to_string()
        }
    );
}

#[test]
fn saturation_target_version_and_environment_enlargement_invalidate_frozen_claim() {
    let frozen_problem = problem();
    let manifest = frozen_problem.freeze();

    let mut changed = spec();
    changed.numeric.saturation = Interval::new(-1.1, 1.0).unwrap();
    let changed = DeploymentRefinementProblem::admit(changed).unwrap();
    assert_eq!(
        manifest.admit_live(&changed).unwrap_err(),
        DeploymentRefinementError::SaturationDrift
    );

    let mut changed = spec();
    changed.target.toolchain.version = "nightly-2026-07-02".to_string();
    let changed = DeploymentRefinementProblem::admit(changed).unwrap();
    assert_eq!(
        manifest.admit_live(&changed).unwrap_err(),
        DeploymentRefinementError::TargetIdentityDrift
    );

    let mut changed = spec();
    changed.environment = BoundedSet::new("drive-environment.v1")
        .with("temperature_k", Interval::new(240.0, 350.0).unwrap())
        .with("bus_voltage_v", Interval::new(42.0, 58.0).unwrap());
    let changed = DeploymentRefinementProblem::admit(changed).unwrap();
    assert_eq!(
        manifest.admit_live(&changed).unwrap_err(),
        DeploymentRefinementError::EnvironmentEnlargement {
            quantity: "temperature_k".to_string()
        }
    );

    let mut narrowed = spec();
    narrowed.environment = BoundedSet::new("drive-environment.v1")
        .with("temperature_k", Interval::new(260.0, 340.0).unwrap())
        .with("bus_voltage_v", Interval::new(44.0, 56.0).unwrap())
        .with("humidity_fraction", Interval::new(0.0, 0.4).unwrap());
    let narrowed = DeploymentRefinementProblem::admit(narrowed).unwrap();
    manifest
        .admit_live(&narrowed)
        .expect("additional/narrower environment constraints remain inside the proof domain");
}

#[test]
fn changed_fault_and_safe_state_contracts_invalidate_frozen_claim() {
    let manifest = problem().freeze();

    let mut changed = spec();
    changed.faults.model_id = "drive-faults.v4".to_string();
    let changed = DeploymentRefinementProblem::admit(changed).unwrap();
    assert_eq!(
        manifest.admit_live(&changed).unwrap_err(),
        DeploymentRefinementError::FaultModelDrift
    );

    let mut changed = spec();
    changed.faults.source_safe_state = "zero-current-latched".to_string();
    changed.faults.target_safe_state = "zero-current-latched".to_string();
    let changed = DeploymentRefinementProblem::admit(changed).unwrap();
    assert_eq!(
        manifest.admit_live(&changed).unwrap_err(),
        DeploymentRefinementError::SafeStateDrift
    );
}

#[test]
fn retained_manifest_replay_is_versioned_self_consistent_and_exact() {
    let first = problem().freeze();
    let second = problem().freeze();
    assert_eq!(first.root(), second.root());
    assert_eq!(first.canonical_bytes(), second.canonical_bytes());
    first
        .admit_retained(
            first.schema_version(),
            first.canonical_bytes(),
            first.root(),
        )
        .unwrap();

    assert!(matches!(
        first.admit_retained(
            first.schema_version() + 1,
            first.canonical_bytes(),
            first.root()
        ),
        Err(DeploymentRefinementError::SchemaVersionMismatch { .. })
    ));

    let mut tampered = first.canonical_bytes().to_vec();
    let last = tampered.last_mut().unwrap();
    *last ^= 1;
    assert!(matches!(
        first.admit_retained(first.schema_version(), &tampered, first.root()),
        Err(DeploymentRefinementError::RetainedRootMismatch { .. })
    ));

    let mut changed = spec();
    changed
        .assumptions
        .insert("rounding".to_string(), "toward-zero".to_string());
    let other = DeploymentRefinementProblem::admit(changed)
        .unwrap()
        .freeze();
    assert_eq!(
        first
            .admit_retained(
                other.schema_version(),
                other.canonical_bytes(),
                other.root()
            )
            .unwrap_err(),
        DeploymentRefinementError::RetainedManifestMismatch
    );
}

#[test]
fn universal_receipt_requires_established_proof_on_every_axis() {
    let live = problem();
    let manifest = live.freeze();
    let receipt = discharge_universal_claim(&manifest, &live, evidence(manifest.root())).unwrap();
    assert_eq!(receipt.manifest_root(), manifest.root());
    assert_eq!(
        receipt.relation(),
        RefinementRelation::ApproximateSimulation
    );
    assert_eq!(
        receipt
            .evidence()
            .iter()
            .map(|item| item.axis)
            .collect::<Vec<_>>(),
        vec![
            ProofAxis::Numeric,
            ProofAxis::Temporal,
            ProofAxis::Functional,
            ProofAxis::Safety,
        ]
    );

    let mut missing = evidence(manifest.root());
    missing.retain(|item| item.axis != ProofAxis::Safety);
    assert_eq!(
        discharge_universal_claim(&manifest, &live, missing).unwrap_err(),
        DeploymentRefinementError::MissingProofAxis {
            axis: ProofAxis::Safety
        }
    );

    let mut duplicate = evidence(manifest.root());
    duplicate.push(duplicate[0].clone());
    assert_eq!(
        discharge_universal_claim(&manifest, &live, duplicate).unwrap_err(),
        DeploymentRefinementError::DuplicateProofAxis {
            axis: ProofAxis::Numeric
        }
    );
}

#[test]
fn measurement_unknown_and_refutation_never_masquerade_as_universal_proof() {
    let live = problem();
    let manifest = live.freeze();

    let mut measured = evidence(manifest.root());
    measured[1].basis = EvidenceBasis::Measurement { runs: 50_000 };
    assert_eq!(
        discharge_universal_claim(&manifest, &live, measured).unwrap_err(),
        DeploymentRefinementError::MeasuredEvidenceIsNotUniversal {
            axis: ProofAxis::Temporal
        }
    );

    let mut unknown = evidence(manifest.root());
    unknown[2].verdict = EvidenceVerdict::Unknown;
    assert_eq!(
        discharge_universal_claim(&manifest, &live, unknown).unwrap_err(),
        DeploymentRefinementError::UnknownProofAxis {
            axis: ProofAxis::Functional
        }
    );

    let mut refuted = evidence(manifest.root());
    refuted[3].verdict = EvidenceVerdict::Refuted;
    assert_eq!(
        discharge_universal_claim(&manifest, &live, refuted).unwrap_err(),
        DeploymentRefinementError::RefutedProofAxis {
            axis: ProofAxis::Safety
        }
    );
}

#[test]
fn relation_strength_and_manifest_identity_cannot_be_replayed() {
    let simulation = problem();
    let simulation_manifest = simulation.freeze();

    let mut bisimulation_spec = spec();
    bisimulation_spec.relation = RefinementRelation::ApproximateBisimulation;
    let bisimulation = DeploymentRefinementProblem::admit(bisimulation_spec).unwrap();
    let bisimulation_manifest = bisimulation.freeze();
    assert_ne!(simulation_manifest.root(), bisimulation_manifest.root());
    assert_eq!(
        simulation_manifest.admit_live(&bisimulation).unwrap_err(),
        DeploymentRefinementError::ClaimAssumptionDrift
    );

    let stale = evidence(bisimulation_manifest.root());
    assert_eq!(
        discharge_universal_claim(&simulation_manifest, &simulation, stale).unwrap_err(),
        DeploymentRefinementError::EvidenceManifestMismatch {
            axis: ProofAxis::Numeric
        }
    );
}

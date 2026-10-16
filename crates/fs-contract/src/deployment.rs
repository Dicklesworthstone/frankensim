//! Floating-design to deployed-target refinement contracts.
//!
//! A deployment claim is admitted only after the source model, deployed
//! target, toolchain, state/observation relations, timing, numeric behavior,
//! environment, faults, and safe state have been frozen into one deterministic
//! manifest.  This module defines proof obligations; it does not prove them.

use crate::Interval;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Current canonical schema for deployment-refinement manifests.
pub const DEPLOYMENT_REFINEMENT_SCHEMA_VERSION: u32 = 1;

/// Maximum UTF-8 bytes across one admitted problem.
pub const MAX_REFINEMENT_TEXT_BYTES: usize = 65_536;

/// Maximum entries in any named set, fault set, capability set, or assumption map.
pub const MAX_REFINEMENT_SET_ENTRIES: usize = 64;

const MANIFEST_DOMAIN: &[u8] = b"fs-contract/deployment-refinement-manifest.v1";

/// A content-addressed, versioned artifact identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactIdentity {
    /// Stable artifact family or logical name.
    pub name: String,
    /// Exact semantic or producer version.
    pub version: String,
    /// Content digest supplied by the artifact producer.
    pub content_hash: [u8; 32],
}

/// A transition-system identity and the schemas its relations must consume.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TransitionSystemIdentity {
    /// Versioned transition-system artifact.
    pub artifact: ArtifactIdentity,
    /// Exact state-schema identity.
    pub state_schema: String,
    /// Exact observable-schema identity.
    pub observation_schema: String,
}

/// Exact deployed-target and toolchain identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeployedTargetIdentity {
    /// Transition system generated for the target.
    pub system: TransitionSystemIdentity,
    /// Rust-style target triple or equally exact target profile.
    pub target_triple: String,
    /// Exact device/board/silicon revision.
    pub device_revision: String,
    /// Exact compiler, linker, and code-generator artifact.
    pub toolchain: ArtifactIdentity,
    /// Capabilities actually present on this target.
    pub capabilities: BTreeSet<String>,
}

/// State, observation, unit, and frame relation between source and target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceRelation {
    /// Source state schema consumed by `state_map_id`.
    pub source_state_schema: String,
    /// Target state schema produced by `state_map_id`.
    pub target_state_schema: String,
    /// Versioned state-map identity. Empty means no relation and is refused.
    pub state_map_id: String,
    /// Source observation schema consumed by `observation_map_id`.
    pub source_observation_schema: String,
    /// Target observation schema produced by `observation_map_id`.
    pub target_observation_schema: String,
    /// Versioned observation-map identity. Empty means no relation and is refused.
    pub observation_map_id: String,
    /// Source-side unit identity.
    pub source_unit: String,
    /// Target-side unit identity.
    pub target_unit: String,
    /// Source coordinate-frame identity.
    pub source_frame: String,
    /// Target coordinate-frame identity.
    pub target_frame: String,
}

/// Sampling, latency, and jitter assumptions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TimingContract {
    /// Floating/source model sample period.
    pub source_period_ns: u64,
    /// Deployed target sample period.
    pub target_period_ns: u64,
    /// Maximum admitted end-to-end latency.
    pub max_latency_ns: u64,
    /// Maximum admitted release/arrival jitter.
    pub max_jitter_ns: u64,
}

/// Quantization and saturation semantics of the deployed controller.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NumericContract {
    /// Quantizer step in the declared control unit.
    pub quantization_step: f64,
    /// Inclusive deployed saturation interval.
    pub saturation: Interval,
}

/// A versioned bounded environment or disturbance set.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedSet {
    identity: String,
    bounds: BTreeMap<String, Interval>,
}

impl BoundedSet {
    /// Start a named bounded set. Empty bounds are allowed and mean the named
    /// model declares no varying coordinates, not an unnamed/unbounded set.
    #[must_use]
    pub fn new(identity: impl Into<String>) -> Self {
        Self {
            identity: identity.into(),
            bounds: BTreeMap::new(),
        }
    }

    /// Add or replace one named interval.
    #[must_use]
    pub fn with(mut self, quantity: impl Into<String>, interval: Interval) -> Self {
        self.bounds.insert(quantity.into(), interval);
        self
    }

    /// Stable set identity.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }

    /// Sorted bounds.
    pub fn bounds(&self) -> impl Iterator<Item = (&str, Interval)> {
        self.bounds
            .iter()
            .map(|(name, interval)| (name.as_str(), *interval))
    }

    fn offending_enlargement_from(&self, frozen: &Self) -> Option<String> {
        if self.identity != frozen.identity {
            return Some("<set-identity>".to_string());
        }
        for (quantity, frozen_bound) in &frozen.bounds {
            let Some(live_bound) = self.bounds.get(quantity) else {
                // Absence removes a constraint and therefore enlarges the set.
                return Some(quantity.clone());
            };
            if !frozen_bound.contains(live_bound) {
                return Some(quantity.clone());
            }
        }
        None
    }
}

/// Fault coverage and request-to-safe-state semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaultContract {
    /// Versioned fault-model identity.
    pub model_id: String,
    /// Faults the refinement claim promises to cover.
    pub required_faults: BTreeSet<String>,
    /// Faults represented by the deployed transition system.
    pub modeled_faults: BTreeSet<String>,
    /// Source/design safe-state identity.
    pub source_safe_state: String,
    /// Deployed safe-state identity.
    pub target_safe_state: String,
}

/// The exact strength of the requested refinement claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RefinementRelation {
    /// Every target trace is included in the permitted source trace set.
    TraceInclusion,
    /// One-way approximate simulation under the permitted-error relation.
    ApproximateSimulation,
    /// Two-way approximate simulation; strictly distinct from simulation.
    ApproximateBisimulation,
    /// Preservation of the named robust invariant.
    RobustInvariant,
    /// Preservation of the declared closed-loop performance bound.
    PerformanceBound,
}

/// The four independent proof axes. Evidence on one axis never discharges another.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ProofAxis {
    /// Quantization, saturation, and permitted numerical error.
    Numeric,
    /// Sampling, latency, jitter, and temporal error.
    Temporal,
    /// State/observation relations and control-objective behavior.
    Functional,
    /// Environment, faults, invariant, and safe-state behavior.
    Safety,
}

/// Permitted error is explicitly split across all four proof axes.
#[derive(Debug, Clone, PartialEq)]
pub struct PermittedError {
    /// Maximum absolute numeric deviation in the declared unit.
    pub numeric_abs: f64,
    /// Maximum temporal deviation.
    pub temporal_ns: u64,
    /// Named functional equivalence/error relation.
    pub functional_relation: String,
    /// Named safety refinement relation.
    pub safety_relation: String,
}

/// Bounded offline proof resources and cancellation polling cadence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OfflineProofBudget {
    /// Maximum checker/model-check work units.
    pub max_work_units: u64,
    /// Maximum resident bytes.
    pub max_memory_bytes: u64,
    /// Maximum wall-clock duration.
    pub max_wall_time_ns: u64,
    /// Maximum work units between cancellation polls.
    pub cancellation_poll_stride: u64,
    /// Capabilities required by the claim and checker.
    pub required_capabilities: BTreeSet<String>,
}

/// Public construction surface. Admission validates and seals it into a
/// [`DeploymentRefinementProblem`].
#[derive(Debug, Clone, PartialEq)]
pub struct DeploymentRefinementSpec {
    /// Floating/source transition system.
    pub source: TransitionSystemIdentity,
    /// Exact deployed target and toolchain.
    pub target: DeployedTargetIdentity,
    /// State, observation, unit, and frame relation.
    pub interface: InterfaceRelation,
    /// Clock, latency, and jitter assumptions.
    pub timing: TimingContract,
    /// Quantization and saturation assumptions.
    pub numeric: NumericContract,
    /// Versioned plant abstraction.
    pub plant: ArtifactIdentity,
    /// Admitted operating environment.
    pub environment: BoundedSet,
    /// Admitted disturbance set.
    pub disturbances: BoundedSet,
    /// Fault coverage and safe-state contract.
    pub faults: FaultContract,
    /// Finite trace horizon in source samples.
    pub horizon_steps: u64,
    /// Robust invariant identity.
    pub invariant: String,
    /// Closed-loop control-objective identity.
    pub control_objective: String,
    /// Exact relation strength requested.
    pub relation: RefinementRelation,
    /// Four-axis permitted error.
    pub permitted_error: PermittedError,
    /// Explicit proof/cancellation/capability budget.
    pub proof_budget: OfflineProofBudget,
    /// Frozen named assumptions not represented by another typed field.
    pub assumptions: BTreeMap<String, String>,
}

/// An admitted, immutable floating-to-target refinement problem.
#[derive(Debug, Clone, PartialEq)]
pub struct DeploymentRefinementProblem {
    spec: DeploymentRefinementSpec,
}

impl DeploymentRefinementProblem {
    /// Validate every typed seam and seal the problem.
    ///
    /// # Errors
    /// [`DeploymentRefinementError`] identifies the first deterministic
    /// refusal. No partially admitted problem is returned.
    pub fn admit(spec: DeploymentRefinementSpec) -> Result<Self, DeploymentRefinementError> {
        validate_spec(&spec)?;
        Ok(Self { spec })
    }

    /// The admitted immutable specification.
    #[must_use]
    pub fn spec(&self) -> &DeploymentRefinementSpec {
        &self.spec
    }

    /// Freeze this problem into a replayable, content-addressed manifest.
    #[must_use]
    pub fn freeze(&self) -> DeploymentRefinementManifest {
        let canonical_bytes = canonical_problem_bytes(&self.spec);
        let root = fnv1a64(&canonical_bytes);
        DeploymentRefinementManifest {
            frozen: self.clone(),
            canonical_bytes,
            root,
        }
    }

    /// The four independent obligations for the exact requested relation.
    #[must_use]
    pub fn proof_obligations(&self) -> [ProofObligation; 4] {
        [
            ProofObligation {
                axis: ProofAxis::Numeric,
                relation: self.spec.relation,
                requirement: "prove quantization and saturation stay within numeric_abs",
            },
            ProofObligation {
                axis: ProofAxis::Temporal,
                relation: self.spec.relation,
                requirement: "prove sampling, latency, and jitter stay within temporal_ns",
            },
            ProofObligation {
                axis: ProofAxis::Functional,
                relation: self.spec.relation,
                requirement: "prove the state/observation relation preserves the control objective",
            },
            ProofObligation {
                axis: ProofAxis::Safety,
                relation: self.spec.relation,
                requirement: "prove environment, faults, invariant, and safe state are preserved",
            },
        ]
    }
}

/// One axis-specific proof obligation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProofObligation {
    /// Axis this obligation alone may discharge.
    pub axis: ProofAxis,
    /// Exact refinement strength requested.
    pub relation: RefinementRelation,
    /// Stable human/agent-readable obligation summary.
    pub requirement: &'static str,
}

/// Frozen identity and exact bytes for one admitted problem.
#[derive(Debug, Clone, PartialEq)]
pub struct DeploymentRefinementManifest {
    frozen: DeploymentRefinementProblem,
    canonical_bytes: Vec<u8>,
    root: u64,
}

impl DeploymentRefinementManifest {
    /// Manifest schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        DEPLOYMENT_REFINEMENT_SCHEMA_VERSION
    }

    /// Deterministic FNV-1a root of [`Self::canonical_bytes`].
    #[must_use]
    pub const fn root(&self) -> u64 {
        self.root
    }

    /// Exact canonical bytes retained for replay and independent checking.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Frozen admitted problem.
    #[must_use]
    pub fn problem(&self) -> &DeploymentRefinementProblem {
        &self.frozen
    }

    /// Verify a retained transport against the frozen manifest.
    pub fn admit_retained(
        &self,
        schema_version: u32,
        canonical_bytes: &[u8],
        declared_root: u64,
    ) -> Result<(), DeploymentRefinementError> {
        if schema_version != DEPLOYMENT_REFINEMENT_SCHEMA_VERSION {
            return Err(DeploymentRefinementError::SchemaVersionMismatch {
                expected: DEPLOYMENT_REFINEMENT_SCHEMA_VERSION,
                found: schema_version,
            });
        }
        let computed = fnv1a64(canonical_bytes);
        if computed != declared_root {
            return Err(DeploymentRefinementError::RetainedRootMismatch {
                declared: declared_root,
                computed,
            });
        }
        if canonical_bytes != self.canonical_bytes {
            return Err(DeploymentRefinementError::RetainedManifestMismatch);
        }
        Ok(())
    }

    /// Admit a live deployment against the frozen proof assumptions.
    /// Environment/disturbance narrowing is sound; enlargement is refused.
    pub fn admit_live(
        &self,
        live: &DeploymentRefinementProblem,
    ) -> Result<(), DeploymentRefinementError> {
        let frozen = &self.frozen.spec;
        let live = &live.spec;

        if live.source != frozen.source {
            return Err(DeploymentRefinementError::SourceIdentityDrift);
        }
        if live.target != frozen.target {
            return Err(DeploymentRefinementError::TargetIdentityDrift);
        }
        if live.interface.state_map_id != frozen.interface.state_map_id
            || live.interface.source_state_schema != frozen.interface.source_state_schema
            || live.interface.target_state_schema != frozen.interface.target_state_schema
        {
            return Err(DeploymentRefinementError::StateRelationDrift);
        }
        if live.interface.observation_map_id != frozen.interface.observation_map_id
            || live.interface.source_observation_schema
                != frozen.interface.source_observation_schema
            || live.interface.target_observation_schema
                != frozen.interface.target_observation_schema
        {
            return Err(DeploymentRefinementError::ObservationRelationDrift);
        }
        if live.interface.source_unit != frozen.interface.source_unit
            || live.interface.target_unit != frozen.interface.target_unit
        {
            return Err(DeploymentRefinementError::UnitRelationDrift);
        }
        if live.interface.source_frame != frozen.interface.source_frame
            || live.interface.target_frame != frozen.interface.target_frame
        {
            return Err(DeploymentRefinementError::FrameRelationDrift);
        }
        if live.timing != frozen.timing {
            return Err(DeploymentRefinementError::TimingContractDrift);
        }
        if !same_f64(live.numeric.saturation.lo, frozen.numeric.saturation.lo)
            || !same_f64(live.numeric.saturation.hi, frozen.numeric.saturation.hi)
        {
            return Err(DeploymentRefinementError::SaturationDrift);
        }
        if !same_f64(
            live.numeric.quantization_step,
            frozen.numeric.quantization_step,
        ) {
            return Err(DeploymentRefinementError::QuantizationDrift);
        }
        if live.plant != frozen.plant {
            return Err(DeploymentRefinementError::PlantIdentityDrift);
        }
        if let Some(quantity) = live
            .environment
            .offending_enlargement_from(&frozen.environment)
        {
            return Err(DeploymentRefinementError::EnvironmentEnlargement { quantity });
        }
        if let Some(quantity) = live
            .disturbances
            .offending_enlargement_from(&frozen.disturbances)
        {
            return Err(DeploymentRefinementError::DisturbanceEnlargement { quantity });
        }
        if live.faults.source_safe_state != frozen.faults.source_safe_state
            || live.faults.target_safe_state != frozen.faults.target_safe_state
        {
            return Err(DeploymentRefinementError::SafeStateDrift);
        }
        if live.faults != frozen.faults {
            return Err(DeploymentRefinementError::FaultModelDrift);
        }
        if live.horizon_steps != frozen.horizon_steps
            || live.invariant != frozen.invariant
            || live.control_objective != frozen.control_objective
            || live.relation != frozen.relation
            || live.permitted_error != frozen.permitted_error
            || live.proof_budget != frozen.proof_budget
            || live.assumptions != frozen.assumptions
        {
            return Err(DeploymentRefinementError::ClaimAssumptionDrift);
        }
        Ok(())
    }
}

/// Authority type behind one axis of evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceBasis {
    /// Machine-checked symbolic/static proof.
    StaticProof,
    /// Exhaustive model checking over the frozen finite abstraction.
    ExhaustiveModelCheck,
    /// Empirical runs. Useful evidence, but never a universal proof.
    Measurement {
        /// Number of retained measured runs.
        runs: u64,
    },
}

/// Scientific terminal state of one proof attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceVerdict {
    /// The named checker established its axis obligation.
    Established,
    /// The checker could not decide within the frozen budget.
    Unknown,
    /// A counterexample refuted the axis obligation.
    Refuted,
}

/// Evidence offered for exactly one proof axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisEvidence {
    /// Axis this evidence may discharge.
    pub axis: ProofAxis,
    /// Frozen manifest root checked by the evidence producer.
    pub manifest_root: u64,
    /// Independent checker/model-checker identity.
    pub checker: ArtifactIdentity,
    /// Content hash of the proof or retained evidence bundle.
    pub evidence_hash: [u8; 32],
    /// Why the evidence has proof authority (or only measurement scope).
    pub basis: EvidenceBasis,
    /// Established, Unknown, or Refuted remain distinct.
    pub verdict: EvidenceVerdict,
}

/// Universal receipt. Construction is private and only the four-axis
/// discharge gate can create one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeploymentRefinementReceipt {
    manifest_root: u64,
    relation: RefinementRelation,
    evidence: Vec<AxisEvidence>,
}

impl DeploymentRefinementReceipt {
    /// Exact frozen problem this receipt proves.
    #[must_use]
    pub const fn manifest_root(&self) -> u64 {
        self.manifest_root
    }

    /// Exact relation strength proved; stronger relations are never inferred.
    #[must_use]
    pub const fn relation(&self) -> RefinementRelation {
        self.relation
    }

    /// Evidence in stable numeric/temporal/functional/safety order.
    #[must_use]
    pub fn evidence(&self) -> &[AxisEvidence] {
        &self.evidence
    }
}

/// Discharge all four proof axes for one exact live deployment.
///
/// Measurement, Unknown, Refuted, duplicate/missing axes, stale artifacts, and
/// assumption drift all fail closed. A successful receipt means only the exact
/// relation in the frozen manifest; it does not silently strengthen simulation
/// into bisimulation or a performance bound into safety.
pub fn discharge_universal_claim(
    manifest: &DeploymentRefinementManifest,
    live: &DeploymentRefinementProblem,
    evidence: impl IntoIterator<Item = AxisEvidence>,
) -> Result<DeploymentRefinementReceipt, DeploymentRefinementError> {
    manifest.admit_live(live)?;

    let mut by_axis = BTreeMap::new();
    for item in evidence {
        let axis = item.axis;
        if by_axis.contains_key(&axis) {
            return Err(DeploymentRefinementError::DuplicateProofAxis { axis });
        }
        if item.manifest_root != manifest.root {
            return Err(DeploymentRefinementError::EvidenceManifestMismatch { axis });
        }
        let mut checker_text_bytes = 0usize;
        validate_artifact(&item.checker, "checker", &mut checker_text_bytes)?;
        if item.evidence_hash.iter().all(|byte| *byte == 0) {
            return Err(DeploymentRefinementError::InvalidEvidenceHash { axis });
        }
        match item.verdict {
            EvidenceVerdict::Unknown => {
                return Err(DeploymentRefinementError::UnknownProofAxis { axis });
            }
            EvidenceVerdict::Refuted => {
                return Err(DeploymentRefinementError::RefutedProofAxis { axis });
            }
            EvidenceVerdict::Established => {}
        }
        if matches!(item.basis, EvidenceBasis::Measurement { .. }) {
            return Err(DeploymentRefinementError::MeasuredEvidenceIsNotUniversal { axis });
        }
        by_axis.insert(axis, item);
    }

    let mut ordered = Vec::with_capacity(4);
    for axis in [
        ProofAxis::Numeric,
        ProofAxis::Temporal,
        ProofAxis::Functional,
        ProofAxis::Safety,
    ] {
        let Some(item) = by_axis.remove(&axis) else {
            return Err(DeploymentRefinementError::MissingProofAxis { axis });
        };
        ordered.push(item);
    }

    Ok(DeploymentRefinementReceipt {
        manifest_root: manifest.root,
        relation: manifest.frozen.spec.relation,
        evidence: ordered,
    })
}

/// Typed, deterministic refusal surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DeploymentRefinementError {
    /// A required string was empty.
    EmptyField { field: &'static str },
    /// A string exceeded the aggregate bounded input contract.
    TextBudgetExceeded { bytes: usize, limit: usize },
    /// A set exceeded its bounded entry count.
    TooManyEntries { field: &'static str, count: usize },
    /// A content identity used the all-zero sentinel.
    MissingContentHash { field: &'static str },
    /// State relation was omitted.
    MissingStateMap,
    /// Observation relation was omitted.
    MissingObservationMap,
    /// A relation names schemas other than its source/target systems.
    RelationSchemaMismatch { relation: &'static str },
    /// v0 has no implicit unit conversion.
    UnitMismatch,
    /// v0 has no implicit frame transform.
    FrameMismatch,
    /// Sampling clocks are zero or not integer-commensurate.
    IncompatibleClocks,
    /// Timing envelope is internally contradictory.
    InvalidTimingEnvelope,
    /// Quantization or saturation is non-finite, unordered, or inconsistent.
    InvalidNumericContract,
    /// A required fault is absent from the modeled fault set.
    FaultOmission { fault: String },
    /// Source and target safe-state identities differ.
    SafeStateMismatch,
    /// A horizon, resource budget, or cancellation stride is zero/invalid.
    InvalidBudget { field: &'static str },
    /// Target lacks a capability required by the proof/claim.
    MissingCapability { capability: String },
    /// Permitted numeric error is non-finite or negative.
    InvalidPermittedError,
    /// An interval in a named set bypassed `Interval::new` invariants.
    InvalidBound { set: &'static str, quantity: String },
    /// Retained transport uses the wrong producer schema.
    SchemaVersionMismatch { expected: u32, found: u32 },
    /// Retained root does not derive from the retained bytes.
    RetainedRootMismatch { declared: u64, computed: u64 },
    /// Retained bytes are self-consistent but belong to another problem.
    RetainedManifestMismatch,
    /// Source transition system changed after freezing.
    SourceIdentityDrift,
    /// Target, target version, device, capability, or toolchain changed.
    TargetIdentityDrift,
    /// State relation changed.
    StateRelationDrift,
    /// Observation relation changed.
    ObservationRelationDrift,
    /// Unit relation changed.
    UnitRelationDrift,
    /// Frame relation changed.
    FrameRelationDrift,
    /// Timing assumptions changed.
    TimingContractDrift,
    /// Saturation changed.
    SaturationDrift,
    /// Quantization changed.
    QuantizationDrift,
    /// Plant abstraction changed.
    PlantIdentityDrift,
    /// Live environment is not a subset of the frozen environment.
    EnvironmentEnlargement { quantity: String },
    /// Live disturbance set is not a subset of the frozen disturbance set.
    DisturbanceEnlargement { quantity: String },
    /// Safe state changed after freezing.
    SafeStateDrift,
    /// Fault model changed after freezing.
    FaultModelDrift,
    /// Horizon, invariant, objective, relation, error, budget, or assumption drifted.
    ClaimAssumptionDrift,
    /// More than one evidence item attempted to discharge an axis.
    DuplicateProofAxis { axis: ProofAxis },
    /// No evidence item discharged an axis.
    MissingProofAxis { axis: ProofAxis },
    /// Evidence belongs to another frozen problem.
    EvidenceManifestMismatch { axis: ProofAxis },
    /// Evidence used an absent all-zero content hash.
    InvalidEvidenceHash { axis: ProofAxis },
    /// Checker returned Unknown.
    UnknownProofAxis { axis: ProofAxis },
    /// Checker returned a counterexample.
    RefutedProofAxis { axis: ProofAxis },
    /// Measurements cannot discharge a universal refinement obligation.
    MeasuredEvidenceIsNotUniversal { axis: ProofAxis },
}

impl fmt::Display for DeploymentRefinementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyField { field } => write!(f, "required refinement field `{field}` is empty"),
            Self::TextBudgetExceeded { bytes, limit } => write!(
                f,
                "refinement text requires {bytes} bytes; bounded limit is {limit}"
            ),
            Self::TooManyEntries { field, count } => write!(
                f,
                "refinement field `{field}` has {count} entries; limit is {MAX_REFINEMENT_SET_ENTRIES}"
            ),
            Self::MissingContentHash { field } => {
                write!(
                    f,
                    "refinement artifact `{field}` has an all-zero content hash"
                )
            }
            Self::MissingStateMap => f.write_str("deployment refinement has no state map"),
            Self::MissingObservationMap => {
                f.write_str("deployment refinement has no observation map")
            }
            Self::RelationSchemaMismatch { relation } => {
                write!(
                    f,
                    "{relation} relation does not match transition-system schemas"
                )
            }
            Self::UnitMismatch => f.write_str("source and target units are incompatible"),
            Self::FrameMismatch => f.write_str("source and target frames are incompatible"),
            Self::IncompatibleClocks => {
                f.write_str("source and target sampling clocks are incompatible")
            }
            Self::InvalidTimingEnvelope => f.write_str("latency/jitter envelope is invalid"),
            Self::InvalidNumericContract => {
                f.write_str("quantization/saturation contract is invalid")
            }
            Self::FaultOmission { fault } => {
                write!(
                    f,
                    "required fault `{fault}` is absent from the target model"
                )
            }
            Self::SafeStateMismatch => {
                f.write_str("source and target safe-state identities differ")
            }
            Self::InvalidBudget { field } => write!(f, "proof budget `{field}` is invalid"),
            Self::MissingCapability { capability } => {
                write!(f, "target lacks required capability `{capability}`")
            }
            Self::InvalidPermittedError => f.write_str("permitted error relation is invalid"),
            Self::InvalidBound { set, quantity } => {
                write!(f, "{set} bound `{quantity}` is non-finite or unordered")
            }
            Self::SchemaVersionMismatch { expected, found } => {
                write!(
                    f,
                    "manifest schema {found} does not match expected {expected}"
                )
            }
            Self::RetainedRootMismatch { declared, computed } => write!(
                f,
                "retained manifest root {declared:016x} does not match bytes root {computed:016x}"
            ),
            Self::RetainedManifestMismatch => {
                f.write_str("retained manifest belongs to another refinement problem")
            }
            Self::SourceIdentityDrift => f.write_str("source transition-system identity drifted"),
            Self::TargetIdentityDrift => f.write_str("deployed target/toolchain identity drifted"),
            Self::StateRelationDrift => f.write_str("state relation drifted"),
            Self::ObservationRelationDrift => f.write_str("observation relation drifted"),
            Self::UnitRelationDrift => f.write_str("unit relation drifted"),
            Self::FrameRelationDrift => f.write_str("frame relation drifted"),
            Self::TimingContractDrift => f.write_str("timing contract drifted"),
            Self::SaturationDrift => f.write_str("deployed saturation drifted"),
            Self::QuantizationDrift => f.write_str("deployed quantization drifted"),
            Self::PlantIdentityDrift => f.write_str("plant abstraction identity drifted"),
            Self::EnvironmentEnlargement { quantity } => {
                write!(f, "live environment enlarged at `{quantity}`")
            }
            Self::DisturbanceEnlargement { quantity } => {
                write!(f, "live disturbance set enlarged at `{quantity}`")
            }
            Self::SafeStateDrift => f.write_str("safe-state contract drifted"),
            Self::FaultModelDrift => f.write_str("fault model drifted"),
            Self::ClaimAssumptionDrift => f.write_str("refinement claim assumptions drifted"),
            Self::DuplicateProofAxis { axis } => write!(f, "duplicate {axis:?} proof evidence"),
            Self::MissingProofAxis { axis } => write!(f, "missing {axis:?} proof evidence"),
            Self::EvidenceManifestMismatch { axis } => {
                write!(f, "{axis:?} evidence belongs to another manifest")
            }
            Self::InvalidEvidenceHash { axis } => {
                write!(f, "{axis:?} evidence has an all-zero hash")
            }
            Self::UnknownProofAxis { axis } => write!(f, "{axis:?} proof is Unknown"),
            Self::RefutedProofAxis { axis } => write!(f, "{axis:?} proof is Refuted"),
            Self::MeasuredEvidenceIsNotUniversal { axis } => write!(
                f,
                "measured {axis:?} evidence cannot establish universal refinement"
            ),
        }
    }
}

impl std::error::Error for DeploymentRefinementError {}

fn validate_spec(spec: &DeploymentRefinementSpec) -> Result<(), DeploymentRefinementError> {
    let mut text_bytes = 0usize;
    validate_transition(&spec.source, "source", &mut text_bytes)?;
    validate_transition(&spec.target.system, "target", &mut text_bytes)?;
    validate_text("target triple", &spec.target.target_triple, &mut text_bytes)?;
    validate_text(
        "device revision",
        &spec.target.device_revision,
        &mut text_bytes,
    )?;
    validate_artifact(&spec.target.toolchain, "target toolchain", &mut text_bytes)?;
    validate_set_count("target capabilities", spec.target.capabilities.len())?;
    for capability in &spec.target.capabilities {
        validate_text("target capability", capability, &mut text_bytes)?;
    }

    validate_text(
        "source state relation schema",
        &spec.interface.source_state_schema,
        &mut text_bytes,
    )?;
    validate_text(
        "target state relation schema",
        &spec.interface.target_state_schema,
        &mut text_bytes,
    )?;
    if spec.interface.state_map_id.is_empty() {
        return Err(DeploymentRefinementError::MissingStateMap);
    }
    validate_text("state map", &spec.interface.state_map_id, &mut text_bytes)?;
    validate_text(
        "source observation relation schema",
        &spec.interface.source_observation_schema,
        &mut text_bytes,
    )?;
    validate_text(
        "target observation relation schema",
        &spec.interface.target_observation_schema,
        &mut text_bytes,
    )?;
    if spec.interface.observation_map_id.is_empty() {
        return Err(DeploymentRefinementError::MissingObservationMap);
    }
    validate_text(
        "observation map",
        &spec.interface.observation_map_id,
        &mut text_bytes,
    )?;
    for (field, value) in [
        ("source unit", &spec.interface.source_unit),
        ("target unit", &spec.interface.target_unit),
        ("source frame", &spec.interface.source_frame),
        ("target frame", &spec.interface.target_frame),
    ] {
        validate_text(field, value, &mut text_bytes)?;
    }
    if spec.interface.source_state_schema != spec.source.state_schema
        || spec.interface.target_state_schema != spec.target.system.state_schema
    {
        return Err(DeploymentRefinementError::RelationSchemaMismatch { relation: "state" });
    }
    if spec.interface.source_observation_schema != spec.source.observation_schema
        || spec.interface.target_observation_schema != spec.target.system.observation_schema
    {
        return Err(DeploymentRefinementError::RelationSchemaMismatch {
            relation: "observation",
        });
    }
    if spec.interface.source_unit != spec.interface.target_unit {
        return Err(DeploymentRefinementError::UnitMismatch);
    }
    if spec.interface.source_frame != spec.interface.target_frame {
        return Err(DeploymentRefinementError::FrameMismatch);
    }

    let source_period = spec.timing.source_period_ns;
    let target_period = spec.timing.target_period_ns;
    if source_period == 0 || target_period == 0 {
        return Err(DeploymentRefinementError::IncompatibleClocks);
    }
    let (larger, smaller) = if source_period >= target_period {
        (source_period, target_period)
    } else {
        (target_period, source_period)
    };
    if larger % smaller != 0 {
        return Err(DeploymentRefinementError::IncompatibleClocks);
    }
    if spec.timing.max_jitter_ns > spec.timing.max_latency_ns {
        return Err(DeploymentRefinementError::InvalidTimingEnvelope);
    }

    let NumericContract {
        quantization_step,
        saturation,
    } = spec.numeric;
    let saturation_width = saturation.hi - saturation.lo;
    if !quantization_step.is_finite()
        || quantization_step <= 0.0
        || !saturation.lo.is_finite()
        || !saturation.hi.is_finite()
        || saturation.lo > saturation.hi
        || !saturation_width.is_finite()
        || quantization_step > saturation_width
    {
        return Err(DeploymentRefinementError::InvalidNumericContract);
    }

    validate_artifact(&spec.plant, "plant", &mut text_bytes)?;
    validate_bounded_set(&spec.environment, "environment", &mut text_bytes)?;
    validate_bounded_set(&spec.disturbances, "disturbances", &mut text_bytes)?;

    validate_text("fault model", &spec.faults.model_id, &mut text_bytes)?;
    validate_set_count("required faults", spec.faults.required_faults.len())?;
    validate_set_count("modeled faults", spec.faults.modeled_faults.len())?;
    for fault in &spec.faults.required_faults {
        validate_text("required fault", fault, &mut text_bytes)?;
        if !spec.faults.modeled_faults.contains(fault) {
            return Err(DeploymentRefinementError::FaultOmission {
                fault: fault.clone(),
            });
        }
    }
    for fault in &spec.faults.modeled_faults {
        validate_text("modeled fault", fault, &mut text_bytes)?;
    }
    validate_text(
        "source safe state",
        &spec.faults.source_safe_state,
        &mut text_bytes,
    )?;
    validate_text(
        "target safe state",
        &spec.faults.target_safe_state,
        &mut text_bytes,
    )?;
    if spec.faults.source_safe_state != spec.faults.target_safe_state {
        return Err(DeploymentRefinementError::SafeStateMismatch);
    }

    if spec.horizon_steps == 0 {
        return Err(DeploymentRefinementError::InvalidBudget {
            field: "horizon steps",
        });
    }
    validate_text("invariant", &spec.invariant, &mut text_bytes)?;
    validate_text(
        "control objective",
        &spec.control_objective,
        &mut text_bytes,
    )?;
    if !spec.permitted_error.numeric_abs.is_finite() || spec.permitted_error.numeric_abs < 0.0 {
        return Err(DeploymentRefinementError::InvalidPermittedError);
    }
    validate_text(
        "functional error relation",
        &spec.permitted_error.functional_relation,
        &mut text_bytes,
    )?;
    validate_text(
        "safety error relation",
        &spec.permitted_error.safety_relation,
        &mut text_bytes,
    )?;

    for (field, value) in [
        ("max work units", spec.proof_budget.max_work_units),
        ("max memory bytes", spec.proof_budget.max_memory_bytes),
        ("max wall time", spec.proof_budget.max_wall_time_ns),
        (
            "cancellation poll stride",
            spec.proof_budget.cancellation_poll_stride,
        ),
    ] {
        if value == 0 {
            return Err(DeploymentRefinementError::InvalidBudget { field });
        }
    }
    if spec.proof_budget.cancellation_poll_stride > spec.proof_budget.max_work_units {
        return Err(DeploymentRefinementError::InvalidBudget {
            field: "cancellation poll stride",
        });
    }
    validate_set_count(
        "required capabilities",
        spec.proof_budget.required_capabilities.len(),
    )?;
    for capability in &spec.proof_budget.required_capabilities {
        validate_text("required capability", capability, &mut text_bytes)?;
        if !spec.target.capabilities.contains(capability) {
            return Err(DeploymentRefinementError::MissingCapability {
                capability: capability.clone(),
            });
        }
    }

    validate_set_count("assumptions", spec.assumptions.len())?;
    if spec.assumptions.is_empty() {
        return Err(DeploymentRefinementError::EmptyField {
            field: "assumptions",
        });
    }
    for (name, value) in &spec.assumptions {
        validate_text("assumption name", name, &mut text_bytes)?;
        validate_text("assumption value", value, &mut text_bytes)?;
    }
    Ok(())
}

fn validate_transition(
    system: &TransitionSystemIdentity,
    field: &'static str,
    text_bytes: &mut usize,
) -> Result<(), DeploymentRefinementError> {
    validate_artifact(&system.artifact, field, text_bytes)?;
    validate_text("state schema", &system.state_schema, text_bytes)?;
    validate_text("observation schema", &system.observation_schema, text_bytes)
}

fn validate_artifact(
    artifact: &ArtifactIdentity,
    field: &'static str,
    text_bytes: &mut usize,
) -> Result<(), DeploymentRefinementError> {
    validate_text(field, &artifact.name, text_bytes)?;
    validate_text("artifact version", &artifact.version, text_bytes)?;
    if artifact.content_hash.iter().all(|byte| *byte == 0) {
        return Err(DeploymentRefinementError::MissingContentHash { field });
    }
    Ok(())
}

fn validate_bounded_set(
    set: &BoundedSet,
    field: &'static str,
    text_bytes: &mut usize,
) -> Result<(), DeploymentRefinementError> {
    validate_text(field, &set.identity, text_bytes)?;
    validate_set_count(field, set.bounds.len())?;
    for (quantity, interval) in &set.bounds {
        validate_text("bounded quantity", quantity, text_bytes)?;
        if !interval.lo.is_finite() || !interval.hi.is_finite() || interval.lo > interval.hi {
            return Err(DeploymentRefinementError::InvalidBound {
                set: field,
                quantity: quantity.clone(),
            });
        }
    }
    Ok(())
}

fn validate_set_count(field: &'static str, count: usize) -> Result<(), DeploymentRefinementError> {
    if count > MAX_REFINEMENT_SET_ENTRIES {
        Err(DeploymentRefinementError::TooManyEntries { field, count })
    } else {
        Ok(())
    }
}

fn validate_text(
    field: &'static str,
    value: &str,
    text_bytes: &mut usize,
) -> Result<(), DeploymentRefinementError> {
    if value.is_empty() {
        return Err(DeploymentRefinementError::EmptyField { field });
    }
    *text_bytes = text_bytes.saturating_add(value.len());
    if *text_bytes > MAX_REFINEMENT_TEXT_BYTES {
        return Err(DeploymentRefinementError::TextBudgetExceeded {
            bytes: *text_bytes,
            limit: MAX_REFINEMENT_TEXT_BYTES,
        });
    }
    Ok(())
}

fn canonical_problem_bytes(spec: &DeploymentRefinementSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(4_096);
    push_bytes(&mut out, MANIFEST_DOMAIN);
    out.extend_from_slice(&DEPLOYMENT_REFINEMENT_SCHEMA_VERSION.to_le_bytes());
    push_transition(&mut out, &spec.source);
    push_transition(&mut out, &spec.target.system);
    push_str(&mut out, &spec.target.target_triple);
    push_str(&mut out, &spec.target.device_revision);
    push_artifact(&mut out, &spec.target.toolchain);
    push_string_set(&mut out, &spec.target.capabilities);
    push_str(&mut out, &spec.interface.source_state_schema);
    push_str(&mut out, &spec.interface.target_state_schema);
    push_str(&mut out, &spec.interface.state_map_id);
    push_str(&mut out, &spec.interface.source_observation_schema);
    push_str(&mut out, &spec.interface.target_observation_schema);
    push_str(&mut out, &spec.interface.observation_map_id);
    push_str(&mut out, &spec.interface.source_unit);
    push_str(&mut out, &spec.interface.target_unit);
    push_str(&mut out, &spec.interface.source_frame);
    push_str(&mut out, &spec.interface.target_frame);
    push_u64(&mut out, spec.timing.source_period_ns);
    push_u64(&mut out, spec.timing.target_period_ns);
    push_u64(&mut out, spec.timing.max_latency_ns);
    push_u64(&mut out, spec.timing.max_jitter_ns);
    push_f64(&mut out, spec.numeric.quantization_step);
    push_interval(&mut out, spec.numeric.saturation);
    push_artifact(&mut out, &spec.plant);
    push_bounded_set(&mut out, &spec.environment);
    push_bounded_set(&mut out, &spec.disturbances);
    push_str(&mut out, &spec.faults.model_id);
    push_string_set(&mut out, &spec.faults.required_faults);
    push_string_set(&mut out, &spec.faults.modeled_faults);
    push_str(&mut out, &spec.faults.source_safe_state);
    push_str(&mut out, &spec.faults.target_safe_state);
    push_u64(&mut out, spec.horizon_steps);
    push_str(&mut out, &spec.invariant);
    push_str(&mut out, &spec.control_objective);
    out.push(relation_tag(spec.relation));
    push_f64(&mut out, spec.permitted_error.numeric_abs);
    push_u64(&mut out, spec.permitted_error.temporal_ns);
    push_str(&mut out, &spec.permitted_error.functional_relation);
    push_str(&mut out, &spec.permitted_error.safety_relation);
    push_u64(&mut out, spec.proof_budget.max_work_units);
    push_u64(&mut out, spec.proof_budget.max_memory_bytes);
    push_u64(&mut out, spec.proof_budget.max_wall_time_ns);
    push_u64(&mut out, spec.proof_budget.cancellation_poll_stride);
    push_string_set(&mut out, &spec.proof_budget.required_capabilities);
    push_u64(&mut out, spec.assumptions.len() as u64);
    for (name, value) in &spec.assumptions {
        push_str(&mut out, name);
        push_str(&mut out, value);
    }
    out
}

fn push_artifact(out: &mut Vec<u8>, artifact: &ArtifactIdentity) {
    push_str(out, &artifact.name);
    push_str(out, &artifact.version);
    push_bytes(out, &artifact.content_hash);
}

fn push_transition(out: &mut Vec<u8>, system: &TransitionSystemIdentity) {
    push_artifact(out, &system.artifact);
    push_str(out, &system.state_schema);
    push_str(out, &system.observation_schema);
}

fn push_bounded_set(out: &mut Vec<u8>, set: &BoundedSet) {
    push_str(out, &set.identity);
    push_u64(out, set.bounds.len() as u64);
    for (name, interval) in &set.bounds {
        push_str(out, name);
        push_interval(out, *interval);
    }
}

fn push_string_set(out: &mut Vec<u8>, set: &BTreeSet<String>) {
    push_u64(out, set.len() as u64);
    for value in set {
        push_str(out, value);
    }
}

fn push_interval(out: &mut Vec<u8>, interval: Interval) {
    push_f64(out, interval.lo);
    push_f64(out, interval.hi);
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn push_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_str(out: &mut Vec<u8>, value: &str) {
    push_bytes(out, value.as_bytes());
}

fn push_bytes(out: &mut Vec<u8>, value: &[u8]) {
    push_u64(out, value.len() as u64);
    out.extend_from_slice(value);
}

const fn relation_tag(relation: RefinementRelation) -> u8 {
    match relation {
        RefinementRelation::TraceInclusion => 0,
        RefinementRelation::ApproximateSimulation => 1,
        RefinementRelation::ApproximateBisimulation => 2,
        RefinementRelation::RobustInvariant => 3,
        RefinementRelation::PerformanceBound => 4,
    }
}

fn same_f64(left: f64, right: f64) -> bool {
    left.to_bits() == right.to_bits()
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

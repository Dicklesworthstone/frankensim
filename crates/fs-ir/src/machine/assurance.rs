//! Admitted PR-4 operational and assurance semantics for Machine IR.
//!
//! This overlay binds one exact [`MachineBehaviorIdV1`] to structural sensor
//! and experiment lineage, decision-specific Context-of-Use/QoI links,
//! hazards and fault coverage, explicit accounting boundaries, and a finite
//! fidelity-escalation policy. External evidence, laws, monitors, budgets,
//! safety cases, balance audits, and fidelity models remain opaque,
//! content/version-bound references owned by their respective domains.
//!
//! Admission proves bounded reference closure, target typing, canonical
//! identity, explicit hazard ownership/monitoring/expiry, and that every
//! fidelity path terminates in refusal rather than cycling or silently
//! extrapolating. It does **not** authenticate experiments, validate physics,
//! certify safety or regulatory conformance, prove an accounting balance or
//! passivity, rank model accuracy, or authorize runtime model promotion.

use core::fmt;
use core::hash::{Hash, Hasher};
use core::num::NonZeroU64;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;
use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EntityId, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_evidence::vv::{
    AdmittedVvCase, ArtifactHeader, ArtifactId, ArtifactKind, ArtifactRef, AssumptionId,
    EvidenceAxisStatus, QoiId, UnitId,
};
use fs_qty::Dims;

use super::semantics::{AdmittedMachineBehavior, MachineBehaviorIdV1, StateSlotContract};
use super::{
    AdmittedMachineGraph, ClockId, FrameBinding, InterfaceBinding, InterfaceId, MachineElementId,
    MachineGraphIdV1, MachineIdError, ModelRef, PortId, RelationId, RelationMode, RelationSpec,
    StateSlotId, SubsystemId, TerminalCausality, TerminalId, TerminalQuantitySpec, TerminalShape,
    TerminalSpec,
};

/// Version of the Machine-IR operational-assurance overlay identity schema.
pub const MACHINE_ASSURANCE_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum sensor declarations in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_SENSORS: usize = 4_096;
/// Maximum external experiment bindings in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_EXPERIMENTS: usize = 4_096;
/// Maximum Context-of-Use bindings in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_CONTEXTS: usize = 1_024;
/// Maximum hazard declarations in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_HAZARDS: usize = 4_096;
/// Maximum fault declarations in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_FAULTS: usize = 4_096;
/// Maximum accounting windows in one assurance draft.
pub const MAX_MACHINE_ASSURANCE_ACCOUNTING_WINDOWS: usize = 4_096;
/// Maximum fidelity rungs in one assurance policy.
pub const MAX_MACHINE_ASSURANCE_FIDELITY_RUNGS: usize = 4_096;
/// Maximum aggregate nested references inspected by one admission.
pub const MAX_MACHINE_ASSURANCE_NESTED_REFERENCES: usize = 65_536;

const MACHINE_ASSURANCE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(32 * 1_024 * 1_024, 8 * 1_024 * 1_024, 16, 128_000, 65_536);

macro_rules! assurance_id {
    ($(#[$meta:meta])* $name:ident, $schema:ident, $identity:ident, $role:literal, $domain:literal, $context:literal) => {
        #[doc = concat!("Canonical schema marker for `", stringify!($name), "`.")]
        pub enum $schema {}

        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str = $domain;
            const NAME: &'static str = $role;
            const VERSION: u32 = MACHINE_ASSURANCE_SCHEMA_VERSION_V1;
            const CONTEXT: &'static str = $context;
            const FIELDS: &'static [FieldSpec] =
                &[FieldSpec::required("canonical-key", WireType::Utf8)];
        }

        #[doc = concat!("Typed durable digest for `", stringify!($name), "`.")]
        pub type $identity = EntityId<$schema>;

        $(#[$meta])*
        #[derive(Clone)]
        pub struct $name {
            canonical_key: Box<str>,
            receipt: IdentityReceipt<$identity>,
        }

        impl $name {
            /// Admit a canonical human-auditable key.
            ///
            /// # Errors
            /// Refuses noncanonical text or bounded identity publication.
            pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
                let key = key.into();
                super::validate_canonical_key($role, &key)?;
                let receipt = CanonicalEncoder::<$identity, _>::new(
                    super::MACHINE_IDENTITY_LIMITS,
                    NeverCancel,
                )?
                .utf8(Field::new(0, "canonical-key"), &key)?
                .finish()?;
                Ok(Self {
                    canonical_key: key.into_boxed_str(),
                    receipt,
                })
            }

            /// Canonical diagnostic and lowering key.
            #[must_use]
            pub fn canonical_key(&self) -> &str {
                &self.canonical_key
            }

            /// Domain-separated durable identity.
            #[must_use]
            pub const fn identity(&self) -> $identity {
                self.receipt.id()
            }

            /// Complete canonical-preimage receipt.
            #[must_use]
            pub const fn identity_receipt(&self) -> IdentityReceipt<$identity> {
                self.receipt
            }

            fn digest_bytes(&self) -> [u8; 32] {
                *self.identity().as_bytes()
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &Self) -> bool {
                self.identity() == other.identity()
            }
        }

        impl Eq for $name {}

        impl PartialOrd for $name {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Ord for $name {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                self.identity().cmp(&other.identity())
            }
        }

        impl Hash for $name {
            fn hash<H: Hasher>(&self, state: &mut H) {
                self.identity().hash(state);
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.canonical_key)
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("canonical_key", &self.canonical_key)
                    .field("identity", &self.identity())
                    .finish()
            }
        }
    };
}

assurance_id!(
    /// Durable identity of one declared sensor.
    SensorId,
    SensorIdSchemaV1,
    SensorEntityIdV1,
    "sensor-id",
    "org.frankensim.fs-ir.machine.sensor-id.v1",
    "one Machine-IR sensor independent of declaration and sampling order"
);
assurance_id!(
    /// Durable identity of one machine-to-experiment binding.
    ExperimentId,
    ExperimentIdSchemaV1,
    ExperimentEntityIdV1,
    "experiment-id",
    "org.frankensim.fs-ir.machine.experiment-id.v1",
    "one Machine-IR experiment binding independent of evidence serialization order"
);
assurance_id!(
    /// Durable safety-case hazard identity.
    HazardId,
    HazardIdSchemaV1,
    HazardEntityIdV1,
    "hazard-id",
    "org.frankensim.fs-ir.machine.hazard-id.v1",
    "one scoped hazard independent of fault-analysis and report order"
);
assurance_id!(
    /// Durable identity of one declared fault mode.
    FaultId,
    FaultIdSchemaV1,
    FaultEntityIdV1,
    "fault-id",
    "org.frankensim.fs-ir.machine.fault-id.v1",
    "one machine fault mode independent of injection and traversal order"
);
assurance_id!(
    /// Durable identity of one accounting window.
    AccountingWindowId,
    AccountingWindowIdSchemaV1,
    AccountingWindowEntityIdV1,
    "accounting-window-id",
    "org.frankensim.fs-ir.machine.accounting-window-id.v1",
    "one signed accounting boundary and audit window"
);
assurance_id!(
    /// Durable identity of one subsystem fidelity rung.
    FidelityRungId,
    FidelityRungIdSchemaV1,
    FidelityRungEntityIdV1,
    "fidelity-rung-id",
    "org.frankensim.fs-ir.machine.fidelity-rung-id.v1",
    "one subsystem implementation and admitted-use declaration"
);

macro_rules! assurance_ref {
    ($(#[$meta:meta])* $name:ident, $role:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            namespace: Box<str>,
            schema_version: NonZeroU64,
            semantic_digest: [u8; 32],
        }

        impl $name {
            /// Construct an opaque versioned semantic reference.
            ///
            /// # Errors
            /// Refuses a noncanonical namespace or all-zero digest.
            pub fn new(
                namespace: impl Into<String>,
                schema_version: NonZeroU64,
                semantic_digest: [u8; 32],
            ) -> Result<Self, super::MachineReferenceError> {
                let namespace = namespace.into();
                super::validate_canonical_key($role, &namespace)
                    .map_err(super::MachineReferenceError::Namespace)?;
                if semantic_digest == [0; 32] {
                    return Err(super::MachineReferenceError::ZeroDigest { role: $role });
                }
                Ok(Self {
                    namespace: namespace.into_boxed_str(),
                    schema_version,
                    semantic_digest,
                })
            }

            /// External schema namespace.
            #[must_use]
            pub fn namespace(&self) -> &str {
                &self.namespace
            }

            /// Explicit external schema version.
            #[must_use]
            pub const fn schema_version(&self) -> NonZeroU64 {
                self.schema_version
            }

            /// Exact semantic digest supplied by the external owner.
            #[must_use]
            pub const fn semantic_digest(&self) -> [u8; 32] {
                self.semantic_digest
            }

            fn append_canonical(&self, out: &mut Vec<u8>) {
                push_len_prefixed(out, self.namespace.as_bytes());
                out.extend_from_slice(&self.schema_version.get().to_le_bytes());
                out.extend_from_slice(&self.semantic_digest);
            }
        }
    };
}

assurance_ref!(
    /// External semantics for one sensor model.
    SensorModelRef,
    "sensor-model-ref"
);
assurance_ref!(
    /// Exact calibration evidence for one sensor.
    CalibrationRef,
    "sensor-calibration-ref"
);
assurance_ref!(
    /// External model bridging distinct observation clocks.
    SamplingBridgeRef,
    "sampling-bridge-ref"
);
assurance_ref!(
    /// Executable or declarative definition of one QoI.
    QoiDefinitionRef,
    "qoi-definition-ref"
);
assurance_ref!(
    /// Evidence relating an fs-evidence unit to a Machine quantity contract.
    UnitQuantityBridgeRef,
    "unit-quantity-bridge-ref"
);
assurance_ref!(
    /// Decision-specific resource and accuracy budget.
    DecisionBudgetRef,
    "decision-budget-ref"
);
assurance_ref!(
    /// External safety requirement governing one hazard.
    SafetyRequirementRef,
    "safety-requirement-ref"
);
assurance_ref!(
    /// Declared operating envelope for one hazard.
    OperatingEnvelopeRef,
    "operating-envelope-ref"
);
assurance_ref!(
    /// External safety-case artifact for one hazard.
    SafetyCaseRef,
    "safety-case-ref"
);
assurance_ref!(
    /// External semantics for one fault model.
    FaultModelRef,
    "fault-model-ref"
);
assurance_ref!(
    /// External containment semantics for one fault.
    FaultContainmentRef,
    "fault-containment-ref"
);
assurance_ref!(
    /// Exact fault-injection procedure or artifact.
    FaultInjectionRef,
    "fault-injection-ref"
);
assurance_ref!(
    /// Explicit evidence that no positive assurance claim is made.
    NoClaimRef,
    "assurance-no-claim-ref"
);
assurance_ref!(
    /// External definition of one accounting boundary.
    AccountingBoundaryRef,
    "accounting-boundary-ref"
);
assurance_ref!(
    /// External definition of one accounting interval.
    AccountingIntervalRef,
    "accounting-interval-ref"
);
assurance_ref!(
    /// Policy governing one accounting contribution or audit.
    AccountingPolicyRef,
    "accounting-policy-ref"
);
assurance_ref!(
    /// Unique external owner of one modeled loss contribution.
    LossOwnershipRef,
    "loss-ownership-ref"
);
assurance_ref!(
    /// External law for a species, element, or custom balance.
    BalanceLawRef,
    "balance-law-ref"
);
assurance_ref!(
    /// Applicability domain for one fidelity rung.
    ValidityDomainRef,
    "validity-domain-ref"
);
assurance_ref!(
    /// Cost and error model for one fidelity rung.
    CostErrorModelRef,
    "cost-error-model-ref"
);
assurance_ref!(
    /// State-transfer semantics for one escalation edge.
    StateTransferRef,
    "state-transfer-ref"
);
assurance_ref!(
    /// Exact semantic crosswalk between two model representations.
    ModelCrosswalkRef,
    "model-crosswalk-ref"
);
assurance_ref!(
    /// One executable applicability falsifier.
    FalsifierRef,
    "falsifier-ref"
);
assurance_ref!(
    /// Trigger semantics for one fidelity-policy edge.
    EscalationTriggerRef,
    "escalation-trigger-ref"
);
assurance_ref!(
    /// Identity-bound fixed-fidelity replay oracle.
    FixedReplayRef,
    "fixed-replay-ref"
);

/// A sensor observes one exact terminal or admitted state contract.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObservationTarget {
    /// Declared subsystem terminal.
    Terminal(TerminalId),
    /// Declared state slot with a PR-3 behavior contract.
    State(StateSlotId),
}

/// Relationship between the sensor clock and the observed target clock.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ObservationTiming {
    /// Sensor and target use the exact same declared logical clock.
    Direct,
    /// An external resampling/synchronization artifact bridges the clocks.
    ModeledResampling {
        /// Exact bridge between sensor and target clocks.
        bridge: SamplingBridgeRef,
    },
}

/// Structural sensor declaration; calibration authority remains external.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SensorSpec {
    /// Durable sensor identity.
    pub id: SensorId,
    /// Subsystem owning the sensor declaration.
    pub owner: SubsystemId,
    /// Exact observed machine object.
    pub target: ObservationTarget,
    /// Quantity contract exposed by the sensor.
    pub quantity: TerminalQuantitySpec,
    /// Value shape exposed by the sensor.
    pub shape: TerminalShape,
    /// Logical observation clock.
    pub clock: ClockId,
    /// Observation frame and orientation.
    pub frame: FrameBinding,
    /// Direct or explicitly bridged timing semantics.
    pub timing: ObservationTiming,
    /// External sensor-model semantics.
    pub model: SensorModelRef,
    /// Exact calibration evidence.
    pub calibration: CalibrationRef,
    /// Plant-visible or experiment-only exposure.
    pub exposure: SensorExposure,
}

/// Whether a sensor produces a plant signal or exists only in an experiment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SensorExposure {
    /// Sensor emits through one graph-visible output terminal.
    PlantSignal {
        /// Exact output terminal carrying the sensor signal.
        output: TerminalId,
    },
    /// Sensor exists only inside declared experiments.
    ExperimentOnly,
}

/// Exact bridge from a Machine sensor to an fs-evidence instrument identity.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SensorInstrumentBinding {
    /// Machine sensor used by the experiment.
    pub sensor: SensorId,
    /// Exact fs-evidence instrument identity.
    pub instrument: ArtifactId,
}

/// One machine-visible experiment artifact and its exact local dependencies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExperimentSpec {
    /// Durable local experiment-binding identity.
    pub id: ExperimentId,
    /// Exact admitted experiment artifact.
    pub artifact: ArtifactRef,
    /// Exact Context-of-Use governing the experiment.
    pub context: ArtifactRef,
    /// One-to-one machine-sensor to evidence-instrument map.
    pub instruments: Vec<SensorInstrumentBinding>,
    /// Exact QoIs observed by the experiment.
    pub qois: Vec<QoiId>,
}

/// Decision QoI target within the admitted machine.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum QoiTarget {
    /// QoI is defined from one declared sensor.
    Sensor(SensorId),
    /// QoI is defined directly from one terminal.
    Terminal(TerminalId),
    /// QoI is defined directly from one state slot.
    State(StateSlotId),
}

/// One typed input to an aggregate or pointwise QoI definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QoiInput {
    /// Exact machine source for this input.
    pub target: QoiTarget,
    /// Input quantity contract.
    pub quantity: TerminalQuantitySpec,
    /// Input value shape.
    pub shape: TerminalShape,
}

/// Machine binding for one externally named Context-of-Use QoI.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QoiBinding {
    /// Exact fs-evidence QoI identity.
    pub id: QoiId,
    /// Nonempty canonical input set; aggregate QoIs may depend on many inputs.
    pub inputs: Vec<QoiInput>,
    /// Unit declared by the admitted Context-of-Use.
    pub unit: UnitId,
    /// Executable/integral/query semantics remain externally owned.
    pub definition: QoiDefinitionRef,
    /// External evidence for any UnitId-to-Machine-quantity interpretation.
    pub unit_bridge: UnitQuantityBridgeRef,
}

/// Context-qualified QoI identity used outside a single Context binding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ContextQoiKey {
    /// Context-of-Use artifact identity.
    pub context: ArtifactId,
    /// QoI identity within that context.
    pub qoi: QoiId,
}

/// Exact Context-of-Use and validation-plan binding for machine decisions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextBinding {
    /// Exact admitted Context-of-Use artifact.
    pub context: ArtifactRef,
    /// Exact admitted validation-plan artifact.
    pub validation_plan: ArtifactRef,
    /// Complete machine bindings for the context QoIs.
    pub qois: Vec<QoiBinding>,
    /// Decision-specific budget semantics.
    pub budget: DecisionBudgetRef,
}

/// Whether fault analysis structurally covers a hazard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaultCoverage {
    /// At least one admitted fault must reference the hazard.
    Modeled,
    /// The missing fault model is explicit and carries no assurance authority.
    Unmodeled(NoClaimRef),
}

/// Scoped safety hazard with all mandatory assumption-lifecycle fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HazardSpec {
    /// Durable hazard identity.
    pub id: HazardId,
    /// Context-of-Use in which the hazard is relevant.
    pub context: ArtifactRef,
    /// Nonempty exact machine scope.
    pub scope: Vec<MachineScope>,
    /// Governing external safety requirement.
    pub requirement: SafetyRequirementRef,
    /// Exact operating envelope.
    pub operating_envelope: OperatingEnvelopeRef,
    /// Exact external safety-case artifact.
    pub safety_case: SafetyCaseRef,
    /// Exact rows in the context case's admitted assumptions ledger.
    pub assumptions: Vec<AssumptionId>,
    /// Modeled or explicit no-claim fault coverage.
    pub fault_coverage: FaultCoverage,
}

/// One fault mode, affected graph elements, and linked hazards.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaultSpec {
    /// Durable fault identity.
    pub id: FaultId,
    /// Nonempty affected machine scope.
    pub affected: Vec<MachineScope>,
    /// Hazards structurally covered by this fault.
    pub hazards: Vec<HazardId>,
    /// External fault-model semantics.
    pub model: FaultModelRef,
    /// External containment semantics.
    pub containment: FaultContainmentRef,
    /// Exact fault-injection semantics.
    pub injection: FaultInjectionRef,
}

/// Scope vocabulary shared by hazards and faults.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineScope {
    /// Entire admitted machine graph.
    WholeMachine,
    /// One exact subsystem.
    Subsystem(SubsystemId),
    /// One durable owned machine element.
    Element(MachineElementId),
    /// One exact relation.
    Relation(RelationId),
    /// One exact interface.
    Interface(InterfaceId),
}

/// Graph object participating in one explicit accounting statement.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AccountingTarget {
    /// Machine relation providing storage, dissipation, or a source term.
    Relation(RelationId),
    /// Effort/flow port crossing an accounting boundary.
    Port(PortId),
    /// Role-oriented interface crossing an accounting boundary.
    Interface(InterfaceId),
    /// State slot carrying an explicit stored quantity.
    State(StateSlotId),
}

/// Meaning of the extensive balance audited by one accounting window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BalanceKind {
    /// Energy balance.
    Energy,
    /// Enthalpy balance.
    Enthalpy,
    /// Linear-momentum balance.
    LinearMomentum,
    /// Angular-momentum balance.
    AngularMomentum,
    /// Mass balance.
    Mass,
    /// Electric-charge balance.
    ElectricCharge,
    /// Amount-of-substance balance.
    AmountOfSubstance,
    /// Species balance under an external species law.
    Species(BalanceLawRef),
    /// Element balance under an external element law.
    Elements(BalanceLawRef),
    /// Entropy balance.
    Entropy,
    /// Exergy balance.
    Exergy,
    /// Custom extensive balance under an external law.
    Custom(BalanceLawRef),
}

/// Explicit sign law for one accounting contribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AccountingOrientation {
    /// Positive values increase stored quantity.
    StoredIncreasePositive,
    /// Loss contributions are nonnegative.
    NonnegativeLoss,
    /// Boundary inflow is positive.
    IntoBoundaryPositive,
    /// Boundary outflow is positive.
    OutOfBoundaryPositive,
}

/// Structural contribution role; this is not itself a conservation proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AccountingRole {
    /// Stored extensive quantity.
    Storage,
    /// Irreversible loss or dissipation.
    Dissipation,
    /// Source explicitly included inside the balance.
    IncludedSource,
    /// Exchange with the external environment.
    ExternalExchange,
    /// Advective or material stream crossing the boundary.
    Stream,
}

/// One explicit contribution in an accounting window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccountingEntry {
    /// Exact graph object contributing to the statement.
    pub target: AccountingTarget,
    /// Structural accounting role.
    pub role: AccountingRole,
    /// Explicit sign convention.
    pub orientation: AccountingOrientation,
    /// External contribution policy.
    pub policy: AccountingPolicyRef,
    /// Unique loss owner, required only for dissipation.
    pub loss_ownership: Option<LossOwnershipRef>,
}

/// Signed accounting boundary and finite logical audit window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccountingWindow {
    /// Durable accounting-window identity.
    pub id: AccountingWindowId,
    /// Context-of-Use governing the statement.
    pub context: ArtifactRef,
    /// Logical clock used by the interval.
    pub clock: ClockId,
    /// Extensive balance being audited.
    pub balance: BalanceKind,
    /// Quantity and dimensions of the balance.
    pub quantity: TerminalQuantitySpec,
    /// Exact external boundary definition.
    pub boundary: AccountingBoundaryRef,
    /// Exact finite audit interval.
    pub interval: AccountingIntervalRef,
    /// Nonempty signed contribution list.
    pub entries: Vec<AccountingEntry>,
    /// External policy governing the complete audit.
    pub audit_policy: AccountingPolicyRef,
}

/// One fidelity implementation with explicit applicability and falsifiers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FidelityRung {
    /// Durable fidelity-rung identity.
    pub id: FidelityRungId,
    /// Subsystem implemented by this rung.
    pub subsystem: SubsystemId,
    /// Same nominal model-reference family used by the base machine graph.
    pub model: ModelRef,
    /// Exact relation between this representation and the graph model.
    pub model_crosswalk: ModelCrosswalkRef,
    /// Declared applicability domain.
    pub validity_domain: ValidityDomainRef,
    /// Declared cost and error model.
    pub cost_error_model: CostErrorModelRef,
    /// Nonempty applicability-falsifier set.
    pub falsifiers: Vec<FalsifierRef>,
    /// Context-qualified QoIs preserved by the rung.
    pub qois: Vec<ContextQoiKey>,
}

/// Required action after a fidelity-rung applicability trigger fires.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EscalationAction {
    /// Move through the admitted acyclic policy graph.
    Escalate {
        /// Higher-fidelity target rung.
        target: FidelityRungId,
        /// Exact state-transfer semantics.
        transfer: StateTransferRef,
        /// Exact source-to-target model crosswalk.
        crosswalk: ModelCrosswalkRef,
    },
    /// Refuse rather than silently extrapolate.
    Refuse(NoClaimRef),
}

/// One total outgoing transition from a fidelity rung.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EscalationSpec {
    /// Source rung with exactly one outgoing action.
    pub from: FidelityRungId,
    /// Applicability trigger selecting the action.
    pub trigger: EscalationTriggerRef,
    /// Escalation or explicit refusal.
    pub action: EscalationAction,
}

/// Complete fixed-baseline and terminating escalation policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FidelityPolicy {
    /// Exactly one baseline rung for every graph subsystem.
    pub baselines: Vec<FidelityRungId>,
    /// Complete bounded rung set.
    pub rungs: Vec<FidelityRung>,
    /// Exactly one outgoing transition for every rung.
    pub escalations: Vec<EscalationSpec>,
    /// External replay oracle retained as the non-adaptive baseline.
    pub fixed_replay: FixedReplayRef,
}

/// Exact fs-evidence schema-admission receipt bound into Machine-IR identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VvCaseBinding {
    /// Exact admitted Context-of-Use artifact.
    pub context: ArtifactRef,
    /// Bound fs-evidence schema version.
    pub schema_version: u32,
    /// Bound fs-evidence ruleset version.
    pub ruleset_version: u32,
    /// Canonical V&V case hash.
    pub case_hash: ContentHash,
    /// Exact schema-admission receipt hash.
    pub receipt_hash: ContentHash,
}

/// Canonical identity schema for one admitted operational-assurance overlay.
pub enum MachineAssuranceIdentitySchemaV1 {}

impl CanonicalSchema for MachineAssuranceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.assurance.v1";
    const NAME: &'static str = "admitted-machine-assurance";
    const VERSION: u32 = MACHINE_ASSURANCE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "base behavior, sensors and experiments, Context-of-Use QoIs, hazards and faults, accounting boundaries, and terminating fidelity policy";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("assurance-schema-version", WireType::U64),
        FieldSpec::required("base-machine-graph", WireType::Bytes),
        FieldSpec::required("base-machine-behavior", WireType::Bytes),
        FieldSpec::required("sensors", WireType::OrderedBytes),
        FieldSpec::required("experiments", WireType::OrderedBytes),
        FieldSpec::required("contexts", WireType::OrderedBytes),
        FieldSpec::required("vv-case-receipts", WireType::OrderedBytes),
        FieldSpec::required("hazards", WireType::OrderedBytes),
        FieldSpec::required("faults", WireType::OrderedBytes),
        FieldSpec::required("accounting-windows", WireType::OrderedBytes),
        FieldSpec::required("fidelity-baselines", WireType::OrderedBytes),
        FieldSpec::required("fidelity-rungs", WireType::OrderedBytes),
        FieldSpec::required("fidelity-escalations", WireType::OrderedBytes),
        FieldSpec::required("fixed-fidelity-replay", WireType::Bytes),
    ];
}

/// Strong semantic identity of an admitted operational-assurance overlay.
pub type MachineAssuranceIdV1 = ProblemSemanticId<MachineAssuranceIdentitySchemaV1>;

/// Closed rule vocabulary for deterministic assurance diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum MachineAssuranceRule {
    /// A public collection or aggregate nested-reference limit was exceeded.
    ResourceLimit = 1,
    /// The behavior overlay names a different base graph.
    BehaviorGraphMismatch = 2,
    /// No sensors were declared.
    EmptySensors = 3,
    /// One sensor identity was declared more than once.
    DuplicateSensor = 4,
    /// A sensor target is absent from the admitted graph or behavior.
    UnknownSensorTarget = 5,
    /// A sensor references an unknown logical clock.
    UnknownSensorClock = 6,
    /// A sensor quantity differs from its target contract.
    SensorQuantityGap = 7,
    /// A sensor shape differs from its target contract.
    SensorShapeGap = 8,
    /// Direct sensor timing differs from the target clock.
    SensorClockGap = 9,
    /// A sensor frame differs from its target contract.
    SensorFrameGap = 10,
    /// A sensor uses an inadmissible quantity form.
    UnsupportedSensorQuantity = 11,
    /// No Context-of-Use bindings were declared.
    EmptyContexts = 12,
    /// One Context-of-Use was bound more than once.
    DuplicateContext = 13,
    /// A Context-of-Use reference is wrong-kind or stale.
    InvalidContextArtifact = 14,
    /// A validation-plan reference is wrong-kind or stale.
    InvalidValidationPlanArtifact = 15,
    /// A Context-of-Use has no machine QoI bindings.
    EmptyContextQois = 16,
    /// One Context-qualified QoI was bound more than once.
    DuplicateQoi = 17,
    /// A QoI target is absent from the machine.
    UnknownQoiTarget = 18,
    /// A QoI input uses an inadmissible quantity form.
    UnsupportedQoiQuantity = 19,
    /// A QoI input quantity differs from its target.
    QoiQuantityGap = 20,
    /// A QoI input shape differs from its target.
    QoiShapeGap = 21,
    /// No experiment bindings were declared.
    EmptyExperiments = 22,
    /// One local experiment identity was declared more than once.
    DuplicateExperiment = 23,
    /// An experiment reference is wrong-kind or stale.
    InvalidExperimentArtifact = 24,
    /// An experiment references an unknown Context-of-Use.
    UnknownExperimentContext = 25,
    /// An experiment declares no sensor/instrument bindings.
    EmptyExperimentSensors = 26,
    /// One experiment binds a machine sensor more than once.
    DuplicateExperimentSensor = 27,
    /// An experiment references an unknown machine sensor.
    UnknownExperimentSensor = 28,
    /// An experiment declares no QoIs.
    EmptyExperimentQois = 29,
    /// An experiment repeats a QoI identity.
    DuplicateExperimentQoi = 30,
    /// An experiment references an unknown Context-qualified QoI.
    UnknownExperimentQoi = 31,
    /// A declared machine QoI has no experiment evidence binding.
    MissingQoiExperiment = 32,
    /// No hazards were declared.
    EmptyHazards = 33,
    /// One hazard identity was declared more than once.
    DuplicateHazard = 34,
    /// A hazard references an unknown Context-of-Use.
    UnknownHazardContext = 35,
    /// A hazard has no explicit machine scope.
    EmptyHazardScope = 36,
    /// A hazard repeats a scope element.
    DuplicateHazardElement = 37,
    /// A hazard scope element is absent from the machine.
    UnknownHazardElement = 38,
    /// One fault identity was declared more than once.
    DuplicateFault = 40,
    /// A fault has no affected machine scope.
    EmptyFaultAffected = 41,
    /// A fault repeats an affected scope element.
    DuplicateFaultElement = 42,
    /// A fault scope element is absent from the machine.
    UnknownFaultElement = 43,
    /// A fault covers no declared hazard.
    EmptyFaultHazards = 44,
    /// A fault repeats a hazard edge.
    DuplicateFaultHazard = 45,
    /// A fault references an unknown hazard.
    UnknownFaultHazard = 46,
    /// A modeled hazard has no covering fault.
    UncoveredModeledHazard = 47,
    /// An explicitly unmodeled hazard also has a fault edge.
    ContradictoryFaultCoverage = 48,
    /// No accounting windows were declared.
    EmptyAccountingWindows = 49,
    /// One accounting-window identity was declared more than once.
    DuplicateAccountingWindow = 50,
    /// An accounting window references an unknown context.
    UnknownAccountingContext = 51,
    /// An accounting window references an unknown clock.
    UnknownAccountingClock = 52,
    /// An accounting window uses an inadmissible quantity form.
    UnsupportedAccountingQuantity = 53,
    /// An accounting window has no signed contributions.
    EmptyAccountingEntries = 54,
    /// One accounting target appears more than once in a window.
    DuplicateAccountingEntry = 55,
    /// An accounting target is absent from the machine.
    UnknownAccountingTarget = 56,
    /// Role, target, and orientation are inconsistent.
    InvalidAccountingRole = 57,
    /// The fidelity policy has no rungs.
    EmptyFidelityRungs = 58,
    /// One fidelity-rung identity was declared more than once.
    DuplicateFidelityRung = 59,
    /// A fidelity rung references an unknown subsystem.
    UnknownFidelitySubsystem = 60,
    /// A fidelity rung declares no applicability falsifier.
    EmptyFidelityFalsifiers = 61,
    /// A fidelity rung repeats a falsifier.
    DuplicateFidelityFalsifier = 62,
    /// A fidelity rung declares no Context-qualified QoIs.
    EmptyFidelityQois = 63,
    /// A fidelity rung repeats a Context-qualified QoI.
    DuplicateFidelityQoi = 64,
    /// A fidelity rung references an unknown QoI.
    UnknownFidelityQoi = 65,
    /// A declared machine QoI appears on no fidelity rung.
    MissingQoiFidelityCoverage = 66,
    /// One baseline rung was listed more than once.
    DuplicateFidelityBaseline = 67,
    /// A baseline references an unknown rung.
    UnknownFidelityBaseline = 68,
    /// One subsystem has multiple baseline rungs.
    DuplicateSubsystemBaseline = 69,
    /// A graph subsystem has no baseline rung.
    MissingSubsystemBaseline = 70,
    /// A fidelity rung has multiple outgoing actions.
    DuplicateEscalation = 71,
    /// A fidelity rung has no outgoing action.
    MissingEscalation = 72,
    /// An escalation references an unknown rung.
    UnknownEscalationRung = 73,
    /// An escalation crosses subsystem ownership.
    CrossSubsystemEscalation = 74,
    /// An escalation targets its source rung.
    SelfEscalation = 75,
    /// The fidelity escalation graph contains a cycle.
    FidelityEscalationCycle = 76,
    /// An exact external artifact reference failed closure.
    InvalidArtifactReference = 77,
    /// Canonical identity publication failed.
    Identity = 78,
    /// A bound Context-of-Use has no supplied admitted V&V case.
    MissingVvCase = 79,
    /// A supplied V&V case has no machine Context binding.
    UnexpectedVvCase = 80,
    /// A supplied V&V receipt does not verify its exact case.
    InvalidVvCaseReceipt = 81,
    /// Machine and admitted-context QoI sets differ.
    ContextQoiSetMismatch = 82,
    /// A QoI definition has no machine inputs.
    EmptyQoiInputs = 83,
    /// A QoI definition repeats an input target.
    DuplicateQoiInput = 84,
    /// A subsystem baseline model differs from its graph model.
    BaselineModelMismatch = 85,
    /// Machine and admitted-experiment QoI sets differ.
    ExperimentQoiSetMismatch = 86,
    /// A sensor owner is absent from the graph.
    UnknownSensorOwner = 87,
    /// A plant-signal sensor names an unknown output terminal.
    UnknownSensorOutput = 88,
    /// A sensor and plant-signal output have different owners.
    SensorOutputOwnerMismatch = 89,
    /// A plant-signal sensor names a non-output terminal.
    SensorOutputCausalityGap = 90,
    /// A plant-signal output contract differs from its sensor.
    SensorOutputContractGap = 91,
    /// An experiment instrument is absent from admitted evidence.
    UnknownExperimentInstrument = 92,
    /// One experiment binds an instrument more than once.
    DuplicateExperimentInstrument = 93,
    /// Machine and admitted-experiment instrument sets differ.
    ExperimentInstrumentSetMismatch = 94,
    /// A hazard has no admitted assumption links.
    EmptyHazardAssumptions = 95,
    /// A hazard repeats an assumption link.
    DuplicateHazardAssumption = 96,
    /// A hazard assumption is absent from its admitted case.
    UnknownHazardAssumption = 97,
    /// A dissipative contribution has no loss owner.
    MissingLossOwnership = 98,
    /// A non-dissipative contribution declares a loss owner.
    UnexpectedLossOwnership = 99,
    /// One loss owner claims multiple contributions in a window.
    DuplicateLossOwnership = 100,
    /// A fidelity rung is unreachable from its subsystem baseline.
    UnreachableFidelityRung = 101,
    /// Multiple sensors claim the same plant output terminal.
    DuplicateSensorOutput = 102,
    /// A machine QoI unit differs from its admitted-context unit.
    QoiUnitMismatch = 103,
    /// An admitted evidence experiment has no local binding.
    MissingExperimentBinding = 104,
    /// A known balance kind has incompatible quantity dimensions.
    AccountingBalanceQuantityGap = 105,
    /// One evidence artifact has multiple local experiment identities.
    DuplicateExperimentArtifactBinding = 106,
    /// Two accounting targets claim the same role-qualified graph atom.
    OverlappingAccountingTarget = 107,
    /// An escalation target drops a source-rung QoI.
    FidelityQoiDrop = 108,
}

impl MachineAssuranceRule {
    /// Stable machine-readable rule code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ResourceLimit => "MachineAssuranceResourceLimit",
            Self::BehaviorGraphMismatch => "MachineAssuranceBehaviorGraphMismatch",
            Self::EmptySensors => "MachineAssuranceEmptySensors",
            Self::DuplicateSensor => "MachineAssuranceDuplicateSensor",
            Self::UnknownSensorTarget => "MachineAssuranceUnknownSensorTarget",
            Self::UnknownSensorClock => "MachineAssuranceUnknownSensorClock",
            Self::SensorQuantityGap => "MachineAssuranceSensorQuantityGap",
            Self::SensorShapeGap => "MachineAssuranceSensorShapeGap",
            Self::SensorClockGap => "MachineAssuranceSensorClockGap",
            Self::SensorFrameGap => "MachineAssuranceSensorFrameGap",
            Self::UnsupportedSensorQuantity => "MachineAssuranceUnsupportedSensorQuantity",
            Self::EmptyContexts => "MachineAssuranceEmptyContexts",
            Self::DuplicateContext => "MachineAssuranceDuplicateContext",
            Self::InvalidContextArtifact => "MachineAssuranceInvalidContextArtifact",
            Self::InvalidValidationPlanArtifact => "MachineAssuranceInvalidValidationPlanArtifact",
            Self::EmptyContextQois => "MachineAssuranceEmptyContextQois",
            Self::DuplicateQoi => "MachineAssuranceDuplicateQoi",
            Self::UnknownQoiTarget => "MachineAssuranceUnknownQoiTarget",
            Self::UnsupportedQoiQuantity => "MachineAssuranceUnsupportedQoiQuantity",
            Self::QoiQuantityGap => "MachineAssuranceQoiQuantityGap",
            Self::QoiShapeGap => "MachineAssuranceQoiShapeGap",
            Self::EmptyExperiments => "MachineAssuranceEmptyExperiments",
            Self::DuplicateExperiment => "MachineAssuranceDuplicateExperiment",
            Self::InvalidExperimentArtifact => "MachineAssuranceInvalidExperimentArtifact",
            Self::UnknownExperimentContext => "MachineAssuranceUnknownExperimentContext",
            Self::EmptyExperimentSensors => "MachineAssuranceEmptyExperimentSensors",
            Self::DuplicateExperimentSensor => "MachineAssuranceDuplicateExperimentSensor",
            Self::UnknownExperimentSensor => "MachineAssuranceUnknownExperimentSensor",
            Self::EmptyExperimentQois => "MachineAssuranceEmptyExperimentQois",
            Self::DuplicateExperimentQoi => "MachineAssuranceDuplicateExperimentQoi",
            Self::UnknownExperimentQoi => "MachineAssuranceUnknownExperimentQoi",
            Self::MissingQoiExperiment => "MachineAssuranceMissingQoiExperiment",
            Self::EmptyHazards => "MachineAssuranceEmptyHazards",
            Self::DuplicateHazard => "MachineAssuranceDuplicateHazard",
            Self::UnknownHazardContext => "MachineAssuranceUnknownHazardContext",
            Self::EmptyHazardScope => "MachineAssuranceEmptyHazardScope",
            Self::DuplicateHazardElement => "MachineAssuranceDuplicateHazardElement",
            Self::UnknownHazardElement => "MachineAssuranceUnknownHazardElement",
            Self::DuplicateFault => "MachineAssuranceDuplicateFault",
            Self::EmptyFaultAffected => "MachineAssuranceEmptyFaultAffected",
            Self::DuplicateFaultElement => "MachineAssuranceDuplicateFaultElement",
            Self::UnknownFaultElement => "MachineAssuranceUnknownFaultElement",
            Self::EmptyFaultHazards => "MachineAssuranceEmptyFaultHazards",
            Self::DuplicateFaultHazard => "MachineAssuranceDuplicateFaultHazard",
            Self::UnknownFaultHazard => "MachineAssuranceUnknownFaultHazard",
            Self::UncoveredModeledHazard => "MachineAssuranceUncoveredModeledHazard",
            Self::ContradictoryFaultCoverage => "MachineAssuranceContradictoryFaultCoverage",
            Self::EmptyAccountingWindows => "MachineAssuranceEmptyAccountingWindows",
            Self::DuplicateAccountingWindow => "MachineAssuranceDuplicateAccountingWindow",
            Self::UnknownAccountingContext => "MachineAssuranceUnknownAccountingContext",
            Self::UnknownAccountingClock => "MachineAssuranceUnknownAccountingClock",
            Self::UnsupportedAccountingQuantity => "MachineAssuranceUnsupportedAccountingQuantity",
            Self::EmptyAccountingEntries => "MachineAssuranceEmptyAccountingEntries",
            Self::DuplicateAccountingEntry => "MachineAssuranceDuplicateAccountingEntry",
            Self::UnknownAccountingTarget => "MachineAssuranceUnknownAccountingTarget",
            Self::InvalidAccountingRole => "MachineAssuranceInvalidAccountingRole",
            Self::EmptyFidelityRungs => "MachineAssuranceEmptyFidelityRungs",
            Self::DuplicateFidelityRung => "MachineAssuranceDuplicateFidelityRung",
            Self::UnknownFidelitySubsystem => "MachineAssuranceUnknownFidelitySubsystem",
            Self::EmptyFidelityFalsifiers => "MachineAssuranceEmptyFidelityFalsifiers",
            Self::DuplicateFidelityFalsifier => "MachineAssuranceDuplicateFidelityFalsifier",
            Self::EmptyFidelityQois => "MachineAssuranceEmptyFidelityQois",
            Self::DuplicateFidelityQoi => "MachineAssuranceDuplicateFidelityQoi",
            Self::UnknownFidelityQoi => "MachineAssuranceUnknownFidelityQoi",
            Self::MissingQoiFidelityCoverage => "MachineAssuranceMissingQoiFidelityCoverage",
            Self::DuplicateFidelityBaseline => "MachineAssuranceDuplicateFidelityBaseline",
            Self::UnknownFidelityBaseline => "MachineAssuranceUnknownFidelityBaseline",
            Self::DuplicateSubsystemBaseline => "MachineAssuranceDuplicateSubsystemBaseline",
            Self::MissingSubsystemBaseline => "MachineAssuranceMissingSubsystemBaseline",
            Self::DuplicateEscalation => "MachineAssuranceDuplicateEscalation",
            Self::MissingEscalation => "MachineAssuranceMissingEscalation",
            Self::UnknownEscalationRung => "MachineAssuranceUnknownEscalationRung",
            Self::CrossSubsystemEscalation => "MachineAssuranceCrossSubsystemEscalation",
            Self::SelfEscalation => "MachineAssuranceSelfEscalation",
            Self::FidelityEscalationCycle => "MachineAssuranceFidelityEscalationCycle",
            Self::InvalidArtifactReference => "MachineAssuranceInvalidArtifactReference",
            Self::Identity => "MachineAssuranceIdentity",
            Self::MissingVvCase => "MachineAssuranceMissingVvCase",
            Self::UnexpectedVvCase => "MachineAssuranceUnexpectedVvCase",
            Self::InvalidVvCaseReceipt => "MachineAssuranceInvalidVvCaseReceipt",
            Self::ContextQoiSetMismatch => "MachineAssuranceContextQoiSetMismatch",
            Self::EmptyQoiInputs => "MachineAssuranceEmptyQoiInputs",
            Self::DuplicateQoiInput => "MachineAssuranceDuplicateQoiInput",
            Self::BaselineModelMismatch => "MachineAssuranceBaselineModelMismatch",
            Self::ExperimentQoiSetMismatch => "MachineAssuranceExperimentQoiSetMismatch",
            Self::UnknownSensorOwner => "MachineAssuranceUnknownSensorOwner",
            Self::UnknownSensorOutput => "MachineAssuranceUnknownSensorOutput",
            Self::SensorOutputOwnerMismatch => "MachineAssuranceSensorOutputOwnerMismatch",
            Self::SensorOutputCausalityGap => "MachineAssuranceSensorOutputCausalityGap",
            Self::SensorOutputContractGap => "MachineAssuranceSensorOutputContractGap",
            Self::UnknownExperimentInstrument => "MachineAssuranceUnknownExperimentInstrument",
            Self::DuplicateExperimentInstrument => "MachineAssuranceDuplicateExperimentInstrument",
            Self::ExperimentInstrumentSetMismatch => {
                "MachineAssuranceExperimentInstrumentSetMismatch"
            }
            Self::EmptyHazardAssumptions => "MachineAssuranceEmptyHazardAssumptions",
            Self::DuplicateHazardAssumption => "MachineAssuranceDuplicateHazardAssumption",
            Self::UnknownHazardAssumption => "MachineAssuranceUnknownHazardAssumption",
            Self::MissingLossOwnership => "MachineAssuranceMissingLossOwnership",
            Self::UnexpectedLossOwnership => "MachineAssuranceUnexpectedLossOwnership",
            Self::DuplicateLossOwnership => "MachineAssuranceDuplicateLossOwnership",
            Self::UnreachableFidelityRung => "MachineAssuranceUnreachableFidelityRung",
            Self::DuplicateSensorOutput => "MachineAssuranceDuplicateSensorOutput",
            Self::QoiUnitMismatch => "MachineAssuranceQoiUnitMismatch",
            Self::MissingExperimentBinding => "MachineAssuranceMissingExperimentBinding",
            Self::AccountingBalanceQuantityGap => "MachineAssuranceAccountingBalanceQuantityGap",
            Self::DuplicateExperimentArtifactBinding => {
                "MachineAssuranceDuplicateExperimentArtifactBinding"
            }
            Self::OverlappingAccountingTarget => "MachineAssuranceOverlappingAccountingTarget",
            Self::FidelityQoiDrop => "MachineAssuranceFidelityQoiDrop",
        }
    }
}

/// Stable diagnostic subject for assurance admission.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum MachineAssuranceSubject {
    /// Entire assurance overlay.
    Assurance,
    /// Base machine graph.
    Graph(MachineGraphIdV1),
    /// Sensor declaration.
    Sensor(SensorId),
    /// Local experiment binding.
    Experiment(ExperimentId),
    /// Exact Context-of-Use reference.
    Context(ArtifactRef),
    /// External artifact identity.
    Artifact(ArtifactId),
    /// Context-qualified QoI.
    Qoi(ContextQoiKey),
    /// V&V assumption-ledger row.
    Assumption(AssumptionId),
    /// Hazard declaration.
    Hazard(HazardId),
    /// Fault declaration.
    Fault(FaultId),
    /// Durable machine element.
    Element(MachineElementId),
    /// Logical machine clock.
    Clock(ClockId),
    /// Accounting-window declaration.
    AccountingWindow(AccountingWindowId),
    /// Accounting contribution target.
    AccountingTarget(AccountingTarget),
    /// Fidelity-rung declaration.
    FidelityRung(FidelityRungId),
    /// Machine subsystem.
    Subsystem(SubsystemId),
}

/// One deterministic assurance-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct MachineAssuranceFinding {
    rule: MachineAssuranceRule,
    subject: MachineAssuranceSubject,
    related: Option<MachineAssuranceSubject>,
}

impl MachineAssuranceFinding {
    fn new(
        rule: MachineAssuranceRule,
        subject: MachineAssuranceSubject,
        related: Option<MachineAssuranceSubject>,
    ) -> Self {
        Self {
            rule,
            subject,
            related,
        }
    }

    /// Stable rule that produced this finding.
    #[must_use]
    pub const fn rule(&self) -> MachineAssuranceRule {
        self.rule
    }

    /// Stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }

    /// Primary subject that failed the rule.
    #[must_use]
    pub const fn subject(&self) -> &MachineAssuranceSubject {
        &self.subject
    }

    /// Optional related subject providing deterministic context.
    #[must_use]
    pub const fn related(&self) -> Option<&MachineAssuranceSubject> {
        self.related.as_ref()
    }
}

/// Complete deterministic refusal from assurance admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssuranceRefusal {
    findings: Vec<MachineAssuranceFinding>,
    identity_error: Option<CanonicalError>,
}

impl MachineAssuranceRefusal {
    /// Sorted, duplicate-free complete finding set.
    #[must_use]
    pub fn findings(&self) -> &[MachineAssuranceFinding] {
        &self.findings
    }

    /// Canonical encoder error when identity publication itself failed.
    #[must_use]
    pub const fn identity_error(&self) -> Option<&CanonicalError> {
        self.identity_error.as_ref()
    }

    /// Stable top-level refusal code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        "MachineAssuranceRefused"
    }
}

impl fmt::Display for MachineAssuranceRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "machine assurance refused with {} finding(s)",
            self.findings.len()
        )
    }
}

impl std::error::Error for MachineAssuranceRefusal {}

/// Collection sizes retained for every assurance-admission attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineAssuranceSubmittedCounts {
    /// Submitted sensor count.
    pub sensors: usize,
    /// Submitted local experiment-binding count.
    pub experiments: usize,
    /// Submitted Context-of-Use binding count.
    pub contexts: usize,
    /// Submitted hazard count.
    pub hazards: usize,
    /// Submitted fault count.
    pub faults: usize,
    /// Submitted accounting-window count.
    pub accounting_windows: usize,
    /// Submitted fidelity-rung count.
    pub fidelity_rungs: usize,
    /// Supplied admitted V&V case count.
    pub vv_cases: usize,
    /// Bounded nested count across supplied V&V cases.
    pub vv_nested_references: usize,
    /// Bounded aggregate nested count across the draft and V&V cases.
    pub nested_references: usize,
}

/// Mutable-by-construction PR-4 overlay with no authority before admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssuranceDraft {
    /// Untrusted sensor declarations.
    pub sensors: Vec<SensorSpec>,
    /// Untrusted local experiment bindings.
    pub experiments: Vec<ExperimentSpec>,
    /// Untrusted Context-of-Use and QoI bindings.
    pub contexts: Vec<ContextBinding>,
    /// Untrusted hazard declarations.
    pub hazards: Vec<HazardSpec>,
    /// Untrusted fault declarations.
    pub faults: Vec<FaultSpec>,
    /// Untrusted signed accounting windows.
    pub accounting_windows: Vec<AccountingWindow>,
    /// Untrusted complete fidelity policy.
    pub fidelity: FidelityPolicy,
}

impl MachineAssuranceDraft {
    /// Admit this overlay against one exact graph/behavior pair.
    ///
    /// # Errors
    /// Returns every deterministic finding discovered within public bounds.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
        behavior: &AdmittedMachineBehavior,
        vv_cases: &[AdmittedVvCase],
    ) -> Result<AdmittedMachineAssurance, MachineAssuranceRefusal> {
        self.admit_with_decision(graph, behavior, vv_cases)
            .into_result()
    }

    /// Attempt admission while retaining submitted collection counts.
    #[must_use]
    pub fn admit_with_decision(
        self,
        graph: &AdmittedMachineGraph,
        behavior: &AdmittedMachineBehavior,
        vv_cases: &[AdmittedVvCase],
    ) -> MachineAssuranceAdmissionDecision {
        let submitted = submitted_counts(&self, vv_cases);
        MachineAssuranceAdmissionDecision {
            submitted,
            result: admit_machine_assurance(self, graph, behavior, vv_cases),
        }
    }
}

/// Canonically ordered admitted PR-4 overlay plus its semantic receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineAssurance {
    base_graph: MachineGraphIdV1,
    base_behavior: MachineBehaviorIdV1,
    sensors: Vec<SensorSpec>,
    experiments: Vec<ExperimentSpec>,
    contexts: Vec<ContextBinding>,
    vv_cases: Vec<VvCaseBinding>,
    hazards: Vec<HazardSpec>,
    faults: Vec<FaultSpec>,
    accounting_windows: Vec<AccountingWindow>,
    fidelity: FidelityPolicy,
    receipt: IdentityReceipt<MachineAssuranceIdV1>,
}

impl AdmittedMachineAssurance {
    /// Exact base graph identity.
    #[must_use]
    pub const fn base_graph(&self) -> MachineGraphIdV1 {
        self.base_graph
    }

    /// Exact base behavior identity.
    #[must_use]
    pub const fn base_behavior(&self) -> MachineBehaviorIdV1 {
        self.base_behavior
    }

    /// Domain-separated semantic identity of the admitted overlay.
    #[must_use]
    pub const fn identity(&self) -> MachineAssuranceIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineAssuranceIdV1> {
        self.receipt
    }

    /// Canonically ordered sensors.
    #[must_use]
    pub fn sensors(&self) -> &[SensorSpec] {
        &self.sensors
    }

    /// Canonically ordered experiment bindings.
    #[must_use]
    pub fn experiments(&self) -> &[ExperimentSpec] {
        &self.experiments
    }

    /// Canonically ordered Context-of-Use bindings.
    #[must_use]
    pub fn contexts(&self) -> &[ContextBinding] {
        &self.contexts
    }

    /// Exact evidence-case admission receipts transitively bound by identity.
    #[must_use]
    pub fn vv_cases(&self) -> &[VvCaseBinding] {
        &self.vv_cases
    }

    /// Canonically ordered hazards.
    #[must_use]
    pub fn hazards(&self) -> &[HazardSpec] {
        &self.hazards
    }

    /// Canonically ordered faults.
    #[must_use]
    pub fn faults(&self) -> &[FaultSpec] {
        &self.faults
    }

    /// Canonically ordered accounting windows.
    #[must_use]
    pub fn accounting_windows(&self) -> &[AccountingWindow] {
        &self.accounting_windows
    }

    /// Complete canonical fidelity policy.
    #[must_use]
    pub const fn fidelity(&self) -> &FidelityPolicy {
        &self.fidelity
    }
}

/// Bounded deterministic outcome summary for one assurance admission.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssuranceAdmissionDecision {
    submitted: MachineAssuranceSubmittedCounts,
    result: Result<AdmittedMachineAssurance, MachineAssuranceRefusal>,
}

impl MachineAssuranceAdmissionDecision {
    /// Submitted collection counts retained even on refusal.
    #[must_use]
    pub const fn submitted_counts(&self) -> MachineAssuranceSubmittedCounts {
        self.submitted
    }

    /// Stable admitted/refused outcome code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "MachineAssuranceAdmitted",
            Err(_) => "MachineAssuranceRefused",
        }
    }

    /// Borrow the admitted overlay or complete refusal.
    #[must_use]
    pub fn result(&self) -> Result<&AdmittedMachineAssurance, &MachineAssuranceRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision into its admitted overlay or refusal.
    #[must_use]
    pub fn into_result(self) -> Result<AdmittedMachineAssurance, MachineAssuranceRefusal> {
        self.result
    }
}

fn submitted_counts(
    draft: &MachineAssuranceDraft,
    vv_cases: &[AdmittedVvCase],
) -> MachineAssuranceSubmittedCounts {
    let draft_nested_references = draft
        .experiments
        .iter()
        .fold(0usize, |count, experiment| {
            count
                .saturating_add(experiment.instruments.len())
                .saturating_add(experiment.qois.len())
        })
        .saturating_add(draft.contexts.iter().fold(0usize, |count, context| {
            count.saturating_add(context.qois.len()).saturating_add(
                context
                    .qois
                    .iter()
                    .fold(0usize, |count, qoi| count.saturating_add(qoi.inputs.len())),
            )
        }))
        .saturating_add(draft.hazards.iter().fold(0usize, |count, hazard| {
            count
                .saturating_add(hazard.scope.len())
                .saturating_add(hazard.assumptions.len())
        }))
        .saturating_add(draft.faults.iter().fold(0usize, |count, fault| {
            count
                .saturating_add(fault.affected.len())
                .saturating_add(fault.hazards.len())
        }))
        .saturating_add(
            draft
                .accounting_windows
                .iter()
                .fold(0usize, |count, window| {
                    count.saturating_add(window.entries.len())
                }),
        )
        .saturating_add(draft.fidelity.baselines.len())
        .saturating_add(draft.fidelity.escalations.len())
        .saturating_add(draft.fidelity.rungs.iter().fold(0usize, |count, rung| {
            count
                .saturating_add(rung.falsifiers.len())
                .saturating_add(rung.qois.len())
        }));
    let vv_nested_references = vv_nested_reference_count(vv_cases);
    MachineAssuranceSubmittedCounts {
        sensors: draft.sensors.len(),
        experiments: draft.experiments.len(),
        contexts: draft.contexts.len(),
        hazards: draft.hazards.len(),
        faults: draft.faults.len(),
        accounting_windows: draft.accounting_windows.len(),
        fidelity_rungs: draft.fidelity.rungs.len(),
        vv_cases: vv_cases.len(),
        vv_nested_references,
        nested_references: draft_nested_references.saturating_add(vv_nested_references),
    }
}

const MACHINE_ASSURANCE_NESTED_OVERFLOW_COUNT: usize = MAX_MACHINE_ASSURANCE_NESTED_REFERENCES + 1;

fn add_nested_count(total: &mut usize, amount: usize) -> bool {
    *total = (*total)
        .saturating_add(amount)
        .min(MACHINE_ASSURANCE_NESTED_OVERFLOW_COUNT);
    *total == MACHINE_ASSURANCE_NESTED_OVERFLOW_COUNT
}

fn header_nested_reference_count(header: &ArtifactHeader) -> usize {
    header
        .units()
        .len()
        .saturating_add(header.versions().len())
        .saturating_add(header.capabilities().len())
}

fn vv_nested_reference_count(vv_cases: &[AdmittedVvCase]) -> usize {
    let mut total = 0usize;
    for admitted in vv_cases {
        if add_nested_count(&mut total, vv_case_nested_reference_count(admitted)) {
            break;
        }
    }
    total
}

#[allow(clippy::too_many_lines)]
fn vv_case_nested_reference_count(admitted: &AdmittedVvCase) -> usize {
    let mut total = 0usize;
    macro_rules! add {
        ($amount:expr) => {
            if add_nested_count(&mut total, $amount) {
                return total;
            }
        };
    }

    let case = admitted.case();
    let context = case.context();
    add!(header_nested_reference_count(context.header()));
    add!(context.qois().len());
    add!(context.applicability().numeric().len());
    add!(context.applicability().categorical().len());
    for axis in context.applicability().categorical().values() {
        add!(axis.allowed().len());
    }

    let plan = case.validation_plan();
    add!(header_nested_reference_count(plan.header()));
    add!(plan.by_qoi().len());
    for row in plan.by_qoi().values() {
        add!(row.experiments().len());
        add!(row.metrics().len());
    }

    add!(case.experiments().len());
    for experiment in case.experiments().values() {
        add!(header_nested_reference_count(experiment.header()));
        add!(experiment.qois().len());
        add!(experiment.observation_ids().len());
        add!(experiment.manifest().rows().len());
        add!(experiment.instruments().len());
        if let fs_evidence::vv::ClockSynchronization::Synchronized { clock_ids, .. } =
            experiment.clocks()
        {
            add!(clock_ids.len());
        }
        add!(
            experiment
                .repeatability()
                .covariance()
                .lower_triangle()
                .len()
        );
    }

    add!(case.splits().len());
    for split in case.splits().values() {
        add!(header_nested_reference_count(split.header()));
        add!(split.calibration_ids().len());
        add!(split.validation_ids().len());
        add!(split.blind_holdout_len());
        add!(split.blind_sources().len());
    }

    add!(case.solution_verification().len());
    for solution in case.solution_verification().values() {
        add!(header_nested_reference_count(solution.header()));
    }

    add!(case.predictions().len());
    for prediction in case.predictions().values() {
        add!(header_nested_reference_count(prediction.header()));
        add!(prediction.dependencies().len());
        for dependency in prediction.dependencies() {
            if let Some(observations) = dependency.observations() {
                add!(observations.ids().len());
            }
        }
        add!(prediction.waterfall().terms().len());
        if let fs_evidence::vv::WaterfallMode::Probabilistic { dependence, .. } =
            prediction.waterfall().mode()
        {
            add!(dependence.values().len());
        }
        add!(prediction.validation_metrics().len());
        for metric in prediction.validation_metrics() {
            add!(metric.observations().ids().len());
        }
        add!(prediction.posterior_checks().len());
        for check in prediction.posterior_checks() {
            add!(check.observations().ids().len());
        }
        add!(prediction.applicability_point().numeric().len());
        add!(prediction.applicability_point().categorical().len());
        match prediction.applicability() {
            fs_evidence::vv::ApplicabilityDecision::InDomain => {}
            fs_evidence::vv::ApplicabilityDecision::Demoted { violations }
            | fs_evidence::vv::ApplicabilityDecision::Refused { violations } => {
                add!(violations.len());
            }
        }
        add!(prediction.evidence_axes().axes().len());
        for status in prediction.evidence_axes().axes().values() {
            if let EvidenceAxisStatus::Present { artifacts } = status {
                add!(artifacts.len());
            }
        }
        add!(prediction.assumption_checks().len());
    }

    add!(header_nested_reference_count(case.assumptions().header()));
    add!(case.assumptions().rows().len());
    add!(admitted.receipt().qois().len());
    add!(admitted.receipt().artifact_hashes().len());
    total
}

fn resource_limit_findings(
    counts: MachineAssuranceSubmittedCounts,
) -> Vec<MachineAssuranceFinding> {
    let over_limit = counts.sensors > MAX_MACHINE_ASSURANCE_SENSORS
        || counts.experiments > MAX_MACHINE_ASSURANCE_EXPERIMENTS
        || counts.contexts > MAX_MACHINE_ASSURANCE_CONTEXTS
        || counts.vv_cases > MAX_MACHINE_ASSURANCE_CONTEXTS
        || counts.hazards > MAX_MACHINE_ASSURANCE_HAZARDS
        || counts.faults > MAX_MACHINE_ASSURANCE_FAULTS
        || counts.accounting_windows > MAX_MACHINE_ASSURANCE_ACCOUNTING_WINDOWS
        || counts.fidelity_rungs > MAX_MACHINE_ASSURANCE_FIDELITY_RUNGS
        || counts.nested_references > MAX_MACHINE_ASSURANCE_NESTED_REFERENCES;
    if over_limit {
        vec![MachineAssuranceFinding::new(
            MachineAssuranceRule::ResourceLimit,
            MachineAssuranceSubject::Assurance,
            None,
        )]
    } else {
        Vec::new()
    }
}

fn assurance_refusal(
    mut findings: Vec<MachineAssuranceFinding>,
    identity_error: Option<CanonicalError>,
) -> MachineAssuranceRefusal {
    findings.sort();
    findings.dedup();
    MachineAssuranceRefusal {
        findings,
        identity_error,
    }
}

fn canonicalize_assurance_draft(draft: &mut MachineAssuranceDraft) {
    draft.sensors.sort_by_key(sensor_row);
    for experiment in &mut draft.experiments {
        experiment.instruments.sort();
        experiment.qois.sort();
    }
    draft.experiments.sort_by_key(experiment_row);
    for context in &mut draft.contexts {
        for qoi in &mut context.qois {
            qoi.inputs.sort_by_key(qoi_input_row);
        }
        context.qois.sort_by_key(qoi_binding_row);
    }
    draft.contexts.sort_by_key(context_row);
    for hazard in &mut draft.hazards {
        hazard.scope.sort();
        hazard.assumptions.sort();
    }
    draft.hazards.sort_by_key(hazard_row);
    for fault in &mut draft.faults {
        fault.affected.sort();
        fault.hazards.sort();
    }
    draft.faults.sort_by_key(fault_row);
    for window in &mut draft.accounting_windows {
        window.entries.sort_by_key(accounting_entry_row);
    }
    draft.accounting_windows.sort_by_key(accounting_window_row);
    draft.fidelity.baselines.sort();
    for rung in &mut draft.fidelity.rungs {
        rung.falsifiers.sort();
        rung.qois.sort();
    }
    draft.fidelity.rungs.sort_by_key(fidelity_rung_row);
    draft.fidelity.escalations.sort_by_key(escalation_row);
}

struct AssuranceIndex<'a> {
    clocks: BTreeSet<ClockId>,
    subsystems: BTreeMap<SubsystemId, &'a super::SubsystemSpec>,
    terminals: BTreeMap<TerminalId, &'a TerminalSpec>,
    states: BTreeMap<StateSlotId, &'a StateSlotContract>,
    ports: BTreeSet<PortId>,
    relations: BTreeMap<RelationId, &'a RelationSpec>,
    interfaces: BTreeMap<InterfaceId, &'a InterfaceBinding>,
    elements: BTreeSet<MachineElementId>,
}

impl<'a> AssuranceIndex<'a> {
    fn new(graph: &'a AdmittedMachineGraph, behavior: &'a AdmittedMachineBehavior) -> Self {
        let clocks = graph
            .clocks()
            .iter()
            .map(|clock| clock.id.clone())
            .collect();
        let subsystems = graph
            .subsystems()
            .iter()
            .map(|subsystem| (subsystem.id.clone(), subsystem))
            .collect();
        let terminals = graph
            .terminals()
            .iter()
            .map(|terminal| (terminal.id.clone(), terminal))
            .collect();
        let states = behavior
            .state_contracts()
            .iter()
            .map(|state| (state.id.clone(), state))
            .collect();
        let ports = graph.ports().iter().map(|port| port.id.clone()).collect();
        let relations = graph
            .relations()
            .iter()
            .map(|relation| (relation.id.clone(), relation))
            .collect();
        let interfaces = graph
            .interfaces()
            .iter()
            .map(|interface| (interface.id.clone(), interface))
            .collect();
        let mut elements = BTreeSet::new();
        for subsystem in graph.subsystems() {
            elements.extend(subsystem.bodies.iter().cloned().map(MachineElementId::from));
            elements.extend(
                subsystem
                    .surface_patches
                    .iter()
                    .cloned()
                    .map(MachineElementId::from),
            );
            elements.extend(
                subsystem
                    .contact_features
                    .iter()
                    .cloned()
                    .map(MachineElementId::from),
            );
            elements.extend(
                subsystem
                    .state_slots
                    .iter()
                    .cloned()
                    .map(MachineElementId::from),
            );
        }
        elements.extend(
            graph
                .terminals()
                .iter()
                .map(|terminal| MachineElementId::from(terminal.id.clone())),
        );
        elements.extend(
            graph
                .ports()
                .iter()
                .map(|port| MachineElementId::from(port.id.clone())),
        );
        Self {
            clocks,
            subsystems,
            terminals,
            states,
            ports,
            relations,
            interfaces,
            elements,
        }
    }

    fn observation_contract(
        &self,
        target: &ObservationTarget,
    ) -> Option<(TerminalQuantitySpec, TerminalShape, &ClockId, &FrameBinding)> {
        match target {
            ObservationTarget::Terminal(id) => self.terminals.get(id).map(|terminal| {
                (
                    terminal.quantity,
                    terminal.shape,
                    &terminal.clock,
                    &terminal.frame,
                )
            }),
            ObservationTarget::State(id) => self
                .states
                .get(id)
                .map(|state| (state.quantity, state.shape, &state.clock, &state.frame)),
        }
    }

    fn scope_exists(&self, scope: &MachineScope) -> bool {
        match scope {
            MachineScope::WholeMachine => true,
            MachineScope::Subsystem(id) => self.subsystems.contains_key(id),
            MachineScope::Element(id) => self.elements.contains(id),
            MachineScope::Relation(id) => self.relations.contains_key(id),
            MachineScope::Interface(id) => self.interfaces.contains_key(id),
        }
    }

    fn accounting_target_exists(&self, target: &AccountingTarget) -> bool {
        match target {
            AccountingTarget::Relation(id) => self.relations.contains_key(id),
            AccountingTarget::Port(id) => self.ports.contains(id),
            AccountingTarget::Interface(id) => self.interfaces.contains_key(id),
            AccountingTarget::State(id) => self.states.contains_key(id),
        }
    }

    fn accounting_atoms(
        &self,
        target: &AccountingTarget,
        role: AccountingRole,
    ) -> BTreeSet<AccountingAtom> {
        let mut atoms = BTreeSet::new();
        match target {
            AccountingTarget::Relation(id) => {
                atoms.insert(AccountingAtom::Relation(id.clone()));
                if role == AccountingRole::Storage {
                    if let Some(relation) = self.relations.get(id) {
                        if let RelationMode::Stateful { state_slot } = &relation.mode {
                            atoms.insert(AccountingAtom::State(state_slot.clone()));
                        }
                    }
                }
            }
            AccountingTarget::Port(id) => {
                atoms.insert(AccountingAtom::Port(id.clone()));
            }
            AccountingTarget::Interface(id) => {
                if let Some(interface) = self.interfaces.get(id) {
                    atoms.insert(AccountingAtom::Port(interface.negative.clone()));
                    atoms.insert(AccountingAtom::Port(interface.positive.clone()));
                }
            }
            AccountingTarget::State(id) => {
                atoms.insert(AccountingAtom::State(id.clone()));
            }
        }
        atoms
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum AccountingAtom {
    Relation(RelationId),
    Port(PortId),
    State(StateSlotId),
}

#[derive(Debug, Clone)]
struct VvExperimentView {
    reference: ArtifactRef,
    qois: BTreeSet<QoiId>,
    instruments: BTreeSet<ArtifactId>,
}

#[derive(Debug, Clone)]
struct VvCaseView {
    context: ArtifactRef,
    validation_plan: ArtifactRef,
    qoi_units: BTreeMap<QoiId, UnitId>,
    experiments: BTreeMap<ArtifactId, VvExperimentView>,
    assumptions: BTreeSet<AssumptionId>,
    binding: VvCaseBinding,
}

fn artifact_ref_from_case(
    case: &AdmittedVvCase,
    kind: ArtifactKind,
    id: &ArtifactId,
) -> Option<ArtifactRef> {
    case.receipt()
        .artifact_hashes()
        .get(&(kind, id.clone()))
        .copied()
        .map(|hash| ArtifactRef::new(kind, id.clone(), hash))
}

fn index_vv_cases(
    cases: &[AdmittedVvCase],
    findings: &mut Vec<MachineAssuranceFinding>,
) -> BTreeMap<ArtifactId, VvCaseView> {
    let mut indexed = BTreeMap::new();
    for admitted in cases {
        let case = admitted.case();
        let context_id = case.context().id().clone();
        let subject = || MachineAssuranceSubject::Artifact(context_id.clone());
        if !admitted.receipt().has_valid_binding() || admitted.receipt().verify_case(case).is_err()
        {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidVvCaseReceipt,
                subject(),
                None,
            ));
            continue;
        }
        let Some(context) =
            artifact_ref_from_case(admitted, ArtifactKind::ContextOfUse, &context_id)
        else {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidVvCaseReceipt,
                subject(),
                None,
            ));
            continue;
        };
        let Some(validation_plan) = artifact_ref_from_case(
            admitted,
            ArtifactKind::ValidationPlan,
            case.validation_plan().id(),
        ) else {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidVvCaseReceipt,
                subject(),
                None,
            ));
            continue;
        };
        let mut experiments = BTreeMap::new();
        let mut complete = true;
        for experiment in case.experiments().values() {
            let Some(reference) =
                artifact_ref_from_case(admitted, ArtifactKind::ExperimentArtifact, experiment.id())
            else {
                complete = false;
                break;
            };
            experiments.insert(
                experiment.id().clone(),
                VvExperimentView {
                    reference,
                    qois: experiment.qois().clone(),
                    instruments: experiment
                        .instruments()
                        .iter()
                        .map(|instrument| instrument.instrument_id().clone())
                        .collect(),
                },
            );
        }
        if !complete {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidVvCaseReceipt,
                subject(),
                None,
            ));
            continue;
        }
        let qoi_units = case
            .context()
            .qois()
            .iter()
            .map(|(id, qoi)| (id.clone(), qoi.unit().clone()))
            .collect();
        let binding = VvCaseBinding {
            context: context.clone(),
            schema_version: admitted.receipt().schema_version(),
            ruleset_version: admitted.receipt().ruleset_version(),
            case_hash: admitted.receipt().case_hash(),
            receipt_hash: admitted.receipt().receipt_hash(),
        };
        let view = VvCaseView {
            context,
            validation_plan,
            qoi_units,
            experiments,
            assumptions: case.assumptions().rows().keys().cloned().collect(),
            binding,
        };
        if indexed.insert(context_id.clone(), view).is_some() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidVvCaseReceipt,
                subject(),
                None,
            ));
        }
    }
    indexed
}

#[allow(clippy::too_many_lines)]
fn admit_machine_assurance(
    mut draft: MachineAssuranceDraft,
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    vv_cases: &[AdmittedVvCase],
) -> Result<AdmittedMachineAssurance, MachineAssuranceRefusal> {
    let counts = submitted_counts(&draft, vv_cases);
    let resource_findings = resource_limit_findings(counts);
    if !resource_findings.is_empty() {
        return Err(assurance_refusal(resource_findings, None));
    }

    canonicalize_assurance_draft(&mut draft);
    let mut findings = Vec::new();
    if behavior.base_graph() != graph.identity() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::BehaviorGraphMismatch,
            MachineAssuranceSubject::Graph(graph.identity()),
            Some(MachineAssuranceSubject::Graph(behavior.base_graph())),
        ));
    }
    let index = AssuranceIndex::new(graph, behavior);
    let vv_index = index_vv_cases(vv_cases, &mut findings);
    let sensors = validate_sensors(&draft.sensors, &index, &mut findings);
    let (contexts, qois) =
        validate_contexts(&draft.contexts, &sensors, &index, &vv_index, &mut findings);
    validate_experiments(
        &draft.experiments,
        &sensors,
        &contexts,
        &qois,
        &vv_index,
        &mut findings,
    );
    validate_hazards_and_faults(
        &draft.hazards,
        &draft.faults,
        &contexts,
        &index,
        &vv_index,
        &mut findings,
    );
    validate_accounting(&draft.accounting_windows, &contexts, &index, &mut findings);
    validate_fidelity(&draft.fidelity, &qois, &index, &mut findings);

    let declared_context_ids: BTreeSet<ArtifactId> = contexts.keys().cloned().collect();
    for context_id in vv_index.keys() {
        if !declared_context_ids.contains(context_id) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnexpectedVvCase,
                MachineAssuranceSubject::Artifact(context_id.clone()),
                None,
            ));
        }
    }

    if !findings.is_empty() {
        return Err(assurance_refusal(findings, None));
    }

    let vv_bindings: Vec<VvCaseBinding> = draft
        .contexts
        .iter()
        .map(|context| {
            vv_index
                .get(context.context.id())
                .expect("context closure checked before identity")
                .binding
                .clone()
        })
        .collect();
    let receipt = match machine_assurance_identity(
        &draft,
        graph.identity(),
        behavior.identity(),
        &vv_bindings,
    ) {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(assurance_refusal(
                vec![MachineAssuranceFinding::new(
                    MachineAssuranceRule::Identity,
                    MachineAssuranceSubject::Assurance,
                    None,
                )],
                Some(error),
            ));
        }
    };
    Ok(AdmittedMachineAssurance {
        base_graph: graph.identity(),
        base_behavior: behavior.identity(),
        sensors: draft.sensors,
        experiments: draft.experiments,
        contexts: draft.contexts,
        vv_cases: vv_bindings,
        hazards: draft.hazards,
        faults: draft.faults,
        accounting_windows: draft.accounting_windows,
        fidelity: draft.fidelity,
        receipt,
    })
}

fn validate_sensors<'a>(
    sensors: &'a [SensorSpec],
    index: &AssuranceIndex<'_>,
    findings: &mut Vec<MachineAssuranceFinding>,
) -> BTreeMap<SensorId, &'a SensorSpec> {
    if sensors.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptySensors,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut by_id = BTreeMap::new();
    let mut outputs = BTreeMap::<TerminalId, SensorId>::new();
    for sensor in sensors {
        let subject = || MachineAssuranceSubject::Sensor(sensor.id.clone());
        if let Some(first) = by_id.insert(sensor.id.clone(), sensor) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateSensor,
                subject(),
                Some(MachineAssuranceSubject::Sensor(first.id.clone())),
            ));
        }
        if !index.subsystems.contains_key(&sensor.owner) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownSensorOwner,
                subject(),
                Some(MachineAssuranceSubject::Subsystem(sensor.owner.clone())),
            ));
        }
        if !sensor.quantity.is_admitted() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnsupportedSensorQuantity,
                subject(),
                None,
            ));
        }
        if !index.clocks.contains(&sensor.clock) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownSensorClock,
                subject(),
                Some(MachineAssuranceSubject::Clock(sensor.clock.clone())),
            ));
        }
        match index.observation_contract(&sensor.target) {
            None => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownSensorTarget,
                subject(),
                None,
            )),
            Some((quantity, shape, clock, frame)) => {
                compare_sensor_contract(sensor, quantity, shape, clock, frame, findings);
            }
        }
        if let SensorExposure::PlantSignal { output } = &sensor.exposure {
            if let Some(first) = outputs.insert(output.clone(), sensor.id.clone()) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateSensorOutput,
                    subject(),
                    Some(MachineAssuranceSubject::Sensor(first)),
                ));
            }
            match index.terminals.get(output) {
                None => findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownSensorOutput,
                    subject(),
                    Some(MachineAssuranceSubject::Element(output.clone().into())),
                )),
                Some(terminal) => {
                    if terminal.owner != sensor.owner {
                        findings.push(MachineAssuranceFinding::new(
                            MachineAssuranceRule::SensorOutputOwnerMismatch,
                            subject(),
                            Some(MachineAssuranceSubject::Subsystem(terminal.owner.clone())),
                        ));
                    }
                    if terminal.causality != TerminalCausality::Output {
                        findings.push(MachineAssuranceFinding::new(
                            MachineAssuranceRule::SensorOutputCausalityGap,
                            subject(),
                            Some(MachineAssuranceSubject::Element(output.clone().into())),
                        ));
                    }
                    if terminal.quantity != sensor.quantity
                        || terminal.shape != sensor.shape
                        || terminal.clock != sensor.clock
                        || terminal.frame != sensor.frame
                    {
                        findings.push(MachineAssuranceFinding::new(
                            MachineAssuranceRule::SensorOutputContractGap,
                            subject(),
                            Some(MachineAssuranceSubject::Element(output.clone().into())),
                        ));
                    }
                }
            }
        }
    }
    by_id
}

fn compare_sensor_contract(
    sensor: &SensorSpec,
    quantity: TerminalQuantitySpec,
    shape: TerminalShape,
    clock: &ClockId,
    frame: &FrameBinding,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    let subject = || MachineAssuranceSubject::Sensor(sensor.id.clone());
    if sensor.quantity != quantity {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::SensorQuantityGap,
            subject(),
            None,
        ));
    }
    if sensor.shape != shape {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::SensorShapeGap,
            subject(),
            None,
        ));
    }
    if matches!(sensor.timing, ObservationTiming::Direct) && &sensor.clock != clock {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::SensorClockGap,
            subject(),
            Some(MachineAssuranceSubject::Clock(clock.clone())),
        ));
    }
    if &sensor.frame != frame {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::SensorFrameGap,
            subject(),
            None,
        ));
    }
}

fn validate_contexts<'a>(
    contexts: &'a [ContextBinding],
    sensors: &BTreeMap<SensorId, &SensorSpec>,
    index: &AssuranceIndex<'_>,
    vv_index: &BTreeMap<ArtifactId, VvCaseView>,
    findings: &mut Vec<MachineAssuranceFinding>,
) -> (
    BTreeMap<ArtifactId, ArtifactRef>,
    BTreeMap<ContextQoiKey, &'a QoiBinding>,
) {
    if contexts.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptyContexts,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut context_refs = BTreeMap::new();
    let mut qois = BTreeMap::new();
    for context in contexts {
        let context_subject = || MachineAssuranceSubject::Context(context.context.clone());
        if context.context.kind() != ArtifactKind::ContextOfUse {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidContextArtifact,
                context_subject(),
                None,
            ));
        }
        if context.validation_plan.kind() != ArtifactKind::ValidationPlan {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidValidationPlanArtifact,
                context_subject(),
                Some(MachineAssuranceSubject::Context(
                    context.validation_plan.clone(),
                )),
            ));
        }
        if let Some(first) =
            context_refs.insert(context.context.id().clone(), context.context.clone())
        {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateContext,
                context_subject(),
                Some(MachineAssuranceSubject::Context(first)),
            ));
        }
        let evidence = vv_index.get(context.context.id());
        match evidence {
            None => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::MissingVvCase,
                context_subject(),
                None,
            )),
            Some(evidence) => {
                if evidence.context != context.context
                    || evidence.validation_plan != context.validation_plan
                {
                    findings.push(MachineAssuranceFinding::new(
                        MachineAssuranceRule::InvalidArtifactReference,
                        context_subject(),
                        None,
                    ));
                }
            }
        }
        if context.qois.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyContextQois,
                context_subject(),
                None,
            ));
        }
        let submitted_qois: BTreeSet<QoiId> =
            context.qois.iter().map(|qoi| qoi.id.clone()).collect();
        if let Some(evidence) = evidence {
            let expected: BTreeSet<QoiId> = evidence.qoi_units.keys().cloned().collect();
            if submitted_qois != expected {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::ContextQoiSetMismatch,
                    context_subject(),
                    None,
                ));
            }
        }
        for qoi in &context.qois {
            let key = ContextQoiKey {
                context: context.context.id().clone(),
                qoi: qoi.id.clone(),
            };
            let qoi_subject = || MachineAssuranceSubject::Qoi(key.clone());
            if qois.insert(key.clone(), qoi).is_some() {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateQoi,
                    qoi_subject(),
                    None,
                ));
            }
            if let Some(expected) = evidence.and_then(|case| case.qoi_units.get(&qoi.id)) {
                if expected != &qoi.unit {
                    findings.push(MachineAssuranceFinding::new(
                        MachineAssuranceRule::QoiUnitMismatch,
                        qoi_subject(),
                        None,
                    ));
                }
            }
            if qoi.inputs.is_empty() {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::EmptyQoiInputs,
                    qoi_subject(),
                    None,
                ));
            }
            for pair in qoi.inputs.windows(2) {
                if qoi_input_row(&pair[0]) == qoi_input_row(&pair[1]) {
                    findings.push(MachineAssuranceFinding::new(
                        MachineAssuranceRule::DuplicateQoiInput,
                        qoi_subject(),
                        None,
                    ));
                }
            }
            for input in &qoi.inputs {
                validate_qoi_input(input, &key, sensors, index, findings);
            }
        }
    }
    (context_refs, qois)
}

fn validate_qoi_input(
    input: &QoiInput,
    key: &ContextQoiKey,
    sensors: &BTreeMap<SensorId, &SensorSpec>,
    index: &AssuranceIndex<'_>,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    let subject = || MachineAssuranceSubject::Qoi(key.clone());
    if !input.quantity.is_admitted() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::UnsupportedQoiQuantity,
            subject(),
            None,
        ));
    }
    let contract = match &input.target {
        QoiTarget::Sensor(id) => sensors
            .get(id)
            .map(|sensor| (sensor.quantity, sensor.shape)),
        QoiTarget::Terminal(id) => index
            .terminals
            .get(id)
            .map(|terminal| (terminal.quantity, terminal.shape)),
        QoiTarget::State(id) => index
            .states
            .get(id)
            .map(|state| (state.quantity, state.shape)),
    };
    let Some((quantity, shape)) = contract else {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::UnknownQoiTarget,
            subject(),
            None,
        ));
        return;
    };
    if input.quantity != quantity {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::QoiQuantityGap,
            subject(),
            None,
        ));
    }
    if input.shape != shape {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::QoiShapeGap,
            subject(),
            None,
        ));
    }
}

fn validate_experiments(
    experiments: &[ExperimentSpec],
    sensors: &BTreeMap<SensorId, &SensorSpec>,
    contexts: &BTreeMap<ArtifactId, ArtifactRef>,
    qois: &BTreeMap<ContextQoiKey, &QoiBinding>,
    vv_index: &BTreeMap<ArtifactId, VvCaseView>,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    if experiments.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptyExperiments,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut by_id = BTreeMap::<ExperimentId, &ExperimentSpec>::new();
    let mut qoi_coverage = BTreeSet::<ContextQoiKey>::new();
    let mut artifact_coverage = BTreeSet::<(ArtifactId, ArtifactId)>::new();
    for experiment in experiments {
        let subject = || MachineAssuranceSubject::Experiment(experiment.id.clone());
        if let Some(first) = by_id.insert(experiment.id.clone(), experiment) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateExperiment,
                subject(),
                Some(MachineAssuranceSubject::Experiment(first.id.clone())),
            ));
        }
        if experiment.artifact.kind() != ArtifactKind::ExperimentArtifact {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidExperimentArtifact,
                subject(),
                Some(MachineAssuranceSubject::Context(
                    experiment.artifact.clone(),
                )),
            ));
        }
        match contexts.get(experiment.context.id()) {
            Some(context) if context == &experiment.context => {}
            _ => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownExperimentContext,
                subject(),
                Some(MachineAssuranceSubject::Context(experiment.context.clone())),
            )),
        }
        let evidence = vv_index
            .get(experiment.context.id())
            .and_then(|case| case.experiments.get(experiment.artifact.id()));
        match evidence {
            Some(evidence) if evidence.reference == experiment.artifact => {
                if !artifact_coverage.insert((
                    experiment.context.id().clone(),
                    experiment.artifact.id().clone(),
                )) {
                    findings.push(MachineAssuranceFinding::new(
                        MachineAssuranceRule::DuplicateExperimentArtifactBinding,
                        subject(),
                        Some(MachineAssuranceSubject::Artifact(
                            experiment.artifact.id().clone(),
                        )),
                    ));
                }
            }
            _ => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::InvalidArtifactReference,
                subject(),
                Some(MachineAssuranceSubject::Artifact(
                    experiment.artifact.id().clone(),
                )),
            )),
        }
        if experiment.instruments.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyExperimentSensors,
                subject(),
                None,
            ));
        }
        let mut seen_sensors = BTreeSet::new();
        let mut seen_instruments = BTreeSet::new();
        for binding in &experiment.instruments {
            if !seen_sensors.insert(binding.sensor.clone()) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateExperimentSensor,
                    subject(),
                    Some(MachineAssuranceSubject::Sensor(binding.sensor.clone())),
                ));
            }
            if !seen_instruments.insert(binding.instrument.clone()) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateExperimentInstrument,
                    subject(),
                    Some(MachineAssuranceSubject::Artifact(
                        binding.instrument.clone(),
                    )),
                ));
            }
            if !sensors.contains_key(&binding.sensor) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownExperimentSensor,
                    subject(),
                    Some(MachineAssuranceSubject::Sensor(binding.sensor.clone())),
                ));
            }
            if evidence.is_some_and(|view| !view.instruments.contains(&binding.instrument)) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownExperimentInstrument,
                    subject(),
                    Some(MachineAssuranceSubject::Artifact(
                        binding.instrument.clone(),
                    )),
                ));
            }
        }
        if let Some(evidence) = evidence {
            if seen_instruments != evidence.instruments {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::ExperimentInstrumentSetMismatch,
                    subject(),
                    None,
                ));
            }
        }
        if experiment.qois.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyExperimentQois,
                subject(),
                None,
            ));
        }
        for pair in experiment.qois.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateExperimentQoi,
                    subject(),
                    Some(MachineAssuranceSubject::Qoi(ContextQoiKey {
                        context: experiment.context.id().clone(),
                        qoi: pair[1].clone(),
                    })),
                ));
            }
        }
        let submitted_qois: BTreeSet<QoiId> = experiment.qois.iter().cloned().collect();
        if let Some(evidence) = evidence {
            if submitted_qois != evidence.qois {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::ExperimentQoiSetMismatch,
                    subject(),
                    None,
                ));
            }
        }
        for qoi in &experiment.qois {
            let key = ContextQoiKey {
                context: experiment.context.id().clone(),
                qoi: qoi.clone(),
            };
            if !qois.contains_key(&key) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownExperimentQoi,
                    subject(),
                    Some(MachineAssuranceSubject::Qoi(key)),
                ));
            } else {
                qoi_coverage.insert(key);
            }
        }
    }
    for key in qois.keys() {
        if !qoi_coverage.contains(key) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::MissingQoiExperiment,
                MachineAssuranceSubject::Qoi(key.clone()),
                None,
            ));
        }
    }
    for (context_id, case) in vv_index {
        for experiment_id in case.experiments.keys() {
            if !artifact_coverage.contains(&(context_id.clone(), experiment_id.clone())) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::MissingExperimentBinding,
                    MachineAssuranceSubject::Artifact(experiment_id.clone()),
                    Some(MachineAssuranceSubject::Context(case.context.clone())),
                ));
            }
        }
    }
}

fn validate_hazards_and_faults(
    hazards: &[HazardSpec],
    faults: &[FaultSpec],
    contexts: &BTreeMap<ArtifactId, ArtifactRef>,
    index: &AssuranceIndex<'_>,
    vv_index: &BTreeMap<ArtifactId, VvCaseView>,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    if hazards.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptyHazards,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut hazard_map = BTreeMap::<HazardId, &HazardSpec>::new();
    for hazard in hazards {
        let subject = || MachineAssuranceSubject::Hazard(hazard.id.clone());
        if let Some(first) = hazard_map.insert(hazard.id.clone(), hazard) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateHazard,
                subject(),
                Some(MachineAssuranceSubject::Hazard(first.id.clone())),
            ));
        }
        match contexts.get(hazard.context.id()) {
            Some(context) if context == &hazard.context => {}
            _ => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownHazardContext,
                subject(),
                Some(MachineAssuranceSubject::Context(hazard.context.clone())),
            )),
        }
        if hazard.scope.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyHazardScope,
                subject(),
                None,
            ));
        }
        for pair in hazard.scope.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateHazardElement,
                    subject(),
                    None,
                ));
            }
        }
        for scope in &hazard.scope {
            if !index.scope_exists(scope) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownHazardElement,
                    subject(),
                    machine_scope_subject(scope),
                ));
            }
        }
        if hazard.assumptions.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyHazardAssumptions,
                subject(),
                None,
            ));
        }
        for pair in hazard.assumptions.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateHazardAssumption,
                    subject(),
                    Some(MachineAssuranceSubject::Assumption(pair[1].clone())),
                ));
            }
        }
        let admitted_assumptions = vv_index
            .get(hazard.context.id())
            .map(|case| &case.assumptions);
        for assumption in &hazard.assumptions {
            if admitted_assumptions.is_none_or(|rows| !rows.contains(assumption)) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownHazardAssumption,
                    subject(),
                    Some(MachineAssuranceSubject::Assumption(assumption.clone())),
                ));
            }
        }
    }

    let mut fault_map = BTreeMap::<FaultId, &FaultSpec>::new();
    let mut hazard_coverage = BTreeSet::<HazardId>::new();
    for fault in faults {
        let subject = || MachineAssuranceSubject::Fault(fault.id.clone());
        if let Some(first) = fault_map.insert(fault.id.clone(), fault) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateFault,
                subject(),
                Some(MachineAssuranceSubject::Fault(first.id.clone())),
            ));
        }
        if fault.affected.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyFaultAffected,
                subject(),
                None,
            ));
        }
        for pair in fault.affected.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateFaultElement,
                    subject(),
                    None,
                ));
            }
        }
        for scope in &fault.affected {
            if !index.scope_exists(scope) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownFaultElement,
                    subject(),
                    machine_scope_subject(scope),
                ));
            }
        }
        if fault.hazards.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyFaultHazards,
                subject(),
                None,
            ));
        }
        for pair in fault.hazards.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateFaultHazard,
                    subject(),
                    Some(MachineAssuranceSubject::Hazard(pair[1].clone())),
                ));
            }
        }
        for hazard_id in &fault.hazards {
            match hazard_map.get(hazard_id) {
                None => findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownFaultHazard,
                    subject(),
                    Some(MachineAssuranceSubject::Hazard(hazard_id.clone())),
                )),
                Some(hazard) => {
                    hazard_coverage.insert(hazard_id.clone());
                    if matches!(hazard.fault_coverage, FaultCoverage::Unmodeled(_)) {
                        findings.push(MachineAssuranceFinding::new(
                            MachineAssuranceRule::ContradictoryFaultCoverage,
                            subject(),
                            Some(MachineAssuranceSubject::Hazard(hazard_id.clone())),
                        ));
                    }
                }
            }
        }
    }
    for hazard in hazards {
        let covered = hazard_coverage.contains(&hazard.id);
        match (&hazard.fault_coverage, covered) {
            (FaultCoverage::Modeled, false) => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UncoveredModeledHazard,
                MachineAssuranceSubject::Hazard(hazard.id.clone()),
                None,
            )),
            (FaultCoverage::Unmodeled(_), true) => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::ContradictoryFaultCoverage,
                MachineAssuranceSubject::Hazard(hazard.id.clone()),
                None,
            )),
            _ => {}
        }
    }
}

fn machine_scope_subject(scope: &MachineScope) -> Option<MachineAssuranceSubject> {
    match scope {
        MachineScope::WholeMachine => None,
        MachineScope::Subsystem(id) => Some(MachineAssuranceSubject::Subsystem(id.clone())),
        MachineScope::Element(id) => Some(MachineAssuranceSubject::Element(id.clone())),
        MachineScope::Relation(id) => Some(MachineAssuranceSubject::AccountingTarget(
            AccountingTarget::Relation(id.clone()),
        )),
        MachineScope::Interface(id) => Some(MachineAssuranceSubject::AccountingTarget(
            AccountingTarget::Interface(id.clone()),
        )),
    }
}

fn validate_accounting(
    windows: &[AccountingWindow],
    contexts: &BTreeMap<ArtifactId, ArtifactRef>,
    index: &AssuranceIndex<'_>,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    if windows.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptyAccountingWindows,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut by_id = BTreeMap::<AccountingWindowId, &AccountingWindow>::new();
    for window in windows {
        let subject = || MachineAssuranceSubject::AccountingWindow(window.id.clone());
        if let Some(first) = by_id.insert(window.id.clone(), window) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateAccountingWindow,
                subject(),
                Some(MachineAssuranceSubject::AccountingWindow(first.id.clone())),
            ));
        }
        match contexts.get(window.context.id()) {
            Some(context) if context == &window.context => {}
            _ => findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownAccountingContext,
                subject(),
                Some(MachineAssuranceSubject::Context(window.context.clone())),
            )),
        }
        if !index.clocks.contains(&window.clock) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownAccountingClock,
                subject(),
                Some(MachineAssuranceSubject::Clock(window.clock.clone())),
            ));
        }
        if !window.quantity.is_admitted() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnsupportedAccountingQuantity,
                subject(),
                None,
            ));
        }
        if !balance_quantity_is_compatible(&window.balance, window.quantity.dims()) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::AccountingBalanceQuantityGap,
                subject(),
                None,
            ));
        }
        if window.entries.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyAccountingEntries,
                subject(),
                None,
            ));
        }
        let mut seen_entries = BTreeSet::new();
        let mut claimed_atoms = BTreeMap::<AccountingAtom, AccountingTarget>::new();
        let mut loss_owners = BTreeSet::new();
        for entry in &window.entries {
            if !seen_entries.insert(entry.target.clone()) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateAccountingEntry,
                    subject(),
                    Some(MachineAssuranceSubject::AccountingTarget(
                        entry.target.clone(),
                    )),
                ));
            }
            if !index.accounting_target_exists(&entry.target) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownAccountingTarget,
                    subject(),
                    Some(MachineAssuranceSubject::AccountingTarget(
                        entry.target.clone(),
                    )),
                ));
            } else {
                for atom in index.accounting_atoms(&entry.target, entry.role) {
                    if let Some(first) = claimed_atoms.insert(atom, entry.target.clone()) {
                        if first != entry.target {
                            findings.push(MachineAssuranceFinding::new(
                                MachineAssuranceRule::OverlappingAccountingTarget,
                                subject(),
                                Some(MachineAssuranceSubject::AccountingTarget(first)),
                            ));
                        }
                    }
                }
            }
            if !accounting_role_is_valid(entry) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::InvalidAccountingRole,
                    subject(),
                    Some(MachineAssuranceSubject::AccountingTarget(
                        entry.target.clone(),
                    )),
                ));
            }
            match (&entry.role, &entry.loss_ownership) {
                (AccountingRole::Dissipation, None) => {
                    findings.push(MachineAssuranceFinding::new(
                        MachineAssuranceRule::MissingLossOwnership,
                        subject(),
                        Some(MachineAssuranceSubject::AccountingTarget(
                            entry.target.clone(),
                        )),
                    ));
                }
                (AccountingRole::Dissipation, Some(owner)) => {
                    if !loss_owners.insert(owner.clone()) {
                        findings.push(MachineAssuranceFinding::new(
                            MachineAssuranceRule::DuplicateLossOwnership,
                            subject(),
                            Some(MachineAssuranceSubject::AccountingTarget(
                                entry.target.clone(),
                            )),
                        ));
                    }
                }
                (_, Some(_)) => findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnexpectedLossOwnership,
                    subject(),
                    Some(MachineAssuranceSubject::AccountingTarget(
                        entry.target.clone(),
                    )),
                )),
                (_, None) => {}
            }
        }
    }
}

fn balance_quantity_is_compatible(balance: &BalanceKind, dims: Dims) -> bool {
    match balance {
        BalanceKind::Energy | BalanceKind::Enthalpy | BalanceKind::Exergy => {
            dims == Dims([2, 1, -2, 0, 0, 0])
        }
        BalanceKind::LinearMomentum => dims == Dims([1, 1, -1, 0, 0, 0]),
        BalanceKind::AngularMomentum => dims == Dims([2, 1, -1, 0, 0, 0]),
        BalanceKind::Mass => dims == Dims([0, 1, 0, 0, 0, 0]),
        BalanceKind::ElectricCharge => dims == Dims([0, 0, 1, 1, 0, 0]),
        BalanceKind::AmountOfSubstance => dims == Dims([0, 0, 0, 0, 0, 1]),
        BalanceKind::Entropy => dims == Dims([2, 1, -2, 0, -1, 0]),
        BalanceKind::Species(_) | BalanceKind::Elements(_) | BalanceKind::Custom(_) => true,
    }
}

fn accounting_role_is_valid(entry: &AccountingEntry) -> bool {
    match (&entry.role, &entry.target, entry.orientation) {
        (
            AccountingRole::Storage,
            AccountingTarget::Relation(_) | AccountingTarget::State(_),
            AccountingOrientation::StoredIncreasePositive,
        ) => true,
        (
            AccountingRole::Dissipation,
            AccountingTarget::Relation(_),
            AccountingOrientation::NonnegativeLoss,
        ) => true,
        (
            AccountingRole::IncludedSource,
            AccountingTarget::Relation(_),
            AccountingOrientation::IntoBoundaryPositive
            | AccountingOrientation::OutOfBoundaryPositive,
        ) => true,
        (
            AccountingRole::ExternalExchange | AccountingRole::Stream,
            AccountingTarget::Port(_) | AccountingTarget::Interface(_),
            AccountingOrientation::IntoBoundaryPositive
            | AccountingOrientation::OutOfBoundaryPositive,
        ) => true,
        _ => false,
    }
}

fn validate_fidelity(
    policy: &FidelityPolicy,
    qois: &BTreeMap<ContextQoiKey, &QoiBinding>,
    index: &AssuranceIndex<'_>,
    findings: &mut Vec<MachineAssuranceFinding>,
) {
    if policy.rungs.is_empty() {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::EmptyFidelityRungs,
            MachineAssuranceSubject::Assurance,
            None,
        ));
    }
    let mut rungs = BTreeMap::<FidelityRungId, &FidelityRung>::new();
    let mut qoi_coverage = BTreeSet::<ContextQoiKey>::new();
    for rung in &policy.rungs {
        let subject = || MachineAssuranceSubject::FidelityRung(rung.id.clone());
        if let Some(first) = rungs.insert(rung.id.clone(), rung) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateFidelityRung,
                subject(),
                Some(MachineAssuranceSubject::FidelityRung(first.id.clone())),
            ));
        }
        if !index.subsystems.contains_key(&rung.subsystem) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownFidelitySubsystem,
                subject(),
                Some(MachineAssuranceSubject::Subsystem(rung.subsystem.clone())),
            ));
        }
        if rung.falsifiers.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyFidelityFalsifiers,
                subject(),
                None,
            ));
        }
        for pair in rung.falsifiers.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateFidelityFalsifier,
                    subject(),
                    None,
                ));
            }
        }
        if rung.qois.is_empty() {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::EmptyFidelityQois,
                subject(),
                None,
            ));
        }
        for pair in rung.qois.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::DuplicateFidelityQoi,
                    subject(),
                    Some(MachineAssuranceSubject::Qoi(pair[1].clone())),
                ));
            }
        }
        for qoi in &rung.qois {
            if !qois.contains_key(qoi) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownFidelityQoi,
                    subject(),
                    Some(MachineAssuranceSubject::Qoi(qoi.clone())),
                ));
            } else {
                qoi_coverage.insert(qoi.clone());
            }
        }
    }
    for qoi in qois.keys() {
        if !qoi_coverage.contains(qoi) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::MissingQoiFidelityCoverage,
                MachineAssuranceSubject::Qoi(qoi.clone()),
                None,
            ));
        }
    }

    let mut baseline_ids = BTreeSet::new();
    let mut baseline_by_subsystem = BTreeMap::<SubsystemId, FidelityRungId>::new();
    for baseline in &policy.baselines {
        if !baseline_ids.insert(baseline.clone()) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateFidelityBaseline,
                MachineAssuranceSubject::FidelityRung(baseline.clone()),
                None,
            ));
        }
        let Some(rung) = rungs.get(baseline) else {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownFidelityBaseline,
                MachineAssuranceSubject::FidelityRung(baseline.clone()),
                None,
            ));
            continue;
        };
        if let Some(first) = baseline_by_subsystem.insert(rung.subsystem.clone(), baseline.clone())
        {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateSubsystemBaseline,
                MachineAssuranceSubject::Subsystem(rung.subsystem.clone()),
                Some(MachineAssuranceSubject::FidelityRung(first)),
            ));
        }
        if let Some(subsystem) = index.subsystems.get(&rung.subsystem) {
            if rung.model != subsystem.model {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::BaselineModelMismatch,
                    MachineAssuranceSubject::FidelityRung(rung.id.clone()),
                    Some(MachineAssuranceSubject::Subsystem(rung.subsystem.clone())),
                ));
            }
        }
    }
    for subsystem in index.subsystems.keys() {
        if !baseline_by_subsystem.contains_key(subsystem) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::MissingSubsystemBaseline,
                MachineAssuranceSubject::Subsystem(subsystem.clone()),
                None,
            ));
        }
    }

    let mut outgoing = BTreeMap::<FidelityRungId, Option<FidelityRungId>>::new();
    for escalation in &policy.escalations {
        let subject = || MachineAssuranceSubject::FidelityRung(escalation.from.clone());
        let next = match &escalation.action {
            EscalationAction::Escalate { target, .. } => Some(target.clone()),
            EscalationAction::Refuse(_) => None,
        };
        if outgoing
            .insert(escalation.from.clone(), next.clone())
            .is_some()
        {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::DuplicateEscalation,
                subject(),
                None,
            ));
        }
        let Some(source) = rungs.get(&escalation.from) else {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnknownEscalationRung,
                subject(),
                None,
            ));
            continue;
        };
        if let Some(target_id) = next {
            let Some(target) = rungs.get(&target_id) else {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::UnknownEscalationRung,
                    subject(),
                    Some(MachineAssuranceSubject::FidelityRung(target_id)),
                ));
                continue;
            };
            if escalation.from == target.id {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::SelfEscalation,
                    subject(),
                    None,
                ));
            }
            if source.subsystem != target.subsystem {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::CrossSubsystemEscalation,
                    subject(),
                    Some(MachineAssuranceSubject::FidelityRung(target.id.clone())),
                ));
            }
            if source.qois.iter().any(|qoi| !target.qois.contains(qoi)) {
                findings.push(MachineAssuranceFinding::new(
                    MachineAssuranceRule::FidelityQoiDrop,
                    subject(),
                    Some(MachineAssuranceSubject::FidelityRung(target.id.clone())),
                ));
            }
        }
    }
    for rung in rungs.keys() {
        if !outgoing.contains_key(rung) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::MissingEscalation,
                MachineAssuranceSubject::FidelityRung(rung.clone()),
                None,
            ));
        }
    }

    let mut cyclic = BTreeSet::new();
    for start in rungs.keys() {
        let mut seen = BTreeSet::new();
        let mut current = start.clone();
        loop {
            if !seen.insert(current.clone()) {
                cyclic.extend(seen);
                break;
            }
            match outgoing.get(&current) {
                Some(Some(next)) if rungs.contains_key(next) => current = next.clone(),
                Some(None) | None | Some(Some(_)) => break,
            }
        }
    }
    for rung in cyclic {
        findings.push(MachineAssuranceFinding::new(
            MachineAssuranceRule::FidelityEscalationCycle,
            MachineAssuranceSubject::FidelityRung(rung),
            None,
        ));
    }

    let mut reachable = BTreeSet::new();
    for baseline in baseline_ids {
        let mut current = baseline;
        while reachable.insert(current.clone()) {
            match outgoing.get(&current) {
                Some(Some(next)) if rungs.contains_key(next) => current = next.clone(),
                _ => break,
            }
        }
    }
    for rung in rungs.keys() {
        if !reachable.contains(rung) {
            findings.push(MachineAssuranceFinding::new(
                MachineAssuranceRule::UnreachableFidelityRung,
                MachineAssuranceSubject::FidelityRung(rung.clone()),
                None,
            ));
        }
    }
}

fn machine_assurance_identity(
    draft: &MachineAssuranceDraft,
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    vv_bindings: &[VvCaseBinding],
) -> Result<IdentityReceipt<MachineAssuranceIdV1>, CanonicalError> {
    let sensor_rows: Vec<Vec<u8>> = draft.sensors.iter().map(sensor_row).collect();
    let experiment_rows: Vec<Vec<u8>> = draft.experiments.iter().map(experiment_row).collect();
    let context_rows: Vec<Vec<u8>> = draft.contexts.iter().map(context_row).collect();
    let vv_rows: Vec<Vec<u8>> = vv_bindings.iter().map(vv_case_row).collect();
    let hazard_rows: Vec<Vec<u8>> = draft.hazards.iter().map(hazard_row).collect();
    let fault_rows: Vec<Vec<u8>> = draft.faults.iter().map(fault_row).collect();
    let accounting_rows: Vec<Vec<u8>> = draft
        .accounting_windows
        .iter()
        .map(accounting_window_row)
        .collect();
    let baseline_rows: Vec<Vec<u8>> = draft
        .fidelity
        .baselines
        .iter()
        .map(|id| id.digest_bytes().to_vec())
        .collect();
    let rung_rows: Vec<Vec<u8>> = draft.fidelity.rungs.iter().map(fidelity_rung_row).collect();
    let escalation_rows: Vec<Vec<u8>> = draft
        .fidelity
        .escalations
        .iter()
        .map(escalation_row)
        .collect();
    let mut fixed_replay = Vec::new();
    draft
        .fidelity
        .fixed_replay
        .append_canonical(&mut fixed_replay);

    CanonicalEncoder::<MachineAssuranceIdV1, _>::new(
        MACHINE_ASSURANCE_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "assurance-schema-version"),
        u64::from(MACHINE_ASSURANCE_SCHEMA_VERSION_V1),
    )?
    .bytes(Field::new(1, "base-machine-graph"), graph.as_bytes())?
    .bytes(Field::new(2, "base-machine-behavior"), behavior.as_bytes())?
    .ordered_bytes(
        Field::new(3, "sensors"),
        sensor_rows.len() as u64,
        sensor_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(4, "experiments"),
        experiment_rows.len() as u64,
        experiment_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(5, "contexts"),
        context_rows.len() as u64,
        context_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(6, "vv-case-receipts"),
        vv_rows.len() as u64,
        vv_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(7, "hazards"),
        hazard_rows.len() as u64,
        hazard_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(8, "faults"),
        fault_rows.len() as u64,
        fault_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(9, "accounting-windows"),
        accounting_rows.len() as u64,
        accounting_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(10, "fidelity-baselines"),
        baseline_rows.len() as u64,
        baseline_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(11, "fidelity-rungs"),
        rung_rows.len() as u64,
        rung_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(12, "fidelity-escalations"),
        escalation_rows.len() as u64,
        escalation_rows.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(13, "fixed-fidelity-replay"), &fixed_replay)?
    .finish()
}

fn sensor_row(sensor: &SensorSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(448);
    push_identity(&mut out, &sensor.id.digest_bytes());
    push_identity(&mut out, &sensor.owner.digest_bytes());
    push_observation_target(&mut out, &sensor.target);
    super::push_terminal_quantity(&mut out, sensor.quantity);
    super::push_terminal_shape(&mut out, sensor.shape);
    push_identity(&mut out, &sensor.clock.digest_bytes());
    push_frame(&mut out, &sensor.frame);
    match &sensor.timing {
        ObservationTiming::Direct => out.push(1),
        ObservationTiming::ModeledResampling { bridge } => {
            out.push(2);
            bridge.append_canonical(&mut out);
        }
    }
    sensor.model.append_canonical(&mut out);
    sensor.calibration.append_canonical(&mut out);
    match &sensor.exposure {
        SensorExposure::PlantSignal { output } => {
            out.push(1);
            push_identity(&mut out, &output.digest_bytes());
        }
        SensorExposure::ExperimentOnly => out.push(2),
    }
    out
}

fn push_observation_target(out: &mut Vec<u8>, target: &ObservationTarget) {
    match target {
        ObservationTarget::Terminal(id) => {
            out.push(1);
            push_identity(out, &id.digest_bytes());
        }
        ObservationTarget::State(id) => {
            out.push(2);
            push_identity(out, &id.digest_bytes());
        }
    }
}

fn experiment_row(experiment: &ExperimentSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(256 + 96 * experiment.instruments.len());
    push_identity(&mut out, &experiment.id.digest_bytes());
    push_artifact_ref(&mut out, &experiment.artifact);
    push_artifact_ref(&mut out, &experiment.context);
    out.extend_from_slice(&(experiment.instruments.len() as u64).to_le_bytes());
    for binding in &experiment.instruments {
        push_identity(&mut out, &binding.sensor.digest_bytes());
        push_len_prefixed(&mut out, binding.instrument.as_str().as_bytes());
    }
    out.extend_from_slice(&(experiment.qois.len() as u64).to_le_bytes());
    for qoi in &experiment.qois {
        push_len_prefixed(&mut out, qoi.as_str().as_bytes());
    }
    out
}

fn qoi_input_row(input: &QoiInput) -> Vec<u8> {
    let mut out = Vec::with_capacity(96);
    match &input.target {
        QoiTarget::Sensor(id) => {
            out.push(1);
            push_identity(&mut out, &id.digest_bytes());
        }
        QoiTarget::Terminal(id) => {
            out.push(2);
            push_identity(&mut out, &id.digest_bytes());
        }
        QoiTarget::State(id) => {
            out.push(3);
            push_identity(&mut out, &id.digest_bytes());
        }
    }
    super::push_terminal_quantity(&mut out, input.quantity);
    super::push_terminal_shape(&mut out, input.shape);
    out
}

fn qoi_binding_row(qoi: &QoiBinding) -> Vec<u8> {
    let input_rows: Vec<Vec<u8>> = qoi.inputs.iter().map(qoi_input_row).collect();
    let mut out = Vec::with_capacity(256 + input_rows.iter().map(Vec::len).sum::<usize>());
    push_len_prefixed(&mut out, qoi.id.as_str().as_bytes());
    push_len_prefixed(&mut out, qoi.unit.as_str().as_bytes());
    out.extend_from_slice(&(input_rows.len() as u64).to_le_bytes());
    for row in input_rows {
        push_len_prefixed(&mut out, &row);
    }
    qoi.definition.append_canonical(&mut out);
    qoi.unit_bridge.append_canonical(&mut out);
    out
}

fn context_row(context: &ContextBinding) -> Vec<u8> {
    let qoi_rows: Vec<Vec<u8>> = context.qois.iter().map(qoi_binding_row).collect();
    let mut out = Vec::with_capacity(256 + qoi_rows.iter().map(Vec::len).sum::<usize>());
    push_artifact_ref(&mut out, &context.context);
    push_artifact_ref(&mut out, &context.validation_plan);
    out.extend_from_slice(&(qoi_rows.len() as u64).to_le_bytes());
    for row in qoi_rows {
        push_len_prefixed(&mut out, &row);
    }
    context.budget.append_canonical(&mut out);
    out
}

fn vv_case_row(binding: &VvCaseBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_artifact_ref(&mut out, &binding.context);
    out.extend_from_slice(&binding.schema_version.to_le_bytes());
    out.extend_from_slice(&binding.ruleset_version.to_le_bytes());
    out.extend_from_slice(binding.case_hash.as_bytes());
    out.extend_from_slice(binding.receipt_hash.as_bytes());
    out
}

fn hazard_row(hazard: &HazardSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(512 + 64 * hazard.scope.len());
    push_identity(&mut out, &hazard.id.digest_bytes());
    push_artifact_ref(&mut out, &hazard.context);
    out.extend_from_slice(&(hazard.scope.len() as u64).to_le_bytes());
    for scope in &hazard.scope {
        push_machine_scope(&mut out, scope);
    }
    hazard.requirement.append_canonical(&mut out);
    hazard.operating_envelope.append_canonical(&mut out);
    hazard.safety_case.append_canonical(&mut out);
    out.extend_from_slice(&(hazard.assumptions.len() as u64).to_le_bytes());
    for assumption in &hazard.assumptions {
        push_len_prefixed(&mut out, assumption.as_str().as_bytes());
    }
    match &hazard.fault_coverage {
        FaultCoverage::Modeled => out.push(1),
        FaultCoverage::Unmodeled(no_claim) => {
            out.push(2);
            no_claim.append_canonical(&mut out);
        }
    }
    out
}

fn fault_row(fault: &FaultSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(384 + 64 * fault.affected.len());
    push_identity(&mut out, &fault.id.digest_bytes());
    out.extend_from_slice(&(fault.affected.len() as u64).to_le_bytes());
    for scope in &fault.affected {
        push_machine_scope(&mut out, scope);
    }
    out.extend_from_slice(&(fault.hazards.len() as u64).to_le_bytes());
    for hazard in &fault.hazards {
        push_identity(&mut out, &hazard.digest_bytes());
    }
    fault.model.append_canonical(&mut out);
    fault.containment.append_canonical(&mut out);
    fault.injection.append_canonical(&mut out);
    out
}

fn push_machine_scope(out: &mut Vec<u8>, scope: &MachineScope) {
    match scope {
        MachineScope::WholeMachine => out.push(1),
        MachineScope::Subsystem(id) => {
            out.push(2);
            push_identity(out, &id.digest_bytes());
        }
        MachineScope::Element(id) => {
            out.push(3);
            out.push(id.kind().tag());
            push_identity(out, &id.digest_bytes());
        }
        MachineScope::Relation(id) => {
            out.push(4);
            push_identity(out, &id.digest_bytes());
        }
        MachineScope::Interface(id) => {
            out.push(5);
            push_identity(out, &id.digest_bytes());
        }
    }
}

fn accounting_entry_row(entry: &AccountingEntry) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    push_accounting_target(&mut out, &entry.target);
    out.push(accounting_role_tag(entry.role));
    out.push(accounting_orientation_tag(entry.orientation));
    entry.policy.append_canonical(&mut out);
    match &entry.loss_ownership {
        Some(owner) => {
            out.push(1);
            owner.append_canonical(&mut out);
        }
        None => out.push(0),
    }
    out
}

fn accounting_window_row(window: &AccountingWindow) -> Vec<u8> {
    let entry_rows: Vec<Vec<u8>> = window.entries.iter().map(accounting_entry_row).collect();
    let mut out = Vec::with_capacity(384 + entry_rows.iter().map(Vec::len).sum::<usize>());
    push_identity(&mut out, &window.id.digest_bytes());
    push_artifact_ref(&mut out, &window.context);
    push_identity(&mut out, &window.clock.digest_bytes());
    push_balance_kind(&mut out, &window.balance);
    super::push_terminal_quantity(&mut out, window.quantity);
    window.boundary.append_canonical(&mut out);
    window.interval.append_canonical(&mut out);
    out.extend_from_slice(&(entry_rows.len() as u64).to_le_bytes());
    for row in entry_rows {
        push_len_prefixed(&mut out, &row);
    }
    window.audit_policy.append_canonical(&mut out);
    out
}

fn push_accounting_target(out: &mut Vec<u8>, target: &AccountingTarget) {
    match target {
        AccountingTarget::Relation(id) => {
            out.push(1);
            push_identity(out, &id.digest_bytes());
        }
        AccountingTarget::Port(id) => {
            out.push(2);
            push_identity(out, &id.digest_bytes());
        }
        AccountingTarget::Interface(id) => {
            out.push(3);
            push_identity(out, &id.digest_bytes());
        }
        AccountingTarget::State(id) => {
            out.push(4);
            push_identity(out, &id.digest_bytes());
        }
    }
}

fn push_balance_kind(out: &mut Vec<u8>, balance: &BalanceKind) {
    match balance {
        BalanceKind::Energy => out.push(1),
        BalanceKind::Enthalpy => out.push(2),
        BalanceKind::LinearMomentum => out.push(3),
        BalanceKind::AngularMomentum => out.push(4),
        BalanceKind::Mass => out.push(5),
        BalanceKind::ElectricCharge => out.push(6),
        BalanceKind::AmountOfSubstance => out.push(7),
        BalanceKind::Species(law) => {
            out.push(8);
            law.append_canonical(out);
        }
        BalanceKind::Elements(law) => {
            out.push(9);
            law.append_canonical(out);
        }
        BalanceKind::Entropy => out.push(10),
        BalanceKind::Exergy => out.push(11),
        BalanceKind::Custom(law) => {
            out.push(12);
            law.append_canonical(out);
        }
    }
}

fn accounting_role_tag(role: AccountingRole) -> u8 {
    match role {
        AccountingRole::Storage => 1,
        AccountingRole::Dissipation => 2,
        AccountingRole::IncludedSource => 3,
        AccountingRole::ExternalExchange => 4,
        AccountingRole::Stream => 5,
    }
}

fn accounting_orientation_tag(orientation: AccountingOrientation) -> u8 {
    match orientation {
        AccountingOrientation::StoredIncreasePositive => 1,
        AccountingOrientation::NonnegativeLoss => 2,
        AccountingOrientation::IntoBoundaryPositive => 3,
        AccountingOrientation::OutOfBoundaryPositive => 4,
    }
}

fn context_qoi_row(key: &ContextQoiKey) -> Vec<u8> {
    let mut out = Vec::with_capacity(96);
    push_len_prefixed(&mut out, key.context.as_str().as_bytes());
    push_len_prefixed(&mut out, key.qoi.as_str().as_bytes());
    out
}

fn fidelity_rung_row(rung: &FidelityRung) -> Vec<u8> {
    let mut out = Vec::with_capacity(512 + 96 * rung.qois.len());
    push_identity(&mut out, &rung.id.digest_bytes());
    push_identity(&mut out, &rung.subsystem.digest_bytes());
    rung.model.append_canonical(&mut out);
    rung.model_crosswalk.append_canonical(&mut out);
    rung.validity_domain.append_canonical(&mut out);
    rung.cost_error_model.append_canonical(&mut out);
    out.extend_from_slice(&(rung.falsifiers.len() as u64).to_le_bytes());
    for falsifier in &rung.falsifiers {
        falsifier.append_canonical(&mut out);
    }
    out.extend_from_slice(&(rung.qois.len() as u64).to_le_bytes());
    for qoi in &rung.qois {
        let row = context_qoi_row(qoi);
        push_len_prefixed(&mut out, &row);
    }
    out
}

fn escalation_row(escalation: &EscalationSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(384);
    push_identity(&mut out, &escalation.from.digest_bytes());
    escalation.trigger.append_canonical(&mut out);
    match &escalation.action {
        EscalationAction::Escalate {
            target,
            transfer,
            crosswalk,
        } => {
            out.push(1);
            push_identity(&mut out, &target.digest_bytes());
            transfer.append_canonical(&mut out);
            crosswalk.append_canonical(&mut out);
        }
        EscalationAction::Refuse(no_claim) => {
            out.push(2);
            no_claim.append_canonical(&mut out);
        }
    }
    out
}

fn push_artifact_ref(out: &mut Vec<u8>, reference: &ArtifactRef) {
    out.push(artifact_kind_tag(reference.kind()));
    push_len_prefixed(out, reference.id().as_str().as_bytes());
    out.extend_from_slice(reference.hash().as_bytes());
}

fn artifact_kind_tag(kind: ArtifactKind) -> u8 {
    match kind {
        ArtifactKind::ContextOfUse => 1,
        ArtifactKind::ValidationPlan => 2,
        ArtifactKind::ExperimentArtifact => 3,
        ArtifactKind::CalibrationSplit => 4,
        ArtifactKind::SolutionVerificationReceipt => 5,
        ArtifactKind::PredictionAssessment => 6,
        ArtifactKind::AssumptionsLedger => 7,
    }
}

fn push_frame(out: &mut Vec<u8>, frame: &FrameBinding) {
    push_len_prefixed(out, frame.canonical_key().as_bytes());
    out.push(super::orientation_parity_tag(frame.orientation()));
}

fn push_len_prefixed(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

fn push_identity(out: &mut Vec<u8>, identity: &[u8; 32]) {
    out.extend_from_slice(identity);
}

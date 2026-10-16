//! Admitted PR-3 behavior semantics layered over an admitted machine graph.
//!
//! This module deliberately does not extend the canonical grammar of
//! `MachineGraphIdV1`. A behavior draft binds that already-published identity
//! and then adds state contracts, initial and boundary sources, explicit body
//! motion, event/reset structure, tolerances, and declared dependence. Actual
//! values, laws, guards, reset maps, and statistical models remain nominal,
//! versioned references owned by their execution or evidence domains.
//!
//! Admission proves bounded structural closure only. It does not prove PDE or
//! DAE well-posedness, swept-volume validity, true-flow event completeness,
//! reset regularity, absence of grazing or Zeno behavior, covariance PSD, or
//! physical validation of a referenced distribution/correlation model.
//! Collision policy is closed within each declared event clock; synchronization
//! or simultaneity between distinct clock domains is deliberately unclaimed.

use core::fmt;
use core::hash::{Hash, Hasher};
use core::num::NonZeroU64;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EntityId, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use super::{
    AdmittedMachineGraph, BodyId, ClockId, FrameBinding, MachineClock, MachineElementId,
    MachineGraphIdV1, MachineIdError, StateSlotId, SubsystemId, TerminalCausality, TerminalId,
    TerminalQuantitySpec, TerminalShape,
};

/// Version of the Machine-IR behavior-overlay identity schema.
pub const MACHINE_BEHAVIOR_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum state contracts in one behavior draft.
pub const MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS: usize = 4_096;
/// Maximum initial plus boundary conditions in one behavior draft.
pub const MAX_MACHINE_BEHAVIOR_CONDITIONS: usize = 8_192;
/// Maximum explicit body-motion bindings in one behavior draft.
pub const MAX_MACHINE_BEHAVIOR_MOTIONS: usize = 4_096;
/// Maximum events in one behavior draft.
pub const MAX_MACHINE_BEHAVIOR_EVENTS: usize = 4_096;
/// Maximum tolerance declarations in one behavior draft.
pub const MAX_MACHINE_BEHAVIOR_TOLERANCES: usize = 8_192;
/// Maximum submitted dependence declarations before semantic single-model checks.
pub const MAX_MACHINE_BEHAVIOR_DEPENDENCES: usize = 4_096;
/// Maximum aggregate nested references inspected by one admission.
pub const MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES: usize = 32_768;

const MACHINE_BEHAVIOR_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(24 * 1_024 * 1_024, 8 * 1_024 * 1_024, 8, 96_000, 32_768);

macro_rules! behavior_id {
    ($(#[$meta:meta])* $name:ident, $schema:ident, $identity:ident, $role:literal, $domain:literal, $context:literal) => {
        #[doc = concat!("Canonical schema marker for `", stringify!($name), "`.")]
        pub enum $schema {}

        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str = $domain;
            const NAME: &'static str = $role;
            const VERSION: u32 = MACHINE_BEHAVIOR_SCHEMA_VERSION_V1;
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

behavior_id!(
    /// Durable identity of one machine event declaration.
    EventId,
    EventIdSchemaV1,
    EventEntityIdV1,
    "event-id",
    "org.frankensim.fs-ir.machine.event-id.v1",
    "one Machine-IR event independent of collection and execution order"
);
behavior_id!(
    /// Durable identity of one tolerance declaration.
    ToleranceId,
    ToleranceIdSchemaV1,
    ToleranceEntityIdV1,
    "tolerance-id",
    "org.frankensim.fs-ir.machine.tolerance-id.v1",
    "one Machine-IR tolerance independent of sampling and serialization order"
);

macro_rules! behavior_ref {
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
            /// The referenced domain retains all execution and validation
            /// authority. Machine-IR binds the supplied identity exactly.
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

behavior_ref!(
    /// Opaque fixed condition-value artifact.
    ConditionValueRef,
    "condition-value-ref"
);
behavior_ref!(
    /// Opaque time-history or signal artifact.
    ConditionHistoryRef,
    "condition-history-ref"
);
behavior_ref!(
    /// Opaque probability-distribution artifact.
    DistributionRef,
    "distribution-ref"
);
behavior_ref!(
    /// Opaque body-motion path or law.
    MotionPathRef,
    "motion-path-ref"
);
behavior_ref!(
    /// Opaque event guard artifact.
    GuardRef,
    "guard-ref"
);
behavior_ref!(
    /// Opaque event/reset witness artifact.
    EventWitnessRef,
    "event-witness-ref"
);
behavior_ref!(
    /// Opaque explicit no-claim artifact.
    NoClaimRef,
    "no-claim-ref"
);
behavior_ref!(
    /// Opaque deterministic or set-valued reset relation.
    ResetMapRef,
    "reset-map-ref"
);
behavior_ref!(
    /// Opaque set-valued reset outcome set.
    OutcomeSetRef,
    "outcome-set-ref"
);
behavior_ref!(
    /// Opaque simultaneity-group semantics.
    SimultaneityGroupRef,
    "simultaneity-group-ref"
);
behavior_ref!(
    /// Opaque model/material parameter identity.
    ParameterRef,
    "parameter-ref"
);
behavior_ref!(
    /// Opaque bounded or random tolerance law.
    ToleranceLawRef,
    "tolerance-law-ref"
);
behavior_ref!(
    /// Opaque externally owned correlation/dependence model.
    CorrelationModelRef,
    "correlation-model-ref"
);

/// Refusal from constructing an identity-safe nonnegative scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteNonNegativeError {
    /// NaN and infinities are not semantic scalar values.
    NonFinite,
    /// A tolerance magnitude cannot be negative.
    Negative,
}

impl fmt::Display for FiniteNonNegativeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite => f.write_str("value must be finite"),
            Self::Negative => f.write_str("value must be nonnegative"),
        }
    }
}

impl std::error::Error for FiniteNonNegativeError {}

/// Finite nonnegative binary64 value with canonical signed-zero semantics.
///
/// `-0.0` is normalized to `+0.0`; all identity and equality operations use
/// the retained canonical bits rather than floating comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FiniteNonNegative(u64);

impl FiniteNonNegative {
    /// Validate and canonicalize a finite nonnegative scalar.
    ///
    /// # Errors
    /// Refuses NaN, infinity, and negative values.
    pub fn new(value: f64) -> Result<Self, FiniteNonNegativeError> {
        if !value.is_finite() {
            return Err(FiniteNonNegativeError::NonFinite);
        }
        if value < 0.0 {
            return Err(FiniteNonNegativeError::Negative);
        }
        let canonical = if value == 0.0 { 0.0 } else { value };
        Ok(Self(canonical.to_bits()))
    }

    /// Canonical binary64 value.
    #[must_use]
    pub fn get(self) -> f64 {
        f64::from_bits(self.0)
    }

    /// Canonical IEEE-754 bits.
    #[must_use]
    pub const fn bits(self) -> u64 {
        self.0
    }

    /// Whether the value is exactly zero.
    #[must_use]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }
}

/// Quantity, shape, clock, and frame contract of one state slot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StateSlotContract {
    /// Contracted state slot.
    pub id: StateSlotId,
    /// Owning subsystem; admission checks this against the graph.
    pub owner: SubsystemId,
    /// Quantity dimensions and optional semantic kind.
    pub quantity: TerminalQuantitySpec,
    /// Scalar/vector/tensor/trace value shape.
    pub shape: TerminalShape,
    /// Logical clock used by the state.
    pub clock: ClockId,
    /// Coordinate frame and orientation.
    pub frame: FrameBinding,
}

/// Unique condition target in a behavior overlay.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConditionTarget {
    /// Initial value of one owned state slot.
    Initial(StateSlotId),
    /// Source closure for one `ExternalInput` terminal.
    Boundary(TerminalId),
}

/// Declared continuity of a condition history.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HistoryContinuity {
    /// No discontinuity is declared in the history's admitted domain.
    Continuous,
    /// Jumps occur only at the listed Machine-IR events.
    ResetAtEvents {
        /// Nonempty event set; admission canonicalizes it.
        events: Vec<EventId>,
    },
}

/// Nominal source of an initial or boundary condition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConditionSource {
    /// One immutable fixed-value artifact.
    Fixed(ConditionValueRef),
    /// One time history with explicit continuity/reset semantics.
    History {
        /// Exact history artifact.
        history: ConditionHistoryRef,
        /// Continuity claim made by this graph.
        continuity: HistoryContinuity,
    },
    /// One declared distribution; this is not a validation claim.
    Distribution(DistributionRef),
}

/// Typed condition value bound to an initial or boundary target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConditionBinding {
    /// Unique state-slot or external-terminal role.
    pub target: ConditionTarget,
    /// Exact quantity declaration.
    pub quantity: TerminalQuantitySpec,
    /// Exact value shape.
    pub shape: TerminalShape,
    /// Exact logical clock.
    pub clock: ClockId,
    /// Exact frame/orientation convention.
    pub frame: FrameBinding,
    /// Fixed, historical, or distribution-valued source.
    pub source: ConditionSource,
}

/// Explicit motion behavior for one owned body.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BodyMotion {
    /// Identity path; static behavior is explicit rather than defaulted.
    Static,
    /// Prescribed path owned by an external motion domain.
    Prescribed {
        /// Exact versioned path/law reference.
        path: MotionPathRef,
    },
}

/// One body's explicit motion declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MotionBinding {
    /// Owned body occurrence.
    pub body: BodyId,
    /// Logical clock of the motion law.
    pub clock: ClockId,
    /// Reference frame and orientation for the path.
    pub reference_frame: FrameBinding,
    /// Static identity path or prescribed path.
    pub motion: BodyMotion,
}

/// State or terminal read by an event guard.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EventDependency {
    /// Typed state dependency.
    State(StateSlotId),
    /// Typed terminal dependency.
    Terminal(TerminalId),
}

/// Oriented crossing convention for a scalar guard artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GuardOrientation {
    /// Negative guard value to positive.
    NegativeToPositive,
    /// Positive guard value to negative.
    PositiveToNegative,
    /// Either direction may trigger; uniqueness is not implied.
    Bidirectional,
}

/// Honesty state for local guard crossing semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CrossingSemantics {
    /// The declaration names an artifact claimed to establish transversality.
    Transverse(EventWitnessRef),
    /// The declaration names an artifact claimed to establish grazing.
    Grazing(EventWitnessRef),
    /// No crossing classification is claimed.
    Unknown(NoClaimRef),
}

/// Reset relation following one event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResetSemantics {
    /// Deterministic map writing the listed state slots.
    Deterministic {
        /// Exact reset map.
        map: ResetMapRef,
        /// Nonempty write set; admission canonicalizes it.
        writes: Vec<StateSlotId>,
    },
    /// Exact set-valued relation writing the listed state slots.
    SetValued {
        /// Exact reset relation.
        relation: ResetMapRef,
        /// Exact outcome-set semantics.
        outcomes: OutcomeSetRef,
        /// Nonempty write set; admission canonicalizes it.
        writes: Vec<StateSlotId>,
    },
    /// Event terminates execution under the named restriction relation.
    Terminal {
        /// Exact terminal relation.
        relation: ResetMapRef,
    },
    /// Reset behavior is deliberately unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: NoClaimRef,
    },
}

/// Superdense-time interpretation of potentially simultaneous events.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EventOrder {
    /// Deterministic microstep priority; values must be unique per clock.
    TotalPriority {
        /// Lower values execute at earlier superdense microsteps.
        microstep: u32,
    },
    /// Simultaneous outcomes remain explicitly set-valued.
    SetValued {
        /// Exact simultaneity-group semantics shared by at least two events.
        group: SimultaneityGroupRef,
    },
}

/// One guarded event with reset and explicit simultaneity semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EventSpec {
    /// Durable event identity.
    pub id: EventId,
    /// Event-driven logical clock.
    pub clock: ClockId,
    /// Exact guard artifact.
    pub guard: GuardRef,
    /// Explicit crossing direction.
    pub orientation: GuardOrientation,
    /// Transverse, grazing, or honest unknown classification.
    pub crossing: CrossingSemantics,
    /// Nonempty state/terminal dependency set.
    pub dependencies: Vec<EventDependency>,
    /// Deterministic, set-valued, terminal, or unresolved reset relation.
    pub reset: ResetSemantics,
    /// Explicit superdense collision semantics.
    pub order: EventOrder,
}

/// Durable graph target of one parameter tolerance.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ToleranceTarget {
    /// Parameter owned by a subsystem model as a whole.
    Subsystem(SubsystemId),
    /// Parameter attached to a durable machine element.
    Element(MachineElementId),
}

/// Bounded or random tolerance semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ToleranceSemantics {
    /// Asymmetric guaranteed envelope around a nominal parameter.
    Bounded {
        /// Nonnegative downward width.
        minus: FiniteNonNegative,
        /// Nonnegative upward width.
        plus: FiniteNonNegative,
        /// Exact external law defining application of the envelope.
        law: ToleranceLawRef,
    },
    /// Random variation with a strictly positive scale.
    Random {
        /// Strictly positive quantity-scale factor.
        scale: FiniteNonNegative,
        /// Exact random-variation law.
        law: ToleranceLawRef,
        /// Exact marginal distribution.
        marginal: DistributionRef,
    },
}

/// One typed tolerance on a stable target parameter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ToleranceSpec {
    /// Durable tolerance identity.
    pub id: ToleranceId,
    /// Stable graph owner or element target.
    pub target: ToleranceTarget,
    /// Exact parameter artifact within the target.
    pub parameter: ParameterRef,
    /// Parameter quantity semantics.
    pub quantity: TerminalQuantitySpec,
    /// Parameter value shape.
    pub shape: TerminalShape,
    /// Bounded or random interpretation.
    pub semantics: ToleranceSemantics,
}

/// Explicit dependence interpretation for random tolerances.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DependenceModel {
    /// All named random axes are explicitly declared mutually independent.
    Independent,
    /// Dependence is delegated to an exact externally owned model artifact,
    /// interpreted over the declaration's canonical member order.
    Correlated(CorrelationModelRef),
}

/// One random source axis in the global dependence declaration.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DependenceMember {
    /// Distribution-valued initial or boundary condition.
    Condition(ConditionTarget),
    /// Random tolerance declaration.
    Tolerance(ToleranceId),
}

/// Complete global dependence closure for every random behavior source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DependenceSpec {
    /// Canonical axis set; every random condition/tolerance appears exactly once.
    pub members: Vec<DependenceMember>,
    /// Global mutual independence or one correlated joint-model semantics.
    pub model: DependenceModel,
}

/// Canonical identity schema for one admitted behavior overlay.
pub enum MachineBehaviorIdentitySchemaV1 {}

impl CanonicalSchema for MachineBehaviorIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.behavior.v1";
    const NAME: &'static str = "admitted-machine-behavior";
    const VERSION: u32 = MACHINE_BEHAVIOR_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "base graph, state contracts, conditions, body motions, events and resets, tolerances, and explicit dependence closure";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("behavior-schema-version", WireType::U64),
        FieldSpec::required("base-machine-graph", WireType::Bytes),
        FieldSpec::required("state-contracts", WireType::OrderedBytes),
        FieldSpec::required("conditions", WireType::OrderedBytes),
        FieldSpec::required("motions", WireType::OrderedBytes),
        FieldSpec::required("events", WireType::OrderedBytes),
        FieldSpec::required("tolerances", WireType::OrderedBytes),
        FieldSpec::required("dependences", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of an admitted behavior overlay.
pub type MachineBehaviorIdV1 = ProblemSemanticId<MachineBehaviorIdentitySchemaV1>;

/// Closed rule vocabulary for deterministic behavior diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum MachineBehaviorRule {
    /// A public collection or aggregate nested-reference limit was exceeded.
    ResourceLimit = 1,
    /// One state slot has multiple contracts.
    DuplicateStateContract = 2,
    /// A state contract or initial condition references an unknown state slot.
    UnknownStateSlot = 3,
    /// A state contract names the wrong owning subsystem.
    StateOwnerMismatch = 4,
    /// A state contract references an unknown clock.
    UnknownStateClock = 5,
    /// A state contract carries an inadmissible semantic scalar form.
    UnsupportedStateQuantity = 6,
    /// An owned state slot has no contract.
    MissingStateContract = 7,
    /// One initial/boundary target has multiple condition sources.
    DuplicateCondition = 8,
    /// A boundary condition references an unknown terminal.
    UnknownBoundaryTerminal = 9,
    /// A history was supplied where only an initial value/distribution is valid.
    InvalidInitialSource = 10,
    /// A condition carries an inadmissible semantic scalar form.
    UnsupportedConditionQuantity = 11,
    /// A condition's quantity differs from its target contract.
    ConditionQuantityGap = 12,
    /// A condition's shape differs from its target contract.
    ConditionShapeGap = 13,
    /// A condition's clock differs from its target contract.
    ConditionClockGap = 14,
    /// A condition's frame/orientation differs from its target contract.
    ConditionFrameGap = 15,
    /// A condition references a clock absent from the graph.
    UnknownConditionClock = 16,
    /// A boundary source targets something other than `ExternalInput`.
    BoundaryCausalityGap = 17,
    /// An owned state slot has no initial condition.
    MissingInitialCondition = 18,
    /// An external-input terminal has no boundary condition.
    MissingBoundaryCondition = 19,
    /// A history-reset marker references an unknown event.
    UnknownHistoryEvent = 20,
    /// A history-reset event is repeated.
    DuplicateHistoryEvent = 21,
    /// A reset-delimited history supplied an empty event set.
    EmptyHistoryEvents = 22,
    /// One body has multiple motion declarations.
    DuplicateMotion = 23,
    /// A motion references a body absent from the graph.
    UnknownMotionBody = 24,
    /// A motion references a clock absent from the graph.
    UnknownMotionClock = 25,
    /// An owned body has no explicit static/prescribed motion.
    MissingMotion = 26,
    /// One event identity is declared more than once.
    DuplicateEvent = 27,
    /// An event references a clock absent from the graph.
    UnknownEventClock = 28,
    /// An event uses a clock not declared `EventDriven`.
    NonEventClock = 29,
    /// An event guard declares no state/terminal dependencies.
    EmptyEventDependencies = 30,
    /// An event guard dependency is repeated.
    DuplicateEventDependency = 31,
    /// An event guard dependency is absent from the graph.
    UnknownEventDependency = 32,
    /// Two priority-ordered events on one clock share a microstep.
    DuplicateEventPriority = 33,
    /// A set-valued simultaneity group contains fewer than two events.
    SingletonSimultaneityGroup = 34,
    /// Members of one simultaneity group name different clocks.
    SimultaneityClockGap = 35,
    /// A deterministic/set-valued reset declares no writes.
    EmptyResetWrites = 36,
    /// A reset write target is repeated.
    DuplicateResetWrite = 37,
    /// A reset writes a state slot absent from the graph.
    UnknownResetState = 38,
    /// One tolerance identity is declared more than once.
    DuplicateTolerance = 39,
    /// A tolerance target is absent from the graph.
    UnknownToleranceTarget = 40,
    /// A tolerance carries an inadmissible semantic scalar form.
    UnsupportedToleranceQuantity = 41,
    /// A state/terminal tolerance has a quantity mismatch.
    ToleranceQuantityGap = 42,
    /// A state/terminal tolerance has a shape mismatch.
    ToleranceShapeGap = 43,
    /// A bounded envelope is zero on both sides or a random scale is zero.
    ZeroTolerance = 44,
    /// A dependence declaration has no members.
    EmptyDependence = 45,
    /// A member is repeated within one dependence declaration.
    DuplicateDependenceMember = 46,
    /// A dependence member names an unknown tolerance.
    UnknownToleranceMember = 47,
    /// A bounded-only tolerance was placed in statistical dependence.
    BoundedToleranceInDependence = 48,
    /// A random condition/tolerance appears in more than one declaration.
    DuplicateDependenceCoverage = 49,
    /// A random condition/tolerance has no explicit dependence declaration.
    MissingDependence = 50,
    /// A correlated model has fewer than two distinct members.
    CorrelatedGroupTooSmall = 51,
    /// More than one global joint-dependence declaration was submitted.
    MultipleDependenceModels = 52,
    /// A dependence member names an unknown condition target.
    UnknownConditionMember = 53,
    /// A fixed/history condition was placed in statistical dependence.
    NonRandomConditionInDependence = 54,
    /// One event clock mixes total-priority and set-valued collision policies.
    MixedEventOrderPolicy = 55,
    /// One event clock names more than one set-valued simultaneity group.
    MultipleSimultaneityGroups = 56,
    /// Distinct tolerance IDs bind the same exact target parameter.
    DuplicateToleranceBinding = 57,
    /// Bounded canonical identity publication failed.
    Identity = 58,
}

impl MachineBehaviorRule {
    /// Stable machine-readable rule code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ResourceLimit => "MachineBehaviorResourceLimit",
            Self::DuplicateStateContract => "MachineBehaviorDuplicateStateContract",
            Self::UnknownStateSlot => "MachineBehaviorUnknownStateSlot",
            Self::StateOwnerMismatch => "MachineBehaviorStateOwnerMismatch",
            Self::UnknownStateClock => "MachineBehaviorUnknownStateClock",
            Self::UnsupportedStateQuantity => "MachineBehaviorUnsupportedStateQuantity",
            Self::MissingStateContract => "MachineBehaviorMissingStateContract",
            Self::DuplicateCondition => "MachineBehaviorDuplicateCondition",
            Self::UnknownBoundaryTerminal => "MachineBehaviorUnknownBoundaryTerminal",
            Self::InvalidInitialSource => "MachineBehaviorInvalidInitialSource",
            Self::UnsupportedConditionQuantity => "MachineBehaviorUnsupportedConditionQuantity",
            Self::ConditionQuantityGap => "MachineBehaviorConditionQuantityGap",
            Self::ConditionShapeGap => "MachineBehaviorConditionShapeGap",
            Self::ConditionClockGap => "MachineBehaviorConditionClockGap",
            Self::ConditionFrameGap => "MachineBehaviorConditionFrameGap",
            Self::UnknownConditionClock => "MachineBehaviorUnknownConditionClock",
            Self::BoundaryCausalityGap => "MachineBehaviorBoundaryCausalityGap",
            Self::MissingInitialCondition => "MachineBehaviorMissingInitialCondition",
            Self::MissingBoundaryCondition => "MachineBehaviorMissingBoundaryCondition",
            Self::UnknownHistoryEvent => "MachineBehaviorUnknownHistoryEvent",
            Self::DuplicateHistoryEvent => "MachineBehaviorDuplicateHistoryEvent",
            Self::EmptyHistoryEvents => "MachineBehaviorEmptyHistoryEvents",
            Self::DuplicateMotion => "MachineBehaviorDuplicateMotion",
            Self::UnknownMotionBody => "MachineBehaviorUnknownMotionBody",
            Self::UnknownMotionClock => "MachineBehaviorUnknownMotionClock",
            Self::MissingMotion => "MachineBehaviorMissingMotion",
            Self::DuplicateEvent => "MachineBehaviorDuplicateEvent",
            Self::UnknownEventClock => "MachineBehaviorUnknownEventClock",
            Self::NonEventClock => "MachineBehaviorNonEventClock",
            Self::EmptyEventDependencies => "MachineBehaviorEmptyEventDependencies",
            Self::DuplicateEventDependency => "MachineBehaviorDuplicateEventDependency",
            Self::UnknownEventDependency => "MachineBehaviorUnknownEventDependency",
            Self::DuplicateEventPriority => "MachineBehaviorDuplicateEventPriority",
            Self::SingletonSimultaneityGroup => "MachineBehaviorSingletonSimultaneityGroup",
            Self::SimultaneityClockGap => "MachineBehaviorSimultaneityClockGap",
            Self::EmptyResetWrites => "MachineBehaviorEmptyResetWrites",
            Self::DuplicateResetWrite => "MachineBehaviorDuplicateResetWrite",
            Self::UnknownResetState => "MachineBehaviorUnknownResetState",
            Self::DuplicateTolerance => "MachineBehaviorDuplicateTolerance",
            Self::UnknownToleranceTarget => "MachineBehaviorUnknownToleranceTarget",
            Self::UnsupportedToleranceQuantity => "MachineBehaviorUnsupportedToleranceQuantity",
            Self::ToleranceQuantityGap => "MachineBehaviorToleranceQuantityGap",
            Self::ToleranceShapeGap => "MachineBehaviorToleranceShapeGap",
            Self::ZeroTolerance => "MachineBehaviorZeroTolerance",
            Self::EmptyDependence => "MachineBehaviorEmptyDependence",
            Self::DuplicateDependenceMember => "MachineBehaviorDuplicateDependenceMember",
            Self::UnknownToleranceMember => "MachineBehaviorUnknownToleranceMember",
            Self::BoundedToleranceInDependence => "MachineBehaviorBoundedToleranceInDependence",
            Self::DuplicateDependenceCoverage => "MachineBehaviorDuplicateDependenceCoverage",
            Self::MissingDependence => "MachineBehaviorMissingDependence",
            Self::CorrelatedGroupTooSmall => "MachineBehaviorCorrelatedGroupTooSmall",
            Self::MultipleDependenceModels => "MachineBehaviorMultipleDependenceModels",
            Self::UnknownConditionMember => "MachineBehaviorUnknownConditionMember",
            Self::NonRandomConditionInDependence => "MachineBehaviorNonRandomConditionInDependence",
            Self::MixedEventOrderPolicy => "MachineBehaviorMixedEventOrderPolicy",
            Self::MultipleSimultaneityGroups => "MachineBehaviorMultipleSimultaneityGroups",
            Self::DuplicateToleranceBinding => "MachineBehaviorDuplicateToleranceBinding",
            Self::Identity => "MachineBehaviorIdentity",
        }
    }
}

/// Stable diagnostic subject for behavior admission.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum MachineBehaviorSubject {
    /// Complete overlay.
    Behavior,
    /// State contract or state-closure obligation.
    State(StateSlotId),
    /// Initial or boundary condition target.
    Condition(ConditionTarget),
    /// Body-motion target.
    Motion(BodyId),
    /// Event declaration.
    Event(EventId),
    /// Clock related to a finding.
    Clock(ClockId),
    /// Event dependency or tolerance element target.
    Element(MachineElementId),
    /// Subsystem-level tolerance target.
    Subsystem(SubsystemId),
    /// Tolerance declaration.
    Tolerance(ToleranceId),
    /// Global dependence declaration identified by its canonical axis set.
    Dependence(Vec<DependenceMember>),
}

/// One deterministic behavior-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct MachineBehaviorFinding {
    rule: MachineBehaviorRule,
    subject: MachineBehaviorSubject,
    related: Option<MachineBehaviorSubject>,
}

impl MachineBehaviorFinding {
    fn new(
        rule: MachineBehaviorRule,
        subject: MachineBehaviorSubject,
        related: Option<MachineBehaviorSubject>,
    ) -> Self {
        Self {
            rule,
            subject,
            related,
        }
    }

    /// Closed rule category.
    #[must_use]
    pub const fn rule(&self) -> MachineBehaviorRule {
        self.rule
    }

    /// Stable rule code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }

    /// Primary offending subject.
    #[must_use]
    pub const fn subject(&self) -> &MachineBehaviorSubject {
        &self.subject
    }

    /// Optional related subject.
    #[must_use]
    pub const fn related(&self) -> Option<&MachineBehaviorSubject> {
        self.related.as_ref()
    }
}

/// Complete deterministic refusal from behavior admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineBehaviorRefusal {
    findings: Vec<MachineBehaviorFinding>,
    identity_error: Option<CanonicalError>,
}

impl MachineBehaviorRefusal {
    /// Sorted, duplicate-free findings.
    #[must_use]
    pub fn findings(&self) -> &[MachineBehaviorFinding] {
        &self.findings
    }

    /// Canonical identity error, only for `MachineBehaviorIdentity`.
    #[must_use]
    pub const fn identity_error(&self) -> Option<&CanonicalError> {
        self.identity_error.as_ref()
    }

    /// Stable top-level refusal code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        "MachineBehaviorRefused"
    }
}

impl fmt::Display for MachineBehaviorRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "machine behavior refused with {} finding(s)",
            self.findings.len()
        )
    }
}

impl std::error::Error for MachineBehaviorRefusal {}

/// Collection sizes retained for every behavior-admission attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineBehaviorSubmittedCounts {
    /// State contracts submitted.
    pub state_contracts: usize,
    /// Initial plus boundary conditions submitted.
    pub conditions: usize,
    /// Motion bindings submitted.
    pub motions: usize,
    /// Events submitted.
    pub events: usize,
    /// Tolerances submitted.
    pub tolerances: usize,
    /// Dependence declarations submitted.
    pub dependences: usize,
    /// Aggregate history events, dependencies, reset writes, and group members.
    pub nested_references: usize,
}

/// Mutable-by-construction behavior overlay with no authority before admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineBehaviorDraft {
    /// State-slot type contracts.
    pub state_contracts: Vec<StateSlotContract>,
    /// Initial and boundary condition sources.
    pub conditions: Vec<ConditionBinding>,
    /// Explicit body-motion declarations.
    pub motions: Vec<MotionBinding>,
    /// Guard/reset event declarations.
    pub events: Vec<EventSpec>,
    /// Typed tolerance declarations.
    pub tolerances: Vec<ToleranceSpec>,
    /// Zero or one global dependence closure for all random sources.
    pub dependences: Vec<DependenceSpec>,
}

impl MachineBehaviorDraft {
    /// Admit this overlay against one exact admitted machine graph.
    ///
    /// # Errors
    /// Returns every deterministic finding discovered within public bounds.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineBehavior, MachineBehaviorRefusal> {
        self.admit_with_decision(graph).into_result()
    }

    /// Attempt admission while retaining submitted collection counts.
    #[must_use]
    pub fn admit_with_decision(
        self,
        graph: &AdmittedMachineGraph,
    ) -> MachineBehaviorAdmissionDecision {
        let submitted = submitted_counts(&self);
        MachineBehaviorAdmissionDecision {
            submitted,
            result: admit_machine_behavior(self, graph),
        }
    }
}

/// Canonically ordered admitted behavior plus its semantic receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineBehavior {
    base_graph: MachineGraphIdV1,
    state_contracts: Vec<StateSlotContract>,
    conditions: Vec<ConditionBinding>,
    motions: Vec<MotionBinding>,
    events: Vec<EventSpec>,
    tolerances: Vec<ToleranceSpec>,
    dependences: Vec<DependenceSpec>,
    receipt: IdentityReceipt<MachineBehaviorIdV1>,
}

impl AdmittedMachineBehavior {
    /// Exact graph identity this overlay extends.
    #[must_use]
    pub const fn base_graph(&self) -> MachineGraphIdV1 {
        self.base_graph
    }

    /// Strong semantic identity of the admitted overlay.
    #[must_use]
    pub const fn identity(&self) -> MachineBehaviorIdV1 {
        self.receipt.id()
    }

    /// Complete canonical identity receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineBehaviorIdV1> {
        self.receipt
    }

    /// Canonically ordered state contracts.
    #[must_use]
    pub fn state_contracts(&self) -> &[StateSlotContract] {
        &self.state_contracts
    }

    /// Canonically ordered conditions.
    #[must_use]
    pub fn conditions(&self) -> &[ConditionBinding] {
        &self.conditions
    }

    /// Canonically ordered body motions.
    #[must_use]
    pub fn motions(&self) -> &[MotionBinding] {
        &self.motions
    }

    /// Canonically ordered events.
    #[must_use]
    pub fn events(&self) -> &[EventSpec] {
        &self.events
    }

    /// Canonically ordered tolerances.
    #[must_use]
    pub fn tolerances(&self) -> &[ToleranceSpec] {
        &self.tolerances
    }

    /// Canonical global dependence declaration (empty or singleton after admission).
    #[must_use]
    pub fn dependences(&self) -> &[DependenceSpec] {
        &self.dependences
    }
}

/// Bounded deterministic outcome summary for one behavior admission.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineBehaviorAdmissionDecision {
    submitted: MachineBehaviorSubmittedCounts,
    result: Result<AdmittedMachineBehavior, MachineBehaviorRefusal>,
}

impl MachineBehaviorAdmissionDecision {
    /// Exact collection sizes observed before canonicalization.
    #[must_use]
    pub const fn submitted_counts(&self) -> MachineBehaviorSubmittedCounts {
        self.submitted
    }

    /// Stable top-level decision code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "MachineBehaviorAdmitted",
            Err(_) => "MachineBehaviorRefused",
        }
    }

    /// Borrow the admitted overlay or complete refusal.
    #[must_use]
    pub fn result(&self) -> Result<&AdmittedMachineBehavior, &MachineBehaviorRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision and recover the conventional result.
    #[must_use]
    pub fn into_result(self) -> Result<AdmittedMachineBehavior, MachineBehaviorRefusal> {
        self.result
    }
}

fn submitted_counts(draft: &MachineBehaviorDraft) -> MachineBehaviorSubmittedCounts {
    let nested_references = draft
        .conditions
        .iter()
        .fold(0usize, |count, condition| {
            count.saturating_add(match &condition.source {
                ConditionSource::History {
                    continuity: HistoryContinuity::ResetAtEvents { events },
                    ..
                } => events.len(),
                ConditionSource::Fixed(_)
                | ConditionSource::Distribution(_)
                | ConditionSource::History {
                    continuity: HistoryContinuity::Continuous,
                    ..
                } => 0,
            })
        })
        .saturating_add(draft.events.iter().fold(0usize, |count, event| {
            count
                .saturating_add(event.dependencies.len())
                .saturating_add(reset_writes(&event.reset).map_or(0, |writes| writes.len()))
        }))
        .saturating_add(draft.dependences.iter().fold(0usize, |count, group| {
            count.saturating_add(group.members.len())
        }));
    MachineBehaviorSubmittedCounts {
        state_contracts: draft.state_contracts.len(),
        conditions: draft.conditions.len(),
        motions: draft.motions.len(),
        events: draft.events.len(),
        tolerances: draft.tolerances.len(),
        dependences: draft.dependences.len(),
        nested_references,
    }
}

fn resource_limit_findings(counts: MachineBehaviorSubmittedCounts) -> Vec<MachineBehaviorFinding> {
    let over_limit = counts.state_contracts > MAX_MACHINE_BEHAVIOR_STATE_CONTRACTS
        || counts.conditions > MAX_MACHINE_BEHAVIOR_CONDITIONS
        || counts.motions > MAX_MACHINE_BEHAVIOR_MOTIONS
        || counts.events > MAX_MACHINE_BEHAVIOR_EVENTS
        || counts.tolerances > MAX_MACHINE_BEHAVIOR_TOLERANCES
        || counts.dependences > MAX_MACHINE_BEHAVIOR_DEPENDENCES
        || counts.nested_references > MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES;
    if over_limit {
        vec![MachineBehaviorFinding::new(
            MachineBehaviorRule::ResourceLimit,
            MachineBehaviorSubject::Behavior,
            None,
        )]
    } else {
        Vec::new()
    }
}

fn behavior_refusal(
    mut findings: Vec<MachineBehaviorFinding>,
    identity_error: Option<CanonicalError>,
) -> MachineBehaviorRefusal {
    findings.sort();
    findings.dedup();
    debug_assert!(!findings.is_empty());
    MachineBehaviorRefusal {
        findings,
        identity_error,
    }
}

struct GraphIndex<'a> {
    clocks: BTreeMap<ClockId, MachineClock>,
    subsystems: BTreeSet<SubsystemId>,
    state_owners: BTreeMap<StateSlotId, SubsystemId>,
    bodies: BTreeSet<BodyId>,
    surface_patches: BTreeSet<super::SurfacePatchId>,
    contact_features: BTreeSet<super::ContactFeatureId>,
    terminals: BTreeMap<TerminalId, &'a super::TerminalSpec>,
    ports: BTreeSet<super::PortId>,
}

impl<'a> GraphIndex<'a> {
    fn new(graph: &'a AdmittedMachineGraph) -> Self {
        let clocks = graph
            .clocks()
            .iter()
            .map(|clock| (clock.id.clone(), clock.clock))
            .collect();
        let subsystems = graph
            .subsystems()
            .iter()
            .map(|subsystem| subsystem.id.clone())
            .collect();
        let mut state_owners = BTreeMap::new();
        let mut bodies = BTreeSet::new();
        let mut surface_patches = BTreeSet::new();
        let mut contact_features = BTreeSet::new();
        for subsystem in graph.subsystems() {
            for state in &subsystem.state_slots {
                state_owners.insert(state.clone(), subsystem.id.clone());
            }
            bodies.extend(subsystem.bodies.iter().cloned());
            surface_patches.extend(subsystem.surface_patches.iter().cloned());
            contact_features.extend(subsystem.contact_features.iter().cloned());
        }
        let terminals = graph
            .terminals()
            .iter()
            .map(|terminal| (terminal.id.clone(), terminal))
            .collect();
        let ports = graph.ports().iter().map(|port| port.id.clone()).collect();
        Self {
            clocks,
            subsystems,
            state_owners,
            bodies,
            surface_patches,
            contact_features,
            terminals,
            ports,
        }
    }

    fn element_exists(&self, element: &MachineElementId) -> bool {
        match element {
            MachineElementId::Body(id) => self.bodies.contains(id),
            MachineElementId::SurfacePatch(id) => self.surface_patches.contains(id),
            MachineElementId::ContactFeature(id) => self.contact_features.contains(id),
            MachineElementId::Terminal(id) => self.terminals.contains_key(id),
            MachineElementId::Port(id) => self.ports.contains(id),
            MachineElementId::StateSlot(id) => self.state_owners.contains_key(id),
        }
    }

    fn dependency_exists(&self, dependency: &EventDependency) -> bool {
        match dependency {
            EventDependency::State(id) => self.state_owners.contains_key(id),
            EventDependency::Terminal(id) => self.terminals.contains_key(id),
        }
    }
}

#[allow(clippy::too_many_lines)]
fn admit_machine_behavior(
    mut draft: MachineBehaviorDraft,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineBehavior, MachineBehaviorRefusal> {
    let counts = submitted_counts(&draft);
    let resource_findings = resource_limit_findings(counts);
    if !resource_findings.is_empty() {
        return Err(behavior_refusal(resource_findings, None));
    }

    canonicalize_behavior_draft(&mut draft);
    let index = GraphIndex::new(graph);
    let mut findings = Vec::new();

    for pair in draft.state_contracts.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateStateContract,
                MachineBehaviorSubject::State(pair[1].id.clone()),
                None,
            ));
        }
    }
    let mut state_contracts = BTreeMap::<StateSlotId, &StateSlotContract>::new();
    for contract in &draft.state_contracts {
        match index.state_owners.get(&contract.id) {
            None => findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownStateSlot,
                MachineBehaviorSubject::State(contract.id.clone()),
                None,
            )),
            Some(owner) if *owner != contract.owner => findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::StateOwnerMismatch,
                MachineBehaviorSubject::State(contract.id.clone()),
                Some(MachineBehaviorSubject::Subsystem(contract.owner.clone())),
            )),
            Some(_) => {}
        }
        if !index.clocks.contains_key(&contract.clock) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownStateClock,
                MachineBehaviorSubject::State(contract.id.clone()),
                Some(MachineBehaviorSubject::Clock(contract.clock.clone())),
            ));
        }
        if !contract.quantity.is_admitted() {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnsupportedStateQuantity,
                MachineBehaviorSubject::State(contract.id.clone()),
                None,
            ));
        }
        state_contracts
            .entry(contract.id.clone())
            .or_insert(contract);
    }
    for state in index.state_owners.keys() {
        if !state_contracts.contains_key(state) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingStateContract,
                MachineBehaviorSubject::State(state.clone()),
                None,
            ));
        }
    }

    for pair in draft.conditions.windows(2) {
        if pair[0].target == pair[1].target {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateCondition,
                MachineBehaviorSubject::Condition(pair[1].target.clone()),
                None,
            ));
        }
    }
    let event_ids: BTreeSet<EventId> = draft.events.iter().map(|event| event.id.clone()).collect();
    let mut initial_targets = BTreeSet::new();
    let mut boundary_targets = BTreeSet::new();
    for condition in &draft.conditions {
        let subject = || MachineBehaviorSubject::Condition(condition.target.clone());
        if !condition.quantity.is_admitted() {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnsupportedConditionQuantity,
                subject(),
                None,
            ));
        }
        if !index.clocks.contains_key(&condition.clock) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownConditionClock,
                subject(),
                Some(MachineBehaviorSubject::Clock(condition.clock.clone())),
            ));
        }
        match &condition.target {
            ConditionTarget::Initial(state) => {
                initial_targets.insert(state.clone());
                if !index.state_owners.contains_key(state) {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::UnknownStateSlot,
                        subject(),
                        Some(MachineBehaviorSubject::State(state.clone())),
                    ));
                }
                if matches!(&condition.source, ConditionSource::History { .. }) {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::InvalidInitialSource,
                        subject(),
                        None,
                    ));
                }
                if let Some(contract) = state_contracts.get(state) {
                    compare_condition_contract(
                        condition,
                        contract.quantity,
                        contract.shape,
                        &contract.clock,
                        &contract.frame,
                        &mut findings,
                    );
                }
            }
            ConditionTarget::Boundary(terminal_id) => {
                boundary_targets.insert(terminal_id.clone());
                match index.terminals.get(terminal_id) {
                    None => findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::UnknownBoundaryTerminal,
                        subject(),
                        Some(MachineBehaviorSubject::Element(terminal_id.clone().into())),
                    )),
                    Some(terminal) => {
                        if terminal.causality != TerminalCausality::ExternalInput {
                            findings.push(MachineBehaviorFinding::new(
                                MachineBehaviorRule::BoundaryCausalityGap,
                                subject(),
                                Some(MachineBehaviorSubject::Element(terminal_id.clone().into())),
                            ));
                        }
                        compare_condition_contract(
                            condition,
                            terminal.quantity,
                            terminal.shape,
                            &terminal.clock,
                            &terminal.frame,
                            &mut findings,
                        );
                    }
                }
            }
        }
        validate_history_events(condition, &event_ids, &mut findings);
    }
    for state in index.state_owners.keys() {
        if !initial_targets.contains(state) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingInitialCondition,
                MachineBehaviorSubject::Condition(ConditionTarget::Initial(state.clone())),
                None,
            ));
        }
    }
    for terminal in index.terminals.values() {
        if terminal.causality == TerminalCausality::ExternalInput
            && !boundary_targets.contains(&terminal.id)
        {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingBoundaryCondition,
                MachineBehaviorSubject::Condition(ConditionTarget::Boundary(terminal.id.clone())),
                None,
            ));
        }
    }

    for pair in draft.motions.windows(2) {
        if pair[0].body == pair[1].body {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateMotion,
                MachineBehaviorSubject::Motion(pair[1].body.clone()),
                None,
            ));
        }
    }
    let mut motion_targets = BTreeSet::new();
    for motion in &draft.motions {
        motion_targets.insert(motion.body.clone());
        if !index.bodies.contains(&motion.body) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownMotionBody,
                MachineBehaviorSubject::Motion(motion.body.clone()),
                None,
            ));
        }
        if !index.clocks.contains_key(&motion.clock) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownMotionClock,
                MachineBehaviorSubject::Motion(motion.body.clone()),
                Some(MachineBehaviorSubject::Clock(motion.clock.clone())),
            ));
        }
    }
    for body in &index.bodies {
        if !motion_targets.contains(body) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingMotion,
                MachineBehaviorSubject::Motion(body.clone()),
                None,
            ));
        }
    }

    validate_events(&draft.events, &index, &mut findings);
    validate_tolerances_and_dependences(
        &draft.conditions,
        &draft.tolerances,
        &draft.dependences,
        &index,
        &state_contracts,
        &mut findings,
    );

    if !findings.is_empty() {
        return Err(behavior_refusal(findings, None));
    }

    let receipt = match machine_behavior_identity(&draft, graph.identity()) {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(behavior_refusal(
                vec![MachineBehaviorFinding::new(
                    MachineBehaviorRule::Identity,
                    MachineBehaviorSubject::Behavior,
                    None,
                )],
                Some(error),
            ));
        }
    };

    Ok(AdmittedMachineBehavior {
        base_graph: graph.identity(),
        state_contracts: draft.state_contracts,
        conditions: draft.conditions,
        motions: draft.motions,
        events: draft.events,
        tolerances: draft.tolerances,
        dependences: draft.dependences,
        receipt,
    })
}

fn compare_condition_contract(
    condition: &ConditionBinding,
    quantity: TerminalQuantitySpec,
    shape: TerminalShape,
    clock: &ClockId,
    frame: &FrameBinding,
    findings: &mut Vec<MachineBehaviorFinding>,
) {
    let subject = || MachineBehaviorSubject::Condition(condition.target.clone());
    if condition.quantity != quantity {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ConditionQuantityGap,
            subject(),
            None,
        ));
    }
    if condition.shape != shape {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ConditionShapeGap,
            subject(),
            None,
        ));
    }
    if condition.clock != *clock {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ConditionClockGap,
            subject(),
            Some(MachineBehaviorSubject::Clock(condition.clock.clone())),
        ));
    }
    if condition.frame != *frame {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ConditionFrameGap,
            subject(),
            None,
        ));
    }
}

fn validate_history_events(
    condition: &ConditionBinding,
    event_ids: &BTreeSet<EventId>,
    findings: &mut Vec<MachineBehaviorFinding>,
) {
    let ConditionSource::History {
        continuity: HistoryContinuity::ResetAtEvents { events },
        ..
    } = &condition.source
    else {
        return;
    };
    let subject = || MachineBehaviorSubject::Condition(condition.target.clone());
    if events.is_empty() {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::EmptyHistoryEvents,
            subject(),
            None,
        ));
    }
    for pair in events.windows(2) {
        if pair[0] == pair[1] {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateHistoryEvent,
                subject(),
                Some(MachineBehaviorSubject::Event(pair[1].clone())),
            ));
        }
    }
    for event in events {
        if !event_ids.contains(event) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownHistoryEvent,
                subject(),
                Some(MachineBehaviorSubject::Event(event.clone())),
            ));
        }
    }
}

fn validate_events(
    events: &[EventSpec],
    index: &GraphIndex<'_>,
    findings: &mut Vec<MachineBehaviorFinding>,
) {
    for pair in events.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateEvent,
                MachineBehaviorSubject::Event(pair[1].id.clone()),
                None,
            ));
        }
    }

    let mut priorities = BTreeMap::<(ClockId, u32), EventId>::new();
    let mut groups = BTreeMap::<SimultaneityGroupRef, Vec<(EventId, ClockId)>>::new();
    let mut priority_clocks = BTreeSet::<ClockId>::new();
    let mut groups_by_clock = BTreeMap::<ClockId, BTreeSet<SimultaneityGroupRef>>::new();
    for event in events {
        let subject = || MachineBehaviorSubject::Event(event.id.clone());
        match index.clocks.get(&event.clock) {
            None => findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownEventClock,
                subject(),
                Some(MachineBehaviorSubject::Clock(event.clock.clone())),
            )),
            Some(clock) if *clock != MachineClock::EventDriven => {
                findings.push(MachineBehaviorFinding::new(
                    MachineBehaviorRule::NonEventClock,
                    subject(),
                    Some(MachineBehaviorSubject::Clock(event.clock.clone())),
                ))
            }
            Some(_) => {}
        }

        if event.dependencies.is_empty() {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::EmptyEventDependencies,
                subject(),
                None,
            ));
        }
        for pair in event.dependencies.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineBehaviorFinding::new(
                    MachineBehaviorRule::DuplicateEventDependency,
                    subject(),
                    Some(dependency_subject(&pair[1])),
                ));
            }
        }
        for dependency in &event.dependencies {
            if !index.dependency_exists(dependency) {
                findings.push(MachineBehaviorFinding::new(
                    MachineBehaviorRule::UnknownEventDependency,
                    subject(),
                    Some(dependency_subject(dependency)),
                ));
            }
        }

        match &event.order {
            EventOrder::TotalPriority { microstep } => {
                priority_clocks.insert(event.clock.clone());
                let key = (event.clock.clone(), *microstep);
                if let Some(first) = priorities.insert(key, event.id.clone()) {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::DuplicateEventPriority,
                        subject(),
                        Some(MachineBehaviorSubject::Event(first)),
                    ));
                }
            }
            EventOrder::SetValued { group } => {
                groups_by_clock
                    .entry(event.clock.clone())
                    .or_default()
                    .insert(group.clone());
                groups
                    .entry(group.clone())
                    .or_default()
                    .push((event.id.clone(), event.clock.clone()));
            }
        }

        if let Some(writes) = reset_writes(&event.reset) {
            if writes.is_empty() {
                findings.push(MachineBehaviorFinding::new(
                    MachineBehaviorRule::EmptyResetWrites,
                    subject(),
                    None,
                ));
            }
            for pair in writes.windows(2) {
                if pair[0] == pair[1] {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::DuplicateResetWrite,
                        subject(),
                        Some(MachineBehaviorSubject::State(pair[1].clone())),
                    ));
                }
            }
            for state in writes {
                if !index.state_owners.contains_key(state) {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::UnknownResetState,
                        subject(),
                        Some(MachineBehaviorSubject::State(state.clone())),
                    ));
                }
            }
        }
    }

    for clock in &priority_clocks {
        if groups_by_clock.contains_key(clock) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MixedEventOrderPolicy,
                MachineBehaviorSubject::Clock(clock.clone()),
                None,
            ));
        }
    }
    for (clock, clock_groups) in &groups_by_clock {
        if clock_groups.len() > 1 {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MultipleSimultaneityGroups,
                MachineBehaviorSubject::Clock(clock.clone()),
                None,
            ));
        }
    }

    for members in groups.values() {
        let unique_events: BTreeSet<&EventId> = members.iter().map(|(event, _)| event).collect();
        if unique_events.len() < 2 {
            let event = members.first().map_or_else(
                || unreachable!("group map entries are nonempty"),
                |value| &value.0,
            );
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::SingletonSimultaneityGroup,
                MachineBehaviorSubject::Event(event.clone()),
                None,
            ));
        }
        let clocks: BTreeSet<&ClockId> = members.iter().map(|(_, clock)| clock).collect();
        if clocks.len() > 1 {
            let event = &members[0].0;
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::SimultaneityClockGap,
                MachineBehaviorSubject::Event(event.clone()),
                None,
            ));
        }
    }
}

fn dependency_subject(dependency: &EventDependency) -> MachineBehaviorSubject {
    match dependency {
        EventDependency::State(id) => MachineBehaviorSubject::Element(id.clone().into()),
        EventDependency::Terminal(id) => MachineBehaviorSubject::Element(id.clone().into()),
    }
}

fn reset_writes(reset: &ResetSemantics) -> Option<&[StateSlotId]> {
    match reset {
        ResetSemantics::Deterministic { writes, .. } | ResetSemantics::SetValued { writes, .. } => {
            Some(writes)
        }
        ResetSemantics::Terminal { .. } | ResetSemantics::Unknown { .. } => None,
    }
}

fn reset_writes_mut(reset: &mut ResetSemantics) -> Option<&mut Vec<StateSlotId>> {
    match reset {
        ResetSemantics::Deterministic { writes, .. } | ResetSemantics::SetValued { writes, .. } => {
            Some(writes)
        }
        ResetSemantics::Terminal { .. } | ResetSemantics::Unknown { .. } => None,
    }
}

fn validate_tolerances_and_dependences(
    conditions: &[ConditionBinding],
    tolerances: &[ToleranceSpec],
    dependences: &[DependenceSpec],
    index: &GraphIndex<'_>,
    state_contracts: &BTreeMap<StateSlotId, &StateSlotContract>,
    findings: &mut Vec<MachineBehaviorFinding>,
) {
    for pair in tolerances.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateTolerance,
                MachineBehaviorSubject::Tolerance(pair[1].id.clone()),
                None,
            ));
        }
    }

    let mut tolerance_random = BTreeMap::<ToleranceId, bool>::new();
    let mut tolerance_bindings = BTreeMap::<(ToleranceTarget, ParameterRef), ToleranceId>::new();
    for tolerance in tolerances {
        let subject = || MachineBehaviorSubject::Tolerance(tolerance.id.clone());
        let binding = (tolerance.target.clone(), tolerance.parameter.clone());
        if let Some(first) = tolerance_bindings.insert(binding, tolerance.id.clone()) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::DuplicateToleranceBinding,
                subject(),
                Some(MachineBehaviorSubject::Tolerance(first)),
            ));
        }
        let target_known = match &tolerance.target {
            ToleranceTarget::Subsystem(subsystem) => index.subsystems.contains(subsystem),
            ToleranceTarget::Element(element) => index.element_exists(element),
        };
        if !target_known {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnknownToleranceTarget,
                subject(),
                Some(tolerance_target_subject(&tolerance.target)),
            ));
        }
        if !tolerance.quantity.is_admitted() {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::UnsupportedToleranceQuantity,
                subject(),
                None,
            ));
        }

        match &tolerance.target {
            ToleranceTarget::Element(MachineElementId::Terminal(id)) => {
                if let Some(terminal) = index.terminals.get(id) {
                    compare_tolerance_contract(
                        tolerance,
                        terminal.quantity,
                        terminal.shape,
                        findings,
                    );
                }
            }
            ToleranceTarget::Element(MachineElementId::StateSlot(id)) => {
                if let Some(contract) = state_contracts.get(id) {
                    compare_tolerance_contract(
                        tolerance,
                        contract.quantity,
                        contract.shape,
                        findings,
                    );
                }
            }
            ToleranceTarget::Subsystem(_)
            | ToleranceTarget::Element(MachineElementId::Body(_))
            | ToleranceTarget::Element(MachineElementId::SurfacePatch(_))
            | ToleranceTarget::Element(MachineElementId::ContactFeature(_))
            | ToleranceTarget::Element(MachineElementId::Port(_)) => {}
        }

        let is_random = match &tolerance.semantics {
            ToleranceSemantics::Bounded { minus, plus, .. } => {
                if minus.is_zero() && plus.is_zero() {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::ZeroTolerance,
                        subject(),
                        None,
                    ));
                }
                false
            }
            ToleranceSemantics::Random { scale, .. } => {
                if scale.is_zero() {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::ZeroTolerance,
                        subject(),
                        None,
                    ));
                }
                true
            }
        };
        tolerance_random
            .entry(tolerance.id.clone())
            .or_insert(is_random);
    }

    let mut condition_random = BTreeMap::<ConditionTarget, bool>::new();
    for condition in conditions {
        condition_random
            .entry(condition.target.clone())
            .or_insert(matches!(
                &condition.source,
                ConditionSource::Distribution(_)
            ));
    }

    if dependences.len() > 1 {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::MultipleDependenceModels,
            MachineBehaviorSubject::Behavior,
            None,
        ));
    }
    let mut coverage = BTreeMap::<DependenceMember, Vec<Vec<DependenceMember>>>::new();
    for dependence in dependences {
        let subject = || MachineBehaviorSubject::Dependence(dependence.members.clone());
        if dependence.members.is_empty() {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::EmptyDependence,
                subject(),
                None,
            ));
        }
        for pair in dependence.members.windows(2) {
            if pair[0] == pair[1] {
                findings.push(MachineBehaviorFinding::new(
                    MachineBehaviorRule::DuplicateDependenceMember,
                    subject(),
                    Some(dependence_member_subject(&pair[1])),
                ));
            }
        }
        let distinct: BTreeSet<&DependenceMember> = dependence.members.iter().collect();
        if matches!(&dependence.model, DependenceModel::Correlated(_)) && distinct.len() < 2 {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::CorrelatedGroupTooSmall,
                subject(),
                None,
            ));
        }
        for member in distinct {
            let random = match member {
                DependenceMember::Condition(target) => match condition_random.get(target) {
                    None => {
                        findings.push(MachineBehaviorFinding::new(
                            MachineBehaviorRule::UnknownConditionMember,
                            subject(),
                            Some(MachineBehaviorSubject::Condition(target.clone())),
                        ));
                        false
                    }
                    Some(false) => {
                        findings.push(MachineBehaviorFinding::new(
                            MachineBehaviorRule::NonRandomConditionInDependence,
                            subject(),
                            Some(MachineBehaviorSubject::Condition(target.clone())),
                        ));
                        false
                    }
                    Some(true) => true,
                },
                DependenceMember::Tolerance(tolerance) => match tolerance_random.get(tolerance) {
                    None => {
                        findings.push(MachineBehaviorFinding::new(
                            MachineBehaviorRule::UnknownToleranceMember,
                            subject(),
                            Some(MachineBehaviorSubject::Tolerance(tolerance.clone())),
                        ));
                        false
                    }
                    Some(false) => {
                        findings.push(MachineBehaviorFinding::new(
                            MachineBehaviorRule::BoundedToleranceInDependence,
                            subject(),
                            Some(MachineBehaviorSubject::Tolerance(tolerance.clone())),
                        ));
                        false
                    }
                    Some(true) => true,
                },
            };
            if random {
                let member = member.clone();
                let groups = coverage.entry(member.clone()).or_default();
                if let Some(first) = groups.first() {
                    findings.push(MachineBehaviorFinding::new(
                        MachineBehaviorRule::DuplicateDependenceCoverage,
                        subject(),
                        Some(MachineBehaviorSubject::Dependence(first.clone())),
                    ));
                }
                groups.push(dependence.members.clone());
            }
        }
    }
    for (target, is_random) in condition_random {
        let member = DependenceMember::Condition(target);
        if is_random && !coverage.contains_key(&member) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingDependence,
                dependence_member_subject(&member),
                None,
            ));
        }
    }
    for (tolerance, is_random) in tolerance_random {
        let member = DependenceMember::Tolerance(tolerance);
        if is_random && !coverage.contains_key(&member) {
            findings.push(MachineBehaviorFinding::new(
                MachineBehaviorRule::MissingDependence,
                dependence_member_subject(&member),
                None,
            ));
        }
    }
}

fn dependence_member_subject(member: &DependenceMember) -> MachineBehaviorSubject {
    match member {
        DependenceMember::Condition(target) => MachineBehaviorSubject::Condition(target.clone()),
        DependenceMember::Tolerance(tolerance) => {
            MachineBehaviorSubject::Tolerance(tolerance.clone())
        }
    }
}

fn tolerance_target_subject(target: &ToleranceTarget) -> MachineBehaviorSubject {
    match target {
        ToleranceTarget::Subsystem(id) => MachineBehaviorSubject::Subsystem(id.clone()),
        ToleranceTarget::Element(id) => MachineBehaviorSubject::Element(id.clone()),
    }
}

fn compare_tolerance_contract(
    tolerance: &ToleranceSpec,
    quantity: TerminalQuantitySpec,
    shape: TerminalShape,
    findings: &mut Vec<MachineBehaviorFinding>,
) {
    if tolerance.quantity != quantity {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ToleranceQuantityGap,
            MachineBehaviorSubject::Tolerance(tolerance.id.clone()),
            None,
        ));
    }
    if tolerance.shape != shape {
        findings.push(MachineBehaviorFinding::new(
            MachineBehaviorRule::ToleranceShapeGap,
            MachineBehaviorSubject::Tolerance(tolerance.id.clone()),
            None,
        ));
    }
}

fn canonicalize_behavior_draft(draft: &mut MachineBehaviorDraft) {
    for condition in &mut draft.conditions {
        if let ConditionSource::History {
            continuity: HistoryContinuity::ResetAtEvents { events },
            ..
        } = &mut condition.source
        {
            events.sort();
        }
    }
    for event in &mut draft.events {
        event.dependencies.sort();
        if let Some(writes) = reset_writes_mut(&mut event.reset) {
            writes.sort();
        }
    }
    for dependence in &mut draft.dependences {
        dependence.members.sort();
    }

    draft.state_contracts.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| state_contract_row(left).cmp(&state_contract_row(right)))
    });
    draft.conditions.sort_by(|left, right| {
        left.target
            .cmp(&right.target)
            .then_with(|| condition_row(left).cmp(&condition_row(right)))
    });
    draft.motions.sort_by(|left, right| {
        left.body
            .cmp(&right.body)
            .then_with(|| motion_row(left).cmp(&motion_row(right)))
    });
    draft.events.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| event_row(left).cmp(&event_row(right)))
    });
    draft.tolerances.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| tolerance_row(left).cmp(&tolerance_row(right)))
    });
    draft
        .dependences
        .sort_by_key(|dependence| dependence_row(dependence));
}

fn machine_behavior_identity(
    draft: &MachineBehaviorDraft,
    graph: MachineGraphIdV1,
) -> Result<IdentityReceipt<MachineBehaviorIdV1>, CanonicalError> {
    let state_rows: Vec<Vec<u8>> = draft
        .state_contracts
        .iter()
        .map(state_contract_row)
        .collect();
    let condition_rows: Vec<Vec<u8>> = draft.conditions.iter().map(condition_row).collect();
    let motion_rows: Vec<Vec<u8>> = draft.motions.iter().map(motion_row).collect();
    let event_rows: Vec<Vec<u8>> = draft.events.iter().map(event_row).collect();
    let tolerance_rows: Vec<Vec<u8>> = draft.tolerances.iter().map(tolerance_row).collect();
    let dependence_rows: Vec<Vec<u8>> = draft.dependences.iter().map(dependence_row).collect();

    CanonicalEncoder::<MachineBehaviorIdV1, _>::new(MACHINE_BEHAVIOR_IDENTITY_LIMITS, NeverCancel)?
        .u64(
            Field::new(0, "behavior-schema-version"),
            u64::from(MACHINE_BEHAVIOR_SCHEMA_VERSION_V1),
        )?
        .bytes(Field::new(1, "base-machine-graph"), graph.as_bytes())?
        .ordered_bytes(
            Field::new(2, "state-contracts"),
            state_rows.len() as u64,
            state_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(3, "conditions"),
            condition_rows.len() as u64,
            condition_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(4, "motions"),
            motion_rows.len() as u64,
            motion_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(5, "events"),
            event_rows.len() as u64,
            event_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(6, "tolerances"),
            tolerance_rows.len() as u64,
            tolerance_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(7, "dependences"),
            dependence_rows.len() as u64,
            dependence_rows.iter().map(Vec::as_slice),
        )?
        .finish()
}

fn state_contract_row(contract: &StateSlotContract) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_identity(&mut out, &contract.id.digest_bytes());
    push_identity(&mut out, &contract.owner.digest_bytes());
    super::push_terminal_quantity(&mut out, contract.quantity);
    super::push_terminal_shape(&mut out, contract.shape);
    push_identity(&mut out, &contract.clock.digest_bytes());
    push_frame(&mut out, &contract.frame);
    out
}

fn condition_row(condition: &ConditionBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    match &condition.target {
        ConditionTarget::Initial(state) => {
            out.push(1);
            push_identity(&mut out, &state.digest_bytes());
        }
        ConditionTarget::Boundary(terminal) => {
            out.push(2);
            push_identity(&mut out, &terminal.digest_bytes());
        }
    }
    super::push_terminal_quantity(&mut out, condition.quantity);
    super::push_terminal_shape(&mut out, condition.shape);
    push_identity(&mut out, &condition.clock.digest_bytes());
    push_frame(&mut out, &condition.frame);
    push_condition_source(&mut out, &condition.source);
    out
}

fn push_condition_source(out: &mut Vec<u8>, source: &ConditionSource) {
    match source {
        ConditionSource::Fixed(value) => {
            out.push(1);
            value.append_canonical(out);
        }
        ConditionSource::History {
            history,
            continuity,
        } => {
            out.push(2);
            history.append_canonical(out);
            match continuity {
                HistoryContinuity::Continuous => out.push(1),
                HistoryContinuity::ResetAtEvents { events } => {
                    out.push(2);
                    out.extend_from_slice(&(events.len() as u64).to_le_bytes());
                    for event in events {
                        push_identity(out, &event.digest_bytes());
                    }
                }
            }
        }
        ConditionSource::Distribution(distribution) => {
            out.push(3);
            distribution.append_canonical(out);
        }
    }
}

fn motion_row(motion: &MotionBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_identity(&mut out, &motion.body.digest_bytes());
    push_identity(&mut out, &motion.clock.digest_bytes());
    push_frame(&mut out, &motion.reference_frame);
    match &motion.motion {
        BodyMotion::Static => out.push(1),
        BodyMotion::Prescribed { path } => {
            out.push(2);
            path.append_canonical(&mut out);
        }
    }
    out
}

fn event_row(event: &EventSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(512);
    push_identity(&mut out, &event.id.digest_bytes());
    push_identity(&mut out, &event.clock.digest_bytes());
    event.guard.append_canonical(&mut out);
    out.push(match event.orientation {
        GuardOrientation::NegativeToPositive => 1,
        GuardOrientation::PositiveToNegative => 2,
        GuardOrientation::Bidirectional => 3,
    });
    match &event.crossing {
        CrossingSemantics::Transverse(witness) => {
            out.push(1);
            witness.append_canonical(&mut out);
        }
        CrossingSemantics::Grazing(witness) => {
            out.push(2);
            witness.append_canonical(&mut out);
        }
        CrossingSemantics::Unknown(no_claim) => {
            out.push(3);
            no_claim.append_canonical(&mut out);
        }
    }
    out.extend_from_slice(&(event.dependencies.len() as u64).to_le_bytes());
    for dependency in &event.dependencies {
        match dependency {
            EventDependency::State(state) => {
                out.push(1);
                push_identity(&mut out, &state.digest_bytes());
            }
            EventDependency::Terminal(terminal) => {
                out.push(2);
                push_identity(&mut out, &terminal.digest_bytes());
            }
        }
    }
    push_reset(&mut out, &event.reset);
    match &event.order {
        EventOrder::TotalPriority { microstep } => {
            out.push(1);
            out.extend_from_slice(&microstep.to_le_bytes());
        }
        EventOrder::SetValued { group } => {
            out.push(2);
            group.append_canonical(&mut out);
        }
    }
    out
}

fn push_reset(out: &mut Vec<u8>, reset: &ResetSemantics) {
    match reset {
        ResetSemantics::Deterministic { map, writes } => {
            out.push(1);
            map.append_canonical(out);
            push_state_set(out, writes);
        }
        ResetSemantics::SetValued {
            relation,
            outcomes,
            writes,
        } => {
            out.push(2);
            relation.append_canonical(out);
            outcomes.append_canonical(out);
            push_state_set(out, writes);
        }
        ResetSemantics::Terminal { relation } => {
            out.push(3);
            relation.append_canonical(out);
        }
        ResetSemantics::Unknown { no_claim } => {
            out.push(4);
            no_claim.append_canonical(out);
        }
    }
}

fn push_state_set(out: &mut Vec<u8>, states: &[StateSlotId]) {
    out.extend_from_slice(&(states.len() as u64).to_le_bytes());
    for state in states {
        push_identity(out, &state.digest_bytes());
    }
}

fn tolerance_row(tolerance: &ToleranceSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(320);
    push_identity(&mut out, &tolerance.id.digest_bytes());
    match &tolerance.target {
        ToleranceTarget::Subsystem(subsystem) => {
            out.push(1);
            push_identity(&mut out, &subsystem.digest_bytes());
        }
        ToleranceTarget::Element(element) => {
            out.push(2);
            out.push(element.kind().tag());
            push_identity(&mut out, &element.digest_bytes());
        }
    }
    tolerance.parameter.append_canonical(&mut out);
    super::push_terminal_quantity(&mut out, tolerance.quantity);
    super::push_terminal_shape(&mut out, tolerance.shape);
    match &tolerance.semantics {
        ToleranceSemantics::Bounded { minus, plus, law } => {
            out.push(1);
            out.extend_from_slice(&minus.bits().to_le_bytes());
            out.extend_from_slice(&plus.bits().to_le_bytes());
            law.append_canonical(&mut out);
        }
        ToleranceSemantics::Random {
            scale,
            law,
            marginal,
        } => {
            out.push(2);
            out.extend_from_slice(&scale.bits().to_le_bytes());
            law.append_canonical(&mut out);
            marginal.append_canonical(&mut out);
        }
    }
    out
}

fn dependence_row(dependence: &DependenceSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(64 + 32 * dependence.members.len());
    out.extend_from_slice(&(dependence.members.len() as u64).to_le_bytes());
    for member in &dependence.members {
        match member {
            DependenceMember::Condition(ConditionTarget::Initial(state)) => {
                out.push(1);
                push_identity(&mut out, &state.digest_bytes());
            }
            DependenceMember::Condition(ConditionTarget::Boundary(terminal)) => {
                out.push(2);
                push_identity(&mut out, &terminal.digest_bytes());
            }
            DependenceMember::Tolerance(tolerance) => {
                out.push(3);
                push_identity(&mut out, &tolerance.digest_bytes());
            }
        }
    }
    match &dependence.model {
        DependenceModel::Independent => out.push(1),
        DependenceModel::Correlated(model) => {
            out.push(2);
            model.append_canonical(&mut out);
        }
    }
    out
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

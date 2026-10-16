//! Versioned hybrid-time, event-order, and Zeno semantics (RE.Z1).
//!
//! This module defines the finite admitted object language used by later
//! interval/event checkers. It does not infer Zeno behavior from an event cap,
//! a dense numerical trace, or an opaque witness identifier. Successful
//! validation proves schema consistency and deterministic identity only; the
//! scientific promotion checker belongs to RE.Z2.

use core::fmt;
use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_exec::Cx;

/// Current hybrid/Zeno problem schema version.
pub const ZENO_PROBLEM_SCHEMA_VERSION_V1: u32 = 1;
/// Current hybrid/Zeno claim-descriptor schema version.
pub const ZENO_CLAIM_SCHEMA_VERSION_V1: u32 = 1;
/// Hard pre-work cap on modes in one problem.
pub const MAX_HYBRID_MODES_V1: usize = 256;
/// Hard pre-work cap on events in one problem.
pub const MAX_HYBRID_EVENTS_V1: usize = 1024;
/// Hard pre-work cap on targets of one set-valued reset.
pub const MAX_RESET_TARGETS_V1: usize = 256;
/// Hard pre-work cap on the sum of all reset targets.
pub const MAX_TOTAL_RESET_TARGETS_V1: usize = 4096;

const HYBRID_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 20, 1 << 20, 16, 8192, 8192);

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
}

macro_rules! opaque_hybrid_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct from exact typed digest bytes. Identity alone is not
            /// scientific authority.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact typed digest bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }

        impl DigestBytes for $name {
            fn digest_bytes(&self) -> &[u8; 32] {
                self.as_bytes()
            }
        }
    };
}

opaque_hybrid_id!(
    /// Content identity of a hybrid model.
    HybridModelIdV1
);
opaque_hybrid_id!(
    /// Immutable version identity of a hybrid model.
    HybridModelVersionIdV1
);
opaque_hybrid_id!(
    /// Identity of one discrete mode.
    HybridModeIdV1
);
opaque_hybrid_id!(
    /// Identity of one mode's continuous dynamics.
    ContinuousDynamicsIdV1
);
opaque_hybrid_id!(
    /// Identity of a DAE constraint manifold.
    DaeConstraintIdV1
);
opaque_hybrid_id!(
    /// Identity of an event guard.
    HybridGuardIdV1
);
opaque_hybrid_id!(
    /// Identity of a reset map or relation.
    ResetRelationIdV1
);
opaque_hybrid_id!(
    /// Identity of a hybrid event symbol.
    HybridEventIdV1
);
opaque_hybrid_id!(
    /// Identity of a contact, relay, or other interaction law.
    InteractionLawIdV1
);
opaque_hybrid_id!(
    /// Identity of an event-word/language specification.
    EventLanguageIdV1
);
opaque_hybrid_id!(
    /// Identity of a simultaneity group.
    SimultaneityGroupIdV1
);
opaque_hybrid_id!(
    /// Identity of a physical state or continuation state.
    HybridStateIdV1
);
opaque_hybrid_id!(
    /// Identity of a set of physical states.
    HybridStateSetIdV1
);
opaque_hybrid_id!(
    /// Identity of a frame convention.
    HybridFrameIdV1
);
opaque_hybrid_id!(
    /// Identity of a state/unit convention.
    HybridUnitSystemIdV1
);
opaque_hybrid_id!(
    /// Identity of a time unit/clock convention.
    HybridTimeUnitIdV1
);
opaque_hybrid_id!(
    /// Identity of retained analytic or checker evidence.
    HybridWitnessIdV1
);
opaque_hybrid_id!(
    /// Identity of a trace that is numerical evidence only.
    HybridEventTraceIdV1
);
opaque_hybrid_id!(
    /// Identity of an explicit no-claim statement.
    HybridNoClaimIdV1
);
opaque_hybrid_id!(
    /// Identity of a compliant/smoothed/regularized model transformation.
    HybridRegularizationIdV1
);
opaque_hybrid_id!(
    /// Identity of a post-accumulation continuation rule.
    ContinuationRuleIdV1
);

/// Domain-separated identity schema for one admitted hybrid problem.
pub enum ZenoProblemIdentitySchemaV1 {}

impl CanonicalSchema for ZenoProblemIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-time.zeno-problem.v1";
    const NAME: &'static str = "hybrid-zeno-problem";
    const VERSION: u32 = ZENO_PROBLEM_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "hybrid model, modes, guards, resets, event order, time domain, accumulation, continuation, and budgets";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("model", WireType::Bytes),
        FieldSpec::required("physical-context", WireType::Bytes),
        FieldSpec::required("modes", WireType::CanonicalSet),
        FieldSpec::required("events", WireType::CanonicalSet),
        FieldSpec::required("event-language", WireType::Bytes),
        FieldSpec::required("simultaneous-policy", WireType::Bytes),
        FieldSpec::required("hybrid-time", WireType::Bytes),
        FieldSpec::required("compactness", WireType::Bytes),
        FieldSpec::required("accumulation-candidate", WireType::Bytes),
        FieldSpec::required("continuation", WireType::Bytes),
        FieldSpec::required("analysis-budget", WireType::Bytes),
    ];
}

/// Typed identity of one canonical admitted hybrid problem.
pub type ZenoProblemIdV1 = ProblemSemanticId<ZenoProblemIdentitySchemaV1>;

/// Domain-separated identity schema for one claim descriptor.
pub enum ZenoClaimIdentitySchemaV1 {}

impl CanonicalSchema for ZenoClaimIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-time.zeno-claim.v1";
    const NAME: &'static str = "hybrid-zeno-claim-descriptor";
    const VERSION: u32 = ZENO_CLAIM_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "problem-bound finite-separation, Zeno, warning, regularization, or unknown descriptor plus post-state semantics";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("problem", WireType::Bytes),
        FieldSpec::required("classification", WireType::Bytes),
        FieldSpec::required("post-zeno", WireType::Bytes),
    ];
}

/// Typed identity of one canonical claim descriptor. This is not a theorem.
pub type ZenoClaimIdV1 = ProblemSemanticId<ZenoClaimIdentitySchemaV1>;

/// Original-versus-regularized model lineage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridModelLineageV1 {
    /// Original physical hybrid model.
    Original,
    /// A distinct model produced by a named regularization. No equivalence to
    /// the original event semantics is implied.
    Regularized {
        /// Source model identity.
        source_model: HybridModelIdV1,
        /// Source model version.
        source_version: HybridModelVersionIdV1,
        /// Exact regularization transformation.
        regularization: HybridRegularizationIdV1,
        /// Explicit original-versus-regularized no-equivalence boundary.
        no_equivalence: HybridNoClaimIdV1,
    },
}

/// Supported finite continuous-dynamics classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuousDynamicsClassV1 {
    /// Deterministic finite-dimensional ODE.
    DeterministicOde,
    /// Finite-dimensional differential inclusion; evolution may be set-valued.
    DifferentialInclusion,
    /// Explicitly admitted finite-dimensional DAE.
    AdmittedDae {
        /// Declared positive differentiation index.
        index: u8,
        /// Exact constraint-manifold identity.
        constraint: DaeConstraintIdV1,
    },
    /// Unsupported infinite-dimensional dynamics, retained only so decoded
    /// input can fail closed.
    UnsupportedInfiniteDimensional,
}

/// One finite hybrid mode and its continuous dynamics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HybridModeSpecV1 {
    /// Mode identity.
    pub mode: HybridModeIdV1,
    /// Exact dynamics artifact.
    pub dynamics: ContinuousDynamicsIdV1,
    /// Continuous state dimension.
    pub state_dimension: u32,
    /// Mathematical dynamics class.
    pub class: ContinuousDynamicsClassV1,
}

/// Oriented crossing convention for a scalar guard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GuardOrientationV1 {
    /// Negative guard value to positive.
    NegativeToPositive,
    /// Positive guard value to negative.
    PositiveToNegative,
    /// Either direction may trigger; unique continuation is not implied.
    Bidirectional,
    /// Direction is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Local guard-crossing semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingSemanticsV1 {
    /// Exact transversality evidence is retained.
    Transverse {
        /// Proposition/checker witness.
        witness: HybridWitnessIdV1,
    },
    /// A grazing/tangent contact is part of the admitted model.
    Grazing {
        /// Grazing-classification witness.
        witness: HybridWitnessIdV1,
    },
    /// Crossing class is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Reset semantics following one event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResetSemanticsV1 {
    /// Deterministic map into one target mode.
    Deterministic {
        /// Exact reset map.
        relation: ResetRelationIdV1,
        /// Target mode.
        target: HybridModeIdV1,
    },
    /// Set-valued reset relation into one or more target modes.
    SetValued {
        /// Exact reset relation.
        relation: ResetRelationIdV1,
        /// Possible target modes. Validation canonicalizes this set.
        targets: Vec<HybridModeIdV1>,
        /// Exact post-reset state-set artifact.
        states: HybridStateSetIdV1,
    },
    /// Event terminates the hybrid execution.
    Terminal {
        /// Exact terminal reset/restriction relation.
        relation: ResetRelationIdV1,
    },
    /// Reset behavior is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Physical law attached to an event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InteractionLawV1 {
    /// No contact/relay interaction applies, with an explicit justification.
    None {
        /// Non-applicability artifact.
        justification: HybridNoClaimIdV1,
    },
    /// Contact/impact/friction law.
    Contact {
        /// Exact law identity.
        law: InteractionLawIdV1,
    },
    /// Relay/switching/control law.
    Relay {
        /// Exact law identity.
        law: InteractionLawIdV1,
    },
    /// Other exact event interaction law.
    Other {
        /// Exact law identity.
        law: InteractionLawIdV1,
    },
}

/// Whether one event is exclusive or belongs to a simultaneous group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventSimultaneityV1 {
    /// Event cannot coincide with another admitted event.
    Exclusive {
        /// Exclusivity witness.
        witness: HybridWitnessIdV1,
    },
    /// Event may fire with every event carrying the same group identity.
    Group {
        /// Simultaneity group.
        group: SimultaneityGroupIdV1,
    },
}

/// Lower-bound semantics for inter-event dwell time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DwellSemanticsV1 {
    /// A strictly positive physical-time lower bound.
    PositiveLowerBound {
        /// Lower bound in the declared hybrid-time unit.
        value: f64,
        /// Exact lower-bound witness.
        witness: HybridWitnessIdV1,
    },
    /// Zero dwell is admitted, so a zero-time cycle may exist.
    ZeroAllowed,
    /// No positive dwell claim is available.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// One guard/reset event in the hybrid automaton.
#[derive(Debug, Clone, PartialEq)]
pub struct HybridEventSpecV1 {
    /// Event-symbol identity.
    pub event: HybridEventIdV1,
    /// Source mode.
    pub source_mode: HybridModeIdV1,
    /// Exact guard identity.
    pub guard: HybridGuardIdV1,
    /// Guard orientation.
    pub orientation: GuardOrientationV1,
    /// Transverse, grazing, or unresolved crossing semantics.
    pub crossing: CrossingSemanticsV1,
    /// Reset map/relation semantics.
    pub reset: ResetSemanticsV1,
    /// Contact, relay, or other physical interaction law.
    pub law: InteractionLawV1,
    /// Exclusive or simultaneous-event group semantics.
    pub simultaneity: EventSimultaneityV1,
    /// Inter-event dwell semantics.
    pub dwell: DwellSemanticsV1,
}

/// Event-word/language semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventLanguageSemanticsV1 {
    /// Only finite words up to an explicit length are in scope.
    FiniteWords {
        /// Maximum admitted event-word length.
        max_events_per_word: u32,
    },
    /// Prefix-closed finite-event language without a fixed semantic horizon.
    PrefixClosed,
    /// Infinite event words are part of the mathematical language.
    OmegaLanguage,
    /// Language semantics are unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Exact event alphabet/language reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EventLanguageSpecV1 {
    /// Language artifact identity.
    pub language: EventLanguageIdV1,
    /// Finite, prefix-closed, omega, or unknown semantics.
    pub semantics: EventLanguageSemanticsV1,
}

/// Global interpretation of simultaneous groups.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SimultaneousEventPolicyV1 {
    /// Model asserts that no simultaneous group can occur.
    NoSimultaneousEvents {
        /// Exact exclusion witness.
        witness: HybridWitnessIdV1,
    },
    /// Exact deterministic priority over all group-member events.
    TotalPriority {
        /// Highest-to-lowest priority. Every grouped event appears once.
        ordered_events: Vec<HybridEventIdV1>,
        /// Priority-policy witness.
        witness: HybridWitnessIdV1,
    },
    /// All simultaneous resets commute under a retained theorem.
    Commuting {
        /// Exact commutation witness.
        witness: HybridWitnessIdV1,
    },
    /// Simultaneous firing remains explicitly set-valued.
    SetValued {
        /// Exact outcome-set semantics.
        outcomes: HybridStateSetIdV1,
    },
    /// No event-order selection is claimed.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Physical time-unit conversion carried by the model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HybridTimeScaleV1 {
    /// Time-unit/clock identity.
    pub unit: HybridTimeUnitIdV1,
    /// Positive SI seconds represented by one declared unit.
    pub seconds_per_unit: f64,
}

/// End of the admitted hybrid-time domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HybridTimeEndV1 {
    /// Finite inclusive upper bound in declared time units.
    Finite(f64),
    /// Unbounded future time.
    Infinite,
}

/// Hybrid-time domain and execution-only event cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HybridTimeDomainV1 {
    /// Inclusive start time in declared units.
    pub start: f64,
    /// Finite or infinite end.
    pub end: HybridTimeEndV1,
    /// Optional execution cap. This is never theorem evidence.
    pub event_cap: Option<u64>,
}

/// Compactness/properness state used by later theorem checkers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompactnessSemanticsV1 {
    /// Relevant state/time tube is compact.
    Compact {
        /// Exact compactness witness.
        witness: HybridWitnessIdV1,
    },
    /// Only a named local compactness region is established.
    LocallyCompact {
        /// Exact local-region witness.
        witness: HybridWitnessIdV1,
    },
    /// Compactness is unestablished.
    Unestablished {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// A numerical candidate window, never a Zeno theorem by itself.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AccumulationCandidateV1 {
    /// No accumulation candidate is asserted.
    None {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
    /// Candidate time/state enclosure from a retained numerical trace.
    Window {
        /// Earliest candidate time in declared units.
        earliest: f64,
        /// Latest candidate time in declared units.
        latest: f64,
        /// Candidate state enclosure.
        states: HybridStateSetIdV1,
        /// Numerical trace identity.
        trace: HybridEventTraceIdV1,
    },
}

/// Admitted post-accumulation continuation category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuationCategoryV1 {
    /// A unique continuation map is claimed by retained evidence.
    Unique {
        /// Exact continuation rule.
        rule: ContinuationRuleIdV1,
        /// Uniqueness witness.
        witness: HybridWitnessIdV1,
    },
    /// Continuation is an exact set-valued relation.
    SetValued {
        /// Exact relation identity.
        rule: ContinuationRuleIdV1,
    },
    /// Every admitted accumulation terminates execution.
    Terminal {
        /// Exact terminal rule.
        rule: ContinuationRuleIdV1,
    },
    /// No post-Zeno selection is claimed.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Explicit semantic-analysis budget. It is not an execution receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HybridAnalysisBudgetV1 {
    /// Maximum event-word length a checker may inspect.
    pub max_event_word_len: u32,
    /// Maximum transition steps a checker may inspect.
    pub max_transitions: u64,
    /// Positive wall-time allowance in seconds.
    pub max_wall_seconds: f64,
}

/// Raw versioned hybrid/Zeno problem. It has no authority until validated.
#[derive(Debug, Clone, PartialEq)]
pub struct ZenoProblemIrV1 {
    schema_version: u32,
    /// Exact model identity.
    pub model: HybridModelIdV1,
    /// Exact immutable model version.
    pub model_version: HybridModelVersionIdV1,
    /// Original or explicitly regularized lineage.
    pub lineage: HybridModelLineageV1,
    /// State frame.
    pub frame: HybridFrameIdV1,
    /// State/parameter unit schema.
    pub state_units: HybridUnitSystemIdV1,
    /// Exact initial-state set to which all event/accumulation semantics apply.
    pub initial_states: HybridStateSetIdV1,
    /// Physical time scale.
    pub time_scale: HybridTimeScaleV1,
    /// Finite mode set.
    pub modes: Vec<HybridModeSpecV1>,
    /// Finite event set.
    pub events: Vec<HybridEventSpecV1>,
    /// Event language.
    pub event_language: EventLanguageSpecV1,
    /// Simultaneous-event policy.
    pub simultaneous_policy: SimultaneousEventPolicyV1,
    /// Hybrid-time domain.
    pub time_domain: HybridTimeDomainV1,
    /// Compactness state.
    pub compactness: CompactnessSemanticsV1,
    /// Optional numerical accumulation candidate.
    pub accumulation_candidate: AccumulationCandidateV1,
    /// Post-accumulation continuation category.
    pub continuation: ContinuationCategoryV1,
    /// Semantic analysis budget.
    pub budget: HybridAnalysisBudgetV1,
}

impl ZenoProblemIrV1 {
    /// Construct current-version raw input.
    #[allow(clippy::too_many_arguments)] // Every orthogonal hybrid convention remains explicit.
    #[must_use]
    pub fn new(
        model: HybridModelIdV1,
        model_version: HybridModelVersionIdV1,
        lineage: HybridModelLineageV1,
        frame: HybridFrameIdV1,
        state_units: HybridUnitSystemIdV1,
        initial_states: HybridStateSetIdV1,
        time_scale: HybridTimeScaleV1,
        modes: Vec<HybridModeSpecV1>,
        events: Vec<HybridEventSpecV1>,
        event_language: EventLanguageSpecV1,
        simultaneous_policy: SimultaneousEventPolicyV1,
        time_domain: HybridTimeDomainV1,
        compactness: CompactnessSemanticsV1,
        accumulation_candidate: AccumulationCandidateV1,
        continuation: ContinuationCategoryV1,
        budget: HybridAnalysisBudgetV1,
    ) -> Self {
        Self::with_schema_version(
            ZENO_PROBLEM_SCHEMA_VERSION_V1,
            model,
            model_version,
            lineage,
            frame,
            state_units,
            initial_states,
            time_scale,
            modes,
            events,
            event_language,
            simultaneous_policy,
            time_domain,
            compactness,
            accumulation_candidate,
            continuation,
            budget,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[allow(clippy::too_many_arguments)] // Decoded version plus every semantic axis.
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        model: HybridModelIdV1,
        model_version: HybridModelVersionIdV1,
        lineage: HybridModelLineageV1,
        frame: HybridFrameIdV1,
        state_units: HybridUnitSystemIdV1,
        initial_states: HybridStateSetIdV1,
        time_scale: HybridTimeScaleV1,
        modes: Vec<HybridModeSpecV1>,
        events: Vec<HybridEventSpecV1>,
        event_language: EventLanguageSpecV1,
        simultaneous_policy: SimultaneousEventPolicyV1,
        time_domain: HybridTimeDomainV1,
        compactness: CompactnessSemanticsV1,
        accumulation_candidate: AccumulationCandidateV1,
        continuation: ContinuationCategoryV1,
        budget: HybridAnalysisBudgetV1,
    ) -> Self {
        Self {
            schema_version,
            model,
            model_version,
            lineage,
            frame,
            state_units,
            initial_states,
            time_scale,
            modes,
            events,
            event_language,
            simultaneous_policy,
            time_domain,
            compactness,
            accumulation_candidate,
            continuation,
            budget,
        }
    }

    /// Declared schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }
}

/// Closed candidate/theorem time interval in declared hybrid-time units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HybridTimeIntervalV1 {
    /// Inclusive lower endpoint.
    pub earliest: f64,
    /// Inclusive upper endpoint.
    pub latest: f64,
}

/// Evidence class named by a theorem-shaped descriptor.
///
/// These values are content references, not authority tokens. RE.Z2 must
/// independently check and promote them before a scientific theorem exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZenoEvidenceReferenceV1 {
    /// Analytic derivation/checker artifact.
    Analytic {
        /// Exact evidence identity.
        witness: HybridWitnessIdV1,
    },
    /// Independently interval-validated event/accumulation artifact.
    IntervalValidated {
        /// Exact evidence identity.
        witness: HybridWitnessIdV1,
    },
    /// Execution stopped at an event-count cap. Never theorem evidence.
    EventCap {
        /// Exact numerical trace.
        trace: HybridEventTraceIdV1,
    },
    /// Dense numerical events without a theorem checker.
    NumericalOnly {
        /// Exact numerical trace.
        trace: HybridEventTraceIdV1,
    },
    /// No favorable evidence is claimed.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Mutually exclusive hybrid/Zeno classifications.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ZenoClassificationV1 {
    /// A positive lower bound separates every consecutive event in scope.
    FiniteEventSeparation {
        /// Positive separation in declared time units.
        minimum_separation: f64,
        /// Theorem-shaped evidence reference.
        evidence: ZenoEvidenceReferenceV1,
    },
    /// Event accumulation is asserted inside a closed time interval.
    CertifiedZeno {
        /// Claimed accumulation-time enclosure.
        interval: HybridTimeIntervalV1,
        /// Claimed accumulation-state enclosure.
        states: HybridStateSetIdV1,
        /// Theorem-shaped evidence reference.
        evidence: ZenoEvidenceReferenceV1,
    },
    /// Numerical event density is concerning but proves no accumulation.
    NumericalEventDensityWarning {
        /// Retained trace.
        trace: HybridEventTraceIdV1,
        /// Events observed in the stated window.
        observed_events: u64,
        /// Positive observation-window width in declared time units.
        window: f64,
    },
    /// A distinct regularized model is being used instead of the source.
    RegularizedModel {
        /// Exact validated regularized-problem identity.
        regularized_problem: ZenoProblemIdV1,
        /// Exact regularization transformation.
        regularization: HybridRegularizationIdV1,
        /// Explicit no-equivalence boundary.
        no_equivalence: HybridNoClaimIdV1,
    },
    /// No finite-separation or Zeno classification is claimed.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Post-accumulation state semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PostZenoStateV1 {
    /// No post-Zeno state applies, as for finite-event separation.
    NotApplicable {
        /// Explicit non-applicability artifact.
        justification: HybridNoClaimIdV1,
    },
    /// One exact continuation state.
    Unique {
        /// State identity.
        state: HybridStateIdV1,
        /// Continuation rule.
        rule: ContinuationRuleIdV1,
        /// Uniqueness evidence reference.
        witness: HybridWitnessIdV1,
    },
    /// Exact set of possible continuation states.
    SetValued {
        /// State-set identity.
        states: HybridStateSetIdV1,
        /// Continuation relation.
        rule: ContinuationRuleIdV1,
        /// Set-containment evidence reference.
        witness: HybridWitnessIdV1,
    },
    /// Hybrid execution terminates at accumulation.
    Terminal {
        /// Exact terminal rule.
        rule: ContinuationRuleIdV1,
        /// Terminality evidence reference.
        witness: HybridWitnessIdV1,
    },
    /// Post-accumulation semantics remain unresolved.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: HybridNoClaimIdV1,
    },
}

/// Raw problem-bound claim descriptor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZenoClaimDraftV1 {
    schema_version: u32,
    problem: ZenoProblemIdV1,
    classification: ZenoClassificationV1,
    post_zeno: PostZenoStateV1,
}

impl ZenoClaimDraftV1 {
    /// Construct a current-version claim descriptor.
    #[must_use]
    pub const fn new(
        problem: ZenoProblemIdV1,
        classification: ZenoClassificationV1,
        post_zeno: PostZenoStateV1,
    ) -> Self {
        Self::with_schema_version(
            ZENO_CLAIM_SCHEMA_VERSION_V1,
            problem,
            classification,
            post_zeno,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[must_use]
    pub const fn with_schema_version(
        schema_version: u32,
        problem: ZenoProblemIdV1,
        classification: ZenoClassificationV1,
        post_zeno: PostZenoStateV1,
    ) -> Self {
        Self {
            schema_version,
            problem,
            classification,
            post_zeno,
        }
    }

    /// Declared schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact target problem identity.
    #[must_use]
    pub const fn problem(&self) -> ZenoProblemIdV1 {
        self.problem
    }

    /// Classification descriptor.
    #[must_use]
    pub const fn classification(&self) -> ZenoClassificationV1 {
        self.classification
    }

    /// Post-accumulation state semantics.
    #[must_use]
    pub const fn post_zeno(&self) -> PostZenoStateV1 {
        self.post_zeno
    }
}

/// Collections protected by hard pre-work caps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridCollectionV1 {
    /// Mode set.
    Modes,
    /// Event set.
    Events,
    /// Targets of one set-valued reset.
    ResetTargets,
    /// Sum of all reset targets.
    TotalResetTargets,
    /// Simultaneous-event priority sequence.
    PriorityEvents,
}

/// Identity-bearing object kind associated with a duplicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridIdKindV1 {
    /// Mode identity.
    Mode,
    /// Event identity.
    Event,
    /// Reset target identity within one set.
    ResetTarget,
    /// Priority event identity.
    PriorityEvent,
}

/// Numeric or semantic field associated with a refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridFieldV1 {
    /// Continuous state dimension.
    StateDimension,
    /// DAE index.
    DaeIndex,
    /// Seconds per time unit.
    SecondsPerUnit,
    /// Hybrid-time start.
    TimeStart,
    /// Hybrid-time end.
    TimeEnd,
    /// Execution event cap.
    EventCap,
    /// Positive dwell bound.
    DwellLowerBound,
    /// Candidate accumulation window.
    AccumulationWindow,
    /// Event-word budget.
    EventWordBudget,
    /// Transition budget.
    TransitionBudget,
    /// Wall-time budget.
    WallTimeBudget,
    /// Finite-separation lower bound.
    FiniteEventSeparation,
    /// Numerical-warning observation count.
    ObservedEvents,
    /// Numerical-warning window.
    ObservationWindow,
    /// Claimed Zeno interval.
    ZenoInterval,
}

/// Deterministic fail-closed issue for problem or claim admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HybridSemanticIssueV1 {
    /// Unsupported problem or claim schema.
    UnsupportedSchemaVersion {
        /// Version supplied.
        found: u32,
        /// Supported version.
        supported: u32,
    },
    /// A hard collection cap was exceeded before sorting/graph work.
    TooMany {
        /// Bounded collection.
        collection: HybridCollectionV1,
        /// Items supplied.
        found: usize,
        /// Maximum admitted.
        limit: usize,
    },
    /// A mandatory finite collection is empty.
    EmptyCollection {
        /// Empty collection.
        collection: HybridCollectionV1,
    },
    /// Duplicate identity in a set.
    DuplicateId {
        /// Identity kind.
        kind: HybridIdKindV1,
    },
    /// A numeric value is non-finite, zero, negative, or otherwise invalid.
    InvalidValue {
        /// Invalid field.
        field: HybridFieldV1,
    },
    /// An event source or reset target names no admitted mode.
    UnknownModeReference,
    /// A priority sequence omits, duplicates, or adds grouped events.
    InvalidPriorityOrder,
    /// A declared simultaneity group contains fewer than two events.
    SingletonSimultaneityGroup,
    /// Simultaneous group data conflicts with the global policy.
    SimultaneousPolicyMismatch,
    /// Unsupported infinite-dimensional continuous dynamics.
    UnsupportedDynamicsClass,
    /// Regularized lineage aliases its own source model/version.
    RegularizationSelfReference,
    /// Declared finite event language exceeds the explicit analysis budget.
    EventLanguageBudgetMismatch,
    /// Candidate accumulation window lies outside the hybrid-time domain.
    AccumulationOutsideTimeDomain,
    /// Unique continuation conflicts with set-valued, grazing, unresolved, or
    /// unresolved-simultaneous local semantics.
    UniqueContinuationUnsupported,
    /// A finite-separation descriptor conflicts with any admitted event that
    /// lacks a positive dwell lower bound.
    FiniteSeparationContradictsEventGraph,
    /// A Zeno descriptor lacks the zero/unknown-dwell cycle required for
    /// infinitely many events in finite time in this finite automaton.
    ZenoAccumulationCycleRequired,
    /// Claim descriptor names a different problem.
    TargetProblemMismatch,
    /// Event-count or numerical-only evidence was used for a theorem-shaped
    /// classification.
    InsufficientTheoremEvidence,
    /// Classification and post-Zeno state semantics conflict.
    PostZenoSemanticsMismatch,
    /// Claimed Zeno interval lies outside the admitted hybrid-time domain.
    ClaimOutsideTimeDomain,
    /// Regularized-model claim lacks an exact matching validated lineage.
    RegularizationLineageMismatch,
    /// Cooperative cancellation occurred before publication.
    Cancelled,
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Deterministic fail-closed refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HybridSemanticReportV1 {
    issues: Vec<HybridSemanticIssueV1>,
}

impl HybridSemanticReportV1 {
    fn new(issues: Vec<HybridSemanticIssueV1>) -> Self {
        Self { issues }
    }

    /// Deterministically ordered issues.
    #[must_use]
    pub fn issues(&self) -> &[HybridSemanticIssueV1] {
        &self.issues
    }
}

impl fmt::Display for HybridSemanticReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "hybrid/Zeno semantics refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for HybridSemanticReportV1 {}

/// Sealed canonical hybrid problem plus conservative graph facts.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedZenoProblemV1 {
    ir: ZenoProblemIrV1,
    receipt: IdentityReceipt<ZenoProblemIdV1>,
    has_zero_time_cycle: bool,
    has_nonunique_local_semantics: bool,
}

impl ValidatedZenoProblemV1 {
    /// Observational canonical problem view; not scientific authority.
    #[must_use]
    pub const fn ir(&self) -> &ZenoProblemIrV1 {
        &self.ir
    }

    /// Typed problem identity.
    #[must_use]
    pub const fn problem_id(&self) -> ZenoProblemIdV1 {
        self.receipt.id()
    }

    /// Exact identity and canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<ZenoProblemIdV1> {
        self.receipt
    }

    /// Whether zero/unknown-dwell reset edges contain a directed cycle.
    #[must_use]
    pub const fn has_zero_time_cycle(&self) -> bool {
        self.has_zero_time_cycle
    }

    /// Whether local dynamics, guards, resets, or event order are nonunique or
    /// unresolved independently of a zero-time cycle.
    #[must_use]
    pub const fn has_nonunique_local_semantics(&self) -> bool {
        self.has_nonunique_local_semantics
    }
}

/// Explicit authority boundary carried by every validated descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZenoScientificAuthorityV1 {
    /// Schema consistency and identity do not prove the claimed dynamics.
    ScientificCorrectnessNotProven,
}

/// Sealed, problem-bound, canonically identified claim descriptor.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedZenoClaimDescriptorV1 {
    draft: ZenoClaimDraftV1,
    receipt: IdentityReceipt<ZenoClaimIdV1>,
}

impl ValidatedZenoClaimDescriptorV1 {
    /// Observational descriptor view.
    #[must_use]
    pub const fn draft(&self) -> &ZenoClaimDraftV1 {
        &self.draft
    }

    /// Typed descriptor identity.
    #[must_use]
    pub const fn claim_id(&self) -> ZenoClaimIdV1 {
        self.receipt.id()
    }

    /// Exact identity and canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<ZenoClaimIdV1> {
        self.receipt
    }

    /// Permanent v1 no-claim boundary. RE.Z2 owns theorem promotion.
    #[must_use]
    pub const fn scientific_authority(&self) -> ZenoScientificAuthorityV1 {
        ZenoScientificAuthorityV1::ScientificCorrectnessNotProven
    }
}

/// Validate and canonically identify a finite hybrid/Zeno problem.
///
/// Collection caps and cancellation are checked before sorting, graph
/// closure, or identity publication.
///
/// # Errors
/// Returns [`HybridSemanticReportV1`] for invalid dimensions, time/budget
/// values, references, event order, regularization lineage, continuation
/// overclaims, cancellation, or canonical identity failure.
#[allow(clippy::too_many_lines)] // One ordered cross-field matrix keeps refusal semantics auditable.
#[must_use = "the hybrid/Zeno admission result must be handled before use"]
pub fn validate_zeno_problem_v1(
    mut ir: ZenoProblemIrV1,
    cx: &Cx<'_>,
) -> Result<ValidatedZenoProblemV1, HybridSemanticReportV1> {
    checkpoint(cx)?;
    let mut issues = Vec::new();
    if ir.modes.len() > MAX_HYBRID_MODES_V1 {
        issues.push(HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::Modes,
            found: ir.modes.len(),
            limit: MAX_HYBRID_MODES_V1,
        });
    }
    if ir.events.len() > MAX_HYBRID_EVENTS_V1 {
        issues.push(HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::Events,
            found: ir.events.len(),
            limit: MAX_HYBRID_EVENTS_V1,
        });
    }
    if let SimultaneousEventPolicyV1::TotalPriority { ordered_events, .. } = &ir.simultaneous_policy
        && ordered_events.len() > MAX_HYBRID_EVENTS_V1
    {
        issues.push(HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::PriorityEvents,
            found: ordered_events.len(),
            limit: MAX_HYBRID_EVENTS_V1,
        });
    }
    if !issues.is_empty() {
        return Err(HybridSemanticReportV1::new(issues));
    }

    let mut max_targets = 0_usize;
    let mut total_targets = 0_usize;
    for (index, event) in ir.events.iter().enumerate() {
        if index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        if let ResetSemanticsV1::SetValued { targets, .. } = &event.reset {
            max_targets = max_targets.max(targets.len());
            total_targets = total_targets.saturating_add(targets.len());
        }
    }
    if max_targets > MAX_RESET_TARGETS_V1 {
        issues.push(HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::ResetTargets,
            found: max_targets,
            limit: MAX_RESET_TARGETS_V1,
        });
    }
    if total_targets > MAX_TOTAL_RESET_TARGETS_V1 {
        issues.push(HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::TotalResetTargets,
            found: total_targets,
            limit: MAX_TOTAL_RESET_TARGETS_V1,
        });
    }
    if !issues.is_empty() {
        return Err(HybridSemanticReportV1::new(issues));
    }

    checkpoint(cx)?;
    canonicalize_problem_zeros(&mut ir);
    if ir.schema_version != ZENO_PROBLEM_SCHEMA_VERSION_V1 {
        issues.push(HybridSemanticIssueV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: ZENO_PROBLEM_SCHEMA_VERSION_V1,
        });
    }
    if ir.modes.is_empty() {
        issues.push(HybridSemanticIssueV1::EmptyCollection {
            collection: HybridCollectionV1::Modes,
        });
    }
    if ir.events.is_empty() {
        issues.push(HybridSemanticIssueV1::EmptyCollection {
            collection: HybridCollectionV1::Events,
        });
    }
    if !(ir.time_scale.seconds_per_unit.is_finite() && ir.time_scale.seconds_per_unit > 0.0) {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::SecondsPerUnit,
        });
    }
    if !ir.time_domain.start.is_finite() {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::TimeStart,
        });
    }
    if let HybridTimeEndV1::Finite(end) = ir.time_domain.end
        && (!end.is_finite() || end <= ir.time_domain.start)
    {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::TimeEnd,
        });
    }
    if matches!(ir.time_domain.event_cap, Some(0)) {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::EventCap,
        });
    }
    if ir.budget.max_event_word_len == 0 {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::EventWordBudget,
        });
    }
    if ir.budget.max_transitions == 0 {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::TransitionBudget,
        });
    }
    if !(ir.budget.max_wall_seconds.is_finite() && ir.budget.max_wall_seconds > 0.0) {
        issues.push(HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::WallTimeBudget,
        });
    }
    if let EventLanguageSemanticsV1::FiniteWords {
        max_events_per_word,
    } = ir.event_language.semantics
    {
        if max_events_per_word == 0 {
            issues.push(HybridSemanticIssueV1::InvalidValue {
                field: HybridFieldV1::EventWordBudget,
            });
        } else if max_events_per_word > ir.budget.max_event_word_len {
            issues.push(HybridSemanticIssueV1::EventLanguageBudgetMismatch);
        }
    }
    if let HybridModelLineageV1::Regularized {
        source_model,
        source_version,
        ..
    } = ir.lineage
        && source_model == ir.model
        && source_version == ir.model_version
    {
        issues.push(HybridSemanticIssueV1::RegularizationSelfReference);
    }

    ir.modes.sort_by_key(|mode| mode.mode);
    if ir.modes.windows(2).any(|pair| pair[0].mode == pair[1].mode) {
        issues.push(HybridSemanticIssueV1::DuplicateId {
            kind: HybridIdKindV1::Mode,
        });
    }
    for (index, mode) in ir.modes.iter().enumerate() {
        if index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        if mode.state_dimension == 0 {
            push_once(
                &mut issues,
                HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::StateDimension,
                },
            );
        }
        match mode.class {
            ContinuousDynamicsClassV1::AdmittedDae { index: 0, .. } => push_once(
                &mut issues,
                HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::DaeIndex,
                },
            ),
            ContinuousDynamicsClassV1::UnsupportedInfiniteDimensional => {
                push_once(&mut issues, HybridSemanticIssueV1::UnsupportedDynamicsClass)
            }
            ContinuousDynamicsClassV1::DeterministicOde
            | ContinuousDynamicsClassV1::DifferentialInclusion
            | ContinuousDynamicsClassV1::AdmittedDae { .. } => {}
        }
    }

    ir.events.sort_by_key(|event| event.event);
    if ir
        .events
        .windows(2)
        .any(|pair| pair[0].event == pair[1].event)
    {
        issues.push(HybridSemanticIssueV1::DuplicateId {
            kind: HybridIdKindV1::Event,
        });
    }
    for (index, event) in ir.events.iter_mut().enumerate() {
        if index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        if let DwellSemanticsV1::PositiveLowerBound { value, .. } = event.dwell
            && !(value.is_finite() && value > 0.0)
        {
            push_once(
                &mut issues,
                HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::DwellLowerBound,
                },
            );
        }
        if let ResetSemanticsV1::SetValued { targets, .. } = &mut event.reset {
            if targets.is_empty() {
                push_once(
                    &mut issues,
                    HybridSemanticIssueV1::EmptyCollection {
                        collection: HybridCollectionV1::ResetTargets,
                    },
                );
            }
            targets.sort_unstable();
            if targets.windows(2).any(|pair| pair[0] == pair[1]) {
                push_once(
                    &mut issues,
                    HybridSemanticIssueV1::DuplicateId {
                        kind: HybridIdKindV1::ResetTarget,
                    },
                );
            }
        }
    }

    let mode_exists = |mode: HybridModeIdV1| {
        ir.modes
            .binary_search_by_key(&mode, |candidate| candidate.mode)
            .is_ok()
    };
    for (index, event) in ir.events.iter().enumerate() {
        if index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        if !mode_exists(event.source_mode) {
            push_once(&mut issues, HybridSemanticIssueV1::UnknownModeReference);
        }
        match &event.reset {
            ResetSemanticsV1::Deterministic { target, .. } => {
                if !mode_exists(*target) {
                    push_once(&mut issues, HybridSemanticIssueV1::UnknownModeReference);
                }
            }
            ResetSemanticsV1::SetValued { targets, .. } => {
                if targets.iter().any(|target| !mode_exists(*target)) {
                    push_once(&mut issues, HybridSemanticIssueV1::UnknownModeReference);
                }
            }
            ResetSemanticsV1::Terminal { .. } | ResetSemanticsV1::Unknown { .. } => {}
        }
    }

    validate_accumulation_candidate(&ir, &mut issues);
    let (grouped_events, has_simultaneous_groups) = validate_simultaneous_policy(&ir, &mut issues);
    if !issues.is_empty() {
        return Err(HybridSemanticReportV1::new(issues));
    }

    let mut has_nonunique_local_semantics = ir
        .modes
        .iter()
        .any(|mode| matches!(mode.class, ContinuousDynamicsClassV1::DifferentialInclusion));
    has_nonunique_local_semantics |= matches!(
        ir.event_language.semantics,
        EventLanguageSemanticsV1::Unknown { .. }
    );
    has_nonunique_local_semantics |= ir.events.iter().any(|event| {
        matches!(
            event.orientation,
            GuardOrientationV1::Bidirectional | GuardOrientationV1::Unknown { .. }
        ) || matches!(
            event.crossing,
            CrossingSemanticsV1::Grazing { .. } | CrossingSemanticsV1::Unknown { .. }
        ) || matches!(
            &event.reset,
            ResetSemanticsV1::SetValued { .. } | ResetSemanticsV1::Unknown { .. }
        )
    });
    if has_simultaneous_groups
        && matches!(
            &ir.simultaneous_policy,
            SimultaneousEventPolicyV1::SetValued { .. } | SimultaneousEventPolicyV1::Unknown { .. }
        )
    {
        has_nonunique_local_semantics = true;
    }
    debug_assert_eq!(has_simultaneous_groups, !grouped_events.is_empty());
    let has_zero_time_cycle = zero_time_graph_has_cycle(&ir, cx)?;
    if matches!(ir.continuation, ContinuationCategoryV1::Unique { .. })
        && has_nonunique_local_semantics
    {
        issues.push(HybridSemanticIssueV1::UniqueContinuationUnsupported);
    }
    if !issues.is_empty() {
        return Err(HybridSemanticReportV1::new(issues));
    }

    let receipt = problem_receipt(&ir, cx).map_err(identity_report)?;
    Ok(ValidatedZenoProblemV1 {
        ir,
        receipt,
        has_zero_time_cycle,
        has_nonunique_local_semantics,
    })
}

/// Validate a problem-bound classification descriptor.
///
/// `regularized` is required only for [`ZenoClassificationV1::RegularizedModel`]
/// and must carry the exact source lineage. The returned value deliberately
/// retains [`ZenoScientificAuthorityV1::ScientificCorrectnessNotProven`].
///
/// # Errors
/// Returns [`HybridSemanticReportV1`] for schema/target mismatch, invalid
/// intervals, event-cap theorem overclaim, post-state/category mismatch,
/// invalid regularization lineage, cancellation, or identity failure.
#[must_use = "the claim-descriptor admission result must be handled before retention"]
#[allow(clippy::too_many_lines)] // One classification/post-state matrix keeps all refusals explicit.
pub fn validate_zeno_claim_descriptor_v1(
    mut draft: ZenoClaimDraftV1,
    problem: &ValidatedZenoProblemV1,
    regularized: Option<&ValidatedZenoProblemV1>,
    cx: &Cx<'_>,
) -> Result<ValidatedZenoClaimDescriptorV1, HybridSemanticReportV1> {
    checkpoint(cx)?;
    canonicalize_claim_zeros(&mut draft);
    let mut issues = Vec::new();
    if draft.schema_version != ZENO_CLAIM_SCHEMA_VERSION_V1 {
        issues.push(HybridSemanticIssueV1::UnsupportedSchemaVersion {
            found: draft.schema_version,
            supported: ZENO_CLAIM_SCHEMA_VERSION_V1,
        });
    }
    if draft.problem != problem.problem_id() {
        issues.push(HybridSemanticIssueV1::TargetProblemMismatch);
    }

    match draft.classification {
        ZenoClassificationV1::FiniteEventSeparation {
            minimum_separation,
            evidence,
        } => {
            if !(minimum_separation.is_finite() && minimum_separation > 0.0) {
                issues.push(HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::FiniteEventSeparation,
                });
            }
            require_theorem_evidence(evidence, &mut issues);
            if problem
                .ir
                .events
                .iter()
                .any(|event| !matches!(event.dwell, DwellSemanticsV1::PositiveLowerBound { .. }))
            {
                issues.push(HybridSemanticIssueV1::FiniteSeparationContradictsEventGraph);
            }
            if !matches!(draft.post_zeno, PostZenoStateV1::NotApplicable { .. }) {
                issues.push(HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        ZenoClassificationV1::CertifiedZeno {
            interval, evidence, ..
        } => {
            if !valid_interval(interval) {
                issues.push(HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::ZenoInterval,
                });
            } else if !interval_in_domain(interval, problem.ir.time_domain) {
                issues.push(HybridSemanticIssueV1::ClaimOutsideTimeDomain);
            }
            require_theorem_evidence(evidence, &mut issues);
            if !problem.has_zero_time_cycle {
                issues.push(HybridSemanticIssueV1::ZenoAccumulationCycleRequired);
            }
            if matches!(draft.post_zeno, PostZenoStateV1::NotApplicable { .. }) {
                issues.push(HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        ZenoClassificationV1::NumericalEventDensityWarning {
            observed_events,
            window,
            ..
        } => {
            if observed_events == 0 {
                issues.push(HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::ObservedEvents,
                });
            }
            if !(window.is_finite() && window > 0.0) {
                issues.push(HybridSemanticIssueV1::InvalidValue {
                    field: HybridFieldV1::ObservationWindow,
                });
            }
            if !matches!(draft.post_zeno, PostZenoStateV1::Unresolved { .. }) {
                issues.push(HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        ZenoClassificationV1::RegularizedModel {
            regularized_problem,
            regularization,
            no_equivalence,
        } => {
            let lineage_matches = regularized.is_some_and(|candidate| {
                candidate.problem_id() == regularized_problem
                    && matches!(
                        candidate.ir.lineage,
                        HybridModelLineageV1::Regularized {
                            source_model,
                            source_version,
                            regularization: candidate_regularization,
                            no_equivalence: candidate_no_equivalence,
                        } if source_model == problem.ir.model
                            && source_version == problem.ir.model_version
                            && candidate_regularization == regularization
                            && candidate_no_equivalence == no_equivalence
                    )
            });
            if !lineage_matches || regularized_problem == problem.problem_id() {
                issues.push(HybridSemanticIssueV1::RegularizationLineageMismatch);
            }
            if !matches!(draft.post_zeno, PostZenoStateV1::Unresolved { .. }) {
                issues.push(HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        ZenoClassificationV1::Unknown { .. } => {
            if !matches!(draft.post_zeno, PostZenoStateV1::Unresolved { .. }) {
                issues.push(HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
    }

    validate_post_zeno(draft.post_zeno, problem, &mut issues);
    if !issues.is_empty() {
        return Err(HybridSemanticReportV1::new(issues));
    }
    let receipt = claim_receipt(&draft, cx).map_err(identity_report)?;
    Ok(ValidatedZenoClaimDescriptorV1 { draft, receipt })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), HybridSemanticReportV1> {
    cx.checkpoint()
        .map_err(|_| HybridSemanticReportV1::new(vec![HybridSemanticIssueV1::Cancelled]))
}

fn identity_report(error: CanonicalError) -> HybridSemanticReportV1 {
    let issue = if matches!(&error, CanonicalError::Cancelled { .. }) {
        HybridSemanticIssueV1::Cancelled
    } else {
        HybridSemanticIssueV1::Identity(error)
    };
    HybridSemanticReportV1::new(vec![issue])
}

fn push_once(issues: &mut Vec<HybridSemanticIssueV1>, issue: HybridSemanticIssueV1) {
    if !issues.contains(&issue) {
        issues.push(issue);
    }
}

fn canonicalize_zero(value: &mut f64) {
    if *value == 0.0 {
        *value = 0.0;
    }
}

fn canonicalize_problem_zeros(ir: &mut ZenoProblemIrV1) {
    canonicalize_zero(&mut ir.time_scale.seconds_per_unit);
    canonicalize_zero(&mut ir.time_domain.start);
    if let HybridTimeEndV1::Finite(end) = &mut ir.time_domain.end {
        canonicalize_zero(end);
    }
    canonicalize_zero(&mut ir.budget.max_wall_seconds);
    if let AccumulationCandidateV1::Window {
        earliest, latest, ..
    } = &mut ir.accumulation_candidate
    {
        canonicalize_zero(earliest);
        canonicalize_zero(latest);
    }
    for event in &mut ir.events {
        if let DwellSemanticsV1::PositiveLowerBound { value, .. } = &mut event.dwell {
            canonicalize_zero(value);
        }
    }
}

fn canonicalize_claim_zeros(draft: &mut ZenoClaimDraftV1) {
    match &mut draft.classification {
        ZenoClassificationV1::FiniteEventSeparation {
            minimum_separation, ..
        } => canonicalize_zero(minimum_separation),
        ZenoClassificationV1::CertifiedZeno { interval, .. } => {
            canonicalize_zero(&mut interval.earliest);
            canonicalize_zero(&mut interval.latest);
        }
        ZenoClassificationV1::NumericalEventDensityWarning { window, .. } => {
            canonicalize_zero(window);
        }
        ZenoClassificationV1::RegularizedModel { .. } | ZenoClassificationV1::Unknown { .. } => {}
    }
}

fn valid_interval(interval: HybridTimeIntervalV1) -> bool {
    interval.earliest.is_finite()
        && interval.latest.is_finite()
        && interval.earliest <= interval.latest
}

fn interval_in_domain(interval: HybridTimeIntervalV1, domain: HybridTimeDomainV1) -> bool {
    interval.earliest >= domain.start
        && match domain.end {
            HybridTimeEndV1::Finite(end) => interval.latest <= end,
            HybridTimeEndV1::Infinite => true,
        }
}

fn validate_accumulation_candidate(ir: &ZenoProblemIrV1, issues: &mut Vec<HybridSemanticIssueV1>) {
    if let AccumulationCandidateV1::Window {
        earliest, latest, ..
    } = ir.accumulation_candidate
    {
        let interval = HybridTimeIntervalV1 { earliest, latest };
        if !valid_interval(interval) {
            issues.push(HybridSemanticIssueV1::InvalidValue {
                field: HybridFieldV1::AccumulationWindow,
            });
        } else if !interval_in_domain(interval, ir.time_domain) {
            issues.push(HybridSemanticIssueV1::AccumulationOutsideTimeDomain);
        }
    }
}

fn validate_simultaneous_policy(
    ir: &ZenoProblemIrV1,
    issues: &mut Vec<HybridSemanticIssueV1>,
) -> (Vec<HybridEventIdV1>, bool) {
    let mut groups: BTreeMap<SimultaneityGroupIdV1, Vec<HybridEventIdV1>> = BTreeMap::new();
    for event in &ir.events {
        if let EventSimultaneityV1::Group { group } = event.simultaneity {
            groups.entry(group).or_default().push(event.event);
        }
    }
    let mut grouped_events = Vec::new();
    for events in groups.values_mut() {
        events.sort_unstable();
        if events.len() < 2 {
            push_once(issues, HybridSemanticIssueV1::SingletonSimultaneityGroup);
        }
        grouped_events.extend_from_slice(events);
    }
    grouped_events.sort_unstable();
    let has_groups = !grouped_events.is_empty();
    match &ir.simultaneous_policy {
        SimultaneousEventPolicyV1::NoSimultaneousEvents { .. } if has_groups => {
            issues.push(HybridSemanticIssueV1::SimultaneousPolicyMismatch);
        }
        SimultaneousEventPolicyV1::TotalPriority { ordered_events, .. } => {
            let mut canonical = ordered_events.clone();
            canonical.sort_unstable();
            if canonical.windows(2).any(|pair| pair[0] == pair[1]) {
                push_once(
                    issues,
                    HybridSemanticIssueV1::DuplicateId {
                        kind: HybridIdKindV1::PriorityEvent,
                    },
                );
            }
            if canonical != grouped_events {
                issues.push(HybridSemanticIssueV1::InvalidPriorityOrder);
            }
        }
        SimultaneousEventPolicyV1::NoSimultaneousEvents { .. }
        | SimultaneousEventPolicyV1::Commuting { .. }
        | SimultaneousEventPolicyV1::SetValued { .. }
        | SimultaneousEventPolicyV1::Unknown { .. } => {}
    }
    (grouped_events, has_groups)
}

fn zero_time_graph_has_cycle(
    ir: &ZenoProblemIrV1,
    cx: &Cx<'_>,
) -> Result<bool, HybridSemanticReportV1> {
    let n = ir.modes.len();
    let mut adjacency = vec![false; n.saturating_mul(n)];
    for (event_index, event) in ir.events.iter().enumerate() {
        if event_index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        if matches!(event.dwell, DwellSemanticsV1::PositiveLowerBound { .. }) {
            continue;
        }
        let source = ir
            .modes
            .binary_search_by_key(&event.source_mode, |mode| mode.mode)
            .expect("mode references were validated before graph construction");
        match &event.reset {
            ResetSemanticsV1::Deterministic { target, .. } => {
                let target = ir
                    .modes
                    .binary_search_by_key(target, |mode| mode.mode)
                    .expect("reset targets were validated before graph construction");
                adjacency[source * n + target] = true;
            }
            ResetSemanticsV1::SetValued { targets, .. } => {
                for target in targets {
                    let target = ir
                        .modes
                        .binary_search_by_key(target, |mode| mode.mode)
                        .expect("reset targets were validated before graph construction");
                    adjacency[source * n + target] = true;
                }
            }
            ResetSemanticsV1::Terminal { .. } | ResetSemanticsV1::Unknown { .. } => {}
        }
    }

    let mut indegree = vec![0_usize; n];
    for source in 0..n {
        for target in 0..n {
            if adjacency[source * n + target] {
                indegree[target] += 1;
            }
        }
    }
    let mut queue: Vec<usize> = indegree
        .iter()
        .enumerate()
        .filter_map(|(index, degree)| (*degree == 0).then_some(index))
        .collect();
    let mut cursor = 0_usize;
    let mut removed = 0_usize;
    while cursor < queue.len() {
        if cursor.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        let source = queue[cursor];
        cursor += 1;
        removed += 1;
        for target in 0..n {
            if adjacency[source * n + target] {
                indegree[target] -= 1;
                if indegree[target] == 0 {
                    queue.push(target);
                }
            }
        }
    }
    Ok(removed != n)
}

fn require_theorem_evidence(
    evidence: ZenoEvidenceReferenceV1,
    issues: &mut Vec<HybridSemanticIssueV1>,
) {
    if !matches!(
        evidence,
        ZenoEvidenceReferenceV1::Analytic { .. }
            | ZenoEvidenceReferenceV1::IntervalValidated { .. }
    ) {
        issues.push(HybridSemanticIssueV1::InsufficientTheoremEvidence);
    }
}

fn validate_post_zeno(
    post: PostZenoStateV1,
    problem: &ValidatedZenoProblemV1,
    issues: &mut Vec<HybridSemanticIssueV1>,
) {
    match post {
        PostZenoStateV1::NotApplicable { .. } | PostZenoStateV1::Unresolved { .. } => {}
        PostZenoStateV1::Unique { rule, witness, .. } => {
            if problem.has_nonunique_local_semantics
                || !matches!(
                    problem.ir.continuation,
                    ContinuationCategoryV1::Unique {
                        rule: declared_rule,
                        witness: declared_witness,
                    } if declared_rule == rule && declared_witness == witness
                )
            {
                push_once(issues, HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        PostZenoStateV1::SetValued { rule, .. } => {
            if !matches!(
                problem.ir.continuation,
                ContinuationCategoryV1::SetValued {
                    rule: declared_rule,
                } if declared_rule == rule
            ) {
                push_once(issues, HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
        PostZenoStateV1::Terminal { rule, .. } => {
            if !matches!(
                problem.ir.continuation,
                ContinuationCategoryV1::Terminal {
                    rule: declared_rule,
                } if declared_rule == rule
            ) {
                push_once(issues, HybridSemanticIssueV1::PostZenoSemanticsMismatch);
            }
        }
    }
}

fn push_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn canonical_f64_bits(value: f64) -> u64 {
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    push_u64(out, canonical_f64_bits(value));
}

fn push_digest<I: DigestBytes>(out: &mut Vec<u8>, id: I) {
    out.extend_from_slice(id.digest_bytes());
}

fn model_bytes(ir: &ZenoProblemIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_digest(&mut out, ir.model);
    push_digest(&mut out, ir.model_version);
    match ir.lineage {
        HybridModelLineageV1::Original => out.push(0),
        HybridModelLineageV1::Regularized {
            source_model,
            source_version,
            regularization,
            no_equivalence,
        } => {
            out.push(1);
            push_digest(&mut out, source_model);
            push_digest(&mut out, source_version);
            push_digest(&mut out, regularization);
            push_digest(&mut out, no_equivalence);
        }
    }
    out
}

fn physical_context_bytes(ir: &ZenoProblemIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(160);
    push_digest(&mut out, ir.frame);
    push_digest(&mut out, ir.state_units);
    push_digest(&mut out, ir.initial_states);
    push_digest(&mut out, ir.time_scale.unit);
    push_f64(&mut out, ir.time_scale.seconds_per_unit);
    out
}

fn mode_bytes(mode: &HybridModeSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(128);
    push_digest(&mut out, mode.mode);
    push_digest(&mut out, mode.dynamics);
    push_u32(&mut out, mode.state_dimension);
    match mode.class {
        ContinuousDynamicsClassV1::DeterministicOde => out.push(0),
        ContinuousDynamicsClassV1::DifferentialInclusion => out.push(1),
        ContinuousDynamicsClassV1::AdmittedDae { index, constraint } => {
            out.push(2);
            out.push(index);
            push_digest(&mut out, constraint);
        }
        ContinuousDynamicsClassV1::UnsupportedInfiniteDimensional => out.push(3),
    }
    out
}

fn push_reset(out: &mut Vec<u8>, reset: &ResetSemanticsV1) {
    match reset {
        ResetSemanticsV1::Deterministic { relation, target } => {
            out.push(0);
            push_digest(out, *relation);
            push_digest(out, *target);
        }
        ResetSemanticsV1::SetValued {
            relation,
            targets,
            states,
        } => {
            out.push(1);
            push_digest(out, *relation);
            push_u64(out, targets.len() as u64);
            for target in targets {
                push_digest(out, *target);
            }
            push_digest(out, *states);
        }
        ResetSemanticsV1::Terminal { relation } => {
            out.push(2);
            push_digest(out, *relation);
        }
        ResetSemanticsV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(out, *no_claim);
        }
    }
}

fn event_bytes(event: &HybridEventSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(384);
    push_digest(&mut out, event.event);
    push_digest(&mut out, event.source_mode);
    push_digest(&mut out, event.guard);
    match event.orientation {
        GuardOrientationV1::NegativeToPositive => out.push(0),
        GuardOrientationV1::PositiveToNegative => out.push(1),
        GuardOrientationV1::Bidirectional => out.push(2),
        GuardOrientationV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(&mut out, no_claim);
        }
    }
    match event.crossing {
        CrossingSemanticsV1::Transverse { witness } => {
            out.push(0);
            push_digest(&mut out, witness);
        }
        CrossingSemanticsV1::Grazing { witness } => {
            out.push(1);
            push_digest(&mut out, witness);
        }
        CrossingSemanticsV1::Unknown { no_claim } => {
            out.push(2);
            push_digest(&mut out, no_claim);
        }
    }
    push_reset(&mut out, &event.reset);
    match event.law {
        InteractionLawV1::None { justification } => {
            out.push(0);
            push_digest(&mut out, justification);
        }
        InteractionLawV1::Contact { law } => {
            out.push(1);
            push_digest(&mut out, law);
        }
        InteractionLawV1::Relay { law } => {
            out.push(2);
            push_digest(&mut out, law);
        }
        InteractionLawV1::Other { law } => {
            out.push(3);
            push_digest(&mut out, law);
        }
    }
    match event.simultaneity {
        EventSimultaneityV1::Exclusive { witness } => {
            out.push(0);
            push_digest(&mut out, witness);
        }
        EventSimultaneityV1::Group { group } => {
            out.push(1);
            push_digest(&mut out, group);
        }
    }
    match event.dwell {
        DwellSemanticsV1::PositiveLowerBound { value, witness } => {
            out.push(0);
            push_f64(&mut out, value);
            push_digest(&mut out, witness);
        }
        DwellSemanticsV1::ZeroAllowed => out.push(1),
        DwellSemanticsV1::Unknown { no_claim } => {
            out.push(2);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn event_language_bytes(language: EventLanguageSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(40);
    push_digest(&mut out, language.language);
    match language.semantics {
        EventLanguageSemanticsV1::FiniteWords {
            max_events_per_word,
        } => {
            out.push(0);
            push_u32(&mut out, max_events_per_word);
        }
        EventLanguageSemanticsV1::PrefixClosed => out.push(1),
        EventLanguageSemanticsV1::OmegaLanguage => out.push(2),
        EventLanguageSemanticsV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn simultaneous_policy_bytes(policy: &SimultaneousEventPolicyV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(128);
    match policy {
        SimultaneousEventPolicyV1::NoSimultaneousEvents { witness } => {
            out.push(0);
            push_digest(&mut out, *witness);
        }
        SimultaneousEventPolicyV1::TotalPriority {
            ordered_events,
            witness,
        } => {
            out.push(1);
            push_u64(&mut out, ordered_events.len() as u64);
            for event in ordered_events {
                push_digest(&mut out, *event);
            }
            push_digest(&mut out, *witness);
        }
        SimultaneousEventPolicyV1::Commuting { witness } => {
            out.push(2);
            push_digest(&mut out, *witness);
        }
        SimultaneousEventPolicyV1::SetValued { outcomes } => {
            out.push(3);
            push_digest(&mut out, *outcomes);
        }
        SimultaneousEventPolicyV1::Unknown { no_claim } => {
            out.push(4);
            push_digest(&mut out, *no_claim);
        }
    }
    out
}

fn hybrid_time_bytes(domain: HybridTimeDomainV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(32);
    push_f64(&mut out, domain.start);
    match domain.end {
        HybridTimeEndV1::Finite(end) => {
            out.push(0);
            push_f64(&mut out, end);
        }
        HybridTimeEndV1::Infinite => out.push(1),
    }
    match domain.event_cap {
        Some(cap) => {
            out.push(1);
            push_u64(&mut out, cap);
        }
        None => out.push(0),
    }
    out
}

fn compactness_bytes(compactness: CompactnessSemanticsV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(40);
    match compactness {
        CompactnessSemanticsV1::Compact { witness } => {
            out.push(0);
            push_digest(&mut out, witness);
        }
        CompactnessSemanticsV1::LocallyCompact { witness } => {
            out.push(1);
            push_digest(&mut out, witness);
        }
        CompactnessSemanticsV1::Unestablished { no_claim } => {
            out.push(2);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn accumulation_bytes(candidate: AccumulationCandidateV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(112);
    match candidate {
        AccumulationCandidateV1::None { no_claim } => {
            out.push(0);
            push_digest(&mut out, no_claim);
        }
        AccumulationCandidateV1::Window {
            earliest,
            latest,
            states,
            trace,
        } => {
            out.push(1);
            push_f64(&mut out, earliest);
            push_f64(&mut out, latest);
            push_digest(&mut out, states);
            push_digest(&mut out, trace);
        }
    }
    out
}

fn continuation_bytes(continuation: ContinuationCategoryV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(72);
    match continuation {
        ContinuationCategoryV1::Unique { rule, witness } => {
            out.push(0);
            push_digest(&mut out, rule);
            push_digest(&mut out, witness);
        }
        ContinuationCategoryV1::SetValued { rule } => {
            out.push(1);
            push_digest(&mut out, rule);
        }
        ContinuationCategoryV1::Terminal { rule } => {
            out.push(2);
            push_digest(&mut out, rule);
        }
        ContinuationCategoryV1::Unresolved { no_claim } => {
            out.push(3);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn budget_bytes(budget: HybridAnalysisBudgetV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(24);
    push_u32(&mut out, budget.max_event_word_len);
    push_u64(&mut out, budget.max_transitions);
    push_f64(&mut out, budget.max_wall_seconds);
    out
}

fn problem_receipt(
    ir: &ZenoProblemIrV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<ZenoProblemIdV1>, CanonicalError> {
    let model = model_bytes(ir);
    let physical_context = physical_context_bytes(ir);
    let modes: Vec<Vec<u8>> = ir.modes.iter().map(mode_bytes).collect();
    let events: Vec<Vec<u8>> = ir.events.iter().map(event_bytes).collect();
    let event_language = event_language_bytes(ir.event_language);
    let simultaneous_policy = simultaneous_policy_bytes(&ir.simultaneous_policy);
    let hybrid_time = hybrid_time_bytes(ir.time_domain);
    let compactness = compactness_bytes(ir.compactness);
    let accumulation = accumulation_bytes(ir.accumulation_candidate);
    let continuation = continuation_bytes(ir.continuation);
    let budget = budget_bytes(ir.budget);

    CanonicalEncoder::<ZenoProblemIdV1, _>::new(HYBRID_IDENTITY_LIMITS, || {
        cx.is_cancel_requested()
    })?
    .bytes(Field::new(0, "model"), &model)?
    .bytes(Field::new(1, "physical-context"), &physical_context)?
    .canonical_set(
        Field::new(2, "modes"),
        modes.len() as u64,
        modes.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(3, "events"),
        events.len() as u64,
        events.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(4, "event-language"), &event_language)?
    .bytes(Field::new(5, "simultaneous-policy"), &simultaneous_policy)?
    .bytes(Field::new(6, "hybrid-time"), &hybrid_time)?
    .bytes(Field::new(7, "compactness"), &compactness)?
    .bytes(Field::new(8, "accumulation-candidate"), &accumulation)?
    .bytes(Field::new(9, "continuation"), &continuation)?
    .bytes(Field::new(10, "analysis-budget"), &budget)?
    .finish()
}

fn push_evidence(out: &mut Vec<u8>, evidence: ZenoEvidenceReferenceV1) {
    match evidence {
        ZenoEvidenceReferenceV1::Analytic { witness } => {
            out.push(0);
            push_digest(out, witness);
        }
        ZenoEvidenceReferenceV1::IntervalValidated { witness } => {
            out.push(1);
            push_digest(out, witness);
        }
        ZenoEvidenceReferenceV1::EventCap { trace } => {
            out.push(2);
            push_digest(out, trace);
        }
        ZenoEvidenceReferenceV1::NumericalOnly { trace } => {
            out.push(3);
            push_digest(out, trace);
        }
        ZenoEvidenceReferenceV1::Unknown { no_claim } => {
            out.push(4);
            push_digest(out, no_claim);
        }
    }
}

fn classification_bytes(classification: ZenoClassificationV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    match classification {
        ZenoClassificationV1::FiniteEventSeparation {
            minimum_separation,
            evidence,
        } => {
            out.push(0);
            push_f64(&mut out, minimum_separation);
            push_evidence(&mut out, evidence);
        }
        ZenoClassificationV1::CertifiedZeno {
            interval,
            states,
            evidence,
        } => {
            out.push(1);
            push_f64(&mut out, interval.earliest);
            push_f64(&mut out, interval.latest);
            push_digest(&mut out, states);
            push_evidence(&mut out, evidence);
        }
        ZenoClassificationV1::NumericalEventDensityWarning {
            trace,
            observed_events,
            window,
        } => {
            out.push(2);
            push_digest(&mut out, trace);
            push_u64(&mut out, observed_events);
            push_f64(&mut out, window);
        }
        ZenoClassificationV1::RegularizedModel {
            regularized_problem,
            regularization,
            no_equivalence,
        } => {
            out.push(3);
            out.extend_from_slice(regularized_problem.as_bytes());
            push_digest(&mut out, regularization);
            push_digest(&mut out, no_equivalence);
        }
        ZenoClassificationV1::Unknown { no_claim } => {
            out.push(4);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn post_zeno_bytes(post: PostZenoStateV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(128);
    match post {
        PostZenoStateV1::NotApplicable { justification } => {
            out.push(0);
            push_digest(&mut out, justification);
        }
        PostZenoStateV1::Unique {
            state,
            rule,
            witness,
        } => {
            out.push(1);
            push_digest(&mut out, state);
            push_digest(&mut out, rule);
            push_digest(&mut out, witness);
        }
        PostZenoStateV1::SetValued {
            states,
            rule,
            witness,
        } => {
            out.push(2);
            push_digest(&mut out, states);
            push_digest(&mut out, rule);
            push_digest(&mut out, witness);
        }
        PostZenoStateV1::Terminal { rule, witness } => {
            out.push(3);
            push_digest(&mut out, rule);
            push_digest(&mut out, witness);
        }
        PostZenoStateV1::Unresolved { no_claim } => {
            out.push(4);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn claim_receipt(
    draft: &ZenoClaimDraftV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<ZenoClaimIdV1>, CanonicalError> {
    let classification = classification_bytes(draft.classification);
    let post_zeno = post_zeno_bytes(draft.post_zeno);
    CanonicalEncoder::<ZenoClaimIdV1, _>::new(HYBRID_IDENTITY_LIMITS, || cx.is_cancel_requested())?
        .bytes(Field::new(0, "problem"), draft.problem.as_bytes())?
        .bytes(Field::new(1, "classification"), &classification)?
        .bytes(Field::new(2, "post-zeno"), &post_zeno)?
        .finish()
}

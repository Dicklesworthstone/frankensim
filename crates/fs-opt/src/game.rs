//! Versioned quantified reach-avoid and viability game semantics (RE.Q1).
//!
//! This module admits finite, explicitly typed game descriptions. It seals
//! quantifier order, information timing, strategy representation, model and
//! set lineage, horizon/stopping semantics, composition, and proof polarity
//! into one deterministic identity. Admission is not a winning-set theorem.

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_exec::Cx;

/// Current quantified-game schema.
pub const GAME_PROBLEM_SCHEMA_VERSION_V1: u32 = 1;
/// Hard cap on quantified variables.
pub const MAX_GAME_QUANTIFIERS_V1: usize = 16;
/// Hard cap on information grants.
pub const MAX_INFORMATION_GRANTS_V1: usize = 64;
/// Hard cap on strategy descriptions.
pub const MAX_GAME_STRATEGIES_V1: usize = 8;
/// Hard cap on dependencies summed over every strategy.
pub const MAX_STRATEGY_DEPENDENCIES_V1: usize = 256;
/// Hard cap on composed problem references.
pub const MAX_GAME_COMPONENTS_V1: usize = 128;

const GAME_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 20, 1 << 20, 16, 8192, 8192);

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
}

macro_rules! opaque_game_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct from exact typed digest bytes.
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

opaque_game_id!(
    /// Content identity of the dynamical model.
    GameModelIdV1
);
opaque_game_id!(
    /// Immutable version identity of the dynamical model.
    GameModelVersionIdV1
);
opaque_game_id!(
    /// State-space identity.
    GameStateSpaceIdV1
);
opaque_game_id!(
    /// State-frame identity.
    GameFrameIdV1
);
opaque_game_id!(
    /// Unit-system identity.
    GameUnitSystemIdV1
);
opaque_game_id!(
    /// Time-unit identity.
    GameTimeUnitIdV1
);
opaque_game_id!(
    /// Initial-state-set identity.
    InitialSetIdV1
);
opaque_game_id!(
    /// Target-set identity.
    TargetSetIdV1
);
opaque_game_id!(
    /// Unsafe-set identity.
    UnsafeSetIdV1
);
opaque_game_id!(
    /// Control-set identity.
    ControlSetIdV1
);
opaque_game_id!(
    /// Disturbance-set identity.
    DisturbanceSetIdV1
);
opaque_game_id!(
    /// Uncertain-parameter-set identity.
    ParameterSetIdV1
);
opaque_game_id!(
    /// Hybrid-mode-set identity.
    GameModeSetIdV1
);
opaque_game_id!(
    /// Hybrid-event-set identity.
    GameEventSetIdV1
);
opaque_game_id!(
    /// DAE constraint-manifold identity.
    GameDaeConstraintIdV1
);
opaque_game_id!(
    /// Observation-map identity.
    GameObservationMapIdV1
);
opaque_game_id!(
    /// Strategy artifact identity.
    GameStrategyArtifactIdV1
);
opaque_game_id!(
    /// State-dependent stopping-rule identity.
    GameStoppingRuleIdV1
);
opaque_game_id!(
    /// Composition-interface identity.
    GameCompositionInterfaceIdV1
);
opaque_game_id!(
    /// Exact analytic or checker witness identity.
    GameWitnessIdV1
);
opaque_game_id!(
    /// Deterministic state/time-cell decomposition identity.
    GameCellDecompositionIdV1
);
opaque_game_id!(
    /// Explicit no-claim artifact identity.
    GameNoClaimIdV1
);

/// Domain-separated identity schema for an admitted game problem.
pub enum GameProblemIdentitySchemaV1 {}

impl CanonicalSchema for GameProblemIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-opt.quantified-game.v1";
    const NAME: &'static str = "quantified-reach-avoid-viability-game";
    const VERSION: u32 = GAME_PROBLEM_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "model, typed sets, objective and polarity, quantifier prefix, information, strategies, horizon, stopping, composition, and budgets";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("model", WireType::Bytes),
        FieldSpec::required("state-sets", WireType::Bytes),
        FieldSpec::required("action-sets", WireType::Bytes),
        FieldSpec::required("claim", WireType::Bytes),
        FieldSpec::required("quantifier-prefix", WireType::OrderedBytes),
        FieldSpec::required("information-pattern", WireType::CanonicalSet),
        FieldSpec::required("strategies", WireType::CanonicalSet),
        FieldSpec::required("horizon", WireType::Bytes),
        FieldSpec::required("stopping", WireType::Bytes),
        FieldSpec::required("composition", WireType::Bytes),
        FieldSpec::required("analysis-budget", WireType::Bytes),
    ];
}

/// Typed semantic identity of one admitted game problem.
pub type GameProblemIdV1 = ProblemSemanticId<GameProblemIdentitySchemaV1>;

/// Common context carried independently by every state-set role.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateSetContextV1 {
    /// State space.
    pub state_space: GameStateSpaceIdV1,
    /// State frame.
    pub frame: GameFrameIdV1,
    /// State units.
    pub units: GameUnitSystemIdV1,
    /// Model version under which the set is meaningful.
    pub model_version: GameModelVersionIdV1,
}

/// Typed initial-state set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InitialSetV1 {
    /// Initial-set artifact.
    pub set: InitialSetIdV1,
    /// Exact physical/model context.
    pub context: StateSetContextV1,
}

/// Typed reach target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TargetSetV1 {
    /// Target-set artifact.
    pub set: TargetSetIdV1,
    /// Exact physical/model context.
    pub context: StateSetContextV1,
}

/// Typed unsafe set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnsafeSetV1 {
    /// Unsafe-set artifact.
    pub set: UnsafeSetIdV1,
    /// Exact physical/model context.
    pub context: StateSetContextV1,
}

/// Typed control action set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ControlSetV1 {
    /// Control-set artifact.
    pub set: ControlSetIdV1,
    /// Control units.
    pub units: GameUnitSystemIdV1,
    /// Owning model version.
    pub model_version: GameModelVersionIdV1,
}

/// Typed disturbance signal/action set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DisturbanceSetV1 {
    /// Disturbance-set artifact.
    pub set: DisturbanceSetIdV1,
    /// Disturbance units.
    pub units: GameUnitSystemIdV1,
    /// Owning model version.
    pub model_version: GameModelVersionIdV1,
}

/// Typed uncertain-parameter set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParameterSetV1 {
    /// Parameter-set artifact.
    pub set: ParameterSetIdV1,
    /// Parameter units.
    pub units: GameUnitSystemIdV1,
    /// Owning model version.
    pub model_version: GameModelVersionIdV1,
}

/// Typed hybrid mode set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameModeSetV1 {
    /// Mode-set artifact.
    pub modes: GameModeSetIdV1,
    /// Owning model version.
    pub model_version: GameModelVersionIdV1,
}

/// Typed hybrid event set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameEventSetV1 {
    /// Event-set artifact.
    pub events: GameEventSetIdV1,
    /// Owning model version.
    pub model_version: GameModelVersionIdV1,
}

/// Hybrid execution scope relevant to game claims.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridZenoScopeV1 {
    /// A retained theorem excludes accumulation on the admitted horizon.
    Excluded {
        /// Exact exclusion witness.
        witness: GameWitnessIdV1,
    },
    /// Accumulation is admitted but unresolved. Reachability remains Unknown.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: GameNoClaimIdV1,
    },
    /// The game uses a distinct named regularized model.
    Regularized {
        /// Regularization artifact.
        regularization: GameWitnessIdV1,
        /// Explicit no-equivalence boundary.
        no_equivalence: GameNoClaimIdV1,
    },
}

/// Supported game dynamics classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameModelClassV1 {
    /// Deterministic controlled dynamics.
    Deterministic,
    /// Continuous differential game.
    DifferentialGame,
    /// Finite hybrid game.
    Hybrid {
        /// Exact mode set.
        modes: GameModeSetV1,
        /// Exact event set.
        events: GameEventSetV1,
        /// Zeno scope.
        zeno: HybridZenoScopeV1,
    },
    /// Explicit finite-dimensional DAE game.
    AdmittedDae {
        /// Positive differentiation index.
        index: u8,
        /// Constraint-manifold artifact.
        constraint: GameDaeConstraintIdV1,
    },
    /// Unsupported infinite-dimensional game, retained only for fail-closed decode.
    UnsupportedInfiniteDimensional,
}

/// Versioned physical model description.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameModelSpecV1 {
    /// Model identity.
    pub model: GameModelIdV1,
    /// Immutable version.
    pub version: GameModelVersionIdV1,
    /// State space.
    pub state_space: GameStateSpaceIdV1,
    /// State frame.
    pub frame: GameFrameIdV1,
    /// State units.
    pub state_units: GameUnitSystemIdV1,
    /// Dynamics class.
    pub class: GameModelClassV1,
}

/// Quantified game variable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GameVariableV1 {
    /// Initial state.
    InitialState,
    /// Controller action/signal.
    Control,
    /// Adversarial disturbance.
    Disturbance,
    /// Uncertain but fixed or evolving parameter.
    Parameter,
}

impl fmt::Display for GameVariableV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::InitialState => "initial-state",
            Self::Control => "control",
            Self::Disturbance => "disturbance",
            Self::Parameter => "parameter",
        };
        f.write_str(name)
    }
}

/// Quantifier kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameQuantifierV1 {
    /// Existential choice.
    Exists,
    /// Universal adversarial/environmental choice.
    ForAll,
}

impl fmt::Display for GameQuantifierV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Exists => "exists",
            Self::ForAll => "forall",
        })
    }
}

/// Typed domain named by one quantifier clause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameQuantifierDomainV1 {
    /// Initial-state domain.
    Initial(InitialSetIdV1),
    /// Control domain.
    Control(ControlSetIdV1),
    /// Disturbance domain.
    Disturbance(DisturbanceSetIdV1),
    /// Parameter domain.
    Parameter(ParameterSetIdV1),
}

/// One ordered quantifier clause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameQuantifierClauseV1 {
    /// Exists or forall.
    pub quantifier: GameQuantifierV1,
    /// Variable bound.
    pub variable: GameVariableV1,
    /// Exact typed domain.
    pub domain: GameQuantifierDomainV1,
}

/// Exact, order-sensitive game quantifier prefix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GameQuantifierPrefixV1 {
    clauses: Vec<GameQuantifierClauseV1>,
}

impl GameQuantifierPrefixV1 {
    /// Construct an order-sensitive prefix.
    #[must_use]
    pub fn new(clauses: Vec<GameQuantifierClauseV1>) -> Self {
        Self { clauses }
    }

    /// Ordered clauses.
    #[must_use]
    pub fn clauses(&self) -> &[GameQuantifierClauseV1] {
        &self.clauses
    }
}

impl fmt::Display for GameQuantifierPrefixV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, clause) in self.clauses.iter().enumerate() {
            if index != 0 {
                f.write_str("; ")?;
            }
            write!(f, "{} {} in ", clause.quantifier, clause.variable)?;
            write_domain(f, clause.domain)?;
        }
        Ok(())
    }
}

/// Player whose information and strategy are described.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GamePlayerV1 {
    /// Controller/protagonist.
    Controller,
    /// Disturbance/adversary.
    Disturbance,
}

impl fmt::Display for GamePlayerV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Controller => "controller",
            Self::Disturbance => "disturbance-player",
        })
    }
}

/// Observable subject.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum InformationSubjectV1 {
    /// Continuous state.
    State,
    /// Hybrid mode.
    Mode,
    /// Hybrid event.
    Event,
    /// Control signal.
    Control,
    /// Disturbance signal.
    Disturbance,
    /// Uncertain parameter.
    Parameter,
}

/// When a subject becomes observable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ObservationAvailabilityV1 {
    /// Subject is hidden.
    Hidden,
    /// Only its initial value is available.
    InitialOnly,
    /// Current value is available.
    Current,
    /// Full history through the current time is available.
    HistoryThroughCurrent,
    /// Value arrives after a nonnegative physical lag.
    Delayed {
        /// Lag in the declared horizon time unit.
        lag: f64,
    },
    /// Future values are exposed; admission refuses this anticipative grant.
    FutureTrajectory,
}

/// One exact information grant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InformationGrantV1 {
    /// Recipient.
    pub player: GamePlayerV1,
    /// Subject.
    pub subject: InformationSubjectV1,
    /// Observation map.
    pub observation: GameObservationMapIdV1,
    /// Timing semantics.
    pub availability: ObservationAvailabilityV1,
}

/// Complete nonanticipative information pattern.
#[derive(Debug, Clone, PartialEq)]
pub struct InformationPatternV1 {
    /// Grants, canonicalized as a set by admission.
    pub grants: Vec<InformationGrantV1>,
}

/// How a strategy reads one subject.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StrategyAccessV1 {
    /// Initial value only.
    InitialOnly,
    /// Current value.
    Current,
    /// History through current time.
    HistoryThroughCurrent,
    /// Delayed value/history.
    Delayed {
        /// Strategy's declared lag.
        lag: f64,
    },
    /// Any future value; admission always refuses.
    FutureTrajectory,
}

/// One strategy dependency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrategyDependencyV1 {
    /// Observed subject.
    pub subject: InformationSubjectV1,
    /// Exact temporal access.
    pub access: StrategyAccessV1,
}

/// Strategy representation and causal class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrategyRepresentationV1 {
    /// One open-loop signal chosen without online observations.
    OpenLoop {
        /// Exact signal/policy artifact.
        artifact: GameStrategyArtifactIdV1,
    },
    /// State-feedback policy.
    StateFeedback {
        /// Exact policy artifact.
        artifact: GameStrategyArtifactIdV1,
    },
    /// Explicit nonanticipative history policy.
    NonanticipativeFeedback {
        /// Exact policy artifact.
        artifact: GameStrategyArtifactIdV1,
        /// Finite memory-state bound.
        memory_states: u32,
    },
    /// Hybrid-mode-aware feedback policy.
    HybridModeFeedback {
        /// Exact policy artifact.
        artifact: GameStrategyArtifactIdV1,
    },
    /// Set-valued strategy relation.
    SetValued {
        /// Exact relation artifact.
        artifact: GameStrategyArtifactIdV1,
    },
    /// Strategy class is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: GameNoClaimIdV1,
    },
}

impl fmt::Display for StrategyRepresentationV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::OpenLoop { .. } => "open-loop",
            Self::StateFeedback { .. } => "state-feedback",
            Self::NonanticipativeFeedback { .. } => "nonanticipative-feedback",
            Self::HybridModeFeedback { .. } => "hybrid-mode-feedback",
            Self::SetValued { .. } => "set-valued",
            Self::Unknown { .. } => "unknown",
        })
    }
}

/// Strategy for one player.
#[derive(Debug, Clone, PartialEq)]
pub struct GameStrategySpecV1 {
    /// Player implementing the strategy.
    pub player: GamePlayerV1,
    /// Strategy class.
    pub representation: StrategyRepresentationV1,
    /// Observation dependencies, canonicalized as a set.
    pub dependencies: Vec<StrategyDependencyV1>,
}

/// Game objective.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameObjectiveV1 {
    /// Reach target before entering unsafe set.
    ReachAvoid,
    /// Remain outside unsafe set throughout the admitted horizon.
    Viability,
}

/// Requested proof-set polarity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameProofPolarityV1 {
    /// Certified subset of the true winning/viability set.
    Inner,
    /// Certified superset of the true losing/nonviable complement boundary.
    Outer,
    /// Exact equality is requested.
    Exact,
    /// No polarity claim is available.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: GameNoClaimIdV1,
    },
}

/// Objective and requested proof polarity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameClaimSpecV1 {
    /// Reach-avoid or viability.
    pub objective: GameObjectiveV1,
    /// Inner, outer, exact, or Unknown.
    pub polarity: GameProofPolarityV1,
}

/// Physical horizon.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GameHorizonV1 {
    /// Closed finite time interval.
    Finite {
        /// Inclusive start.
        start: f64,
        /// Inclusive end.
        end: f64,
        /// Physical time unit.
        unit: GameTimeUnitIdV1,
        /// SI seconds per unit.
        seconds_per_unit: f64,
    },
    /// Infinite horizon. V1 admits the object but exposes only Unknown.
    Infinite {
        /// Inclusive start.
        start: f64,
        /// Physical time unit.
        unit: GameTimeUnitIdV1,
        /// SI seconds per unit.
        seconds_per_unit: f64,
        /// Explicit no-claim artifact.
        no_claim: GameNoClaimIdV1,
    },
}

/// Stopping semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameStoppingSemanticsV1 {
    /// Stop only at the finite horizon.
    FixedHorizon,
    /// Stop on first target hit.
    FirstTarget,
    /// Stop on first target or unsafe hit.
    FirstTargetOrUnsafe,
    /// Exact state-dependent stopping rule.
    StateDependent {
        /// Rule artifact.
        rule: GameStoppingRuleIdV1,
    },
    /// Stop on an admitted hybrid terminal event.
    HybridTerminalEvent {
        /// Event-set containing terminal events.
        events: GameEventSetIdV1,
    },
    /// Never stop; valid only for an infinite horizon and still Unknown in v1.
    Never,
}

/// Context carried by a referenced composed game.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GameComponentRefV1 {
    /// Component problem identity.
    pub problem: GameProblemIdV1,
    /// Component model version.
    pub model_version: GameModelVersionIdV1,
    /// State space.
    pub state_space: GameStateSpaceIdV1,
    /// Frame.
    pub frame: GameFrameIdV1,
    /// State units.
    pub units: GameUnitSystemIdV1,
}

/// Composition semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GameCompositionV1 {
    /// Atomic game.
    Atomic,
    /// Ordered sequential composition.
    Sequential {
        /// Ordered components.
        components: Vec<GameComponentRefV1>,
        /// Exact handoff/interface.
        interface: GameCompositionInterfaceIdV1,
    },
    /// Order-independent parallel product.
    Parallel {
        /// Canonical component set.
        components: Vec<GameComponentRefV1>,
        /// Exact coupling/interface.
        interface: GameCompositionInterfaceIdV1,
    },
    /// Order-independent adversarial product.
    AdversarialProduct {
        /// Canonical component set.
        components: Vec<GameComponentRefV1>,
        /// Exact game interface.
        interface: GameCompositionInterfaceIdV1,
    },
}

/// Explicit bounded analysis budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GameAnalysisBudgetV1 {
    /// Exact deterministic cell decomposition.
    pub cell_decomposition: GameCellDecompositionIdV1,
    /// Maximum deterministic cells.
    pub max_cells: u64,
    /// Maximum transition/event expansions.
    pub max_transitions: u64,
    /// Maximum retained strategy nodes.
    pub max_strategy_nodes: u64,
    /// Positive wall-time allowance in seconds.
    pub max_wall_seconds: f64,
}

/// Raw versioned game problem. It has no authority until admitted.
#[derive(Debug, Clone, PartialEq)]
pub struct GameProblemIrV1 {
    schema_version: u32,
    /// Physical model.
    pub model: GameModelSpecV1,
    /// Initial set.
    pub initial: InitialSetV1,
    /// Target set.
    pub target: TargetSetV1,
    /// Unsafe set.
    pub unsafe_set: UnsafeSetV1,
    /// Control domain.
    pub controls: ControlSetV1,
    /// Disturbance domain.
    pub disturbances: DisturbanceSetV1,
    /// Parameter domain.
    pub parameters: ParameterSetV1,
    /// Objective and requested polarity.
    pub claim: GameClaimSpecV1,
    /// Exact ordered quantifiers.
    pub quantifiers: GameQuantifierPrefixV1,
    /// Nonanticipative information grants.
    pub information: InformationPatternV1,
    /// Strategy representations.
    pub strategies: Vec<GameStrategySpecV1>,
    /// Horizon.
    pub horizon: GameHorizonV1,
    /// Stopping semantics.
    pub stopping: GameStoppingSemanticsV1,
    /// Atomic or composed problem.
    pub composition: GameCompositionV1,
    /// Analysis budget.
    pub budget: GameAnalysisBudgetV1,
}

impl GameProblemIrV1 {
    /// Construct a current-version raw game.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        model: GameModelSpecV1,
        initial: InitialSetV1,
        target: TargetSetV1,
        unsafe_set: UnsafeSetV1,
        controls: ControlSetV1,
        disturbances: DisturbanceSetV1,
        parameters: ParameterSetV1,
        claim: GameClaimSpecV1,
        quantifiers: GameQuantifierPrefixV1,
        information: InformationPatternV1,
        strategies: Vec<GameStrategySpecV1>,
        horizon: GameHorizonV1,
        stopping: GameStoppingSemanticsV1,
        composition: GameCompositionV1,
        budget: GameAnalysisBudgetV1,
    ) -> Self {
        Self::with_schema_version(
            GAME_PROBLEM_SCHEMA_VERSION_V1,
            model,
            initial,
            target,
            unsafe_set,
            controls,
            disturbances,
            parameters,
            claim,
            quantifiers,
            information,
            strategies,
            horizon,
            stopping,
            composition,
            budget,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        model: GameModelSpecV1,
        initial: InitialSetV1,
        target: TargetSetV1,
        unsafe_set: UnsafeSetV1,
        controls: ControlSetV1,
        disturbances: DisturbanceSetV1,
        parameters: ParameterSetV1,
        claim: GameClaimSpecV1,
        quantifiers: GameQuantifierPrefixV1,
        information: InformationPatternV1,
        strategies: Vec<GameStrategySpecV1>,
        horizon: GameHorizonV1,
        stopping: GameStoppingSemanticsV1,
        composition: GameCompositionV1,
        budget: GameAnalysisBudgetV1,
    ) -> Self {
        Self {
            schema_version,
            model,
            initial,
            target,
            unsafe_set,
            controls,
            disturbances,
            parameters,
            claim,
            quantifiers,
            information,
            strategies,
            horizon,
            stopping,
            composition,
            budget,
        }
    }

    /// Declared schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }
}

impl fmt::Display for GameProblemIrV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("quantifiers=[")?;
        fmt::Display::fmt(&self.quantifiers, f)?;
        f.write_str("]; strategies=[")?;
        for (index, strategy) in self.strategies.iter().enumerate() {
            if index != 0 {
                f.write_str("; ")?;
            }
            write!(f, "{}:{}", strategy.player, strategy.representation)?;
        }
        f.write_str("]")
    }
}

/// Collection protected by a hard pre-work cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameCollectionV1 {
    /// Quantifier clauses.
    Quantifiers,
    /// Information grants.
    InformationGrants,
    /// Strategies.
    Strategies,
    /// Total strategy dependencies.
    StrategyDependencies,
    /// Composition components.
    Components,
}

/// Set/context role associated with a mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameContextRoleV1 {
    /// Initial set.
    Initial,
    /// Target set.
    Target,
    /// Unsafe set.
    Unsafe,
    /// Control set.
    Control,
    /// Disturbance set.
    Disturbance,
    /// Parameter set.
    Parameter,
    /// Hybrid modes.
    Modes,
    /// Hybrid events.
    Events,
    /// Composition component.
    Component,
}

/// Invalid scalar field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameFieldV1 {
    /// DAE differentiation index.
    DaeIndex,
    /// Horizon start.
    HorizonStart,
    /// Horizon end.
    HorizonEnd,
    /// Seconds per time unit.
    SecondsPerUnit,
    /// Observation delay.
    ObservationDelay,
    /// Strategy delay.
    StrategyDelay,
    /// Nonanticipative memory bound.
    StrategyMemory,
    /// Cell budget.
    CellBudget,
    /// Transition budget.
    TransitionBudget,
    /// Strategy-node budget.
    StrategyNodeBudget,
    /// Wall-time budget.
    WallTimeBudget,
}

/// Why v1 cannot expose a theorem-eligible game claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameUnknownReasonV1 {
    /// Infinite-horizon theorem modules are not admitted in v1.
    InfiniteHorizon,
    /// Hybrid Zeno behavior is unresolved.
    ZenoUnresolved,
    /// A regularized model is distinct from the original.
    RegularizedHybridModel,
    /// At least one strategy is unresolved.
    StrategyUnresolved,
    /// Requested proof polarity is Unknown.
    PolarityUnresolved,
    /// Set-valued strategy semantics need a later theorem module.
    SetValuedStrategy,
}

/// Derived scientific availability. Identity admission alone proves no game theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GameClaimAvailabilityV1 {
    /// Structurally eligible for a later theorem checker at this polarity.
    Eligible {
        /// Exact requested polarity.
        polarity: GameProofPolarityV1,
    },
    /// Structurally admitted but theorem availability is Unknown.
    Unknown {
        /// Exact fail-closed reason.
        reason: GameUnknownReasonV1,
    },
}

/// Deterministic fail-closed admission issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GameSemanticIssueV1 {
    /// Unsupported decoded schema.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Supported version.
        supported: u32,
    },
    /// Hard cap exceeded before proportional work.
    TooMany {
        /// Capped collection.
        collection: GameCollectionV1,
        /// Items supplied.
        found: usize,
        /// Maximum admitted.
        limit: usize,
    },
    /// A required collection is empty.
    EmptyCollection {
        /// Empty collection.
        collection: GameCollectionV1,
    },
    /// A scalar is non-finite, zero, negative, or otherwise invalid.
    InvalidValue {
        /// Invalid field.
        field: GameFieldV1,
    },
    /// Typed set or component context does not match the model.
    ContextMismatch {
        /// Mismatching role.
        role: GameContextRoleV1,
    },
    /// Unsupported infinite-dimensional dynamics.
    UnsupportedModelClass,
    /// A variable is bound more than once.
    DuplicateQuantifiedVariable {
        /// Duplicated variable.
        variable: GameVariableV1,
    },
    /// Required control or disturbance variable is unbound.
    MissingQuantifiedVariable {
        /// Missing variable.
        variable: GameVariableV1,
    },
    /// Variable and typed domain disagree, or the domain artifact is wrong.
    QuantifierDomainMismatch {
        /// Offending variable.
        variable: GameVariableV1,
    },
    /// Quantifier polarity contradicts the fixed controller/adversary roles.
    QuantifierPolarityMismatch {
        /// Offending variable.
        variable: GameVariableV1,
        /// Supplied polarity.
        found: GameQuantifierV1,
        /// Polarity required by the v1 player model.
        required: GameQuantifierV1,
    },
    /// Duplicate information grant for one player/subject.
    DuplicateInformationGrant,
    /// An information grant itself exposes future values.
    AnticipativeInformation,
    /// Strategy dependency exposes future values.
    AnticipativeStrategy,
    /// Strategy names a hidden or absent observation.
    HiddenOrMissingObservation,
    /// Strategy access is stronger/earlier than the grant permits.
    ObservationTimingMismatch,
    /// Open-loop strategy names an online dependency.
    OpenLoopHasOnlineDependency,
    /// A mode/event observation is used outside a hybrid model.
    HybridInformationForNonHybridModel,
    /// Required controller strategy is absent.
    MissingControllerStrategy,
    /// More than one strategy is supplied for a player.
    DuplicatePlayerStrategy {
        /// Duplicated player.
        player: GamePlayerV1,
    },
    /// A strategy repeats the same dependency subject.
    DuplicateStrategyDependency,
    /// A universal disturbance followed by an existential open-loop control
    /// would lower the prefix into a trajectory-clairvoyant controller.
    ClairvoyantQuantifierLowering,
    /// Stopping semantics conflict with horizon or model class.
    StoppingSemanticsMismatch,
    /// Composition needs at least two components.
    TooFewCompositionComponents,
    /// A component identity is repeated.
    DuplicateCompositionComponent,
    /// Cooperative cancellation occurred before publication.
    Cancelled,
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Complete deterministic refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GameSemanticReportV1 {
    issues: Vec<GameSemanticIssueV1>,
}

impl GameSemanticReportV1 {
    fn new(issues: Vec<GameSemanticIssueV1>) -> Self {
        Self { issues }
    }

    /// Ordered issues.
    #[must_use]
    pub fn issues(&self) -> &[GameSemanticIssueV1] {
        &self.issues
    }
}

impl fmt::Display for GameSemanticReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "quantified game admission refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for GameSemanticReportV1 {}

/// Sealed canonical game problem.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedGameProblemV1 {
    ir: GameProblemIrV1,
    receipt: IdentityReceipt<GameProblemIdV1>,
    availability: GameClaimAvailabilityV1,
}

impl ValidatedGameProblemV1 {
    /// Canonical admitted problem view.
    #[must_use]
    pub const fn ir(&self) -> &GameProblemIrV1 {
        &self.ir
    }

    /// Typed semantic identity.
    #[must_use]
    pub const fn problem_id(&self) -> GameProblemIdV1 {
        self.receipt.id()
    }

    /// Exact identity/preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<GameProblemIdV1> {
        self.receipt
    }

    /// Later-theorem eligibility or exact Unknown reason.
    #[must_use]
    pub const fn claim_availability(&self) -> GameClaimAvailabilityV1 {
        self.availability
    }
}

/// Validate and canonically identify a quantified game problem.
///
/// Hard caps and cooperative cancellation are checked before sorting or
/// identity publication. Successful admission proves schema consistency and
/// deterministic semantics only; it does not prove a winning or viability set.
///
/// # Errors
/// Returns a deterministic report for cap, context, quantifier, information,
/// strategy, horizon, stopping, composition, cancellation, or identity
/// failures.
#[allow(clippy::too_many_lines)]
#[must_use = "game admission must be handled before the problem is used"]
pub fn validate_game_problem_v1(
    mut ir: GameProblemIrV1,
    cx: &Cx<'_>,
) -> Result<ValidatedGameProblemV1, GameSemanticReportV1> {
    checkpoint(cx)?;
    let mut issues = Vec::new();
    let component_count = composition_components(&ir.composition).len();
    check_cap(
        &mut issues,
        GameCollectionV1::Quantifiers,
        ir.quantifiers.clauses.len(),
        MAX_GAME_QUANTIFIERS_V1,
    );
    check_cap(
        &mut issues,
        GameCollectionV1::InformationGrants,
        ir.information.grants.len(),
        MAX_INFORMATION_GRANTS_V1,
    );
    check_cap(
        &mut issues,
        GameCollectionV1::Strategies,
        ir.strategies.len(),
        MAX_GAME_STRATEGIES_V1,
    );
    check_cap(
        &mut issues,
        GameCollectionV1::Components,
        component_count,
        MAX_GAME_COMPONENTS_V1,
    );
    if !issues.is_empty() {
        return Err(GameSemanticReportV1::new(issues));
    }
    let dependency_count = ir.strategies.iter().fold(0_usize, |total, strategy| {
        total.saturating_add(strategy.dependencies.len())
    });
    check_cap(
        &mut issues,
        GameCollectionV1::StrategyDependencies,
        dependency_count,
        MAX_STRATEGY_DEPENDENCIES_V1,
    );
    if !issues.is_empty() {
        return Err(GameSemanticReportV1::new(issues));
    }

    canonicalize_problem_zeros(&mut ir);
    if ir.schema_version != GAME_PROBLEM_SCHEMA_VERSION_V1 {
        issues.push(GameSemanticIssueV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: GAME_PROBLEM_SCHEMA_VERSION_V1,
        });
    }
    if ir.quantifiers.clauses.is_empty() {
        issues.push(GameSemanticIssueV1::EmptyCollection {
            collection: GameCollectionV1::Quantifiers,
        });
    }
    if ir.information.grants.is_empty() {
        issues.push(GameSemanticIssueV1::EmptyCollection {
            collection: GameCollectionV1::InformationGrants,
        });
    }
    if ir.strategies.is_empty() {
        issues.push(GameSemanticIssueV1::EmptyCollection {
            collection: GameCollectionV1::Strategies,
        });
    }

    validate_model_and_contexts(&ir, &mut issues);
    validate_horizon_budget_and_stopping(&ir, &mut issues);
    validate_quantifiers(&ir, &mut issues);

    ir.information
        .grants
        .sort_by_cached_key(information_grant_bytes);
    for pair in ir.information.grants.windows(2) {
        if pair[0].player == pair[1].player && pair[0].subject == pair[1].subject {
            push_once(&mut issues, GameSemanticIssueV1::DuplicateInformationGrant);
        }
    }
    for (index, grant) in ir.information.grants.iter().enumerate() {
        if index.is_multiple_of(32) {
            checkpoint(cx)?;
        }
        validate_information_grant(grant, &ir.model, &mut issues);
    }

    for strategy in &mut ir.strategies {
        strategy
            .dependencies
            .sort_by_cached_key(strategy_dependency_bytes);
        if strategy
            .dependencies
            .windows(2)
            .any(|pair| pair[0].subject == pair[1].subject)
        {
            push_once(
                &mut issues,
                GameSemanticIssueV1::DuplicateStrategyDependency,
            );
        }
    }
    ir.strategies.sort_by_cached_key(strategy_bytes);
    for pair in ir.strategies.windows(2) {
        if pair[0].player == pair[1].player {
            push_once(
                &mut issues,
                GameSemanticIssueV1::DuplicatePlayerStrategy {
                    player: pair[0].player,
                },
            );
        }
    }
    for (index, strategy) in ir.strategies.iter().enumerate() {
        if index.is_multiple_of(16) {
            checkpoint(cx)?;
        }
        validate_strategy(strategy, &ir, &mut issues);
    }
    if !ir
        .strategies
        .iter()
        .any(|strategy| strategy.player == GamePlayerV1::Controller)
    {
        issues.push(GameSemanticIssueV1::MissingControllerStrategy);
    }
    validate_clairvoyant_lowering(&ir, &mut issues);
    canonicalize_and_validate_composition(&mut ir, &mut issues);
    if !issues.is_empty() {
        return Err(GameSemanticReportV1::new(issues));
    }

    let availability = derive_claim_availability(&ir);
    let receipt = game_problem_receipt(&ir, cx).map_err(identity_report)?;
    Ok(ValidatedGameProblemV1 {
        ir,
        receipt,
        availability,
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), GameSemanticReportV1> {
    cx.checkpoint()
        .map_err(|_| GameSemanticReportV1::new(vec![GameSemanticIssueV1::Cancelled]))
}

fn identity_report(error: CanonicalError) -> GameSemanticReportV1 {
    let issue = if matches!(&error, CanonicalError::Cancelled { .. }) {
        GameSemanticIssueV1::Cancelled
    } else {
        GameSemanticIssueV1::Identity(error)
    };
    GameSemanticReportV1::new(vec![issue])
}

fn check_cap(
    issues: &mut Vec<GameSemanticIssueV1>,
    collection: GameCollectionV1,
    found: usize,
    limit: usize,
) {
    if found > limit {
        issues.push(GameSemanticIssueV1::TooMany {
            collection,
            found,
            limit,
        });
    }
}

fn push_once(issues: &mut Vec<GameSemanticIssueV1>, issue: GameSemanticIssueV1) {
    if !issues.contains(&issue) {
        issues.push(issue);
    }
}

fn canonicalize_zero(value: &mut f64) {
    if *value == 0.0 {
        *value = 0.0;
    }
}

fn canonicalize_problem_zeros(ir: &mut GameProblemIrV1) {
    match &mut ir.horizon {
        GameHorizonV1::Finite {
            start,
            end,
            seconds_per_unit,
            ..
        } => {
            canonicalize_zero(start);
            canonicalize_zero(end);
            canonicalize_zero(seconds_per_unit);
        }
        GameHorizonV1::Infinite {
            start,
            seconds_per_unit,
            ..
        } => {
            canonicalize_zero(start);
            canonicalize_zero(seconds_per_unit);
        }
    }
    canonicalize_zero(&mut ir.budget.max_wall_seconds);
    for grant in &mut ir.information.grants {
        if let ObservationAvailabilityV1::Delayed { lag } = &mut grant.availability {
            canonicalize_zero(lag);
        }
    }
    for strategy in &mut ir.strategies {
        for dependency in &mut strategy.dependencies {
            if let StrategyAccessV1::Delayed { lag } = &mut dependency.access {
                canonicalize_zero(lag);
            }
        }
    }
}

fn validate_model_and_contexts(ir: &GameProblemIrV1, issues: &mut Vec<GameSemanticIssueV1>) {
    let expected = StateSetContextV1 {
        state_space: ir.model.state_space,
        frame: ir.model.frame,
        units: ir.model.state_units,
        model_version: ir.model.version,
    };
    for (role, context) in [
        (GameContextRoleV1::Initial, ir.initial.context),
        (GameContextRoleV1::Target, ir.target.context),
        (GameContextRoleV1::Unsafe, ir.unsafe_set.context),
    ] {
        if context != expected {
            issues.push(GameSemanticIssueV1::ContextMismatch { role });
        }
    }
    for (role, version) in [
        (GameContextRoleV1::Control, ir.controls.model_version),
        (
            GameContextRoleV1::Disturbance,
            ir.disturbances.model_version,
        ),
        (GameContextRoleV1::Parameter, ir.parameters.model_version),
    ] {
        if version != ir.model.version {
            issues.push(GameSemanticIssueV1::ContextMismatch { role });
        }
    }
    match ir.model.class {
        GameModelClassV1::AdmittedDae { index: 0, .. } => {
            issues.push(GameSemanticIssueV1::InvalidValue {
                field: GameFieldV1::DaeIndex,
            });
        }
        GameModelClassV1::Hybrid { modes, events, .. } => {
            if modes.model_version != ir.model.version {
                issues.push(GameSemanticIssueV1::ContextMismatch {
                    role: GameContextRoleV1::Modes,
                });
            }
            if events.model_version != ir.model.version {
                issues.push(GameSemanticIssueV1::ContextMismatch {
                    role: GameContextRoleV1::Events,
                });
            }
        }
        GameModelClassV1::UnsupportedInfiniteDimensional => {
            issues.push(GameSemanticIssueV1::UnsupportedModelClass);
        }
        GameModelClassV1::Deterministic
        | GameModelClassV1::DifferentialGame
        | GameModelClassV1::AdmittedDae { .. } => {}
    }
}

fn validate_horizon_budget_and_stopping(
    ir: &GameProblemIrV1,
    issues: &mut Vec<GameSemanticIssueV1>,
) {
    match ir.horizon {
        GameHorizonV1::Finite {
            start,
            end,
            seconds_per_unit,
            ..
        } => {
            if !start.is_finite() {
                issues.push(GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::HorizonStart,
                });
            }
            if !end.is_finite() || end <= start {
                issues.push(GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::HorizonEnd,
                });
            }
            if !(seconds_per_unit.is_finite() && seconds_per_unit > 0.0) {
                issues.push(GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::SecondsPerUnit,
                });
            }
            if matches!(ir.stopping, GameStoppingSemanticsV1::Never) {
                issues.push(GameSemanticIssueV1::StoppingSemanticsMismatch);
            }
        }
        GameHorizonV1::Infinite {
            start,
            seconds_per_unit,
            ..
        } => {
            if !start.is_finite() {
                issues.push(GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::HorizonStart,
                });
            }
            if !(seconds_per_unit.is_finite() && seconds_per_unit > 0.0) {
                issues.push(GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::SecondsPerUnit,
                });
            }
            if matches!(ir.stopping, GameStoppingSemanticsV1::FixedHorizon) {
                issues.push(GameSemanticIssueV1::StoppingSemanticsMismatch);
            }
        }
    }
    if ir.budget.max_cells == 0 {
        issues.push(GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::CellBudget,
        });
    }
    if ir.budget.max_transitions == 0 {
        issues.push(GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::TransitionBudget,
        });
    }
    if ir.budget.max_strategy_nodes == 0 {
        issues.push(GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::StrategyNodeBudget,
        });
    }
    if !(ir.budget.max_wall_seconds.is_finite() && ir.budget.max_wall_seconds > 0.0) {
        issues.push(GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::WallTimeBudget,
        });
    }
    if let GameStoppingSemanticsV1::HybridTerminalEvent { events } = ir.stopping {
        if !matches!(
            ir.model.class,
            GameModelClassV1::Hybrid {
                events: GameEventSetV1 {
                    events: admitted,
                    ..
                },
                ..
            } if admitted == events
        ) {
            issues.push(GameSemanticIssueV1::StoppingSemanticsMismatch);
        }
    }
}

fn validate_quantifiers(ir: &GameProblemIrV1, issues: &mut Vec<GameSemanticIssueV1>) {
    let mut seen = Vec::new();
    for clause in &ir.quantifiers.clauses {
        if seen.contains(&clause.variable) {
            issues.push(GameSemanticIssueV1::DuplicateQuantifiedVariable {
                variable: clause.variable,
            });
        } else {
            seen.push(clause.variable);
        }
        if !domain_matches(clause.variable, clause.domain, ir) {
            issues.push(GameSemanticIssueV1::QuantifierDomainMismatch {
                variable: clause.variable,
            });
        }
        let required = match clause.variable {
            GameVariableV1::Control => Some(GameQuantifierV1::Exists),
            GameVariableV1::Disturbance => Some(GameQuantifierV1::ForAll),
            GameVariableV1::InitialState | GameVariableV1::Parameter => None,
        };
        if let Some(required) = required
            && clause.quantifier != required
        {
            issues.push(GameSemanticIssueV1::QuantifierPolarityMismatch {
                variable: clause.variable,
                found: clause.quantifier,
                required,
            });
        }
    }
    for required in [GameVariableV1::Control, GameVariableV1::Disturbance] {
        if !seen.contains(&required) {
            issues.push(GameSemanticIssueV1::MissingQuantifiedVariable { variable: required });
        }
    }
}

fn domain_matches(
    variable: GameVariableV1,
    domain: GameQuantifierDomainV1,
    ir: &GameProblemIrV1,
) -> bool {
    match (variable, domain) {
        (GameVariableV1::InitialState, GameQuantifierDomainV1::Initial(id)) => id == ir.initial.set,
        (GameVariableV1::Control, GameQuantifierDomainV1::Control(id)) => id == ir.controls.set,
        (GameVariableV1::Disturbance, GameQuantifierDomainV1::Disturbance(id)) => {
            id == ir.disturbances.set
        }
        (GameVariableV1::Parameter, GameQuantifierDomainV1::Parameter(id)) => {
            id == ir.parameters.set
        }
        _ => false,
    }
}

fn validate_information_grant(
    grant: &InformationGrantV1,
    model: &GameModelSpecV1,
    issues: &mut Vec<GameSemanticIssueV1>,
) {
    match grant.availability {
        ObservationAvailabilityV1::Delayed { lag } => {
            if !(lag.is_finite() && lag >= 0.0) {
                push_once(
                    issues,
                    GameSemanticIssueV1::InvalidValue {
                        field: GameFieldV1::ObservationDelay,
                    },
                );
            }
        }
        ObservationAvailabilityV1::FutureTrajectory => {
            push_once(issues, GameSemanticIssueV1::AnticipativeInformation);
        }
        ObservationAvailabilityV1::Hidden
        | ObservationAvailabilityV1::InitialOnly
        | ObservationAvailabilityV1::Current
        | ObservationAvailabilityV1::HistoryThroughCurrent => {}
    }
    if matches!(
        grant.subject,
        InformationSubjectV1::Mode | InformationSubjectV1::Event
    ) && !matches!(model.class, GameModelClassV1::Hybrid { .. })
    {
        push_once(
            issues,
            GameSemanticIssueV1::HybridInformationForNonHybridModel,
        );
    }
}

fn validate_strategy(
    strategy: &GameStrategySpecV1,
    ir: &GameProblemIrV1,
    issues: &mut Vec<GameSemanticIssueV1>,
) {
    if let StrategyRepresentationV1::NonanticipativeFeedback {
        memory_states: 0, ..
    } = strategy.representation
    {
        push_once(
            issues,
            GameSemanticIssueV1::InvalidValue {
                field: GameFieldV1::StrategyMemory,
            },
        );
    }
    if matches!(
        strategy.representation,
        StrategyRepresentationV1::HybridModeFeedback { .. }
    ) && !matches!(ir.model.class, GameModelClassV1::Hybrid { .. })
    {
        push_once(
            issues,
            GameSemanticIssueV1::HybridInformationForNonHybridModel,
        );
    }
    for dependency in &strategy.dependencies {
        if matches!(dependency.access, StrategyAccessV1::FutureTrajectory) {
            push_once(issues, GameSemanticIssueV1::AnticipativeStrategy);
        }
        if let StrategyAccessV1::Delayed { lag } = dependency.access
            && !(lag.is_finite() && lag >= 0.0)
        {
            push_once(
                issues,
                GameSemanticIssueV1::InvalidValue {
                    field: GameFieldV1::StrategyDelay,
                },
            );
        }
        if matches!(
            dependency.subject,
            InformationSubjectV1::Mode | InformationSubjectV1::Event
        ) && !matches!(ir.model.class, GameModelClassV1::Hybrid { .. })
        {
            push_once(
                issues,
                GameSemanticIssueV1::HybridInformationForNonHybridModel,
            );
        }
        let grant =
            ir.information.grants.iter().find(|grant| {
                grant.player == strategy.player && grant.subject == dependency.subject
            });
        match grant {
            None
            | Some(InformationGrantV1 {
                availability: ObservationAvailabilityV1::Hidden,
                ..
            }) => {
                push_once(issues, GameSemanticIssueV1::HiddenOrMissingObservation);
            }
            Some(grant) if !access_is_permitted(dependency.access, grant.availability) => {
                push_once(issues, GameSemanticIssueV1::ObservationTimingMismatch);
            }
            Some(_) => {}
        }
        if matches!(
            strategy.representation,
            StrategyRepresentationV1::OpenLoop { .. }
        ) && !matches!(dependency.access, StrategyAccessV1::InitialOnly)
        {
            push_once(issues, GameSemanticIssueV1::OpenLoopHasOnlineDependency);
        }
    }
}

fn access_is_permitted(access: StrategyAccessV1, availability: ObservationAvailabilityV1) -> bool {
    match availability {
        ObservationAvailabilityV1::Hidden | ObservationAvailabilityV1::FutureTrajectory => false,
        ObservationAvailabilityV1::InitialOnly => {
            matches!(access, StrategyAccessV1::InitialOnly)
        }
        ObservationAvailabilityV1::Current => matches!(
            access,
            StrategyAccessV1::InitialOnly | StrategyAccessV1::Current
        ),
        ObservationAvailabilityV1::HistoryThroughCurrent => {
            !matches!(access, StrategyAccessV1::FutureTrajectory)
        }
        ObservationAvailabilityV1::Delayed {
            lag: observation_lag,
        } => match access {
            StrategyAccessV1::InitialOnly => observation_lag == 0.0,
            StrategyAccessV1::Delayed { lag: strategy_lag } => strategy_lag >= observation_lag,
            StrategyAccessV1::Current
            | StrategyAccessV1::HistoryThroughCurrent
            | StrategyAccessV1::FutureTrajectory => false,
        },
    }
}

fn validate_clairvoyant_lowering(ir: &GameProblemIrV1, issues: &mut Vec<GameSemanticIssueV1>) {
    let disturbance_forall = ir.quantifiers.clauses.iter().position(|clause| {
        clause.variable == GameVariableV1::Disturbance
            && clause.quantifier == GameQuantifierV1::ForAll
    });
    let control_exists = ir.quantifiers.clauses.iter().position(|clause| {
        clause.variable == GameVariableV1::Control && clause.quantifier == GameQuantifierV1::Exists
    });
    let controller_is_open_loop = ir.strategies.iter().any(|strategy| {
        strategy.player == GamePlayerV1::Controller
            && matches!(
                strategy.representation,
                StrategyRepresentationV1::OpenLoop { .. }
            )
    });
    if disturbance_forall
        .zip(control_exists)
        .is_some_and(|(disturbance, control)| disturbance < control)
        && controller_is_open_loop
    {
        issues.push(GameSemanticIssueV1::ClairvoyantQuantifierLowering);
    }
}

fn canonicalize_and_validate_composition(
    ir: &mut GameProblemIrV1,
    issues: &mut Vec<GameSemanticIssueV1>,
) {
    let components = match &mut ir.composition {
        GameCompositionV1::Atomic => return,
        GameCompositionV1::Sequential { components, .. } => components,
        GameCompositionV1::Parallel { components, .. }
        | GameCompositionV1::AdversarialProduct { components, .. } => {
            components.sort_by(|left, right| left.problem.as_bytes().cmp(right.problem.as_bytes()));
            components
        }
    };
    if components.len() < 2 {
        issues.push(GameSemanticIssueV1::TooFewCompositionComponents);
    }
    if components
        .windows(2)
        .any(|pair| pair[0].problem == pair[1].problem)
    {
        issues.push(GameSemanticIssueV1::DuplicateCompositionComponent);
    } else {
        let mut identities: Vec<_> = components.iter().map(|item| item.problem).collect();
        identities.sort_by(|left, right| left.as_bytes().cmp(right.as_bytes()));
        if identities.windows(2).any(|pair| pair[0] == pair[1]) {
            issues.push(GameSemanticIssueV1::DuplicateCompositionComponent);
        }
    }
    for component in components {
        if component.model_version != ir.model.version
            || component.state_space != ir.model.state_space
            || component.frame != ir.model.frame
            || component.units != ir.model.state_units
        {
            push_once(
                issues,
                GameSemanticIssueV1::ContextMismatch {
                    role: GameContextRoleV1::Component,
                },
            );
        }
    }
}

fn derive_claim_availability(ir: &GameProblemIrV1) -> GameClaimAvailabilityV1 {
    if matches!(ir.horizon, GameHorizonV1::Infinite { .. }) {
        return GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::InfiniteHorizon,
        };
    }
    if let GameModelClassV1::Hybrid { zeno, .. } = ir.model.class {
        match zeno {
            HybridZenoScopeV1::Unresolved { .. } => {
                return GameClaimAvailabilityV1::Unknown {
                    reason: GameUnknownReasonV1::ZenoUnresolved,
                };
            }
            HybridZenoScopeV1::Regularized { .. } => {
                return GameClaimAvailabilityV1::Unknown {
                    reason: GameUnknownReasonV1::RegularizedHybridModel,
                };
            }
            HybridZenoScopeV1::Excluded { .. } => {}
        }
    }
    if ir.strategies.iter().any(|strategy| {
        matches!(
            strategy.representation,
            StrategyRepresentationV1::Unknown { .. }
        )
    }) {
        return GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::StrategyUnresolved,
        };
    }
    if ir.strategies.iter().any(|strategy| {
        matches!(
            strategy.representation,
            StrategyRepresentationV1::SetValued { .. }
        )
    }) {
        return GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::SetValuedStrategy,
        };
    }
    if matches!(ir.claim.polarity, GameProofPolarityV1::Unknown { .. }) {
        return GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::PolarityUnresolved,
        };
    }
    GameClaimAvailabilityV1::Eligible {
        polarity: ir.claim.polarity,
    }
}

fn composition_components(composition: &GameCompositionV1) -> &[GameComponentRefV1] {
    match composition {
        GameCompositionV1::Atomic => &[],
        GameCompositionV1::Sequential { components, .. }
        | GameCompositionV1::Parallel { components, .. }
        | GameCompositionV1::AdversarialProduct { components, .. } => components,
    }
}

fn domain_name(domain: GameQuantifierDomainV1) -> &'static str {
    match domain {
        GameQuantifierDomainV1::Initial(_) => "initial-set",
        GameQuantifierDomainV1::Control(_) => "control-set",
        GameQuantifierDomainV1::Disturbance(_) => "disturbance-set",
        GameQuantifierDomainV1::Parameter(_) => "parameter-set",
    }
}

fn write_domain(f: &mut fmt::Formatter<'_>, domain: GameQuantifierDomainV1) -> fmt::Result {
    write!(f, "{}(", domain_name(domain))?;
    let bytes = match &domain {
        GameQuantifierDomainV1::Initial(id) => id.as_bytes(),
        GameQuantifierDomainV1::Control(id) => id.as_bytes(),
        GameQuantifierDomainV1::Disturbance(id) => id.as_bytes(),
        GameQuantifierDomainV1::Parameter(id) => id.as_bytes(),
    };
    for byte in bytes {
        write!(f, "{byte:02x}")?;
    }
    f.write_str(")")
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

fn push_state_context(out: &mut Vec<u8>, context: StateSetContextV1) {
    push_digest(out, context.state_space);
    push_digest(out, context.frame);
    push_digest(out, context.units);
    push_digest(out, context.model_version);
}

fn model_bytes(model: GameModelSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(352);
    push_digest(&mut out, model.model);
    push_digest(&mut out, model.version);
    push_digest(&mut out, model.state_space);
    push_digest(&mut out, model.frame);
    push_digest(&mut out, model.state_units);
    match model.class {
        GameModelClassV1::Deterministic => out.push(0),
        GameModelClassV1::DifferentialGame => out.push(1),
        GameModelClassV1::Hybrid {
            modes,
            events,
            zeno,
        } => {
            out.push(2);
            push_digest(&mut out, modes.modes);
            push_digest(&mut out, modes.model_version);
            push_digest(&mut out, events.events);
            push_digest(&mut out, events.model_version);
            match zeno {
                HybridZenoScopeV1::Excluded { witness } => {
                    out.push(0);
                    push_digest(&mut out, witness);
                }
                HybridZenoScopeV1::Unresolved { no_claim } => {
                    out.push(1);
                    push_digest(&mut out, no_claim);
                }
                HybridZenoScopeV1::Regularized {
                    regularization,
                    no_equivalence,
                } => {
                    out.push(2);
                    push_digest(&mut out, regularization);
                    push_digest(&mut out, no_equivalence);
                }
            }
        }
        GameModelClassV1::AdmittedDae { index, constraint } => {
            out.push(3);
            out.push(index);
            push_digest(&mut out, constraint);
        }
        GameModelClassV1::UnsupportedInfiniteDimensional => out.push(4),
    }
    out
}

fn state_sets_bytes(ir: &GameProblemIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(480);
    push_digest(&mut out, ir.initial.set);
    push_state_context(&mut out, ir.initial.context);
    push_digest(&mut out, ir.target.set);
    push_state_context(&mut out, ir.target.context);
    push_digest(&mut out, ir.unsafe_set.set);
    push_state_context(&mut out, ir.unsafe_set.context);
    out
}

fn action_sets_bytes(ir: &GameProblemIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(288);
    push_digest(&mut out, ir.controls.set);
    push_digest(&mut out, ir.controls.units);
    push_digest(&mut out, ir.controls.model_version);
    push_digest(&mut out, ir.disturbances.set);
    push_digest(&mut out, ir.disturbances.units);
    push_digest(&mut out, ir.disturbances.model_version);
    push_digest(&mut out, ir.parameters.set);
    push_digest(&mut out, ir.parameters.units);
    push_digest(&mut out, ir.parameters.model_version);
    out
}

fn claim_bytes(claim: GameClaimSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(40);
    out.push(match claim.objective {
        GameObjectiveV1::ReachAvoid => 0,
        GameObjectiveV1::Viability => 1,
    });
    match claim.polarity {
        GameProofPolarityV1::Inner => out.push(0),
        GameProofPolarityV1::Outer => out.push(1),
        GameProofPolarityV1::Exact => out.push(2),
        GameProofPolarityV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn variable_tag(variable: GameVariableV1) -> u8 {
    match variable {
        GameVariableV1::InitialState => 0,
        GameVariableV1::Control => 1,
        GameVariableV1::Disturbance => 2,
        GameVariableV1::Parameter => 3,
    }
}

fn player_tag(player: GamePlayerV1) -> u8 {
    match player {
        GamePlayerV1::Controller => 0,
        GamePlayerV1::Disturbance => 1,
    }
}

fn subject_tag(subject: InformationSubjectV1) -> u8 {
    match subject {
        InformationSubjectV1::State => 0,
        InformationSubjectV1::Mode => 1,
        InformationSubjectV1::Event => 2,
        InformationSubjectV1::Control => 3,
        InformationSubjectV1::Disturbance => 4,
        InformationSubjectV1::Parameter => 5,
    }
}

fn quantifier_clause_bytes(clause: &GameQuantifierClauseV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(35);
    out.push(match clause.quantifier {
        GameQuantifierV1::Exists => 0,
        GameQuantifierV1::ForAll => 1,
    });
    out.push(variable_tag(clause.variable));
    match clause.domain {
        GameQuantifierDomainV1::Initial(id) => {
            out.push(0);
            push_digest(&mut out, id);
        }
        GameQuantifierDomainV1::Control(id) => {
            out.push(1);
            push_digest(&mut out, id);
        }
        GameQuantifierDomainV1::Disturbance(id) => {
            out.push(2);
            push_digest(&mut out, id);
        }
        GameQuantifierDomainV1::Parameter(id) => {
            out.push(3);
            push_digest(&mut out, id);
        }
    }
    out
}

fn information_grant_bytes(grant: &InformationGrantV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(44);
    out.push(player_tag(grant.player));
    out.push(subject_tag(grant.subject));
    push_digest(&mut out, grant.observation);
    match grant.availability {
        ObservationAvailabilityV1::Hidden => out.push(0),
        ObservationAvailabilityV1::InitialOnly => out.push(1),
        ObservationAvailabilityV1::Current => out.push(2),
        ObservationAvailabilityV1::HistoryThroughCurrent => out.push(3),
        ObservationAvailabilityV1::Delayed { lag } => {
            out.push(4);
            push_f64(&mut out, lag);
        }
        ObservationAvailabilityV1::FutureTrajectory => out.push(5),
    }
    out
}

fn strategy_dependency_bytes(dependency: &StrategyDependencyV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(10);
    out.push(subject_tag(dependency.subject));
    match dependency.access {
        StrategyAccessV1::InitialOnly => out.push(0),
        StrategyAccessV1::Current => out.push(1),
        StrategyAccessV1::HistoryThroughCurrent => out.push(2),
        StrategyAccessV1::Delayed { lag } => {
            out.push(3);
            push_f64(&mut out, lag);
        }
        StrategyAccessV1::FutureTrajectory => out.push(4),
    }
    out
}

fn strategy_bytes(strategy: &GameStrategySpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    out.push(player_tag(strategy.player));
    match strategy.representation {
        StrategyRepresentationV1::OpenLoop { artifact } => {
            out.push(0);
            push_digest(&mut out, artifact);
        }
        StrategyRepresentationV1::StateFeedback { artifact } => {
            out.push(1);
            push_digest(&mut out, artifact);
        }
        StrategyRepresentationV1::NonanticipativeFeedback {
            artifact,
            memory_states,
        } => {
            out.push(2);
            push_digest(&mut out, artifact);
            push_u32(&mut out, memory_states);
        }
        StrategyRepresentationV1::HybridModeFeedback { artifact } => {
            out.push(3);
            push_digest(&mut out, artifact);
        }
        StrategyRepresentationV1::SetValued { artifact } => {
            out.push(4);
            push_digest(&mut out, artifact);
        }
        StrategyRepresentationV1::Unknown { no_claim } => {
            out.push(5);
            push_digest(&mut out, no_claim);
        }
    }
    push_u64(&mut out, strategy.dependencies.len() as u64);
    for dependency in &strategy.dependencies {
        let encoded = strategy_dependency_bytes(dependency);
        push_u64(&mut out, encoded.len() as u64);
        out.extend_from_slice(&encoded);
    }
    out
}

fn horizon_bytes(horizon: GameHorizonV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(96);
    match horizon {
        GameHorizonV1::Finite {
            start,
            end,
            unit,
            seconds_per_unit,
        } => {
            out.push(0);
            push_f64(&mut out, start);
            push_f64(&mut out, end);
            push_digest(&mut out, unit);
            push_f64(&mut out, seconds_per_unit);
        }
        GameHorizonV1::Infinite {
            start,
            unit,
            seconds_per_unit,
            no_claim,
        } => {
            out.push(1);
            push_f64(&mut out, start);
            push_digest(&mut out, unit);
            push_f64(&mut out, seconds_per_unit);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn stopping_bytes(stopping: GameStoppingSemanticsV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(40);
    match stopping {
        GameStoppingSemanticsV1::FixedHorizon => out.push(0),
        GameStoppingSemanticsV1::FirstTarget => out.push(1),
        GameStoppingSemanticsV1::FirstTargetOrUnsafe => out.push(2),
        GameStoppingSemanticsV1::StateDependent { rule } => {
            out.push(3);
            push_digest(&mut out, rule);
        }
        GameStoppingSemanticsV1::HybridTerminalEvent { events } => {
            out.push(4);
            push_digest(&mut out, events);
        }
        GameStoppingSemanticsV1::Never => out.push(5),
    }
    out
}

fn push_component(out: &mut Vec<u8>, component: GameComponentRefV1) {
    out.extend_from_slice(component.problem.as_bytes());
    push_digest(out, component.model_version);
    push_digest(out, component.state_space);
    push_digest(out, component.frame);
    push_digest(out, component.units);
}

fn composition_bytes(composition: &GameCompositionV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(512);
    match composition {
        GameCompositionV1::Atomic => out.push(0),
        GameCompositionV1::Sequential {
            components,
            interface,
        } => {
            out.push(1);
            push_digest(&mut out, *interface);
            push_u64(&mut out, components.len() as u64);
            for component in components {
                push_component(&mut out, *component);
            }
        }
        GameCompositionV1::Parallel {
            components,
            interface,
        } => {
            out.push(2);
            push_digest(&mut out, *interface);
            push_u64(&mut out, components.len() as u64);
            for component in components {
                push_component(&mut out, *component);
            }
        }
        GameCompositionV1::AdversarialProduct {
            components,
            interface,
        } => {
            out.push(3);
            push_digest(&mut out, *interface);
            push_u64(&mut out, components.len() as u64);
            for component in components {
                push_component(&mut out, *component);
            }
        }
    }
    out
}

fn budget_bytes(budget: GameAnalysisBudgetV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(64);
    push_digest(&mut out, budget.cell_decomposition);
    push_u64(&mut out, budget.max_cells);
    push_u64(&mut out, budget.max_transitions);
    push_u64(&mut out, budget.max_strategy_nodes);
    push_f64(&mut out, budget.max_wall_seconds);
    out
}

fn game_problem_receipt(
    ir: &GameProblemIrV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<GameProblemIdV1>, CanonicalError> {
    let model = model_bytes(ir.model);
    let state_sets = state_sets_bytes(ir);
    let action_sets = action_sets_bytes(ir);
    let claim = claim_bytes(ir.claim);
    let quantifiers: Vec<_> = ir
        .quantifiers
        .clauses
        .iter()
        .map(quantifier_clause_bytes)
        .collect();
    let information: Vec<_> = ir
        .information
        .grants
        .iter()
        .map(information_grant_bytes)
        .collect();
    let strategies: Vec<_> = ir.strategies.iter().map(strategy_bytes).collect();
    let horizon = horizon_bytes(ir.horizon);
    let stopping = stopping_bytes(ir.stopping);
    let composition = composition_bytes(&ir.composition);
    let budget = budget_bytes(ir.budget);

    CanonicalEncoder::<GameProblemIdV1, _>::new(GAME_IDENTITY_LIMITS, || cx.is_cancel_requested())?
        .bytes(Field::new(0, "model"), &model)?
        .bytes(Field::new(1, "state-sets"), &state_sets)?
        .bytes(Field::new(2, "action-sets"), &action_sets)?
        .bytes(Field::new(3, "claim"), &claim)?
        .ordered_bytes(
            Field::new(4, "quantifier-prefix"),
            quantifiers.len() as u64,
            quantifiers.iter().map(Vec::as_slice),
        )?
        .canonical_set(
            Field::new(5, "information-pattern"),
            information.len() as u64,
            information.iter().map(Vec::as_slice),
        )?
        .canonical_set(
            Field::new(6, "strategies"),
            strategies.len() as u64,
            strategies.iter().map(Vec::as_slice),
        )?
        .bytes(Field::new(7, "horizon"), &horizon)?
        .bytes(Field::new(8, "stopping"), &stopping)?
        .bytes(Field::new(9, "composition"), &composition)?
        .bytes(Field::new(10, "analysis-budget"), &budget)?
        .finish()
}

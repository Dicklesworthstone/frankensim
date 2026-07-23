//! fs-scenario: the boundary-condition and load-case ALGEBRA (plan
//! patch Rev D). A Region answers "what is the domain?"; a PhysicsModel
//! answers "which equations?"; a `Scenario` answers "WHAT IS BEING DONE
//! TO IT?" — as a typed value with dimensional analysis (fs-qty),
//! provenance (seed + canonical IR), and validity checks. Boundary
//! conditions are where simulations quietly become invalid; this crate
//! makes that class of mistake a structured, fixable refusal instead.
//!
//! Scenario objects name geometry through the [`entity`] module's persistent
//! identities (`Assembly -> Part -> Region | Surface | Interface`) rather than
//! through strings a rename can silently orphan.
//!
//! Layer: L3 (FLUX support). Runtime deps: `std`, fs-blake3, fs-qty,
//! fs-rand, fs-cheb, fs-exec, fs-ga, fs-ivl, fs-motion, fs-math. The Design Ledger stores
//! scenarios as canonical-IR artifacts — that integration lives ABOVE this layer
//! (exercised here via a dev-dependency in conformance tests).

pub mod bc;
pub mod ensemble;
pub mod entity;
pub mod frame;
pub mod ir;
pub mod payload;
pub mod scenario;
pub mod sensor;
pub mod signal;

pub use bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
pub use ensemble::{
    DEFAULT_REALIZATION_BUDGET, Realization, RealizationBudget, SpectrumModel, StochasticEnsemble,
};
pub use entity::{
    Binding, BindingRow, BindingTable, ContactSide, Correspondence, DEFAULT_ENTITY_BUDGET, Datum,
    DatumFeature, DatumId, Entity, EntityBudget, EntityCatalog, EntityDeclaration, EntityError,
    EntityId, EntityKind, EntityRef, EntityStatus, EvidenceTier, GeometryFingerprint,
    IdentityReceipt, ImportOutcome, ImportRevision, ImportScope, ImportStep, ImportedEntity,
    InterfacePair, InterfacePairing, KindExpectation, LegacyMigration, MatchBasis, NameLookup,
    Placement, PlacementBasis, RebindEvent, ReferenceSite, Resolution, ResolutionFault, Tolerance,
    ToleranceKind, ToleranceSource, migrate_legacy_scenario, scenario_reference_sites,
    validate_bindings,
};
pub use frame::{
    Frame, FrameId, FrameMotion, FrameMotionKind, FrameMotionLoweringError, FrameTree,
    FrameTreeMotorPath, WORLD,
};
pub use scenario::{
    Combination, ContactLaw, ContactModel, DEFAULT_VALIDATION_BUDGET, Environment, LoadCase,
    Scenario, ValidationBudget, ValidationError, ValidationPlan, Violation,
};
pub use sensor::{
    CompiledSensorOperator, MAX_SENSOR_STATE_DIMENSION, MAX_SENSOR_SUPPORT_TERMS,
    MAX_SENSOR_TEXT_BYTES, ObservationSupport, ObservationTerm, PlacementUncertainty,
    SENSOR_IDENTITY_DOMAIN, SENSOR_SCHEMA_VERSION, ScenarioSensor, SensorCalibration,
    SensorComparison, SensorDynamics, SensorError, SensorKind, SensorLocation, SensorMount,
    SensorObservationParts, SensorQuantity,
};
pub use signal::{ChebProfile, Interp, TimeSignal};

use core::fmt;

/// Half-open byte span in one supplied scenario-IR source artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct IrSourceSpan {
    /// First byte belonging to the offending form or token.
    pub start: usize,
    /// First byte after the offending form or token.
    pub end: usize,
}

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured scenario failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum ScenarioError {
    /// A value's dimensions disagree with what the context demands.
    Dimensions {
        /// What was being evaluated.
        context: String,
        /// Expected SI exponents (m, kg, s, K, A, mol).
        expected: [i8; 6],
        /// Supplied exponents.
        got: [i8; 6],
    },
    /// A frame reference could not be resolved.
    Frame {
        /// Diagnosis.
        what: String,
    },
    /// An evaluation was structurally impossible (empty table, bad time).
    Evaluate {
        /// Diagnosis.
        what: String,
    },
    /// Canonical-IR text failed to parse.
    Parse {
        /// Exact half-open byte span of the offending form or token.
        span: IrSourceSpan,
        /// Deterministic structural path within the scenario form.
        path: String,
        /// Diagnosis.
        what: String,
    },
    /// A Machine-IR graph role was presented in a boundary-kind slot.
    ReservedBoundaryRole {
        /// Exact reserved wire role that was refused.
        role: &'static str,
        /// Exact half-open byte span of the reserved role token.
        span: IrSourceSpan,
        /// Deterministic structural path within the scenario form.
        path: String,
    },
    /// IR codec resource admission or fallible allocation refused before any
    /// scenario or canonical byte string was published.
    Resource {
        /// Stable codec operation (`decode` or `encode`).
        operation: &'static str,
        /// Deterministic phase within the operation.
        phase: &'static str,
        /// Stable resource name (`heap_bytes`, `output_bytes`, or `work`).
        resource: &'static str,
        /// Requested resource units.
        requested: u128,
        /// Admitted resource units.
        limit: u128,
        /// Conservatively completed logical work units.
        completed: u128,
        /// Preflighted logical work units, or zero before preflight.
        planned: u128,
    },
}

impl fmt::Display for ScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScenarioError::Dimensions {
                context,
                expected,
                got,
            } => write!(
                f,
                "{context}: dimension mismatch (expected {expected:?}, got {got:?})"
            ),
            ScenarioError::Frame { what } => write!(f, "frame error: {what}"),
            ScenarioError::Evaluate { what } => write!(f, "evaluation error: {what}"),
            ScenarioError::Parse { span, path, what } => write!(
                f,
                "IR parse error at {path} (bytes {}..{}): {what}",
                span.start, span.end
            ),
            ScenarioError::ReservedBoundaryRole { role, span, path } => write!(
                f,
                "reserved machine-graph role {role:?} at {path} (bytes {}..{}) cannot be encoded as a boundary-condition kind; declare it as a typed fs-ir machine relation",
                span.start, span.end
            ),
            ScenarioError::Resource {
                operation,
                phase,
                resource,
                requested,
                limit,
                completed,
                planned,
            } => write!(
                f,
                "scenario IR {operation} refused during {phase}: {resource} request {requested} exceeds or could not satisfy limit {limit} after {completed}/{planned} planned work units"
            ),
        }
    }
}

impl std::error::Error for ScenarioError {}

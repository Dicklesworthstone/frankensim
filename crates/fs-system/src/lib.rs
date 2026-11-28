//! Composition kernel: components joined by connections, stepped as one.
//!
//! A coupled device is a set of [`fs_ports::component::Component`]s and the
//! connections between their ports. This crate owns what every pairing used
//! to write for itself: the state store with atomic publication, the
//! equations of the connections, the order of evaluation, the interface
//! solve, step negotiation, and the energy, entropy and declared stream-mass ledger.
//!
//! ```text
//! SystemBuilder::add / add_with_tangent      components and initial states
//!              ::bond / common_effort / common_flow / heat_route / signal
//!              ::transformer / power_map    lossless power conversion
//!              ::power_map_members           one map from several ports
//!              ::common_flow_ratios / common_effort_ratios
//!                                            a junction whose members each
//!                                            sit behind their own ratio
//!              ::signal_map                 select or mix vector signal coordinates
//!              ::boundary / scheduled_boundary  prescribed inputs and events
//!              ::seal -> System              admission and plan
//! System::step                               negotiate, solve, audit, admit, publish
//! System::ledger                             window accounts
//! System::checkpoint / restore               accepted-state rewind and replay
//! System::step_adaptive                      state/account error-controlled steps
//! System::step_to_event                      locate a state-triggered threshold
//! System::step_to_events                     locate the earliest threshold in a set
//! ```
//!
//! [`elements`] holds ideal lumped elements (mass, spring, damper, force,
//! heat capacity, conductance, heater), their power-domain variants, the
//! dual storage, admittance and flow source for common-effort networks, and
//! the field in a gap of variable width that joins a circuit to a mass.
//! Mechanical, electrical and thermal devices can be described as data and
//! run through the same kernel. With feature `phs`, `circuit` runs any
//! Kirchhoff netlist of fs-phs as one Component whose resistors heat what
//! their heat ports are joined to, and `battery` does the same for the
//! equivalent-circuit battery of fs-phs and its two resistances.
//! `commutated_drive` runs the battery-fed commutated motor drive of fs-phs,
//! one gate program per step, its resistances and switches heating their
//! ports.
//!
//! # Version 0
//!
//! The contract in fs-ports is v0 and so is this driver. It solves looped
//! blocks by a Newton iteration on their interface inputs, using each
//! component's port tangent when it offers one and a forward difference
//! otherwise; the inputs that are copies are substituted and the ones that
//! close the loops are eliminated together. Blocks whose owners all offer
//! `StageResidual` can instead solve their internal and interface equations
//! together, including owners with no internal unknowns. Multirate stepping
//! uses explicitly configured rate groups. Durable snapshots require every
//! component to supply a native state codec.
//! Step-doubling error control uses an explicitly supplied state
//! norm and convergence order. The remaining work is in beads `.5.3`, `.5.5` and
//! `.5.7` to `.5.10` of `frankensim-rc-1109-root-3gc12`.

#![deny(unsafe_code)]
#![warn(missing_docs)]

pub mod accelerator;
#[cfg(feature = "plate-slices")]
pub mod acoustic;
mod adaptive;
#[cfg(feature = "phs")]
pub mod battery;
#[cfg(feature = "phs")]
pub mod circuit;
#[cfg(feature = "phs")]
pub mod commutated_drive;
#[cfg(feature = "constitutive")]
pub mod constitutive;
pub mod diagnostics;
pub mod elements;
mod erased;
pub mod ledger;
pub mod lowering;
#[cfg(feature = "phs")]
pub mod phs;
#[cfg(feature = "phs")]
pub mod phs_exact;
mod solve;
mod system;
mod wiring;

pub use accelerator::{
    AcceptedStateCache, AitkenConfig, BlockAcceleratorState, ControllerState, IqnIlsConfig,
    MAX_IQN_HISTORY, StateCacheKey,
};
pub use adaptive::{AdaptiveConfig, AdaptiveError, AdaptiveReport, ErrorNorm, ErrorTolerance};
pub use diagnostics::{
    BlockSolveMode, BlockSolverPolicy, BlockStepDiagnostic, BudgetKind, FixedPointMethod,
    MonolithicConfig, PartialStepOutcome, ReceiptDiagnostic, StepBudgets, StepDiagnostic,
    StepVerdict,
};
pub use erased::{
    LineNotOffered, OwnerError, StageResidualNotOffered, StateTypeMismatch, TangentNotOffered,
};
pub use ledger::{
    CarriedInventoryOutcome, EnergyLeakAttribution, Ledger, LedgerExport, WindowLedgerOutcome,
};
pub use lowering::{
    AcousticAssemblySpec, AcousticFamilyConstructor, AcousticFamilyContext, AcousticFamilyKey,
    AcousticMemberLoweringReceipt, AcousticMemberSpec, CardStore, ComponentLoweringReceipt,
    ConstitutiveModelRef, FidelityRungKind, FidelitySelectionPolicy, GearMeshConstructor,
    GearMeshLoweringContext, GearMeshLoweringKey, GearMeshLoweringReceipt, GroupConstructor,
    GroupLoweringContext, GroupSelection, GroupSelector, InterfaceCardEntry, InterfaceConstructor,
    InterfaceKindTag, InterfaceLoweringContext, InterfaceLoweringKey, InterfaceLoweringReceipt,
    JointConstructor, JointLoweringContext, JointLoweringKey, JointLoweringReceipt,
    LOWERING_GEOMETRY_HASH_DOMAIN, LOWERING_RECEIPT_HASH_DOMAIN, LOWERING_RECEIPT_SCHEMA,
    LoweredAcousticMemberOutput, LoweredGearMeshOutput, LoweredGroupOutput, LoweredInterfaceOutput,
    LoweredJointOutput, LoweredPartOutput, LoweredSystem, LoweringError, LoweringGraph,
    LoweringReceipt, LoweringRegistry, MatDbQueryEnv, POWER_PORT_SCHEMA_HASH_DOMAIN,
    ParameterProvenance, ParameterSource, PartConstructor, PartLoweringContext, PartLoweringKey,
    PowerPortSchemaCard, ProvenanceRecorder, ResolvedInterfaceCard, SourceLoweringReceipt,
    geometry_axial_area_over_length_m, geometry_generator_label, geometry_polar_inertia_over_rho,
    geometry_volume_m3, register_kernel_defaults, resolve_geometry_fingerprint,
    spiral_spring_stiffness_factor_m3,
};
pub use system::{
    AlgebraicLoopAudit, AnalysisError, AnalysisKind, AnalysisOutcome, AnalysisRequest, BlockClass,
    BlockCostEstimate, BlockPlan, BoundaryChannel, BoundaryFn, CAUSALIZATION_AUDIT_SCHEMA,
    CausalDetermination, CausalMatchingEntry, CausalizationAudit, Checkpoint, CheckpointError,
    CrossRateChannel, EXECUTION_PLAN_SCHEMA, EventBracket, EventConfig, EventDirection, EventError,
    EventReport, EventSetConfig, EventSetError, EventSetReport, EventSetStatus, EventStatus,
    ExecutionPlan, HarmonicResponse, LinearizeConfig, LocatedEvent, MAX_CROSS_RATE_ITERATIONS,
    MAX_MULTIRATE_SUBSTEPS, ModalAnalysis, ModalConfig, ModalMode, MultirateConfig, MultirateError,
    MultirateEventReport, MultiratePlan, MultirateStepReport, OutputChannel, PeriodicHbConfig,
    PeriodicHbSolution, PeriodicShootingConfig, PeriodicShootingSolution, Phase, RateGroup,
    RequiredTier, SensitivityConfig, SlowInterpolation, SolveConfig, StateMap,
    StateSpaceLinearization, StationaryConfig, StationarySolution, StepAdjointSensitivity,
    StepObjective, StepReport, StepTangentSensitivity, System, SystemBuilder, SystemError,
    SystemImage, TearStrategy, TimeScaleGroup, TimeScaleTier, TrajectoryAdjointResult,
    TrajectoryStepPartials, heat_sink,
};
pub use wiring::{AdmissionError, HeatSink, NodeId, PortRef, RouteTemperature};

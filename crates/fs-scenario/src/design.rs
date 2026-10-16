//! Design variables, manufacturing constraints, objectives, and failure modes.
//!
//! These are the declarations an optimization phase will consume, landed as
//! SCHEMA NOW so the 0.1 project format does not need a breaking change when
//! optimization arrives (bead `f85xj.17.8`; execution is `f85xj.6.14` scope).
//! Nothing here optimizes anything.
//!
//! Doctrine, all structural rather than advisory:
//! - **Nothing is implied.** An optimizer's target, its robustness posture,
//!   and the constraints it must respect are DECLARED. [`RobustnessPosture`]
//!   is a required field of [`OptimizationObjective`], so an objective with an
//!   unstated posture is unrepresentable rather than merely discouraged.
//! - **There is no unsourced constraint.** Every manufacturing constraint,
//!   objective, and failure mode carries a [`DesignSource`], mirroring the
//!   "no unsourced tolerance" rule this crate already applies to geometry.
//! - **An optimization request without manufacturing constraints REFUSES.**
//!   A design change nobody can build is not an optimum, so
//!   [`DesignBlock::optimization_request`] fails closed.
//! - **Bounds carry units.** Each variable kind and constraint kind declares
//!   its expected [`Dims`], and a mismatch is a named violation rather than a
//!   silent reinterpretation.
//!
//! Persistence is a CONTENT-ADDRESSED SIDECAR ([`ScenarioDesignExtension`]),
//! bound to the hash of the scenario's canonical IR. This is deliberate: the
//! scenario form is strictly positional with fixed arity, so adding a section
//! to it would move every pinned canonical hash and force an IR version bump.
//! The sidecar leaves `SCENARIO_IR_VERSION` untouched, exactly as the sensor
//! extension does.
//!
//! No-claim boundaries:
//! - A design variable NAMES what may change. It does not prove the named
//!   parameter exists in any geometry kernel, that a material card resolves,
//!   or that a fan is purchasable. Selections are recorded as citations
//!   ([`ContentHash`]), never as resolutions — this crate is L3 and cannot
//!   reach the material or catalogue layers that would resolve them.
//! - A manufacturing constraint is a DECLARATION with provenance, not a
//!   manufacturability proof. Satisfying every declared constraint does not
//!   mean a part is buildable; it means no declared rule was violated.
//! - A failure mode is a named threshold, not a physics model. Severity
//!   orders consequences; it does not predict occurrence.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher, hash_bytes};
use fs_qty::{Dims, QtyAny};

use crate::entity::{EntityCatalog, EntityRef, KindExpectation};
use crate::ir::write_ir;
use crate::scenario::{Scenario, Violation};

/// Schema version of the design sidecar wire format.
pub const DESIGN_EXTENSION_SCHEMA_VERSION: u16 = 1;

/// Domain separation for the design-extension identity.
pub const DESIGN_EXTENSION_IDENTITY_DOMAIN: &str = "org.frankensim.scenario.design-extension.v1";

/// Ledger artifact kind for a retained design sidecar.
pub const DESIGN_EXTENSION_ARTIFACT_KIND: &str = "scenario-design-extension-v1";

/// Upper bound on an admitted design-extension wire payload.
pub const MAX_DESIGN_EXTENSION_WIRE_BYTES: usize = 16 * 1024 * 1024;

const DESIGN_EXTENSION_MAGIC: &[u8; 8] = b"FSDSGN1\0";

/// Dimensionless, for discrete and count-valued declarations.
const DIMENSIONLESS: Dims = Dims([0, 0, 0, 0, 0, 0]);

/// Length, for thickness/spacing/clearance constraints and position bounds.
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Temperature, for thermal failure thresholds.
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);

/// Bounds on the size of a declared design block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DesignBudget {
    /// Maximum declared design variables.
    pub max_variables: usize,
    /// Maximum declared manufacturing constraints.
    pub max_constraints: usize,
    /// Maximum declared objectives.
    pub max_objectives: usize,
    /// Maximum declared failure modes.
    pub max_failure_modes: usize,
    /// Maximum discrete choices in one variable domain.
    pub max_choices: usize,
    /// Maximum bytes in any declared free-text field.
    pub max_text_bytes: usize,
}

/// Default design-block bounds.
pub const DEFAULT_DESIGN_BUDGET: DesignBudget = DesignBudget {
    max_variables: 4_096,
    max_constraints: 4_096,
    max_objectives: 256,
    max_failure_modes: 4_096,
    max_choices: 4_096,
    max_text_bytes: 4_096,
};

impl Default for DesignBudget {
    fn default() -> Self {
        DEFAULT_DESIGN_BUDGET
    }
}

/// Where a declared design object came from. There is no unsourced constraint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DesignSource {
    /// A manufacturing standard or process specification clause.
    Standard {
        /// Clause identifier, e.g. `IPC-2221B §6.2`.
        clause: String,
    },
    /// A supplier or process capability document.
    ProcessCapability {
        /// Supplier or process owner.
        supplier: String,
        /// Document identity and revision.
        document: String,
    },
    /// A component datasheet limit.
    Datasheet {
        /// Part identity.
        part: String,
        /// Datasheet revision.
        revision: String,
    },
    /// An internal engineering policy.
    InternalPolicy {
        /// Policy identifier.
        policy: String,
    },
    /// An engineering assumption, with the rationale that justifies it.
    Assumed {
        /// Why this value was assumed.
        rationale: String,
    },
}

impl DesignSource {
    /// Stable label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        match self {
            Self::Standard { .. } => "standard",
            Self::ProcessCapability { .. } => "process-capability",
            Self::Datasheet { .. } => "datasheet",
            Self::InternalPolicy { .. } => "internal-policy",
            Self::Assumed { .. } => "assumed",
        }
    }

    /// The free-text fields that must be nonempty for this source.
    fn nonempty_fields(&self) -> Vec<(&'static str, &str)> {
        match self {
            Self::Standard { clause } => vec![("clause", clause.as_str())],
            Self::ProcessCapability { supplier, document } => {
                vec![("supplier", supplier.as_str()), ("document", document)]
            }
            Self::Datasheet { part, revision } => {
                vec![("part", part.as_str()), ("revision", revision)]
            }
            Self::InternalPolicy { policy } => vec![("policy", policy.as_str())],
            Self::Assumed { rationale } => vec![("rationale", rationale.as_str())],
        }
    }

    fn absorb(&self, hasher: &mut DomainHasher) {
        hasher.update(self.label().as_bytes());
        for (name, value) in self.nonempty_fields() {
            hasher.update(name.as_bytes());
            hasher.update(&(value.len() as u64).to_le_bytes());
            hasher.update(value.as_bytes());
        }
    }
}

/// What a design variable is allowed to change.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignVariableKind {
    /// A geometry parameter driven by a declared parameterization handle.
    GeometryParameter,
    /// A choice among declared material cards.
    MaterialSelection,
    /// A choice among declared fan operating points.
    FanSelection,
    /// A layout position offset.
    LayoutPosition,
}

impl DesignVariableKind {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::GeometryParameter => "geometry-parameter",
            Self::MaterialSelection => "material-selection",
            Self::FanSelection => "fan-selection",
            Self::LayoutPosition => "layout-position",
        }
    }

    /// Whether this kind varies over a continuous range or a discrete set.
    #[must_use]
    pub const fn is_discrete(self) -> bool {
        matches!(self, Self::MaterialSelection | Self::FanSelection)
    }

    /// The entity kinds this variable may target.
    #[must_use]
    pub const fn expectation(self) -> KindExpectation {
        match self {
            Self::MaterialSelection => KindExpectation::Domain,
            _ => KindExpectation::Any,
        }
    }

    /// Expected dimensions for a continuous domain of this kind.
    #[must_use]
    pub const fn expected_dims(self) -> Dims {
        match self {
            Self::GeometryParameter | Self::LayoutPosition => LENGTH,
            Self::MaterialSelection | Self::FanSelection => DIMENSIONLESS,
        }
    }
}

/// The set of values a design variable may take.
#[derive(Debug, Clone, PartialEq)]
pub enum VariableDomain {
    /// A closed continuous interval, in coherent SI base units.
    Continuous {
        /// Inclusive lower bound.
        lo: QtyAny,
        /// Inclusive upper bound.
        hi: QtyAny,
    },
    /// A discrete choice set, recorded as CITATIONS rather than resolutions.
    ///
    /// A choice is a content address of the card or catalogue entry it names.
    /// This crate is L3 and cannot reach the layers that would resolve one, so
    /// a citation here proves only that a specific identity was declared.
    Discrete {
        /// Declared choices, in caller order.
        choices: Vec<ContentHash>,
    },
}

impl VariableDomain {
    /// Stable label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        match self {
            Self::Continuous { .. } => "continuous",
            Self::Discrete { .. } => "discrete",
        }
    }
}

/// A quantity in the scenario that an optimizer may change.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignVariable {
    name: String,
    target: EntityRef,
    kind: DesignVariableKind,
    domain: VariableDomain,
    source: DesignSource,
    row: usize,
}

impl DesignVariable {
    /// Declare a design variable.
    ///
    /// Structural admission happens in [`DesignBlock::declare_variable`];
    /// dimensional, bound-order and provenance checks are reported by
    /// [`DesignBlock::validate`] so a caller gets the whole repair list at
    /// once rather than one fault per attempt.
    #[must_use]
    pub fn new(
        name: &str,
        target: EntityRef,
        kind: DesignVariableKind,
        domain: VariableDomain,
        source: DesignSource,
    ) -> Self {
        Self {
            name: name.to_string(),
            target,
            kind,
            domain,
            source,
            row: 0,
        }
    }

    /// Declared variable name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The entity this variable modifies.
    #[must_use]
    pub const fn target(&self) -> EntityRef {
        self.target
    }

    /// What kind of change this variable expresses.
    #[must_use]
    pub const fn kind(&self) -> DesignVariableKind {
        self.kind
    }

    /// The admitted value set.
    #[must_use]
    pub const fn domain(&self) -> &VariableDomain {
        &self.domain
    }

    /// Declared provenance.
    #[must_use]
    pub const fn source(&self) -> &DesignSource {
        &self.source
    }

    /// Row assigned by the owning block.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// A manufacturing rule an optimizer must respect and a human can audit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ManufacturingConstraintKind {
    /// Minimum wall thickness.
    MinimumWallThickness,
    /// Minimum spacing between fins or features.
    MinimumFeatureSpacing,
    /// Minimum assembly clearance.
    MinimumAssemblyClearance,
    /// Available stock thickness.
    StockThickness,
}

impl ManufacturingConstraintKind {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::MinimumWallThickness => "minimum-wall-thickness",
            Self::MinimumFeatureSpacing => "minimum-feature-spacing",
            Self::MinimumAssemblyClearance => "minimum-assembly-clearance",
            Self::StockThickness => "stock-thickness",
        }
    }

    /// Expected dimensions. Every current kind is a length.
    #[must_use]
    pub const fn expected_dims(self) -> Dims {
        LENGTH
    }
}

/// A declared, sourced manufacturing rule.
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturingConstraint {
    subject: EntityRef,
    kind: ManufacturingConstraintKind,
    magnitude: QtyAny,
    source: DesignSource,
    row: usize,
}

impl ManufacturingConstraint {
    /// Declare a manufacturing constraint.
    #[must_use]
    pub fn new(
        subject: EntityRef,
        kind: ManufacturingConstraintKind,
        magnitude: QtyAny,
        source: DesignSource,
    ) -> Self {
        Self {
            subject,
            kind,
            magnitude,
            source,
            row: 0,
        }
    }

    /// The entity this constraint applies to.
    #[must_use]
    pub const fn subject(&self) -> EntityRef {
        self.subject
    }

    /// Which rule this is.
    #[must_use]
    pub const fn kind(&self) -> ManufacturingConstraintKind {
        self.kind
    }

    /// The declared magnitude.
    #[must_use]
    pub const fn magnitude(&self) -> QtyAny {
        self.magnitude
    }

    /// Declared provenance.
    #[must_use]
    pub const fn source(&self) -> &DesignSource {
        &self.source
    }

    /// Row assigned by the owning block.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// Which direction an objective improves in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectiveSense {
    /// Smaller is better.
    Minimize,
    /// Larger is better.
    Maximize,
}

impl ObjectiveSense {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Minimize => "minimize",
            Self::Maximize => "maximize",
        }
    }
}

/// How an objective treats uncertainty. Required, never inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RobustnessPosture {
    /// Optimize the nominal value only, with no robustness claim.
    Nominal,
    /// Optimize a conditional value at risk at the declared tail level.
    ConditionalValueAtRisk {
        /// Tail level in `(0, 1)`; larger is more tail-averse.
        beta: f64,
    },
    /// Optimize the worst case over the declared operating envelope.
    WorstCase,
}

impl RobustnessPosture {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Nominal => "nominal",
            Self::ConditionalValueAtRisk { .. } => "cvar",
            Self::WorstCase => "worst-case",
        }
    }

    /// Whether this posture makes any robustness claim at all.
    #[must_use]
    pub const fn is_robust(self) -> bool {
        !matches!(self, Self::Nominal)
    }
}

/// A declared optimization target.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationObjective {
    qoi: String,
    sense: ObjectiveSense,
    posture: RobustnessPosture,
    source: DesignSource,
    row: usize,
}

impl OptimizationObjective {
    /// Declare an objective. The robustness posture is a required argument, so
    /// an objective with an unstated posture cannot be constructed.
    #[must_use]
    pub fn new(
        qoi: &str,
        sense: ObjectiveSense,
        posture: RobustnessPosture,
        source: DesignSource,
    ) -> Self {
        Self {
            qoi: qoi.to_string(),
            sense,
            posture,
            source,
            row: 0,
        }
    }

    /// The quantity of interest being optimized, by declared name.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Optimization direction.
    #[must_use]
    pub const fn sense(&self) -> ObjectiveSense {
        self.sense
    }

    /// Declared robustness posture.
    #[must_use]
    pub const fn posture(&self) -> RobustnessPosture {
        self.posture
    }

    /// Declared provenance.
    #[must_use]
    pub const fn source(&self) -> &DesignSource {
        &self.source
    }

    /// Row assigned by the owning block.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// Engineering consequence of reaching a failure mode. Ordered weakest first.
///
/// Mirrors the requirement severity classes used by the project layer. This
/// crate is L3 and cannot depend on that layer, so the ordering is declared
/// here and the higher layer maps onto it rather than the reverse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DesignSeverity {
    /// Preferred derating or reliability margin.
    ReliabilityDerating,
    /// Component or assembly damage limit.
    DamageLimit,
    /// Safety or regulatory limit.
    SafetyCritical,
}

impl DesignSeverity {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::ReliabilityDerating => "reliability-derating",
            Self::DamageLimit => "damage-limit",
            Self::SafetyCritical => "safety-critical",
        }
    }
}

/// A named failure threshold for a component class.
#[derive(Debug, Clone, PartialEq)]
pub struct FailureMode {
    name: String,
    component_class: String,
    severity: DesignSeverity,
    threshold: QtyAny,
    source: DesignSource,
    row: usize,
}

impl FailureMode {
    /// Declare a failure mode.
    #[must_use]
    pub fn new(
        name: &str,
        component_class: &str,
        severity: DesignSeverity,
        threshold: QtyAny,
        source: DesignSource,
    ) -> Self {
        Self {
            name: name.to_string(),
            component_class: component_class.to_string(),
            severity,
            threshold,
            source,
            row: 0,
        }
    }

    /// Declared failure-mode name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Component class this mode applies to.
    #[must_use]
    pub fn component_class(&self) -> &str {
        &self.component_class
    }

    /// Consequence class.
    #[must_use]
    pub const fn severity(&self) -> DesignSeverity {
        self.severity
    }

    /// The declared threshold.
    #[must_use]
    pub const fn threshold(&self) -> QtyAny {
        self.threshold
    }

    /// Declared provenance.
    #[must_use]
    pub const fn source(&self) -> &DesignSource {
        &self.source
    }

    /// Row assigned by the owning block.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// Refusal while declaring into a design block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DesignError {
    /// A collection reached its declared bound.
    CapacityExceeded {
        /// Which collection.
        resource: &'static str,
        /// Requested size.
        requested: usize,
        /// Admitted bound.
        limit: usize,
    },
    /// The process could not reserve memory for the declaration.
    AllocationRefused {
        /// Which collection.
        resource: &'static str,
    },
}

impl fmt::Display for DesignError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CapacityExceeded {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "design {resource} capacity exceeded: requested {requested}, limit {limit}"
            ),
            Self::AllocationRefused { resource } => {
                write!(formatter, "could not reserve memory for design {resource}")
            }
        }
    }
}

impl std::error::Error for DesignError {}

impl DesignError {
    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::CapacityExceeded { .. } => "design-capacity-exceeded",
            Self::AllocationRefused { .. } => "design-allocation-refused",
        }
    }
}

/// Why an optimization request was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OptimizationRefusal {
    /// No objective was declared, so there is nothing to optimize.
    NoObjective,
    /// No manufacturing constraint was declared.
    ///
    /// A design change nobody can build is not an optimum. This refusal is the
    /// no-silent-defaults doctrine applied to optimization: the absence of
    /// constraints is treated as an unanswered question, never as "none apply".
    NoManufacturingConstraint,
    /// No design variable was declared, so nothing may change.
    NoDesignVariable,
}

impl fmt::Display for OptimizationRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoObjective => write!(
                formatter,
                "no objective is declared; an optimization target is never inferred"
            ),
            Self::NoManufacturingConstraint => write!(
                formatter,
                "no manufacturing constraint is declared; a design change nobody can build is not \
                 an optimum, so an unconstrained request is refused rather than treated as \
                 unconstrained-by-choice"
            ),
            Self::NoDesignVariable => write!(
                formatter,
                "no design variable is declared; nothing is permitted to change"
            ),
        }
    }
}

impl OptimizationRefusal {
    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoObjective => "design-optimization-no-objective",
            Self::NoManufacturingConstraint => "design-optimization-unconstrained",
            Self::NoDesignVariable => "design-optimization-no-variable",
        }
    }

    /// How to repair the request.
    #[must_use]
    pub fn fix(&self) -> String {
        match self {
            Self::NoObjective => {
                "declare an OptimizationObjective naming the QoI, sense, and robustness posture"
                    .to_string()
            }
            Self::NoManufacturingConstraint => {
                "declare the manufacturing constraints the design must respect, each with its \
                 DesignSource; if none genuinely apply, that is itself a sourced declaration"
                    .to_string()
            }
            Self::NoDesignVariable => {
                "declare at least one DesignVariable naming what may change".to_string()
            }
        }
    }
}

/// An admitted optimization request: objectives, variables, and the
/// constraints they must respect.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationRequest<'a> {
    objectives: &'a [OptimizationObjective],
    variables: &'a [DesignVariable],
    constraints: &'a [ManufacturingConstraint],
}

impl<'a> OptimizationRequest<'a> {
    /// Declared objectives.
    #[must_use]
    pub const fn objectives(&self) -> &'a [OptimizationObjective] {
        self.objectives
    }

    /// Declared design variables.
    #[must_use]
    pub const fn variables(&self) -> &'a [DesignVariable] {
        self.variables
    }

    /// Declared manufacturing constraints.
    #[must_use]
    pub const fn constraints(&self) -> &'a [ManufacturingConstraint] {
        self.constraints
    }

    /// Whether every objective makes a robustness claim.
    #[must_use]
    pub fn fully_robust(&self) -> bool {
        self.objectives
            .iter()
            .all(|objective| objective.posture().is_robust())
    }
}

/// The declared design block for one scenario.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignBlock {
    variables: Vec<DesignVariable>,
    constraints: Vec<ManufacturingConstraint>,
    objectives: Vec<OptimizationObjective>,
    failure_modes: Vec<FailureMode>,
    budget: DesignBudget,
}

impl Default for DesignBlock {
    fn default() -> Self {
        Self::new()
    }
}

fn push_checked<T>(
    values: &mut Vec<T>,
    value: T,
    resource: &'static str,
) -> Result<(), DesignError> {
    values
        .try_reserve(1)
        .map_err(|_| DesignError::AllocationRefused { resource })?;
    values.push(value);
    Ok(())
}

impl DesignBlock {
    /// An empty block under the default budget.
    #[must_use]
    pub fn new() -> Self {
        Self::with_budget(DEFAULT_DESIGN_BUDGET)
    }

    /// An empty block under a caller-declared budget.
    #[must_use]
    pub fn with_budget(budget: DesignBudget) -> Self {
        Self {
            variables: Vec::new(),
            constraints: Vec::new(),
            objectives: Vec::new(),
            failure_modes: Vec::new(),
            budget,
        }
    }

    /// The declared budget.
    #[must_use]
    pub const fn budget(&self) -> DesignBudget {
        self.budget
    }

    /// Declared design variables.
    #[must_use]
    pub fn variables(&self) -> &[DesignVariable] {
        &self.variables
    }

    /// Declared manufacturing constraints.
    #[must_use]
    pub fn constraints(&self) -> &[ManufacturingConstraint] {
        &self.constraints
    }

    /// Declared objectives.
    #[must_use]
    pub fn objectives(&self) -> &[OptimizationObjective] {
        &self.objectives
    }

    /// Declared failure modes.
    #[must_use]
    pub fn failure_modes(&self) -> &[FailureMode] {
        &self.failure_modes
    }

    /// Declare a design variable, returning its row.
    ///
    /// # Errors
    ///
    /// [`DesignError::CapacityExceeded`] beyond the budget, and
    /// [`DesignError::AllocationRefused`] when memory cannot be reserved.
    pub fn declare_variable(&mut self, variable: DesignVariable) -> Result<usize, DesignError> {
        if self.variables.len() >= self.budget.max_variables {
            return Err(DesignError::CapacityExceeded {
                resource: "variables",
                requested: self.variables.len() + 1,
                limit: self.budget.max_variables,
            });
        }
        let row = self.variables.len();
        let mut variable = variable;
        variable.row = row;
        push_checked(&mut self.variables, variable, "variables")?;
        Ok(row)
    }

    /// Declare a manufacturing constraint, returning its row.
    ///
    /// # Errors
    ///
    /// As [`DesignBlock::declare_variable`].
    pub fn declare_constraint(
        &mut self,
        constraint: ManufacturingConstraint,
    ) -> Result<usize, DesignError> {
        if self.constraints.len() >= self.budget.max_constraints {
            return Err(DesignError::CapacityExceeded {
                resource: "constraints",
                requested: self.constraints.len() + 1,
                limit: self.budget.max_constraints,
            });
        }
        let row = self.constraints.len();
        let mut constraint = constraint;
        constraint.row = row;
        push_checked(&mut self.constraints, constraint, "constraints")?;
        Ok(row)
    }

    /// Declare an objective, returning its row.
    ///
    /// # Errors
    ///
    /// As [`DesignBlock::declare_variable`].
    pub fn declare_objective(
        &mut self,
        objective: OptimizationObjective,
    ) -> Result<usize, DesignError> {
        if self.objectives.len() >= self.budget.max_objectives {
            return Err(DesignError::CapacityExceeded {
                resource: "objectives",
                requested: self.objectives.len() + 1,
                limit: self.budget.max_objectives,
            });
        }
        let row = self.objectives.len();
        let mut objective = objective;
        objective.row = row;
        push_checked(&mut self.objectives, objective, "objectives")?;
        Ok(row)
    }

    /// Declare a failure mode, returning its row.
    ///
    /// # Errors
    ///
    /// As [`DesignBlock::declare_variable`].
    pub fn declare_failure_mode(&mut self, mode: FailureMode) -> Result<usize, DesignError> {
        if self.failure_modes.len() >= self.budget.max_failure_modes {
            return Err(DesignError::CapacityExceeded {
                resource: "failure modes",
                requested: self.failure_modes.len() + 1,
                limit: self.budget.max_failure_modes,
            });
        }
        let row = self.failure_modes.len();
        let mut mode = mode;
        mode.row = row;
        push_checked(&mut self.failure_modes, mode, "failure modes")?;
        Ok(row)
    }

    /// Look up a failure mode by name and component class.
    #[must_use]
    pub fn failure_mode(&self, name: &str, component_class: &str) -> Option<&FailureMode> {
        self.failure_modes
            .iter()
            .find(|mode| mode.name == name && mode.component_class == component_class)
    }

    /// Admit an optimization request, or refuse with every reason.
    ///
    /// Refusal is total and lists all reasons, so a caller repairs the whole
    /// set rather than rediscovering one problem per attempt.
    ///
    /// # Errors
    ///
    /// [`OptimizationRefusal`] rows for a missing objective, missing
    /// manufacturing constraints, or no declared variable.
    pub fn optimization_request(
        &self,
    ) -> Result<OptimizationRequest<'_>, Vec<OptimizationRefusal>> {
        let mut refusals = Vec::new();
        if self.objectives.is_empty() {
            refusals.push(OptimizationRefusal::NoObjective);
        }
        if self.variables.is_empty() {
            refusals.push(OptimizationRefusal::NoDesignVariable);
        }
        if self.constraints.is_empty() {
            refusals.push(OptimizationRefusal::NoManufacturingConstraint);
        }
        if refusals.is_empty() {
            Ok(OptimizationRequest {
                objectives: &self.objectives,
                variables: &self.variables,
                constraints: &self.constraints,
            })
        } else {
            Err(refusals)
        }
    }

    /// Report every structural, dimensional and provenance fault at once.
    #[must_use]
    #[allow(clippy::too_many_lines)] // One linear pass per collection; splitting it would scatter the codes.
    pub fn validate(&self, catalog: &EntityCatalog) -> Vec<Violation> {
        let mut out = Vec::new();
        let limit = self.budget.max_text_bytes;

        let mut check_source = |source: &DesignSource, what: &str, out: &mut Vec<Violation>| {
            for (field, value) in source.nonempty_fields() {
                if value.trim().is_empty() {
                    out.push(Violation {
                        code: "design-source-empty",
                        what: format!("{what}: {} source field `{field}` is empty", source.label()),
                        fix: format!("give `{field}` the exact document, clause or rationale text"),
                    });
                }
                if value.len() > limit {
                    out.push(Violation {
                        code: "design-text-too-long",
                        what: format!(
                            "{what}: source field `{field}` is {} bytes, limit {limit}",
                            value.len()
                        ),
                        fix: "shorten the field to an identifier and cite the document separately"
                            .to_string(),
                    });
                }
            }
        };

        for variable in &self.variables {
            let what = format!("design variable row {} (`{}`)", variable.row, variable.name);
            if variable.name.trim().is_empty() {
                out.push(Violation {
                    code: "design-variable-name-empty",
                    what: format!("design variable row {} has an empty name", variable.row),
                    fix: "give the variable a stable declared name".to_string(),
                });
            }
            check_source(&variable.source, &what, &mut out);

            if !variable
                .kind
                .expectation()
                .admits(variable.target.target().kind())
            {
                out.push(Violation {
                    code: "design-variable-target-kind",
                    what: format!(
                        "{what}: a {} variable may not target a {} entity",
                        variable.kind.label(),
                        variable.target.target().kind().label()
                    ),
                    fix: format!(
                        "retarget the variable at a {} entity",
                        variable.kind.expectation().label()
                    ),
                });
            }
            if catalog.resolve(variable.target).is_err() {
                out.push(Violation {
                    code: "design-variable-dangling-target",
                    what: format!(
                        "{what}: target {} does not resolve in the catalog",
                        variable.target.target().short_token()
                    ),
                    fix: "declare the entity, or retarget the variable at a live one".to_string(),
                });
            }

            match &variable.domain {
                VariableDomain::Continuous { lo, hi } => {
                    if variable.kind.is_discrete() {
                        out.push(Violation {
                            code: "design-variable-domain-kind",
                            what: format!(
                                "{what}: a {} variable needs a discrete choice set",
                                variable.kind.label()
                            ),
                            fix: "declare the admitted choices as content-addressed citations"
                                .to_string(),
                        });
                    }
                    let expected = variable.kind.expected_dims();
                    for (label, qty) in [("lo", lo), ("hi", hi)] {
                        if qty.dims != expected {
                            out.push(Violation {
                                code: "design-variable-dims",
                                what: format!(
                                    "{what}: {label} bound has dims {:?}, expected {expected:?}",
                                    qty.dims
                                ),
                                fix: format!(
                                    "express the {label} bound in the units a {} requires",
                                    variable.kind.label()
                                ),
                            });
                        }
                        if !qty.value.is_finite() {
                            out.push(Violation {
                                code: "design-variable-nonfinite",
                                what: format!("{what}: {label} bound is not finite"),
                                fix: "declare a finite bound; an unbounded variable is not a \
                                      declaration"
                                    .to_string(),
                            });
                        }
                    }
                    if lo.value.is_finite() && hi.value.is_finite() && lo.value > hi.value {
                        out.push(Violation {
                            code: "design-variable-bound-order",
                            what: format!(
                                "{what}: lower bound {} exceeds upper bound {}",
                                lo.value, hi.value
                            ),
                            fix: "order the bounds so lo <= hi".to_string(),
                        });
                    }
                }
                VariableDomain::Discrete { choices } => {
                    if !variable.kind.is_discrete() {
                        out.push(Violation {
                            code: "design-variable-domain-kind",
                            what: format!(
                                "{what}: a {} variable needs a continuous range",
                                variable.kind.label()
                            ),
                            fix: "declare lo/hi bounds with the variable's units".to_string(),
                        });
                    }
                    if choices.is_empty() {
                        out.push(Violation {
                            code: "design-variable-domain-empty",
                            what: format!("{what}: the discrete choice set is empty"),
                            fix: "declare at least one admitted choice, or remove the variable"
                                .to_string(),
                        });
                    }
                    if choices.len() > self.budget.max_choices {
                        out.push(Violation {
                            code: "design-variable-choices-exceeded",
                            what: format!(
                                "{what}: {} choices exceed the budget of {}",
                                choices.len(),
                                self.budget.max_choices
                            ),
                            fix: "reduce the declared choice set".to_string(),
                        });
                    }
                    let mut seen: Vec<&ContentHash> = Vec::new();
                    for choice in choices {
                        if seen.contains(&choice) {
                            out.push(Violation {
                                code: "design-variable-choice-duplicate",
                                what: format!(
                                    "{what}: choice {} is declared twice",
                                    choice.to_hex()
                                ),
                                fix: "remove the duplicate citation".to_string(),
                            });
                        } else {
                            seen.push(choice);
                        }
                    }
                }
            }
        }

        for constraint in &self.constraints {
            let what = format!(
                "manufacturing constraint row {} ({})",
                constraint.row,
                constraint.kind.label()
            );
            check_source(&constraint.source, &what, &mut out);
            let expected = constraint.kind.expected_dims();
            if constraint.magnitude.dims != expected {
                out.push(Violation {
                    code: "design-constraint-dims",
                    what: format!(
                        "{what}: magnitude has dims {:?}, expected {expected:?}",
                        constraint.magnitude.dims
                    ),
                    fix: "express the constraint in coherent SI base units".to_string(),
                });
            }
            if !constraint.magnitude.value.is_finite() || constraint.magnitude.value <= 0.0 {
                out.push(Violation {
                    code: "design-constraint-magnitude",
                    what: format!(
                        "{what}: magnitude {} is not a positive finite length",
                        constraint.magnitude.value
                    ),
                    fix: "declare a positive finite magnitude".to_string(),
                });
            }
            if catalog.resolve(constraint.subject).is_err() {
                out.push(Violation {
                    code: "design-constraint-dangling-subject",
                    what: format!(
                        "{what}: subject {} does not resolve in the catalog",
                        constraint.subject.target().short_token()
                    ),
                    fix: "declare the entity, or retarget the constraint".to_string(),
                });
            }
        }

        for objective in &self.objectives {
            let what = format!("objective row {} (`{}`)", objective.row, objective.qoi);
            check_source(&objective.source, &what, &mut out);
            if objective.qoi.trim().is_empty() {
                out.push(Violation {
                    code: "design-objective-qoi-empty",
                    what: format!("objective row {} names no QoI", objective.row),
                    fix: "name the quantity of interest this objective optimizes".to_string(),
                });
            }
            if let RobustnessPosture::ConditionalValueAtRisk { beta } = objective.posture
                && !(beta.is_finite() && beta > 0.0 && beta < 1.0)
            {
                out.push(Violation {
                    code: "design-objective-cvar-beta",
                    what: format!("{what}: CVaR beta {beta} is not in the open interval (0, 1)"),
                    fix: "declare a tail level strictly between 0 and 1".to_string(),
                });
            }
        }

        for mode in &self.failure_modes {
            let what = format!("failure mode row {} (`{}`)", mode.row, mode.name);
            check_source(&mode.source, &what, &mut out);
            if mode.name.trim().is_empty() {
                out.push(Violation {
                    code: "design-failure-mode-name-empty",
                    what: format!("failure mode row {} has an empty name", mode.row),
                    fix: "give the failure mode a stable declared name".to_string(),
                });
            }
            if mode.component_class.trim().is_empty() {
                out.push(Violation {
                    code: "design-failure-mode-class-empty",
                    what: format!("{what}: component class is empty"),
                    fix: "name the component class this mode applies to".to_string(),
                });
            }
            if !mode.threshold.value.is_finite() {
                out.push(Violation {
                    code: "design-failure-mode-threshold",
                    what: format!("{what}: threshold is not finite"),
                    fix: "declare a finite threshold".to_string(),
                });
            }
            if mode.threshold.dims != TEMPERATURE && mode.threshold.dims != DIMENSIONLESS {
                out.push(Violation {
                    code: "design-failure-mode-dims",
                    what: format!(
                        "{what}: threshold has dims {:?}; thermal modes are K and proxy modes are \
                         dimensionless",
                        mode.threshold.dims
                    ),
                    fix: "express the threshold in K, or as a dimensionless proxy".to_string(),
                });
            }
        }

        out
    }
}

/// A content-addressed design sidecar bound to one scenario.
///
/// Bound to the hash of the scenario's canonical IR rather than folded into
/// it, so `SCENARIO_IR_VERSION` and every pinned canonical hash stay unmoved.
#[derive(Debug, Clone, PartialEq)]
pub struct ScenarioDesignExtension {
    scenario_ir_hash: ContentHash,
    block: DesignBlock,
    identity: ContentHash,
}

impl ScenarioDesignExtension {
    /// Bind a design block to a scenario.
    #[must_use]
    pub fn new(scenario: &Scenario, block: DesignBlock) -> Self {
        let scenario_ir_hash = hash_bytes(write_ir(scenario).as_bytes());
        let identity = design_identity(scenario_ir_hash, &block);
        Self {
            scenario_ir_hash,
            block,
            identity,
        }
    }

    /// The bound scenario's canonical IR hash.
    #[must_use]
    pub const fn scenario_ir_hash(&self) -> ContentHash {
        self.scenario_ir_hash
    }

    /// The declared design block.
    #[must_use]
    pub const fn block(&self) -> &DesignBlock {
        &self.block
    }

    /// Domain-separated identity over the typed declaration fields.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Whether this sidecar belongs to the given scenario.
    ///
    /// A sidecar that does not verify is not evidence about that scenario, so
    /// callers must check rather than assume.
    #[must_use]
    pub fn verifies_scenario(&self, scenario: &Scenario) -> bool {
        hash_bytes(write_ir(scenario).as_bytes()) == self.scenario_ir_hash
    }
}

fn absorb_qty(hasher: &mut DomainHasher, qty: QtyAny) {
    hasher.update(&qty.value.to_bits().to_le_bytes());
    for exponent in qty.dims.0 {
        hasher.update(&exponent.to_le_bytes());
    }
}

fn absorb_text(hasher: &mut DomainHasher, value: &str) {
    hasher.update(&(value.len() as u64).to_le_bytes());
    hasher.update(value.as_bytes());
}

/// Domain-separated identity over a bound design block.
///
/// Binds the typed fields, never a rendered string, and length-prefixes every
/// text field so two different declarations cannot share a preimage.
#[must_use]
pub fn design_identity(scenario_ir_hash: ContentHash, block: &DesignBlock) -> ContentHash {
    let mut hasher = DomainHasher::new(DESIGN_EXTENSION_IDENTITY_DOMAIN);
    hasher.update(&DESIGN_EXTENSION_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&scenario_ir_hash.0);

    hasher.update(&(block.variables.len() as u64).to_le_bytes());
    for variable in &block.variables {
        absorb_text(&mut hasher, &variable.name);
        hasher.update(&variable.target.target().digest().0);
        hasher.update(variable.kind.label().as_bytes());
        hasher.update(variable.domain.label().as_bytes());
        match &variable.domain {
            VariableDomain::Continuous { lo, hi } => {
                absorb_qty(&mut hasher, *lo);
                absorb_qty(&mut hasher, *hi);
            }
            VariableDomain::Discrete { choices } => {
                hasher.update(&(choices.len() as u64).to_le_bytes());
                for choice in choices {
                    hasher.update(&choice.0);
                }
            }
        }
        variable.source.absorb(&mut hasher);
    }

    hasher.update(&(block.constraints.len() as u64).to_le_bytes());
    for constraint in &block.constraints {
        hasher.update(&constraint.subject.target().digest().0);
        hasher.update(constraint.kind.label().as_bytes());
        absorb_qty(&mut hasher, constraint.magnitude);
        constraint.source.absorb(&mut hasher);
    }

    hasher.update(&(block.objectives.len() as u64).to_le_bytes());
    for objective in &block.objectives {
        absorb_text(&mut hasher, &objective.qoi);
        hasher.update(objective.sense.label().as_bytes());
        hasher.update(objective.posture.label().as_bytes());
        if let RobustnessPosture::ConditionalValueAtRisk { beta } = objective.posture {
            hasher.update(&beta.to_bits().to_le_bytes());
        }
        objective.source.absorb(&mut hasher);
    }

    hasher.update(&(block.failure_modes.len() as u64).to_le_bytes());
    for mode in &block.failure_modes {
        absorb_text(&mut hasher, &mode.name);
        absorb_text(&mut hasher, &mode.component_class);
        hasher.update(mode.severity.label().as_bytes());
        absorb_qty(&mut hasher, mode.threshold);
        mode.source.absorb(&mut hasher);
    }

    hasher.finalize()
}

/// The sidecar's stable wire prefix, for callers writing their own envelopes.
#[must_use]
pub const fn design_extension_magic() -> &'static [u8; 8] {
    DESIGN_EXTENSION_MAGIC
}

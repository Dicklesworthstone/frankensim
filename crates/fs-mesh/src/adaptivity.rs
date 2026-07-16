//! Goal-oriented adaptivity accounting.
//!
//! This module is the dependency-clean seam between L3 error estimators and
//! L2 mesh evolution.  It does not import `fs-adjoint` or `fs-ir`: callers pass
//! their retained identities as opaque 32-byte digests.  Admission validates
//! the complete numerical accounting and publishes a canonical,
//! declaration-only receipt; it does not promote caller-supplied bounds into a
//! certificate.

use fs_exec::Cx;

/// Canonical schema name for [`AdaptivityReceipt::to_json`].
pub const ADAPTIVITY_RECEIPT_SCHEMA_V1: &str = "fs-mesh-adaptivity-receipt-v1";

/// Canonical schema name for [`ConservativeRemapReport::to_json`].
pub const CONSERVATIVE_REMAP_REPORT_SCHEMA_V1: &str = "fs-mesh-conservative-cell-remap-report-v1";

/// Maximum source cells in one conservative remap request.
pub const MAX_REMAP_SOURCE_CELLS: usize = 1_000_000;

/// Maximum target cells in one conservative remap request.
pub const MAX_REMAP_TARGET_CELLS: usize = 1_000_000;

/// Maximum source-to-target overlap contributions in one remap request.
pub const MAX_REMAP_CONTRIBUTIONS: usize = 4_000_000;

/// Maximum simultaneously retained auxiliary bytes for one remap request.
pub const MAX_REMAP_AUXILIARY_BYTES: usize = 256 * 1024 * 1024;

/// Largest admitted absolute defect in a source row's unity partition.
pub const MAX_REMAP_SOURCE_COVERAGE_TOLERANCE: f64 = 1.0e-6;

const ID_BYTES: usize = 32;
const REMAP_CANCELLATION_STRIDE: usize = 4_096;

/// Why adaptivity accounting could not be admitted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdaptivityError {
    /// A mesh-evolution record did not change the state identity.
    UnchangedMeshState,
    /// Before and after snapshots named different quantities of interest.
    QoiMismatch,
    /// The before snapshot did not name the lineage source state.
    BeforeStateMismatch,
    /// The after snapshot did not name the lineage target state.
    AfterStateMismatch,
    /// Physical-topology change was declared without a connectivity change.
    PhysicalTopologyWithoutConnectivity,
    /// An adaptivity step attempted to suppress its discontinuity flag without
    /// retained continuity evidence.
    GradientContinuityUnproven,
    /// Declared effects contradicted the selected action's fixed semantics.
    ActionEffectsMismatch {
        /// Action whose fixed semantics were contradicted.
        action: AdaptivityAction,
    },
    /// A non-negative finite accounting field was negative or non-finite.
    InvalidNonnegative {
        /// Stable field name.
        field: &'static str,
        /// Exact rejected floating-point representation.
        value_bits: u64,
    },
    /// A signed accounting field was non-finite.
    InvalidFinite {
        /// Stable field name.
        field: &'static str,
        /// Exact rejected floating-point representation.
        value_bits: u64,
    },
    /// Adding the admitted estimator and conversion bounds overflowed.
    QoiBoundOverflow,
}

impl core::fmt::Display for AdaptivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnchangedMeshState => {
                f.write_str("mesh evolution requires distinct source and target state IDs")
            }
            Self::QoiMismatch => {
                f.write_str("before and after snapshots must name the same declared QoI")
            }
            Self::BeforeStateMismatch => {
                f.write_str("before snapshot must name the lineage source state")
            }
            Self::AfterStateMismatch => {
                f.write_str("after snapshot must name the lineage target state")
            }
            Self::PhysicalTopologyWithoutConnectivity => f.write_str(
                "a declared physical-topology change also requires a connectivity change",
            ),
            Self::GradientContinuityUnproven => f.write_str(
                "adaptivity requires a gradient-discontinuity flag until continuity evidence exists",
            ),
            Self::ActionEffectsMismatch { action } => write!(
                f,
                "declared effects contradict the fixed semantics of {}",
                action.as_str()
            ),
            Self::InvalidNonnegative { field, value_bits } => write!(
                f,
                "{field} must be finite and non-negative (rejected bits {value_bits:#018x})"
            ),
            Self::InvalidFinite { field, value_bits } => write!(
                f,
                "{field} must be finite (rejected bits {value_bits:#018x})"
            ),
            Self::QoiBoundOverflow => {
                f.write_str("estimator plus conversion QoI bounds must remain finite")
            }
        }
    }
}

impl std::error::Error for AdaptivityError {}

fn hex(bytes: &[u8; ID_BYTES]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(ID_BYTES * 2);
    for byte in bytes {
        let byte = usize::from(*byte);
        out.push(char::from(DIGITS[byte >> 4]));
        out.push(char::from(DIGITS[byte & 0x0f]));
    }
    out
}

fn is_zero(value: f64) -> bool {
    value.to_bits() << 1 == 0
}

fn admit_nonnegative(field: &'static str, value: f64) -> Result<f64, AdaptivityError> {
    if !value.is_finite() || value < 0.0 {
        return Err(AdaptivityError::InvalidNonnegative {
            field,
            value_bits: value.to_bits(),
        });
    }
    Ok(if is_zero(value) { 0.0 } else { value })
}

fn admit_finite(field: &'static str, value: f64) -> Result<f64, AdaptivityError> {
    if !value.is_finite() {
        return Err(AdaptivityError::InvalidFinite {
            field,
            value_bits: value.to_bits(),
        });
    }
    Ok(if is_zero(value) { 0.0 } else { value })
}

/// Stable identity of the declared quantity of interest.
///
/// The bytes are an adapter boundary for the owning L3/L6 identity type; this
/// type assigns no semantic or certificate authority to them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QoiId([u8; ID_BYTES]);

impl QoiId {
    /// Retain exact digest bytes from the owning identity system. Parsing does
    /// not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// Stable identity of a mesh/discretization state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MeshStateId([u8; ID_BYTES]);

impl MeshStateId {
    /// Retain exact digest bytes from the owning identity system. Parsing does
    /// not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// Stable identity of an upstream QoI-error evidence artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QoiEvidenceId([u8; ID_BYTES]);

impl QoiEvidenceId {
    /// Retain exact digest bytes from the owning identity system. Parsing does
    /// not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// Stable identity of the upstream split/merge/remesh lineage record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LineageRecordId([u8; ID_BYTES]);

impl LineageRecordId {
    /// Retain exact digest bytes from the owning identity system. Parsing does
    /// not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// Stable identity of the declared conserved quantity, units, and balance
/// convention used by a remap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RemapInvariantId([u8; ID_BYTES]);

impl RemapInvariantId {
    /// Retain exact digest bytes from the owning quantity/state system.
    /// Parsing does not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// Stable identity of the upstream transfer/remap evidence artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RemapEvidenceId([u8; ID_BYTES]);

impl RemapEvidenceId {
    /// Retain exact digest bytes from the owning transfer implementation.
    /// Parsing does not add semantic authority.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; ID_BYTES]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact retained bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

/// The mesh/discretization action performed by an adaptivity step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptivityAction {
    /// Subdivide cells or elements.
    HRefine,
    /// Remove admitted spatial refinement.
    HCoarsen,
    /// Raise local approximation order.
    PEnrich,
    /// Lower local approximation order.
    PReduce,
    /// Reconnect or relocate a mesh under an anisotropic metric.
    AnisotropicRemesh,
    /// Relocate a moving mesh while preserving its declared connectivity.
    Untangle,
    /// Split one durable topological entity into multiple entities.
    Split,
    /// Merge multiple durable topological entities into one entity.
    Merge,
}

impl AdaptivityAction {
    const fn as_str(self) -> &'static str {
        match self {
            Self::HRefine => "h-refine",
            Self::HCoarsen => "h-coarsen",
            Self::PEnrich => "p-enrich",
            Self::PReduce => "p-reduce",
            Self::AnisotropicRemesh => "anisotropic-remesh",
            Self::Untangle => "untangle",
            Self::Split => "split",
            Self::Merge => "merge",
        }
    }
}

/// Caller-declared effects of an adaptivity action.
///
/// These values are retained rather than inferred from [`AdaptivityAction`]: a
/// remesh may relocate without reconnecting, a split may affect a discrete
/// patch without changing body topology, and an untangling solve may still be
/// nonsmooth. Receipt authority remains declaration-only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdaptivityEffects {
    connectivity_changed: bool,
    physical_topology_changed: bool,
    gradient_discontinuity: bool,
}

impl AdaptivityEffects {
    /// Admit declared action effects. A physical-topology change without any
    /// connectivity change is internally inconsistent and refuses.
    pub fn new(
        connectivity_changed: bool,
        physical_topology_changed: bool,
        gradient_discontinuity: bool,
    ) -> Result<Self, AdaptivityError> {
        if physical_topology_changed && !connectivity_changed {
            return Err(AdaptivityError::PhysicalTopologyWithoutConnectivity);
        }
        if !gradient_discontinuity {
            return Err(AdaptivityError::GradientContinuityUnproven);
        }
        Ok(Self {
            connectivity_changed,
            physical_topology_changed,
            gradient_discontinuity,
        })
    }

    /// Whether connectivity was declared changed.
    #[must_use]
    pub const fn connectivity_changed(&self) -> bool {
        self.connectivity_changed
    }

    /// Whether physical topology was declared changed.
    #[must_use]
    pub const fn physical_topology_changed(&self) -> bool {
        self.physical_topology_changed
    }

    /// Whether downstream gradients were declared discontinuous.
    #[must_use]
    pub const fn gradient_discontinuity(&self) -> bool {
        self.gradient_discontinuity
    }
}

/// The declared source of the refinement request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptivityTrigger {
    /// A declared QoI-error estimate requested the change.
    GoalOriented,
    /// A contact feature requested local resolution.
    Contact,
    /// A wear model requested local resolution or state transfer.
    Wear,
    /// A fracture model requested local resolution or a topology event.
    Fracture,
    /// Moving-domain mesh quality requested relocation or untangling.
    MovingMesh,
}

impl AdaptivityTrigger {
    const fn as_str(self) -> &'static str {
        match self {
            Self::GoalOriented => "goal-oriented",
            Self::Contact => "contact",
            Self::Wear => "wear",
            Self::Fracture => "fracture",
            Self::MovingMesh => "moving-mesh",
        }
    }
}

/// Retained source-to-target lineage for one adaptivity step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopologyLineage {
    action: AdaptivityAction,
    effects: AdaptivityEffects,
    source: MeshStateId,
    target: MeshStateId,
    record: LineageRecordId,
}

impl TopologyLineage {
    /// Admit a source-to-target transition backed by an upstream lineage
    /// record. A no-op state transition refuses rather than publishing a fake
    /// evolution event.
    pub fn new(
        action: AdaptivityAction,
        effects: AdaptivityEffects,
        source: MeshStateId,
        target: MeshStateId,
        record: LineageRecordId,
    ) -> Result<Self, AdaptivityError> {
        if source == target {
            return Err(AdaptivityError::UnchangedMeshState);
        }
        let effects_match_action = match action {
            AdaptivityAction::HRefine | AdaptivityAction::HCoarsen => {
                effects.connectivity_changed && !effects.physical_topology_changed
            }
            AdaptivityAction::PEnrich | AdaptivityAction::PReduce | AdaptivityAction::Untangle => {
                !effects.connectivity_changed && !effects.physical_topology_changed
            }
            AdaptivityAction::AnisotropicRemesh => true,
            AdaptivityAction::Split | AdaptivityAction::Merge => effects.connectivity_changed,
        };
        if !effects_match_action {
            return Err(AdaptivityError::ActionEffectsMismatch { action });
        }
        Ok(Self {
            action,
            effects,
            source,
            target,
            record,
        })
    }

    /// The declared mesh/discretization action.
    #[must_use]
    pub const fn action(&self) -> AdaptivityAction {
        self.action
    }

    /// Caller-declared effects retained for downstream routing.
    #[must_use]
    pub const fn effects(&self) -> AdaptivityEffects {
        self.effects
    }

    /// Source state identity.
    #[must_use]
    pub const fn source(&self) -> MeshStateId {
        self.source
    }

    /// Target state identity.
    #[must_use]
    pub const fn target(&self) -> MeshStateId {
        self.target
    }

    /// Upstream lineage-record identity.
    #[must_use]
    pub const fn record(&self) -> LineageRecordId {
        self.record
    }
}

/// QoI-error accounting at one mesh/discretization state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QoiBoundSnapshot {
    state: MeshStateId,
    qoi: QoiId,
    evidence: QoiEvidenceId,
    estimator_upper_bound: f64,
    conversion_upper_bound: f64,
    total_upper_bound: f64,
}

impl QoiBoundSnapshot {
    /// Admit non-negative finite estimator and representation-conversion error
    /// bounds. Their conservative sum must also remain finite.
    pub fn new(
        state: MeshStateId,
        qoi: QoiId,
        evidence: QoiEvidenceId,
        estimator_upper_bound: f64,
        conversion_upper_bound: f64,
    ) -> Result<Self, AdaptivityError> {
        let estimator_upper_bound =
            admit_nonnegative("estimator_upper_bound", estimator_upper_bound)?;
        let conversion_upper_bound =
            admit_nonnegative("conversion_upper_bound", conversion_upper_bound)?;
        let rounded_sum = estimator_upper_bound + conversion_upper_bound;
        if !rounded_sum.is_finite() {
            return Err(AdaptivityError::QoiBoundOverflow);
        }
        let (rounded_sum, residual) =
            fs_math::eft::two_sum(estimator_upper_bound, conversion_upper_bound);
        let total_upper_bound = if residual > 0.0 {
            fs_math::next_up(rounded_sum)
        } else {
            rounded_sum
        };
        if !total_upper_bound.is_finite() {
            return Err(AdaptivityError::QoiBoundOverflow);
        }
        Ok(Self {
            state,
            qoi,
            evidence,
            estimator_upper_bound,
            conversion_upper_bound,
            total_upper_bound,
        })
    }

    /// Mesh/discretization state to which this accounting applies.
    #[must_use]
    pub const fn state(&self) -> MeshStateId {
        self.state
    }

    /// Declared QoI identity.
    #[must_use]
    pub const fn qoi(&self) -> QoiId {
        self.qoi
    }

    /// Identity of the upstream evidence artifact that supplied the bounds.
    #[must_use]
    pub const fn evidence(&self) -> QoiEvidenceId {
        self.evidence
    }

    /// Estimator contribution to the declared QoI upper bound.
    #[must_use]
    pub const fn estimator_upper_bound(&self) -> f64 {
        self.estimator_upper_bound
    }

    /// Representation-conversion contribution to the declared QoI upper
    /// bound.
    #[must_use]
    pub const fn conversion_upper_bound(&self) -> f64 {
        self.conversion_upper_bound
    }

    /// Conservative sum of estimator and conversion contributions.
    #[must_use]
    pub const fn total_upper_bound(&self) -> f64 {
        self.total_upper_bound
    }
}

/// Whether a declared remap balance defect met its declared tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BalanceStatus {
    /// Absolute balance defect is at or below the declared tolerance.
    WithinDeclaredTolerance,
    /// Absolute balance defect exceeds the declared tolerance.
    ExceededDeclaredTolerance,
}

impl BalanceStatus {
    const fn as_str(self) -> &'static str {
        match self {
            Self::WithinDeclaredTolerance => "within-declared-tolerance",
            Self::ExceededDeclaredTolerance => "exceeded-declared-tolerance",
        }
    }
}

/// Field/internal-state transfer accounting for one adaptivity step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemapAccounting {
    invariant: RemapInvariantId,
    evidence: RemapEvidenceId,
    balance_defect: f64,
    balance_tolerance: f64,
    projection_error: f64,
    balance_status: BalanceStatus,
}

impl RemapAccounting {
    /// Admit a signed balance defect, a non-negative balance tolerance, and a
    /// non-negative projection-error upper bound, retaining both the declared
    /// invariant and the upstream transfer-evidence identity.
    pub fn new(
        invariant: RemapInvariantId,
        evidence: RemapEvidenceId,
        balance_defect: f64,
        balance_tolerance: f64,
        projection_error: f64,
    ) -> Result<Self, AdaptivityError> {
        let balance_defect = admit_finite("balance_defect", balance_defect)?;
        let balance_tolerance = admit_nonnegative("balance_tolerance", balance_tolerance)?;
        let projection_error = admit_nonnegative("projection_error", projection_error)?;
        let balance_status = if balance_defect.abs() <= balance_tolerance {
            BalanceStatus::WithinDeclaredTolerance
        } else {
            BalanceStatus::ExceededDeclaredTolerance
        };
        Ok(Self {
            invariant,
            evidence,
            balance_defect,
            balance_tolerance,
            projection_error,
            balance_status,
        })
    }

    /// Declared conserved quantity, units, and balance convention.
    #[must_use]
    pub const fn invariant(&self) -> RemapInvariantId {
        self.invariant
    }

    /// Identity of the upstream transfer/remap evidence artifact.
    #[must_use]
    pub const fn evidence(&self) -> RemapEvidenceId {
        self.evidence
    }

    /// Signed declared balance defect after transfer.
    #[must_use]
    pub const fn balance_defect(&self) -> f64 {
        self.balance_defect
    }

    /// Declared tolerance used to classify the balance defect.
    #[must_use]
    pub const fn balance_tolerance(&self) -> f64 {
        self.balance_tolerance
    }

    /// Declared projection-error upper bound.
    #[must_use]
    pub const fn projection_error(&self) -> f64 {
        self.projection_error
    }

    /// Classification against the declared balance tolerance.
    #[must_use]
    pub const fn balance_status(&self) -> BalanceStatus {
        self.balance_status
    }
}

/// One canonical source-cell to target-cell overlap fraction.
///
/// Descriptors are validated by [`conservative_cell_remap`]. The caller must
/// provide them in strictly increasing `(source, target)` order, with every
/// source represented and no duplicate pair. Requiring canonical input avoids
/// an uninterruptible allocation-bearing sort in this cancellation-sensitive
/// kernel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemapContribution {
    source: usize,
    target: usize,
    fraction: f64,
}

impl RemapContribution {
    /// Describe one source-to-target overlap fraction. Admission occurs only
    /// when the complete request reaches [`conservative_cell_remap`].
    #[must_use]
    pub const fn new(source: usize, target: usize, fraction: f64) -> Self {
        Self {
            source,
            target,
            fraction,
        }
    }

    /// Source-cell ordinal.
    #[must_use]
    pub const fn source(&self) -> usize {
        self.source
    }

    /// Target-cell ordinal.
    #[must_use]
    pub const fn target(&self) -> usize {
        self.target
    }

    /// Fraction of the source-cell extensive value sent to the target.
    #[must_use]
    pub const fn fraction(&self) -> f64 {
        self.fraction
    }
}

/// Fail-closed refusal from the piecewise-constant conservative remap kernel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConservativeRemapError {
    /// No source-cell extensive values were supplied.
    EmptySourceValues,
    /// No target cells were requested.
    EmptyTargetSet,
    /// No overlap contributions were supplied.
    EmptyContributions,
    /// A request dimension exceeded its static ceiling.
    TooMany {
        /// Stable dimension name.
        field: &'static str,
        /// Requested count.
        count: usize,
        /// Static ceiling.
        cap: usize,
    },
    /// The conservative auxiliary-memory estimate overflowed or exceeded its
    /// static ceiling.
    AuxiliaryMemoryExceeded {
        /// Estimated bytes, or `usize::MAX` when the estimate overflowed.
        bytes: usize,
        /// Static byte ceiling.
        cap: usize,
    },
    /// A tolerance was negative, non-finite, or outside its admitted domain.
    InvalidTolerance {
        /// Stable tolerance field name.
        field: &'static str,
        /// Exact rejected representation.
        value_bits: u64,
    },
    /// A source-cell extensive value was non-finite.
    NonFiniteSourceValue {
        /// Source-cell ordinal.
        source: usize,
        /// Exact rejected representation.
        value_bits: u64,
    },
    /// An overlap fraction was not finite and in `(0, 1]`.
    InvalidFraction {
        /// Canonical contribution ordinal.
        contribution: usize,
        /// Exact rejected representation.
        value_bits: u64,
    },
    /// A contribution named a missing source cell.
    SourceOutOfRange {
        /// Canonical contribution ordinal.
        contribution: usize,
        /// Rejected source ordinal.
        source: usize,
        /// Number of source cells.
        source_count: usize,
    },
    /// A contribution named a missing target cell.
    TargetOutOfRange {
        /// Canonical contribution ordinal.
        contribution: usize,
        /// Rejected target ordinal.
        target: usize,
        /// Number of target cells.
        target_count: usize,
    },
    /// Contributions were duplicated or not in canonical source/target order.
    NonCanonicalContributions {
        /// Canonical contribution ordinal at which ordering failed.
        contribution: usize,
        /// Previous `(source, target)` pair.
        previous: [usize; 2],
        /// Current `(source, target)` pair.
        current: [usize; 2],
    },
    /// A source cell had no overlap contribution.
    MissingSourceCoverage {
        /// First uncovered source-cell ordinal.
        source: usize,
    },
    /// A target cell had no incoming overlap contribution.
    MissingTargetCoverage {
        /// First uncovered target-cell ordinal.
        target: usize,
    },
    /// One source row did not partition unity within the declared tolerance.
    SourceCoverageExceeded {
        /// Source-cell ordinal.
        source: usize,
        /// Measured sum of its overlap fractions.
        coverage_bits: u64,
        /// Absolute measured distance from one.
        defect_bits: u64,
        /// Declared source-coverage tolerance.
        tolerance_bits: u64,
    },
    /// Finite source data overflowed a deterministic arithmetic stage.
    ArithmeticOverflow {
        /// Stable arithmetic stage.
        stage: &'static str,
        /// Source, target, or contribution ordinal for the stage.
        index: usize,
    },
    /// The final extensive balance exceeded its declared absolute tolerance.
    BalanceExceeded {
        /// Signed target-minus-source balance defect.
        defect_bits: u64,
        /// Declared absolute tolerance.
        tolerance_bits: u64,
    },
    /// A bounded output or scratch allocation was refused.
    Resource {
        /// Allocation stage.
        stage: &'static str,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

impl core::fmt::Display for ConservativeRemapError {
    #[allow(clippy::too_many_lines)] // every structured refusal stays actionable
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptySourceValues => f.write_str("conservative remap has no source cells"),
            Self::EmptyTargetSet => f.write_str("conservative remap has no target cells"),
            Self::EmptyContributions => {
                f.write_str("conservative remap has no overlap contributions")
            }
            Self::TooMany { field, count, cap } => {
                write!(
                    f,
                    "conservative remap {field} count {count} exceeds cap {cap}"
                )
            }
            Self::AuxiliaryMemoryExceeded { bytes, cap } => write!(
                f,
                "conservative remap auxiliary-memory estimate {bytes} bytes exceeds cap {cap}"
            ),
            Self::InvalidTolerance { field, value_bits } => write!(
                f,
                "conservative remap {field} is outside its finite non-negative domain (bits {value_bits:#018x})"
            ),
            Self::NonFiniteSourceValue { source, value_bits } => write!(
                f,
                "conservative remap source cell {source} has non-finite value bits {value_bits:#018x}"
            ),
            Self::InvalidFraction {
                contribution,
                value_bits,
            } => write!(
                f,
                "conservative remap contribution {contribution} fraction is not finite and in (0, 1] (bits {value_bits:#018x})"
            ),
            Self::SourceOutOfRange {
                contribution,
                source,
                source_count,
            } => write!(
                f,
                "conservative remap contribution {contribution} names source {source}, but source count is {source_count}"
            ),
            Self::TargetOutOfRange {
                contribution,
                target,
                target_count,
            } => write!(
                f,
                "conservative remap contribution {contribution} names target {target}, but target count is {target_count}"
            ),
            Self::NonCanonicalContributions {
                contribution,
                previous,
                current,
            } => write!(
                f,
                "conservative remap contribution {contribution} pair [{}, {}] does not follow [{}, {}] in strict canonical order",
                current[0], current[1], previous[0], previous[1]
            ),
            Self::MissingSourceCoverage { source } => write!(
                f,
                "conservative remap source cell {source} has no overlap contribution"
            ),
            Self::MissingTargetCoverage { target } => write!(
                f,
                "conservative remap target cell {target} has no incoming overlap contribution"
            ),
            Self::SourceCoverageExceeded {
                source,
                coverage_bits,
                defect_bits,
                tolerance_bits,
            } => write!(
                f,
                "conservative remap source cell {source} overlap sum {} has unity defect {} above tolerance {}",
                f64::from_bits(*coverage_bits),
                f64::from_bits(*defect_bits),
                f64::from_bits(*tolerance_bits)
            ),
            Self::ArithmeticOverflow { stage, index } => write!(
                f,
                "conservative remap arithmetic became non-finite at {stage} index {index}"
            ),
            Self::BalanceExceeded {
                defect_bits,
                tolerance_bits,
            } => write!(
                f,
                "conservative remap target-minus-source balance defect {} exceeds tolerance {}",
                f64::from_bits(*defect_bits),
                f64::from_bits(*tolerance_bits)
            ),
            Self::Resource { stage } => {
                write!(f, "conservative remap allocation refused at {stage}")
            }
            Self::Cancelled => f.write_str("conservative remap cancelled before publication"),
        }
    }
}

impl std::error::Error for ConservativeRemapError {}

#[derive(Debug, Clone, Copy, Default)]
struct DeterministicSum {
    value: f64,
    correction: f64,
}

impl DeterministicSum {
    fn add(&mut self, addend: f64) -> bool {
        let (value, residual) = fs_math::eft::two_sum(self.value, addend);
        let (correction, carry) = fs_math::eft::two_sum(self.correction, residual);
        let (value, residual) = fs_math::eft::two_sum(value, correction);
        let correction = residual + carry;
        if !value.is_finite() || !correction.is_finite() {
            return false;
        }
        self.value = value;
        self.correction = correction;
        true
    }

    fn finish(self) -> Option<f64> {
        let (value, residual) = fs_math::eft::two_sum(self.value, self.correction);
        let value = value + residual;
        value
            .is_finite()
            .then_some(if is_zero(value) { 0.0 } else { value })
    }
}

fn remap_tolerance(
    field: &'static str,
    value: f64,
    maximum: Option<f64>,
) -> Result<f64, ConservativeRemapError> {
    if !value.is_finite() || value < 0.0 || maximum.is_some_and(|maximum| value > maximum) {
        return Err(ConservativeRemapError::InvalidTolerance {
            field,
            value_bits: value.to_bits(),
        });
    }
    Ok(if is_zero(value) { 0.0 } else { value })
}

fn finish_source_coverage(
    source: usize,
    coverage: DeterministicSum,
    tolerance: f64,
) -> Result<f64, ConservativeRemapError> {
    let coverage = coverage
        .finish()
        .ok_or(ConservativeRemapError::ArithmeticOverflow {
            stage: "source-coverage-sum",
            index: source,
        })?;
    let defect = (coverage - 1.0).abs();
    if !defect.is_finite() {
        return Err(ConservativeRemapError::ArithmeticOverflow {
            stage: "source-coverage-defect",
            index: source,
        });
    }
    if defect > tolerance {
        return Err(ConservativeRemapError::SourceCoverageExceeded {
            source,
            coverage_bits: coverage.to_bits(),
            defect_bits: defect.to_bits(),
            tolerance_bits: tolerance.to_bits(),
        });
    }
    Ok(defect)
}

/// Measured accounting for one successful extensive cell remap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConservativeRemapReport {
    source_cells: usize,
    target_cells: usize,
    contributions: usize,
    source_coverage_tolerance: f64,
    maximum_source_coverage_defect: f64,
    balance_tolerance: f64,
    source_total: f64,
    target_total: f64,
    balance_defect: f64,
    preserves_nonnegative_inputs: bool,
}

impl ConservativeRemapReport {
    /// Number of source cells.
    #[must_use]
    pub const fn source_cells(&self) -> usize {
        self.source_cells
    }

    /// Number of target cells.
    #[must_use]
    pub const fn target_cells(&self) -> usize {
        self.target_cells
    }

    /// Number of admitted overlap contributions.
    #[must_use]
    pub const fn contributions(&self) -> usize {
        self.contributions
    }

    /// Declared absolute tolerance for each source row's unity partition.
    #[must_use]
    pub const fn source_coverage_tolerance(&self) -> f64 {
        self.source_coverage_tolerance
    }

    /// Largest measured absolute source-row deviation from one.
    #[must_use]
    pub const fn maximum_source_coverage_defect(&self) -> f64 {
        self.maximum_source_coverage_defect
    }

    /// Declared absolute extensive-balance tolerance.
    #[must_use]
    pub const fn balance_tolerance(&self) -> f64 {
        self.balance_tolerance
    }

    /// Deterministic measured source extensive total.
    #[must_use]
    pub const fn source_total(&self) -> f64 {
        self.source_total
    }

    /// Deterministic measured target extensive total.
    #[must_use]
    pub const fn target_total(&self) -> f64 {
        self.target_total
    }

    /// Signed target-minus-source balance defect.
    #[must_use]
    pub const fn balance_defect(&self) -> f64 {
        self.balance_defect
    }

    /// Whether non-negative input values were observed and therefore every
    /// admitted non-negative weighted target remained non-negative.
    #[must_use]
    pub const fn preserves_nonnegative_inputs(&self) -> bool {
        self.preserves_nonnegative_inputs
    }

    /// Convert measured kernel accounting into the existing adaptivity seam.
    /// The caller supplies only the owning invariant/evidence identities; this
    /// method cannot promote the measured f64 result into a certificate.
    pub fn accounting(
        &self,
        invariant: RemapInvariantId,
        evidence: RemapEvidenceId,
        projection_error: f64,
    ) -> Result<RemapAccounting, AdaptivityError> {
        RemapAccounting::new(
            invariant,
            evidence,
            self.balance_defect,
            self.balance_tolerance,
            projection_error,
        )
    }

    /// Canonical measured report for an owning evidence artifact.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"{}\",\"authority\":\"measured-f64\",\
             \"source_cells\":{},\"target_cells\":{},\"contributions\":{},\
             \"source_coverage_tolerance\":{:.17e},\
             \"maximum_source_coverage_defect\":{:.17e},\
             \"balance_tolerance\":{:.17e},\"source_total\":{:.17e},\
             \"target_total\":{:.17e},\"balance_defect\":{:.17e},\
             \"preserves_nonnegative_inputs\":{},\
             \"no_claim\":\"overlap fractions, units, geometry, and continuum conservation are caller-owned\"}}",
            CONSERVATIVE_REMAP_REPORT_SCHEMA_V1,
            self.source_cells,
            self.target_cells,
            self.contributions,
            self.source_coverage_tolerance,
            self.maximum_source_coverage_defect,
            self.balance_tolerance,
            self.source_total,
            self.target_total,
            self.balance_defect,
            self.preserves_nonnegative_inputs,
        )
    }
}

/// Published target values and measured accounting from a complete remap.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservativeRemapOutcome {
    target_values: Vec<f64>,
    report: ConservativeRemapReport,
}

impl ConservativeRemapOutcome {
    /// Borrow target-cell extensive values in target ordinal order.
    #[must_use]
    pub fn target_values(&self) -> &[f64] {
        &self.target_values
    }

    /// Borrow the measured remap report.
    #[must_use]
    pub const fn report(&self) -> &ConservativeRemapReport {
        &self.report
    }

    /// Consume the outcome and return its target values.
    #[must_use]
    pub fn into_target_values(self) -> Vec<f64> {
        self.target_values
    }
}

fn remap_checkpoint(cx: &Cx<'_>) -> Result<(), ConservativeRemapError> {
    cx.checkpoint()
        .map_err(|_| ConservativeRemapError::Cancelled)
}

fn remap_auxiliary_bytes(target_count: usize) -> Result<usize, ConservativeRemapError> {
    let per_target = core::mem::size_of::<DeterministicSum>()
        .checked_add(core::mem::size_of::<bool>())
        .and_then(|bytes| bytes.checked_add(core::mem::size_of::<f64>()))
        .ok_or(ConservativeRemapError::AuxiliaryMemoryExceeded {
            bytes: usize::MAX,
            cap: MAX_REMAP_AUXILIARY_BYTES,
        })?;
    let bytes = target_count.checked_mul(per_target).ok_or(
        ConservativeRemapError::AuxiliaryMemoryExceeded {
            bytes: usize::MAX,
            cap: MAX_REMAP_AUXILIARY_BYTES,
        },
    )?;
    if bytes > MAX_REMAP_AUXILIARY_BYTES {
        return Err(ConservativeRemapError::AuxiliaryMemoryExceeded {
            bytes,
            cap: MAX_REMAP_AUXILIARY_BYTES,
        });
    }
    Ok(bytes)
}

/// Remap one piecewise-constant extensive scalar through canonical overlap
/// fractions.
///
/// Every source row must partition unity within `source_coverage_tolerance`,
/// every target must receive at least one contribution, and the measured final
/// target-minus-source total must fit `balance_tolerance`. Values and overlap
/// products remain ordinary f64 measurements; the report does not certify the
/// caller's geometry, units, overlap construction, or continuum conservation.
/// Cancellation is checked at bounded scan/publication strides and never
/// exposes a partially filled target vector.
#[allow(clippy::too_many_lines)] // staged admission and publication stay visibly transactional
pub fn conservative_cell_remap(
    source_values: &[f64],
    target_count: usize,
    contributions: &[RemapContribution],
    source_coverage_tolerance: f64,
    balance_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<ConservativeRemapOutcome, ConservativeRemapError> {
    remap_checkpoint(cx)?;
    if source_values.is_empty() {
        return Err(ConservativeRemapError::EmptySourceValues);
    }
    if target_count == 0 {
        return Err(ConservativeRemapError::EmptyTargetSet);
    }
    if contributions.is_empty() {
        return Err(ConservativeRemapError::EmptyContributions);
    }
    for (field, count, cap) in [
        ("source_cells", source_values.len(), MAX_REMAP_SOURCE_CELLS),
        ("target_cells", target_count, MAX_REMAP_TARGET_CELLS),
        (
            "contributions",
            contributions.len(),
            MAX_REMAP_CONTRIBUTIONS,
        ),
    ] {
        if count > cap {
            return Err(ConservativeRemapError::TooMany { field, count, cap });
        }
    }
    let source_coverage_tolerance = remap_tolerance(
        "source_coverage_tolerance",
        source_coverage_tolerance,
        Some(MAX_REMAP_SOURCE_COVERAGE_TOLERANCE),
    )?;
    let balance_tolerance = remap_tolerance("balance_tolerance", balance_tolerance, None)?;
    let _auxiliary_bytes = remap_auxiliary_bytes(target_count)?;

    let mut source_total = DeterministicSum::default();
    let mut preserves_nonnegative_inputs = true;
    for (source, value) in source_values.iter().copied().enumerate() {
        if source % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        if !value.is_finite() {
            return Err(ConservativeRemapError::NonFiniteSourceValue {
                source,
                value_bits: value.to_bits(),
            });
        }
        preserves_nonnegative_inputs &= value >= 0.0;
        if !source_total.add(value) {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "source-total",
                index: source,
            });
        }
    }
    let source_total = source_total
        .finish()
        .ok_or(ConservativeRemapError::ArithmeticOverflow {
            stage: "source-total-publication",
            index: source_values.len() - 1,
        })?;

    let mut target_sums = Vec::new();
    target_sums
        .try_reserve_exact(target_count)
        .map_err(|_| ConservativeRemapError::Resource {
            stage: "target-sums",
        })?;
    for target in 0..target_count {
        if target % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        target_sums.push(DeterministicSum::default());
    }
    let mut target_received = Vec::new();
    target_received
        .try_reserve_exact(target_count)
        .map_err(|_| ConservativeRemapError::Resource {
            stage: "target-coverage",
        })?;
    for target in 0..target_count {
        if target % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        target_received.push(false);
    }

    let mut previous_pair = None;
    let mut active_source = None;
    let mut active_coverage = DeterministicSum::default();
    let mut maximum_source_coverage_defect = 0.0f64;
    for (index, contribution) in contributions.iter().copied().enumerate() {
        if index % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        if contribution.source >= source_values.len() {
            return Err(ConservativeRemapError::SourceOutOfRange {
                contribution: index,
                source: contribution.source,
                source_count: source_values.len(),
            });
        }
        if contribution.target >= target_count {
            return Err(ConservativeRemapError::TargetOutOfRange {
                contribution: index,
                target: contribution.target,
                target_count,
            });
        }
        if !contribution.fraction.is_finite()
            || contribution.fraction <= 0.0
            || contribution.fraction > 1.0
        {
            return Err(ConservativeRemapError::InvalidFraction {
                contribution: index,
                value_bits: contribution.fraction.to_bits(),
            });
        }
        let pair = (contribution.source, contribution.target);
        if let Some(previous) = previous_pair
            && pair <= previous
        {
            return Err(ConservativeRemapError::NonCanonicalContributions {
                contribution: index,
                previous: [previous.0, previous.1],
                current: [pair.0, pair.1],
            });
        }
        previous_pair = Some(pair);

        match active_source {
            None => {
                if contribution.source != 0 {
                    return Err(ConservativeRemapError::MissingSourceCoverage { source: 0 });
                }
                active_source = Some(contribution.source);
            }
            Some(source) if source != contribution.source => {
                let defect =
                    finish_source_coverage(source, active_coverage, source_coverage_tolerance)?;
                maximum_source_coverage_defect = maximum_source_coverage_defect.max(defect);
                let expected = source + 1;
                if contribution.source != expected {
                    return Err(ConservativeRemapError::MissingSourceCoverage { source: expected });
                }
                active_source = Some(contribution.source);
                active_coverage = DeterministicSum::default();
            }
            Some(_) => {}
        }
        if !active_coverage.add(contribution.fraction) {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "source-coverage-accumulation",
                index,
            });
        }
        let value = source_values[contribution.source] * contribution.fraction;
        if !value.is_finite() {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "weighted-contribution",
                index,
            });
        }
        if !target_sums[contribution.target].add(value) {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "target-accumulation",
                index: contribution.target,
            });
        }
        target_received[contribution.target] = true;
    }

    let Some(last_source) = active_source else {
        return Err(ConservativeRemapError::EmptyContributions);
    };
    let defect = finish_source_coverage(last_source, active_coverage, source_coverage_tolerance)?;
    maximum_source_coverage_defect = maximum_source_coverage_defect.max(defect);
    if last_source + 1 != source_values.len() {
        return Err(ConservativeRemapError::MissingSourceCoverage {
            source: last_source + 1,
        });
    }
    for (target, &received) in target_received.iter().enumerate() {
        if target % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        if !received {
            return Err(ConservativeRemapError::MissingTargetCoverage { target });
        }
    }

    let mut target_values = Vec::new();
    target_values.try_reserve_exact(target_count).map_err(|_| {
        ConservativeRemapError::Resource {
            stage: "target-values",
        }
    })?;
    let mut target_total = DeterministicSum::default();
    for (target, sum) in target_sums.into_iter().enumerate() {
        if target % REMAP_CANCELLATION_STRIDE == 0 {
            remap_checkpoint(cx)?;
        }
        let value = sum
            .finish()
            .ok_or(ConservativeRemapError::ArithmeticOverflow {
                stage: "target-value-publication",
                index: target,
            })?;
        if preserves_nonnegative_inputs && value < 0.0 {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "nonnegative-target-publication",
                index: target,
            });
        }
        if !target_total.add(value) {
            return Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "target-total",
                index: target,
            });
        }
        target_values.push(value);
    }
    let target_total = target_total
        .finish()
        .ok_or(ConservativeRemapError::ArithmeticOverflow {
            stage: "target-total-publication",
            index: target_count - 1,
        })?;
    let mut balance = DeterministicSum::default();
    if !balance.add(target_total) || !balance.add(-source_total) {
        return Err(ConservativeRemapError::ArithmeticOverflow {
            stage: "balance-defect",
            index: 0,
        });
    }
    let balance_defect = balance
        .finish()
        .ok_or(ConservativeRemapError::ArithmeticOverflow {
            stage: "balance-defect-publication",
            index: 0,
        })?;
    if balance_defect.abs() > balance_tolerance {
        return Err(ConservativeRemapError::BalanceExceeded {
            defect_bits: balance_defect.to_bits(),
            tolerance_bits: balance_tolerance.to_bits(),
        });
    }
    remap_checkpoint(cx)?;

    Ok(ConservativeRemapOutcome {
        target_values,
        report: ConservativeRemapReport {
            source_cells: source_values.len(),
            target_cells: target_count,
            contributions: contributions.len(),
            source_coverage_tolerance,
            maximum_source_coverage_defect,
            balance_tolerance,
            source_total,
            target_total,
            balance_defect,
            preserves_nonnegative_inputs,
        },
    })
}

/// Strict before/after trend of the declared total QoI upper bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QoiBoundTrend {
    /// The after bound is strictly smaller than the before bound.
    Decreased,
    /// The after bound has the identical admitted `f64` value.
    Unchanged,
    /// The after bound is strictly larger than the before bound.
    Increased,
}

impl QoiBoundTrend {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Decreased => "decreased",
            Self::Unchanged => "unchanged",
            Self::Increased => "increased",
        }
    }
}

/// Authority carried by this accounting-only receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptivityReceiptAuthority {
    /// Inputs are retained and validated, but their scientific authority is
    /// neither recreated nor promoted by `fs-mesh`.
    DeclarationOnly,
}

impl AdaptivityReceiptAuthority {
    const fn as_str(self) -> &'static str {
        match self {
            Self::DeclarationOnly => "declaration-only",
        }
    }
}

/// Complete accounting receipt for one mesh/discretization change.
///
/// Construction is the only publication path. It binds the requested trigger,
/// retained topology lineage, before/after evidence identities, error-ledger
/// components, remap defects, and the strict QoI-bound trend. A failed
/// refinement therefore remains visible as [`QoiBoundTrend::Unchanged`] or
/// [`QoiBoundTrend::Increased`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptivityReceipt {
    trigger: AdaptivityTrigger,
    lineage: TopologyLineage,
    before: QoiBoundSnapshot,
    after: QoiBoundSnapshot,
    remap: RemapAccounting,
    qoi_trend: QoiBoundTrend,
}

impl AdaptivityReceipt {
    /// Admit one complete adaptivity accounting record.
    pub fn admit(
        trigger: AdaptivityTrigger,
        lineage: TopologyLineage,
        before: QoiBoundSnapshot,
        after: QoiBoundSnapshot,
        remap: RemapAccounting,
    ) -> Result<Self, AdaptivityError> {
        if before.qoi != after.qoi {
            return Err(AdaptivityError::QoiMismatch);
        }
        if before.state != lineage.source {
            return Err(AdaptivityError::BeforeStateMismatch);
        }
        if after.state != lineage.target {
            return Err(AdaptivityError::AfterStateMismatch);
        }
        let qoi_trend = if after.total_upper_bound < before.total_upper_bound {
            QoiBoundTrend::Decreased
        } else if after.total_upper_bound.to_bits() == before.total_upper_bound.to_bits() {
            QoiBoundTrend::Unchanged
        } else {
            QoiBoundTrend::Increased
        };
        Ok(Self {
            trigger,
            lineage,
            before,
            after,
            remap,
            qoi_trend,
        })
    }

    /// Declared source of the adaptivity request.
    #[must_use]
    pub const fn trigger(&self) -> AdaptivityTrigger {
        self.trigger
    }

    /// Retained source-to-target topology lineage.
    #[must_use]
    pub const fn lineage(&self) -> &TopologyLineage {
        &self.lineage
    }

    /// QoI-error accounting before the step.
    #[must_use]
    pub const fn before(&self) -> &QoiBoundSnapshot {
        &self.before
    }

    /// QoI-error accounting after the step.
    #[must_use]
    pub const fn after(&self) -> &QoiBoundSnapshot {
        &self.after
    }

    /// Field/internal-state transfer accounting.
    #[must_use]
    pub const fn remap(&self) -> &RemapAccounting {
        &self.remap
    }

    /// Strict before/after QoI-bound trend.
    #[must_use]
    pub const fn qoi_trend(&self) -> QoiBoundTrend {
        self.qoi_trend
    }

    /// Whether the requested QoI upper bound strictly decreased.
    #[must_use]
    pub const fn qoi_bound_decreased(&self) -> bool {
        matches!(self.qoi_trend, QoiBoundTrend::Decreased)
    }

    /// Receipt authority. This cannot be selected by the caller.
    #[must_use]
    pub const fn authority(&self) -> AdaptivityReceiptAuthority {
        AdaptivityReceiptAuthority::DeclarationOnly
    }

    /// Canonical JSON suitable for content-addressing by an owning ledger.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"{}\",\"authority\":\"{}\",\"trigger\":\"{}\",\
             \"action\":\"{}\",\"declared_connectivity_changed\":{},\
             \"declared_physical_topology_changed\":{},\
             \"declared_gradient_discontinuity\":{},\
             \"lineage_record_id\":\"{}\",\"source_mesh_state_id\":\"{}\",\
             \"target_mesh_state_id\":\"{}\",\"qoi_id\":\"{}\",\
             \"before_evidence_id\":\"{}\",\"before_estimator_upper_bound\":{:.17e},\
             \"before_conversion_upper_bound\":{:.17e},\"before_total_upper_bound\":{:.17e},\
             \"after_evidence_id\":\"{}\",\"after_estimator_upper_bound\":{:.17e},\
             \"after_conversion_upper_bound\":{:.17e},\"after_total_upper_bound\":{:.17e},\
             \"qoi_bound_trend\":\"{}\",\"qoi_bound_decreased\":{},\
             \"remap_invariant_id\":\"{}\",\"remap_evidence_id\":\"{}\",\
             \"balance_defect\":{:.17e},\"balance_tolerance\":{:.17e},\
             \"balance_status\":\"{}\",\"projection_error\":{:.17e}}}",
            ADAPTIVITY_RECEIPT_SCHEMA_V1,
            self.authority().as_str(),
            self.trigger.as_str(),
            self.lineage.action.as_str(),
            self.lineage.effects.connectivity_changed,
            self.lineage.effects.physical_topology_changed,
            self.lineage.effects.gradient_discontinuity,
            hex(self.lineage.record.as_bytes()),
            hex(self.lineage.source.as_bytes()),
            hex(self.lineage.target.as_bytes()),
            hex(self.before.qoi.as_bytes()),
            hex(self.before.evidence.as_bytes()),
            self.before.estimator_upper_bound,
            self.before.conversion_upper_bound,
            self.before.total_upper_bound,
            hex(self.after.evidence.as_bytes()),
            self.after.estimator_upper_bound,
            self.after.conversion_upper_bound,
            self.after.total_upper_bound,
            self.qoi_trend.as_str(),
            self.qoi_bound_decreased(),
            hex(self.remap.invariant.as_bytes()),
            hex(self.remap.evidence.as_bytes()),
            self.remap.balance_defect,
            self.remap.balance_tolerance,
            self.remap.balance_status.as_str(),
            self.remap.projection_error,
        )
    }
}

//! Goal-oriented adaptivity accounting.
//!
//! This module is the dependency-clean seam between L3 error estimators and
//! L2 mesh evolution.  It does not import `fs-adjoint` or `fs-ir`: callers pass
//! their retained identities as opaque 32-byte digests.  Admission validates
//! the complete numerical accounting and publishes a canonical,
//! declaration-only receipt; it does not promote caller-supplied bounds into a
//! certificate.

/// Canonical schema name for [`AdaptivityReceipt::to_json`].
pub const ADAPTIVITY_RECEIPT_SCHEMA_V1: &str = "fs-mesh-adaptivity-receipt-v1";

const ID_BYTES: usize = 32;

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

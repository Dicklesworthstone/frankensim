//! Bounded finite-population propagation through hierarchical manufacturing
//! structure and closed piecewise-quadratic clearance laws.
//!
//! This module complements the first-order correlated-stack lane. Correlation
//! and hierarchy are represented by an explicit weighted finite population;
//! nonlinear and mode-switching behavior is represented by a bounded,
//! callback-free piecewise-quadratic law. The result is descriptive binary64
//! evidence for the supplied population, not a distributional or reliability
//! claim.

use std::{collections::BTreeMap, num::NonZeroU64};

/// Receipt/evaluation schema for this fixed-order structured lane.
pub const STRUCTURED_PROPAGATION_SCHEMA_V1: u32 = 1;

/// Maximum number of hierarchy nodes retained by schema version one.
pub const MAX_STRUCTURED_NODES_V1: usize = 8_192;

/// Maximum number of weighted leaves retained by schema version one.
pub const MAX_STRUCTURED_LEAVES_V1: usize = 4_096;

/// Maximum root-to-leaf edge depth admitted by schema version one.
pub const MAX_STRUCTURED_DEPTH_V1: usize = 16;

/// Maximum number of piecewise response laws retained by schema version one.
pub const MAX_STRUCTURED_LAWS_V1: usize = 64;

/// Maximum number of quadratic pieces in one response law.
pub const MAX_STRUCTURED_PIECES_PER_LAW_V1: usize = 64;

/// Maximum byte length of one stable node, law, or mode key.
pub const MAX_STRUCTURED_KEY_BYTES_V1: usize = 128;

/// Maximum byte length of the external structured-model namespace.
pub const MAX_STRUCTURED_NAMESPACE_BYTES_V1: usize = 256;

/// Largest integer multiplicity sum exactly representable by binary64.
pub const MAX_EXACT_STRUCTURED_WEIGHT_V1: u64 = 1_u64 << 53;

/// Stable ordinal identity of one node in a structured population model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StructuredNodeId(
    /// Zero-based node ordinal in the retained model.
    pub u32,
);

impl StructuredNodeId {
    /// Zero-based node index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

/// Stable ordinal identity of one piecewise-quadratic law.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StructuredLawId(
    /// Zero-based law ordinal in the retained model.
    pub u16,
);

impl StructuredLawId {
    /// Zero-based law index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

/// Caller-supplied external identity retained with a structured receipt.
///
/// Admission checks grammar and a nonzero digest; it does not compute the
/// digest from model bytes or authenticate the external owner.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuredModelIdentity {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    semantic_digest: [u8; 32],
}

impl StructuredModelIdentity {
    /// Admit a bounded slash-separated external model identity.
    ///
    /// # Errors
    ///
    /// Refuses a noncanonical namespace or an all-zero semantic digest.
    pub fn try_new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        semantic_digest: [u8; 32],
    ) -> Result<Self, StructuredPropagationError> {
        let namespace = namespace.into();
        validate_namespace(&namespace)?;
        if semantic_digest == [0; 32] {
            return Err(StructuredPropagationError::ZeroSemanticDigest);
        }
        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            semantic_digest,
        })
    }

    /// External model namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Explicit external schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact semantic digest supplied by the external model owner.
    #[must_use]
    pub const fn semantic_digest(&self) -> [u8; 32] {
        self.semantic_digest
    }
}

/// One node in a parent-before-child manufacturing hierarchy.
#[derive(Debug, Clone, PartialEq)]
pub enum StructuredNodeSpec {
    /// A grouping node such as a process, lot, or assembly batch.
    Branch {
        /// Globally unique stable key.
        key: String,
        /// Parent branch, or none for the one root at ordinal zero.
        parent: Option<StructuredNodeId>,
    },
    /// A weighted terminal population observation.
    Leaf {
        /// Globally unique stable key.
        key: String,
        /// Parent branch.
        parent: StructuredNodeId,
        /// Exact finite-population multiplicity.
        relative_weight: NonZeroU64,
        /// Supplied raw clearance before law-domain clamping.
        raw_clearance: f64,
        /// Piecewise law used by this observation.
        law: StructuredLawId,
    },
}

impl StructuredNodeSpec {
    fn key(&self) -> &str {
        match self {
            Self::Branch { key, .. } | Self::Leaf { key, .. } => key,
        }
    }

    fn parent(&self) -> Option<StructuredNodeId> {
        match self {
            Self::Branch { parent, .. } => *parent,
            Self::Leaf { parent, .. } => Some(*parent),
        }
    }
}

/// One quadratic response piece evaluated as (a x + b) x + c.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticResponsePiece {
    /// Stable mode key for this interval.
    pub mode_key: String,
    /// Quadratic coefficient.
    pub a: f64,
    /// Linear coefficient.
    pub b: f64,
    /// Constant coefficient.
    pub c: f64,
}

/// A bounded, globally single-valued piecewise-quadratic clearance law.
#[derive(Debug, Clone, PartialEq)]
pub struct PiecewiseQuadraticLaw {
    /// Stable law key.
    pub key: String,
    /// Finite lower clamp bound.
    pub lower_bound: f64,
    /// Finite upper clamp bound.
    pub upper_bound: f64,
    /// Strictly increasing knots, including both clamp bounds.
    ///
    /// The number of knots must equal the number of pieces plus one.
    pub knots: Vec<f64>,
    /// Owner of every interior knot, in ascending knot order.
    ///
    /// The number of owners must equal the number of pieces minus one. Making
    /// boundary ownership explicit admits deadbands such as minus one through
    /// plus one inclusive without overlapping predicates or order-dependent
    /// callback behavior.
    pub interior_knot_owners: Vec<InteriorKnotOwner>,
    /// Quadratic responses in ascending interval order.
    pub pieces: Vec<QuadraticResponsePiece>,
}

/// Which adjacent response piece owns an exact interior-knot value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InteriorKnotOwner {
    /// The piece below the knot owns equality.
    LowerPiece,
    /// The piece above the knot owns equality.
    UpperPiece,
}

/// Exact structured finite-population model supplied for one evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredPopulationModel {
    /// Caller-supplied external identity retained alongside this model.
    pub identity: StructuredModelIdentity,
    /// Parent-before-child hierarchy. Ordinal zero must be the one root branch.
    pub nodes: Vec<StructuredNodeSpec>,
    /// Piecewise laws referenced by leaf ordinals.
    pub laws: Vec<PiecewiseQuadraticLaw>,
}

/// Stable category of a bounded retained resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredResource {
    /// Hierarchy nodes.
    Nodes,
    /// Terminal weighted observations.
    Leaves,
    /// Root-to-leaf hierarchy depth.
    HierarchyDepth,
    /// Piecewise response laws.
    Laws,
    /// Pieces in one response law.
    PiecesPerLaw,
    /// Exact total finite-population multiplicity.
    TotalWeight,
}

/// Stable key category used by a refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredKeyRole {
    /// Hierarchy node key.
    Node,
    /// Piecewise law key.
    Law,
    /// Law-local response-mode key.
    Mode,
}

/// Why one hierarchy relationship was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredTopologyIssue {
    /// Ordinal zero must be a branch.
    RootMustBeBranch,
    /// The root branch must not declare a parent.
    RootHasParent,
    /// Every non-root branch must declare a parent.
    MissingParent,
    /// A parent must have a lower ordinal than its child.
    ParentNotBeforeChild,
    /// A parent ordinal must identify a branch rather than a leaf.
    ParentIsLeaf,
    /// A branch without children does not describe a hierarchy.
    EmptyBranch,
}

/// Why one piecewise law layout was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredLawIssue {
    /// The lower bound is not strictly below the upper bound.
    BoundsNotIncreasing,
    /// Knots do not number exactly one more than pieces.
    KnotPieceCountMismatch,
    /// Interior-knot owners do not number exactly one fewer than pieces.
    InteriorOwnerCountMismatch,
    /// The first knot is not bit-identical to the lower bound.
    FirstKnotDoesNotMatchLowerBound,
    /// The last knot is not bit-identical to the upper bound.
    LastKnotDoesNotMatchUpperBound,
    /// Adjacent knots are not strictly increasing.
    KnotsNotIncreasing,
}

/// One quadratic coefficient named by a scalar refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredCoefficient {
    /// Quadratic coefficient a.
    Quadratic,
    /// Linear coefficient b.
    Linear,
    /// Constant coefficient c.
    Constant,
}

/// Exact location of an invalid structured-model scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredScalarLocation {
    /// Raw clearance on one leaf node.
    LeafRawClearance { node: StructuredNodeId },
    /// Lower clamp bound on one law.
    LawLowerBound { law: StructuredLawId },
    /// Upper clamp bound on one law.
    LawUpperBound { law: StructuredLawId },
    /// One knot in a law.
    LawKnot {
        /// Law ordinal.
        law: StructuredLawId,
        /// Knot ordinal.
        knot: usize,
    },
    /// One coefficient in one response piece.
    PieceCoefficient {
        /// Law ordinal.
        law: StructuredLawId,
        /// Piece ordinal.
        piece: usize,
        /// Coefficient role.
        coefficient: StructuredCoefficient,
    },
}

/// Stable binary64 refusal class for the structured lane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredNumericIssue {
    /// NaN or infinity was supplied or derived.
    NonFinite,
    /// Exact zero used the noncanonical negative-zero encoding.
    NonCanonicalNegativeZero,
    /// A nonzero product or rescaling rounded to zero.
    Underflow,
    /// A normalized mean escaped the convex hull scale because of arithmetic.
    OutsideNormalizedRange,
    /// A variance-like quantity became negative.
    Negative,
}

/// Fixed operation whose derived value could not be represented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredEvaluationStage {
    /// Product a x in the first Horner stage.
    QuadraticProduct,
    /// Sum a x + b in the first Horner stage.
    LinearSum,
    /// Product (a x + b) x in the second Horner stage.
    LinearProduct,
    /// Final sum (a x + b) x + c.
    OutputSum,
}

/// Scope of a moment-arithmetic refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredMomentScope {
    /// Direct whole-population audit.
    Population,
    /// One hierarchy node.
    Node(StructuredNodeId),
    /// One law-local response mode.
    Mode {
        /// Law ordinal.
        law: StructuredLawId,
        /// Piece ordinal.
        piece: usize,
    },
    /// Whole-population decomposition by selected response mode.
    ModeDecomposition,
}

/// Derived moment quantity that could not be represented safely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuredMomentQuantity {
    /// One nonzero observation normalized by the global output scale.
    NormalizedObservation,
    /// A weighted normalized mean.
    Mean,
    /// Immediate-child or within-mode variance.
    WithinVariance,
    /// Immediate-child or between-mode variance.
    BetweenVariance,
    /// Total population variance.
    TotalVariance,
    /// Standard deviation obtained from a positive normalized variance.
    StandardDeviation,
    /// A fixed-order compensated sum.
    Accumulation,
}

/// Refusal from identity admission, model admission, evaluation, or moments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StructuredPropagationError {
    /// The external namespace is not a bounded canonical slash-separated key.
    InvalidNamespace {
        /// Bounded prefix of the rejected namespace.
        namespace: String,
        /// Stable grammar explanation.
        reason: &'static str,
    },
    /// An all-zero digest cannot identify an external model.
    ZeroSemanticDigest,
    /// One bounded resource is empty or exceeds its schema limit.
    ResourceLimit {
        /// Resource class.
        resource: StructuredResource,
        /// Supplied or derived count.
        actual: u64,
        /// Inclusive versioned maximum.
        max: u64,
    },
    /// One stable key is empty, oversized, or noncanonical.
    InvalidKey {
        /// Key role.
        role: StructuredKeyRole,
        /// Owning law for a law-local mode key.
        owner_index: Option<usize>,
        /// Node, law, or piece ordinal.
        index: usize,
        /// Bounded prefix of the rejected key.
        key: String,
        /// Stable grammar explanation.
        reason: &'static str,
    },
    /// Two stable keys collide in their declared scope.
    DuplicateKey {
        /// Key role.
        role: StructuredKeyRole,
        /// Owning law for law-local mode keys.
        owner_index: Option<usize>,
        /// First ordinal.
        first_index: usize,
        /// Colliding ordinal.
        duplicate_index: usize,
        /// Exact canonical key.
        key: String,
    },
    /// One hierarchy edge or branch violates the parent-before-child grammar.
    InvalidTopology {
        /// Node being admitted.
        node: StructuredNodeId,
        /// Related parent ordinal when one was supplied.
        related: Option<StructuredNodeId>,
        /// Stable topology class.
        issue: StructuredTopologyIssue,
    },
    /// A leaf references a law ordinal outside the supplied law table.
    InvalidLawReference {
        /// Leaf node.
        node: StructuredNodeId,
        /// Supplied law ordinal.
        law: StructuredLawId,
        /// Number of available laws.
        available: usize,
    },
    /// One piecewise law has malformed bounds, knots, or piece cardinality.
    InvalidLawLayout {
        /// Law ordinal.
        law: StructuredLawId,
        /// Stable layout class.
        issue: StructuredLawIssue,
        /// Relevant knot ordinal when applicable.
        knot: Option<usize>,
    },
    /// One supplied scalar is not finite or uses noncanonical negative zero.
    InvalidScalar {
        /// Exact scalar location.
        location: StructuredScalarLocation,
        /// Stable numeric class.
        issue: StructuredNumericIssue,
    },
    /// A finite admitted leaf and law produced an unrepresentable result.
    InvalidEvaluation {
        /// Leaf node.
        node: StructuredNodeId,
        /// Selected law.
        law: StructuredLawId,
        /// Selected response piece.
        piece: usize,
        /// Fixed Horner operation.
        stage: StructuredEvaluationStage,
        /// Stable numeric class.
        issue: StructuredNumericIssue,
    },
    /// Moment accumulation or rescaling was not representable safely.
    InvalidMoment {
        /// Failing scope.
        scope: StructuredMomentScope,
        /// Failing quantity.
        quantity: StructuredMomentQuantity,
        /// Stable numeric class.
        issue: StructuredNumericIssue,
    },
}

/// Whether law-domain clamping changed one raw clearance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClampDisposition {
    /// Raw clearance was below the lower law bound.
    LowerClamped,
    /// Raw clearance already lay inside the closed law domain.
    Unchanged,
    /// Raw clearance was above the upper law bound.
    UpperClamped,
}

/// Retained evaluation of one weighted population leaf.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredLeafReceipt {
    node: StructuredNodeId,
    relative_weight: NonZeroU64,
    raw_clearance: f64,
    clamped_clearance: f64,
    clamp_disposition: ClampDisposition,
    law: StructuredLawId,
    selected_piece: usize,
    output: f64,
}

impl StructuredLeafReceipt {
    /// Leaf node identity.
    #[must_use]
    pub const fn node(&self) -> StructuredNodeId {
        self.node
    }

    /// Exact finite-population multiplicity.
    #[must_use]
    pub const fn relative_weight(&self) -> NonZeroU64 {
        self.relative_weight
    }

    /// Supplied raw clearance.
    #[must_use]
    pub const fn raw_clearance(&self) -> f64 {
        self.raw_clearance
    }

    /// Clearance after deterministic closed-domain clamping.
    #[must_use]
    pub const fn clamped_clearance(&self) -> f64 {
        self.clamped_clearance
    }

    /// Whether clamping changed the raw clearance.
    #[must_use]
    pub const fn clamp_disposition(&self) -> ClampDisposition {
        self.clamp_disposition
    }

    /// Selected law identity.
    #[must_use]
    pub const fn law(&self) -> StructuredLawId {
        self.law
    }

    /// Selected law-local piece ordinal.
    #[must_use]
    pub const fn selected_piece(&self) -> usize {
        self.selected_piece
    }

    /// Evaluated piecewise-quadratic response.
    #[must_use]
    pub const fn output(&self) -> f64 {
        self.output
    }
}

/// Population moments for one hierarchy node.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredNodeMomentReceipt {
    node: StructuredNodeId,
    relative_weight: u64,
    descendant_leaf_count: usize,
    mean: f64,
    within_child_variance: f64,
    between_child_variance: f64,
    total_variance: f64,
    standard_deviation: f64,
    decomposition_residual: f64,
}

impl StructuredNodeMomentReceipt {
    /// Node identity.
    #[must_use]
    pub const fn node(&self) -> StructuredNodeId {
        self.node
    }

    /// Exact total multiplicity of descendant leaves.
    #[must_use]
    pub const fn relative_weight(&self) -> u64 {
        self.relative_weight
    }

    /// Number of distinct descendant leaf records.
    #[must_use]
    pub const fn descendant_leaf_count(&self) -> usize {
        self.descendant_leaf_count
    }

    /// Weighted population mean.
    #[must_use]
    pub const fn mean(&self) -> f64 {
        self.mean
    }

    /// Immediate-child within variance.
    #[must_use]
    pub const fn within_child_variance(&self) -> f64 {
        self.within_child_variance
    }

    /// Immediate-child between variance.
    #[must_use]
    pub const fn between_child_variance(&self) -> f64 {
        self.between_child_variance
    }

    /// Total population variance below this node.
    #[must_use]
    pub const fn total_variance(&self) -> f64 {
        self.total_variance
    }

    /// Square root of the total population variance.
    #[must_use]
    pub const fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }

    /// Binary64 diagnostic total minus (within plus between).
    #[must_use]
    pub const fn decomposition_residual(&self) -> f64 {
        self.decomposition_residual
    }
}

/// Moments for one law-local response mode.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredModeMomentReceipt {
    law: StructuredLawId,
    piece: usize,
    relative_weight: u64,
    leaf_count: usize,
    mean: Option<f64>,
    variance: Option<f64>,
    standard_deviation: Option<f64>,
}

impl StructuredModeMomentReceipt {
    /// Law containing this mode.
    #[must_use]
    pub const fn law(&self) -> StructuredLawId {
        self.law
    }

    /// Law-local piece ordinal.
    #[must_use]
    pub const fn piece(&self) -> usize {
        self.piece
    }

    /// Exact total multiplicity selecting this mode.
    #[must_use]
    pub const fn relative_weight(&self) -> u64 {
        self.relative_weight
    }

    /// Number of distinct leaves selecting this mode.
    #[must_use]
    pub const fn leaf_count(&self) -> usize {
        self.leaf_count
    }

    /// Weighted mode mean, or none when the supplied population never selects it.
    #[must_use]
    pub const fn mean(&self) -> Option<f64> {
        self.mean
    }

    /// Population variance within this mode, or none when unobserved.
    #[must_use]
    pub const fn variance(&self) -> Option<f64> {
        self.variance
    }

    /// Population standard deviation within this mode, or none when unobserved.
    #[must_use]
    pub const fn standard_deviation(&self) -> Option<f64> {
        self.standard_deviation
    }
}

/// Privately constructed result of one admitted structured-population evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredPropagationReceipt {
    schema_version: u32,
    model: StructuredPopulationModel,
    leaves: Box<[StructuredLeafReceipt]>,
    nodes: Box<[StructuredNodeMomentReceipt]>,
    modes: Box<[StructuredModeMomentReceipt]>,
    total_weight: u64,
    mean: f64,
    variance: f64,
    standard_deviation: f64,
    direct_mean: f64,
    direct_variance: f64,
    hierarchy_mean_residual: f64,
    hierarchy_variance_residual: f64,
    within_mode_variance: f64,
    between_mode_variance: f64,
    mode_decomposed_variance: f64,
    mode_decomposition_residual: f64,
}

impl StructuredPropagationReceipt {
    /// Fixed structured-propagation receipt/evaluation schema.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact admitted model, including hierarchy, weights, laws, and identity.
    #[must_use]
    pub const fn model(&self) -> &StructuredPopulationModel {
        &self.model
    }

    /// Leaf evaluations in node-ordinal order.
    #[must_use]
    pub fn leaves(&self) -> &[StructuredLeafReceipt] {
        &self.leaves
    }

    /// Per-node moments in node-ordinal order; element zero is the root.
    #[must_use]
    pub fn nodes(&self) -> &[StructuredNodeMomentReceipt] {
        &self.nodes
    }

    /// Law-major, piece-minor response-mode moments.
    #[must_use]
    pub fn modes(&self) -> &[StructuredModeMomentReceipt] {
        &self.modes
    }

    /// Exact sum of finite-population leaf multiplicities.
    #[must_use]
    pub const fn total_weight(&self) -> u64 {
        self.total_weight
    }

    /// Hierarchically accumulated root mean.
    #[must_use]
    pub const fn mean(&self) -> f64 {
        self.mean
    }

    /// Hierarchically accumulated root population variance.
    #[must_use]
    pub const fn variance(&self) -> f64 {
        self.variance
    }

    /// Square root of the root population variance.
    #[must_use]
    pub const fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }

    /// Independent direct two-pass audit mean over all leaves.
    #[must_use]
    pub const fn direct_mean(&self) -> f64 {
        self.direct_mean
    }

    /// Independent direct two-pass audit variance over all leaves.
    #[must_use]
    pub const fn direct_variance(&self) -> f64 {
        self.direct_variance
    }

    /// Direct audit mean minus hierarchical root mean.
    #[must_use]
    pub const fn hierarchy_mean_residual(&self) -> f64 {
        self.hierarchy_mean_residual
    }

    /// Direct audit variance minus hierarchical root variance.
    #[must_use]
    pub const fn hierarchy_variance_residual(&self) -> f64 {
        self.hierarchy_variance_residual
    }

    /// Weighted variance retained inside selected response modes.
    #[must_use]
    pub const fn within_mode_variance(&self) -> f64 {
        self.within_mode_variance
    }

    /// Weighted variance between selected response-mode means.
    #[must_use]
    pub const fn between_mode_variance(&self) -> f64 {
        self.between_mode_variance
    }

    /// Binary64 sum of within-mode and between-mode variance.
    #[must_use]
    pub const fn mode_decomposed_variance(&self) -> f64 {
        self.mode_decomposed_variance
    }

    /// Root variance minus the mode-decomposed variance.
    #[must_use]
    pub const fn mode_decomposition_residual(&self) -> f64 {
        self.mode_decomposition_residual
    }
}

#[derive(Debug, Clone, Copy)]
struct InternalMoment {
    weight: u64,
    leaf_count: usize,
    mean: f64,
    within_sum_squares: f64,
    between_sum_squares: f64,
    total_sum_squares: f64,
}

impl Default for InternalMoment {
    fn default() -> Self {
        Self {
            weight: 0,
            leaf_count: 0,
            mean: 0.0,
            within_sum_squares: 0.0,
            between_sum_squares: 0.0,
            total_sum_squares: 0.0,
        }
    }
}

struct ValidatedTopology {
    children: Vec<Vec<usize>>,
    leaf_indices: Vec<usize>,
    total_weight: u64,
}

fn bounded_utf8_prefix(value: &str, max_bytes: usize) -> String {
    let mut end = value.len().min(max_bytes);
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    value[..end].to_string()
}

fn validate_namespace(namespace: &str) -> Result<(), StructuredPropagationError> {
    let reason = if namespace.is_empty() {
        Some("namespace must not be empty")
    } else if namespace.len() > MAX_STRUCTURED_NAMESPACE_BYTES_V1 {
        Some("namespace exceeds the versioned byte cap")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| segment.is_empty())
    {
        Some("namespace segments must not be empty")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| !is_canonical_key_bytes(segment))
    {
        Some(
            "namespace segments must start lowercase and use lowercase ASCII letters, digits, or single interior hyphens",
        )
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(StructuredPropagationError::InvalidNamespace {
            namespace: bounded_utf8_prefix(namespace, MAX_STRUCTURED_NAMESPACE_BYTES_V1),
            reason,
        })
    } else {
        Ok(())
    }
}

fn is_canonical_key_bytes(key: &[u8]) -> bool {
    key.first().is_some_and(u8::is_ascii_lowercase)
        && key.last() != Some(&b'-')
        && !key.windows(2).any(|pair| pair == b"--")
        && key
            .iter()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'-')
}

fn validate_key(
    role: StructuredKeyRole,
    owner_index: Option<usize>,
    index: usize,
    key: &str,
) -> Result<(), StructuredPropagationError> {
    let reason = if key.is_empty() {
        Some("key must not be empty")
    } else if key.len() > MAX_STRUCTURED_KEY_BYTES_V1 {
        Some("key exceeds the versioned byte cap")
    } else if !key.is_ascii() || !is_canonical_key_bytes(key.as_bytes()) {
        Some(
            "key must start lowercase and use lowercase ASCII letters, digits, or single interior hyphens",
        )
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(StructuredPropagationError::InvalidKey {
            role,
            owner_index,
            index,
            key: bounded_utf8_prefix(key, MAX_STRUCTURED_KEY_BYTES_V1),
            reason,
        })
    } else {
        Ok(())
    }
}

fn scalar_issue(value: f64) -> Option<StructuredNumericIssue> {
    if !value.is_finite() {
        Some(StructuredNumericIssue::NonFinite)
    } else if value == 0.0 && value.is_sign_negative() {
        Some(StructuredNumericIssue::NonCanonicalNegativeZero)
    } else {
        None
    }
}

fn law_id(index: usize) -> StructuredLawId {
    StructuredLawId(index as u16)
}

fn node_id(index: usize) -> StructuredNodeId {
    StructuredNodeId(index as u32)
}

fn validate_laws(model: &StructuredPopulationModel) -> Result<(), StructuredPropagationError> {
    if model.laws.is_empty() || model.laws.len() > MAX_STRUCTURED_LAWS_V1 {
        return Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::Laws,
            actual: model.laws.len() as u64,
            max: MAX_STRUCTURED_LAWS_V1 as u64,
        });
    }

    let mut law_keys = BTreeMap::new();
    for (law_index, law) in model.laws.iter().enumerate() {
        validate_key(StructuredKeyRole::Law, None, law_index, &law.key)?;
        if let Some(&first_index) = law_keys.get(&law.key) {
            return Err(StructuredPropagationError::DuplicateKey {
                role: StructuredKeyRole::Law,
                owner_index: None,
                first_index,
                duplicate_index: law_index,
                key: law.key.clone(),
            });
        }
        law_keys.insert(law.key.clone(), law_index);

        let id = law_id(law_index);
        for (location, value) in [
            (
                StructuredScalarLocation::LawLowerBound { law: id },
                law.lower_bound,
            ),
            (
                StructuredScalarLocation::LawUpperBound { law: id },
                law.upper_bound,
            ),
        ] {
            if let Some(issue) = scalar_issue(value) {
                return Err(StructuredPropagationError::InvalidScalar { location, issue });
            }
        }
        if law.lower_bound >= law.upper_bound {
            return Err(StructuredPropagationError::InvalidLawLayout {
                law: id,
                issue: StructuredLawIssue::BoundsNotIncreasing,
                knot: None,
            });
        }
        if law.pieces.is_empty() || law.pieces.len() > MAX_STRUCTURED_PIECES_PER_LAW_V1 {
            return Err(StructuredPropagationError::ResourceLimit {
                resource: StructuredResource::PiecesPerLaw,
                actual: law.pieces.len() as u64,
                max: MAX_STRUCTURED_PIECES_PER_LAW_V1 as u64,
            });
        }
        if law.knots.len() != law.pieces.len() + 1 {
            return Err(StructuredPropagationError::InvalidLawLayout {
                law: id,
                issue: StructuredLawIssue::KnotPieceCountMismatch,
                knot: None,
            });
        }
        if law.interior_knot_owners.len() + 1 != law.pieces.len() {
            return Err(StructuredPropagationError::InvalidLawLayout {
                law: id,
                issue: StructuredLawIssue::InteriorOwnerCountMismatch,
                knot: None,
            });
        }
        for (knot_index, &knot) in law.knots.iter().enumerate() {
            if let Some(issue) = scalar_issue(knot) {
                return Err(StructuredPropagationError::InvalidScalar {
                    location: StructuredScalarLocation::LawKnot {
                        law: id,
                        knot: knot_index,
                    },
                    issue,
                });
            }
        }
        if law.knots[0].to_bits() != law.lower_bound.to_bits() {
            return Err(StructuredPropagationError::InvalidLawLayout {
                law: id,
                issue: StructuredLawIssue::FirstKnotDoesNotMatchLowerBound,
                knot: Some(0),
            });
        }
        let last_knot = law.knots.len() - 1;
        if law.knots[last_knot].to_bits() != law.upper_bound.to_bits() {
            return Err(StructuredPropagationError::InvalidLawLayout {
                law: id,
                issue: StructuredLawIssue::LastKnotDoesNotMatchUpperBound,
                knot: Some(last_knot),
            });
        }
        for knot_index in 1..law.knots.len() {
            if law.knots[knot_index - 1] >= law.knots[knot_index] {
                return Err(StructuredPropagationError::InvalidLawLayout {
                    law: id,
                    issue: StructuredLawIssue::KnotsNotIncreasing,
                    knot: Some(knot_index),
                });
            }
        }

        let mut mode_keys = BTreeMap::new();
        for (piece_index, piece) in law.pieces.iter().enumerate() {
            validate_key(
                StructuredKeyRole::Mode,
                Some(law_index),
                piece_index,
                &piece.mode_key,
            )?;
            if let Some(&first_index) = mode_keys.get(&piece.mode_key) {
                return Err(StructuredPropagationError::DuplicateKey {
                    role: StructuredKeyRole::Mode,
                    owner_index: Some(law_index),
                    first_index,
                    duplicate_index: piece_index,
                    key: piece.mode_key.clone(),
                });
            }
            mode_keys.insert(piece.mode_key.clone(), piece_index);
            for (coefficient, value) in [
                (StructuredCoefficient::Quadratic, piece.a),
                (StructuredCoefficient::Linear, piece.b),
                (StructuredCoefficient::Constant, piece.c),
            ] {
                if let Some(issue) = scalar_issue(value) {
                    return Err(StructuredPropagationError::InvalidScalar {
                        location: StructuredScalarLocation::PieceCoefficient {
                            law: id,
                            piece: piece_index,
                            coefficient,
                        },
                        issue,
                    });
                }
            }
        }
    }
    Ok(())
}

fn validate_topology(
    model: &StructuredPopulationModel,
) -> Result<ValidatedTopology, StructuredPropagationError> {
    if model.nodes.is_empty() || model.nodes.len() > MAX_STRUCTURED_NODES_V1 {
        return Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::Nodes,
            actual: model.nodes.len() as u64,
            max: MAX_STRUCTURED_NODES_V1 as u64,
        });
    }
    match &model.nodes[0] {
        StructuredNodeSpec::Branch { parent: None, .. } => {}
        StructuredNodeSpec::Branch {
            parent: Some(parent),
            ..
        } => {
            return Err(StructuredPropagationError::InvalidTopology {
                node: node_id(0),
                related: Some(*parent),
                issue: StructuredTopologyIssue::RootHasParent,
            });
        }
        StructuredNodeSpec::Leaf { .. } => {
            return Err(StructuredPropagationError::InvalidTopology {
                node: node_id(0),
                related: None,
                issue: StructuredTopologyIssue::RootMustBeBranch,
            });
        }
    }

    let mut node_keys = BTreeMap::new();
    let mut children = vec![Vec::new(); model.nodes.len()];
    let mut depths = vec![0_usize; model.nodes.len()];
    let mut leaf_indices = Vec::new();
    let mut total_weight = 0_u64;
    for (index, node) in model.nodes.iter().enumerate() {
        validate_key(StructuredKeyRole::Node, None, index, node.key())?;
        if let Some(&first_index) = node_keys.get(node.key()) {
            return Err(StructuredPropagationError::DuplicateKey {
                role: StructuredKeyRole::Node,
                owner_index: None,
                first_index,
                duplicate_index: index,
                key: node.key().to_string(),
            });
        }
        node_keys.insert(node.key().to_string(), index);

        if index > 0 {
            let parent = node
                .parent()
                .ok_or(StructuredPropagationError::InvalidTopology {
                    node: node_id(index),
                    related: None,
                    issue: StructuredTopologyIssue::MissingParent,
                })?;
            let parent_index = parent.index();
            if parent_index >= index {
                return Err(StructuredPropagationError::InvalidTopology {
                    node: node_id(index),
                    related: Some(parent),
                    issue: StructuredTopologyIssue::ParentNotBeforeChild,
                });
            }
            if !matches!(model.nodes[parent_index], StructuredNodeSpec::Branch { .. }) {
                return Err(StructuredPropagationError::InvalidTopology {
                    node: node_id(index),
                    related: Some(parent),
                    issue: StructuredTopologyIssue::ParentIsLeaf,
                });
            }
            children[parent_index].push(index);
            depths[index] = depths[parent_index] + 1;
            if depths[index] > MAX_STRUCTURED_DEPTH_V1 {
                return Err(StructuredPropagationError::ResourceLimit {
                    resource: StructuredResource::HierarchyDepth,
                    actual: depths[index] as u64,
                    max: MAX_STRUCTURED_DEPTH_V1 as u64,
                });
            }
        }

        if let StructuredNodeSpec::Leaf {
            relative_weight,
            raw_clearance,
            law,
            ..
        } = node
        {
            if law.index() >= model.laws.len() {
                return Err(StructuredPropagationError::InvalidLawReference {
                    node: node_id(index),
                    law: *law,
                    available: model.laws.len(),
                });
            }
            if let Some(issue) = scalar_issue(*raw_clearance) {
                return Err(StructuredPropagationError::InvalidScalar {
                    location: StructuredScalarLocation::LeafRawClearance {
                        node: node_id(index),
                    },
                    issue,
                });
            }
            total_weight = total_weight.checked_add(relative_weight.get()).ok_or(
                StructuredPropagationError::ResourceLimit {
                    resource: StructuredResource::TotalWeight,
                    actual: u64::MAX,
                    max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
                },
            )?;
            if total_weight > MAX_EXACT_STRUCTURED_WEIGHT_V1 {
                return Err(StructuredPropagationError::ResourceLimit {
                    resource: StructuredResource::TotalWeight,
                    actual: total_weight,
                    max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
                });
            }
            leaf_indices.push(index);
        }
    }

    if leaf_indices.is_empty() || leaf_indices.len() > MAX_STRUCTURED_LEAVES_V1 {
        return Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::Leaves,
            actual: leaf_indices.len() as u64,
            max: MAX_STRUCTURED_LEAVES_V1 as u64,
        });
    }
    for (index, node) in model.nodes.iter().enumerate() {
        if matches!(node, StructuredNodeSpec::Branch { .. }) && children[index].is_empty() {
            return Err(StructuredPropagationError::InvalidTopology {
                node: node_id(index),
                related: None,
                issue: StructuredTopologyIssue::EmptyBranch,
            });
        }
    }

    Ok(ValidatedTopology {
        children,
        leaf_indices,
        total_weight,
    })
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn select_piece(law: &PiecewiseQuadraticLaw, value: f64) -> usize {
    for piece in 0..law.pieces.len() {
        let upper = law.knots[piece + 1];
        if value < upper || piece + 1 == law.pieces.len() {
            return piece;
        }
        if value == upper {
            return match law.interior_knot_owners[piece] {
                InteriorKnotOwner::LowerPiece => piece,
                InteriorKnotOwner::UpperPiece => piece + 1,
            };
        }
    }
    law.pieces.len().saturating_sub(1)
}

fn evaluate_piece(
    node: StructuredNodeId,
    law: StructuredLawId,
    piece_index: usize,
    piece: &QuadraticResponsePiece,
    value: f64,
) -> Result<f64, StructuredPropagationError> {
    let quadratic_product = piece.a * value;
    if !quadratic_product.is_finite() {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::QuadraticProduct,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if piece.a != 0.0 && value != 0.0 && quadratic_product == 0.0 {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::QuadraticProduct,
            issue: StructuredNumericIssue::Underflow,
        });
    }

    let linear_sum = quadratic_product + piece.b;
    if !linear_sum.is_finite() {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::LinearSum,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    let linear_product = linear_sum * value;
    if !linear_product.is_finite() {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::LinearProduct,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if linear_sum != 0.0 && value != 0.0 && linear_product == 0.0 {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::LinearProduct,
            issue: StructuredNumericIssue::Underflow,
        });
    }

    let output = linear_product + piece.c;
    if !output.is_finite() {
        return Err(StructuredPropagationError::InvalidEvaluation {
            node,
            law,
            piece: piece_index,
            stage: StructuredEvaluationStage::OutputSum,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    Ok(canonical_zero(output))
}

fn rescale_mean(
    scale: f64,
    normalized_mean: f64,
    scope: StructuredMomentScope,
) -> Result<f64, StructuredPropagationError> {
    if !normalized_mean.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if normalized_mean.abs() > 1.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::OutsideNormalizedRange,
        });
    }
    if scale == 0.0 || normalized_mean == 0.0 {
        return Ok(0.0);
    }
    let mean = scale * normalized_mean;
    if !mean.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if mean == 0.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::Underflow,
        });
    }
    Ok(mean)
}

fn rescale_variance(
    scale: f64,
    normalized_variance: f64,
    scope: StructuredMomentScope,
    quantity: StructuredMomentQuantity,
) -> Result<(f64, f64), StructuredPropagationError> {
    if !normalized_variance.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if normalized_variance < 0.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity,
            issue: StructuredNumericIssue::Negative,
        });
    }
    if scale == 0.0 || normalized_variance == 0.0 {
        return Ok((0.0, 0.0));
    }
    let normalized_standard_deviation = normalized_variance.sqrt();
    let standard_deviation = scale * normalized_standard_deviation;
    if !standard_deviation.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::StandardDeviation,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if standard_deviation == 0.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::StandardDeviation,
            issue: StructuredNumericIssue::Underflow,
        });
    }
    let variance = standard_deviation * standard_deviation;
    if !variance.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if variance == 0.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity,
            issue: StructuredNumericIssue::Underflow,
        });
    }
    Ok((standard_deviation, variance))
}

fn checked_accumulation(
    value: f64,
    scope: StructuredMomentScope,
) -> Result<f64, StructuredPropagationError> {
    if value.is_finite() {
        Ok(canonical_zero(value))
    } else {
        Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Accumulation,
            issue: StructuredNumericIssue::NonFinite,
        })
    }
}

fn direct_normalized_moment(
    observations: &[(u64, f64)],
    scope: StructuredMomentScope,
) -> Result<InternalMoment, StructuredPropagationError> {
    let weight = observations.iter().try_fold(0_u64, |total, (weight, _)| {
        total
            .checked_add(*weight)
            .ok_or(StructuredPropagationError::ResourceLimit {
                resource: StructuredResource::TotalWeight,
                actual: u64::MAX,
                max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
            })
    })?;
    if weight == 0 || weight > MAX_EXACT_STRUCTURED_WEIGHT_V1 {
        return Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::TotalWeight,
            actual: weight,
            max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
        });
    }
    let weight_as_f64 = weight as f64;
    let weighted_sum = checked_accumulation(
        compensated_sum(
            observations
                .iter()
                .map(|(weight, value)| *weight as f64 * *value),
        ),
        scope,
    )?;
    let mean = canonical_zero(weighted_sum / weight_as_f64);
    if !mean.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    if mean.abs() > 1.0 {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::Mean,
            issue: StructuredNumericIssue::OutsideNormalizedRange,
        });
    }

    let mut terms = Vec::with_capacity(observations.len());
    for (observation_weight, value) in observations {
        let difference = *value - mean;
        let square = difference * difference;
        if difference != 0.0 && square == 0.0 {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::TotalVariance,
                issue: StructuredNumericIssue::Underflow,
            });
        }
        let term = *observation_weight as f64 * square;
        if !term.is_finite() {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::TotalVariance,
                issue: StructuredNumericIssue::NonFinite,
            });
        }
        terms.push(term);
    }
    let total_sum_squares = checked_accumulation(compensated_sum(terms), scope)?;
    Ok(InternalMoment {
        weight,
        leaf_count: observations.len(),
        mean,
        within_sum_squares: total_sum_squares,
        between_sum_squares: 0.0,
        total_sum_squares,
    })
}

fn build_hierarchy_moments(
    model: &StructuredPopulationModel,
    topology: &ValidatedTopology,
    normalized_by_node: &[Option<f64>],
) -> Result<Vec<InternalMoment>, StructuredPropagationError> {
    let mut moments = vec![InternalMoment::default(); model.nodes.len()];
    for (index, node) in model.nodes.iter().enumerate() {
        let StructuredNodeSpec::Leaf {
            relative_weight, ..
        } = node
        else {
            continue;
        };
        let mean = normalized_by_node[index].ok_or(StructuredPropagationError::InvalidMoment {
            scope: StructuredMomentScope::Node(node_id(index)),
            quantity: StructuredMomentQuantity::NormalizedObservation,
            issue: StructuredNumericIssue::NonFinite,
        })?;
        moments[index] = InternalMoment {
            weight: relative_weight.get(),
            leaf_count: 1,
            mean,
            within_sum_squares: 0.0,
            between_sum_squares: 0.0,
            total_sum_squares: 0.0,
        };
    }

    for index in (0..model.nodes.len()).rev() {
        if !matches!(model.nodes[index], StructuredNodeSpec::Branch { .. }) {
            continue;
        }
        let scope = StructuredMomentScope::Node(node_id(index));
        let child_indices = &topology.children[index];
        let weight = child_indices.iter().try_fold(0_u64, |total, child| {
            total.checked_add(moments[*child].weight).ok_or(
                StructuredPropagationError::ResourceLimit {
                    resource: StructuredResource::TotalWeight,
                    actual: u64::MAX,
                    max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
                },
            )
        })?;
        if weight == 0 || weight > MAX_EXACT_STRUCTURED_WEIGHT_V1 {
            return Err(StructuredPropagationError::ResourceLimit {
                resource: StructuredResource::TotalWeight,
                actual: weight,
                max: MAX_EXACT_STRUCTURED_WEIGHT_V1,
            });
        }
        let leaf_count = child_indices
            .iter()
            .map(|child| moments[*child].leaf_count)
            .sum();
        let weighted_sum = checked_accumulation(
            compensated_sum(
                child_indices
                    .iter()
                    .map(|child| moments[*child].weight as f64 * moments[*child].mean),
            ),
            scope,
        )?;
        let mean = canonical_zero(weighted_sum / weight as f64);
        if !mean.is_finite() {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::Mean,
                issue: StructuredNumericIssue::NonFinite,
            });
        }
        if mean.abs() > 1.0 {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::Mean,
                issue: StructuredNumericIssue::OutsideNormalizedRange,
            });
        }
        let within_sum_squares = checked_accumulation(
            compensated_sum(
                child_indices
                    .iter()
                    .map(|child| moments[*child].total_sum_squares),
            ),
            scope,
        )?;
        let mut between_terms = Vec::with_capacity(child_indices.len());
        for child in child_indices {
            let difference = moments[*child].mean - mean;
            let square = difference * difference;
            if difference != 0.0 && square == 0.0 {
                return Err(StructuredPropagationError::InvalidMoment {
                    scope,
                    quantity: StructuredMomentQuantity::BetweenVariance,
                    issue: StructuredNumericIssue::Underflow,
                });
            }
            let term = moments[*child].weight as f64 * square;
            if !term.is_finite() {
                return Err(StructuredPropagationError::InvalidMoment {
                    scope,
                    quantity: StructuredMomentQuantity::BetweenVariance,
                    issue: StructuredNumericIssue::NonFinite,
                });
            }
            between_terms.push(term);
        }
        let between_sum_squares = checked_accumulation(compensated_sum(between_terms), scope)?;
        let total_sum_squares =
            checked_accumulation(within_sum_squares + between_sum_squares, scope)?;
        moments[index] = InternalMoment {
            weight,
            leaf_count,
            mean,
            within_sum_squares,
            between_sum_squares,
            total_sum_squares,
        };
    }
    Ok(moments)
}

fn publish_node_moments(
    moments: &[InternalMoment],
    scale: f64,
) -> Result<Vec<StructuredNodeMomentReceipt>, StructuredPropagationError> {
    let mut receipts = Vec::with_capacity(moments.len());
    for (index, moment) in moments.iter().enumerate() {
        let scope = StructuredMomentScope::Node(node_id(index));
        let denominator = moment.weight as f64;
        let mean = rescale_mean(scale, moment.mean, scope)?;
        let (_, within_child_variance) = rescale_variance(
            scale,
            moment.within_sum_squares / denominator,
            scope,
            StructuredMomentQuantity::WithinVariance,
        )?;
        let (_, between_child_variance) = rescale_variance(
            scale,
            moment.between_sum_squares / denominator,
            scope,
            StructuredMomentQuantity::BetweenVariance,
        )?;
        let (standard_deviation, total_variance) = rescale_variance(
            scale,
            moment.total_sum_squares / denominator,
            scope,
            StructuredMomentQuantity::TotalVariance,
        )?;
        let decomposed = within_child_variance + between_child_variance;
        if !decomposed.is_finite() {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::TotalVariance,
                issue: StructuredNumericIssue::NonFinite,
            });
        }
        receipts.push(StructuredNodeMomentReceipt {
            node: node_id(index),
            relative_weight: moment.weight,
            descendant_leaf_count: moment.leaf_count,
            mean,
            within_child_variance,
            between_child_variance,
            total_variance,
            standard_deviation,
            decomposition_residual: canonical_zero(total_variance - decomposed),
        });
    }
    Ok(receipts)
}

struct PublishedModes {
    receipts: Vec<StructuredModeMomentReceipt>,
    within_variance: f64,
    between_variance: f64,
    decomposed_variance: f64,
}

fn publish_mode_moments(
    model: &StructuredPopulationModel,
    leaves: &[StructuredLeafReceipt],
    normalized_outputs: &[f64],
    root_mean: f64,
    total_weight: u64,
    scale: f64,
) -> Result<PublishedModes, StructuredPropagationError> {
    let mut offsets = Vec::with_capacity(model.laws.len());
    let mut total_modes = 0_usize;
    for law in &model.laws {
        offsets.push(total_modes);
        total_modes += law.pieces.len();
    }
    let mut buckets = vec![Vec::<(u64, f64)>::new(); total_modes];
    for (leaf_index, leaf) in leaves.iter().enumerate() {
        let bucket = offsets[leaf.law.index()] + leaf.selected_piece;
        buckets[bucket].push((leaf.relative_weight.get(), normalized_outputs[leaf_index]));
    }

    let mut receipts = Vec::with_capacity(total_modes);
    let mut summaries = Vec::with_capacity(total_modes);
    let mut bucket_index = 0_usize;
    for (law_index, law) in model.laws.iter().enumerate() {
        let id = law_id(law_index);
        for piece_index in 0..law.pieces.len() {
            let scope = StructuredMomentScope::Mode {
                law: id,
                piece: piece_index,
            };
            let observations = &buckets[bucket_index];
            if observations.is_empty() {
                receipts.push(StructuredModeMomentReceipt {
                    law: id,
                    piece: piece_index,
                    relative_weight: 0,
                    leaf_count: 0,
                    mean: None,
                    variance: None,
                    standard_deviation: None,
                });
                summaries.push(None);
            } else {
                let summary = direct_normalized_moment(observations, scope)?;
                let mean = rescale_mean(scale, summary.mean, scope)?;
                let (standard_deviation, variance) = rescale_variance(
                    scale,
                    summary.total_sum_squares / summary.weight as f64,
                    scope,
                    StructuredMomentQuantity::TotalVariance,
                )?;
                receipts.push(StructuredModeMomentReceipt {
                    law: id,
                    piece: piece_index,
                    relative_weight: summary.weight,
                    leaf_count: summary.leaf_count,
                    mean: Some(mean),
                    variance: Some(variance),
                    standard_deviation: Some(standard_deviation),
                });
                summaries.push(Some(summary));
            }
            bucket_index += 1;
        }
    }

    let scope = StructuredMomentScope::ModeDecomposition;
    let within_sum_squares = checked_accumulation(
        compensated_sum(
            summaries
                .iter()
                .flatten()
                .map(|summary| summary.total_sum_squares),
        ),
        scope,
    )?;
    let mut between_terms = Vec::new();
    for summary in summaries.iter().flatten() {
        let difference = summary.mean - root_mean;
        let square = difference * difference;
        if difference != 0.0 && square == 0.0 {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::BetweenVariance,
                issue: StructuredNumericIssue::Underflow,
            });
        }
        let term = summary.weight as f64 * square;
        if !term.is_finite() {
            return Err(StructuredPropagationError::InvalidMoment {
                scope,
                quantity: StructuredMomentQuantity::BetweenVariance,
                issue: StructuredNumericIssue::NonFinite,
            });
        }
        between_terms.push(term);
    }
    let between_sum_squares = checked_accumulation(compensated_sum(between_terms), scope)?;
    let (_, within_variance) = rescale_variance(
        scale,
        within_sum_squares / total_weight as f64,
        scope,
        StructuredMomentQuantity::WithinVariance,
    )?;
    let (_, between_variance) = rescale_variance(
        scale,
        between_sum_squares / total_weight as f64,
        scope,
        StructuredMomentQuantity::BetweenVariance,
    )?;
    let decomposed_variance = within_variance + between_variance;
    if !decomposed_variance.is_finite() {
        return Err(StructuredPropagationError::InvalidMoment {
            scope,
            quantity: StructuredMomentQuantity::TotalVariance,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    Ok(PublishedModes {
        receipts,
        within_variance,
        between_variance,
        decomposed_variance,
    })
}

/// Admit and evaluate one bounded hierarchical finite population.
///
/// Every leaf is clamped to its selected law domain, selects exactly one
/// explicitly owned knot interval, and evaluates a fixed two-multiply,
/// two-add Horner response. Integer multiplicities are accumulated exactly up
/// to 2^53 before conversion to binary64. Hierarchy and mode moments use fixed
/// input order, compensated sums, a global output scale, and the population
/// variance denominator W.
///
/// # Errors
///
/// Refuses unstable identity or keys, an empty or oversized model, malformed
/// topology or law tables, invalid scalars, unrepresentable response
/// arithmetic, or unsafe moment accumulation/rescaling. No partial receipt is
/// returned.
pub fn propagate_structured_population(
    model: &StructuredPopulationModel,
) -> Result<StructuredPropagationReceipt, StructuredPropagationError> {
    validate_laws(model)?;
    let topology = validate_topology(model)?;

    let mut leaves = Vec::with_capacity(topology.leaf_indices.len());
    let mut scale = 0.0_f64;
    for (index, node) in model.nodes.iter().enumerate() {
        let StructuredNodeSpec::Leaf {
            relative_weight,
            raw_clearance,
            law,
            ..
        } = node
        else {
            continue;
        };
        let response_law = &model.laws[law.index()];
        let (clamped_clearance, clamp_disposition) = if *raw_clearance < response_law.lower_bound {
            (response_law.lower_bound, ClampDisposition::LowerClamped)
        } else if *raw_clearance > response_law.upper_bound {
            (response_law.upper_bound, ClampDisposition::UpperClamped)
        } else {
            (*raw_clearance, ClampDisposition::Unchanged)
        };
        let selected_piece = select_piece(response_law, clamped_clearance);
        let output = evaluate_piece(
            node_id(index),
            *law,
            selected_piece,
            &response_law.pieces[selected_piece],
            clamped_clearance,
        )?;
        scale = scale.max(output.abs());
        leaves.push(StructuredLeafReceipt {
            node: node_id(index),
            relative_weight: *relative_weight,
            raw_clearance: *raw_clearance,
            clamped_clearance,
            clamp_disposition,
            law: *law,
            selected_piece,
            output,
        });
    }

    let mut normalized_outputs = Vec::with_capacity(leaves.len());
    let mut normalized_by_node = vec![None; model.nodes.len()];
    for leaf in &leaves {
        let normalized = if scale == 0.0 {
            0.0
        } else {
            leaf.output / scale
        };
        if !normalized.is_finite() {
            return Err(StructuredPropagationError::InvalidMoment {
                scope: StructuredMomentScope::Population,
                quantity: StructuredMomentQuantity::NormalizedObservation,
                issue: StructuredNumericIssue::NonFinite,
            });
        }
        if leaf.output != 0.0 && normalized == 0.0 {
            return Err(StructuredPropagationError::InvalidMoment {
                scope: StructuredMomentScope::Population,
                quantity: StructuredMomentQuantity::NormalizedObservation,
                issue: StructuredNumericIssue::Underflow,
            });
        }
        let normalized = canonical_zero(normalized);
        normalized_by_node[leaf.node.index()] = Some(normalized);
        normalized_outputs.push(normalized);
    }

    let hierarchy = build_hierarchy_moments(model, &topology, &normalized_by_node)?;
    let node_receipts = publish_node_moments(&hierarchy, scale)?;
    let root = hierarchy[0];
    debug_assert_eq!(root.weight, topology.total_weight);

    let direct_observations = leaves
        .iter()
        .zip(&normalized_outputs)
        .map(|(leaf, normalized)| (leaf.relative_weight.get(), *normalized))
        .collect::<Vec<_>>();
    let direct = direct_normalized_moment(&direct_observations, StructuredMomentScope::Population)?;
    let direct_mean = rescale_mean(scale, direct.mean, StructuredMomentScope::Population)?;
    let (_, direct_variance) = rescale_variance(
        scale,
        direct.total_sum_squares / direct.weight as f64,
        StructuredMomentScope::Population,
        StructuredMomentQuantity::TotalVariance,
    )?;

    let published_modes = publish_mode_moments(
        model,
        &leaves,
        &normalized_outputs,
        root.mean,
        topology.total_weight,
        scale,
    )?;
    let root_receipt = &node_receipts[0];
    let hierarchy_mean_residual = canonical_zero(direct_mean - root_receipt.mean);
    let hierarchy_variance_residual = canonical_zero(direct_variance - root_receipt.total_variance);
    let mode_decomposition_residual =
        canonical_zero(root_receipt.total_variance - published_modes.decomposed_variance);
    if !hierarchy_mean_residual.is_finite()
        || !hierarchy_variance_residual.is_finite()
        || !mode_decomposition_residual.is_finite()
    {
        return Err(StructuredPropagationError::InvalidMoment {
            scope: StructuredMomentScope::Population,
            quantity: StructuredMomentQuantity::Accumulation,
            issue: StructuredNumericIssue::NonFinite,
        });
    }
    let mean = root_receipt.mean;
    let variance = root_receipt.total_variance;
    let standard_deviation = root_receipt.standard_deviation;

    Ok(StructuredPropagationReceipt {
        schema_version: STRUCTURED_PROPAGATION_SCHEMA_V1,
        model: model.clone(),
        leaves: leaves.into_boxed_slice(),
        nodes: node_receipts.into_boxed_slice(),
        modes: published_modes.receipts.into_boxed_slice(),
        total_weight: topology.total_weight,
        mean,
        variance,
        standard_deviation,
        direct_mean,
        direct_variance,
        hierarchy_mean_residual,
        hierarchy_variance_residual,
        within_mode_variance: published_modes.within_variance,
        between_mode_variance: published_modes.between_variance,
        mode_decomposed_variance: published_modes.decomposed_variance,
        mode_decomposition_residual,
    })
}

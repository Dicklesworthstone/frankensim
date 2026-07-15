//! The typed problem graph: expression nodes over manifold variables,
//! hash-consed (common subexpressions are IDENTICAL node ids), with
//! per-step validation — shape rules, dimension rules (fs-qty `Dims`),
//! and differentiability-CLASS propagation, so "this objective is
//! non-smooth through that min()" is knowable at BUILD time and
//! optimizer routing can refuse with the offending node named.

use crate::admission::{self, AdmissionCaps, AdmissionReport, ProblemAdmission};
use crate::serial::{LegacyProblemHash, ProblemSemanticId};
use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::Dims;
use std::collections::{BTreeMap, btree_map::Entry};

/// A design variable handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VarId(pub u32);

/// An expression node handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct NodeId(pub u32);

/// The manifold a variable lives on — with the retraction metadata the
/// gradient stack consumes ("optimize an orientation" never becomes
/// "optimize 9 numbers and renormalize when it explodes").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Manifold {
    /// Flat Euclidean space.
    Rn {
        /// Dimension.
        dim: u32,
    },
    /// The unit sphere in `ambient` dimensions (points stored ambient).
    Sphere {
        /// Ambient dimension (≥ 2).
        ambient: u32,
    },
    /// Rotations, stored as unit quaternions (w, x, y, z).
    So3,
    /// Orthonormal `p`-frames in `n` dimensions (column-major storage).
    Stiefel {
        /// Ambient dimension.
        n: u32,
        /// Frame size (≤ n).
        p: u32,
    },
}

impl Manifold {
    /// Storage length of one point, computed with CHECKED arithmetic.
    /// `None` means the descriptor's formulas leave the `u32` domain
    /// (e.g. `Stiefel` with `n * p` overflow) — such a descriptor can
    /// never enter a [`Problem`] because [`ProblemBuilder::var`]
    /// validates through [`Manifold::validate`] first.
    #[must_use]
    pub fn point_dim(&self) -> Option<u32> {
        match *self {
            Manifold::Rn { dim } => Some(dim),
            Manifold::Sphere { ambient } => Some(ambient),
            Manifold::So3 => Some(4),
            Manifold::Stiefel { n, p } => n.checked_mul(p),
        }
    }

    /// Tangent-space dimension (what a Riemannian gradient has),
    /// computed with CHECKED arithmetic; `None` when the formula is
    /// not representable or meaningful for the raw descriptor.
    #[must_use]
    pub fn tangent_dim(&self) -> Option<u32> {
        match *self {
            Manifold::Rn { dim } => Some(dim),
            Manifold::Sphere { ambient } => ambient.checked_sub(1),
            Manifold::So3 => Some(3),
            Manifold::Stiefel { n, p } => {
                let np = u64::from(n).checked_mul(u64::from(p))?;
                let correction = u64::from(p).checked_mul(u64::from(p) + 1)? / 2;
                u32::try_from(np.checked_sub(correction)?).ok()
            }
        }
    }

    /// Validate this descriptor against the versioned admission policy:
    /// `Rn` needs `dim >= 1` (a zero-storage variable cannot bind a
    /// point), `Sphere` needs `ambient >= 2` (the 0/1-dimensional
    /// "spheres" have empty/degenerate tangent spaces), `Stiefel` needs
    /// `1 <= p <= n`, and every point/tangent formula must stay inside
    /// the checked `u32` domain and the per-variable dimension cap.
    ///
    /// # Errors
    /// [`OptError::ManifoldInvalid`] naming the violated rule.
    pub fn validate(&self, caps: &AdmissionCaps) -> Result<(), OptError> {
        admission::validate_manifold(self, caps)
    }
}

/// Value shape of an expression node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    /// One number.
    Scalar,
    /// A fixed-length vector.
    Vector(u32),
}

/// Differentiability class, propagated bottom-up (the minimum over a
/// node's children and its own contribution).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Class {
    /// Non-smooth (kinks: min/max/abs, CVaR, quantiles).
    C0,
    /// Once differentiable.
    C1,
    /// Smooth on its domain.
    Smooth,
}

/// One expression node. `f64` payloads serialize by BIT PATTERN, so
/// identity (hash-consing, problem hashes, round-trips) is exact.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// A variable reference (shape = its manifold's point storage).
    Var(VarId),
    /// One component of a vector-valued node.
    Component {
        /// Source vector node.
        of: NodeId,
        /// Component index.
        index: u32,
    },
    /// A dimensioned constant.
    Const {
        /// Value in coherent SI units.
        value: f64,
        /// Dimension vector.
        dims: Dims,
    },
    /// Sum (shapes and dimensions must match).
    Add(NodeId, NodeId),
    /// Difference (shapes and dimensions must match).
    Sub(NodeId, NodeId),
    /// Product (scalar×scalar or scalar×vector; dimensions add).
    Mul(NodeId, NodeId),
    /// Quotient (scalars; dimensions subtract).
    Div(NodeId, NodeId),
    /// Negation.
    Neg(NodeId),
    /// Integer power (scalar; dimensions scale by the exponent).
    Powi {
        /// Base node.
        base: NodeId,
        /// Exponent.
        exp: i32,
    },
    /// Square root (scalar; even dimension exponents halve).
    Sqrt(NodeId),
    /// Exponential (dimensionless scalar).
    Exp(NodeId),
    /// Natural log (dimensionless scalar).
    Ln(NodeId),
    /// Hyperbolic tangent (dimensionless scalar).
    Tanh(NodeId),
    /// Inner product of same-length vectors (dimensions add).
    Dot(NodeId, NodeId),
    /// Squared Euclidean norm of a vector (dimensions double).
    NormSq(NodeId),
    /// Pointwise minimum — C0: POISONS smooth-optimizer routing.
    Min(NodeId, NodeId),
    /// Pointwise maximum — C0.
    Max(NodeId, NodeId),
    /// Absolute value — C0.
    Abs(NodeId),
    /// A PDE residual node `physics(u, θ) = 0` referencing a FLUX study
    /// (first-class, with adjoint availability metadata).
    PdeResidual {
        /// Study identifier (FLUX side).
        study: String,
        /// The design variable the physics depends on.
        over: VarId,
        /// Whether an adjoint gradient path exists.
        adjoint_available: bool,
        /// Declared dimensions of the residual.
        dims: Dims,
    },
    /// Expectation over a UQ configuration (preserves the child class).
    Expectation {
        /// Integrand.
        of: NodeId,
        /// UQ configuration identifier.
        uq_config: String,
    },
    /// Conditional value-at-risk — C0 (kink at the VaR).
    Cvar {
        /// Integrand.
        of: NodeId,
        /// Tail level in (0, 1).
        alpha: f64,
        /// UQ configuration identifier.
        uq_config: String,
    },
    /// Quantile — C0.
    Quantile {
        /// Integrand.
        of: NodeId,
        /// Quantile level in (0, 1).
        q: f64,
        /// UQ configuration identifier.
        uq_config: String,
    },
}

/// Allocation-free hash-consing lookup key. External strings contribute
/// a full-width content hash; exact [`Expr`] comparison inside the bucket
/// resolves the theoretical collision case without treating a hash as
/// semantic equality.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct ExprFingerprint {
    kind: u8,
    words: [u64; 3],
    dims: [i8; 6],
    text: ContentHash,
}

impl ExprFingerprint {
    const EMPTY_TEXT: ContentHash = ContentHash([0; 32]);

    fn plain(kind: u8, words: [u64; 3], dims: [i8; 6]) -> Self {
        Self {
            kind,
            words,
            dims,
            text: Self::EMPTY_TEXT,
        }
    }

    fn with_text(kind: u8, words: [u64; 3], dims: [i8; 6], text: &str) -> Self {
        Self {
            kind,
            words,
            dims,
            text: hash_bytes(text.as_bytes()),
        }
    }

    fn of(expression: &Expr) -> Self {
        match expression {
            Expr::Var(variable) => Self::plain(0, [u64::from(variable.0), 0, 0], [0; 6]),
            Expr::Component { of, index } => {
                Self::plain(1, [u64::from(of.0), u64::from(*index), 0], [0; 6])
            }
            Expr::Const { value, dims } => Self::plain(2, [value.to_bits(), 0, 0], dims.0),
            Expr::Add(a, b) => Self::binary(3, *a, *b),
            Expr::Sub(a, b) => Self::binary(4, *a, *b),
            Expr::Mul(a, b) => Self::binary(5, *a, *b),
            Expr::Div(a, b) => Self::binary(6, *a, *b),
            Expr::Neg(node) => Self::unary(7, *node),
            Expr::Powi { base, exp } => {
                Self::plain(8, [u64::from(base.0), i64::from(*exp) as u64, 0], [0; 6])
            }
            Expr::Sqrt(node) => Self::unary(9, *node),
            Expr::Exp(node) => Self::unary(10, *node),
            Expr::Ln(node) => Self::unary(11, *node),
            Expr::Tanh(node) => Self::unary(12, *node),
            Expr::Dot(a, b) => Self::binary(13, *a, *b),
            Expr::NormSq(node) => Self::unary(14, *node),
            Expr::Min(a, b) => Self::binary(15, *a, *b),
            Expr::Max(a, b) => Self::binary(16, *a, *b),
            Expr::Abs(node) => Self::unary(17, *node),
            Expr::PdeResidual {
                study,
                over,
                adjoint_available,
                dims,
            } => Self::pde(study, *over, *adjoint_available, *dims),
            Expr::Expectation { of, uq_config } => Self::expectation(*of, uq_config),
            Expr::Cvar {
                of,
                alpha,
                uq_config,
            } => Self::cvar(*of, *alpha, uq_config),
            Expr::Quantile { of, q, uq_config } => Self::quantile(*of, *q, uq_config),
        }
    }

    fn unary(kind: u8, node: NodeId) -> Self {
        Self::plain(kind, [u64::from(node.0), 0, 0], [0; 6])
    }

    fn binary(kind: u8, left: NodeId, right: NodeId) -> Self {
        Self::plain(kind, [u64::from(left.0), u64::from(right.0), 0], [0; 6])
    }

    fn pde(study: &str, over: VarId, adjoint_available: bool, dims: Dims) -> Self {
        Self::with_text(
            18,
            [u64::from(over.0), u64::from(u8::from(adjoint_available)), 0],
            dims.0,
            study,
        )
    }

    fn expectation(of: NodeId, uq_config: &str) -> Self {
        Self::with_text(19, [u64::from(of.0), 0, 0], [0; 6], uq_config)
    }

    fn cvar(of: NodeId, alpha: f64, uq_config: &str) -> Self {
        Self::with_text(20, [u64::from(of.0), alpha.to_bits(), 0], [0; 6], uq_config)
    }

    fn quantile(of: NodeId, q: f64, uq_config: &str) -> Self {
        Self::with_text(21, [u64::from(of.0), q.to_bits(), 0], [0; 6], uq_config)
    }
}

#[derive(Debug)]
enum InternEntry {
    One(NodeId),
    Collision(Vec<NodeId>),
}

/// Objective sense.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sense {
    /// Minimize.
    Minimize,
    /// Maximize.
    Maximize,
}

/// Constraint kind (semantics/repair live in fs-constraint; the IR
/// owns the graph substrate only).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintKind {
    /// `g(x) = 0`.
    EqZero,
    /// `g(x) ≤ 0`.
    LeZero,
}

/// A reference to a bilevel inner problem. The two variants are
/// deliberately NON-INTERCHANGEABLE: a full-width BLAKE3-backed
/// [`ProblemSemanticId`] is the only strong spelling, while a legacy
/// 64-bit FNV identity parsed from a historical artifact stays
/// QUARANTINED — it is inspectable provenance with no execution or
/// certificate authority, and no API widens or reinterprets it as
/// strong (admission lists it in
/// [`ProblemAdmission::quarantined_legacy_identities`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BilevelRef {
    /// Full-width semantic identity of the inner problem.
    Semantic(ProblemSemanticId),
    /// Quarantined legacy FNV-64 identity from a v1/v2 artifact.
    LegacyFnv(LegacyProblemHash),
}

/// Problem-structure annotations (representable, not yet executed).
#[derive(Debug, Clone, PartialEq)]
pub enum ProblemTag {
    /// Multiple fidelity levels are available.
    MultiFidelity {
        /// Number of levels (validated: `1..=` the admission cap).
        levels: u32,
    },
    /// Chance constraint: `P(g ≤ 0) ≥ prob`.
    ChanceConstrained {
        /// Required probability (validated: finite, strictly in (0, 1)).
        prob: f64,
    },
    /// Bilevel structure (inner problem referenced by typed identity).
    Bilevel {
        /// Inner problem reference.
        inner: BilevelRef,
    },
}

/// Evaluation budget (P4: attached to the problem, enforced by
/// consumers like the toy descent).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EvalBudget {
    /// Maximum objective evaluations (0 = unlimited).
    pub max_evals: u64,
}

/// One declared variable.
#[derive(Debug, Clone, PartialEq)]
pub struct Variable {
    /// Human name (diagnostics).
    pub name: String,
    /// Where it lives.
    pub manifold: Manifold,
    /// Dimensions of its components.
    pub dims: Dims,
}

/// One objective entry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Objective {
    /// Root node (Scalar).
    pub node: NodeId,
    /// Direction.
    pub sense: Sense,
    /// Multi-objective weight (1.0 default).
    pub weight: f64,
}

/// One constraint entry.
#[derive(Debug, Clone, PartialEq)]
pub struct Constraint {
    /// Root node (Scalar).
    pub node: NodeId,
    /// Kind.
    pub kind: ConstraintKind,
    /// Human name (diagnostics).
    pub name: String,
}

/// Teaching errors for graph construction and routing.
#[derive(Debug, Clone, PartialEq)]
pub enum OptError {
    /// A referenced node id does not exist.
    UnknownNode {
        /// The offending id.
        id: u32,
    },
    /// A referenced variable id does not exist.
    UnknownVar {
        /// The offending id.
        id: u32,
    },
    /// Operand shapes are incompatible.
    ShapeMismatch {
        /// Operation name.
        op: &'static str,
        /// Left/first shape.
        left: Shape,
        /// Right/second shape (or expected).
        right: Shape,
    },
    /// Operand dimensions are incompatible.
    DimMismatch {
        /// Operation name.
        op: &'static str,
        /// Left dims (fs-qty exponent vector).
        left: [i8; 6],
        /// Right dims.
        right: [i8; 6],
    },
    /// Combining two operand dimension vectors would leave the
    /// representable i8 exponent domain (mul/div/dot/norm_sq).
    DimSumOverflow {
        /// Operation that combined the dimensions.
        op: &'static str,
        /// Left operand dims.
        left: [i8; 6],
        /// Right operand dims.
        right: [i8; 6],
    },
    /// Dimension exponents would leave the representable runtime domain.
    DimOverflow {
        /// Operation that attempted the dimension arithmetic.
        op: &'static str,
        /// Base dimensions before scaling.
        dims: [i8; 6],
        /// Requested integer exponent.
        exponent: i32,
    },
    /// A transcendental was fed a dimensioned quantity.
    NonDimensionless {
        /// Operation name.
        op: &'static str,
        /// The offending dims.
        dims: [i8; 6],
    },
    /// `sqrt` of odd dimension exponents.
    OddDims {
        /// The offending dims.
        dims: [i8; 6],
    },
    /// A parameter left its valid range.
    BadParam {
        /// What.
        what: &'static str,
        /// Value received.
        value: f64,
    },
    /// A vector component index out of range.
    IndexOut {
        /// Index asked.
        index: u32,
        /// Vector length.
        len: u32,
    },
    /// An objective/constraint root must be scalar.
    NotScalar {
        /// The offending node.
        node: u32,
    },
    /// The problem is non-smooth for the requested optimizer family.
    NonsmoothForFamily {
        /// Requested family.
        family: &'static str,
        /// The node that poisons smoothness.
        node: u32,
        /// What that node is.
        kind: String,
        /// Its propagated class.
        class: Class,
    },
    /// A PDE node lacks an adjoint path for a gradient-based family.
    NoAdjoint {
        /// The PDE node.
        node: u32,
        /// The study it references.
        study: String,
    },
    /// Serialized text failed to parse.
    Parse {
        /// 1-based line number.
        line: usize,
        /// What went wrong.
        what: String,
    },
    /// A problem cannot be represented by a requested historical wire
    /// version without weakening or changing its typed semantics.
    WireIncompatible {
        /// Historical wire version requested.
        version: &'static str,
        /// The incompatible construct.
        what: String,
    },
    /// A node the IR carries but cannot execute (PDE/stochastic).
    Unevaluable {
        /// The node.
        node: u32,
        /// What it is and who executes it.
        kind: String,
    },
    /// Cancelled mid-run (descent).
    Cancelled,
    /// Budget exhausted (P4 receipt, not a failure of the math).
    BudgetExhausted {
        /// Evaluations spent.
        spent: u64,
    },
    /// A manifold descriptor violates the admission policy.
    ManifoldInvalid {
        /// The violated rule, teaching text.
        what: String,
    },
    /// A payload that must be finite is NaN or infinite. The exact bit
    /// pattern is retained because NaN payloads do not survive Display.
    NonFinite {
        /// Which payload.
        what: &'static str,
        /// The offending value's IEEE-754 bit pattern.
        bits: u64,
    },
    /// A runtime variable binding contains a NaN or infinity.
    BindingNonFinite {
        /// The variable whose point is malformed.
        var: u32,
        /// The component within that manifold point.
        component: u32,
        /// The offending value's exact IEEE-754 bit pattern.
        bits: u64,
    },
    /// Runtime arithmetic produced a NaN or infinity at a graph node.
    EvalNonFinite {
        /// The node whose result is malformed.
        node: u32,
        /// Vector component, or `None` for a scalar result.
        component: Option<u32>,
        /// The offending value's exact IEEE-754 bit pattern.
        bits: u64,
    },
    /// A raw retraction point or step has the wrong storage length.
    RetractionLen {
        /// Which retraction input was malformed.
        input: &'static str,
        /// Required storage length for this manifold.
        expected: u32,
        /// Storage length supplied by the caller.
        got: u64,
    },
    /// A per-item or aggregate admission cap was exceeded.
    CapExceeded {
        /// Which cap.
        what: &'static str,
        /// Count that was attempted.
        count: u64,
        /// The versioned cap.
        cap: u64,
    },
    /// More variable bindings were supplied than the problem declares.
    BindingCount {
        /// Declared variable count.
        vars: u32,
        /// Bindings supplied.
        got: u64,
    },
    /// A supplied binding's length does not match its manifold storage.
    BindingLen {
        /// The variable.
        var: u32,
        /// Its manifold point storage length.
        expected: u32,
        /// The binding length supplied.
        got: u64,
    },
}

impl core::fmt::Display for OptError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            OptError::UnknownNode { id } => write!(
                f,
                "node id {id} does not exist in this problem; use ids returned by the builder"
            ),
            OptError::UnknownVar { id } => write!(f, "variable id {id} does not exist"),
            OptError::ShapeMismatch { op, left, right } => write!(
                f,
                "`{op}` got incompatible shapes {left:?} and {right:?}; scalars and \
                 vectors do not mix here"
            ),
            OptError::DimMismatch { op, left, right } => write!(
                f,
                "`{op}` got incompatible dimensions {left:?} vs {right:?}; only \
                 same-dimension quantities add or compare"
            ),
            OptError::DimSumOverflow { op, left, right } => write!(
                f,
                "`{op}` would combine dimensions {left:?} and {right:?} outside the \
                 supported i8 exponent domain; rescale the formulation"
            ),
            OptError::DimOverflow { op, dims, exponent } => write!(
                f,
                "`{op}` would scale dimensions {dims:?} by {exponent} outside the supported i8 exponent domain; rescale the formulation or use a dimensionless base"
            ),
            OptError::NonDimensionless { op, dims } => write!(
                f,
                "`{op}` needs a dimensionless argument, got exponents {dims:?}; divide \
                 by a reference quantity first"
            ),
            OptError::OddDims { dims } => write!(
                f,
                "sqrt of odd dimension exponents {dims:?} has no dimensional meaning"
            ),
            OptError::BadParam { what, value } => {
                write!(f, "`{what}` = {value} is outside its valid range")
            }
            OptError::IndexOut { index, len } => {
                write!(
                    f,
                    "component {index} of a length-{len} vector does not exist"
                )
            }
            OptError::NotScalar { node } => write!(
                f,
                "objective/constraint roots must be SCALAR; node {node} is a vector — \
                 reduce it (dot, norm_sq, component) first"
            ),
            OptError::NonsmoothForFamily {
                family,
                node,
                kind,
                class,
            } => write!(
                f,
                "this problem is non-smooth for {family}: node {node} ({kind}) has \
                 class {class:?} — route to a subgradient/gradient-free family or \
                 replace the kink with a smooth surrogate"
            ),
            OptError::NoAdjoint { node, study } => write!(
                f,
                "PDE node {node} (study `{study}`) has no adjoint path; a \
                 gradient-based family cannot differentiate through it"
            ),
            OptError::Unevaluable { node, kind } => write!(
                f,
                "node {node} is carried by the IR but not evaluable here: {kind}"
            ),
            OptError::Parse { line, what } => write!(f, "parse error at line {line}: {what}"),
            OptError::WireIncompatible { version, what } => {
                write!(f, "problem cannot be encoded as {version}: {what}")
            }
            OptError::Cancelled => write!(f, "cancelled between descent steps"),
            OptError::BudgetExhausted { spent } => write!(
                f,
                "evaluation budget exhausted after {spent} evaluations (P4 receipt)"
            ),
            OptError::ManifoldInvalid { what } => {
                write!(f, "manifold descriptor rejected: {what}")
            }
            OptError::NonFinite { what, bits } => write!(
                f,
                "`{what}` must be finite; got {} (bits {bits:016X}) — non-finite \
                 payloads cannot carry graph authority",
                f64::from_bits(*bits)
            ),
            OptError::BindingNonFinite {
                var,
                component,
                bits,
            } => write!(
                f,
                "binding for variable {var} has non-finite component {component}: {} \
                 (bits {bits:016X}); runtime points must be finite before evaluation",
                f64::from_bits(*bits)
            ),
            OptError::EvalNonFinite {
                node,
                component,
                bits,
            } => match component {
                Some(component) => write!(
                    f,
                    "node {node} produced non-finite vector component {component}: {} \
                     (bits {bits:016X}); the result has no evaluation authority",
                    f64::from_bits(*bits)
                ),
                None => write!(
                    f,
                    "node {node} produced a non-finite scalar: {} (bits {bits:016X}); \
                     the result has no evaluation authority",
                    f64::from_bits(*bits)
                ),
            },
            OptError::RetractionLen {
                input,
                expected,
                got,
            } => write!(
                f,
                "{input} has length {got}, but this manifold requires length {expected}; \
                 retraction refuses malformed storage"
            ),
            OptError::CapExceeded { what, count, cap } => write!(
                f,
                "admission cap exceeded: {what} = {count} > {cap}; split the problem \
                 or raise the cap through an explicit AdmissionCaps"
            ),
            OptError::BindingCount { vars, got } => write!(
                f,
                "{got} bindings supplied for {vars} declared variable(s); bindings \
                 are indexed by VarId and must not exceed the declaration list"
            ),
            OptError::BindingLen { var, expected, got } => write!(
                f,
                "binding for variable {var} has length {got}, but its manifold \
                 stores points of length {expected}"
            ),
        }
    }
}

impl std::error::Error for OptError {}

/// Optimizer families for routing checks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimizerFamily {
    /// Quasi-Newton (needs C1+ and adjoints through PDE nodes).
    Lbfgs,
    /// Newton-type (needs Smooth and adjoints).
    Newton,
    /// Subgradient/bundle methods (accepts C0).
    SubgradientBundle,
    /// Gradient-free (accepts anything).
    GradientFree,
}

/// The finished problem: graph + roots + metadata. SEALED — fields are
/// crate-private, so the only public paths into a `Problem` are the
/// validating [`ProblemBuilder`] and the parser (which rebuilds through
/// the builder). Accessors are read-only and cheap; ID-indexed
/// accessors are CHECKED and refuse unknown ids instead of panicking.
#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    pub(crate) vars: Vec<Variable>,
    pub(crate) exprs: Vec<Expr>,
    pub(crate) objectives: Vec<Objective>,
    pub(crate) constraints: Vec<Constraint>,
    pub(crate) tags: Vec<ProblemTag>,
    pub(crate) budget: EvalBudget,
    pub(crate) shapes: Vec<Shape>,
    pub(crate) dims: Vec<Dims>,
    pub(crate) classes: Vec<Class>,
    pub(crate) depths: Vec<u32>,
    pub(crate) total_work: u64,
}

impl Problem {
    /// Declared variables (read-only).
    #[must_use]
    pub fn vars(&self) -> &[Variable] {
        &self.vars
    }

    /// Expression nodes (read-only; ids index this list).
    #[must_use]
    pub fn exprs(&self) -> &[Expr] {
        &self.exprs
    }

    /// Objectives (read-only).
    #[must_use]
    pub fn objectives(&self) -> &[Objective] {
        &self.objectives
    }

    /// Constraints (read-only).
    #[must_use]
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// Structure annotations (read-only).
    #[must_use]
    pub fn tags(&self) -> &[ProblemTag] {
        &self.tags
    }

    /// Evaluation budget (P4).
    #[must_use]
    pub fn budget(&self) -> EvalBudget {
        self.budget
    }

    /// Checked expression accessor.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn expr(&self, n: NodeId) -> Result<&Expr, OptError> {
        self.exprs
            .get(n.0 as usize)
            .ok_or(OptError::UnknownNode { id: n.0 })
    }

    /// Checked variable accessor.
    ///
    /// # Errors
    /// [`OptError::UnknownVar`].
    pub fn variable(&self, v: VarId) -> Result<&Variable, OptError> {
        self.vars
            .get(v.0 as usize)
            .ok_or(OptError::UnknownVar { id: v.0 })
    }

    /// Shape of a node (checked).
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn shape(&self, n: NodeId) -> Result<Shape, OptError> {
        self.shapes
            .get(n.0 as usize)
            .copied()
            .ok_or(OptError::UnknownNode { id: n.0 })
    }

    /// Dimensions of a node (checked).
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn node_dims(&self, n: NodeId) -> Result<Dims, OptError> {
        self.dims
            .get(n.0 as usize)
            .copied()
            .ok_or(OptError::UnknownNode { id: n.0 })
    }

    /// Propagated differentiability class of a node (checked).
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn class(&self, n: NodeId) -> Result<Class, OptError> {
        self.classes
            .get(n.0 as usize)
            .copied()
            .ok_or(OptError::UnknownNode { id: n.0 })
    }

    /// Root-to-leaf depth of a node, derived and capped before the node
    /// entered the sealed graph.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn node_depth(&self, n: NodeId) -> Result<u32, OptError> {
        self.depths
            .get(n.0 as usize)
            .copied()
            .ok_or(OptError::UnknownNode { id: n.0 })
    }

    /// Deterministic admission-work charge retained with the sealed
    /// graph. Re-admission re-derives this value instead of trusting it.
    #[must_use]
    pub fn total_admission_work(&self) -> u64 {
        self.total_work
    }

    /// Re-validate the complete problem through the versioned common
    /// admission validator and mint its [`ProblemSemanticId`]. Builder
    /// output always admits (the builder enforces the same leaf rules);
    /// this is the defense-in-depth chokepoint for deserialized,
    /// migrated, or future foreign constructions.
    ///
    /// # Errors
    /// A complete, deterministically ordered [`AdmissionReport`].
    pub fn admit(&self) -> Result<ProblemAdmission, AdmissionReport> {
        admission::admit_with_caps(self, &AdmissionCaps::default())
    }

    /// [`Problem::admit`] under explicit caps.
    ///
    /// # Errors
    /// A complete, deterministically ordered [`AdmissionReport`].
    pub fn admit_with_caps(
        &self,
        caps: &AdmissionCaps,
    ) -> Result<ProblemAdmission, AdmissionReport> {
        admission::admit_with_caps(self, caps)
    }

    /// The class-propagation trace: one line per node (build order),
    /// naming each node's own contribution and resulting class.
    #[must_use]
    pub fn class_trace(&self) -> Vec<String> {
        self.exprs
            .iter()
            .enumerate()
            .map(|(i, e)| format!("node {i}: {} -> {:?}", expr_kind_name(e), self.classes[i]))
            .collect()
    }

    /// Route the problem to an optimizer family, refusing with the
    /// OFFENDING NODE named when the class or adjoint metadata does not
    /// support it.
    ///
    /// # Errors
    /// [`OptError::NonsmoothForFamily`] / [`OptError::NoAdjoint`].
    pub fn route(&self, family: OptimizerFamily) -> Result<(), OptError> {
        let min_class = match family {
            OptimizerFamily::Lbfgs => Class::C1,
            OptimizerFamily::Newton => Class::Smooth,
            OptimizerFamily::SubgradientBundle | OptimizerFamily::GradientFree => Class::C0,
        };
        let needs_adjoint = matches!(family, OptimizerFamily::Lbfgs | OptimizerFamily::Newton);
        let family_name = match family {
            OptimizerFamily::Lbfgs => "L-BFGS",
            OptimizerFamily::Newton => "Newton",
            OptimizerFamily::SubgradientBundle => "subgradient/bundle",
            OptimizerFamily::GradientFree => "gradient-free",
        };
        let roots: Vec<NodeId> = self
            .objectives
            .iter()
            .map(|o| o.node)
            .chain(self.constraints.iter().map(|c| c.node))
            .collect();
        for root in roots {
            for n in self.reachable(root)? {
                let i = n.0 as usize;
                if self.classes[i] < min_class && own_class(&self.exprs[i]) < min_class {
                    return Err(OptError::NonsmoothForFamily {
                        family: family_name,
                        node: n.0,
                        kind: expr_kind_name(&self.exprs[i]).to_string(),
                        class: self.classes[i],
                    });
                }
                if needs_adjoint
                    && let Expr::PdeResidual {
                        study,
                        adjoint_available: false,
                        ..
                    } = &self.exprs[i]
                {
                    return Err(OptError::NoAdjoint {
                        node: n.0,
                        study: study.clone(),
                    });
                }
            }
        }
        Ok(())
    }

    /// Nodes reachable from `root` (build order). The root is CHECKED;
    /// interior child ids are valid by the sealed-arena invariant.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`] when `root` is not in this problem.
    pub fn reachable(&self, root: NodeId) -> Result<Vec<NodeId>, OptError> {
        if root.0 as usize >= self.exprs.len() {
            return Err(OptError::UnknownNode { id: root.0 });
        }
        let mut seen = vec![false; self.exprs.len()];
        let mut stack = vec![root];
        let mut out = Vec::new();
        while let Some(n) = stack.pop() {
            if std::mem::replace(&mut seen[n.0 as usize], true) {
                continue;
            }
            out.push(n);
            for c in children(&self.exprs[n.0 as usize]) {
                stack.push(c);
            }
        }
        out.sort_unstable();
        Ok(out)
    }
}

/// Child node ids of an expression. The fixed-size representation keeps
/// validation and cap checks allocation-free. Public so downstream
/// evaluators (fs-constraint's interval walker, future consumers) can
/// drive EXPLICIT-STACK traversals instead of recursing (bead
/// frankensim-xf8v7: recursion over admitted graphs is a stack-overflow
/// hazard at the depth cap; iterate instead).
pub fn children(e: &Expr) -> impl Iterator<Item = NodeId> {
    let children = match *e {
        Expr::Var(_) | Expr::Const { .. } | Expr::PdeResidual { .. } => [None, None],
        Expr::Component { of, .. }
        | Expr::Neg(of)
        | Expr::Powi { base: of, .. }
        | Expr::Sqrt(of)
        | Expr::Exp(of)
        | Expr::Ln(of)
        | Expr::Tanh(of)
        | Expr::NormSq(of)
        | Expr::Abs(of)
        | Expr::Expectation { of, .. }
        | Expr::Cvar { of, .. }
        | Expr::Quantile { of, .. } => [Some(of), None],
        Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Dot(a, b)
        | Expr::Min(a, b)
        | Expr::Max(a, b) => [Some(a), Some(b)],
    };
    children.into_iter().flatten()
}

/// A node's OWN class contribution (children aside).
pub(crate) fn own_class(e: &Expr) -> Class {
    match e {
        Expr::Min(..)
        | Expr::Max(..)
        | Expr::Abs(_)
        | Expr::Cvar { .. }
        | Expr::Quantile { .. } => Class::C0,
        _ => Class::Smooth,
    }
}

/// Exact hash-consing equality. Float payloads compare by bit pattern,
/// matching canonical wire identity (`-0.0` and `0.0` remain distinct).
fn expr_identity_eq(left: &Expr, right: &Expr) -> bool {
    match (left, right) {
        (Expr::Var(a), Expr::Var(b)) => a == b,
        (
            Expr::Component {
                of: a_of,
                index: a_index,
            },
            Expr::Component {
                of: b_of,
                index: b_index,
            },
        ) => a_of == b_of && a_index == b_index,
        (
            Expr::Const {
                value: a_value,
                dims: a_dims,
            },
            Expr::Const {
                value: b_value,
                dims: b_dims,
            },
        ) => a_value.to_bits() == b_value.to_bits() && a_dims == b_dims,
        (Expr::Add(a1, a2), Expr::Add(b1, b2))
        | (Expr::Sub(a1, a2), Expr::Sub(b1, b2))
        | (Expr::Mul(a1, a2), Expr::Mul(b1, b2))
        | (Expr::Div(a1, a2), Expr::Div(b1, b2))
        | (Expr::Dot(a1, a2), Expr::Dot(b1, b2))
        | (Expr::Min(a1, a2), Expr::Min(b1, b2))
        | (Expr::Max(a1, a2), Expr::Max(b1, b2)) => a1 == b1 && a2 == b2,
        (Expr::Neg(a), Expr::Neg(b))
        | (Expr::Sqrt(a), Expr::Sqrt(b))
        | (Expr::Exp(a), Expr::Exp(b))
        | (Expr::Ln(a), Expr::Ln(b))
        | (Expr::Tanh(a), Expr::Tanh(b))
        | (Expr::NormSq(a), Expr::NormSq(b))
        | (Expr::Abs(a), Expr::Abs(b)) => a == b,
        (
            Expr::Powi {
                base: a_base,
                exp: a_exp,
            },
            Expr::Powi {
                base: b_base,
                exp: b_exp,
            },
        ) => a_base == b_base && a_exp == b_exp,
        (
            Expr::PdeResidual {
                study: a_study,
                over: a_over,
                adjoint_available: a_adjoint,
                dims: a_dims,
            },
            Expr::PdeResidual {
                study: b_study,
                over: b_over,
                adjoint_available: b_adjoint,
                dims: b_dims,
            },
        ) => a_study == b_study && a_over == b_over && a_adjoint == b_adjoint && a_dims == b_dims,
        (
            Expr::Expectation {
                of: a_of,
                uq_config: a_config,
            },
            Expr::Expectation {
                of: b_of,
                uq_config: b_config,
            },
        ) => a_of == b_of && a_config == b_config,
        (
            Expr::Cvar {
                of: a_of,
                alpha: a_alpha,
                uq_config: a_config,
            },
            Expr::Cvar {
                of: b_of,
                alpha: b_alpha,
                uq_config: b_config,
            },
        ) => a_of == b_of && a_alpha.to_bits() == b_alpha.to_bits() && a_config == b_config,
        (
            Expr::Quantile {
                of: a_of,
                q: a_q,
                uq_config: a_config,
            },
            Expr::Quantile {
                of: b_of,
                q: b_q,
                uq_config: b_config,
            },
        ) => a_of == b_of && a_q.to_bits() == b_q.to_bits() && a_config == b_config,
        _ => false,
    }
}

/// Stable kind name (diagnostics, traces, serialization).
pub(crate) fn expr_kind_name(e: &Expr) -> &'static str {
    match e {
        Expr::Var(_) => "var",
        Expr::Component { .. } => "component",
        Expr::Const { .. } => "const",
        Expr::Add(..) => "add",
        Expr::Sub(..) => "sub",
        Expr::Mul(..) => "mul",
        Expr::Div(..) => "div",
        Expr::Neg(_) => "neg",
        Expr::Powi { .. } => "powi",
        Expr::Sqrt(_) => "sqrt",
        Expr::Exp(_) => "exp",
        Expr::Ln(_) => "ln",
        Expr::Tanh(_) => "tanh",
        Expr::Dot(..) => "dot",
        Expr::NormSq(_) => "norm_sq",
        Expr::Min(..) => "min",
        Expr::Max(..) => "max",
        Expr::Abs(_) => "abs",
        Expr::PdeResidual { .. } => "pde_residual",
        Expr::Expectation { .. } => "expectation",
        Expr::Cvar { .. } => "cvar",
        Expr::Quantile { .. } => "quantile",
    }
}

/// Incremental, validating builder. Every `Result` is a teaching error;
/// hash-consing makes repeated subexpressions return the SAME id (the
/// G0 common-subexpression identity). Every constructor validates
/// through the SAME versioned leaf rules the admission validator uses
/// (`admission::derive_expr` and friends), and every rejection leaves
/// the intern table, ids, storage/byte/work accounting, budget, and
/// ordering unchanged.
#[derive(Debug)]
pub struct ProblemBuilder {
    vars: Vec<Variable>,
    exprs: Vec<Expr>,
    shapes: Vec<Shape>,
    dims: Vec<Dims>,
    classes: Vec<Class>,
    depths: Vec<u32>,
    intern: BTreeMap<ExprFingerprint, InternEntry>,
    objectives: Vec<Objective>,
    constraints: Vec<Constraint>,
    tags: Vec<ProblemTag>,
    budget: EvalBudget,
    caps: AdmissionCaps,
    total_point_storage: u64,
    total_retained_bytes: u64,
    total_work: u64,
}

impl Default for ProblemBuilder {
    fn default() -> Self {
        ProblemBuilder::new()
    }
}

impl ProblemBuilder {
    /// Empty builder (unlimited budget, versioned default caps).
    #[must_use]
    pub fn new() -> Self {
        ProblemBuilder::with_caps(AdmissionCaps::default())
    }

    /// Empty builder under explicit admission caps (tests, sandboxes).
    #[must_use]
    pub fn with_caps(caps: AdmissionCaps) -> Self {
        ProblemBuilder {
            vars: Vec::new(),
            exprs: Vec::new(),
            shapes: Vec::new(),
            dims: Vec::new(),
            classes: Vec::new(),
            depths: Vec::new(),
            intern: BTreeMap::new(),
            objectives: Vec::new(),
            constraints: Vec::new(),
            tags: Vec::new(),
            budget: EvalBudget { max_evals: 0 },
            caps,
            total_point_storage: 0,
            total_retained_bytes: 0,
            total_work: 0,
        }
    }

    /// Declare a variable. Validates the manifold descriptor (checked
    /// point/tangent formulas, per-variable dimension cap), the name
    /// length, and the variable-count, total-storage, target-address,
    /// retained-byte, and aggregate-work caps BEFORE assigning a `VarId`.
    ///
    /// # Errors
    /// [`OptError::ManifoldInvalid`] / [`OptError::CapExceeded`].
    pub fn var(&mut self, name: &str, manifold: Manifold, dims: Dims) -> Result<VarId, OptError> {
        admission::validate_name("variable name", name, &self.caps)?;
        admission::validate_manifold(&manifold, &self.caps)?;
        if self.vars.len() as u64 >= u64::from(self.caps.max_vars) {
            return Err(OptError::CapExceeded {
                what: "variables",
                count: self.vars.len() as u64 + 1,
                cap: u64::from(self.caps.max_vars),
            });
        }
        let storage = u64::from(
            manifold
                .point_dim()
                .expect("validate_manifold proved the point formula representable"),
        );
        let total = self
            .total_point_storage
            .checked_add(storage)
            .unwrap_or(u64::MAX);
        if total > self.caps.max_total_point_storage {
            return Err(OptError::CapExceeded {
                what: "total point storage",
                count: total,
                cap: self.caps.max_total_point_storage,
            });
        }
        if usize::try_from(total).is_err() {
            return Err(OptError::CapExceeded {
                what: "target packed point storage",
                count: total,
                cap: usize::MAX as u64,
            });
        }
        let total_retained_bytes = admission::checked_retained_total(
            self.total_retained_bytes,
            admission::name_retained_bytes(name),
            &self.caps,
        )?;
        let total_work = admission::checked_work_total(self.total_work, 1, &self.caps)?;
        self.total_point_storage = total;
        self.total_retained_bytes = total_retained_bytes;
        self.total_work = total_work;
        self.vars.push(Variable {
            name: name.to_string(),
            manifold,
            dims,
        });
        Ok(VarId((self.vars.len() - 1) as u32))
    }

    /// Attach the evaluation budget (P4). Any `u64` is valid (0 =
    /// unlimited), so this stays infallible.
    pub fn set_budget(&mut self, max_evals: u64) {
        self.budget.max_evals = max_evals;
    }

    /// Attach a structure tag (validated: fidelity levels in `1..=cap`,
    /// chance probability finite and strictly inside (0, 1), bilevel
    /// references typed).
    ///
    /// # Errors
    /// [`OptError::BadParam`] / [`OptError::CapExceeded`].
    pub fn tag(&mut self, tag: ProblemTag) -> Result<(), OptError> {
        admission::validate_tag(&tag, &self.caps)?;
        if self.tags.len() as u64 >= u64::from(self.caps.max_tags) {
            return Err(OptError::CapExceeded {
                what: "tags",
                count: self.tags.len() as u64 + 1,
                cap: u64::from(self.caps.max_tags),
            });
        }
        let total_work = admission::checked_work_total(self.total_work, 1, &self.caps)?;
        self.tags.push(tag);
        self.total_work = total_work;
        Ok(())
    }

    fn require_scalar_root(&self, node: NodeId) -> Result<(), OptError> {
        let shape = self
            .shapes
            .get(node.0 as usize)
            .ok_or(OptError::UnknownNode { id: node.0 })?;
        if *shape != Shape::Scalar {
            return Err(OptError::NotScalar { node: node.0 });
        }
        Ok(())
    }

    fn unary_depth(&self, node: NodeId) -> Result<u32, OptError> {
        let probe = Expr::Neg(node);
        let lookup = |child: NodeId| self.depths.get(child.0 as usize).copied();
        admission::derive_depth(&probe, &lookup, &self.caps)
    }

    fn new_expr_totals(
        &self,
        retained_charge: u64,
        work_charge: u64,
    ) -> Result<(u64, u64), OptError> {
        if self.exprs.len() as u64 >= u64::from(self.caps.max_nodes) {
            return Err(OptError::CapExceeded {
                what: "expression nodes",
                count: self.exprs.len() as u64 + 1,
                cap: u64::from(self.caps.max_nodes),
            });
        }
        let total_work = admission::checked_work_total(self.total_work, work_charge, &self.caps)?;
        let total_retained_bytes = admission::checked_retained_total(
            self.total_retained_bytes,
            retained_charge,
            &self.caps,
        )?;
        Ok((total_retained_bytes, total_work))
    }

    fn commit_expr(
        &mut self,
        expression: Expr,
        fingerprint: ExprFingerprint,
        derived: (Shape, Dims, Class),
        depth: u32,
        total_retained_bytes: u64,
        total_work: u64,
    ) -> NodeId {
        let (shape, dims, class) = derived;
        self.exprs.push(expression);
        self.shapes.push(shape);
        self.dims.push(dims);
        self.classes.push(class);
        self.depths.push(depth);
        self.total_retained_bytes = total_retained_bytes;
        self.total_work = total_work;
        let id = NodeId((self.exprs.len() - 1) as u32);
        match self.intern.entry(fingerprint) {
            Entry::Vacant(slot) => {
                slot.insert(InternEntry::One(id));
            }
            Entry::Occupied(mut slot) => {
                let entry = slot.get_mut();
                match entry {
                    InternEntry::One(previous) => {
                        let previous = *previous;
                        *entry = InternEntry::Collision(vec![previous, id]);
                    }
                    InternEntry::Collision(ids) => ids.push(id),
                }
            }
        }
        id
    }

    fn interned_id(
        &self,
        fingerprint: &ExprFingerprint,
        mut matches: impl FnMut(&Expr) -> bool,
    ) -> Option<NodeId> {
        let entry = self.intern.get(fingerprint)?;
        let mut check = |id: NodeId| {
            let expression = self.exprs.get(id.0 as usize)?;
            matches(expression).then_some(id)
        };
        match entry {
            InternEntry::One(id) => check(*id),
            InternEntry::Collision(ids) => ids.iter().find_map(|id| check(*id)),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn push_borrowed_string_expr(
        &mut self,
        value: &str,
        derived: (Shape, Dims, Class),
        depth: u32,
        work_charge: u64,
        fingerprint: ExprFingerprint,
        mut matches: impl FnMut(&Expr) -> bool,
        make: impl FnOnce(String) -> Expr,
    ) -> Result<NodeId, OptError> {
        if let Some(id) = self.interned_id(&fingerprint, &mut matches) {
            return Ok(id);
        }
        let (total_retained_bytes, total_work) =
            self.new_expr_totals(admission::string_expr_retained_bytes(value), work_charge)?;
        // All leaf, depth, count, work, and retained-byte checks have
        // passed. This is the first point at which the external string
        // becomes owned.
        let expression = make(value.to_owned());
        Ok(self.commit_expr(
            expression,
            fingerprint,
            derived,
            depth,
            total_retained_bytes,
            total_work,
        ))
    }

    /// Validate a candidate expression through the shared admission
    /// rules, then intern it. Interned hits return the EXISTING id
    /// (identical expressions were valid when first admitted); every
    /// rejection happens before any mutation.
    fn push_checked(&mut self, e: Expr) -> Result<NodeId, OptError> {
        let derived = {
            let lookup = |n: NodeId| {
                let i = n.0 as usize;
                (i < self.exprs.len()).then(|| (self.shapes[i], self.dims[i], self.classes[i]))
            };
            admission::derive_expr(&e, &lookup, &self.vars, &self.caps)?
        };
        let depth = {
            let lookup = |n: NodeId| self.depths.get(n.0 as usize).copied();
            admission::derive_depth(&e, &lookup, &self.caps)?
        };
        let fingerprint = ExprFingerprint::of(&e);
        if let Some(id) = self.interned_id(&fingerprint, |existing| expr_identity_eq(existing, &e))
        {
            return Ok(id);
        }
        let (total_retained_bytes, total_work) =
            self.new_expr_totals(admission::expr_retained_bytes(&e), admission::expr_work(&e))?;
        Ok(self.commit_expr(
            e,
            fingerprint,
            derived,
            depth,
            total_retained_bytes,
            total_work,
        ))
    }

    /// Variable reference node.
    ///
    /// # Errors
    /// [`OptError::UnknownVar`].
    pub fn var_ref(&mut self, v: VarId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Var(v))
    }

    /// Component extraction.
    ///
    /// # Errors
    /// Shape/index teaching errors.
    pub fn component(&mut self, of: NodeId, index: u32) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Component { of, index })
    }

    /// Dimensioned constant. The value must be FINITE — non-finite
    /// constants cannot acquire graph authority.
    ///
    /// # Errors
    /// [`OptError::NonFinite`].
    pub fn konst(&mut self, value: f64, dims: Dims) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Const { value, dims })
    }

    /// Sum.
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn add(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Add(a, b))
    }

    /// Difference.
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn sub(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Sub(a, b))
    }

    /// Pointwise minimum (C0 — poisons smooth routing, on purpose).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn min_of(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Min(a, b))
    }

    /// Pointwise maximum (C0).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn max_of(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Max(a, b))
    }

    /// Product (scalar×scalar or scalar×vector; dimensions add,
    /// CHECKED against the i8 exponent domain).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn mul(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Mul(a, b))
    }

    /// Quotient (scalars; dimensions subtract, CHECKED).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn div(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Div(a, b))
    }

    /// Negation.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn neg(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Neg(a))
    }

    /// Integer power (scalar).
    ///
    /// # Errors
    /// Shape or dimension-overflow teaching errors.
    pub fn powi(&mut self, base: NodeId, exp: i32) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Powi { base, exp })
    }

    /// Square root (even dimension exponents halve).
    ///
    /// # Errors
    /// Shape/[`OptError::OddDims`] teaching errors.
    pub fn sqrt(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Sqrt(a))
    }

    /// Exponential (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn exp(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Exp(a))
    }

    /// Natural log (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn ln(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Ln(a))
    }

    /// Hyperbolic tangent (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn tanh(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Tanh(a))
    }

    /// Inner product of same-length vectors.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn dot(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Dot(a, b))
    }

    /// Squared norm of a vector.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn norm_sq(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::NormSq(a))
    }

    /// Absolute value (C0).
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn abs(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.push_checked(Expr::Abs(a))
    }

    /// PDE residual node (FLUX study reference + adjoint metadata).
    ///
    /// # Errors
    /// [`OptError::UnknownVar`] / [`OptError::CapExceeded`].
    pub fn pde_residual(
        &mut self,
        study: &str,
        over: VarId,
        adjoint_available: bool,
        dims: Dims,
    ) -> Result<NodeId, OptError> {
        let derived = admission::derive_pde_leaf(study, over, dims, &self.vars, &self.caps)?;
        self.push_borrowed_string_expr(
            study,
            derived,
            1,
            1,
            ExprFingerprint::pde(study, over, adjoint_available, dims),
            |existing| {
                matches!(
                    existing,
                    Expr::PdeResidual {
                        study: existing_study,
                        over: existing_over,
                        adjoint_available: existing_adjoint,
                        dims: existing_dims,
                    } if existing_study == study
                        && *existing_over == over
                        && *existing_adjoint == adjoint_available
                        && *existing_dims == dims
                )
            },
            |owned| Expr::PdeResidual {
                study: owned,
                over,
                adjoint_available,
                dims,
            },
        )
    }

    /// Expectation over a UQ configuration.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn expectation(&mut self, of: NodeId, uq_config: &str) -> Result<NodeId, OptError> {
        let lookup = |node: NodeId| {
            let index = node.0 as usize;
            (index < self.exprs.len())
                .then(|| (self.shapes[index], self.dims[index], self.classes[index]))
        };
        let derived = admission::derive_expectation_leaf(of, uq_config, &lookup, &self.caps)?;
        let depth = self.unary_depth(of)?;
        self.push_borrowed_string_expr(
            uq_config,
            derived,
            depth,
            2,
            ExprFingerprint::expectation(of, uq_config),
            |existing| {
                matches!(
                    existing,
                    Expr::Expectation {
                        of: existing_of,
                        uq_config: existing_config,
                    } if *existing_of == of && existing_config == uq_config
                )
            },
            |owned| Expr::Expectation {
                of,
                uq_config: owned,
            },
        )
    }

    /// CVaR at tail level `alpha` (C0).
    ///
    /// # Errors
    /// Shape/[`OptError::BadParam`] teaching errors.
    pub fn cvar(&mut self, of: NodeId, alpha: f64, uq_config: &str) -> Result<NodeId, OptError> {
        let lookup = |node: NodeId| {
            let index = node.0 as usize;
            (index < self.exprs.len())
                .then(|| (self.shapes[index], self.dims[index], self.classes[index]))
        };
        let derived = admission::derive_cvar_leaf(of, alpha, uq_config, &lookup, &self.caps)?;
        let depth = self.unary_depth(of)?;
        self.push_borrowed_string_expr(
            uq_config,
            derived,
            depth,
            2,
            ExprFingerprint::cvar(of, alpha, uq_config),
            |existing| {
                matches!(
                    existing,
                    Expr::Cvar {
                        of: existing_of,
                        alpha: existing_alpha,
                        uq_config: existing_config,
                    } if *existing_of == of
                        && existing_alpha.to_bits() == alpha.to_bits()
                        && existing_config == uq_config
                )
            },
            |owned| Expr::Cvar {
                of,
                alpha,
                uq_config: owned,
            },
        )
    }

    /// Quantile at level `q` (C0).
    ///
    /// # Errors
    /// Shape/[`OptError::BadParam`] teaching errors.
    pub fn quantile(&mut self, of: NodeId, q: f64, uq_config: &str) -> Result<NodeId, OptError> {
        let lookup = |node: NodeId| {
            let index = node.0 as usize;
            (index < self.exprs.len())
                .then(|| (self.shapes[index], self.dims[index], self.classes[index]))
        };
        let derived = admission::derive_quantile_leaf(of, q, uq_config, &lookup, &self.caps)?;
        let depth = self.unary_depth(of)?;
        self.push_borrowed_string_expr(
            uq_config,
            derived,
            depth,
            2,
            ExprFingerprint::quantile(of, q, uq_config),
            |existing| {
                matches!(
                    existing,
                    Expr::Quantile {
                        of: existing_of,
                        q: existing_q,
                        uq_config: existing_config,
                    } if *existing_of == of
                        && existing_q.to_bits() == q.to_bits()
                        && existing_config == uq_config
                )
            },
            |owned| Expr::Quantile {
                of,
                q,
                uq_config: owned,
            },
        )
    }

    /// Declare an objective (scalar root; weight FINITE and
    /// NONNEGATIVE — `Sense` already carries direction, so signed
    /// weights are refused rather than silently flipping it; `-0.0` is
    /// refused because bit-pattern serialization would give two wire
    /// identities to one meaning).
    ///
    /// # Errors
    /// [`OptError::NotScalar`] / [`OptError::BadParam`] /
    /// [`OptError::CapExceeded`].
    pub fn objective(&mut self, node: NodeId, sense: Sense, weight: f64) -> Result<(), OptError> {
        self.require_scalar_root(node)?;
        admission::validate_weight(weight)?;
        if self.objectives.len() as u64 >= u64::from(self.caps.max_objectives) {
            return Err(OptError::CapExceeded {
                what: "objectives",
                count: self.objectives.len() as u64 + 1,
                cap: u64::from(self.caps.max_objectives),
            });
        }
        let total_work = admission::checked_work_total(self.total_work, 1, &self.caps)?;
        self.objectives.push(Objective {
            node,
            sense,
            weight,
        });
        self.total_work = total_work;
        Ok(())
    }

    /// Declare a constraint (scalar root; semantics live in
    /// fs-constraint).
    ///
    /// # Errors
    /// [`OptError::NotScalar`] / [`OptError::CapExceeded`].
    pub fn constraint(
        &mut self,
        node: NodeId,
        kind: ConstraintKind,
        name: &str,
    ) -> Result<(), OptError> {
        self.require_scalar_root(node)?;
        admission::validate_name("constraint name", name, &self.caps)?;
        if self.constraints.len() as u64 >= u64::from(self.caps.max_constraints) {
            return Err(OptError::CapExceeded {
                what: "constraints",
                count: self.constraints.len() as u64 + 1,
                cap: u64::from(self.caps.max_constraints),
            });
        }
        let total_retained_bytes = admission::checked_retained_total(
            self.total_retained_bytes,
            admission::name_retained_bytes(name),
            &self.caps,
        )?;
        let total_work = admission::checked_work_total(self.total_work, 1, &self.caps)?;
        self.constraints.push(Constraint {
            node,
            kind,
            name: name.to_string(),
        });
        self.total_retained_bytes = total_retained_bytes;
        self.total_work = total_work;
        Ok(())
    }

    /// Finish. The graph is valid by construction: fields are sealed
    /// and every mutating path above validated through the same
    /// versioned rules [`Problem::admit`] re-checks, so builder output
    /// always admits (the conformance suite pins that agreement).
    #[must_use]
    pub fn finish(self) -> Problem {
        Problem {
            vars: self.vars,
            exprs: self.exprs,
            objectives: self.objectives,
            constraints: self.constraints,
            tags: self.tags,
            budget: self.budget,
            shapes: self.shapes,
            dims: self.dims,
            classes: self.classes,
            depths: self.depths,
            total_work: self.total_work,
        }
    }
}

#[cfg(test)]
mod cap_allocation_tests {
    use super::{Dims, Manifold, ProblemBuilder, ProblemTag};

    #[derive(Debug, PartialEq, Eq)]
    struct StorageSnapshot {
        vars: (usize, usize),
        exprs: (usize, usize),
        shapes: (usize, usize),
        dims: (usize, usize),
        classes: (usize, usize),
        depths: (usize, usize),
        interned: usize,
        objectives: (usize, usize),
        constraints: (usize, usize),
        tags: (usize, usize),
        owned_string_capacity: usize,
        retained_bytes: u64,
        work: u64,
    }

    fn snapshot(builder: &ProblemBuilder) -> StorageSnapshot {
        let owned_string_capacity = builder
            .vars
            .iter()
            .map(|variable| variable.name.capacity())
            .chain(builder.constraints.iter().map(|item| item.name.capacity()))
            .chain(
                builder
                    .exprs
                    .iter()
                    .filter_map(|expression| match expression {
                        super::Expr::PdeResidual { study, .. } => Some(study.capacity()),
                        super::Expr::Expectation { uq_config, .. }
                        | super::Expr::Cvar { uq_config, .. }
                        | super::Expr::Quantile { uq_config, .. } => Some(uq_config.capacity()),
                        _ => None,
                    }),
            )
            .sum();
        StorageSnapshot {
            vars: (builder.vars.len(), builder.vars.capacity()),
            exprs: (builder.exprs.len(), builder.exprs.capacity()),
            shapes: (builder.shapes.len(), builder.shapes.capacity()),
            dims: (builder.dims.len(), builder.dims.capacity()),
            classes: (builder.classes.len(), builder.classes.capacity()),
            depths: (builder.depths.len(), builder.depths.capacity()),
            interned: builder.intern.len(),
            objectives: (builder.objectives.len(), builder.objectives.capacity()),
            constraints: (builder.constraints.len(), builder.constraints.capacity()),
            tags: (builder.tags.len(), builder.tags.capacity()),
            owned_string_capacity,
            retained_bytes: builder.total_retained_bytes,
            work: builder.total_work,
        }
    }

    #[test]
    fn cap_plus_one_paths_do_not_grow_builder_storage() {
        let mut node_caps = super::AdmissionCaps::default();
        node_caps.max_nodes = 1;
        let mut nodes = ProblemBuilder::with_caps(node_caps);
        let root = nodes.konst(1.0, Dims::NONE).expect("sole node");
        let before = snapshot(&nodes);
        assert!(nodes.neg(root).is_err());
        assert_eq!(snapshot(&nodes), before, "push_checked cap+1 is atomic");

        let mut uq = ProblemBuilder::new();
        let scalar = uq.konst(1.0, Dims::NONE).expect("scalar");
        uq.caps.max_total_retained_bytes = uq
            .total_retained_bytes
            .saturating_add(super::admission::string_expr_retained_bytes("uq01"))
            - 1;
        let before = snapshot(&uq);
        assert!(uq.expectation(scalar, "uq01").is_err());
        assert_eq!(snapshot(&uq), before, "UQ cap+1 does not clone input");

        let mut pde = ProblemBuilder::new();
        let variable = pde
            .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("variable");
        pde.caps.max_total_retained_bytes = pde
            .total_retained_bytes
            .saturating_add(super::admission::string_expr_retained_bytes("flow"))
            - 1;
        let before = snapshot(&pde);
        assert!(
            pde.pde_residual("flow", variable, false, Dims::NONE)
                .is_err()
        );
        assert_eq!(snapshot(&pde), before, "PDE cap+1 does not clone input");

        let mut work = ProblemBuilder::new();
        work.caps.max_total_work = work.total_work;
        let before = snapshot(&work);
        assert!(work.tag(ProblemTag::MultiFidelity { levels: 1 }).is_err());
        assert_eq!(snapshot(&work), before, "work cap+1 is atomic");
    }

    #[test]
    fn fixed_fingerprints_preserve_exact_identity_and_resolve_collisions() {
        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("variable");
        let unrelated = builder.konst(1.0, Dims::NONE).expect("constant");

        // Inject the theoretical full-width-hash collision state: this
        // bucket points at an expression that is not the PDE candidate.
        // Exact expression comparison must decline it, retain the PDE,
        // and find that PDE (not the colliding entry) on the next call.
        let fingerprint = super::ExprFingerprint::pde("flow", variable, false, Dims::NONE);
        let replaced = builder
            .intern
            .insert(fingerprint, super::InternEntry::One(unrelated));
        assert!(replaced.is_none(), "synthetic collision key is fresh");
        let pde = builder
            .pde_residual("flow", variable, false, Dims::NONE)
            .expect("colliding fingerprint still admits distinct expression");
        assert_ne!(pde, unrelated);
        assert_eq!(
            builder
                .pde_residual("flow", variable, false, Dims::NONE)
                .expect("exact PDE duplicate"),
            pde
        );
        assert!(matches!(
            builder.intern.get(&fingerprint),
            Some(super::InternEntry::Collision(ids))
                if ids.len() == 2 && ids[0] == unrelated && ids[1] == pde
        ));
        assert_ne!(
            builder
                .pde_residual("flow-2", variable, false, Dims::NONE)
                .expect("different external identity"),
            pde
        );

        let expectation = builder.expectation(unrelated, "uq").expect("expectation");
        assert_eq!(
            builder
                .expectation(unrelated, "uq")
                .expect("exact expectation duplicate"),
            expectation
        );
        assert_ne!(
            builder
                .expectation(unrelated, "uq-2")
                .expect("different UQ identity"),
            expectation
        );

        let positive_zero = builder.konst(0.0, Dims::NONE).expect("positive zero");
        let negative_zero = builder.konst(-0.0, Dims::NONE).expect("negative zero");
        assert_ne!(
            positive_zero, negative_zero,
            "float payload identity is bitwise"
        );
    }
}

//! The typed problem graph: expression nodes over manifold variables,
//! hash-consed (common subexpressions are IDENTICAL node ids), with
//! per-step validation — shape rules, dimension rules (fs-qty `Dims`),
//! and differentiability-CLASS propagation, so "this objective is
//! non-smooth through that min()" is knowable at BUILD time and
//! optimizer routing can refuse with the offending node named.

use fs_qty::Dims;
use std::collections::BTreeMap;

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
    /// Storage length of one point.
    #[must_use]
    pub fn point_dim(&self) -> u32 {
        match *self {
            Manifold::Rn { dim } => dim,
            Manifold::Sphere { ambient } => ambient,
            Manifold::So3 => 4,
            Manifold::Stiefel { n, p } => n * p,
        }
    }

    /// Tangent-space dimension (what a Riemannian gradient has).
    #[must_use]
    pub fn tangent_dim(&self) -> u32 {
        match *self {
            Manifold::Rn { dim } => dim,
            Manifold::Sphere { ambient } => ambient - 1,
            Manifold::So3 => 3,
            Manifold::Stiefel { n, p } => n * p - p * (p + 1) / 2,
        }
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

/// Problem-structure annotations (representable, not yet executed).
#[derive(Debug, Clone, PartialEq)]
pub enum ProblemTag {
    /// Multiple fidelity levels are available.
    MultiFidelity {
        /// Number of levels.
        levels: u32,
    },
    /// Chance constraint: `P(g ≤ 0) ≥ prob`.
    ChanceConstrained {
        /// Required probability.
        prob: f64,
    },
    /// Bilevel structure (inner problem referenced by hash).
    Bilevel {
        /// Inner problem hash.
        inner_hash: u64,
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
        left: [i8; 5],
        /// Right dims.
        right: [i8; 5],
    },
    /// Dimension exponents would leave the representable runtime domain.
    DimOverflow {
        /// Operation that attempted the dimension arithmetic.
        op: &'static str,
        /// Base dimensions before scaling.
        dims: [i8; 5],
        /// Requested integer exponent.
        exponent: i32,
    },
    /// A transcendental was fed a dimensioned quantity.
    NonDimensionless {
        /// Operation name.
        op: &'static str,
        /// The offending dims.
        dims: [i8; 5],
    },
    /// `sqrt` of odd dimension exponents.
    OddDims {
        /// The offending dims.
        dims: [i8; 5],
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
            OptError::Cancelled => write!(f, "cancelled between descent steps"),
            OptError::BudgetExhausted { spent } => write!(
                f,
                "evaluation budget exhausted after {spent} evaluations (P4 receipt)"
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

/// The finished problem: graph + roots + metadata. Construct through
/// [`ProblemBuilder`]; every accessor is cheap.
#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    /// Declared variables.
    pub vars: Vec<Variable>,
    /// Expression nodes (hash-consed; ids index this list).
    pub exprs: Vec<Expr>,
    /// Objectives (multi-objective representable).
    pub objectives: Vec<Objective>,
    /// Constraints.
    pub constraints: Vec<Constraint>,
    /// Structure annotations.
    pub tags: Vec<ProblemTag>,
    /// Evaluation budget (P4).
    pub budget: EvalBudget,
    pub(crate) shapes: Vec<Shape>,
    pub(crate) dims: Vec<Dims>,
    pub(crate) classes: Vec<Class>,
}

impl Problem {
    /// Shape of a node.
    #[must_use]
    pub fn shape(&self, n: NodeId) -> Shape {
        self.shapes[n.0 as usize]
    }

    /// Dimensions of a node.
    #[must_use]
    pub fn node_dims(&self, n: NodeId) -> Dims {
        self.dims[n.0 as usize]
    }

    /// Propagated differentiability class of a node.
    #[must_use]
    pub fn class(&self, n: NodeId) -> Class {
        self.classes[n.0 as usize]
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
            for n in self.reachable(root) {
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

    /// Nodes reachable from `root` (build order).
    #[must_use]
    pub fn reachable(&self, root: NodeId) -> Vec<NodeId> {
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
        out
    }
}

/// Child node ids of an expression.
pub(crate) fn children(e: &Expr) -> Vec<NodeId> {
    match *e {
        Expr::Var(_) | Expr::Const { .. } | Expr::PdeResidual { .. } => vec![],
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
        | Expr::Quantile { of, .. } => vec![of],
        Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Dot(a, b)
        | Expr::Min(a, b)
        | Expr::Max(a, b) => vec![a, b],
    }
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
/// G0 common-subexpression identity).
#[derive(Debug, Default)]
pub struct ProblemBuilder {
    vars: Vec<Variable>,
    exprs: Vec<Expr>,
    shapes: Vec<Shape>,
    dims: Vec<Dims>,
    classes: Vec<Class>,
    intern: BTreeMap<String, NodeId>,
    objectives: Vec<Objective>,
    constraints: Vec<Constraint>,
    tags: Vec<ProblemTag>,
    budget: EvalBudget,
}

impl ProblemBuilder {
    /// Empty builder (unlimited budget).
    #[must_use]
    pub fn new() -> Self {
        ProblemBuilder {
            budget: EvalBudget { max_evals: 0 },
            ..ProblemBuilder::default()
        }
    }

    /// Declare a variable.
    pub fn var(&mut self, name: &str, manifold: Manifold, dims: Dims) -> VarId {
        self.vars.push(Variable {
            name: name.to_string(),
            manifold,
            dims,
        });
        VarId((self.vars.len() - 1) as u32)
    }

    /// Attach the evaluation budget (P4).
    pub fn set_budget(&mut self, max_evals: u64) {
        self.budget.max_evals = max_evals;
    }

    /// Attach a structure tag.
    pub fn tag(&mut self, tag: ProblemTag) {
        self.tags.push(tag);
    }

    fn check_node(&self, n: NodeId) -> Result<(), OptError> {
        if (n.0 as usize) < self.exprs.len() {
            Ok(())
        } else {
            Err(OptError::UnknownNode { id: n.0 })
        }
    }

    fn scalar(&self, op: &'static str, n: NodeId) -> Result<(), OptError> {
        match self.shapes[n.0 as usize] {
            Shape::Scalar => Ok(()),
            v @ Shape::Vector(_) => Err(OptError::ShapeMismatch {
                op,
                left: v,
                right: Shape::Scalar,
            }),
        }
    }

    fn push(&mut self, e: Expr, shape: Shape, dims: Dims) -> NodeId {
        let key = crate::serial::expr_key(&e);
        if let Some(&id) = self.intern.get(&key) {
            return id;
        }
        let class = children(&e)
            .iter()
            .map(|c| self.classes[c.0 as usize])
            .chain([own_class(&e)])
            .min()
            .unwrap_or(Class::Smooth);
        self.exprs.push(e);
        self.shapes.push(shape);
        self.dims.push(dims);
        self.classes.push(class);
        let id = NodeId((self.exprs.len() - 1) as u32);
        self.intern.insert(key, id);
        id
    }

    /// Variable reference node.
    ///
    /// # Errors
    /// [`OptError::UnknownVar`].
    pub fn var_ref(&mut self, v: VarId) -> Result<NodeId, OptError> {
        let var = self
            .vars
            .get(v.0 as usize)
            .ok_or(OptError::UnknownVar { id: v.0 })?;
        let (shape, dims) = (Shape::Vector(var.manifold.point_dim()), var.dims);
        Ok(self.push(Expr::Var(v), shape, dims))
    }

    /// Component extraction.
    ///
    /// # Errors
    /// Shape/index teaching errors.
    pub fn component(&mut self, of: NodeId, index: u32) -> Result<NodeId, OptError> {
        self.check_node(of)?;
        match self.shapes[of.0 as usize] {
            Shape::Vector(len) if index < len => {
                let dims = self.dims[of.0 as usize];
                Ok(self.push(Expr::Component { of, index }, Shape::Scalar, dims))
            }
            Shape::Vector(len) => Err(OptError::IndexOut { index, len }),
            Shape::Scalar => Err(OptError::ShapeMismatch {
                op: "component",
                left: Shape::Scalar,
                right: Shape::Vector(index + 1),
            }),
        }
    }

    /// Dimensioned constant.
    pub fn konst(&mut self, value: f64, dims: Dims) -> NodeId {
        self.push(Expr::Const { value, dims }, Shape::Scalar, dims)
    }

    fn binary_same(
        &mut self,
        op: &'static str,
        make: fn(NodeId, NodeId) -> Expr,
        a: NodeId,
        b: NodeId,
    ) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.check_node(b)?;
        let (sa, sb) = (self.shapes[a.0 as usize], self.shapes[b.0 as usize]);
        if sa != sb {
            return Err(OptError::ShapeMismatch {
                op,
                left: sa,
                right: sb,
            });
        }
        let (da, db) = (self.dims[a.0 as usize], self.dims[b.0 as usize]);
        if da != db {
            return Err(OptError::DimMismatch {
                op,
                left: da.0,
                right: db.0,
            });
        }
        Ok(self.push(make(a, b), sa, da))
    }

    /// Sum.
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn add(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.binary_same("add", Expr::Add, a, b)
    }

    /// Difference.
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn sub(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.binary_same("sub", Expr::Sub, a, b)
    }

    /// Pointwise minimum (C0 — poisons smooth routing, on purpose).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn min_of(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.binary_same("min", Expr::Min, a, b)
    }

    /// Pointwise maximum (C0).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn max_of(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.binary_same("max", Expr::Max, a, b)
    }

    /// Product (scalar×scalar or scalar×vector).
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn mul(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.check_node(b)?;
        let (sa, sb) = (self.shapes[a.0 as usize], self.shapes[b.0 as usize]);
        let dims = self.dims[a.0 as usize].plus(self.dims[b.0 as usize]);
        let shape = match (sa, sb) {
            (Shape::Scalar, Shape::Scalar) => Shape::Scalar,
            (Shape::Scalar, Shape::Vector(n)) | (Shape::Vector(n), Shape::Scalar) => {
                Shape::Vector(n)
            }
            (l, r) => {
                return Err(OptError::ShapeMismatch {
                    op: "mul",
                    left: l,
                    right: r,
                });
            }
        };
        Ok(self.push(Expr::Mul(a, b), shape, dims))
    }

    /// Quotient (scalars).
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn div(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.check_node(b)?;
        self.scalar("div", a)?;
        self.scalar("div", b)?;
        let dims = Dims([
            self.dims[a.0 as usize].0[0].saturating_sub(self.dims[b.0 as usize].0[0]),
            self.dims[a.0 as usize].0[1].saturating_sub(self.dims[b.0 as usize].0[1]),
            self.dims[a.0 as usize].0[2].saturating_sub(self.dims[b.0 as usize].0[2]),
            self.dims[a.0 as usize].0[3].saturating_sub(self.dims[b.0 as usize].0[3]),
            self.dims[a.0 as usize].0[4].saturating_sub(self.dims[b.0 as usize].0[4]),
        ]);
        Ok(self.push(Expr::Div(a, b), Shape::Scalar, dims))
    }

    /// Negation.
    ///
    /// # Errors
    /// [`OptError::UnknownNode`].
    pub fn neg(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        let (s, d) = (self.shapes[a.0 as usize], self.dims[a.0 as usize]);
        Ok(self.push(Expr::Neg(a), s, d))
    }

    /// Integer power (scalar).
    ///
    /// # Errors
    /// Shape or dimension-overflow teaching errors.
    pub fn powi(&mut self, base: NodeId, exp: i32) -> Result<NodeId, OptError> {
        self.check_node(base)?;
        self.scalar("powi", base)?;
        let d = self.dims[base.0 as usize].0;
        let mut scaled = [0i8; 5];
        for (out, &e) in scaled.iter_mut().zip(&d) {
            let product = i32::from(e).checked_mul(exp).ok_or(OptError::DimOverflow {
                op: "powi",
                dims: d,
                exponent: exp,
            })?;
            *out = i8::try_from(product).map_err(|_| OptError::DimOverflow {
                op: "powi",
                dims: d,
                exponent: exp,
            })?;
        }
        let dims = Dims(scaled);
        Ok(self.push(Expr::Powi { base, exp }, Shape::Scalar, dims))
    }

    /// Square root (even dimension exponents halve).
    ///
    /// # Errors
    /// Shape/[`OptError::OddDims`] teaching errors.
    pub fn sqrt(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.scalar("sqrt", a)?;
        let d = self.dims[a.0 as usize].0;
        if d.iter().any(|e| e % 2 != 0) {
            return Err(OptError::OddDims { dims: d });
        }
        Ok(self.push(Expr::Sqrt(a), Shape::Scalar, Dims(d.map(|e| e / 2))))
    }

    fn transcendental(
        &mut self,
        op: &'static str,
        make: fn(NodeId) -> Expr,
        a: NodeId,
    ) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.scalar(op, a)?;
        let d = self.dims[a.0 as usize];
        if d != Dims::NONE {
            return Err(OptError::NonDimensionless { op, dims: d.0 });
        }
        Ok(self.push(make(a), Shape::Scalar, Dims::NONE))
    }

    /// Exponential (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn exp(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.transcendental("exp", Expr::Exp, a)
    }

    /// Natural log (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn ln(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.transcendental("ln", Expr::Ln, a)
    }

    /// Hyperbolic tangent (dimensionless).
    ///
    /// # Errors
    /// Shape/dimension teaching errors.
    pub fn tanh(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.transcendental("tanh", Expr::Tanh, a)
    }

    /// Inner product of same-length vectors.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn dot(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.check_node(b)?;
        match (self.shapes[a.0 as usize], self.shapes[b.0 as usize]) {
            (Shape::Vector(n), Shape::Vector(m)) if n == m => {
                let dims = self.dims[a.0 as usize].plus(self.dims[b.0 as usize]);
                Ok(self.push(Expr::Dot(a, b), Shape::Scalar, dims))
            }
            (l, r) => Err(OptError::ShapeMismatch {
                op: "dot",
                left: l,
                right: r,
            }),
        }
    }

    /// Squared norm of a vector.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn norm_sq(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        match self.shapes[a.0 as usize] {
            Shape::Vector(_) => {
                let d = self.dims[a.0 as usize];
                Ok(self.push(Expr::NormSq(a), Shape::Scalar, d.plus(d)))
            }
            s @ Shape::Scalar => Err(OptError::ShapeMismatch {
                op: "norm_sq",
                left: s,
                right: Shape::Vector(1),
            }),
        }
    }

    /// Absolute value (C0).
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn abs(&mut self, a: NodeId) -> Result<NodeId, OptError> {
        self.check_node(a)?;
        self.scalar("abs", a)?;
        let d = self.dims[a.0 as usize];
        Ok(self.push(Expr::Abs(a), Shape::Scalar, d))
    }

    /// PDE residual node (FLUX study reference + adjoint metadata).
    ///
    /// # Errors
    /// [`OptError::UnknownVar`].
    pub fn pde_residual(
        &mut self,
        study: &str,
        over: VarId,
        adjoint_available: bool,
        dims: Dims,
    ) -> Result<NodeId, OptError> {
        if (over.0 as usize) >= self.vars.len() {
            return Err(OptError::UnknownVar { id: over.0 });
        }
        Ok(self.push(
            Expr::PdeResidual {
                study: study.to_string(),
                over,
                adjoint_available,
                dims,
            },
            Shape::Scalar,
            dims,
        ))
    }

    /// Expectation over a UQ configuration.
    ///
    /// # Errors
    /// Shape teaching errors.
    pub fn expectation(&mut self, of: NodeId, uq_config: &str) -> Result<NodeId, OptError> {
        self.check_node(of)?;
        self.scalar("expectation", of)?;
        let d = self.dims[of.0 as usize];
        Ok(self.push(
            Expr::Expectation {
                of,
                uq_config: uq_config.to_string(),
            },
            Shape::Scalar,
            d,
        ))
    }

    /// CVaR at tail level `alpha` (C0).
    ///
    /// # Errors
    /// Shape/[`OptError::BadParam`] teaching errors.
    pub fn cvar(&mut self, of: NodeId, alpha: f64, uq_config: &str) -> Result<NodeId, OptError> {
        self.check_node(of)?;
        self.scalar("cvar", of)?;
        if !(alpha > 0.0 && alpha < 1.0) {
            return Err(OptError::BadParam {
                what: "cvar alpha",
                value: alpha,
            });
        }
        let d = self.dims[of.0 as usize];
        Ok(self.push(
            Expr::Cvar {
                of,
                alpha,
                uq_config: uq_config.to_string(),
            },
            Shape::Scalar,
            d,
        ))
    }

    /// Quantile at level `q` (C0).
    ///
    /// # Errors
    /// Shape/[`OptError::BadParam`] teaching errors.
    pub fn quantile(&mut self, of: NodeId, q: f64, uq_config: &str) -> Result<NodeId, OptError> {
        self.check_node(of)?;
        self.scalar("quantile", of)?;
        if !(q > 0.0 && q < 1.0) {
            return Err(OptError::BadParam {
                what: "quantile q",
                value: q,
            });
        }
        let d = self.dims[of.0 as usize];
        Ok(self.push(
            Expr::Quantile {
                of,
                q,
                uq_config: uq_config.to_string(),
            },
            Shape::Scalar,
            d,
        ))
    }

    /// Declare an objective (scalar root).
    ///
    /// # Errors
    /// [`OptError::NotScalar`].
    pub fn objective(&mut self, node: NodeId, sense: Sense, weight: f64) -> Result<(), OptError> {
        self.check_node(node)?;
        if self.shapes[node.0 as usize] != Shape::Scalar {
            return Err(OptError::NotScalar { node: node.0 });
        }
        self.objectives.push(Objective {
            node,
            sense,
            weight,
        });
        Ok(())
    }

    /// Declare a constraint (scalar root; semantics live in
    /// fs-constraint).
    ///
    /// # Errors
    /// [`OptError::NotScalar`].
    pub fn constraint(
        &mut self,
        node: NodeId,
        kind: ConstraintKind,
        name: &str,
    ) -> Result<(), OptError> {
        self.check_node(node)?;
        if self.shapes[node.0 as usize] != Shape::Scalar {
            return Err(OptError::NotScalar { node: node.0 });
        }
        self.constraints.push(Constraint {
            node,
            kind,
            name: name.to_string(),
        });
        Ok(())
    }

    /// Finish (the graph is valid by construction).
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
        }
    }
}

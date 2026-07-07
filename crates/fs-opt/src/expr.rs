//! The typed expression graph (plan §9.1): objectives and constraints as
//! nodes over design variables, with TWO propagated static properties —
//! physical dimensions (fs-qty, so a pressure never adds to a stress) and
//! the DIFFERENTIABILITY CLASS ("this objective is non-smooth through that
//! min()" is known at build time and routed to the right optimizer
//! family, never discovered at iteration 400).
//!
//! The graph is HASH-CONSED: structurally identical nodes get the same id
//! (the common-subexpression identity law, G0-tested), which also makes
//! problem hashes meaningful.
//!
//! Determinism: arena order is insertion order; hash-consing keys are
//! structural; serialization is canonical — identical builds give
//! identical graphs, ids, and hashes (P2).

use std::collections::BTreeMap;

use crate::{OptDiag, VarId};
use fs_qty::Dims;

/// Node id in the expression arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ExprId(pub u32);

/// How differentiable a value is along every path that produced it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Smoothness {
    /// C¹ (or better) with available derivatives.
    Smooth,
    /// Continuous with kinks (min/max/abs/CVaR): subgradients exist.
    Kinked,
    /// No derivative claim (opaque physics without adjoints, black boxes).
    Opaque,
}

impl Smoothness {
    /// Combining values: the WEAKEST class wins (conservative).
    #[must_use]
    pub fn join(self, other: Smoothness) -> Smoothness {
        self.max(other)
    }

    /// Stable lowercase name.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Smoothness::Smooth => "smooth",
            Smoothness::Kinked => "kinked",
            Smoothness::Opaque => "opaque",
        }
    }
}

/// One expression node. Children are ids INTO the same arena.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ExprNode {
    /// A scalar design variable component (reduction of a variable is a
    /// higher-level concern; the graph is scalar-typed).
    Var(VarId),
    /// A dimensioned constant (value bits + dims — canonical).
    Const {
        /// Value bits (total order for hash-consing).
        bits: u64,
        /// Physical dimensions.
        dims: Dims,
    },
    /// Σ of ≥ 2 addends (same dims).
    Add(Vec<ExprId>),
    /// Product of two factors (dims add).
    Mul(ExprId, ExprId),
    /// Negation.
    Neg(ExprId),
    /// Integer power (dims scale).
    Powi(ExprId, i32),
    /// Pointwise minimum (KINK).
    Min(Vec<ExprId>),
    /// Pointwise maximum (KINK).
    Max(Vec<ExprId>),
    /// Absolute value (KINK at 0).
    Abs(ExprId),
    /// A PDE residual/QoI from a FLUX study: opaque physics with declared
    /// adjoint availability (plan: adjoint metadata is first-class).
    PdeQoi {
        /// Study identity (ledger reference).
        study: String,
        /// Output dims.
        dims: Dims,
        /// Whether a discrete adjoint is available (smooth if so).
        adjoint: bool,
    },
    /// Expectation over a UQ configuration (class-preserving).
    Expect {
        /// The integrand.
        of: ExprId,
        /// UQ configuration identity.
        uq: String,
    },
    /// Conditional value-at-risk at level β (kinked — RU reformulation is
    /// downstream's smoothing choice, the raw operator is kinked).
    CVaR {
        /// The tail quantity.
        of: ExprId,
        /// Tail level bits (canonical).
        beta_bits: u64,
    },
}

/// The arena: nodes + per-node propagated (dims, smoothness).
#[derive(Debug, Clone, Default)]
pub struct ExprGraph {
    nodes: Vec<ExprNode>,
    dims: Vec<Dims>,
    class: Vec<Smoothness>,
    intern: BTreeMap<ExprNode, ExprId>,
}

impl ExprGraph {
    /// An empty graph.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Node count.
    #[must_use]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// True when empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// The node behind an id.
    #[must_use]
    pub fn node(&self, id: ExprId) -> &ExprNode {
        &self.nodes[id.0 as usize]
    }

    /// Propagated dimensions of a node.
    #[must_use]
    pub fn dims(&self, id: ExprId) -> Dims {
        self.dims[id.0 as usize]
    }

    /// Propagated differentiability class of a node.
    #[must_use]
    pub fn class(&self, id: ExprId) -> Smoothness {
        self.class[id.0 as usize]
    }

    fn check(&self, id: ExprId) -> Result<(), OptDiag> {
        if (id.0 as usize) < self.nodes.len() {
            Ok(())
        } else {
            Err(OptDiag::UnknownNode { id: id.0 })
        }
    }

    /// Intern a node with its propagated properties (hash-consing: an
    /// identical structure returns the EXISTING id — the
    /// common-subexpression law).
    fn intern(&mut self, node: ExprNode, dims: Dims, class: Smoothness) -> ExprId {
        if let Some(&id) = self.intern.get(&node) {
            return id;
        }
        let id = ExprId(u32::try_from(self.nodes.len()).expect("graph < 2^32 nodes"));
        self.intern.insert(node.clone(), id);
        self.nodes.push(node);
        self.dims.push(dims);
        self.class.push(class);
        id
    }

    /// A variable-component leaf.
    pub(crate) fn var(&mut self, v: VarId, dims: Dims) -> ExprId {
        self.intern(ExprNode::Var(v), dims, Smoothness::Smooth)
    }

    /// A dimensioned constant.
    ///
    /// # Errors
    /// [`OptDiag::NonFinite`] for NaN/∞.
    pub fn constant(&mut self, value: f64, dims: Dims) -> Result<ExprId, OptDiag> {
        if !value.is_finite() {
            return Err(OptDiag::NonFinite { value });
        }
        Ok(self.intern(ExprNode::Const { bits: value.to_bits(), dims }, dims, Smoothness::Smooth))
    }

    /// Σ addends (≥ 2, same dims).
    ///
    /// # Errors
    /// [`OptDiag`] naming the offending node on dims mismatch.
    pub fn add(&mut self, terms: &[ExprId]) -> Result<ExprId, OptDiag> {
        if terms.len() < 2 {
            return Err(OptDiag::Arity { op: "add", need: 2, got: terms.len() });
        }
        for &t in terms {
            self.check(t)?;
        }
        let d0 = self.dims(terms[0]);
        for &t in &terms[1..] {
            if self.dims(t) != d0 {
                return Err(OptDiag::DimsMismatch {
                    op: "add",
                    node: t.0,
                    left: d0.unit_string(),
                    right: self.dims(t).unit_string(),
                });
            }
        }
        let class =
            terms.iter().map(|&t| self.class(t)).fold(Smoothness::Smooth, Smoothness::join);
        let mut sorted = terms.to_vec();
        sorted.sort_unstable(); // canonical: addition commutes
        Ok(self.intern(ExprNode::Add(sorted), d0, class))
    }

    /// Product (dims add).
    ///
    /// # Errors
    /// [`OptDiag::UnknownNode`] for dangling ids.
    pub fn mul(&mut self, a: ExprId, b: ExprId) -> Result<ExprId, OptDiag> {
        self.check(a)?;
        self.check(b)?;
        let dims = self.dims(a).plus(self.dims(b));
        let class = self.class(a).join(self.class(b));
        let (a, b) = if a <= b { (a, b) } else { (b, a) }; // canonical
        Ok(self.intern(ExprNode::Mul(a, b), dims, class))
    }

    /// Negation.
    ///
    /// # Errors
    /// [`OptDiag::UnknownNode`] for dangling ids.
    pub fn neg(&mut self, a: ExprId) -> Result<ExprId, OptDiag> {
        self.check(a)?;
        Ok(self.intern(ExprNode::Neg(a), self.dims(a), self.class(a)))
    }

    /// Integer power (dims scale by the exponent).
    ///
    /// # Errors
    /// [`OptDiag::UnknownNode`] for dangling ids.
    pub fn powi(&mut self, a: ExprId, e: i32) -> Result<ExprId, OptDiag> {
        self.check(a)?;
        let mut dims = Dims([0; 5]);
        let step = self.dims(a);
        for _ in 0..e.unsigned_abs() {
            dims = if e >= 0 { dims.plus(step) } else { dims.minus(step) };
        }
        Ok(self.intern(ExprNode::Powi(a, e), dims, self.class(a)))
    }

    /// Pointwise min (introduces a KINK).
    ///
    /// # Errors
    /// [`OptDiag`] on arity/dims defects.
    pub fn min(&mut self, terms: &[ExprId]) -> Result<ExprId, OptDiag> {
        self.kink("min", terms)
    }

    /// Pointwise max (introduces a KINK).
    ///
    /// # Errors
    /// [`OptDiag`] on arity/dims defects.
    pub fn max(&mut self, terms: &[ExprId]) -> Result<ExprId, OptDiag> {
        self.kink("max", terms)
    }

    fn kink(&mut self, op: &'static str, terms: &[ExprId]) -> Result<ExprId, OptDiag> {
        if terms.len() < 2 {
            return Err(OptDiag::Arity { op, need: 2, got: terms.len() });
        }
        for &t in terms {
            self.check(t)?;
        }
        let d0 = self.dims(terms[0]);
        for &t in &terms[1..] {
            if self.dims(t) != d0 {
                return Err(OptDiag::DimsMismatch {
                    op,
                    node: t.0,
                    left: d0.unit_string(),
                    right: self.dims(t).unit_string(),
                });
            }
        }
        let class = terms
            .iter()
            .map(|&t| self.class(t))
            .fold(Smoothness::Kinked, Smoothness::join);
        let mut sorted = terms.to_vec();
        sorted.sort_unstable();
        let node =
            if op == "min" { ExprNode::Min(sorted) } else { ExprNode::Max(sorted) };
        Ok(self.intern(node, d0, class))
    }

    /// |a| (kink at 0).
    ///
    /// # Errors
    /// [`OptDiag::UnknownNode`] for dangling ids.
    pub fn abs(&mut self, a: ExprId) -> Result<ExprId, OptDiag> {
        self.check(a)?;
        let class = self.class(a).join(Smoothness::Kinked);
        Ok(self.intern(ExprNode::Abs(a), self.dims(a), class))
    }

    /// A FLUX study QoI node: smooth iff a discrete adjoint is declared,
    /// otherwise OPAQUE (routed to derivative-free families).
    pub fn pde_qoi(&mut self, study: &str, dims: Dims, adjoint: bool) -> ExprId {
        let class = if adjoint { Smoothness::Smooth } else { Smoothness::Opaque };
        self.intern(
            ExprNode::PdeQoi { study: study.to_string(), dims, adjoint },
            dims,
            class,
        )
    }

    /// Expectation over a UQ configuration (class-preserving; the
    /// stochastic estimator's own noise is the e-process layer's concern).
    ///
    /// # Errors
    /// [`OptDiag::UnknownNode`] for dangling ids.
    pub fn expect(&mut self, of: ExprId, uq: &str) -> Result<ExprId, OptDiag> {
        self.check(of)?;
        let (dims, class) = (self.dims(of), self.class(of));
        Ok(self.intern(ExprNode::Expect { of, uq: uq.to_string() }, dims, class))
    }

    /// CVaR_β (kinked; the smooth RU reformulation is an optimizer-side
    /// choice, not a property of the operator).
    ///
    /// # Errors
    /// [`OptDiag`] on dangling ids or β outside (0, 1).
    pub fn cvar(&mut self, of: ExprId, beta: f64) -> Result<ExprId, OptDiag> {
        self.check(of)?;
        if !(beta.is_finite() && beta > 0.0 && beta < 1.0) {
            return Err(OptDiag::BadBeta { beta });
        }
        let dims = self.dims(of);
        let class = self.class(of).join(Smoothness::Kinked);
        Ok(self.intern(ExprNode::CVaR { of, beta_bits: beta.to_bits() }, dims, class))
    }

    /// Substitute every occurrence of variable `v` with expression `e`
    /// (rebuilds bottom-up; the substitution law of the G0 suite).
    ///
    /// # Errors
    /// Propagates rebuild diagnostics (cannot fail on well-formed graphs).
    pub fn substitute(&mut self, root: ExprId, v: VarId, e: ExprId) -> Result<ExprId, OptDiag> {
        self.check(root)?;
        self.check(e)?;
        let mut memo: BTreeMap<ExprId, ExprId> = BTreeMap::new();
        self.subst_inner(root, v, e, &mut memo)
    }

    fn subst_inner(
        &mut self,
        id: ExprId,
        v: VarId,
        e: ExprId,
        memo: &mut BTreeMap<ExprId, ExprId>,
    ) -> Result<ExprId, OptDiag> {
        if let Some(&m) = memo.get(&id) {
            return Ok(m);
        }
        let node = self.node(id).clone();
        let out = match node {
            ExprNode::Var(w) if w == v => e,
            ExprNode::Var(_) | ExprNode::Const { .. } | ExprNode::PdeQoi { .. } => id,
            ExprNode::Add(ts) => {
                let mapped: Result<Vec<ExprId>, OptDiag> =
                    ts.iter().map(|&t| self.subst_inner(t, v, e, memo)).collect();
                self.add(&mapped?)?
            }
            ExprNode::Mul(a, b) => {
                let (a, b) =
                    (self.subst_inner(a, v, e, memo)?, self.subst_inner(b, v, e, memo)?);
                self.mul(a, b)?
            }
            ExprNode::Neg(a) => {
                let a = self.subst_inner(a, v, e, memo)?;
                self.neg(a)?
            }
            ExprNode::Powi(a, p) => {
                let a = self.subst_inner(a, v, e, memo)?;
                self.powi(a, p)?
            }
            ExprNode::Min(ts) => {
                let mapped: Result<Vec<ExprId>, OptDiag> =
                    ts.iter().map(|&t| self.subst_inner(t, v, e, memo)).collect();
                self.min(&mapped?)?
            }
            ExprNode::Max(ts) => {
                let mapped: Result<Vec<ExprId>, OptDiag> =
                    ts.iter().map(|&t| self.subst_inner(t, v, e, memo)).collect();
                self.max(&mapped?)?
            }
            ExprNode::Abs(a) => {
                let a = self.subst_inner(a, v, e, memo)?;
                self.abs(a)?
            }
            ExprNode::Expect { of, uq } => {
                let of = self.subst_inner(of, v, e, memo)?;
                self.expect(of, &uq)?
            }
            ExprNode::CVaR { of, beta_bits } => {
                let of = self.subst_inner(of, v, e, memo)?;
                self.cvar(of, f64::from_bits(beta_bits))?
            }
        };
        memo.insert(id, out);
        Ok(out)
    }
}

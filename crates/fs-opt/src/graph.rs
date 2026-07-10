//! The typed expression graph: arena nodes with hash-consing (structural
//! deduplication makes the common-subexpression identity hold BY
//! CONSTRUCTION), dimension checking on every constructor, bottom-up
//! differentiability-class propagation, exact reverse-mode gradients with
//! deterministic accumulation order, and structure-only nodes (PDE
//! residuals, stochastic operators) that are representable and validated
//! today and evaluable when FLUX/UQ land.

use crate::manifold::Manifold;
use fs_qty::Dims;
use std::collections::HashMap;

/// Node handle (arena index).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub u32);

/// Variable handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VarId(pub u32);

/// Differentiability class, propagated bottom-up (worst-of-children
/// composed with the op's own class).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DiffClass {
    /// C¹ everywhere on the evaluated domain.
    Smooth,
    /// Continuous with kinks (min/max/abs) — subgradient-safe only.
    NonSmooth,
    /// Not differentiable / not yet evaluable (structure placeholders).
    NonDiff,
}

/// Optimizer families for routing validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimizerFamily {
    /// Smooth-only line-search methods (L-BFGS, TR-Newton).
    SmoothGradient,
    /// Subgradient/bundle/proximal: tolerates kinks.
    Subgradient,
    /// Derivative-free (CMA-ES, BO): tolerates NonDiff evaluable nodes.
    DerivativeFree,
}

/// Structured validation failure naming the offending node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidationError {
    /// Dimensions disagree where they must match.
    DimensionMismatch {
        /// The node being constructed/validated.
        node: NodeId,
        /// Human-readable op name.
        op: &'static str,
        /// Left/first operand dims (unit string).
        left: String,
        /// Right/second operand dims (unit string).
        right: String,
    },
    /// A transcendental applied to a dimensioned quantity.
    DimensionedTranscendental {
        /// Offending node.
        node: NodeId,
        /// Op name.
        op: &'static str,
        /// The argument's unit string.
        arg: String,
    },
    /// The graph's class exceeds what the optimizer family accepts.
    ClassTooRough {
        /// First offending node on a worst-class path.
        node: NodeId,
        /// Op that introduced the roughness.
        op: &'static str,
        /// The class found.
        found: DiffClass,
        /// The family that rejected it.
        family: OptimizerFamily,
    },
}

/// Node payload. Structure-only variants carry metadata, not values.
#[derive(Debug, Clone, PartialEq)]
pub enum Node {
    /// Constant with dimensions.
    Const(f64, Dims),
    /// One component of a variable.
    Component(VarId, u32),
    /// Σᵢ vᵢ² of a whole variable (dimensionless variables only in v1).
    NormSq(VarId),
    /// ⟨u, v⟩ of two same-length variables.
    Dot(VarId, VarId),
    /// Binary arithmetic.
    Add(NodeId, NodeId),
    /// Subtraction.
    Sub(NodeId, NodeId),
    /// Multiplication.
    Mul(NodeId, NodeId),
    /// Division.
    Div(NodeId, NodeId),
    /// Negation.
    Neg(NodeId),
    /// Integer power.
    Powi(NodeId, i8),
    /// Pointwise minimum (NonSmooth).
    Min(NodeId, NodeId),
    /// Pointwise maximum (NonSmooth).
    Max(NodeId, NodeId),
    /// Absolute value (NonSmooth).
    Abs(NodeId),
    /// sin (dimensionless).
    Sin(NodeId),
    /// cos (dimensionless).
    Cos(NodeId),
    /// exp (dimensionless).
    Exp(NodeId),
    /// ln (dimensionless).
    Ln(NodeId),
    /// sqrt (dims must halve evenly).
    Sqrt(NodeId),
    /// tanh (dimensionless).
    Tanh(NodeId),
    /// PDE residual placeholder: physics(u, θ) = 0 as a first-class node.
    PdeResidual {
        /// FLUX study identity (content hash once FLUX exists).
        study: u64,
        /// Whether an adjoint is available (routing metadata).
        adjoint_available: bool,
    },
    /// Expectation of an inner expression under a UQ configuration.
    Expectation {
        /// The integrand.
        inner: NodeId,
        /// UQ configuration identity.
        uq_config: u64,
    },
}

/// A typed optimization problem under incremental construction.
#[derive(Debug, Clone, Default)]
pub struct Problem {
    nodes: Vec<Node>,
    dims: Vec<Dims>,
    class: Vec<DiffClass>,
    vars: Vec<(String, Manifold)>,
    /// Hash-consing table: structural key → existing node.
    dedup: HashMap<String, NodeId>,
    /// The designated objective (last set).
    objective: Option<NodeId>,
}

/// Evaluation failure (structure-only nodes).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Unevaluable {
    /// The node that cannot be evaluated yet.
    pub node: NodeId,
}

const DIMLESS: Dims = Dims([0; 5]);

impl Problem {
    /// Empty problem.
    #[must_use]
    pub fn new() -> Problem {
        Problem::default()
    }

    /// Declare a variable on a manifold; the name is for diagnostics and
    /// serialization (must be unique — panics otherwise; a modeling
    /// error).
    pub fn variable(&mut self, name: &str, manifold: Manifold) -> VarId {
        assert!(
            self.vars.iter().all(|(n, _)| n != name),
            "variable name '{name}' already declared"
        );
        self.vars.push((name.to_string(), manifold));
        VarId(u32::try_from(self.vars.len() - 1).expect("var count"))
    }

    /// Variable metadata.
    #[must_use]
    pub fn manifold(&self, v: VarId) -> &Manifold {
        &self.vars[v.0 as usize].1
    }

    /// Number of variables.
    #[must_use]
    pub fn var_count(&self) -> usize {
        self.vars.len()
    }

    /// Set the objective node.
    pub fn set_objective(&mut self, node: NodeId) {
        self.objective = Some(node);
    }

    /// The objective node (panics if unset — a modeling error).
    #[must_use]
    pub fn objective(&self) -> NodeId {
        self.objective.expect("objective not set")
    }

    fn push(&mut self, node: Node, dims: Dims, class: DiffClass) -> NodeId {
        // Structural hash-consing: identical payloads collapse to one id
        // (CSE identity BY CONSTRUCTION).
        let key = format!("{node:?}");
        if let Some(&id) = self.dedup.get(&key) {
            return id;
        }
        let id = NodeId(u32::try_from(self.nodes.len()).expect("node count"));
        self.nodes.push(node);
        self.dims.push(dims);
        self.class.push(class);
        self.dedup.insert(key, id);
        id
    }

    /// Dimensioned constant.
    pub fn constant(&mut self, v: f64, dims: Dims) -> NodeId {
        self.push(Node::Const(v, dims), dims, DiffClass::Smooth)
    }

    /// Dimensionless constant.
    pub fn scalar(&mut self, v: f64) -> NodeId {
        self.constant(v, DIMLESS)
    }

    /// Component leaf (dimensionless variables in v1 — Qty-valued
    /// variables are recorded follow-up scope).
    pub fn component(&mut self, v: VarId, idx: u32) -> NodeId {
        let n = self.vars[v.0 as usize].1.ambient_dim();
        assert!((idx as usize) < n, "component {idx} out of range for dim {n}");
        self.push(Node::Component(v, idx), DIMLESS, DiffClass::Smooth)
    }

    /// ‖v‖² convenience leaf.
    pub fn norm_sq(&mut self, v: VarId) -> NodeId {
        self.push(Node::NormSq(v), DIMLESS, DiffClass::Smooth)
    }

    /// ⟨u, v⟩ convenience leaf (same ambient dimension required).
    pub fn dot(&mut self, u: VarId, v: VarId) -> NodeId {
        assert_eq!(
            self.vars[u.0 as usize].1.ambient_dim(),
            self.vars[v.0 as usize].1.ambient_dim(),
            "dot requires equal ambient dimensions"
        );
        self.push(Node::Dot(u, v), DIMLESS, DiffClass::Smooth)
    }

    fn class2(&self, a: NodeId, b: NodeId) -> DiffClass {
        self.class[a.0 as usize].max(self.class[b.0 as usize])
    }

    /// a + b (dims must match).
    ///
    /// # Errors
    /// [`ValidationError::DimensionMismatch`] naming the node.
    pub fn add(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, ValidationError> {
        self.dim_match("add", a, b)?;
        let d = self.dims[a.0 as usize];
        let c = self.class2(a, b);
        Ok(self.push(Node::Add(a, b), d, c))
    }

    /// a − b (dims must match).
    ///
    /// # Errors
    /// [`ValidationError::DimensionMismatch`].
    pub fn sub(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, ValidationError> {
        self.dim_match("sub", a, b)?;
        let d = self.dims[a.0 as usize];
        let c = self.class2(a, b);
        Ok(self.push(Node::Sub(a, b), d, c))
    }

    /// a·b (dims compose).
    pub fn mul(&mut self, a: NodeId, b: NodeId) -> NodeId {
        let d = self.dims[a.0 as usize].plus(self.dims[b.0 as usize]);
        let c = self.class2(a, b);
        self.push(Node::Mul(a, b), d, c)
    }

    /// a/b (dims compose).
    pub fn div(&mut self, a: NodeId, b: NodeId) -> NodeId {
        let d = self.dims[a.0 as usize].minus(self.dims[b.0 as usize]);
        let c = self.class2(a, b);
        self.push(Node::Div(a, b), d, c)
    }

    /// −a.
    pub fn neg(&mut self, a: NodeId) -> NodeId {
        let d = self.dims[a.0 as usize];
        let c = self.class[a.0 as usize];
        self.push(Node::Neg(a), d, c)
    }

    /// aⁿ (dims scale by n).
    pub fn powi(&mut self, a: NodeId, n: i8) -> NodeId {
        let d = self.dims[a.0 as usize].times(n);
        let c = self.class[a.0 as usize];
        self.push(Node::Powi(a, n), d, c)
    }

    /// min(a, b): NonSmooth by construction.
    ///
    /// # Errors
    /// [`ValidationError::DimensionMismatch`].
    pub fn min(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, ValidationError> {
        self.dim_match("min", a, b)?;
        let d = self.dims[a.0 as usize];
        let c = self.class2(a, b).max(DiffClass::NonSmooth);
        Ok(self.push(Node::Min(a, b), d, c))
    }

    /// max(a, b): NonSmooth by construction.
    ///
    /// # Errors
    /// [`ValidationError::DimensionMismatch`].
    pub fn max(&mut self, a: NodeId, b: NodeId) -> Result<NodeId, ValidationError> {
        self.dim_match("max", a, b)?;
        let d = self.dims[a.0 as usize];
        let c = self.class2(a, b).max(DiffClass::NonSmooth);
        Ok(self.push(Node::Max(a, b), d, c))
    }

    /// |a|: NonSmooth by construction.
    pub fn abs(&mut self, a: NodeId) -> NodeId {
        let d = self.dims[a.0 as usize];
        let c = self.class[a.0 as usize].max(DiffClass::NonSmooth);
        self.push(Node::Abs(a), d, c)
    }

    /// Transcendentals (dimensionless argument enforced).
    ///
    /// # Errors
    /// [`ValidationError::DimensionedTranscendental`].
    pub fn unary(
        &mut self,
        op: &'static str,
        a: NodeId,
    ) -> Result<NodeId, ValidationError> {
        let d = self.dims[a.0 as usize];
        if d != DIMLESS && op != "sqrt" {
            return Err(ValidationError::DimensionedTranscendental {
                node: a,
                op,
                arg: d.unit_string(),
            });
        }
        let c = self.class[a.0 as usize];
        let (node, out_d) = match op {
            "sin" => (Node::Sin(a), DIMLESS),
            "cos" => (Node::Cos(a), DIMLESS),
            "exp" => (Node::Exp(a), DIMLESS),
            "ln" => (Node::Ln(a), DIMLESS),
            "tanh" => (Node::Tanh(a), DIMLESS),
            "sqrt" => {
                // Dims must halve evenly.
                let ok = d.0.iter().all(|&e| e % 2 == 0);
                if !ok {
                    return Err(ValidationError::DimensionedTranscendental {
                        node: a,
                        op,
                        arg: d.unit_string(),
                    });
                }
                let half = Dims([d.0[0] / 2, d.0[1] / 2, d.0[2] / 2, d.0[3] / 2, d.0[4] / 2]);
                (Node::Sqrt(a), half)
            }
            _ => panic!("unknown unary op {op}"),
        };
        Ok(self.push(node, out_d, c))
    }

    /// PDE-residual structure node (representable now, evaluable when
    /// FLUX lands). Class: NonDiff unless an adjoint is declared.
    pub fn pde_residual(&mut self, study: u64, adjoint_available: bool) -> NodeId {
        let c = if adjoint_available { DiffClass::Smooth } else { DiffClass::NonDiff };
        self.push(Node::PdeResidual { study, adjoint_available }, DIMLESS, c)
    }

    /// Expectation structure node over an inner expression.
    pub fn expectation(&mut self, inner: NodeId, uq_config: u64) -> NodeId {
        let c = self.class[inner.0 as usize];
        let d = self.dims[inner.0 as usize];
        self.push(Node::Expectation { inner, uq_config }, d, c)
    }

    fn dim_match(
        &self,
        op: &'static str,
        a: NodeId,
        b: NodeId,
    ) -> Result<(), ValidationError> {
        let (da, db) = (self.dims[a.0 as usize], self.dims[b.0 as usize]);
        if da == db {
            Ok(())
        } else {
            Err(ValidationError::DimensionMismatch {
                node: NodeId(u32::try_from(self.nodes.len()).expect("node count")),
                op,
                left: da.unit_string(),
                right: db.unit_string(),
            })
        }
    }

    /// The propagated differentiability class of a node.
    #[must_use]
    pub fn diff_class(&self, n: NodeId) -> DiffClass {
        self.class[n.0 as usize]
    }

    /// Validate the objective for an optimizer family: structured
    /// diagnostics name the first offending node on a worst-class path.
    ///
    /// # Errors
    /// [`ValidationError::ClassTooRough`].
    pub fn validate_for(&self, family: OptimizerFamily) -> Result<(), ValidationError> {
        let obj = self.objective();
        let worst = self.class[obj.0 as usize];
        let limit = match family {
            OptimizerFamily::SmoothGradient => DiffClass::Smooth,
            OptimizerFamily::Subgradient => DiffClass::NonSmooth,
            OptimizerFamily::DerivativeFree => DiffClass::NonDiff,
        };
        if worst <= limit {
            return Ok(());
        }
        // Find the first node that introduced the offending class.
        for (i, c) in self.class.iter().enumerate() {
            if *c == worst {
                let op = match &self.nodes[i] {
                    Node::Min(..) => "min",
                    Node::Max(..) => "max",
                    Node::Abs(..) => "abs",
                    Node::PdeResidual { .. } => "pde_residual(no adjoint)",
                    _ => "node",
                };
                return Err(ValidationError::ClassTooRough {
                    node: NodeId(u32::try_from(i).expect("idx")),
                    op,
                    found: worst,
                    family,
                });
            }
        }
        unreachable!("worst class must exist in the table");
    }

    /// Evaluate the objective at a point (map var → component values).
    ///
    /// # Errors
    /// [`Unevaluable`] on structure-only nodes.
    pub fn eval(&self, x: &HashMap<VarId, Vec<f64>>) -> Result<f64, Unevaluable> {
        let obj = self.objective();
        let mut vals = vec![0.0f64; self.nodes.len()];
        for i in 0..=obj.0 as usize {
            vals[i] = self.eval_node(NodeId(u32::try_from(i).expect("idx")), &vals, x)?;
        }
        Ok(vals[obj.0 as usize])
    }

    #[allow(clippy::too_many_lines)]
    fn eval_node(
        &self,
        id: NodeId,
        vals: &[f64],
        x: &HashMap<VarId, Vec<f64>>,
    ) -> Result<f64, Unevaluable> {
        use fs_math::det;
        let v = |n: NodeId| vals[n.0 as usize];
        Ok(match &self.nodes[id.0 as usize] {
            Node::Const(c, _) => *c,
            Node::Component(var, idx) => x[var][*idx as usize],
            Node::NormSq(var) => x[var].iter().map(|t| t * t).sum(),
            Node::Dot(u, w) => x[u].iter().zip(&x[w]).map(|(a, b)| a * b).sum(),
            Node::Add(a, b) => v(*a) + v(*b),
            Node::Sub(a, b) => v(*a) - v(*b),
            Node::Mul(a, b) => v(*a) * v(*b),
            Node::Div(a, b) => v(*a) / v(*b),
            Node::Neg(a) => -v(*a),
            Node::Powi(a, n) => det::powi(v(*a), i32::from(*n)),
            Node::Min(a, b) => v(*a).min(v(*b)),
            Node::Max(a, b) => v(*a).max(v(*b)),
            Node::Abs(a) => v(*a).abs(),
            Node::Sin(a) => det::sin(v(*a)),
            Node::Cos(a) => det::cos(v(*a)),
            Node::Exp(a) => det::exp(v(*a)),
            Node::Ln(a) => det::ln(v(*a)),
            Node::Sqrt(a) => det::sqrt(v(*a)),
            Node::Tanh(a) => det::tanh(v(*a)),
            Node::PdeResidual { .. } | Node::Expectation { .. } => {
                return Err(Unevaluable { node: id });
            }
        })
    }

    /// Exact reverse-mode gradient of the objective w.r.t. every variable
    /// component (deterministic accumulation: forward pass by node index,
    /// reverse pass in exact reverse order).
    ///
    /// # Errors
    /// [`Unevaluable`] on structure-only nodes.
    pub fn gradient(
        &self,
        x: &HashMap<VarId, Vec<f64>>,
    ) -> Result<HashMap<VarId, Vec<f64>>, Unevaluable> {
        use fs_math::det;
        let obj = self.objective();
        let n = obj.0 as usize + 1;
        let mut vals = vec![0.0f64; self.nodes.len()];
        for i in 0..n {
            vals[i] = self.eval_node(NodeId(u32::try_from(i).expect("idx")), &vals, x)?;
        }
        let mut bar = vec![0.0f64; n];
        bar[obj.0 as usize] = 1.0;
        let mut grads: HashMap<VarId, Vec<f64>> = x
            .iter()
            .map(|(k, v)| (*k, vec![0.0; v.len()]))
            .collect();
        for i in (0..n).rev() {
            let b = bar[i];
            if b == 0.0 {
                continue;
            }
            let v = |n: NodeId| vals[n.0 as usize];
            match &self.nodes[i] {
                Node::Const(..) => {}
                Node::Component(var, idx) => {
                    grads.get_mut(var).expect("var present")[*idx as usize] += b;
                }
                Node::NormSq(var) => {
                    let g = grads.get_mut(var).expect("var present");
                    for (gi, &xi) in g.iter_mut().zip(&x[var]) {
                        *gi = (2.0 * b).mul_add(xi, *gi);
                    }
                }
                Node::Dot(u, w) => {
                    if u == w {
                        let g = grads.get_mut(u).expect("var");
                        for (gi, &xi) in g.iter_mut().zip(&x[u]) {
                            *gi = (2.0 * b).mul_add(xi, *gi);
                        }
                    } else {
                        let xw: Vec<f64> = x[w].clone();
                        let xu: Vec<f64> = x[u].clone();
                        let gu = grads.get_mut(u).expect("var");
                        for (gi, &wi) in gu.iter_mut().zip(&xw) {
                            *gi = b.mul_add(wi, *gi);
                        }
                        let gw = grads.get_mut(w).expect("var");
                        for (gi, &ui) in gw.iter_mut().zip(&xu) {
                            *gi = b.mul_add(ui, *gi);
                        }
                    }
                }
                Node::Add(a, c) => {
                    bar[a.0 as usize] += b;
                    bar[c.0 as usize] += b;
                }
                Node::Sub(a, c) => {
                    bar[a.0 as usize] += b;
                    bar[c.0 as usize] -= b;
                }
                Node::Mul(a, c) => {
                    bar[a.0 as usize] = b.mul_add(v(*c), bar[a.0 as usize]);
                    bar[c.0 as usize] = b.mul_add(v(*a), bar[c.0 as usize]);
                }
                Node::Div(a, c) => {
                    let vc = v(*c);
                    bar[a.0 as usize] += b / vc;
                    bar[c.0 as usize] -= b * v(*a) / (vc * vc);
                }
                Node::Neg(a) => bar[a.0 as usize] -= b,
                Node::Powi(a, k) => {
                    let va = v(*a);
                    bar[a.0 as usize] = (b * f64::from(*k))
                        .mul_add(det::powi(va, i32::from(*k) - 1), bar[a.0 as usize]);
                }
                Node::Min(a, c) => {
                    // Subgradient convention: the SMALLER branch gets the
                    // pull; exact tie → lower node id (deterministic).
                    let pick = if v(*a) < v(*c) || (v(*a) == v(*c) && a <= c) { a } else { c };
                    bar[pick.0 as usize] += b;
                }
                Node::Max(a, c) => {
                    let pick = if v(*a) > v(*c) || (v(*a) == v(*c) && a <= c) { a } else { c };
                    bar[pick.0 as usize] += b;
                }
                Node::Abs(a) => {
                    let s = if v(*a) >= 0.0 { 1.0 } else { -1.0 };
                    bar[a.0 as usize] = (b * s).mul_add(1.0, bar[a.0 as usize]);
                }
                Node::Sin(a) => {
                    bar[a.0 as usize] = b.mul_add(det::cos(v(*a)), bar[a.0 as usize]);
                }
                Node::Cos(a) => {
                    bar[a.0 as usize] = (-b).mul_add(det::sin(v(*a)), bar[a.0 as usize]);
                }
                Node::Exp(a) => {
                    bar[a.0 as usize] = b.mul_add(vals[i], bar[a.0 as usize]);
                }
                Node::Ln(a) => bar[a.0 as usize] += b / v(*a),
                Node::Sqrt(a) => bar[a.0 as usize] += b / (2.0 * vals[i]),
                Node::Tanh(a) => {
                    let t = vals[i];
                    bar[a.0 as usize] = b.mul_add(t.mul_add(-t, 1.0), bar[a.0 as usize]);
                }
                Node::PdeResidual { .. } | Node::Expectation { .. } => {
                    return Err(Unevaluable {
                        node: NodeId(u32::try_from(i).expect("idx")),
                    });
                }
            }
        }
        Ok(grads)
    }

    /// Raw node table access for serialization.
    #[must_use]
    pub(crate) fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// Variable table access for serialization.
    #[must_use]
    pub(crate) fn vars(&self) -> &[(String, Manifold)] {
        &self.vars
    }

    /// Re-import a node during parsing (validation re-runs).
    pub(crate) fn push_raw(&mut self, node: Node, dims: Dims, class: DiffClass) -> NodeId {
        self.push(node, dims, class)
    }

    /// Node dims (for serialization).
    #[must_use]
    pub(crate) fn dims_of(&self, n: NodeId) -> Dims {
        self.dims[n.0 as usize]
    }
}

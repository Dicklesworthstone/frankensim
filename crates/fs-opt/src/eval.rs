//! Evaluation + the TOY Riemannian descent that proves the manifold
//! metadata is consumable: retractions keep iterates ON their
//! manifolds ("optimize an orientation" never becomes "optimize 9
//! numbers and renormalize when it explodes"). Gradients here are
//! finite differences through retractions — a deliberately simple
//! consumer; exact adjoints are the gradient-stack bead's.

use crate::admission::AdmissionCaps;
use crate::ir::{Expr, Manifold, NodeId, OptError, Problem};
use fs_exec::Cx;

/// Squared-norm/dot tolerance for deciding whether a finite stored
/// point is already on its declared manifold.
const MANIFOLD_MEMBERSHIP_TOL: f64 = 1e-10;
/// Candidates below this squared norm are numerically rank-deficient;
/// normalizing them would fabricate a direction.
const RETRACTION_MIN_NORM_SQ: f64 = 1e-24;

/// An evaluated node value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Scalar.
    S(f64),
    /// Vector.
    V(Vec<f64>),
}

impl Value {
    /// The scalar, if it is one.
    #[must_use]
    pub fn scalar(&self) -> Option<f64> {
        match self {
            Value::S(v) => Some(*v),
            Value::V(_) => None,
        }
    }
}

/// Evaluate `node` with variable `bindings` (one point per variable,
/// stored per its manifold; bindings are indexed by `VarId` and may be
/// a PREFIX of the declaration list — referencing an unbound variable
/// teaches with `UnknownVar`). PDE and stochastic nodes are NOT
/// evaluable here — they refuse with a teaching error (their execution
/// belongs to FLUX/UQ runners; the IR only carries them).
///
/// # Errors
/// [`OptError::Unevaluable`] / [`OptError::UnknownVar`] /
/// [`OptError::UnknownNode`] / [`OptError::BindingCount`] /
/// [`OptError::BindingLen`] / [`OptError::BindingNonFinite`] /
/// [`OptError::BindingDomain`] / [`OptError::EvalNonFinite`] /
/// [`OptError::CapExceeded`].
pub fn eval(problem: &Problem, node: NodeId, bindings: &[Vec<f64>]) -> Result<Value, OptError> {
    if node.0 as usize >= problem.exprs.len() {
        return Err(OptError::UnknownNode { id: node.0 });
    }
    let caps = AdmissionCaps::default();
    let depth = problem.node_depth(node)?;
    if depth > caps.max_graph_depth {
        return Err(OptError::CapExceeded {
            what: "graph depth",
            count: u64::from(depth),
            cap: u64::from(caps.max_graph_depth),
        });
    }
    let work = problem.total_admission_work();
    if work > caps.max_total_work {
        return Err(OptError::CapExceeded {
            what: "total admission work",
            count: work,
            cap: caps.max_total_work,
        });
    }
    if bindings.len() > problem.vars.len() {
        return Err(OptError::BindingCount {
            vars: problem.vars.len() as u32,
            got: bindings.len() as u64,
        });
    }
    for (i, (binding, var)) in bindings.iter().zip(&problem.vars).enumerate() {
        let expected = var
            .manifold
            .point_dim()
            .expect("sealed problems carry validated manifolds");
        if binding.len() as u64 != u64::from(expected) {
            return Err(OptError::BindingLen {
                var: i as u32,
                expected,
                got: binding.len() as u64,
            });
        }
        for (component_index, component) in binding.iter().enumerate() {
            if !component.is_finite() {
                return Err(OptError::BindingNonFinite {
                    var: i as u32,
                    component: component_index as u32,
                    bits: component.to_bits(),
                });
            }
        }
        var.manifold.validate_binding_domain(binding, i as u32)?;
    }
    // Arena order guarantees every dependency has a lower id, so a
    // prefix is sufficient; unrelated later nodes cannot force memo
    // allocation for this evaluation.
    //
    // EXPLICIT-STACK EVALUATION (bead frankensim-xf8v7): the walk is
    // physically iterative so no admitted graph — at the depth cap or
    // otherwise — can overflow the call stack. Reachability first (a
    // worklist; children always have lower ids), then one bottom-up
    // sweep in arena order, which is a topological order by
    // construction. Only reachable nodes are evaluated, so unreachable
    // Unevaluable nodes still cannot poison an evaluation, exactly as
    // under the recursive walk.
    let root = node.0 as usize;
    let mut memo: Vec<Option<Value>> = vec![None; root + 1];
    let mut reachable = vec![false; root + 1];
    let mut worklist = vec![node];
    while let Some(n) = worklist.pop() {
        let i = n.0 as usize;
        if !reachable[i] {
            reachable[i] = true;
            worklist.extend(crate::ir::children(&problem.exprs[i]));
        }
    }
    for i in 0..=root {
        if reachable[i] {
            let value = eval_node(problem, NodeId(i as u32), bindings, &memo)?;
            memo[i] = Some(value);
        }
    }
    Ok(memo[root]
        .take()
        .expect("the root is reachable from itself"))
}

/// Evaluate ONE node whose children are already in `memo` (guaranteed
/// by the bottom-up arena-order sweep in [`eval`]). Never recurses.
#[allow(clippy::too_many_lines)] // one arm per node kind: the evaluator IS the semantics
fn eval_node(
    problem: &Problem,
    node: NodeId,
    bindings: &[Vec<f64>],
    memo: &[Option<Value>],
) -> Result<Value, OptError> {
    let ev = |n: NodeId| -> Value {
        memo[n.0 as usize]
            .clone()
            .expect("arena order: children are evaluated before their parents")
    };
    let scalar = |v: Value| -> f64 {
        match v {
            Value::S(x) => x,
            Value::V(_) => unreachable!("builder enforced scalar shape"),
        }
    };
    let out = match &problem.exprs[node.0 as usize] {
        Expr::Var(v) => {
            let x = bindings
                .get(v.0 as usize)
                .ok_or(OptError::UnknownVar { id: v.0 })?;
            Value::V(x.clone())
        }
        Expr::Component { of, index } => {
            let v = ev(*of);
            match v {
                Value::V(xs) => Value::S(xs[*index as usize]),
                Value::S(_) => unreachable!("builder enforced vector shape"),
            }
        }
        Expr::Const { value, .. } => Value::S(*value),
        Expr::Add(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x + y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p + q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Sub(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x - y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p - q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Mul(a, b) => match (ev(*a), ev(*b)) {
            (Value::S(x), Value::S(y)) => Value::S(x * y),
            (Value::S(s), Value::V(v)) | (Value::V(v), Value::S(s)) => {
                Value::V(v.iter().map(|p| p * s).collect())
            }
            _ => unreachable!("builder rejected vector*vector"),
        },
        Expr::Div(a, b) => {
            let (x, y) = (scalar(ev(*a)), scalar(ev(*b)));
            Value::S(x / y)
        }
        Expr::Neg(a) => match ev(*a) {
            Value::S(x) => Value::S(-x),
            Value::V(v) => Value::V(v.iter().map(|p| -p).collect()),
        },
        Expr::Powi { base, exp } => Value::S(fs_math::det::powi(scalar(ev(*base)), *exp)),
        Expr::Sqrt(a) => Value::S(fs_math::det::sqrt(scalar(ev(*a)))),
        Expr::Exp(a) => Value::S(fs_math::det::exp(scalar(ev(*a)))),
        Expr::Ln(a) => Value::S(fs_math::det::ln(scalar(ev(*a)))),
        Expr::Tanh(a) => Value::S(fs_math::det::tanh(scalar(ev(*a)))),
        Expr::Dot(a, b) => match (ev(*a), ev(*b)) {
            (Value::V(x), Value::V(y)) => Value::S(x.iter().zip(&y).map(|(p, q)| p * q).sum()),
            _ => unreachable!("builder enforced vectors"),
        },
        Expr::NormSq(a) => match ev(*a) {
            Value::V(x) => Value::S(x.iter().map(|p| p * p).sum()),
            Value::S(_) => unreachable!("builder enforced vector"),
        },
        Expr::Min(a, b) => Value::S(scalar(ev(*a)).min(scalar(ev(*b)))),
        Expr::Max(a, b) => Value::S(scalar(ev(*a)).max(scalar(ev(*b)))),
        Expr::Abs(a) => Value::S(scalar(ev(*a)).abs()),
        Expr::PdeResidual { study, .. } => {
            return Err(OptError::Unevaluable {
                node: node.0,
                kind: format!("pde_residual `{study}` (FLUX executes physics, not the IR)"),
            });
        }
        Expr::Expectation { uq_config, .. }
        | Expr::Cvar { uq_config, .. }
        | Expr::Quantile { uq_config, .. } => {
            return Err(OptError::Unevaluable {
                node: node.0,
                kind: format!("stochastic node over `{uq_config}` (UQ runners execute these)"),
            });
        }
    };
    match &out {
        Value::S(value) if !value.is_finite() => {
            return Err(OptError::EvalNonFinite {
                node: node.0,
                component: None,
                bits: value.to_bits(),
            });
        }
        Value::V(values) => {
            if let Some((component, value)) = values
                .iter()
                .enumerate()
                .find(|(_, value)| !value.is_finite())
            {
                return Err(OptError::EvalNonFinite {
                    node: node.0,
                    component: Some(component as u32),
                    bits: value.to_bits(),
                });
            }
        }
        Value::S(_) => {}
    }
    Ok(out)
}

impl Manifold {
    /// Descent parameter dimension (what the FD gradient has): ambient
    /// storage for Rn/Sphere/Stiefel (projection happens inside the
    /// retraction), axis-angle 3 for SO(3). CHECKED like
    /// [`Manifold::point_dim`]; `None` only for descriptors a sealed
    /// problem can never contain.
    #[must_use]
    pub fn param_dim(&self) -> Option<u32> {
        match *self {
            Manifold::So3 => Some(3),
            m => m.point_dim(),
        }
    }

    fn label(&self) -> &'static str {
        match self {
            Manifold::Rn { .. } => "Rn",
            Manifold::Sphere { .. } => "Sphere",
            Manifold::So3 => "SO(3)",
            Manifold::Stiefel { .. } => "Stiefel",
        }
    }

    fn domain_error(
        &self,
        what: &'static str,
        location: Option<(u32, u32)>,
        measurement: f64,
    ) -> OptError {
        OptError::RetractionDomain {
            manifold: self.label(),
            what,
            location,
            measurement_bits: measurement.to_bits(),
        }
    }

    fn validate_retraction_finite(input: &'static str, values: &[f64]) -> Result<(), OptError> {
        for (component, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(OptError::RetractionNonFinite {
                    input,
                    component: component as u32,
                    bits: value.to_bits(),
                });
            }
        }
        Ok(())
    }

    fn validate_point_domain(&self, x: &[f64]) -> Result<(), OptError> {
        match *self {
            Manifold::Rn { .. } => Ok(()),
            Manifold::Sphere { .. } | Manifold::So3 => {
                let norm_sq = x.iter().map(|value| value * value).sum::<f64>();
                if !norm_sq.is_finite() {
                    return Err(self.domain_error(
                        "point norm squared must be finite",
                        None,
                        norm_sq,
                    ));
                }
                if (norm_sq - 1.0).abs() > MANIFOLD_MEMBERSHIP_TOL {
                    return Err(self.domain_error(
                        "point must have unit norm before retraction",
                        None,
                        norm_sq,
                    ));
                }
                Ok(())
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                for column in 0..p {
                    for against in 0..=column {
                        let dot = (0..n)
                            .map(|row| x[column * n + row] * x[against * n + row])
                            .sum::<f64>();
                        let location = Some((column as u32, against as u32));
                        if !dot.is_finite() {
                            return Err(self.domain_error(
                                "point Gram entry must be finite",
                                location,
                                dot,
                            ));
                        }
                        let expected = if column == against { 1.0 } else { 0.0 };
                        if (dot - expected).abs() > MANIFOLD_MEMBERSHIP_TOL {
                            return Err(self.domain_error(
                                "point columns must be orthonormal",
                                location,
                                dot,
                            ));
                        }
                    }
                }
                Ok(())
            }
        }
    }

    fn validate_binding_domain(&self, x: &[f64], var: u32) -> Result<(), OptError> {
        match self.validate_point_domain(x) {
            Ok(()) => Ok(()),
            Err(OptError::RetractionDomain {
                manifold,
                what,
                location,
                measurement_bits,
            }) => Err(OptError::BindingDomain {
                var,
                manifold,
                what,
                location,
                measurement_bits,
            }),
            Err(other) => Err(other),
        }
    }

    fn validate_retraction_output(&self, output: Vec<f64>) -> Result<Vec<f64>, OptError> {
        Self::validate_retraction_finite("retraction output", &output)?;
        self.validate_point_domain(&output)?;
        Ok(output)
    }

    /// Retract: move from `x` along parameter vector `t`, landing ON
    /// the manifold. Rn: translation. Sphere: normalize(x+t). SO(3):
    /// right-multiply by `exp(ω/2)` (unit quaternion). Stiefel:
    /// Gram-Schmidt of `X+T` (QR retraction). Raw point and parameter
    /// storage must exactly match this manifold and contain finite
    /// components. The base point must already belong to the manifold;
    /// zero-norm and rank-deficient candidates are refused rather than
    /// normalized into fabricated points.
    ///
    /// # Errors
    /// [`OptError::ManifoldInvalid`], [`OptError::RetractionLen`],
    /// [`OptError::RetractionNonFinite`], or [`OptError::RetractionDomain`].
    #[must_use]
    #[allow(clippy::too_many_lines)] // one auditable branch per manifold and refusal phase
    pub fn retract(&self, x: &[f64], t: &[f64]) -> Result<Vec<f64>, OptError> {
        self.validate(&AdmissionCaps::default())?;
        let point_dim = self.point_dim().ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{self:?} has no representable point dimension"),
        })?;
        if x.len() as u64 != u64::from(point_dim) {
            return Err(OptError::RetractionLen {
                input: "retraction point",
                expected: point_dim,
                got: x.len() as u64,
            });
        }
        Self::validate_retraction_finite("retraction point", x)?;
        self.validate_point_domain(x)?;
        let param_dim = self.param_dim().ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{self:?} has no representable retraction parameter dimension"),
        })?;
        if t.len() as u64 != u64::from(param_dim) {
            return Err(OptError::RetractionLen {
                input: "retraction step",
                expected: param_dim,
                got: t.len() as u64,
            });
        }
        Self::validate_retraction_finite("retraction step", t)?;
        match *self {
            Manifold::Rn { .. } => {
                let output = x.iter().zip(t).map(|(a, b)| a + b).collect();
                self.validate_retraction_output(output)
            }
            Manifold::Sphere { .. } => {
                let y: Vec<f64> = x.iter().zip(t).map(|(a, b)| a + b).collect();
                Self::validate_retraction_finite("retraction candidate", &y)?;
                let norm_sq = y.iter().map(|value| value * value).sum::<f64>();
                if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                    return Err(self.domain_error(
                        "candidate norm squared must be finite and nonsingular",
                        None,
                        norm_sq,
                    ));
                }
                let norm = norm_sq.sqrt();
                self.validate_retraction_output(y.iter().map(|value| value / norm).collect())
            }
            Manifold::So3 => {
                let half = [t[0] * 0.5, t[1] * 0.5, t[2] * 0.5];
                let angle_sq = half[0] * half[0] + half[1] * half[1] + half[2] * half[2];
                if !angle_sq.is_finite() {
                    return Err(self.domain_error(
                        "step half-angle norm squared must be finite",
                        None,
                        angle_sq,
                    ));
                }
                let ang = angle_sq.sqrt();
                let (s, c) = if ang > 1e-12 {
                    (ang.sin() / ang, ang.cos())
                } else {
                    (1.0, 1.0)
                };
                let dq = [c, half[0] * s, half[1] * s, half[2] * s];
                let q = [x[0], x[1], x[2], x[3]];
                let mut out = [
                    q[0] * dq[0] - q[1] * dq[1] - q[2] * dq[2] - q[3] * dq[3],
                    q[0] * dq[1] + q[1] * dq[0] + q[2] * dq[3] - q[3] * dq[2],
                    q[0] * dq[2] - q[1] * dq[3] + q[2] * dq[0] + q[3] * dq[1],
                    q[0] * dq[3] + q[1] * dq[2] - q[2] * dq[1] + q[3] * dq[0],
                ];
                Self::validate_retraction_finite("retraction candidate", &out)?;
                let norm_sq = out.iter().map(|value| value * value).sum::<f64>();
                if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                    return Err(self.domain_error(
                        "candidate norm squared must be finite and nonsingular",
                        None,
                        norm_sq,
                    ));
                }
                let norm = norm_sq.sqrt();
                for v in &mut out {
                    *v /= norm;
                }
                self.validate_retraction_output(out.to_vec())
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let mut cols = Vec::with_capacity(p);
                for column in 0..p {
                    let candidate = (0..n)
                        .map(|row| x[column * n + row] + t[column * n + row])
                        .collect::<Vec<f64>>();
                    Self::validate_retraction_finite("retraction candidate", &candidate)?;
                    cols.push(candidate);
                }
                // Deterministic Gram-Schmidt (QR retraction).
                for j in 0..p {
                    for k in 0..j {
                        let d: f64 = (0..n).map(|i| cols[j][i] * cols[k][i]).sum();
                        if !d.is_finite() {
                            return Err(self.domain_error(
                                "candidate Gram projection must be finite",
                                Some((j as u32, k as u32)),
                                d,
                            ));
                        }
                        let prior = cols[k].clone();
                        for (cj, ck) in cols[j].iter_mut().zip(&prior) {
                            *cj -= d * ck;
                        }
                        Self::validate_retraction_finite("retraction candidate", &cols[j])?;
                    }
                    let norm_sq = cols[j].iter().map(|value| value * value).sum::<f64>();
                    if !norm_sq.is_finite() || norm_sq <= RETRACTION_MIN_NORM_SQ {
                        return Err(self.domain_error(
                            "candidate column is rank-deficient",
                            Some((j as u32, j as u32)),
                            norm_sq,
                        ));
                    }
                    let norm = norm_sq.sqrt();
                    for v in &mut cols[j] {
                        *v /= norm;
                    }
                }
                self.validate_retraction_output(cols.concat())
            }
        }
    }
}

/// Descent policy.
#[derive(Debug, Clone, Copy)]
pub struct DescentOptions {
    /// Maximum steps.
    pub steps: u32,
    /// Learning rate.
    pub lr: f64,
    /// Finite-difference step.
    pub fd_h: f64,
}

impl Default for DescentOptions {
    fn default() -> Self {
        DescentOptions {
            steps: 200,
            lr: 0.2,
            fd_h: 1e-6,
        }
    }
}

/// What the descent did.
#[derive(Debug, Clone)]
pub struct DescentReport {
    /// Final point.
    pub x: Vec<f64>,
    /// Initial objective value.
    pub f0: f64,
    /// Final objective value.
    pub f_final: f64,
    /// Objective evaluations spent.
    pub evals: u64,
    /// Steps taken.
    pub steps_taken: u32,
    /// True when the P4 budget stopped the run (a receipt, not an
    /// error — the point is still valid).
    pub budget_stopped: bool,
}

/// Toy Riemannian gradient descent of a closure over ONE manifold
/// variable: FD gradient through the retraction, fixed step. Proves
/// retraction metadata is consumable; polls cancellation each step and
/// honors `max_evals` (0 = unlimited). The manifold, start point
/// (length AND finite components), and step policy (`fd_h`/`lr`
/// finite, positive) are leaf-gated BEFORE any descent arithmetic.
///
/// # Errors
/// [`OptError::Cancelled`] / [`OptError::ManifoldInvalid`] /
/// [`OptError::BindingLen`] / [`OptError::NonFinite`] /
/// [`OptError::BadParam`].
pub fn descend_fn(
    manifold: Manifold,
    f: &dyn Fn(&[f64]) -> f64,
    x0: &[f64],
    opts: DescentOptions,
    max_evals: u64,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    let checked_f = |x: &[f64]| finite_descent_value(f(x), "descent objective result");
    descend_fn_checked(manifold, &checked_f, x0, opts, max_evals, cx)
}

fn finite_descent_value(value: f64, what: &'static str) -> Result<f64, OptError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OptError::NonFinite {
            what,
            bits: value.to_bits(),
        })
    }
}

fn descent_checkpoint(cx: &Cx<'_>) -> Result<(), OptError> {
    cx.checkpoint().map_err(|_| OptError::Cancelled)
}

fn descend_fn_checked(
    manifold: Manifold,
    f: &dyn Fn(&[f64]) -> Result<f64, OptError>,
    x0: &[f64],
    opts: DescentOptions,
    max_evals: u64,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    manifold.validate(&AdmissionCaps::default())?;
    let point_dim = manifold
        .point_dim()
        .ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{manifold:?} has no representable point dimension"),
        })?;
    if x0.len() as u64 != u64::from(point_dim) {
        return Err(OptError::BindingLen {
            var: 0,
            expected: point_dim,
            got: x0.len() as u64,
        });
    }
    // Leaf gating (review High #6, bead j3vb5): a non-finite start
    // point or degenerate step policy must refuse BEFORE any descent
    // arithmetic — NaN would otherwise propagate through retractions
    // and finite differences as plausible-looking garbage.
    for component in x0 {
        if !component.is_finite() {
            return Err(OptError::NonFinite {
                what: "descent initial point component",
                bits: component.to_bits(),
            });
        }
    }
    manifold.validate_point_domain(x0)?;
    if !(opts.fd_h.is_finite() && opts.fd_h > 0.0) {
        return Err(OptError::BadParam {
            what: "descent finite-difference step fd_h (finite, > 0)",
            value: opts.fd_h,
        });
    }
    if !(opts.lr.is_finite() && opts.lr > 0.0) {
        return Err(OptError::BadParam {
            what: "descent learning rate lr (finite, > 0; descent, not ascent)",
            value: opts.lr,
        });
    }
    let param_dim = manifold
        .param_dim()
        .ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{manifold:?} has no representable descent parameter dimension"),
        })?;
    let pd = usize::try_from(param_dim).map_err(|_| OptError::ManifoldInvalid {
        what: format!("{manifold:?} descent parameter dimension does not fit usize"),
    })?;
    let atomic_step_evals = u64::from(param_dim).saturating_mul(2).saturating_add(1);
    descent_checkpoint(cx)?;
    let mut x = x0.to_vec();
    let mut evals = 0u64;
    let mut budget_stopped = false;
    let f0 = f(&x)?;
    evals += 1;
    descent_checkpoint(cx)?;
    let mut steps_taken = 0;
    'outer: for _ in 0..opts.steps {
        descent_checkpoint(cx)?;
        // One landed step is atomic: reserve its complete central-FD
        // gradient (two probes per parameter) plus the terminal value
        // before spending the first probe. A one-short cap therefore
        // leaves both x and the evaluation ledger at the prior receipt.
        if max_evals > 0 && evals.saturating_add(atomic_step_evals) > max_evals {
            budget_stopped = true;
            break 'outer;
        }
        let mut g = vec![0.0; pd];
        for (i, gi) in g.iter_mut().enumerate() {
            descent_checkpoint(cx)?;
            let mut t = vec![0.0; pd];
            t[i] = opts.fd_h;
            descent_checkpoint(cx)?;
            let xp = manifold.retract(&x, &t)?;
            descent_checkpoint(cx)?;
            let fp = f(&xp)?;
            evals += 1;
            descent_checkpoint(cx)?;
            t[i] = -opts.fd_h;
            let xm = manifold.retract(&x, &t)?;
            descent_checkpoint(cx)?;
            let fm = f(&xm)?;
            evals += 1;
            descent_checkpoint(cx)?;
            *gi = finite_descent_value(
                (fp - fm) / (2.0 * opts.fd_h),
                "finite-difference gradient component",
            )?;
        }
        descent_checkpoint(cx)?;
        let step: Vec<f64> = g.iter().map(|v| -opts.lr * v).collect();
        descent_checkpoint(cx)?;
        x = manifold.retract(&x, &step)?;
        descent_checkpoint(cx)?;
        steps_taken += 1;
    }
    let f_final = if steps_taken == 0 {
        f0
    } else {
        descent_checkpoint(cx)?;
        let value = f(&x)?;
        evals += 1;
        descent_checkpoint(cx)?;
        value
    };
    descent_checkpoint(cx)?;
    Ok(DescentReport {
        x,
        f0,
        f_final,
        evals,
        steps_taken,
        budget_stopped,
    })
}

/// Toy descent of a problem's FIRST objective over its FIRST variable
/// (the IR-driven variant; enforces `problem.budget` per P4).
///
/// # Errors
/// [`OptError::Cancelled`] / evaluation teaching errors.
pub fn descend_ir(
    problem: &Problem,
    x0: &[f64],
    opts: DescentOptions,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    // A problem with no objective or no variable is unsolvable — return a
    // structured error, never an index panic (`ProblemBuilder` does not
    // require either, and `descend_ir` is public).
    let obj = *problem
        .objectives
        .first()
        .ok_or(OptError::IndexOut { index: 0, len: 0 })?;
    let sign = match obj.sense {
        crate::ir::Sense::Minimize => 1.0,
        crate::ir::Sense::Maximize => -1.0,
    };
    let manifold = problem
        .vars
        .first()
        .ok_or(OptError::IndexOut { index: 0, len: 0 })?
        .manifold;
    // The fallible objective runs only inside the shared descent seam,
    // after manifold/start/options leaf validation. Its first successful
    // invocation is f0 and is counted exactly once, including under a
    // one-evaluation budget.
    let f = |x: &[f64]| -> Result<f64, OptError> {
        let scalar = eval(problem, obj.node, &[x.to_vec()])?
            .scalar()
            .ok_or(OptError::NotScalar { node: obj.node.0 })?;
        finite_descent_value(sign * obj.weight * scalar, "weighted IR objective")
    };
    descend_fn_checked(manifold, &f, x0, opts, problem.budget.max_evals, cx)
}

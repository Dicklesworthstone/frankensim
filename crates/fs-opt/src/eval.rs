//! Evaluation + the TOY Riemannian descent that proves the manifold
//! metadata is consumable: retractions keep iterates ON their
//! manifolds ("optimize an orientation" never becomes "optimize 9
//! numbers and renormalize when it explodes"). Gradients here are
//! finite differences through retractions — a deliberately simple
//! consumer; exact adjoints are the gradient-stack bead's.

use crate::ir::{Expr, Manifold, NodeId, OptError, Problem};
use fs_exec::Cx;

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
/// [`OptError::BindingLen`].
pub fn eval(problem: &Problem, node: NodeId, bindings: &[Vec<f64>]) -> Result<Value, OptError> {
    if node.0 as usize >= problem.exprs.len() {
        return Err(OptError::UnknownNode { id: node.0 });
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
    }
    let mut memo: Vec<Option<Value>> = vec![None; problem.exprs.len()];
    eval_at(problem, node, bindings, &mut memo)
}

#[allow(clippy::too_many_lines)] // one arm per node kind: the evaluator IS the semantics
fn eval_at(
    problem: &Problem,
    node: NodeId,
    bindings: &[Vec<f64>],
    memo: &mut Vec<Option<Value>>,
) -> Result<Value, OptError> {
    if let Some(v) = &memo[node.0 as usize] {
        return Ok(v.clone());
    }
    let ev = |n: NodeId, memo: &mut Vec<Option<Value>>| eval_at(problem, n, bindings, memo);
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
            let v = ev(*of, memo)?;
            match v {
                Value::V(xs) => Value::S(xs[*index as usize]),
                Value::S(_) => unreachable!("builder enforced vector shape"),
            }
        }
        Expr::Const { value, .. } => Value::S(*value),
        Expr::Add(a, b) => match (ev(*a, memo)?, ev(*b, memo)?) {
            (Value::S(x), Value::S(y)) => Value::S(x + y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p + q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Sub(a, b) => match (ev(*a, memo)?, ev(*b, memo)?) {
            (Value::S(x), Value::S(y)) => Value::S(x - y),
            (Value::V(x), Value::V(y)) => Value::V(x.iter().zip(&y).map(|(p, q)| p - q).collect()),
            _ => unreachable!("builder enforced matching shapes"),
        },
        Expr::Mul(a, b) => match (ev(*a, memo)?, ev(*b, memo)?) {
            (Value::S(x), Value::S(y)) => Value::S(x * y),
            (Value::S(s), Value::V(v)) | (Value::V(v), Value::S(s)) => {
                Value::V(v.iter().map(|p| p * s).collect())
            }
            _ => unreachable!("builder rejected vector*vector"),
        },
        Expr::Div(a, b) => {
            let (x, y) = (scalar(ev(*a, memo)?), scalar(ev(*b, memo)?));
            Value::S(x / y)
        }
        Expr::Neg(a) => match ev(*a, memo)? {
            Value::S(x) => Value::S(-x),
            Value::V(v) => Value::V(v.iter().map(|p| -p).collect()),
        },
        Expr::Powi { base, exp } => Value::S(fs_math::det::powi(scalar(ev(*base, memo)?), *exp)),
        Expr::Sqrt(a) => Value::S(fs_math::det::sqrt(scalar(ev(*a, memo)?))),
        Expr::Exp(a) => Value::S(fs_math::det::exp(scalar(ev(*a, memo)?))),
        Expr::Ln(a) => Value::S(fs_math::det::ln(scalar(ev(*a, memo)?))),
        Expr::Tanh(a) => Value::S(fs_math::det::tanh(scalar(ev(*a, memo)?))),
        Expr::Dot(a, b) => match (ev(*a, memo)?, ev(*b, memo)?) {
            (Value::V(x), Value::V(y)) => Value::S(x.iter().zip(&y).map(|(p, q)| p * q).sum()),
            _ => unreachable!("builder enforced vectors"),
        },
        Expr::NormSq(a) => match ev(*a, memo)? {
            Value::V(x) => Value::S(x.iter().map(|p| p * p).sum()),
            Value::S(_) => unreachable!("builder enforced vector"),
        },
        Expr::Min(a, b) => Value::S(scalar(ev(*a, memo)?).min(scalar(ev(*b, memo)?))),
        Expr::Max(a, b) => Value::S(scalar(ev(*a, memo)?).max(scalar(ev(*b, memo)?))),
        Expr::Abs(a) => Value::S(scalar(ev(*a, memo)?).abs()),
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
    memo[node.0 as usize] = Some(out.clone());
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

    /// Retract: move from `x` along parameter vector `t`, landing ON
    /// the manifold. Rn: translation. Sphere: normalize(x+t). SO(3):
    /// right-multiply by `exp(ω/2)` (unit quaternion). Stiefel:
    /// Gram-Schmidt of `X+T` (QR retraction).
    #[must_use]
    pub fn retract(&self, x: &[f64], t: &[f64]) -> Vec<f64> {
        match *self {
            Manifold::Rn { .. } => x.iter().zip(t).map(|(a, b)| a + b).collect(),
            Manifold::Sphere { .. } => {
                let y: Vec<f64> = x.iter().zip(t).map(|(a, b)| a + b).collect();
                let n = y.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
                y.iter().map(|v| v / n).collect()
            }
            Manifold::So3 => {
                let half = [t[0] * 0.5, t[1] * 0.5, t[2] * 0.5];
                let ang = (half[0] * half[0] + half[1] * half[1] + half[2] * half[2]).sqrt();
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
                let n = out.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
                for v in &mut out {
                    *v /= n;
                }
                out.to_vec()
            }
            Manifold::Stiefel { n, p } => {
                let (n, p) = (n as usize, p as usize);
                let mut cols: Vec<Vec<f64>> = (0..p)
                    .map(|j| {
                        (0..n)
                            .map(|i| x[j * n + i] + t[j * n + i])
                            .collect::<Vec<f64>>()
                    })
                    .collect();
                // Deterministic Gram-Schmidt (QR retraction).
                for j in 0..p {
                    for k in 0..j {
                        let d: f64 = (0..n).map(|i| cols[j][i] * cols[k][i]).sum();
                        let prior = cols[k].clone();
                        for (cj, ck) in cols[j].iter_mut().zip(&prior) {
                            *cj -= d * ck;
                        }
                    }
                    let nn = cols[j]
                        .iter()
                        .map(|v| v * v)
                        .sum::<f64>()
                        .sqrt()
                        .max(1e-300);
                    for v in &mut cols[j] {
                        *v /= nn;
                    }
                }
                cols.concat()
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
/// honors `max_evals` (0 = unlimited).
///
/// # Errors
/// [`OptError::Cancelled`].
pub fn descend_fn(
    manifold: Manifold,
    f: &dyn Fn(&[f64]) -> f64,
    x0: &[f64],
    opts: DescentOptions,
    max_evals: u64,
    cx: &Cx<'_>,
) -> Result<DescentReport, OptError> {
    let mut x = x0.to_vec();
    let mut evals = 0u64;
    let mut budget_stopped = false;
    let f0 = f(&x);
    evals += 1;
    let pd = manifold
        .param_dim()
        .and_then(|d| usize::try_from(d).ok())
        .ok_or_else(|| OptError::ManifoldInvalid {
            what: format!("{manifold:?} has no representable descent parameter dimension"),
        })?;
    let mut steps_taken = 0;
    'outer: for _ in 0..opts.steps {
        if cx.checkpoint().is_err() {
            return Err(OptError::Cancelled);
        }
        let mut g = vec![0.0; pd];
        for (i, gi) in g.iter_mut().enumerate() {
            if max_evals > 0 && evals + 2 > max_evals {
                budget_stopped = true;
                break 'outer;
            }
            let mut t = vec![0.0; pd];
            t[i] = opts.fd_h;
            let fp = f(&manifold.retract(&x, &t));
            t[i] = -opts.fd_h;
            let fm = f(&manifold.retract(&x, &t));
            evals += 2;
            *gi = (fp - fm) / (2.0 * opts.fd_h);
        }
        let step: Vec<f64> = g.iter().map(|v| -opts.lr * v).collect();
        x = manifold.retract(&x, &step);
        steps_taken += 1;
    }
    let f_final = f(&x);
    evals += 1;
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
    // Surface evaluation errors (PDE/stochastic nodes) BEFORE descending.
    eval(problem, obj.node, &[x0.to_vec()])?;
    let manifold = problem
        .vars
        .first()
        .ok_or(OptError::IndexOut { index: 0, len: 0 })?
        .manifold;
    let f = |x: &[f64]| -> f64 {
        let v = eval(problem, obj.node, &[x.to_vec()]).expect("checked evaluable above");
        sign * v.scalar().expect("objective roots are scalar")
    };
    descend_fn(manifold, &f, x0, opts, problem.budget.max_evals, cx)
}

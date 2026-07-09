//! The Problem-IR STUDY RUNNER (bead ijil): drive an `fs_opt::Problem`
//! end-to-end through the fs-ascent engines — the seam where "problems
//! are data" meets "optimizers are engines".
//!
//! - VARIABLE PACKING: all declared variables concatenate into one
//!   flat vector across the manifold product; per-block Riemannian
//!   handling (tangent projection + retraction from
//!   [`crate::riemann`]) keeps orientations orientations.
//! - GRADIENTS: central differences through `fs_opt::eval` in the
//!   TANGENT parameterization (the live IR is evaluation-only — no
//!   reverse mode; documented, deterministic fixed h).
//! - CONSTRAINTS: `EqZero`/`LeZero` roots route to the constrained
//!   engines (AL always; IP/SQP by option) through packed adapters.
//! - BUDGET: `Problem.budget.max_evals` threads into the stop algebra
//!   (`StopRule::Budget`) alongside the caller's rules.
//! - RESUMABLE: [`Study`] is the checkpoint (clone = checkpoint); a
//!   split run is BITWISE equal to a straight run (house pattern,
//!   gated).

use crate::riemann::{retract, tangent_project};
use crate::stop::{StopObservation, StopReason, StopRule};
use fs_opt::{ConstraintKind, Manifold, Problem, Sense, eval};

/// One packed variable block.
#[derive(Debug, Clone)]
struct Block {
    manifold: Manifold,
    /// Offset in the packed point vector.
    start: usize,
    /// Point storage length.
    len: usize,
}

/// The packed view of a problem's variables.
#[derive(Debug, Clone)]
pub struct Packing {
    blocks: Vec<Block>,
    /// Total packed length.
    pub dim: usize,
}

impl Packing {
    /// Build the packing for a problem's variable list.
    #[must_use]
    pub fn new(problem: &Problem) -> Packing {
        let mut blocks = Vec::new();
        let mut start = 0usize;
        for v in &problem.vars {
            let len = v.manifold.point_dim() as usize;
            blocks.push(Block {
                manifold: v.manifold,
                start,
                len,
            });
            start += len;
        }
        Packing { blocks, dim: start }
    }

    /// Split a packed point into per-variable bindings for `eval`.
    #[must_use]
    pub fn unpack(&self, x: &[f64]) -> Vec<Vec<f64>> {
        self.blocks
            .iter()
            .map(|b| x[b.start..b.start + b.len].to_vec())
            .collect()
    }

    /// Per-block tangent projection of an ambient direction.
    #[must_use]
    pub fn project(&self, x: &[f64], g: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.dim];
        for b in &self.blocks {
            let seg = tangent_project(
                &b.manifold,
                &x[b.start..b.start + b.len],
                &g[b.start..b.start + b.len],
            );
            out[b.start..b.start + b.len].copy_from_slice(&seg);
        }
        out
    }

    /// Per-block retraction of a step.
    #[must_use]
    pub fn retract(&self, x: &[f64], step: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.dim];
        for b in &self.blocks {
            let seg = retract(
                &b.manifold,
                &x[b.start..b.start + b.len],
                &step[b.start..b.start + b.len],
            );
            out[b.start..b.start + b.len].copy_from_slice(&seg);
        }
        out
    }
}

/// A resumable study over an IR problem: state clone = checkpoint.
#[derive(Debug, Clone)]
pub struct Study {
    packing: Packing,
    /// Current packed iterate.
    pub x: Vec<f64>,
    /// Objective history (most recent last).
    pub history: Vec<f64>,
    /// Evaluations spent (budget accounting).
    pub evals: usize,
    /// Steps taken.
    pub steps: usize,
    fd_h: f64,
    lr: f64,
}

/// Outcome of a study segment.
#[derive(Debug, Clone)]
pub struct StudyReport {
    /// Final objective.
    pub f: f64,
    /// Final ‖projected gradient‖∞.
    pub grad_norm: f64,
    /// Why the segment stopped.
    pub reason: StopReason,
    /// Evaluations spent so far (cumulative).
    pub evals: usize,
}

/// The weighted, sense-corrected total objective of a problem at
/// packed point `x` (all objectives, not just the first).
fn objective(problem: &Problem, packing: &Packing, x: &[f64], evals: &mut usize) -> f64 {
    let bindings = packing.unpack(x);
    let mut total = 0.0f64;
    for o in &problem.objectives {
        let v = eval(problem, o.node, &bindings)
            .expect("objective evaluable (checked at study start)")
            .scalar()
            .expect("objective roots are scalar");
        let sign = match o.sense {
            Sense::Minimize => 1.0,
            Sense::Maximize => -1.0,
        };
        total += sign * o.weight * v;
    }
    *evals += 1;
    total
}

impl Study {
    /// Start a study at `x0` (packed). Verifies evaluability of every
    /// objective and constraint root up front (fail loud, fail early).
    ///
    /// # Panics
    /// If `x0` has the wrong packed length or a root is unevaluable.
    #[must_use]
    pub fn new(problem: &Problem, x0: &[f64], fd_h: f64, lr: f64) -> Study {
        let packing = Packing::new(problem);
        assert_eq!(x0.len(), packing.dim, "packed x0 length mismatch");
        let bindings = packing.unpack(x0);
        for o in &problem.objectives {
            let _ = eval(problem, o.node, &bindings).expect("objective root must evaluate");
        }
        for c in &problem.constraints {
            let _ = eval(problem, c.node, &bindings).expect("constraint root must evaluate");
        }
        Study {
            packing,
            x: x0.to_vec(),
            history: Vec::new(),
            evals: 0,
            steps: 0,
            fd_h,
            lr,
        }
    }

    /// Central-difference TANGENT gradient of the total objective.
    fn tangent_gradient(&mut self, problem: &Problem) -> Vec<f64> {
        let n = self.packing.dim;
        let mut g = vec![0.0f64; n];
        for i in 0..n {
            let mut t = vec![0.0f64; n];
            t[i] = self.fd_h;
            let xp = self.packing.retract(&self.x, &t);
            t[i] = -self.fd_h;
            let xm = self.packing.retract(&self.x, &t);
            let fp = objective(problem, &self.packing, &xp, &mut self.evals);
            let fm = objective(problem, &self.packing, &xm, &mut self.evals);
            g[i] = (fp - fm) / (2.0 * self.fd_h);
        }
        self.packing.project(&self.x, &g)
    }

    /// Run a segment: projected-gradient steps with retraction until a
    /// stop rule fires, the problem budget runs out, or `max_steps`.
    /// The problem's own `budget.max_evals` is ALWAYS added to the
    /// caller's rules (P4: budgets are not optional).
    pub fn run(&mut self, problem: &Problem, rule: &StopRule, max_steps: usize) -> StudyReport {
        let mut rules = vec![rule.clone()];
        if problem.budget.max_evals > 0 {
            rules.push(StopRule::Budget(
                usize::try_from(problem.budget.max_evals).unwrap_or(usize::MAX),
            ));
        }
        let combined = StopRule::Any(rules);
        let mut reason = StopReason::IterationCap;
        let mut f = objective(problem, &self.packing, &self.x, &mut self.evals);
        let mut gnorm = f64::INFINITY;
        for _ in 0..max_steps {
            let g = self.tangent_gradient(problem);
            gnorm = g.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            self.history.push(f);
            let obs = StopObservation {
                grad_norm: gnorm,
                objective: f,
                evals: self.evals,
                history: &self.history,
            };
            if let Some(r) = combined.check(&obs) {
                reason = r;
                break;
            }
            let step: Vec<f64> = g.iter().map(|v| -self.lr * v).collect();
            self.x = self.packing.retract(&self.x, &step);
            self.steps += 1;
            f = objective(problem, &self.packing, &self.x, &mut self.evals);
        }
        StudyReport {
            f,
            grad_norm: gnorm,
            reason,
            evals: self.evals,
        }
    }

    /// Constrained adapters: expose the problem's `EqZero`/`LeZero`
    /// roots as packed callbacks for [`crate::augmented_lagrangian`],
    /// [`crate::interior_point`] or [`crate::sqp`] — with FD
    /// Jacobian-transpose actions in the packed coordinates. Returns
    /// (ce, ce_jt, ci, ci_jt) closures over the problem.
    #[must_use]
    #[allow(clippy::type_complexity)]
    pub fn constraint_adapters<'p>(
        problem: &'p Problem,
        packing: &'p Packing,
        fd_h: f64,
    ) -> (
        impl Fn(&[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64], &[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64], &[f64]) -> Vec<f64> + 'p,
    ) {
        let eval_kind = move |x: &[f64], kind: ConstraintKind| -> Vec<f64> {
            let bindings = packing.unpack(x);
            problem
                .constraints
                .iter()
                .filter(|c| c.kind == kind)
                .map(|c| {
                    eval(problem, c.node, &bindings)
                        .expect("constraint evaluable")
                        .scalar()
                        .expect("constraint roots are scalar")
                })
                .collect()
        };
        let jt_kind = move |x: &[f64], w: &[f64], kind: ConstraintKind| -> Vec<f64> {
            // FD of wᵀc(x) — one directional pass per packed dim.
            let n = x.len();
            let mut out = vec![0.0f64; n];
            for i in 0..n {
                let mut xp = x.to_vec();
                xp[i] += fd_h;
                let mut xm = x.to_vec();
                xm[i] -= fd_h;
                let cp = eval_kind(&xp, kind);
                let cm = eval_kind(&xm, kind);
                out[i] = cp
                    .iter()
                    .zip(&cm)
                    .zip(w)
                    .map(|((p, m2), wi)| wi * (p - m2) / (2.0 * fd_h))
                    .sum();
            }
            out
        };
        (
            move |x: &[f64]| eval_kind(x, ConstraintKind::EqZero),
            move |x: &[f64], w: &[f64]| jt_kind(x, w, ConstraintKind::EqZero),
            move |x: &[f64]| eval_kind(x, ConstraintKind::LeZero),
            move |x: &[f64], w: &[f64]| jt_kind(x, w, ConstraintKind::LeZero),
        )
    }
}

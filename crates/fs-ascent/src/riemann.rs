//! Riemannian L-BFGS over fs-opt's manifold METADATA: ambient
//! gradient → tangent projection → search along RETRACTED curves,
//! with memory pairs transported by projection (the standard
//! projection-based vector transport for embedded manifolds).
//! fs-opt's `Manifold` is metadata (dimensions, kinds); the
//! OPERATIONS (projection, retraction) live here with the engines
//! that consume them. Manifold invariants (unit norms) hold to
//! roundoff along the whole iterate path — tested, not assumed.

use crate::stop::{StopObservation, StopReason, StopRule};
use fs_opt::Manifold;
use std::collections::VecDeque;

/// Project an ambient vector onto the tangent space at `x`.
///
/// # Panics
/// For metadata-only manifolds (Stiefel — its retraction lands with
/// its consumer bead; using it in a descent today is a modeling
/// error, surfaced loudly).
#[must_use]
pub fn tangent_project(man: &Manifold, x: &[f64], v: &[f64]) -> Vec<f64> {
    match man {
        Manifold::Rn { .. } => v.to_vec(),
        Manifold::Sphere { .. } | Manifold::So3 => {
            // v − (v·x)·x (x assumed unit).
            let dot: f64 = v.iter().zip(x).map(|(a, b)| a * b).sum();
            v.iter()
                .zip(x)
                .map(|(vi, xi)| (-dot).mul_add(*xi, *vi))
                .collect()
        }
        Manifold::Stiefel { .. } => {
            panic!("Stiefel is metadata-only until its consumer bead supplies the retraction")
        }
    }
}

/// Retract from `x` along `step` (metric-projection retractions).
///
/// # Panics
/// Same metadata-only policy as [`tangent_project`].
#[must_use]
pub fn retract(man: &Manifold, x: &[f64], step: &[f64]) -> Vec<f64> {
    match man {
        Manifold::Rn { .. } => x.iter().zip(step).map(|(a, b)| a + b).collect(),
        Manifold::Sphere { .. } | Manifold::So3 => {
            let moved: Vec<f64> = x.iter().zip(step).map(|(a, b)| a + b).collect();
            let nrm = fs_math::det::sqrt(moved.iter().map(|t| t * t).sum());
            moved.iter().map(|t| t / nrm).collect()
        }
        Manifold::Stiefel { .. } => {
            panic!("Stiefel is metadata-only until its consumer bead supplies the retraction")
        }
    }
}

/// Outcome of a Riemannian L-BFGS run.
#[derive(Debug, Clone)]
pub struct RiemannianReport {
    /// Why the run stopped.
    pub reason: StopReason,
    /// Final Riemannian gradient norm (∞).
    pub grad_norm: f64,
    /// Final objective.
    pub f: f64,
    /// Iterations.
    pub iters: usize,
    /// Evaluations.
    pub evals: usize,
    /// Worst manifold-constraint violation observed along the path
    /// (the invariant certificate; ~1e−15 for sphere/SO3).
    pub worst_violation: f64,
}

/// Resumable Riemannian L-BFGS on a single-manifold variable.
#[derive(Debug, Clone)]
pub struct RiemannianLbfgs {
    /// The manifold.
    pub manifold: Manifold,
    /// Current point (on the manifold).
    pub x: Vec<f64>,
    /// Current objective.
    pub f: f64,
    /// Current RIEMANNIAN gradient (tangent).
    pub g: Vec<f64>,
    pairs: VecDeque<(Vec<f64>, Vec<f64>, f64)>,
    memory: usize,
    /// Iterations performed.
    pub iters: usize,
    /// Evaluations spent.
    pub evals: usize,
    /// Objective history.
    pub history: Vec<f64>,
    worst_violation: f64,
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0f64, f64::max)
}

fn violation(man: &Manifold, x: &[f64]) -> f64 {
    match man {
        Manifold::Sphere { .. } | Manifold::So3 => {
            let n2: f64 = x.iter().map(|v| v * v).sum();
            (fs_math::det::sqrt(n2) - 1.0).abs()
        }
        _ => 0.0,
    }
}

impl RiemannianLbfgs {
    /// Start at a point ON the manifold.
    #[must_use]
    pub fn new(manifold: Manifold, x0: &[f64], memory: usize, fg: crate::FnGrad<'_>) -> Self {
        let (f, g_amb) = fg(x0);
        let g = tangent_project(&manifold, x0, &g_amb);
        let worst = violation(&manifold, x0);
        RiemannianLbfgs {
            manifold,
            x: x0.to_vec(),
            f,
            g,
            pairs: VecDeque::new(),
            memory,
            iters: 0,
            evals: 1,
            history: vec![f],
            worst_violation: worst,
        }
    }

    fn direction(&self) -> Vec<f64> {
        let mut q = self.g.clone();
        let mut alphas = Vec::with_capacity(self.pairs.len());
        for (s, y, rho) in self.pairs.iter().rev() {
            let a = rho * s.iter().zip(&q).map(|(si, qi)| si * qi).sum::<f64>();
            for (qi, yi) in q.iter_mut().zip(y) {
                *qi = a.mul_add(-yi, *qi);
            }
            alphas.push(a);
        }
        if let Some((s, y, _)) = self.pairs.back() {
            let sy: f64 = s.iter().zip(y).map(|(a, b)| a * b).sum();
            let yy: f64 = y.iter().map(|v| v * v).sum();
            for qi in &mut q {
                *qi *= sy / yy;
            }
        }
        for ((s, y, rho), a) in self.pairs.iter().zip(alphas.iter().rev()) {
            let b = rho * y.iter().zip(&q).map(|(yi, qi)| yi * qi).sum::<f64>();
            let coeff = a - b;
            for (qi, si) in q.iter_mut().zip(s) {
                *qi = coeff.mul_add(*si, *qi);
            }
        }
        for qi in &mut q {
            *qi = -*qi;
        }
        q
    }

    /// Run with backtracking-Armijo along retracted curves (curve
    /// search; strong Wolfe on manifolds needs transported derivative
    /// bookkeeping — recorded follow-up, Armijo suffices for the
    /// convergence gates here).
    pub fn run(
        &mut self,
        fg: crate::FnGrad<'_>,
        rule: &StopRule,
        max_iters: usize,
    ) -> RiemannianReport {
        let mut reason = StopReason::IterationCap;
        for _ in 0..max_iters {
            let obs = StopObservation {
                grad_norm: inf_norm(&self.g),
                objective: self.f,
                evals: self.evals,
                history: &self.history,
            };
            if let Some(r) = rule.check(&obs) {
                reason = r;
                break;
            }
            let mut d = self.direction();
            let mut dphi0: f64 = d.iter().zip(&self.g).map(|(a, b)| a * b).sum();
            if dphi0 >= 0.0 {
                d = self.g.iter().map(|g| -g).collect();
                dphi0 = -self.g.iter().map(|g| g * g).sum::<f64>();
            }
            // Backtracking Armijo along the retraction.
            let mut alpha = 1.0f64;
            let mut accepted = None;
            for _ in 0..40 {
                let step: Vec<f64> = d.iter().map(|di| alpha * di).collect();
                let x_new = retract(&self.manifold, &self.x, &step);
                let (f_new, g_amb) = fg(&x_new);
                self.evals += 1;
                if f_new <= (1e-4 * alpha).mul_add(dphi0, self.f) {
                    accepted = Some((x_new, f_new, g_amb));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((x_new, f_new, g_amb)) = accepted else {
                reason = StopReason::Stall;
                break;
            };
            self.worst_violation = self.worst_violation.max(violation(&self.manifold, &x_new));
            let g_new = tangent_project(&self.manifold, &x_new, &g_amb);
            // Transport old quantities to the new tangent space by
            // projection, then form the curvature pair.
            let g_old_t = tangent_project(&self.manifold, &x_new, &self.g);
            let s_raw: Vec<f64> = d.iter().map(|di| alpha * di).collect();
            let s = tangent_project(&self.manifold, &x_new, &s_raw);
            let y: Vec<f64> = g_new.iter().zip(&g_old_t).map(|(a, b)| a - b).collect();
            let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
            if sy > 1e-14 {
                if self.pairs.len() == self.memory {
                    self.pairs.pop_front();
                }
                // Transport EXISTING pairs by projection too.
                for (ps, py, _) in &mut self.pairs {
                    *ps = tangent_project(&self.manifold, &x_new, ps);
                    *py = tangent_project(&self.manifold, &x_new, py);
                }
                self.pairs.push_back((s, y, 1.0 / sy));
            }
            self.x = x_new;
            self.f = f_new;
            self.g = g_new;
            self.iters += 1;
            self.history.push(self.f);
        }
        RiemannianReport {
            reason,
            grad_norm: inf_norm(&self.g),
            f: self.f,
            iters: self.iters,
            evals: self.evals,
            worst_violation: self.worst_violation,
        }
    }
}

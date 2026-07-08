//! L-BFGS — THE default first-order engine: limited-memory two-loop
//! recursion (m ≈ 17-class), strong-Wolfe line search, RESUMABLE
//! state (clone = checkpoint; split runs bitwise-equal to straight
//! runs — the house pattern), stopping through the condition algebra,
//! and gradient-norm certificates in every report.

use crate::stop::{StopObservation, StopReason, StopRule};
use crate::wolfe::strong_wolfe;
use std::collections::VecDeque;

/// Resumable L-BFGS state. Plain data: `clone()` is a checkpoint.
#[derive(Debug, Clone)]
pub struct LbfgsState {
    /// Current iterate.
    pub x: Vec<f64>,
    /// Current objective.
    pub f: f64,
    /// Current gradient.
    pub g: Vec<f64>,
    /// Memory pairs (s, y, 1/yᵀs), oldest first.
    pairs: VecDeque<(Vec<f64>, Vec<f64>, f64)>,
    /// Memory length m.
    pub memory: usize,
    /// Iterations performed.
    pub iters: usize,
    /// Function+gradient evaluations spent.
    pub evals: usize,
    /// Objective history (per accepted iterate).
    pub history: Vec<f64>,
}

/// Outcome of an L-BFGS run.
#[derive(Debug, Clone)]
pub struct LbfgsReport {
    /// Why the run stopped.
    pub reason: StopReason,
    /// Final ‖g‖∞ (the unconstrained certificate).
    pub grad_norm: f64,
    /// Final objective.
    pub f: f64,
    /// Iterations performed (cumulative across resumes).
    pub iters: usize,
    /// Evaluations spent (cumulative).
    pub evals: usize,
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0f64, f64::max)
}

impl LbfgsState {
    /// Start at `x0` with memory `m` (one evaluation is spent here).
    #[must_use]
    pub fn new(x0: &[f64], memory: usize, fg: crate::FnGrad<'_>) -> Self {
        let (f, g) = fg(x0);
        LbfgsState {
            x: x0.to_vec(),
            f,
            g,
            pairs: VecDeque::new(),
            memory,
            iters: 0,
            evals: 1,
            history: vec![f],
        }
    }

    /// Two-loop recursion: d = −H·g from the memory pairs.
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
        // Initial scaling γ = sᵀy/yᵀy from the most recent pair.
        if let Some((s, y, _)) = self.pairs.back() {
            let sy: f64 = s.iter().zip(y).map(|(a, b)| a * b).sum();
            let yy: f64 = y.iter().map(|v| v * v).sum();
            let gamma = sy / yy;
            for qi in &mut q {
                *qi *= gamma;
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

    /// Run until the stop rule fires or `max_iters` ADDITIONAL
    /// iterations; resumable (call again to continue bitwise).
    pub fn run(&mut self, fg: crate::FnGrad<'_>, rule: &StopRule, max_iters: usize) -> LbfgsReport {
        let n = self.x.len();
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
            let mut dphi0: f64 = d.iter().zip(&self.g).map(|(di, gi)| di * gi).sum();
            if dphi0 >= 0.0 {
                // Memory produced a non-descent direction (can happen
                // right after a resume with stale curvature on
                // nonconvex terrain): fall back to steepest descent.
                d = self.g.iter().map(|g| -g).collect();
                dphi0 = -self.g.iter().map(|g| g * g).sum::<f64>();
            }
            let x0 = self.x.clone();
            let f0 = self.f;
            let mut evals_here = 0usize;
            let mut last: Option<(f64, Vec<f64>, Vec<f64>)> = None;
            let outcome = {
                let mut phi = |alpha: f64| -> (f64, f64) {
                    let xt: Vec<f64> = x0
                        .iter()
                        .zip(&d)
                        .map(|(xi, di)| alpha.mul_add(*di, *xi))
                        .collect();
                    let (f, g) = fg(&xt);
                    evals_here += 1;
                    let dphi: f64 = g.iter().zip(&d).map(|(gi, di)| gi * di).sum();
                    last = Some((f, g, xt));
                    (f, dphi)
                };
                strong_wolfe(&mut phi, f0, dphi0, 1.0, 1e-4, 0.9)
            };
            self.evals += evals_here;
            if !outcome.success {
                reason = StopReason::Stall;
                break;
            }
            let (f_new, g_new, x_new) = last.expect("line search evaluated at least once");
            // Curvature pair.
            let s: Vec<f64> = x_new.iter().zip(&x0).map(|(a, b)| a - b).collect();
            let y: Vec<f64> = g_new.iter().zip(&self.g).map(|(a, b)| a - b).collect();
            let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
            if sy > 1e-14 * inf_norm(&s) * inf_norm(&y).max(1e-30) {
                if self.pairs.len() == self.memory {
                    self.pairs.pop_front();
                }
                self.pairs.push_back((s, y, 1.0 / sy));
            }
            self.x = x_new;
            self.f = f_new;
            self.g = g_new;
            self.iters += 1;
            self.history.push(self.f);
            let _ = n;
        }
        LbfgsReport {
            reason,
            grad_norm: inf_norm(&self.g),
            f: self.f,
            iters: self.iters,
            evals: self.evals,
        }
    }
}

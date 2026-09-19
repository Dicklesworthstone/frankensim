//! Trust-region Newton–Krylov: Steihaug-CG on the quadratic model
//! with NEGATIVE-CURVATURE handling (follow the direction to the
//! boundary — the feature that separates TR from line-search Newton
//! on nonconvex terrain), classical radius update laws (G0-tested),
//! and matrix-free Hessian-vector products via caller-supplied
//! closures. Second-order adjoints are recorded follow-up; the
//! finite-difference-of-gradients Hv helper carries its O(√ε)
//! accuracy in its name rather than hiding it.

use crate::stop::{StopObservation, StopReason, StopRule};
use fs_exec::Cx;

/// Outcome of a trust-region run.
#[derive(Debug, Clone)]
pub struct TrustRegionReport {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Final objective.
    pub f: f64,
    /// Final ‖g‖∞.
    pub grad_norm: f64,
    /// Outer iterations.
    pub iters: usize,
    /// Function+gradient evaluations.
    pub evals: usize,
    /// Hessian-vector products spent.
    pub hv_evals: usize,
    /// Steps that hit the boundary via negative curvature.
    pub negative_curvature_hits: usize,
}

/// Why control returned to the caller. A cancellation pause is not convergence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrustRegionProgress {
    /// Cancellation was observed at a complete outer-iteration boundary.
    Paused,
    /// A numerical stopping rule, resource limit, or iteration cap fired.
    Stopped(StopReason),
}

/// Attributed outcome of a resumable run, with cumulative solver accounting.
#[derive(Debug, Clone)]
pub struct TrustRegionRunReport {
    /// Pause or numerical/resource stop attribution.
    pub progress: TrustRegionProgress,
    /// Last accepted point and cumulative callback counts.
    pub solution: TrustRegionReport,
}

/// Resumable Newton–Krylov state. `clone()` retains the radius as well as the
/// accepted point, gradient, history, and counters; resumption performs no
/// extra initial evaluation. Resume with the SAME objective and Hessian
/// callbacks. This closure-based engine cannot identify a changed problem.
#[derive(Debug, Clone)]
pub struct TrustRegionState {
    x: Vec<f64>,
    f: f64,
    g: Vec<f64>,
    delta: f64,
    iters: usize,
    evals: usize,
    hv_evals: usize,
    negative_curvature_hits: usize,
    history: Vec<f64>,
    stalled: bool,
}

fn inf_norm(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0.0f64, f64::max)
}

/// Validate before any callback; every budget leaf is a hard ceiling even
/// under `All`, where ordinary boolean composition would permit overspending.
fn admit_rule(rule: &StopRule, evals: usize) -> Option<usize> {
    match rule {
        StopRule::GradNorm(t) => {
            assert!(t.is_finite() && *t >= 0.0, "invalid trust-region gradient tolerance");
            None
        }
        StopRule::ObjectiveBelow(t) => {
            assert!(t.is_finite(), "invalid trust-region objective target");
            None
        }
        StopRule::Budget(b) => {
            assert!(*b >= evals, "trust-region budget is below retained evaluations");
            Some(*b)
        }
        StopRule::Stall { rel, window } => {
            assert!(rel.is_finite() && *rel >= 0.0, "invalid trust-region stall tolerance");
            assert!(*window > 0 && *window < usize::MAX, "invalid trust-region stall window");
            None
        }
        StopRule::Any(rules) | StopRule::All(rules) => {
            assert!(!rules.is_empty(), "empty trust-region stopping rule");
            rules.iter().filter_map(|rule| admit_rule(rule, evals)).min()
        }
    }
}

impl TrustRegionState {
    /// Start from a finite point, spending one objective/gradient evaluation.
    ///
    /// # Panics
    /// Panics for a non-finite point/objective/gradient or a gradient of the
    /// wrong dimension, matching the other ASCENT engines' admission policy.
    #[must_use]
    pub fn new(x0: &[f64], fg: crate::FnGrad<'_>) -> Self {
        assert!(x0.iter().all(|v| v.is_finite()), "non-finite trust-region start");
        let (f, g) = fg(x0);
        assert_eq!(g.len(), x0.len(), "trust-region gradient dimension mismatch");
        assert!(
            f.is_finite() && g.iter().all(|v| v.is_finite()),
            "non-finite trust-region initial evaluation"
        );
        Self {
            x: x0.to_vec(),
            f,
            g,
            delta: 1.0,
            iters: 0,
            evals: 1,
            hv_evals: 0,
            negative_curvature_hits: 0,
            history: vec![f],
            stalled: false,
        }
    }

    /// Current trust radius, retained across checkpoints and rejected trials.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.delta
    }

    /// Gradient at the last accepted point.
    #[must_use]
    pub fn gradient(&self) -> &[f64] {
        &self.g
    }

    /// Initial objective followed by one entry per completed outer iteration,
    /// including rejected trials so stagnation rules observe lack of progress.
    #[must_use]
    pub fn history(&self) -> &[f64] {
        &self.history
    }

    /// Snapshot of the last accepted point; this does not evaluate the model.
    #[must_use]
    pub fn report(&self) -> TrustRegionReport {
        TrustRegionReport {
            x: self.x.clone(),
            f: self.f,
            grad_norm: inf_norm(&self.g),
            iters: self.iters,
            evals: self.evals,
            hv_evals: self.hv_evals,
            negative_curvature_hits: self.negative_curvature_hits,
        }
    }

    /// Run at most `max_iters` ADDITIONAL outer iterations. Budget leaves
    /// count cumulative objective/gradient callbacks (not Hessian products)
    /// and take precedence over simultaneous convergence. The constructor's
    /// initial callback is already included. Invalid rules panic before work.
    pub fn run(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
    ) -> TrustRegionRunReport {
        self.run_with_pause(fg, hv_at, rule, max_iters, &mut || false)
    }

    /// As [`Self::run`], polling cancellation only at complete, resume-safe
    /// outer-iteration boundaries. Callback execution and an in-flight
    /// Steihaug solve are not interruptible here. A pause leaves the accepted
    /// state untouched and can be resumed with a fresh cancellation context.
    pub fn run_cancellable(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        cx: &Cx,
    ) -> TrustRegionRunReport {
        self.run_with_pause(fg, hv_at, rule, max_iters, &mut || {
            cx.checkpoint().is_err()
        })
    }

    fn run_with_pause(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        paused: &mut dyn FnMut() -> bool,
    ) -> TrustRegionRunReport {
        let budget = admit_rule(rule, self.evals);
        let mut completed = 0usize;
        let progress = loop {
            if paused() {
                break TrustRegionProgress::Paused;
            }
            if budget.is_some_and(|cap| self.evals >= cap) {
                break TrustRegionProgress::Stopped(StopReason::Budget);
            }
            let obs = StopObservation {
                grad_norm: inf_norm(&self.g),
                objective: self.f,
                evals: self.evals,
                history: &self.history,
            };
            if let Some(reason) = rule.check(&obs) {
                break TrustRegionProgress::Stopped(reason);
            }
            if self.stalled || self.delta < 1e-14 {
                break TrustRegionProgress::Stopped(StopReason::Stall);
            }
            if completed == max_iters {
                break TrustRegionProgress::Stopped(StopReason::IterationCap);
            }
            self.stalled = !self.iterate(fg, hv_at);
            self.iters += 1;
            self.history.push(self.f);
            completed += 1;
        };
        TrustRegionRunReport {
            progress,
            solution: self.report(),
        }
    }

    // false means an unusable model/precision stall; no candidate is accepted.
    fn iterate(&mut self, fg: crate::FnGrad<'_>, hv_at: crate::FnHv<'_>) -> bool {
        let step = {
            let xc = self.x.clone();
            let mut hv = |v: &[f64]| hv_at(&xc, v);
            steihaug(&self.g, &mut hv, self.delta, 1e-8)
        };
        let (p, _hit, neg, hv_count) = match step {
            Ok(step) => step,
            Err(hv_count) => {
                self.hv_evals += hv_count;
                return false;
            }
        };
        self.hv_evals += hv_count;
        if neg {
            self.negative_curvature_hits += 1;
        }
        // Keep the established floating-point operation order on valid runs.
        let hp = hv_at(&self.x, &p);
        self.hv_evals += 1;
        assert_eq!(hp.len(), self.x.len(), "trust-region Hessian dimension mismatch");
        if hp.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let gp: f64 = self.g.iter().zip(&p).map(|(a, b)| a * b).sum();
        let php: f64 = p.iter().zip(&hp).map(|(a, b)| a * b).sum();
        let model_decrease = -gp - 0.5 * php;
        if !model_decrease.is_finite() || model_decrease <= 0.0 {
            return false;
        }
        let x_new: Vec<f64> = self.x.iter().zip(&p).map(|(a, b)| a + b).collect();
        if x_new == self.x {
            return false;
        }
        if x_new.iter().any(|v| !v.is_finite()) {
            self.delta *= 0.25;
            return true;
        }
        let (f_new, g_new) = fg(&x_new);
        self.evals += 1;
        assert_eq!(g_new.len(), self.x.len(), "trust-region gradient dimension mismatch");
        // A trial outside the objective's domain is a rejected step, not a
        // new incumbent. In particular, -infinity must not look like an
        // infinite improvement, nor may NaN freeze the radius through rho.
        if !f_new.is_finite() || g_new.iter().any(|v| !v.is_finite()) {
            self.delta *= 0.25;
            return true;
        }
        let actual = self.f - f_new;
        let rho = if model_decrease.abs() < 1e-300 {
            0.0
        } else {
            actual / model_decrease
        };
        if !actual.is_finite() || !rho.is_finite() {
            self.delta *= 0.25;
            return true;
        }
        let p_norm: f64 = p.iter().map(|v| v * v).sum::<f64>().sqrt();
        if rho < 0.25 {
            self.delta *= 0.25;
        } else if rho > 0.75 && (p_norm - self.delta).abs() < 1e-10 * self.delta {
            self.delta = (2.0 * self.delta).min(1e8);
        }
        if rho > 1e-4 {
            self.x = x_new;
            self.f = f_new;
            self.g = g_new;
        }
        true
    }
}

/// Steihaug-CG: approximately minimize m(p) = gᵀp + ½pᵀHp within
/// ‖p‖ ≤ Δ. Returns (step, hit_boundary, negative_curvature, hv_count),
/// or the spent Hessian count if the model/arithmetic is unusable.
fn steihaug(
    g: &[f64],
    hv: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    delta: f64,
    tol: f64,
) -> Result<(Vec<f64>, bool, bool, usize), usize> {
    let n = g.len();
    let mut p = vec![0.0f64; n];
    let mut r: Vec<f64> = g.iter().map(|v| -v).collect();
    let mut d = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let g_norm = rr.sqrt();
    let mut hv_count = 0usize;
    if !rr.is_finite() {
        return Err(hv_count);
    }
    for _ in 0..n.saturating_mul(2) {
        if rr.sqrt() < tol * g_norm.max(1e-30) {
            return Ok((p, false, false, hv_count));
        }
        if d.iter().any(|v| !v.is_finite()) {
            return Err(hv_count);
        }
        let hd = hv(&d);
        hv_count += 1;
        assert_eq!(hd.len(), n, "trust-region Hessian dimension mismatch");
        if hd.iter().any(|v| !v.is_finite()) {
            return Err(hv_count);
        }
        let dhd: f64 = d.iter().zip(&hd).map(|(a, b)| a * b).sum();
        if !dhd.is_finite() {
            return Err(hv_count);
        }
        if dhd <= 0.0 {
            // Negative curvature: follow d to the boundary.
            let tau = boundary_tau(&p, &d, delta);
            for i in 0..n {
                p[i] = tau.mul_add(d[i], p[i]);
            }
            return if p.iter().all(|v| v.is_finite()) {
                Ok((p, true, true, hv_count))
            } else {
                Err(hv_count)
            };
        }
        let alpha = rr / dhd;
        if !alpha.is_finite() {
            return Err(hv_count);
        }
        let mut p_next = p.clone();
        for i in 0..n {
            p_next[i] = alpha.mul_add(d[i], p_next[i]);
        }
        let pn_norm: f64 = p_next.iter().map(|v| v * v).sum::<f64>().sqrt();
        if pn_norm >= delta {
            let tau = boundary_tau(&p, &d, delta);
            for i in 0..n {
                p[i] = tau.mul_add(d[i], p[i]);
            }
            return if p.iter().all(|v| v.is_finite()) {
                Ok((p, true, false, hv_count))
            } else {
                Err(hv_count)
            };
        }
        p = p_next;
        for i in 0..n {
            r[i] = alpha.mul_add(-hd[i], r[i]);
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        if !rr_new.is_finite() {
            return Err(hv_count);
        }
        let beta = rr_new / rr;
        if !beta.is_finite() {
            return Err(hv_count);
        }
        rr = rr_new;
        for i in 0..n {
            d[i] = beta.mul_add(d[i], r[i]);
        }
    }
    if p.iter().all(|v| v.is_finite()) {
        Ok((p, false, false, hv_count))
    } else {
        Err(hv_count)
    }
}

/// Positive τ with ‖p + τ·d‖ = Δ.
fn boundary_tau(p: &[f64], d: &[f64], delta: f64) -> f64 {
    let pd: f64 = p.iter().zip(d).map(|(a, b)| a * b).sum();
    let dd: f64 = d.iter().map(|v| v * v).sum();
    let pp: f64 = p.iter().map(|v| v * v).sum();
    let disc = pd.mul_add(pd, dd * (delta * delta - pp));
    (-pd + fs_math::det::sqrt(disc.max(0.0))) / dd
}

/// Trust-region Newton–Krylov with the classical radius laws
/// (shrink ×¼ below ρ = ¼, grow ×2 above ρ = ¾ at the boundary).
/// `fg` returns (f, gradient); `hv` is the Hessian-vector product at
/// the CURRENT iterate (the driver re-binds it per iterate).
pub fn trust_region_newton(
    x0: &[f64],
    fg: crate::FnGrad<'_>,
    hv_at: crate::FnHv<'_>,
    grad_tol: f64,
    max_iters: usize,
) -> TrustRegionReport {
    let rule = StopRule::GradNorm(grad_tol);
    let _ = admit_rule(&rule, 0);
    let mut state = TrustRegionState::new(x0, fg);
    state.run(fg, hv_at, &rule, max_iters).solution
}

/// Finite-difference-of-gradients Hessian-vector product: the interim
/// path until second-order adjoints land. Accuracy is O(√ε)·‖H‖ — in
/// the NAME and the doc, not hidden: prefer exact Hv (duals, or
/// adjoint-of-adjoint when it ships) wherever it reaches.
pub fn hv_fd_of_gradients(fg: crate::FnGrad<'_>, x: &[f64], v: &[f64], eps: f64) -> Vec<f64> {
    let xp: Vec<f64> = x
        .iter()
        .zip(v)
        .map(|(xi, vi)| eps.mul_add(*vi, *xi))
        .collect();
    let xm: Vec<f64> = x
        .iter()
        .zip(v)
        .map(|(xi, vi)| eps.mul_add(-vi, *xi))
        .collect();
    let (_, gp) = fg(&xp);
    let (_, gm) = fg(&xm);
    gp.iter()
        .zip(&gm)
        .map(|(a, b)| (a - b) / (2.0 * eps))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
        let residual = x[1] - x[0] * x[0];
        (
            100.0 * residual * residual + (1.0 - x[0]).powi(2),
            vec![-400.0 * x[0] * residual - 2.0 * (1.0 - x[0]), 200.0 * residual],
        )
    }

    fn rosenbrock_hv(x: &[f64], v: &[f64]) -> Vec<f64> {
        vec![
            (1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0) * v[0]
                - 400.0 * x[0] * v[1],
            -400.0 * x[0] * v[0] + 200.0 * v[1],
        ]
    }

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    fn assert_same_state(a: &TrustRegionState, b: &TrustRegionState) {
        assert_eq!(bits(&a.x), bits(&b.x));
        assert_eq!(a.f.to_bits(), b.f.to_bits());
        assert_eq!(bits(&a.g), bits(&b.g));
        assert_eq!(a.delta.to_bits(), b.delta.to_bits());
        assert_eq!(bits(&a.history), bits(&b.history));
        assert_eq!(a.iters, b.iters);
        assert_eq!(a.evals, b.evals);
        assert_eq!(a.hv_evals, b.hv_evals);
        assert_eq!(a.negative_curvature_hits, b.negative_curvature_hits);
        assert_eq!(a.stalled, b.stalled);
    }

    #[test]
    fn trust_region_every_split_preserves_complete_state() {
        let initial = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
        let rule = StopRule::GradNorm(1e-12);
        let mut straight = initial.clone();
        straight.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12);
        for split in 0..=12 {
            let mut segmented = initial.clone();
            segmented.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, split);
            let mut resumed = segmented.clone();
            resumed.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12 - split);
            assert_same_state(&straight, &resumed);
        }
    }

    #[test]
    fn trust_region_pause_is_mutation_free_and_resumable() {
        let initial = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
        let rule = StopRule::GradNorm(1e-12);
        let mut paused = initial.clone();
        let mut boundaries = 0;
        let outcome = paused.run_with_pause(
            &mut rosenbrock, &mut rosenbrock_hv, &rule, 10,
            &mut || { boundaries += 1; boundaries > 2 },
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Paused);
        assert_eq!(outcome.solution.iters, 2);
        let checkpoint = paused.clone();
        let repeated = paused.run_with_pause(
            &mut |_| panic!("objective evaluated during pause"),
            &mut |_, _| panic!("Hessian evaluated during pause"),
            &rule, 10, &mut || true,
        );
        assert_eq!(repeated.progress, TrustRegionProgress::Paused);
        assert_same_state(&checkpoint, &paused);
        paused.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 10);
        let mut straight = initial;
        straight.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12);
        assert_same_state(&straight, &paused);
    }

    #[test]
    fn trust_region_nested_all_budget_is_a_hard_callback_ceiling() {
        use std::cell::Cell;
        let calls = Cell::new(0);
        let mut fg = |x: &[f64]| { calls.set(calls.get() + 1); rosenbrock(x) };
        let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut fg);
        let rule = StopRule::All(vec![
            StopRule::GradNorm(0.0),
            StopRule::Any(vec![StopRule::Budget(3), StopRule::Budget(5)]),
        ]);
        let outcome = state.run(&mut fg, &mut rosenbrock_hv, &rule, 100);
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Budget));
        assert_eq!(outcome.solution.evals, 3);
        assert_eq!(calls.get(), 3);
        let checkpoint = state.clone();
        state.run(&mut fg, &mut rosenbrock_hv, &rule, 100);
        assert_eq!(calls.get(), 3);
        assert_same_state(&checkpoint, &state);
    }

    #[test]
    fn trust_region_budget_has_priority_over_convergence() {
        let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
        let mut state = TrustRegionState::new(&[0.0], &mut fg);
        let outcome = state.run(
            &mut fg, &mut |_, _| panic!("unexpected Hessian evaluation"),
            &StopRule::Any(vec![StopRule::GradNorm(1e-8), StopRule::Budget(1)]), 100,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Budget));
        assert_eq!(outcome.solution.iters, 0);
        assert_eq!(outcome.solution.evals, 1);
    }

    #[test]
    fn trust_region_zero_iteration_segment_spends_no_callbacks() {
        let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
        let checkpoint = state.clone();
        let outcome = state.run(
            &mut |_| panic!("unexpected objective evaluation"),
            &mut |_, _| panic!("unexpected Hessian evaluation"),
            &StopRule::GradNorm(1e-12), 0,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::IterationCap));
        assert_same_state(&checkpoint, &state);
    }

    #[test]
    fn trust_region_wrapper_and_state_return_the_same_solution() {
        let legacy = trust_region_newton(
            &[-1.2, 1.0], &mut rosenbrock, &mut rosenbrock_hv, 1e-8, 100,
        );
        let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
        let run = state.run(&mut rosenbrock, &mut rosenbrock_hv, &StopRule::GradNorm(1e-8), 100);
        assert_eq!(run.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
        assert!(run.solution.grad_norm <= 1e-8);
        assert_eq!(bits(&legacy.x), bits(&run.solution.x));
        assert_eq!(legacy.f.to_bits(), run.solution.f.to_bits());
        assert_eq!(legacy.grad_norm.to_bits(), run.solution.grad_norm.to_bits());
        assert_eq!(legacy.iters, run.solution.iters);
        assert_eq!(legacy.evals, run.solution.evals);
        assert_eq!(legacy.hv_evals, run.solution.hv_evals);
        assert_eq!(legacy.negative_curvature_hits, run.solution.negative_curvature_hits);
    }

    #[test]
    #[should_panic(expected = "budget is below retained evaluations")]
    fn trust_region_refuses_underwritten_resume_budget() {
        let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
        state.run(
            &mut |_| panic!("unexpected objective evaluation"),
            &mut |_, _| panic!("unexpected Hessian evaluation"),
            &StopRule::Budget(0), 1,
        );
    }

    #[test]
    fn trust_region_rejects_nonfinite_objectives_without_poisoning_state() {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut fg = |x: &[f64]| {
                let f = if x[0] == 1.0 { invalid } else { (x[0] - 2.0).powi(2) };
                (f, vec![2.0 * (x[0] - 2.0)])
            };
            let mut hv = |_: &[f64], v: &[f64]| vec![2.0 * v[0]];
            let mut state = TrustRegionState::new(&[0.0], &mut fg);
            let run = state.run(&mut fg, &mut hv, &StopRule::GradNorm(1e-10), 1);
            assert_eq!(run.solution.x, vec![0.0]);
            assert_eq!(run.solution.f, 4.0);
            assert_eq!(state.gradient(), &[-4.0]);
            assert_eq!(state.radius(), 0.25);
            assert_eq!(state.history(), &[4.0, 4.0]);
            assert_eq!(run.solution.evals, 2);
            let recovered = state.run(&mut fg, &mut hv, &StopRule::GradNorm(1e-10), 100);
            assert_eq!(recovered.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
            assert!((recovered.solution.x[0] - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn trust_region_finite_objective_with_nan_gradient_is_not_accepted() {
        let mut fg = |x: &[f64]| {
            let g = if x[0] == 1.0 { f64::NAN } else { 2.0 * (x[0] - 2.0) };
            ((x[0] - 2.0).powi(2), vec![g])
        };
        let mut state = TrustRegionState::new(&[0.0], &mut fg);
        let run = state.run(
            &mut fg, &mut |_, v| vec![2.0 * v[0]], &StopRule::GradNorm(1e-10), 1,
        );
        assert_eq!(run.solution.x, vec![0.0]);
        assert_eq!(run.solution.grad_norm, 4.0);
        assert_eq!(state.radius(), 0.25);
    }

    #[test]
    fn trust_region_nonfinite_hessian_stalls_without_objective_trials() {
        let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
        let mut state = TrustRegionState::new(&[1.0], &mut fg);
        let outcome = state.run(
            &mut |_| panic!("invalid Hessian must not reach the objective"),
            &mut |_, _| vec![f64::NAN], &StopRule::GradNorm(1e-8), 100,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
        assert_eq!(outcome.solution.x, vec![1.0]);
        assert_eq!(outcome.solution.grad_norm, 2.0);
        assert_eq!(outcome.solution.evals, 1);
        assert_eq!(outcome.solution.hv_evals, 1);
        let checkpoint = state.clone();
        state.run(
            &mut |_| panic!("stalled state evaluated the objective"),
            &mut |_, _| panic!("stalled state evaluated the Hessian"),
            &StopRule::GradNorm(1e-8), 100,
        );
        assert_same_state(&checkpoint, &state);
    }

    #[test]
    fn trust_region_nonpositive_model_cannot_accept_an_uphill_step() {
        let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
        let mut state = TrustRegionState::new(&[1.0], &mut fg);
        let mut calls = 0;
        let mut inconsistent_hv = |_: &[f64], v: &[f64]| {
            calls += 1;
            let scale = if calls == 1 { 2.0 } else { 10.0 };
            vec![scale * v[0]]
        };
        let outcome = state.run(
            &mut |_| panic!("nonpositive model must not reach the objective"),
            &mut inconsistent_hv, &StopRule::GradNorm(1e-8), 100,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
        assert_eq!(outcome.solution.x, vec![1.0]);
        assert_eq!(outcome.solution.evals, 1);
        assert_eq!(outcome.solution.hv_evals, 2);
    }

    #[test]
    fn trust_region_precision_noop_stops_instead_of_reevaluating() {
        let mut state = TrustRegionState::new(&[1e100], &mut |_| (0.0, vec![-1.0]));
        let outcome = state.run(
            &mut |_| panic!("unchanged point must not be reevaluated"),
            &mut |_, _| vec![0.0], &StopRule::GradNorm(1e-8), 100,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
        assert_eq!(outcome.solution.evals, 1);
        assert_eq!(outcome.solution.x, vec![1e100]);
    }

    #[test]
    fn trust_region_internal_overflow_is_not_convergence() {
        let mut state = TrustRegionState::new(&[0.0], &mut |_| (1.0, vec![1e308]));
        let outcome = state.run(
            &mut |_| panic!("overflowed model reached objective"),
            &mut |_, _| panic!("overflowed residual reached Hessian"),
            &StopRule::GradNorm(1e-8), 10,
        );
        assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
        assert_eq!(outcome.solution.evals, 1);
        assert_eq!(outcome.solution.hv_evals, 0);
        assert_eq!(outcome.solution.grad_norm, 1e308);
    }

    #[test]
    #[should_panic(expected = "Hessian dimension mismatch")]
    fn trust_region_refuses_malformed_hessian_products() {
        let mut state = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![2.0]));
        state.run(
            &mut |_| panic!("bad Hessian reached objective"),
            &mut |_, _| vec![], &StopRule::GradNorm(1e-8), 1,
        );
    }

    #[test]
    #[should_panic(expected = "gradient dimension mismatch")]
    fn trust_region_refuses_malformed_trial_gradients() {
        let mut state = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![2.0]));
        state.run(
            &mut |_| (0.0, vec![]), &mut |_, v| vec![2.0 * v[0]],
            &StopRule::GradNorm(1e-8), 1,
        );
    }

    #[test]
    #[should_panic(expected = "initial evaluation")]
    fn trust_region_refuses_nan_initial_gradients() {
        let _ = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![f64::NAN]));
    }
}

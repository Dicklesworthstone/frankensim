//! Trust-region Newton–Krylov: Steihaug-CG on the quadratic model
//! with NEGATIVE-CURVATURE handling (follow the direction to the
//! boundary — the feature that separates TR from line-search Newton
//! on nonconvex terrain), classical radius update laws (G0-tested),
//! and matrix-free Hessian-vector products via caller-supplied
//! closures. Second-order adjoints are recorded follow-up; the
//! finite-difference-of-gradients Hv helper carries its O(√ε)
//! accuracy in its name rather than hiding it.

use crate::stop::{StopObservation, StopReason, StopRule};

#[path = "engine/range.rs"]
mod range;

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
    /// Cancellation was observed before publishing the next complete iteration.
    Paused,
    /// The next Krylov or final-model Hessian product exceeds its cumulative
    /// ceiling. The incomplete trial is discarded, but spent work is retained.
    HessianBudget,
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

enum IterationOutcome {
    Complete { usable: bool, negative_curvature: bool },
    Paused,
    HessianBudget,
}

enum SteihaugStop {
    Unusable(usize),
    Paused(usize),
    Budget(usize),
}

fn inf_norm(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0.0f64, f64::max)
}

// Multiplying an objective by a positive constant must not turn its Newton
// step into an overflow or a zero residual. Outside the ordinary exponent
// range, solve the equivalent model (g/s, H/s) with an exact power-of-two s.
// Keeping s = 1 in the ordinary range preserves the established trajectories.
// The range leaves headroom for the squared residual and cubic curvature
// products in CG; it is not an admission guarantee for arbitrary conditioning.
fn quadratic_scale(g: &[f64]) -> f64 {
    let magnitude = inf_norm(g);
    let bits = magnitude.to_bits();
    let exponent = bits >> 52;
    if magnitude == 0.0 || (1023 - 64..=1023 + 64).contains(&exponent) {
        1.0
    } else if exponent != 0 {
        f64::from_bits(exponent << 52)
    } else {
        // Subnormal gradients have no implicit leading bit. Their highest
        // stored mantissa bit still gives a nonzero, exact binary scale.
        f64::from_bits(1u64 << (63 - bits.leading_zeros()))
    }
}

fn scaled_difference(a: f64, b: f64, scale: f64) -> f64 {
    let difference = a - b;
    if difference.is_finite() {
        difference / scale
    } else {
        // Opposite finite objective values can overflow on subtraction even
        // when the scaled actual decrease and agreement ratio are finite.
        a / scale - b / scale
    }
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

    /// As [`Self::run`], with an additional cumulative ceiling on ALL Hessian
    /// products, including the final quadratic-model evaluation. No placeholder
    /// product is sent to the solver on exhaustion. A last permitted product
    /// may complete an iteration; a limit hit midway preserves the accepted
    /// point, radius and history, without refunding products already attempted.
    /// Raising the ceiling resumes from that point and restarts its Krylov solve.
    /// Objective/gradient limits remain the hard `Budget` leaves of `rule`.
    pub fn run_bounded(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        max_hessian_products: usize,
    ) -> TrustRegionRunReport {
        self.run_with_limits(fg, hv_at, rule, max_iters, max_hessian_products, &mut || false)
    }

    pub(super) fn run_with_pause(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        paused: &mut dyn FnMut() -> bool,
    ) -> TrustRegionRunReport {
        self.run_with_limits(fg, hv_at, rule, max_iters, usize::MAX, paused)
    }

    pub(super) fn run_with_limits(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        rule: &StopRule,
        max_iters: usize,
        max_hessian_products: usize,
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
            match self.iterate(fg, hv_at, max_hessian_products, paused) {
                IterationOutcome::Paused => break TrustRegionProgress::Paused,
                IterationOutcome::HessianBudget => break TrustRegionProgress::HessianBudget,
                IterationOutcome::Complete { usable, negative_curvature } => {
                    self.stalled = !usable;
                    self.negative_curvature_hits += usize::from(negative_curvature);
                }
            }
            self.iters += 1;
            self.history.push(self.f);
            completed += 1;
        };
        TrustRegionRunReport {
            progress,
            solution: self.report(),
        }
    }

    // Only callback counters change before the final cancellation checkpoint.
    // An unusable model/precision stall never publishes a candidate.
    fn iterate(
        &mut self,
        fg: crate::FnGrad<'_>,
        hv_at: crate::FnHv<'_>,
        max_hessian_products: usize,
        paused: &mut dyn FnMut() -> bool,
    ) -> IterationOutcome {
        let scale = quadratic_scale(&self.g);
        let scaled_gradient =
            (scale != 1.0).then(|| self.g.iter().map(|v| v / scale).collect::<Vec<_>>());
        let gradient = scaled_gradient.as_deref().unwrap_or(&self.g);
        let step = {
            let xc = self.x.clone();
            let mut hv = |v: &[f64]| {
                let mut product = hv_at(&xc, v);
                if scale != 1.0 {
                    for value in &mut product {
                        *value /= scale;
                    }
                }
                product
            };
            steihaug(gradient, &mut hv, self.delta, 1e-8,
                max_hessian_products.saturating_sub(self.hv_evals), paused)
        };
        let (p, _hit, neg, hv_count) = match step {
            Ok(step) => step,
            Err(SteihaugStop::Unusable(hv_count)) => {
                self.hv_evals += hv_count;
                return IterationOutcome::Complete { usable: false, negative_curvature: false };
            }
            Err(SteihaugStop::Paused(hv_count)) => {
                self.hv_evals += hv_count;
                return IterationOutcome::Paused;
            }
            Err(SteihaugStop::Budget(hv_count)) => {
                self.hv_evals += hv_count;
                return IterationOutcome::HessianBudget;
            }
        };
        self.hv_evals += hv_count;
        let complete = |usable| IterationOutcome::Complete { usable, negative_curvature: neg };
        if paused() {
            return IterationOutcome::Paused;
        }
        if self.hv_evals >= max_hessian_products {
            return IterationOutcome::HessianBudget;
        }
        // Keep the established floating-point operation order on valid runs.
        let mut hp = hv_at(&self.x, &p);
        self.hv_evals += 1;
        if paused() {
            return IterationOutcome::Paused;
        }
        assert_eq!(hp.len(), self.x.len(), "trust-region Hessian dimension mismatch");
        if hp.iter().any(|v| !v.is_finite()) {
            return complete(false);
        }
        if scale != 1.0 {
            for value in &mut hp {
                *value /= scale;
            }
        }
        let gp: f64 = gradient.iter().zip(&p).map(|(a, b)| a * b).sum();
        let php: f64 = p.iter().zip(&hp).map(|(a, b)| a * b).sum();
        let model_decrease = -gp - 0.5 * php;
        if !model_decrease.is_finite() || model_decrease <= 0.0 {
            return complete(false);
        }
        let x_new: Vec<f64> = self.x.iter().zip(&p).map(|(a, b)| a + b).collect();
        if x_new == self.x {
            return complete(false);
        }
        if x_new.iter().any(|v| !v.is_finite()) {
            self.delta *= 0.25;
            return complete(true);
        }
        if paused() {
            return IterationOutcome::Paused;
        }
        let (f_new, g_new) = fg(&x_new);
        self.evals += 1;
        if paused() {
            return IterationOutcome::Paused;
        }
        assert_eq!(g_new.len(), self.x.len(), "trust-region gradient dimension mismatch");
        // A trial outside the objective's domain is a rejected step, not a
        // new incumbent. In particular, -infinity must not look like an
        // infinite improvement, nor may NaN freeze the radius through rho.
        if !f_new.is_finite() || g_new.iter().any(|v| !v.is_finite()) {
            self.delta *= 0.25;
            return complete(true);
        }
        let actual = scaled_difference(self.f, f_new, scale);
        // Positivity and finiteness were checked above. An absolute floor
        // would make identical models disagree after changing objective units.
        let rho = actual / model_decrease;
        if !actual.is_finite() || !rho.is_finite() {
            self.delta *= 0.25;
            return complete(true);
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
        complete(true)
    }
}

/// Steihaug-CG: approximately minimize m(p) = gᵀp + ½pᵀHp within
/// ‖p‖ ≤ Δ. Returns (step, hit_boundary, negative_curvature, hv_count),
/// or the stop reason with the spent Hessian count.
fn steihaug(
    g: &[f64],
    hv: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    delta: f64,
    tol: f64,
    max_hv: usize,
    paused: &mut dyn FnMut() -> bool,
) -> Result<(Vec<f64>, bool, bool, usize), SteihaugStop> {
    let n = g.len();
    let mut p = vec![0.0f64; n];
    let mut r: Vec<f64> = g.iter().map(|v| -v).collect();
    let mut d = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let g_norm = rr.sqrt();
    let mut hv_count = 0usize;
    if !rr.is_finite() {
        return Err(SteihaugStop::Unusable(hv_count));
    }
    for _ in 0..n.saturating_mul(2) {
        if rr.sqrt() < tol * g_norm.max(1e-30) {
            return Ok((p, false, false, hv_count));
        }
        if d.iter().any(|v| !v.is_finite()) {
            return Err(SteihaugStop::Unusable(hv_count));
        }
        if paused() {
            return Err(SteihaugStop::Paused(hv_count));
        }
        if hv_count >= max_hv {
            return Err(SteihaugStop::Budget(hv_count));
        }
        let hd = hv(&d);
        hv_count += 1;
        if paused() {
            return Err(SteihaugStop::Paused(hv_count));
        }
        assert_eq!(hd.len(), n, "trust-region Hessian dimension mismatch");
        if hd.iter().any(|v| !v.is_finite()) {
            return Err(SteihaugStop::Unusable(hv_count));
        }
        let (negative_curvature, alpha) = range::curvature(rr, &d, &hd);
        if negative_curvature {
            // Negative curvature: follow d to the boundary.
            return if range::to_boundary(&mut p, &d, delta) {
                Ok((p, true, true, hv_count))
            } else {
                Err(SteihaugStop::Unusable(hv_count))
            };
        }
        let mut p_next = p.clone();
        for i in 0..n {
            p_next[i] = alpha.mul_add(d[i], p_next[i]);
        }
        let pn_norm: f64 = p_next.iter().map(|v| v * v).sum::<f64>().sqrt();
        // Positive curvature can be so small that the unrestricted CG step
        // overflows. Its finite trust-boundary truncation is still usable.
        if !alpha.is_finite() || pn_norm >= delta {
            return if range::to_boundary(&mut p, &d, delta) {
                Ok((p, true, false, hv_count))
            } else {
                Err(SteihaugStop::Unusable(hv_count))
            };
        }
        p = p_next;
        for i in 0..n {
            r[i] = alpha.mul_add(-hd[i], r[i]);
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        if !rr_new.is_finite() {
            return Err(SteihaugStop::Unusable(hv_count));
        }
        let beta = rr_new / rr;
        if !beta.is_finite() {
            return Err(SteihaugStop::Unusable(hv_count));
        }
        rr = rr_new;
        for i in 0..n {
            d[i] = beta.mul_add(d[i], r[i]);
        }
    }
    if p.iter().all(|v| v.is_finite()) {
        Ok((p, false, false, hv_count))
    } else {
        Err(SteihaugStop::Unusable(hv_count))
    }
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
#[path = "engine/tests.rs"]
mod tests;
#[cfg(test)]
#[path = "engine/budget_tests.rs"]
mod budget_tests;

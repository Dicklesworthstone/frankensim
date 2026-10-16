//! Riemannian L-BFGS over fs-opt's authoritative manifold runtime:
//! ambient point gradient → parameter gradient → strong-Wolfe search along
//! retracted curves, with memory pairs moved by the transport associated with
//! each manifold. Point storage and optimizer-parameter storage remain
//! distinct end to end (notably SO(3): four quaternion lanes versus three
//! right/body lanes). Manifold invariants hold to roundoff along the whole
//! iterate path — tested, not assumed.

use crate::stop::{StopObservation, StopReason, StopRule};
use crate::wolfe::strong_wolfe_with_budget;
use fs_opt::Manifold;
use std::collections::VecDeque;

type CurvaturePair = (Vec<f64>, Vec<f64>, f64);

const CURVATURE_RELATIVE_FLOOR: f64 = 1e-14;

fn point_dim(man: &Manifold) -> usize {
    usize::try_from(
        man.layout()
            .expect("Riemannian solver requires a valid manifold descriptor")
            .point_dim()
            .get(),
    )
    .expect("manifold point dimension must fit usize")
}

fn parameter_dim(man: &Manifold) -> usize {
    usize::try_from(
        man.layout()
            .expect("Riemannian solver requires a valid manifold descriptor")
            .param_dim()
            .get(),
    )
    .expect("manifold parameter dimension must fit usize")
}

fn assert_finite(label: &str, values: &[f64]) {
    assert!(
        values.iter().all(|value| value.is_finite()),
        "{label} entries must be finite"
    );
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "dot-product dimensions must match");
    a.iter().zip(b).map(|(left, right)| left * right).sum()
}

fn l2_norm(values: &[f64]) -> f64 {
    assert!(!values.is_empty(), "norm input must not be empty");
    assert_finite("norm input", values);
    // Deterministic scaled sum-of-squares (the classic LASSQ recurrence)
    // avoids both fabricating zero for subnormal squares and overflowing the
    // intermediate sum for representable inputs.
    let mut scale = 0.0f64;
    let mut sumsq = 1.0f64;
    for value in values {
        let magnitude = value.abs();
        if magnitude == 0.0 {
            continue;
        }
        if scale < magnitude {
            let ratio = scale / magnitude;
            sumsq = 1.0 + sumsq * ratio * ratio;
            scale = magnitude;
        } else {
            let ratio = magnitude / scale;
            sumsq += ratio * ratio;
        }
    }
    let norm = if scale == 0.0 {
        0.0
    } else {
        scale * fs_math::det::sqrt(sumsq)
    };
    assert!(norm.is_finite(), "norm must remain finite");
    norm
}

fn inf_norm(values: &[f64]) -> f64 {
    assert!(!values.is_empty(), "norm input must not be empty");
    assert_finite("norm input", values);
    values
        .iter()
        .map(|value| value.abs())
        .fold(0.0f64, f64::max)
}

fn violation(man: &Manifold, x: &[f64]) -> f64 {
    match man {
        Manifold::Sphere { .. } | Manifold::So3 => (l2_norm(x) - 1.0).abs(),
        Manifold::Stiefel { n, p } => {
            let (n, p) = (*n as usize, *p as usize);
            let mut residual = 0.0_f64;
            for column in 0..p {
                for against in 0..=column {
                    let gram: f64 = (0..n)
                        .map(|row| x[column * n + row] * x[against * n + row])
                        .sum();
                    let expected = if column == against { 1.0 } else { 0.0 };
                    residual = residual.max((gram - expected).abs());
                }
            }
            residual
        }
        Manifold::Rn { .. } => 0.0,
    }
}

fn authoritative_point(man: &Manifold, x: &[f64]) -> Vec<f64> {
    assert_eq!(
        x.len(),
        point_dim(man),
        "manifold point dimension does not match its descriptor"
    );
    assert_finite("manifold point", x);
    let zero_step = vec![0.0; parameter_dim(man)];
    let validated = man
        .retract(x, &zero_step)
        .unwrap_or_else(|error| panic!("manifold point failed fs-opt authority: {error}"));
    // SO(3) has an explicit antipodal representative contract. The other
    // manifolds have no quotient-representative ambiguity, so retain the
    // caller's already-authority-validated point bits instead of silently
    // moving a Sphere/Stiefel start by another normalization/QR roundoff.
    if matches!(man, Manifold::So3) {
        validated
    } else {
        x.to_vec()
    }
}

fn assert_point(man: &Manifold, x: &[f64]) {
    let _canonical = authoritative_point(man, x);
}

fn assert_tangent(man: &Manifold, x: &[f64], vector: &[f64]) {
    assert_eq!(
        vector.len(),
        parameter_dim(man),
        "tangent-vector dimension must match the manifold parameter layout"
    );
    assert_finite("tangent vector", vector);
    man.validate_parameter_tangent(x, vector)
        .unwrap_or_else(|error| panic!("tangent vector failed fs-opt authority: {error}"));
}

fn checked_fg(man: &Manifold, fg: crate::FnGrad<'_>, x: &[f64]) -> (f64, Vec<f64>) {
    let (f, ambient_gradient) = fg(x);
    assert!(f.is_finite(), "objective value must be finite");
    assert_eq!(
        ambient_gradient.len(),
        x.len(),
        "objective callback gradient dimension must match point storage"
    );
    assert_finite("objective ambient gradient", &ambient_gradient);
    let parameter_gradient = man
        .parameter_gradient(x, &ambient_gradient)
        .unwrap_or_else(|error| panic!("objective gradient pullback failed: {error}"));
    (f, parameter_gradient)
}

/// Pull an ambient point-storage vector back to optimizer-parameter storage at
/// `x` through the authoritative fs-opt manifold operation.
///
/// For SO(3), the returned vector has three right/body lanes even though `x`
/// and `v` have four quaternion lanes. This public wrapper is retained for
/// compatibility; new code may call [`Manifold::parameter_gradient`] directly.
///
/// # Panics
/// If the descriptor, point, ambient vector, or resulting tangent is invalid.
#[must_use]
pub fn tangent_project(man: &Manifold, x: &[f64], v: &[f64]) -> Vec<f64> {
    man.parameter_gradient(x, v)
        .unwrap_or_else(|error| panic!("ambient vector failed fs-opt manifold pullback: {error}"))
}

/// Retract from point-storage `x` along optimizer-parameter `step` through the
/// authoritative fs-opt manifold operation.
///
/// # Panics
/// If the descriptor, point, step, or landing is invalid.
#[must_use]
pub fn retract(man: &Manifold, x: &[f64], step: &[f64]) -> Vec<f64> {
    man.retract(x, step)
        .unwrap_or_else(|error| panic!("retraction failed fs-opt manifold authority: {error}"))
}

fn project_parameter_direction(man: &Manifold, point: &[f64], direction: &[f64]) -> Vec<f64> {
    match man {
        Manifold::Rn { .. } | Manifold::So3 => {
            man.validate_parameter_tangent(point, direction)
                .unwrap_or_else(|error| panic!("parameter direction is invalid: {error}"));
            direction.to_vec()
        }
        Manifold::Sphere { .. } | Manifold::Stiefel { .. } => man
            .parameter_gradient(point, direction)
            .unwrap_or_else(|error| panic!("parameter direction projection failed: {error}")),
    }
}

fn vector_transport(
    man: &Manifold,
    from: &[f64],
    step: &[f64],
    to: &[f64],
    vector: &[f64],
) -> Vec<f64> {
    man.transport_parameter(from, step, to, vector)
        .unwrap_or_else(|error| panic!("vector transport failed fs-opt authority: {error}"))
}

fn curvature_rho(s: &[f64], y: &[f64]) -> Option<f64> {
    assert_eq!(s.len(), y.len(), "curvature-pair dimensions must match");
    let s_norm = l2_norm(s);
    let y_norm = l2_norm(y);
    if s_norm == 0.0 || y_norm == 0.0 {
        return None;
    }
    let sy = dot(s, y);
    // Compare the dimensionless curvature cosine rather than adding an
    // absolute floor to either vector norm. This keeps admission symmetric
    // and invariant under reciprocal unit rescaling of s and y.
    let relative_curvature: f64 = s
        .iter()
        .zip(y)
        .map(|(step, gradient_delta)| (step / s_norm) * (gradient_delta / y_norm))
        .sum();
    if !sy.is_finite()
        || sy <= 0.0
        || !relative_curvature.is_finite()
        || relative_curvature <= CURVATURE_RELATIVE_FLOOR
    {
        return None;
    }
    let rho = 1.0 / sy;
    rho.is_finite().then_some(rho)
}

fn transport_pairs(
    man: &Manifold,
    from: &[f64],
    step: &[f64],
    to: &[f64],
    pairs: &mut VecDeque<CurvaturePair>,
) {
    let mut transported = VecDeque::with_capacity(pairs.len());
    while let Some((s, y, _old_rho)) = pairs.pop_front() {
        let s = vector_transport(man, from, step, to, &s);
        let y = vector_transport(man, from, step, to, &y);
        if let Some(rho) = curvature_rho(&s, &y) {
            transported.push_back((s, y, rho));
        }
    }
    *pairs = transported;
}

fn evaluation_budget(rule: &StopRule) -> Option<usize> {
    match rule {
        StopRule::Budget(budget) => Some(*budget),
        StopRule::Any(children) | StopRule::All(children) => {
            children.iter().filter_map(evaluation_budget).min()
        }
        StopRule::GradNorm(_) | StopRule::ObjectiveBelow(_) | StopRule::Stall { .. } => None,
    }
}

fn assert_valid_stop_rule(rule: &StopRule) {
    match rule {
        StopRule::GradNorm(tolerance) => assert!(
            tolerance.is_finite() && *tolerance >= 0.0,
            "gradient-norm tolerance must be finite and nonnegative"
        ),
        StopRule::ObjectiveBelow(target) => {
            assert!(target.is_finite(), "objective target must be finite");
        }
        StopRule::Budget(_) => {}
        StopRule::Stall { rel, window } => {
            assert!(
                rel.is_finite() && *rel >= 0.0,
                "stagnation tolerance must be finite and nonnegative"
            );
            assert!(*window >= 1, "stagnation window must be at least one");
        }
        StopRule::Any(children) | StopRule::All(children) => {
            assert!(
                !children.is_empty(),
                "composite stop rules must contain at least one child"
            );
            for child in children {
                assert_valid_stop_rule(child);
            }
        }
    }
}

fn stop_reason(
    rule: &StopRule,
    hard_budget: Option<usize>,
    observation: &StopObservation<'_>,
) -> Option<StopReason> {
    if hard_budget.is_some_and(|budget| observation.evals >= budget) {
        Some(StopReason::Budget)
    } else {
        rule.check(observation)
    }
}

/// Outcome of a Riemannian L-BFGS run.
#[derive(Debug, Clone)]
pub struct RiemannianReport {
    /// Why the run stopped.
    pub reason: StopReason,
    /// Final Riemannian gradient norm (∞).
    pub grad_norm: f64,
    /// Final intrinsic Riemannian gradient norm (L2 under the induced metric).
    pub grad_l2_norm: f64,
    /// Final objective.
    pub f: f64,
    /// Iterations.
    pub iters: usize,
    /// Evaluations.
    pub evals: usize,
    /// Worst manifold-constraint violation observed along the path (unit-norm
    /// residual for Sphere/SO(3), Gram residual for Stiefel).
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
    /// Current Riemannian gradient in optimizer-parameter storage (SO(3):
    /// three body lanes while [`RiemannianLbfgs::x`] has four quaternion lanes).
    pub g: Vec<f64>,
    pairs: VecDeque<CurvaturePair>,
    memory: usize,
    /// Iterations performed.
    pub iters: usize,
    /// Evaluations spent.
    pub evals: usize,
    /// Objective history.
    pub history: Vec<f64>,
    worst_violation: f64,
}

impl RiemannianLbfgs {
    /// Start at a point ON the manifold.
    #[must_use]
    pub fn new(manifold: Manifold, x0: &[f64], memory: usize, fg: crate::FnGrad<'_>) -> Self {
        let x = authoritative_point(&manifold, x0);
        let (f, g) = checked_fg(&manifold, fg, &x);
        let worst = violation(&manifold, &x);
        RiemannianLbfgs {
            manifold,
            x,
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

    fn assert_valid_state(&self) {
        assert_point(&self.manifold, &self.x);
        assert!(self.f.is_finite(), "current objective must be finite");
        assert_eq!(
            self.g.len(),
            parameter_dim(&self.manifold),
            "current gradient dimension must match the manifold parameter layout"
        );
        assert_tangent(&self.manifold, &self.x, &self.g);
        assert!(
            self.evals >= 1,
            "Riemannian state must retain its initial evaluation"
        );
        let expected_history_len = self
            .iters
            .checked_add(1)
            .expect("iteration count must admit a history tail");
        assert_eq!(
            self.history.len(),
            expected_history_len,
            "objective history must contain the start plus every accepted iterate"
        );
        assert_finite("objective history", &self.history);
        assert_eq!(
            self.history.last().map(|value| value.to_bits()),
            Some(self.f.to_bits()),
            "objective history tail must equal the current objective"
        );
        assert!(
            self.worst_violation.is_finite()
                && self.worst_violation >= violation(&self.manifold, &self.x),
            "stored manifold-violation certificate must cover the current point"
        );
        assert!(
            self.pairs.len() <= self.memory,
            "curvature memory exceeds its configured capacity"
        );
        for (s, y, rho) in &self.pairs {
            assert_tangent(&self.manifold, &self.x, s);
            assert_tangent(&self.manifold, &self.x, y);
            let expected = curvature_rho(s, y)
                .expect("retained curvature pair must remain finite and positive");
            assert_eq!(
                rho.to_bits(),
                expected.to_bits(),
                "retained curvature reciprocal must match transported vectors"
            );
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

    fn report(&self, reason: StopReason) -> RiemannianReport {
        RiemannianReport {
            reason,
            grad_norm: inf_norm(&self.g),
            grad_l2_norm: l2_norm(&self.g),
            f: self.f,
            iters: self.iters,
            evals: self.evals,
            worst_violation: self.worst_violation,
        }
    }

    /// Run a budget-exact strong-Wolfe search along retracted curves.
    ///
    /// Every budget leaf in the stop algebra is a hard callback cap, even
    /// below `All`: no line-search phase may overshoot it. Directional
    /// derivatives pair the fs-opt parameter gradient with the authoritative
    /// retraction-curve parameter velocity.
    pub fn run(
        &mut self,
        fg: crate::FnGrad<'_>,
        rule: &StopRule,
        max_iters: usize,
    ) -> RiemannianReport {
        assert_valid_stop_rule(rule);
        self.assert_valid_state();
        let hard_budget = evaluation_budget(rule);
        if let Some(budget) = hard_budget {
            assert!(
                budget >= self.evals,
                "evaluation budget {budget} cannot underwrite {} callbacks already retained by the Riemannian state",
                self.evals,
            );
        }
        let mut reason = StopReason::IterationCap;
        let initial = StopObservation {
            grad_norm: inf_norm(&self.g),
            objective: self.f,
            evals: self.evals,
            history: &self.history,
        };
        if let Some(stopped) = stop_reason(rule, hard_budget, &initial) {
            return self.report(stopped);
        }

        'iterations: for _ in 0..max_iters {
            let raw_direction = self.direction();
            let raw_is_finite = raw_direction.iter().all(|value| value.is_finite());
            let mut direction = if raw_is_finite {
                project_parameter_direction(&self.manifold, &self.x, &raw_direction)
            } else {
                Vec::new()
            };
            let mut dphi0 = if raw_is_finite {
                dot(&direction, &self.g)
            } else {
                f64::NAN
            };
            if !dphi0.is_finite() || dphi0 >= 0.0 {
                direction = self.g.iter().map(|gradient| -gradient).collect();
                direction = project_parameter_direction(&self.manifold, &self.x, &direction);
                dphi0 = dot(&direction, &self.g);
            }
            if !dphi0.is_finite() || dphi0 >= 0.0 {
                // The configured stop algebra was checked above. If it did
                // not accept this state, a numerically unusable direction is
                // a stall; it must not mint an unrequested GradNorm reason.
                reason = StopReason::Stall;
                break;
            }

            let manifold = self.manifold;
            let origin = self.x.clone();
            let origin_gradient = self.g.clone();
            let remaining = hard_budget
                .map(|budget| budget.saturating_sub(self.evals))
                .unwrap_or(usize::MAX);
            let mut accepted = None;
            let mut trial_worst_violation = self.worst_violation;
            let outcome = {
                let mut curve = |alpha: f64| {
                    let raw_step: Vec<f64> = direction.iter().map(|value| alpha * value).collect();
                    let curve = manifold
                        .retract_curve(&origin, &direction, alpha)
                        .unwrap_or_else(|error| {
                            panic!("line-search curve failed fs-opt authority: {error}")
                        });
                    let trial = curve.point;
                    trial_worst_violation = trial_worst_violation.max(violation(&manifold, &trial));
                    let (f_trial, gradient) = checked_fg(&manifold, fg, &trial);
                    let derivative = dot(&gradient, &curve.velocity);
                    assert!(
                        derivative.is_finite(),
                        "retracted curve derivative must remain finite"
                    );
                    accepted = Some((alpha, raw_step, trial, f_trial, gradient));
                    (f_trial, derivative)
                };
                strong_wolfe_with_budget(&mut curve, self.f, dphi0, 1.0, 1e-4, 0.9, remaining)
            };
            // The invariant certificate covers every point exposed to the
            // objective, including rejected line-search trials.
            self.worst_violation = trial_worst_violation;
            self.evals = self
                .evals
                .checked_add(outcome.evals)
                .expect("evaluation accounting must not overflow");
            if !outcome.success {
                reason = if hard_budget.is_some() && outcome.evals == remaining {
                    StopReason::Budget
                } else {
                    StopReason::Stall
                };
                break;
            }
            let (alpha, raw_step, x_new, f_new, g_new) =
                accepted.expect("successful line search must retain its accepted trial");
            assert_eq!(
                alpha.to_bits(),
                outcome.alpha.to_bits(),
                "line-search outcome must name the retained trial"
            );
            assert_eq!(
                f_new.to_bits(),
                outcome.f_new.to_bits(),
                "line-search outcome value must match the retained trial"
            );

            let decision_noop = origin
                .iter()
                .zip(&x_new)
                .all(|(old, new)| old.to_bits() == new.to_bits());
            if decision_noop {
                let observation = StopObservation {
                    grad_norm: inf_norm(&self.g),
                    objective: self.f,
                    evals: self.evals,
                    history: &self.history,
                };
                reason = stop_reason(rule, hard_budget, &observation).unwrap_or(StopReason::Stall);
                break;
            }

            // Every retained pair moves to the new tangent space after every
            // accepted iterate, even when the new secant pair is rejected.
            transport_pairs(&manifold, &origin, &raw_step, &x_new, &mut self.pairs);
            let transported_old_gradient =
                vector_transport(&manifold, &origin, &raw_step, &x_new, &origin_gradient);
            let s = vector_transport(&manifold, &origin, &raw_step, &x_new, &raw_step);
            let y: Vec<f64> = g_new
                .iter()
                .zip(&transported_old_gradient)
                .map(|(new, old)| new - old)
                .collect();
            if self.memory == 0 {
                self.pairs.clear();
            } else if let Some(rho) = curvature_rho(&s, &y) {
                while self.pairs.len() >= self.memory {
                    self.pairs.pop_front();
                }
                self.pairs.push_back((s, y, rho));
            }

            self.worst_violation = self.worst_violation.max(violation(&manifold, &x_new));
            self.x = x_new;
            self.f = f_new;
            self.g = g_new;
            self.iters = self
                .iters
                .checked_add(1)
                .expect("iteration accounting must not overflow");
            self.history.push(self.f);
            self.assert_valid_state();

            let observation = StopObservation {
                grad_norm: inf_norm(&self.g),
                objective: self.f,
                evals: self.evals,
                history: &self.history,
            };
            if let Some(stopped) = stop_reason(rule, hard_budget, &observation) {
                reason = stopped;
                break 'iterations;
            }
        }
        self.report(reason)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::panic::{AssertUnwindSafe, catch_unwind};

    #[test]
    fn sphere_transport_preserves_tangency_inner_products_and_round_trip() {
        let manifold = Manifold::Sphere { ambient: 3 };
        let from = [1.0, 0.0, 0.0];
        let step = [0.0, 0.5, 0.0];
        let to = retract(&manifold, &from, &step);
        let first = [0.0, 1.0, 0.0];
        let second = [0.0, 0.0, 1.0];
        let moved_first = vector_transport(&manifold, &from, &step, &to, &first);
        let moved_second = vector_transport(&manifold, &from, &step, &to, &second);
        assert!(dot(&to, &moved_first).abs() <= 8.0 * f64::EPSILON);
        assert!(dot(&to, &moved_second).abs() <= 8.0 * f64::EPSILON);
        assert!(
            (dot(&first, &second) - dot(&moved_first, &moved_second)).abs() <= 8.0 * f64::EPSILON
        );
        assert!((l2_norm(&first) - l2_norm(&moved_first)).abs() <= 8.0 * f64::EPSILON);
        let scale = dot(&to, &from);
        let reverse_step: Vec<f64> = from
            .iter()
            .zip(&to)
            .map(|(target, current)| target / scale - current)
            .collect();
        let returned = retract(&manifold, &to, &reverse_step);
        let round_trip = vector_transport(&manifold, &to, &reverse_step, &returned, &moved_first);
        assert!(
            from.iter()
                .zip(&returned)
                .all(|(expected, actual)| (expected - actual).abs() <= 8.0 * f64::EPSILON)
        );
        assert!(
            first
                .iter()
                .zip(&round_trip)
                .all(|(expected, actual)| (expected - actual).abs() <= 8.0 * f64::EPSILON)
        );
    }

    #[test]
    fn transported_memory_recomputes_rho_and_drops_invalid_curvature() {
        let manifold = Manifold::Sphere { ambient: 3 };
        let from = [1.0, 0.0, 0.0];
        let step = [0.0, 0.5, 0.0];
        let to = retract(&manifold, &from, &step);
        let mut pairs = VecDeque::from([
            (vec![0.0, 1.0, 0.0], vec![0.0, 2.0, 0.0], 99.0),
            (vec![0.0, 0.0, 1.0], vec![0.0, 0.0, -1.0], 1.0),
        ]);
        transport_pairs(&manifold, &from, &step, &to, &mut pairs);
        assert_eq!(pairs.len(), 1, "negative transported curvature must drop");
        let (s, y, rho) = pairs.front().expect("positive pair retained");
        assert_eq!(
            rho.to_bits(),
            curvature_rho(s, y)
                .expect("retained pair remains valid")
                .to_bits()
        );
        assert_ne!(rho.to_bits(), 99.0f64.to_bits());
    }

    #[test]
    fn so3_wrappers_keep_quaternion_points_and_body_parameters_distinct() {
        let manifold = Manifold::So3;
        let point = [1.0, 0.0, 0.0, 0.0];
        let ambient_gradient = [0.0, 2.0, 0.0, 0.0];
        let gradient = tangent_project(&manifold, &point, &ambient_gradient);
        assert_eq!(gradient, vec![1.0, 0.0, 0.0]);

        let step = [0.2, -0.1, 0.3];
        let landed = retract(&manifold, &point, &step);
        assert_eq!(landed.len(), 4);
        let transported = vector_transport(&manifold, &point, &step, &landed, &gradient);
        assert_eq!(transported, gradient);
        assert_tangent(&manifold, &landed, &transported);
    }

    #[test]
    fn stiefel_wrappers_use_qr_retraction_and_its_differential_transport() {
        let manifold = Manifold::Stiefel { n: 4, p: 2 };
        let point = [
            0.5, 0.5, 0.5, 0.5, // first column
            0.5, -0.5, 0.5, -0.5, // second column
        ];
        let ambient = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75];
        let tangent = tangent_project(&manifold, &point, &ambient);
        assert_tangent(&manifold, &point, &tangent);
        let step: Vec<f64> = tangent.iter().map(|value| 0.125 * value).collect();
        let landed = retract(&manifold, &point, &step);
        let transported = vector_transport(&manifold, &point, &step, &landed, &tangent);
        let authoritative = manifold
            .transport_parameter(&point, &step, &landed, &tangent)
            .expect("valid Stiefel transport");
        assert_eq!(transported, authoritative);
        assert!(violation(&manifold, &landed) <= 2.0e-15);
        assert_tangent(&manifold, &landed, &transported);
    }

    #[test]
    fn curvature_admission_is_invariant_to_reciprocal_unit_rescaling() {
        let step = [1.0, 2.0, -0.5];
        let gradient_delta = [3.0, 4.0, 1.0];
        let step_rescaled = step.map(|value| value * 1.0e100);
        let gradient_rescaled = gradient_delta.map(|value| value * 1.0e-100);
        let original = curvature_rho(&step, &gradient_delta).expect("positive curvature");
        let rescaled = curvature_rho(&step_rescaled, &gradient_rescaled)
            .expect("reciprocal unit rescaling must preserve admission");
        assert!(
            (original - rescaled).abs() <= 16.0 * f64::EPSILON * original.abs().max(1.0),
            "reciprocal unit rescaling changed curvature rho: {original:e} versus {rescaled:e}"
        );
    }

    #[test]
    fn zero_memory_is_exact_steepest_descent_without_a_hidden_pair() {
        let mut objective = |x: &[f64]| (0.5 * x[0] * x[0], vec![x[0]]);
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[1.0], 0, &mut objective);
        let report = state.run(&mut objective, &StopRule::GradNorm(0.0), 1);
        assert_eq!(report.reason, StopReason::GradNorm);
        assert_eq!(report.evals, 2);
        assert!(state.pairs.is_empty());
    }

    #[test]
    fn nested_budget_is_a_hard_cap_inside_line_search() {
        let calls = Cell::new(0usize);
        let mut objective = |x: &[f64]| {
            calls.set(calls.get() + 1);
            (0.0, vec![x[0]])
        };
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[1.0], 4, &mut objective);
        let rule = StopRule::All(vec![StopRule::Budget(3), StopRule::GradNorm(0.0)]);
        let report = state.run(&mut objective, &rule, 100);
        assert_eq!(report.reason, StopReason::Budget);
        assert_eq!(report.evals, 3);
        assert_eq!(calls.get(), 3);
        assert_eq!(
            report.iters, 0,
            "failed search must not partially mutate state"
        );
    }

    #[test]
    fn absolute_budget_cannot_underwrite_existing_state() {
        let calls = Cell::new(0usize);
        let mut objective = |x: &[f64]| {
            calls.set(calls.get() + 1);
            (0.5 * x[0] * x[0], vec![x[0]])
        };
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[1.0], 4, &mut objective);
        assert_eq!(calls.get(), 1, "construction retains its initial callback");

        let refusal = catch_unwind(AssertUnwindSafe(|| {
            state.run(&mut objective, &StopRule::Budget(0), 10)
        }));
        assert!(
            refusal.is_err(),
            "an absolute budget below already-spent state must refuse"
        );
        assert_eq!(calls.get(), 1, "refusal must occur before another callback");
        assert_eq!(state.evals, 1);
        assert_eq!(state.iters, 0);
        assert_eq!(state.history, vec![0.5]);
    }

    #[test]
    fn accepted_precision_noop_stalls_without_state_transition() {
        let initial = Cell::new(true);
        let mut objective = |_x: &[f64]| {
            let gradient = if initial.replace(false) { 1.0 } else { 0.0 };
            (1.0e308, vec![gradient])
        };
        let mut state =
            RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[1.0e308], 4, &mut objective);
        let report = state.run(
            &mut objective,
            &StopRule::Any(vec![StopRule::Budget(100), StopRule::GradNorm(0.0)]),
            10,
        );
        assert_eq!(report.reason, StopReason::Stall);
        assert_eq!(report.grad_norm.to_bits(), 1.0f64.to_bits());
        assert_eq!(report.iters, 0);
        assert_eq!(state.history.len(), 1);
        assert_eq!(state.x[0].to_bits(), 1.0e308f64.to_bits());
    }

    #[test]
    fn invalid_initial_dimension_is_rejected_before_callback() {
        let calls = Cell::new(0usize);
        let result = catch_unwind(AssertUnwindSafe(|| {
            let mut objective = |_x: &[f64]| {
                calls.set(calls.get() + 1);
                (0.0, vec![0.0, 0.0])
            };
            let _ = RiemannianLbfgs::new(
                Manifold::Sphere { ambient: 3 },
                &[1.0, 0.0],
                4,
                &mut objective,
            );
        }));
        assert!(result.is_err());
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn malformed_callback_cannot_mint_a_zero_norm_certificate() {
        let result = catch_unwind(AssertUnwindSafe(|| {
            let mut objective = |_x: &[f64]| (0.0, vec![f64::NAN, 0.0]);
            let _ = RiemannianLbfgs::new(
                Manifold::Sphere { ambient: 2 },
                &[1.0, 0.0],
                4,
                &mut objective,
            );
        }));
        assert!(result.is_err());
    }

    #[test]
    fn underflowed_slope_does_not_fabricate_gradient_convergence() {
        let mut objective = |_x: &[f64]| (0.0, vec![1.0e-200]);
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[0.0], 4, &mut objective);
        let report = state.run(&mut objective, &StopRule::Budget(100), 10);
        assert_eq!(report.reason, StopReason::Stall);
        assert_eq!(report.grad_norm.to_bits(), 1.0e-200f64.to_bits());
        assert_eq!(report.grad_l2_norm.to_bits(), 1.0e-200f64.to_bits());
        assert_eq!(report.evals, 1);
    }

    #[test]
    fn stationary_state_does_not_bypass_the_configured_stop_algebra() {
        let mut objective = |_x: &[f64]| (0.0, vec![0.0]);
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[0.0], 4, &mut objective);
        let report = state.run(&mut objective, &StopRule::Budget(10), 10);
        assert_eq!(report.reason, StopReason::Stall);
        assert_eq!(report.evals, 1);
    }

    #[test]
    fn malformed_stop_rule_is_refused_before_an_additional_callback() {
        let calls = Cell::new(0usize);
        let mut objective = |_x: &[f64]| {
            calls.set(calls.get() + 1);
            (0.0, vec![1.0])
        };
        let mut state = RiemannianLbfgs::new(Manifold::Rn { dim: 1 }, &[0.0], 4, &mut objective);
        let result = catch_unwind(AssertUnwindSafe(|| {
            let _ = state.run(&mut objective, &StopRule::Any(Vec::new()), 10);
        }));
        assert!(result.is_err());
        assert_eq!(calls.get(), 1, "only construction may evaluate");
    }
}

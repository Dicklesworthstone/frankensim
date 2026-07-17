//! Riemannian L-BFGS over fs-opt's manifold METADATA: ambient
//! gradient → tangent projection → strong-Wolfe search along RETRACTED
//! curves, with memory pairs moved by deterministic isometric transport.
//! fs-opt's `Manifold` is metadata (dimensions, kinds); the
//! OPERATIONS (projection, retraction) live here with the engines
//! that consume them. Manifold invariants (unit norms) hold to
//! roundoff along the whole iterate path — tested, not assumed.

use crate::stop::{StopObservation, StopReason, StopRule};
use crate::wolfe::strong_wolfe_with_budget;
use fs_opt::Manifold;
use std::collections::VecDeque;

type CurvaturePair = (Vec<f64>, Vec<f64>, f64);

const MANIFOLD_ADMISSION_TOLERANCE: f64 = 1e-12;
const TANGENCY_TOLERANCE: f64 = 1e-10;
const CURVATURE_RELATIVE_FLOOR: f64 = 1e-14;

fn point_dim(man: &Manifold) -> usize {
    let raw = match *man {
        Manifold::Rn { dim } => {
            assert!(dim >= 1, "Rn dimension must be at least one");
            dim
        }
        Manifold::Sphere { ambient } => {
            assert!(
                ambient >= 2,
                "sphere ambient dimension must be at least two"
            );
            ambient
        }
        Manifold::So3 => 4,
        Manifold::Stiefel { .. } => {
            panic!("Stiefel is metadata-only until its consumer bead supplies the retraction")
        }
    };
    usize::try_from(raw).expect("manifold point dimension must fit usize")
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
        _ => 0.0,
    }
}

fn assert_point(man: &Manifold, x: &[f64]) {
    assert_eq!(
        x.len(),
        point_dim(man),
        "manifold point dimension does not match its descriptor"
    );
    assert_finite("manifold point", x);
    assert!(
        violation(man, x) <= MANIFOLD_ADMISSION_TOLERANCE,
        "sphere/SO(3) point must have unit norm within {MANIFOLD_ADMISSION_TOLERANCE:e}"
    );
}

fn assert_tangent(man: &Manifold, x: &[f64], vector: &[f64]) {
    assert_eq!(
        vector.len(),
        x.len(),
        "tangent-vector dimension must match the manifold point"
    );
    assert_finite("tangent vector", vector);
    if matches!(man, Manifold::Sphere { .. } | Manifold::So3) {
        let residual = dot(x, vector).abs();
        let scale = l2_norm(vector).max(f64::MIN_POSITIVE);
        assert!(
            residual <= TANGENCY_TOLERANCE * scale,
            "sphere/SO(3) vector is not tangent: residual {residual:e}"
        );
    }
}

fn checked_fg(fg: crate::FnGrad<'_>, x: &[f64]) -> (f64, Vec<f64>) {
    let (f, gradient) = fg(x);
    assert!(f.is_finite(), "objective value must be finite");
    assert_eq!(
        gradient.len(),
        x.len(),
        "objective gradient dimension must match the decision vector"
    );
    assert_finite("objective gradient", &gradient);
    (f, gradient)
}

/// Project an ambient vector onto the tangent space at `x`.
///
/// # Panics
/// For metadata-only manifolds (Stiefel — its retraction lands with
/// its consumer bead; using it in a descent today is a modeling
/// error, surfaced loudly).
#[must_use]
pub fn tangent_project(man: &Manifold, x: &[f64], v: &[f64]) -> Vec<f64> {
    assert_point(man, x);
    assert_eq!(
        v.len(),
        x.len(),
        "ambient-vector dimension must match the manifold point"
    );
    assert_finite("ambient vector", v);
    let projected = match man {
        Manifold::Rn { .. } => v.to_vec(),
        Manifold::Sphere { .. } | Manifold::So3 => {
            // v − (v·x)·x.
            let normal_component = dot(v, x);
            v.iter()
                .zip(x)
                .map(|(vi, xi)| (-normal_component).mul_add(*xi, *vi))
                .collect()
        }
        Manifold::Stiefel { .. } => unreachable!("point admission rejects Stiefel"),
    };
    assert_finite("projected tangent vector", &projected);
    projected
}

/// Retract from `x` along `step` (metric-projection retractions).
///
/// # Panics
/// Same metadata-only policy as [`tangent_project`].
#[must_use]
pub fn retract(man: &Manifold, x: &[f64], step: &[f64]) -> Vec<f64> {
    assert_point(man, x);
    assert_eq!(
        step.len(),
        x.len(),
        "retraction-step dimension must match the manifold point"
    );
    assert_finite("retraction step", step);
    let retracted: Vec<f64> = match man {
        Manifold::Rn { .. } => x.iter().zip(step).map(|(a, b)| a + b).collect(),
        Manifold::Sphere { .. } | Manifold::So3 => {
            let moved: Vec<f64> = x.iter().zip(step).map(|(a, b)| a + b).collect();
            assert_finite("retraction trial", &moved);
            let nrm = l2_norm(&moved);
            assert!(nrm > 0.0, "retraction trial must have nonzero norm");
            moved.iter().map(|t| t / nrm).collect()
        }
        Manifold::Stiefel { .. } => unreachable!("point admission rejects Stiefel"),
    };
    assert_point(man, &retracted);
    retracted
}

/// Isometric vector transport on the sphere/SO(3) shortest geodesic.
/// Rn transport is the identity. The normalized retraction of a tangent step
/// cannot reach the antipode, but the denominator is still checked so corrupt
/// resumed state fails closed instead of contaminating curvature memory.
fn vector_transport(man: &Manifold, from: &[f64], to: &[f64], vector: &[f64]) -> Vec<f64> {
    assert_point(man, from);
    assert_point(man, to);
    let vector = tangent_project(man, from, vector);
    let transported = match man {
        Manifold::Rn { .. } => vector,
        Manifold::Sphere { .. } | Manifold::So3 => {
            let denominator = 1.0 + dot(from, to);
            assert!(
                denominator.is_finite() && denominator > 64.0 * f64::EPSILON,
                "sphere/SO(3) vector transport is undefined at antipodal points"
            );
            let coefficient = dot(&vector, to) / denominator;
            vector
                .iter()
                .zip(from.iter().zip(to))
                .map(|(value, (old, new))| (-coefficient).mul_add(old + new, *value))
                .collect()
        }
        Manifold::Stiefel { .. } => unreachable!("point admission rejects Stiefel"),
    };
    // Remove the final few ulps of normal component introduced by arithmetic.
    tangent_project(man, to, &transported)
}

fn retraction_velocity(
    man: &Manifold,
    origin: &[f64],
    direction: &[f64],
    alpha: f64,
    trial: &[f64],
) -> Vec<f64> {
    match man {
        Manifold::Rn { .. } => direction.to_vec(),
        Manifold::Sphere { .. } | Manifold::So3 => {
            let moved: Vec<f64> = origin
                .iter()
                .zip(direction)
                .map(|(x, d)| alpha.mul_add(*d, *x))
                .collect();
            let scale = l2_norm(&moved);
            tangent_project(man, trial, direction)
                .into_iter()
                .map(|value| value / scale)
                .collect()
        }
        Manifold::Stiefel { .. } => unreachable!("point admission rejects Stiefel"),
    }
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

fn transport_pairs(man: &Manifold, from: &[f64], to: &[f64], pairs: &mut VecDeque<CurvaturePair>) {
    let mut transported = VecDeque::with_capacity(pairs.len());
    while let Some((s, y, _old_rho)) = pairs.pop_front() {
        let s = vector_transport(man, from, to, &s);
        let y = vector_transport(man, from, to, &y);
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
        assert_point(&manifold, x0);
        let (f, g_amb) = checked_fg(fg, x0);
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

    fn assert_valid_state(&self) {
        assert_point(&self.manifold, &self.x);
        assert!(self.f.is_finite(), "current objective must be finite");
        assert_eq!(
            self.g.len(),
            self.x.len(),
            "current gradient dimension must match the decision vector"
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
    /// below `All`: no line-search phase may overshoot it. Sphere/SO(3)
    /// derivatives use the exact velocity of the normalized retraction.
    pub fn run(
        &mut self,
        fg: crate::FnGrad<'_>,
        rule: &StopRule,
        max_iters: usize,
    ) -> RiemannianReport {
        assert_valid_stop_rule(rule);
        self.assert_valid_state();
        let hard_budget = evaluation_budget(rule);
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
                tangent_project(&self.manifold, &self.x, &raw_direction)
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
                direction = tangent_project(&self.manifold, &self.x, &direction);
                dphi0 = -dot(&self.g, &self.g);
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
                    let step: Vec<f64> = direction.iter().map(|value| alpha * value).collect();
                    let trial = retract(&manifold, &origin, &step);
                    trial_worst_violation = trial_worst_violation.max(violation(&manifold, &trial));
                    let (f_trial, ambient_gradient) = checked_fg(fg, &trial);
                    let gradient = tangent_project(&manifold, &trial, &ambient_gradient);
                    let velocity =
                        retraction_velocity(&manifold, &origin, &direction, alpha, &trial);
                    let derivative = dot(&gradient, &velocity);
                    assert!(
                        derivative.is_finite(),
                        "retracted curve derivative must remain finite"
                    );
                    accepted = Some((alpha, trial, f_trial, gradient));
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
            let (alpha, x_new, f_new, g_new) =
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
            transport_pairs(&manifold, &origin, &x_new, &mut self.pairs);
            let transported_old_gradient =
                vector_transport(&manifold, &origin, &x_new, &origin_gradient);
            let raw_step: Vec<f64> = direction.iter().map(|value| alpha * value).collect();
            let s = vector_transport(&manifold, &origin, &x_new, &raw_step);
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
        let to = [0.0, 1.0, 0.0];
        let first = [0.0, 1.0, 0.0];
        let second = [0.0, 0.0, 1.0];
        let moved_first = vector_transport(&manifold, &from, &to, &first);
        let moved_second = vector_transport(&manifold, &from, &to, &second);
        assert!(dot(&to, &moved_first).abs() <= 8.0 * f64::EPSILON);
        assert!(dot(&to, &moved_second).abs() <= 8.0 * f64::EPSILON);
        assert!(
            (dot(&first, &second) - dot(&moved_first, &moved_second)).abs() <= 8.0 * f64::EPSILON
        );
        assert!((l2_norm(&first) - l2_norm(&moved_first)).abs() <= 8.0 * f64::EPSILON);
        let round_trip = vector_transport(&manifold, &to, &from, &moved_first);
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
        let to = [0.0, 1.0, 0.0];
        let mut pairs = VecDeque::from([
            (vec![0.0, 1.0, 0.0], vec![0.0, 2.0, 0.0], 99.0),
            (vec![0.0, 0.0, 1.0], vec![0.0, 0.0, -1.0], 1.0),
        ]);
        transport_pairs(&manifold, &from, &to, &mut pairs);
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

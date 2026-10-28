//! Local extended-Kalman prediction with caller-scaled vector-model probes.
use super::{PredictionMean, invalid, predict_prefix_factors, prediction_work};
use crate::factor::{FactorBelief, try_f64_vec};
use crate::{
    AssimError, Cx, WorkPlan, WorkProgress, canonical_f64_bits, canonicalize_zero,
    validate_leaf_identity, validate_state_dimension,
};
use fs_blake3::{ContentHash, DomainHasher};

/// A discrete nonlinear physical transition with independent additive noise.
///
/// Each positive `steps[i]` is an ABSOLUTE central-difference perturbation in
/// the caller's units for state component i. The callback implements the entire
/// time interval; neither it nor its Jacobian is multiplied by dt here.
/// Covariance uses the local approximation J P J^T + Q at the prior mean.
/// This is an extended-Kalman approximation, not a nonlinear moment integral,
/// automatic differentiation, a derivative certificate, or a flow error bound.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearDynamics {
    noise: FactorBelief,
    steps: Vec<f64>,
    identity: String,
    content_id: ContentHash,
}

impl NonlinearDynamics {
    /// Admit zero-mean process noise and one finite positive step per physical
    /// component. Singular process noise needs no jitter. Its independence
    /// from the physical and nuisance priors is a caller model assumption.
    pub fn new(
        noise: FactorBelief,
        mut steps: Vec<f64>,
        identity: impl Into<String>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        let n = noise.dim();
        validate_state_dimension(n)?;
        if steps.len() != n || steps.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return Err(invalid("nonlinear dynamics require one positive finite step per component"));
        }
        if noise.mean().iter().any(|m| *m != 0.0) {
            return Err(invalid("nonlinear process noise must have zero mean"));
        }
        let identity = identity.into();
        validate_leaf_identity("nonlinear dynamics", &identity)?;
        let nu = n as u128;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(
            0, 0, 0, 0, 4 * nu * nu + 8 * nu, 0,
            8 * nu * nu + 32 * nu + identity.len() as u128 + 128,
        )?)?;
        progress.checkpoint("nonlinear dynamics admission")?;
        let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.nonlinear-dynamics.v1");
        hash.update(&(identity.len() as u64).to_le_bytes());
        hash.update(identity.as_bytes());
        hash.update(&(n as u64).to_le_bytes());
        progress.hash_bytes("nonlinear dynamics identity", identity.len() as u128 + 16)?;
        for value in &mut steps { *value = canonicalize_zero(*value); }
        for value in steps.iter().chain(&noise.upper).chain(&noise.diag) {
            hash_number(&mut hash, *value, &mut progress)?;
        }
        progress.checkpoint("nonlinear dynamics publication")?;
        Ok(Self { noise, steps, identity, content_id: hash.finalize() })
    }

    /// Number of physical state variables.
    #[must_use]
    pub fn dimension(&self) -> usize { self.steps.len() }

    /// Absolute requested perturbations, ordered like the physical state.
    #[must_use]
    pub fn steps(&self) -> &[f64] { &self.steps }

    /// Caller-declared model identity; callback code is not authenticated.
    #[must_use]
    pub fn identity(&self) -> &str { &self.identity }

    /// Binds the declaration, finite-difference policy, and process noise.
    #[must_use]
    pub const fn content_id(&self) -> ContentHash { self.content_id }

    /// Predict mean f(mean) and factor covariance J P J^T + Q.
    ///
    /// Exactly 1+2*n callback evaluations on success: the mean, then plus/minus
    /// probes for each component in order. Every call receives the same Cx for
    /// finer cancellation and model-owned budgeting. Callback cost itself is
    /// not inferred; the wrapper charges its numerical/hash work and all probe
    /// checkpoints cumulatively with the square-root propagation.
    ///
    /// The callback must be deterministic for its declared identity and return
    /// n finite values or an explicit domain error. No failed probe is skipped.
    /// The unmodified prior remains available after every refusal. Finite
    /// differences do not certify smoothness or detect every poor step choice.
    pub fn predict_fd<F>(
        &self,
        prior: &FactorBelief,
        predict: F,
        cx: &Cx<'_>,
    ) -> Result<FactorBelief, AssimError>
    where
        F: FnMut(&[f64], &Cx<'_>) -> Result<Vec<f64>, AssimError>,
    {
        if prior.dim() != self.dimension() {
            return Err(AssimError::DimMismatch { state: prior.dim(), operator: self.dimension() });
        }
        self.predict_prefix_fd(prior, predict, cx).map(|(belief, _)| belief)
    }

    pub(crate) fn predict_prefix_fd<F>(
        &self,
        prior: &FactorBelief,
        mut predict: F,
        cx: &Cx<'_>,
    ) -> Result<(FactorBelief, ContentHash), AssimError>
    where
        F: FnMut(&[f64], &Cx<'_>) -> Result<Vec<f64>, AssimError>,
    {
        let n = self.dimension();
        let d = prior.dim();
        validate_state_dimension(d)?;
        if d < n {
            return Err(AssimError::DimMismatch { state: d, operator: n });
        }
        let nu = n as u128;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(
            0, 0, 0, 0, prediction_work(n, d) + 64 * nu * nu + 64 * nu + 256,
            0, 32 * nu * nu + 64 * nu + 128,
        )?)?;
        progress.checkpoint("nonlinear prediction admission")?;
        let mut point = try_f64_vec(n, "nonlinear prediction point", &mut progress)?;
        point.copy_from_slice(&prior.mean[..n]);
        let mut plus = try_f64_vec(n, "nonlinear plus states", &mut progress)?;
        let mut minus = try_f64_vec(n, "nonlinear minus states", &mut progress)?;
        for i in 0..n {
            plus[i] = point[i] + self.steps[i];
            minus[i] = point[i] - self.steps[i];
            let span = plus[i] - minus[i];
            if !plus[i].is_finite() || !minus[i].is_finite()
                || plus[i] <= point[i] || minus[i] >= point[i]
                || !span.is_finite() || span <= 0.0
            {
                return Err(AssimError::FiniteDifferenceStepUnrepresentable { component: i });
            }
            progress.scalar("nonlinear perturbations", 8)?;
        }
        let mut jacobian = Vec::new();
        jacobian.try_reserve_exact(n).map_err(|_| AssimError::AllocationRefused {
            stage: "nonlinear dynamics Jacobian", requested_bytes: nu * nu * 8,
        })?;
        for _ in 0..n {
            jacobian.push(try_f64_vec(n, "nonlinear dynamics Jacobian", &mut progress)?);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.nonlinear-forecast.v1");
        hash.update(self.content_id.as_bytes());
        progress.hash_bytes("nonlinear prediction identity", 32)?;
        for value in &point { hash_number(&mut hash, *value, &mut progress)?; }
        let value = evaluate(&mut predict, &point, n, cx, &mut progress)?;
        for x in &value { hash_number(&mut hash, *x, &mut progress)?; }
        for column in 0..n {
            let original = point[column];
            point[column] = plus[column];
            let high = evaluate(&mut predict, &point, n, cx, &mut progress)?;
            point[column] = minus[column];
            let low = evaluate(&mut predict, &point, n, cx, &mut progress)?;
            point[column] = original;
            hash_number(&mut hash, plus[column], &mut progress)?;
            hash_number(&mut hash, minus[column], &mut progress)?;
            let span = plus[column] - minus[column];
            for row in 0..n {
                let derivative = (high[row] - low[row]) / span;
                if !derivative.is_finite() {
                    return Err(AssimError::NonFiniteComputation {
                        stage: "nonlinear dynamics Jacobian",
                    });
                }
                jacobian[row][column] = canonicalize_zero(derivative);
                progress.scalar("nonlinear dynamics Jacobian", 3)?;
                hash_number(&mut hash, high[row], &mut progress)?;
                hash_number(&mut hash, low[row], &mut progress)?;
                hash_number(&mut hash, derivative, &mut progress)?;
            }
        }
        let belief = predict_prefix_factors(
            &jacobian, PredictionMean::Evaluated(&value), &self.noise, prior, &mut progress,
        )?;
        progress.checkpoint("nonlinear prediction publication")?;
        Ok((belief, hash.finalize()))
    }
}

fn evaluate<F>(
    predict: &mut F,
    state: &[f64],
    n: usize,
    cx: &Cx<'_>,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Vec<f64>, AssimError>
where
    F: FnMut(&[f64], &Cx<'_>) -> Result<Vec<f64>, AssimError>,
{
    progress.checkpoint("nonlinear model evaluation")?;
    let mut value = predict(state, cx)?;
    progress.checkpoint("nonlinear model completion")?;
    if value.len() != n {
        return Err(invalid("nonlinear model output dimension differs from its physical state"));
    }
    for entry in &mut value {
        if !entry.is_finite() {
            return Err(AssimError::NonFiniteComputation { stage: "nonlinear dynamics output" });
        }
        *entry = canonicalize_zero(*entry);
        progress.scalar("nonlinear model output", 1)?;
    }
    Ok(value)
}

fn hash_number(
    hash: &mut DomainHasher,
    value: f64,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    hash.update(&canonical_f64_bits(value).to_le_bytes());
    progress.hash_bytes("nonlinear prediction identity", 8)
}

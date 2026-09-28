//! Gaussian-copula Monte Carlo with explicitly bounded uniform marginals.
//!
//! The existing joint-Gaussian producer owns PSD admission, Philox addressing,
//! cancellation, refusal, observations and restart. This adapter maps its
//! correlated standard normals through Phi, then through each declared uniform
//! inverse CDF. The supplied matrix describes the LATENT NORMALS, not Pearson
//! correlations of the physical uniform values. A matrix without this joint-law
//! declaration is still refused by the ordinary product executor.
//!
//! Phi uses the production deterministic erfc, including the lower tail rather
//! than subtracting an almost-one erf. Rounded endpoints are possible in extreme
//! tails. No clipping, redrawing, tail rejection or fitted correlation is used.
//! This is floating-point sampling of a declared model, not a distributional,
//! PSD, physical-validity or continuum-error certificate.

use core::fmt::{self, Display};

use fs_blake3::{ContentHash, DomainHasher};

use crate::{CorrelationModel, ParameterUncertainty, UncertaintyKind, UqCheckpointError,
    UqExecution, UqPlan, UqResult};

/// A joint sampler whose evaluator always receives physical uniform values.
///
/// The marginal plan must explicitly use `Independent`: it supplies marginal
/// laws and execution policy only, and the separately supplied copula matrix
/// replaces that independence. Another declared joint law is never discarded.
/// Between sample ordinals the existing Monte Carlo independence assumptions
/// remain unchanged; dependence is WITHIN each physical parameter vector.
#[derive(Debug, Clone)]
pub struct GaussianCopulaExecution {
    marginals: Vec<ParameterUncertainty>,
    execution: UqExecution,
}

impl GaussianCopulaExecution {
    /// Admit finite uniform supports and a latent-normal correlation matrix.
    ///
    /// Reuses the existing numerical PSD factor, including singular matrices.
    /// The matrix uses parameter declaration order and has a unit diagonal.
    /// Zero-width marginals remain exactly constant; their latent coordinates
    /// still participate in the declared matrix and consume ordinary draws.
    /// Only the existing Monte Carlo method is admitted, not QMC or intervals.
    /// No physical model is evaluated during admission.
    pub fn new(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>])
        -> Result<Self, &'static str>
    {
        let latent = latent_plan(marginal_plan, latent_correlation)?;
        let execution = UqExecution::new(&latent)?;
        Ok(Self { marginals: marginal_plan.parameters.clone(), execution })
    }

    /// Physical marginal declarations; their units are not the latent units.
    #[must_use]
    pub fn marginals(&self) -> &[ParameterUncertainty] { &self.marginals }

    /// Read-only Monte Carlo observations and inference machinery.
    ///
    /// Observations are actual physical QoIs, so compliance assessments act on
    /// the original threshold. Its `plan()` describes LATENT NORMAL coordinates,
    /// not physical inputs. Use this adapter's checkpoint methods: they also bind
    /// the physical supports, units and transform semantics.
    #[must_use]
    pub const fn monte_carlo(&self) -> &UqExecution { &self.execution }

    /// Physical-QoI report with the original execution/refusal semantics.
    #[must_use]
    pub fn report(&self) -> UqResult { self.execution.report() }

    /// At most this many further physical evaluations. Failed samples are terminal.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluator: F) -> UqResult
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |parameters| evaluator(parameters).map(Some))
    }

    /// Interrupted evaluations retry the SAME physical vector after restoration.
    /// `None` may only describe unfinished work, never a discarded inconvenient
    /// draw. Mapping or evaluator errors refuse the existing execution permanently.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, cancelled: C, mut evaluator: F,
    ) -> UqResult
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        let marginals = &self.marginals;
        self.execution.advance_interruptible(allowance, cancelled, |latent| {
            let physical = physical_parameters(marginals, latent).map_err(EvaluationError::Transform)?;
            evaluator(&physical).map_err(EvaluationError::Producer)
        })
    }

    /// Persist through the existing checkpoint owner, additionally binding every
    /// physical support and unit. `model` must bind the actual evaluator and all
    /// fixed physical inputs; an integrity checksum is not authentication.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        self.execution.checkpoint(self.bound_model(model))
    }

    /// Re-admit both the latent law and physical supports before restoring the
    /// observation prefix. Changing even only a physical endpoint or unit refuses.
    pub fn restore(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>],
        model: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError>
    {
        let mut restored = Self::new(marginal_plan, latent_correlation)
            .map_err(UqCheckpointError::InvalidPlan)?;
        restored.execution = UqExecution::restore(restored.execution.plan(), restored.bound_model(model), bytes)?;
        Ok(restored)
    }

    fn bound_model(&self, model: ContentHash) -> ContentHash {
        bound_model(&self.marginals, model)
    }
}

// Share admission and the exact transformation with randomized QMC. The
// method-specific owner still admits the method, PSD matrix and work policy.
pub(super) fn latent_plan(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>])
    -> Result<UqPlan, &'static str>
{
    let dim = marginal_plan.parameters.len();
    if dim == 0 || dim > 256 || latent_correlation.len() != dim
        || latent_correlation.iter().any(|row| row.len() != dim)
    { return Err("copula requires 1..=256 uniform marginals and a matching latent matrix"); }
    if !matches!(&marginal_plan.correlation, CorrelationModel::Independent) {
        return Err("copula construction requires an Independent marginal plan; another joint law cannot be overwritten");
    }
    for parameter in &marginal_plan.parameters {
        if parameter.unit.is_empty() { return Err("uniform marginals require explicit physical units"); }
        match parameter.kind {
            UncertaintyKind::AleatoryUniform { lo, hi }
                if lo.is_finite() && hi.is_finite() && lo <= hi => {},
            _ => return Err("Gaussian copula requires explicit finite uniform supports"),
        }
    }
    let mut latent = marginal_plan.clone();
    latent.parameters = marginal_plan.parameters.iter().map(|parameter|
        ParameterUncertainty::gaussian(&parameter.name, 0.0, 1.0, "1")).collect();
    latent.correlation = CorrelationModel::JointGaussian { matrix: latent_correlation.to_vec() };
    Ok(latent)
}

pub(super) fn bound_model(marginals: &[ParameterUncertainty], model: ContentHash) -> ContentHash {
    // The inner checkpoint already binds the exact latent matrix, parameter
    // order/names, sampler, QoI, seed, threshold and lifetime budget.
    let mut hash = DomainHasher::new("org.frankensim.uq.uniform-gaussian-copula.v1");
    hash.update(model.as_bytes());
    hash.update(&(marginals.len() as u64).to_le_bytes());
    for parameter in marginals {
        hash.update(&(parameter.unit.len() as u64).to_le_bytes());
        hash.update(parameter.unit.as_bytes());
        let UncertaintyKind::AleatoryUniform { lo, hi } = parameter.kind else {
            unreachable!("private admitted marginals");
        };
        hash.update(&lo.to_bits().to_le_bytes());
        hash.update(&hi.to_bits().to_le_bytes());
    }
    hash.finalize()
}

#[derive(Debug)]
pub(super) enum EvaluationError<E> { Transform(&'static str), Producer(E) }
impl<E: Display> Display for EvaluationError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self { Self::Transform(reason) => write!(f, "copula transform: {reason}"),
            Self::Producer(error) => error.fmt(f) }
    }
}

pub(super) fn physical_parameters(marginals: &[ParameterUncertainty], latent: &[f64]) -> Result<Vec<f64>, &'static str> {
    if marginals.len() != latent.len() { return Err("latent sample arity differs"); }
    marginals.iter().zip(latent).map(|(parameter, &z)| {
        if !z.is_finite() { return Err("nonfinite latent normal"); }
        let UncertaintyKind::AleatoryUniform { lo, hi } = parameter.kind else {
            return Err("expected uniform marginal");
        };
        // This operation tree is part of the checkpoint transform version.
        let u = 0.5 * fs_math::det::erfc(-z * core::f64::consts::FRAC_1_SQRT_2);
        if !(u.is_finite() && (0.0..=1.0).contains(&u)) { return Err("normal CDF outside [0,1]"); }
        let value = if lo == hi { lo } else { (1.0 - u) * lo + u * hi };
        if !(value.is_finite() && value >= lo && value <= hi) {
            return Err("mapped uniform value outside finite support");
        }
        Ok(value)
    }).collect()
}

#[cfg(test)]
mod tests;

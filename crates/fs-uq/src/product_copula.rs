//! Gaussian-copula Monte Carlo with declared uniform and empirical marginals.
//!
//! The existing joint-Gaussian producer owns PSD admission, Philox addressing,
//! cancellation, refusal, observations and restart. This adapter maps its
//! correlated standard normals through Phi, then through each declared inverse
//! CDF. The supplied matrix describes the LATENT NORMALS, not Pearson
//! correlations of the physical values. A matrix without this joint-law
//! declaration is still refused by the ordinary product executor.
//!
//! Phi uses the production deterministic erfc, including the lower tail rather
//! than subtracting an almost-one erf. Rounded endpoints are possible in extreme
//! tails. No clipping, redrawing, tail rejection or fitted correlation is used.
//! This is floating-point sampling of a declared model, not a distributional,
//! PSD, physical-validity or continuum-error certificate.

use core::fmt::{self, Display};

use fs_blake3::{ContentHash, DomainHasher};

use crate::{CorrelationModel, EmpiricalMarginal, ParameterUncertainty, UncertaintyKind, UqCheckpointError,
    UqExecution, UqPlan, UqResult};

/// Frozen physical-coordinate mean controls, distinct from latent controls.
pub mod mean_control;

/// A joint sampler whose evaluator receives physical values in declaration order.
///
/// The marginal plan must explicitly use `Independent`: it supplies marginal
/// laws and execution policy only, and the separately supplied copula matrix
/// replaces that independence. Another declared joint law is never discarded.
/// Between sample ordinals the existing Monte Carlo independence assumptions
/// remain unchanged; dependence is WITHIN each physical parameter vector.
#[derive(Debug, Clone)]
pub struct GaussianCopulaExecution {
    marginals: CopulaMarginals,
    execution: UqExecution,
}

impl GaussianCopulaExecution {
    /// Admit uniform or empirical marginals and a latent-normal correlation matrix.
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
        let (marginals, latent) = admit_marginals(marginal_plan, latent_correlation)?;
        let execution = UqExecution::new(&latent)?;
        Ok(Self { marginals, execution })
    }

    /// Physical marginal declarations; their units are not the latent units.
    #[must_use]
    pub fn marginals(&self) -> &[ParameterUncertainty] { &self.marginals.parameters }

    /// Empirical declarations in evaluator order, including source counts and
    /// observed ranges. Draws are confined to the observed values.
    pub fn empirical_marginals(&self) -> impl Iterator<Item = (&ParameterUncertainty, &EmpiricalMarginal)> {
        self.marginals.empirical()
    }

    /// Read-only Monte Carlo observations and inference machinery.
    ///
    /// Observations are actual physical QoIs, so compliance assessments act on
    /// the original threshold. Its `plan()` describes LATENT NORMAL coordinates,
    /// ordered by name for plans containing an empirical marginal. The adapter's
    /// checkpoint methods also bind physical supports, units and transform semantics.
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
            let physical = marginals.transform(latent).map_err(EvaluationError::Transform)?;
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
    /// observation prefix. Changes to a sample, unit or support refuse restoration.
    pub fn restore(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>],
        model: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError>
    {
        let (marginals, latent) = admit_marginals(marginal_plan, latent_correlation)
            .map_err(UqCheckpointError::InvalidPlan)?;
        let execution = UqExecution::restore(&latent, marginals.bound_model(model), bytes)?;
        Ok(Self { marginals, execution })
    }

    fn bound_model(&self, model: ContentHash) -> ContentHash {
        self.marginals.bound_model(model)
    }
}

/// The indices map physical declaration order to the latent sampler's coordinates.
#[derive(Debug, Clone)]
pub(super) struct CopulaMarginals {
    pub(super) parameters: Vec<ParameterUncertainty>,
    latent_indices: Vec<usize>,
}

// The method-specific owner admits the latent PSD matrix and execution policy.
pub(super) fn admit_marginals(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>])
    -> Result<(CopulaMarginals, UqPlan), &'static str>
{
    let parameters = &marginal_plan.parameters;
    let dim = parameters.len();
    if dim == 0 || dim > 256 || latent_correlation.len() != dim
        || latent_correlation.iter().any(|row| row.len() != dim)
    { return Err("copula requires 1..=256 marginals and a matching latent matrix"); }
    if !matches!(&marginal_plan.correlation, CorrelationModel::Independent) {
        return Err("copula construction requires an Independent marginal plan; another joint law cannot be overwritten");
    }
    let has_empirical = parameters.iter().any(|p| matches!(p.kind, UncertaintyKind::AleatoryEmpirical(_)));
    for parameter in parameters {
        if parameter.unit.is_empty() || (has_empirical && parameter.unit.trim().is_empty()) {
            return Err("copula marginals require explicit physical units");
        }
        match parameter.kind {
            UncertaintyKind::AleatoryUniform { lo, hi }
                if lo.is_finite() && hi.is_finite() && lo <= hi => {},
            UncertaintyKind::AleatoryEmpirical(_) => {},
            _ => return Err("Gaussian copula requires a finite uniform support or a declared empirical sample"),
        }
    }
    let mut latent = marginal_plan.clone();
    // Assign the same random stream and factorization to each named input when
    // the caller permutes both declarations and matrix. Uniform-only v1 draws
    // used declaration order and must retain that addressing for replay.
    let mut order: Vec<_> = (0..dim).collect();
    if has_empirical {
        order.sort_unstable_by(|&i, &j| parameters[i].name.cmp(&parameters[j].name));
    }
    latent.parameters = order.iter()
        .map(|&i| ParameterUncertainty::gaussian(&parameters[i].name, 0.0, 1.0, "1"))
        .collect();
    let matrix = order.iter()
        .map(|&i| order.iter().map(|&j| latent_correlation[i][j]).collect())
        .collect();
    latent.correlation = CorrelationModel::JointGaussian { matrix };
    let mut latent_indices = vec![0; dim];
    for (latent_index, &physical_index) in order.iter().enumerate() {
        latent_indices[physical_index] = latent_index;
    }
    Ok((CopulaMarginals { parameters: parameters.clone(), latent_indices }, latent))
}

impl CopulaMarginals {
    pub(super) fn empirical(&self) -> impl Iterator<Item = (&ParameterUncertainty, &EmpiricalMarginal)> {
        self.parameters.iter().filter_map(|p| match &p.kind {
            UncertaintyKind::AleatoryEmpirical(m) => Some((p, m)),
            _ => None,
        })
    }

    pub(super) fn bound_model(&self, model: ContentHash) -> ContentHash {
        if self.empirical().next().is_some() {
            // Bind physical declaration order independently of the canonical latent order.
            let mut bytes = model.as_bytes().to_vec();
            bytes.extend_from_slice(&(self.parameters.len() as u64).to_le_bytes());
            for parameter in &self.parameters {
                crate::product_checkpoint::encode_parameter(&mut bytes, parameter);
            }
            return fs_blake3::hash_domain("org.frankensim.uq.empirical-gaussian-copula.v1", &bytes);
        }
        uniform_bound_model(&self.parameters, model)
    }

    pub(super) fn transform(&self, latent: &[f64]) -> Result<Vec<f64>, &'static str> {
        if self.parameters.len() != latent.len() { return Err("latent sample arity differs"); }
        self.parameters.iter().zip(&self.latent_indices)
            .map(|(parameter, &i)| physical_value(parameter, latent[i])).collect()
    }
}

fn uniform_bound_model(marginals: &[ParameterUncertainty], model: ContentHash) -> ContentHash {
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

fn physical_value(parameter: &ParameterUncertainty, z: f64) -> Result<f64, &'static str> {
    if !z.is_finite() { return Err("nonfinite latent normal"); }
    // This operation tree is part of the checkpoint transform version.
    let u = 0.5 * fs_math::det::erfc(-z * core::f64::consts::FRAC_1_SQRT_2);
    if !(u.is_finite() && (0.0..=1.0).contains(&u)) { return Err("normal CDF outside [0,1]"); }
    match parameter.kind {
        UncertaintyKind::AleatoryUniform { lo, hi } => {
            let value = if lo == hi { lo } else { (1.0 - u) * lo + u * hi };
            if !(value.is_finite() && value >= lo && value <= hi) {
                return Err("mapped uniform value outside finite support");
            }
            Ok(value)
        }
        UncertaintyKind::AleatoryEmpirical(ref marginal) => Ok(marginal.quantile(u)),
        _ => Err("unsupported copula marginal"),
    }
}

#[cfg(test)]
mod tests;

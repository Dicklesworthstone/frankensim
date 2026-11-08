//! Joint passive calibration on authored tensor directions and relaxation times.
//!
//! Each nonnegative fitted strength multiplies an original dimensionless Gram
//! tensor `B^T B`. Equilibrium/Maxwell strengths have units Pa; Kelvin strengths
//! have units Pa s. Signed normal/shear couplings are retained. Every experiment
//! starts with fully relaxed Maxwell memory and zero strain rate at its first
//! observation, including at nonzero initial strain. Subsequent stresses are
//! the left-hand limits at the end of the incoming linear strain ramps, which
//! matters when Kelvin strain rates change at an observation.
//!
//! This fits the supplied tensor cone, not arbitrary tensor directions,
//! relaxation times, hidden prehistory or noise covariance. It does not prove
//! parameter uniqueness, predict an experimental validity band, or validate
//! unsampled histories. A fit is published only after the original admitted
//! constitutive history independently reproduces every authored stress bound.

use super::{
    TensorProny, TensorPronyBranch, TensorPronyError, TensorPronyHistory, WEIGHTS, apply, dot,
    gram, mandel, transpose_add,
};
use crate::visco::{PronyRamp, ViscoError, nnls};
use crate::{Tangent6, Voigt};
use fs_math::det;
use std::sync::Arc;

/// A caller-authored tensor direction with one fitted nonnegative strength.
/// Gram rows use Mandel coordinates `[xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx]`.
/// Their dimensionless normalization is retained; no eigenbasis is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TensorPronyHistoryTerm {
    /// An equilibrium spring contribution `strength_pa * B^T B`.
    Equilibrium {
        /// Original dimensionless Gram rows; rank deficiency is allowed.
        gram_rows: Tangent6,
    },
    /// One original Maxwell arm `strength_pa * B^T B` at a fixed time.
    Maxwell {
        /// Original dimensionless Gram rows; rank deficiency is allowed.
        gram_rows: Tangent6,
        /// Positive fixed relaxation time [s]. Equal times with different
        /// tensor directions are allowed, and original branch order is retained.
        relaxation_time_s: f64,
    },
    /// A Kelvin dashpot contribution `strength_pa_s * B^T B`.
    Kelvin {
        /// Original dimensionless Gram rows; rank deficiency is allowed.
        gram_rows: Tangent6,
    },
}

impl TensorPronyHistoryTerm {
    fn factor(&self) -> &Tangent6 {
        match self {
            Self::Equilibrium { gram_rows }
            | Self::Maxwell { gram_rows, .. }
            | Self::Kelvin { gram_rows } => gram_rows,
        }
    }
}

/// One six-component endpoint observation in an original material experiment.
#[derive(Debug, Clone, Copy)]
pub struct TensorPronyHistorySample {
    /// Nonnegative time [s], strictly increasing within this experiment.
    pub time_s: f64,
    /// Dimensionless SMALL strain in tensor-shear Voigt coordinates, not
    /// engineering shear. Consecutive observations define linear ramps.
    pub strain: Voigt,
    /// Observed physical stress [Pa], in the same tensor-shear Voigt order.
    pub stress_pa: Voigt,
    /// Positive absolute error allowed for EACH observed stress component [Pa].
    /// Reciprocals weight independent scalar residuals; no noise floor is added.
    pub stress_error_pa: Voigt,
}

/// Explicit capacities for the joint small dense fitting problem.
#[derive(Debug, Clone, Copy)]
pub struct TensorPronyHistoryFitBudget {
    /// Maximum number of independent, separately relaxed experiments.
    pub max_experiments: usize,
    /// Maximum TOTAL endpoint observations across all experiments.
    pub max_samples: usize,
    /// Maximum fitted tensor strengths; the existing NNLS owner supports 64.
    pub max_terms: usize,
}

/// Admitted passive law and independent replay of every calibration experiment.
#[derive(Debug, Clone)]
pub struct TensorPronyHistoryFit {
    /// Canonical material with SPD equilibrium and nonnegative damping channels.
    pub model: TensorProny,
    /// Original term order, including zeros. Units are Pa for equilibrium and
    /// Maxwell terms, Pa s for Kelvin terms, at the authored factor normalization.
    pub strengths: Vec<f64>,
    /// Original history replay, indexed by experiment then observation [Pa].
    pub predicted_stress_pa: Vec<Vec<Voigt>>,
    /// Largest component error divided by its own supplied bound (<= 1).
    pub max_error_ratio: f64,
    /// RMS of those scalar ratios; this is not an uncertainty estimate.
    pub rms_error_ratio: f64,
}

/// No partial fitted law is returned on refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum TensorPronyHistoryFitError {
    /// Invalid experiment, tensor basis, time, strain limit or arithmetic.
    Invalid(&'static str),
    /// An explicit input capacity was exceeded before matrix allocation.
    Capacity {
        /// Capacity name.
        quantity: &'static str,
        /// Requested count.
        requested: usize,
        /// Effective limit, including the fixed 64-term dense ceiling.
        limit: usize,
    },
    /// This original term has exactly zero stress on all supplied histories.
    /// Its strength cannot be fitted from these experiments.
    UnexcitedTerm(usize),
    /// Original equilibrium or Maxwell admission/history refused the fitted law.
    Material(TensorPronyError),
    /// Original Kelvin tensor admission refused the fitted viscosity.
    Viscosity(ViscoError),
    /// Independent replay missed an authored component error bound.
    StressResidual {
        /// Original experiment index.
        experiment: usize,
        /// Original observation index within that experiment.
        sample: usize,
        /// Tensor-shear stress component in 0..6.
        component: usize,
        /// Absolute error divided by this component's allowed error.
        ratio: f64,
    },
    /// A numerical payload could not be reserved.
    Allocation,
    /// Cooperative cancellation before publication of a complete fitted law.
    Cancelled,
}

impl core::fmt::Display for TensorPronyHistoryFitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MAT-TENSOR-HISTORY-FIT: {self:?}")
    }
}
impl std::error::Error for TensorPronyHistoryFitError {}
type Error = TensorPronyHistoryFitError;
type Result<T> = core::result::Result<T, Error>;

/// Jointly fit a passive tensor law to original multiaxial loading experiments.
///
/// The objective sums squared stress errors divided by their authored squared
/// component bounds. Each experiment starts independently with relaxed memory
/// and zero strain rate. At least one equilibrium term and one observation per
/// experiment are required; six times the total sample count must be at least
/// the number of fitted strengths. The equilibrium SUM must be strictly SPD
/// after fitting; no stiffness floor or damping repair is introduced.
///
/// # Errors
/// Invalid data/budgets, unexcited terms, failed passive admission, an unmet
/// observation bound, allocation failure or unrepresentable arithmetic.
pub fn fit_tensor_prony_histories(
    experiments: &[&[TensorPronyHistorySample]],
    terms: &[TensorPronyHistoryTerm],
    strain_limit: f64,
    budget: TensorPronyHistoryFitBudget,
) -> Result<TensorPronyHistoryFit> {
    fit_tensor_prony_histories_with_cancel(experiments, terms, strain_limit, budget, || false)
}

/// Cancellable joint fitting with the same original-history semantics.
///
/// Polls admission, matrix assembly, column equilibration, material construction
/// and independent replay, and immediately before/after the existing bounded
/// small dense NNLS solve. It does not interrupt inside that solver.
///
/// # Errors
/// [`fit_tensor_prony_histories`]'s refusals, or cancellation.
#[allow(clippy::too_many_lines)] // One bounded calibration transaction; the original solver and constitutive owners remain separate.
pub fn fit_tensor_prony_histories_with_cancel(
    experiments: &[&[TensorPronyHistorySample]],
    terms: &[TensorPronyHistoryTerm],
    strain_limit: f64,
    budget: TensorPronyHistoryFitBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<TensorPronyHistoryFit> {
    poll(&mut cancelled)?;
    capacity("experiments", experiments.len(), budget.max_experiments)?;
    capacity("tensor terms", terms.len(), budget.max_terms.min(64))?;
    let samples = experiments.iter().try_fold(0_usize, |count, e| {
        count.checked_add(e.len()).ok_or(Error::Allocation)
    })?;
    capacity("total observations", samples, budget.max_samples)?;
    let rows = samples.checked_mul(6).ok_or(Error::Allocation)?;
    let cols = terms.len();
    if experiments.is_empty()
        || experiments.iter().any(|e| e.is_empty())
        || cols == 0
        || rows < cols
        || !strain_limit.is_finite()
        || strain_limit <= 0.0
        || !terms
            .iter()
            .any(|t| matches!(t, TensorPronyHistoryTerm::Equilibrium { .. }))
    {
        return Err(Error::Invalid(
            "need nonempty experiments, equilibrium terms, enough observations and a positive strain limit",
        ));
    }
    let mut tensors = reserved(cols)?;
    for term in terms {
        poll(&mut cancelled)?;
        if let TensorPronyHistoryTerm::Maxwell {
            relaxation_time_s: tau,
            ..
        } = term
            && (!tau.is_finite() || *tau <= 0.0)
        {
            return Err(Error::Invalid(
                "Maxwell relaxation times must be finite and positive",
            ));
        }
        tensors.push(gram(term.factor()).map_err(Error::Material)?);
    }
    for experiment in experiments {
        for (i, sample) in experiment.iter().enumerate() {
            poll(&mut cancelled)?;
            let strain = mandel(&sample.strain);
            let norm = det::sqrt(dot(&strain, &strain));
            if !sample.time_s.is_finite()
                || sample.time_s < 0.0
                || !norm.is_finite()
                || norm > strain_limit
                || sample.strain.iter().any(|v| !v.is_finite())
                || sample.stress_pa.iter().any(|v| !v.is_finite())
                || sample
                    .stress_error_pa
                    .iter()
                    .any(|v| !v.is_finite() || *v <= 0.0)
                || (i > 0 && sample.time_s <= experiment[i - 1].time_s)
            {
                return Err(Error::Invalid(
                    "observations need increasing finite times, admitted strains, finite stresses and positive errors",
                ));
            }
        }
    }
    let extent = rows.checked_mul(cols).ok_or(Error::Allocation)?;
    let mut matrix = filled(extent, 0.0)?;
    let mut rhs = filled(rows, 0.0)?;
    let mut memory = filled(cols, [0.0; 6])?;
    let mut row = 0;
    for experiment in experiments {
        memory.fill([0.0; 6]);
        for (i, sample) in experiment.iter().enumerate() {
            poll(&mut cancelled)?;
            let strain = mandel(&sample.strain);
            let (dt, delta) = if i == 0 {
                (0.0, [0.0; 6])
            } else {
                let previous = mandel(&experiment[i - 1].strain);
                (
                    sample.time_s - experiment[i - 1].time_s,
                    core::array::from_fn(|k| strain[k] - previous[k]),
                )
            };
            for (j, term) in terms.iter().enumerate() {
                let argument = match term {
                    TensorPronyHistoryTerm::Equilibrium { .. } => strain,
                    TensorPronyHistoryTerm::Maxwell {
                        relaxation_time_s, ..
                    } => {
                        if i > 0 {
                            let ramp = PronyRamp::new(dt / relaxation_time_s);
                            for (q, d) in memory[j].iter_mut().zip(delta) {
                                *q = ramp.apply(*q, d).0;
                            }
                        }
                        memory[j]
                    }
                    TensorPronyHistoryTerm::Kelvin { .. } => {
                        if i == 0 {
                            [0.0; 6]
                        } else {
                            delta.map(|d| d / dt)
                        }
                    }
                };
                let projected = apply(term.factor(), &argument);
                let mut stress = [0.0; 6];
                transpose_add(term.factor(), &projected, &mut stress);
                for k in 0..6 {
                    matrix[(row + k) * cols + j] =
                        stress[k] / WEIGHTS[k] / sample.stress_error_pa[k];
                }
            }
            for k in 0..6 {
                rhs[row + k] = sample.stress_pa[k] / sample.stress_error_pa[k];
            }
            if rhs[row..row + 6].iter().any(|v| !v.is_finite())
                || matrix[row * cols..(row + 6) * cols]
                    .iter()
                    .any(|v| !v.is_finite())
            {
                return Err(Error::Invalid(
                    "weighted tensor experiment is unrepresentable",
                ));
            }
            row += 6;
        }
    }
    let mut scales = filled(cols, (0.0, 0.0))?;
    for j in 0..cols {
        poll(&mut cancelled)?;
        let scale = (0..rows)
            .map(|i| matrix[i * cols + j].abs())
            .fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(Error::UnexcitedTerm(j));
        }
        let norm = det::sqrt(
            (0..rows)
                .map(|i| {
                    let value = matrix[i * cols + j] / scale;
                    value * value
                })
                .sum(),
        );
        // Factored scaling avoids overflow of scale*norm and retains the
        // authored tensor normalization when physical strengths are recovered.
        scales[j] = (scale, norm);
        for i in 0..rows {
            matrix[i * cols + j] = (matrix[i * cols + j] / scale) / norm;
        }
    }
    poll(&mut cancelled)?;
    let mut strengths = nnls(&matrix, &rhs, rows, cols);
    poll(&mut cancelled)?;
    for (strength, &(scale, norm)) in strengths.iter_mut().zip(&scales) {
        *strength = (*strength / norm) / scale;
        if !strength.is_finite() || *strength < 0.0 {
            return Err(Error::Invalid(
                "nonnegative tensor strength is unrepresentable",
            ));
        }
    }
    let mut equilibrium = [[0.0; 6]; 6];
    let mut viscosity = [[0.0; 6]; 6];
    let mut has_viscosity = false;
    let mut branches = reserved(cols)?;
    for ((term, original), &strength) in terms.iter().zip(&tensors).zip(&strengths) {
        poll(&mut cancelled)?;
        match term {
            TensorPronyHistoryTerm::Maxwell {
                gram_rows,
                relaxation_time_s,
            } => {
                let root = det::sqrt(strength);
                let factor = gram_rows.map(|r| r.map(|x| root * x));
                branches.push(
                    TensorPronyBranch::from_gram_rows(factor, *relaxation_time_s)
                        .map_err(Error::Material)?,
                );
            }
            TensorPronyHistoryTerm::Equilibrium { .. } | TensorPronyHistoryTerm::Kelvin { .. } => {
                let matrix = if matches!(term, TensorPronyHistoryTerm::Equilibrium { .. }) {
                    &mut equilibrium
                } else {
                    has_viscosity = true;
                    &mut viscosity
                };
                for i in 0..6 {
                    for j in 0..6 {
                        matrix[i][j] += strength * original[i][j];
                    }
                }
            }
        }
    }
    let mut model =
        TensorProny::new_mandel(equilibrium, branches, strain_limit).map_err(Error::Material)?;
    if has_viscosity {
        model = model
            .with_kelvin_viscosity_mandel(viscosity)
            .map_err(Error::Viscosity)?;
    }
    let shared = Arc::new(model.clone());
    let mut predicted_stress_pa = reserved(experiments.len())?;
    let mut max_error_ratio = 0.0_f64;
    let mut squared_error = 0.0;
    for (experiment_index, experiment) in experiments.iter().enumerate() {
        poll(&mut cancelled)?;
        let mut history = relaxed_history(&shared, experiment[0].strain)?;
        let mut predictions = reserved(experiment.len())?;
        for (i, sample) in experiment.iter().enumerate() {
            poll(&mut cancelled)?;
            let stress = if i == 0 {
                history.stress()
            } else {
                history
                    .advance(
                        sample.strain,
                        sample.time_s - experiment[i - 1].time_s,
                        &mut cancelled,
                    )
                    .map_err(material_error)?
                    .stress
            };
            for (k, &value) in stress.iter().enumerate() {
                let ratio = ((value - sample.stress_pa[k]) / sample.stress_error_pa[k]).abs();
                if !ratio.is_finite() || ratio > 1.0 {
                    return Err(Error::StressResidual {
                        experiment: experiment_index,
                        sample: i,
                        component: k,
                        ratio,
                    });
                }
                max_error_ratio = max_error_ratio.max(ratio);
                squared_error = ratio.mul_add(ratio, squared_error);
            }
            predictions.push(stress);
        }
        predicted_stress_pa.push(predictions);
    }
    poll(&mut cancelled)?;
    Ok(TensorPronyHistoryFit {
        model,
        strengths,
        predicted_stress_pa,
        max_error_ratio,
        rms_error_ratio: det::sqrt(squared_error / rows as f64),
    })
}

fn relaxed_history(model: &Arc<TensorProny>, strain: Voigt) -> Result<TensorPronyHistory> {
    let strain = mandel(&strain);
    let elastic = apply(&model.equilibrium_gram, &strain);
    let mut stress = [0.0; 6];
    transpose_add(&model.equilibrium_gram, &elastic, &mut stress);
    let stored = 0.5 * dot(&elastic, &elastic);
    if !stored.is_finite() || stress.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid(
            "relaxed initial stress or storage is unrepresentable",
        ));
    }
    Ok(TensorPronyHistory {
        model: Arc::clone(model),
        strain,
        stress,
        q: filled(model.branches().len(), [0.0; 6])?,
        stored,
        work: 0.0,
        heat: 0.0,
    })
}

fn material_error(error: TensorPronyError) -> Error {
    if error == TensorPronyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Material(error)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn capacity(quantity: &'static str, requested: usize, limit: usize) -> Result<()> {
    if requested > limit {
        Err(Error::Capacity {
            quantity,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}
fn reserved<T>(len: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| Error::Allocation)?;
    Ok(values)
}
fn filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>> {
    let mut values = reserved(len)?;
    values.resize(len, value);
    Ok(values)
}

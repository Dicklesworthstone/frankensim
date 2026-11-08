//! Prescribed strain/stress histories to a passive Maxwell material.
//!
//! The caller supplies the relaxation-time ladder and the absolute stress
//! error allowed at each observation. Original piecewise-linear strain ramps
//! generate every design-matrix column through [`PronyRamp`]; the existing
//! nonnegative least-squares owner fits equilibrium and branch moduli. No
//! sampled modal damping, finite-difference dynamics or new optimizer enters.
//!
//! The FIRST observation declares fully relaxed branch memory at its strain.
//! This is a physical initial condition, not an inferred prehistory. Fitting
//! does not identify relaxation times, prove parameter uniqueness, validate
//! between observations, or establish a microscopic dissipation mechanism.

use super::{GeneralizedMaxwell, PronyRamp, ViscoError, nnls};
use fs_math::det;

/// An endpoint observation of the original strain-controlled experiment.
#[derive(Debug, Clone, Copy)]
pub struct PronyHistorySample {
    /// Nonnegative time [s], strictly increasing through the record.
    pub time_s: f64,
    /// Total dimensionless strain. Consecutive values define linear ramps.
    pub strain: f64,
    /// Observed signed stress [Pa], including compression and reversal.
    pub stress_pa: f64,
    /// Positive absolute stress error allowed at THIS observation [Pa].
    /// Its reciprocal weights the residual, and the final replay must meet
    /// every individual bound. No noise floor or uncertainty is invented.
    pub stress_error_pa: f64,
}

/// Explicit data and parameter capacities for the dense fitting problem.
#[derive(Debug, Clone, Copy)]
pub struct PronyHistoryFitBudget {
    /// Maximum number of observations admitted before allocation.
    pub max_samples: usize,
    /// Maximum number of supplied relaxation times. At most 63 are supported
    /// by this small dense path, leaving one of 64 unknowns for equilibrium.
    pub max_branches: usize,
}

/// A passive fitted model and replay diagnostics on the supplied observations.
#[derive(Debug, Clone)]
pub struct PronyHistoryFit {
    /// Original caller-ordered relaxation times with nonnegative fitted moduli.
    /// Zero coefficients are retained, so parameter indices do not change.
    pub model: GeneralizedMaxwell,
    /// Stress from replaying the fitted material's ORIGINAL `step` method [Pa].
    pub predicted_stress_pa: Vec<f64>,
    /// Largest absolute stress error divided by its own authored bound (<= 1).
    pub max_error_ratio: f64,
    /// RMS of the same dimensionless residuals; not an uncertainty estimate.
    pub rms_error_ratio: f64,
}

/// No partial fitted material is returned on refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PronyHistoryFitError {
    /// Invalid experiment, parameter ladder or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Explicit sample or parameter capacity was exceeded.
    Capacity {
        /// Capacity name.
        quantity: &'static str,
        /// Requested count.
        requested: usize,
        /// Admitted count, including the dense owner's fixed parameter cap.
        limit: usize,
    },
    /// A design column is exactly zero on this experiment: equilibrium is
    /// column 0, and branch j is column j+1. Its modulus cannot be fitted here.
    UnexcitedColumn(usize),
    /// The passive solve did not resolve a positive equilibrium modulus.
    NonPositiveEquilibrium,
    /// A fitted replay exceeds an observation's authored absolute error.
    StressResidual {
        /// First failing observation index.
        sample: usize,
        /// Observed error divided by the supplied absolute bound.
        ratio: f64,
    },
    /// The fitted parameters failed the material owner's admission.
    Material(ViscoError),
    /// Numerical payload could not be reserved.
    Allocation,
    /// Cancellation observed before a complete fit was published.
    Cancelled,
}

impl core::fmt::Display for PronyHistoryFitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MAT-PRONY-HISTORY-FIT: {self:?}")
    }
}
impl std::error::Error for PronyHistoryFitError {}

/// Fit positive equilibrium and nonnegative branch moduli on authored times.
///
/// The first row has fully relaxed branch memory, including for nonzero
/// initial strain. Every later interval is the exact constitutive response
/// to its original linear strain ramp. At least two observations and at least
/// as many observations as fitted coefficients are required. An empty ladder
/// fits a purely elastic material; duplicate relaxation times refuse.
///
/// The objective is the sum of squared stress errors divided by the squared
/// observation-specific error bounds. Column equilibration prevents a change
/// of stress units from changing the nonnegative solver's active-set threshold.
/// Acceptance independently replays the public material state and requires
/// EVERY stress observation to meet its bound. This is a sampled calibration
/// result, not independent validation or an identifiability certificate.
///
/// # Errors
/// Invalid or unexcited data, insufficient capacities, failed passive fitting,
/// an unmet observation error, or unrepresentable material/replay arithmetic.
pub fn fit_prony_history(
    samples: &[PronyHistorySample],
    relaxation_times_s: &[f64],
    budget: PronyHistoryFitBudget,
) -> Result<PronyHistoryFit, PronyHistoryFitError> {
    fit_prony_history_with_cancel(samples, relaxation_times_s, budget, || false)
}

/// Cancellable fitting with the same physical and numerical semantics.
///
/// Polls occur during admission, column assembly, equilibration and replay,
/// and immediately around the existing small dense NNLS solve. This does not
/// promise interruption inside that owner's bounded active-set solve.
///
/// # Errors
/// [`fit_prony_history`]'s refusals, plus cancellation before publication.
#[allow(clippy::too_many_lines)] // One admitted experiment, constrained fit and independent constitutive replay.
pub fn fit_prony_history_with_cancel(
    samples: &[PronyHistorySample],
    relaxation_times_s: &[f64],
    budget: PronyHistoryFitBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PronyHistoryFit, PronyHistoryFitError> {
    use PronyHistoryFitError as Error;
    let poll = |cancelled: &mut dyn FnMut() -> bool| {
        if cancelled() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    };
    poll(&mut cancelled)?;
    let rows = samples.len();
    let branches = relaxation_times_s.len();
    for (quantity, requested, limit) in [
        ("observations", rows, budget.max_samples),
        ("Maxwell branches", branches, budget.max_branches.min(63)),
    ] {
        if requested > limit {
            return Err(Error::Capacity {
                quantity,
                requested,
                limit,
            });
        }
    }
    let cols = branches + 1;
    if rows < 2 || rows < cols {
        return Err(Error::Invalid(
            "need two observations and at least one per fitted modulus",
        ));
    }
    for (i, &tau) in relaxation_times_s.iter().enumerate() {
        if !tau.is_finite() || tau <= 0.0 || relaxation_times_s[..i].contains(&tau) {
            return Err(Error::Invalid(
                "relaxation times must be finite, positive and distinct",
            ));
        }
    }
    for (i, sample) in samples.iter().enumerate() {
        poll(&mut cancelled)?;
        if ![
            sample.time_s,
            sample.strain,
            sample.stress_pa,
            sample.stress_error_pa,
        ]
        .into_iter()
        .all(f64::is_finite)
            || sample.time_s < 0.0
            || sample.stress_error_pa <= 0.0
        {
            return Err(Error::Invalid(
                "observations need finite time/strain/stress and positive stress error",
            ));
        }
        if i > 0
            && (sample.time_s <= samples[i - 1].time_s
                || !(sample.strain - samples[i - 1].strain).is_finite())
        {
            return Err(Error::Invalid(
                "times must increase and strain increments must be representable",
            ));
        }
    }
    let extent = rows.checked_mul(cols).ok_or(Error::Allocation)?;
    let mut matrix = zeros(extent)?;
    let mut rhs = zeros(rows)?;
    let mut memory = vec![0.0; branches];
    for (i, sample) in samples.iter().enumerate() {
        poll(&mut cancelled)?;
        if i > 0 {
            let dt = sample.time_s - samples[i - 1].time_s;
            let delta = sample.strain - samples[i - 1].strain;
            for (q, &tau) in memory.iter_mut().zip(relaxation_times_s) {
                *q = PronyRamp::new(dt / tau).apply(*q, delta).0;
            }
        }
        matrix[i * cols] = sample.strain / sample.stress_error_pa;
        for (j, &q) in memory.iter().enumerate() {
            matrix[i * cols + j + 1] = q / sample.stress_error_pa;
        }
        rhs[i] = sample.stress_pa / sample.stress_error_pa;
        if !rhs[i].is_finite()
            || matrix[i * cols..(i + 1) * cols]
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(Error::Invalid("weighted experiment cannot be represented"));
        }
    }
    let mut scales = vec![(0.0, 0.0); cols];
    for j in 0..cols {
        poll(&mut cancelled)?;
        let scale = (0..rows)
            .map(|i| matrix[i * cols + j].abs())
            .fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(Error::UnexcitedColumn(j));
        }
        let norm = det::sqrt(
            (0..rows)
                .map(|i| {
                    let value = matrix[i * cols + j] / scale;
                    value * value
                })
                .sum(),
        );
        // Keep the norm factored: scale*norm can overflow although the
        // normalized column and final physical modulus remain representable.
        scales[j] = (scale, norm);
        for i in 0..rows {
            matrix[i * cols + j] = (matrix[i * cols + j] / scale) / norm;
        }
    }
    poll(&mut cancelled)?;
    let mut coefficients = nnls(&matrix, &rhs, rows, cols);
    poll(&mut cancelled)?;
    for (coefficient, &(scale, norm)) in coefficients.iter_mut().zip(&scales) {
        *coefficient = (*coefficient / norm) / scale;
        if !coefficient.is_finite() || *coefficient < 0.0 {
            return Err(Error::Invalid(
                "nonnegative fit produced an unrepresentable modulus",
            ));
        }
    }
    if coefficients[0] <= 0.0 {
        return Err(Error::NonPositiveEquilibrium);
    }
    let model = GeneralizedMaxwell::new(
        coefficients[0],
        coefficients[1..]
            .iter()
            .copied()
            .zip(relaxation_times_s.iter().copied())
            .collect(),
    )
    .map_err(Error::Material)?;
    let mut state = model.state();
    state.strain = samples[0].strain;
    state.stress = model.e_inf * state.strain;
    let mut predicted_stress_pa = zeros(rows)?;
    let mut max_error_ratio = 0.0_f64;
    let mut squared_error = 0.0;
    for (i, sample) in samples.iter().enumerate() {
        poll(&mut cancelled)?;
        let stress = if i == 0 {
            state.stress
        } else {
            model.step(
                &mut state,
                sample.strain,
                sample.time_s - samples[i - 1].time_s,
            )
        };
        if ![stress, state.work, state.dissipated]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid(
                "fitted constitutive stress or energy overflowed on replay",
            ));
        }
        let ratio = ((stress - sample.stress_pa) / sample.stress_error_pa).abs();
        if ratio > 1.0 {
            return Err(Error::StressResidual { sample: i, ratio });
        }
        predicted_stress_pa[i] = stress;
        max_error_ratio = max_error_ratio.max(ratio);
        squared_error = ratio.mul_add(ratio, squared_error);
    }
    poll(&mut cancelled)?;
    Ok(PronyHistoryFit {
        model,
        predicted_stress_pa,
        max_error_ratio,
        rms_error_ratio: det::sqrt(squared_error / rows as f64),
    })
}

fn zeros(len: usize) -> Result<Vec<f64>, PronyHistoryFitError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| PronyHistoryFitError::Allocation)?;
    values.resize(len, 0.0);
    Ok(values)
}

//! Nonlinear Gaussian-prior estimation with one retained measurement batch.
//!
//! Write `x = m + L z`, where the original prior is `P = L L^T`, and minimize
//! `0.5 |z|² + 0.5 sum_i ((y_i - h_i(x))/sigma_i)²`. Each Gauss–Newton
//! proposal assimilates the current affine observations into the SAME standard
//! normal prior. Iteration never assimilates the measurements a second time.
//! This is the iterated Kalman/Gauss–Newton relationship of Bell and Cathey,
//! IEEE TAC 38(2), 1993, DOI: 10.1109/9.250476, with an Armijo line search.
//!
//! The existing UD measurement and square-root affine prediction owners solve
//! the linearized problem and transport its covariance. Singular prior factors
//! retain their null directions without jitter. Only a stationary result
//! publishes a Gaussian approximation, whose covariance is the Gauss–Newton
//! curvature approximation at the accepted state, not an exact nonlinear
//! posterior, a certified minimum, or a validated physical-model claim.

use fs_blake3::DomainHasher;
use fs_exec::Cx;

use super::{hash_belief, hash_bytes, hash_f64, hash_u64};
use crate::factor::{
    FactorBelief, batch::assimilate_sequence, forecast::LinearDynamics, scalar_factor_work_estimate,
};
use crate::{
    AssimError, Belief, Observation, WorkPlan, WorkProgress, canonicalize_zero,
    validate_leaf_identity, validate_noise,
};

const MAX_ROWS: usize = 256;
const MAX_ITERATIONS: usize = 256;
const MAX_EVALUATIONS: usize = 65_536;
const ARMIJO: f64 = 1.0e-4;

/// A fixed set of scalar readings with explicitly independent Gaussian noise.
/// Rows and model outputs use declaration order; independence is a caller
/// assumption. Shared calibration must be represented in the state or handled
/// by the correlated-noise APIs, never silently discarded here.
#[derive(Debug, Clone, PartialEq)]
pub struct IndependentMeasurements {
    readings: Vec<f64>,
    variances: Vec<f64>,
    instruments: Vec<String>,
    model_id: String,
}

impl IndependentMeasurements {
    /// Admit finite readings, strictly positive variances, and unique instrument
    /// identities. `model_id` identifies the caller's pure prediction/Jacobian
    /// implementation; it does not authenticate arbitrary callback code.
    pub fn new(
        readings: Vec<f64>,
        variances: Vec<f64>,
        instruments: Vec<String>,
        model_id: impl Into<String>,
    ) -> Result<Self, AssimError> {
        if readings.is_empty() {
            return Err(AssimError::EmptyObservations);
        }
        if readings.len() > MAX_ROWS {
            return Err(AssimError::ObservationCountLimit {
                count: readings.len(),
                max: MAX_ROWS,
            });
        }
        if readings.len() != variances.len() || readings.len() != instruments.len() {
            return Err(invalid("nonlinear measurement row counts differ"));
        }
        let model_id = model_id.into();
        validate_leaf_identity("nonlinear model", &model_id)?;
        for (i, ((&value, &variance), instrument)) in readings
            .iter()
            .zip(&variances)
            .zip(&instruments)
            .enumerate()
        {
            if !value.is_finite() {
                return Err(AssimError::NonFiniteObservationValue);
            }
            validate_noise(variance)?;
            validate_leaf_identity("instrument", instrument)?;
            if instruments[..i].contains(instrument) {
                return Err(AssimError::DuplicateObservationInstrument {
                    instrument: instrument.clone(),
                });
            }
        }
        Ok(Self {
            readings,
            variances,
            instruments,
            model_id,
        })
    }

    /// Number of observed scalars, also the required model/Jacobian row count.
    #[must_use]
    pub fn len(&self) -> usize {
        self.readings.len()
    }
    /// Admitted batches are never empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }
    /// Caller-declared model identity retained in the result identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }
}

/// One model evaluation. The caller may use analytic derivatives, `fs-ad`, or
/// a production adjoint. The implementation checks dimensions and finiteness,
/// not derivative correctness. Use independent derivative checks for the model.
#[derive(Debug, Clone, PartialEq)]
pub enum ModelEvaluation {
    /// Complete model values and physical-coordinate Jacobian rows.
    Value {
        /// One prediction in each declared reading's units.
        predictions: Vec<f64>,
        /// `jacobian[i][j] = d predictions[i] / d state[j]`.
        jacobian: Vec<Vec<f64>>,
    },
    /// A trial left the model's declared physical domain. The line search may
    /// shorten its step. An out-of-domain prior refuses before iteration.
    OutsideDomain,
}

/// Explicit numerical allowances and tolerances in prior-whitened coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IteratedSettings {
    /// Maximum accepted Gauss–Newton steps, in `0..=256`.
    pub max_iterations: usize,
    /// Complete prediction/Jacobian evaluations, including rejected trials.
    pub max_model_evaluations: usize,
    /// Halvings after the full-step trial, in `0..=64`.
    pub max_backtracks: usize,
    /// Absolute infinity-norm stationarity threshold for the whitened gradient.
    pub gradient_tolerance: f64,
    /// Relative infinity-norm step floor. A small step alone is not convergence.
    pub step_tolerance: f64,
}

impl IteratedSettings {
    fn validate(self) -> Result<(), AssimError> {
        if self.max_iterations > MAX_ITERATIONS
            || self.max_model_evaluations == 0
            || self.max_model_evaluations > MAX_EVALUATIONS
            || self.max_backtracks > 64
            || !self.gradient_tolerance.is_finite()
            || self.gradient_tolerance <= 0.0
            || !self.step_tolerance.is_finite()
            || self.step_tolerance <= 0.0
        {
            return Err(invalid("iterated nonlinear settings"));
        }
        Ok(())
    }
}

/// The measured stopping condition; only `Stationary` carries a posterior.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IteratedOutcome {
    /// The accepted state's whitened gradient met the declared threshold.
    Stationary,
    /// The accepted-update allowance ended before stationarity.
    IterationLimit,
    /// The callback allowance ended, counting rejected/domain-invalid trials.
    EvaluationLimit,
    /// No tested step met sufficient decrease within the backtracking allowance.
    LineSearchStalled,
    /// The proposed/halved step was too small, without proving stationarity.
    StepStalled,
}

/// One accepted state; the first row is the original prior at step fraction zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IteratedRecord {
    /// Gaussian-prior plus measurement negative log likelihood, up to constants.
    pub objective: f64,
    /// Infinity norm of the gradient in prior-whitened coordinates.
    pub gradient_norm: f64,
    /// Accepted fraction of the Gauss–Newton proposal (zero for the prior).
    pub step_fraction: f64,
}

/// Atomic nonlinear result. Limits retain the best accepted state and trace,
/// but no posterior covariance is presented as a converged update.
#[derive(Debug, Clone, PartialEq)]
pub struct IteratedAssimilation {
    outcome: IteratedOutcome,
    state: Vec<f64>,
    posterior: Option<Belief>,
    trace: Vec<IteratedRecord>,
    evaluations: usize,
    rejected: usize,
    receipt_id: String,
}

impl IteratedAssimilation {
    /// Actual stopping condition.
    #[must_use]
    pub const fn outcome(&self) -> IteratedOutcome {
        self.outcome
    }
    /// Last accepted physical state, including for an unresolved numerical stop.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    /// Local Gaussian approximation only when the stationarity check passed.
    #[must_use]
    pub fn posterior(&self) -> Option<&Belief> {
        self.posterior.as_ref()
    }
    /// Objective/stationarity evidence at every accepted state.
    #[must_use]
    pub fn trace(&self) -> &[IteratedRecord] {
        &self.trace
    }
    /// All callbacks consumed, including unsuccessful line-search trials.
    #[must_use]
    pub const fn model_evaluations(&self) -> usize {
        self.evaluations
    }
    /// Domain-invalid and insufficient-decrease trials.
    #[must_use]
    pub const fn rejected_trials(&self) -> usize {
        self.rejected
    }
    /// Content identity binding inputs, settings, callback trace and output.
    #[must_use]
    pub fn receipt_id(&self) -> &str {
        &self.receipt_id
    }
}

struct Point {
    state: Vec<f64>,
    residual: Vec<f64>,
    jacobian: Vec<Vec<f64>>,
    gradient: Vec<f64>,
    objective: f64,
}

/// Assimilate one independent nonlinear measurement batch by damped iteration.
///
/// The model returns predictions and a Jacobian together and receives the same
/// cancellation context. It owns intra-evaluation polling and its numerical
/// work/memory; this driver caps the number of calls and checks cancellation
/// before/after each. Rejected calls consume that allowance. Driver work uses
/// the existing cumulative `Cx` accountant and bounded dense factor kernels.
/// A callback error, cancellation or invalid arithmetic publishes no result.
/// Numerical limits return an explicit unresolved outcome, never a posterior.
pub fn assimilate_iterated<F>(
    prior: &Belief,
    measurements: &IndependentMeasurements,
    settings: IteratedSettings,
    mut model: F,
    cx: &Cx<'_>,
) -> Result<IteratedAssimilation, AssimError>
where
    F: FnMut(&[f64], &Cx<'_>) -> Result<ModelEvaluation, AssimError>,
{
    settings.validate()?;
    let n = prior.dim();
    let m = measurements.len();
    let (nu, mu, ev) = (n as u128, m as u128, settings.max_model_evaluations as u128);
    // Charge each reused factor kernel's conservative envelope to this ONE
    // accountant before calling it. Its own local admission must not reissue
    // the workflow's remaining cost/poll allowance at every nonlinear iterate.
    // The caller's model still owns its internal work; calls are separately capped.
    let construction = scalar_factor_work_estimate(n)?.construction + 64;
    let publication = 64 * nu * nu * nu + 128 * nu * nu + 128 * nu + 2048;
    let kernels = 2 * construction
        + publication
        + (settings.max_iterations as u128 + 1) * sequence_work(n, m)?;
    let work = ev * (8 * mu * nu * nu + 64 * mu * nu + 64 * nu + 256)
        + (settings.max_iterations as u128 + 2) * (32 * mu * nu + 64 * nu + 256)
        + kernels;
    let bytes = ev * (32 * mu * nu + 32 * mu + 32 * nu + 64)
        + 32 * nu * nu
        + 4096
        + measurements
            .instruments
            .iter()
            .map(|s| s.len() as u128 + 8)
            .sum::<u128>()
        + measurements.model_id.len() as u128;
    let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, work, 0, 0, bytes)?)?;
    progress.checkpoint("iterated nonlinear admission")?;
    progress.scalar("nonlinear prior factor allowance", construction)?;
    let factors = FactorBelief::from_belief(prior, cx)?;
    let mut root = vec![vec![0.0; n]; n];
    for (i, row) in root.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate().skip(i) {
            let u = if i == j {
                1.0
            } else {
                factors.upper_entry(i, j).expect("valid factor index")
            };
            *value = finite(
                u * factors.diag(j).expect("valid pivot").sqrt(),
                "nonlinear prior root",
            )?;
            progress.scalar("nonlinear prior root", 2)?;
        }
    }
    progress.scalar("nonlinear standard prior allowance", construction)?;
    let standard = FactorBelief::diagonal(vec![0.0; n], &vec![1.0; n], cx)?;
    let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.iterated-nonlinear.v1");
    hash_belief(&mut hash, prior, &mut progress)?;
    hash_bytes(&mut hash, measurements.model_id.as_bytes(), &mut progress)?;
    hash_u64(&mut hash, m as u64, &mut progress)?;
    for ((&reading, &variance), instrument) in measurements
        .readings
        .iter()
        .zip(&measurements.variances)
        .zip(&measurements.instruments)
    {
        hash_f64(&mut hash, reading, &mut progress)?;
        hash_f64(&mut hash, variance, &mut progress)?;
        hash_bytes(&mut hash, instrument.as_bytes(), &mut progress)?;
    }
    for count in [
        settings.max_iterations,
        settings.max_model_evaluations,
        settings.max_backtracks,
    ] {
        hash_u64(&mut hash, count as u64, &mut progress)?;
    }
    hash_f64(&mut hash, settings.gradient_tolerance, &mut progress)?;
    hash_f64(&mut hash, settings.step_tolerance, &mut progress)?;
    let mut z = vec![0.0; n];
    let mut evaluations = 0;
    let mut rejected = 0;
    let mut point = evaluate(
        prior,
        &root,
        &z,
        measurements,
        &mut model,
        cx,
        &mut progress,
        &mut hash,
        &mut evaluations,
    )?
    .ok_or_else(|| invalid("nonlinear model excludes the prior mean"))?;
    let mut trace = vec![record(&point, 0.0)];
    let mut posterior = None;
    let outcome = 'iterations: loop {
        progress.checkpoint("iterated nonlinear stationarity")?;
        let gradient_norm = norm(&point.gradient);
        if gradient_norm <= settings.gradient_tolerance {
            // Recompute the covariance at the accepted Jacobian with the
            // ORIGINAL standard prior, then transport with the existing root owner.
            let local = affine_posterior(&standard, &z, &point, measurements, cx, &mut progress)?;
            progress.scalar("nonlinear covariance publication allowance", publication)?;
            let zero_noise = FactorBelief::diagonal(vec![0.0; n], &vec![0.0; n], cx)?;
            let transform = LinearDynamics::new(
                root.clone(),
                prior.mean().to_vec(),
                zero_noise,
                "iterated-nonlinear-prior-root",
                cx,
            )?;
            let mut physical = transform.predict(&local, cx)?.to_belief(cx)?;
            // Curvature belongs to this accepted state, not to a tiny, untested
            // extra Gauss–Newton step implied by the affine subproblem mean.
            physical.mean.clone_from(&point.state);
            posterior = Some(physical);
            break IteratedOutcome::Stationary;
        }
        if trace.len() - 1 >= settings.max_iterations {
            break IteratedOutcome::IterationLimit;
        }
        if evaluations >= settings.max_model_evaluations {
            break IteratedOutcome::EvaluationLimit;
        }
        let proposal = affine_posterior(&standard, &z, &point, measurements, cx, &mut progress)?;
        let direction = proposal
            .mean()
            .iter()
            .zip(&z)
            .map(|(&to, &from)| finite(to - from, "nonlinear direction"))
            .collect::<Result<Vec<_>, _>>()?;
        let slope = dot(&point.gradient, &direction, &mut progress)?;
        if slope >= 0.0 {
            break IteratedOutcome::StepStalled;
        }
        let floor = settings.step_tolerance * (1.0 + norm(&z));
        let mut alpha = 1.0;
        let mut accepted = None;
        for _ in 0..=settings.max_backtracks {
            progress.checkpoint("iterated nonlinear line search")?;
            if alpha * norm(&direction) <= floor {
                break 'iterations IteratedOutcome::StepStalled;
            }
            if evaluations >= settings.max_model_evaluations {
                break 'iterations IteratedOutcome::EvaluationLimit;
            }
            let trial = z
                .iter()
                .zip(&direction)
                .map(|(&a, &d)| finite(a + alpha * d, "nonlinear trial"))
                .collect::<Result<Vec<_>, _>>()?;
            progress.scalar("iterated nonlinear trial", 2 * nu)?;
            if let Some(next) = evaluate(
                prior,
                &root,
                &trial,
                measurements,
                &mut model,
                cx,
                &mut progress,
                &mut hash,
                &mut evaluations,
            )? {
                if next.objective <= point.objective + ARMIJO * alpha * slope {
                    accepted = Some((trial, next));
                    break;
                }
            }
            rejected += 1;
            alpha *= 0.5;
        }
        let Some((next_z, next)) = accepted else {
            break IteratedOutcome::LineSearchStalled;
        };
        trace.push(record(&next, alpha));
        z = next_z;
        point = next;
    };
    hash_bytes(&mut hash, b"result", &mut progress)?;
    hash_u64(&mut hash, outcome as u64, &mut progress)?;
    hash_u64(&mut hash, evaluations as u64, &mut progress)?;
    hash_u64(&mut hash, rejected as u64, &mut progress)?;
    for value in &point.state {
        hash_f64(&mut hash, *value, &mut progress)?;
    }
    if let Some(belief) = &posterior {
        hash_belief(&mut hash, belief, &mut progress)?;
    }
    hash_u64(&mut hash, trace.len() as u64, &mut progress)?;
    for row in &trace {
        for value in [row.objective, row.gradient_norm, row.step_fraction] {
            hash_f64(&mut hash, value, &mut progress)?;
        }
    }
    progress.checkpoint("iterated nonlinear publication")?;
    Ok(IteratedAssimilation {
        outcome,
        state: point.state,
        posterior,
        trace,
        evaluations,
        rejected,
        receipt_id: format!("iterated-nonlinear:v1:{}", hash.finalize()),
    })
}

#[allow(clippy::too_many_arguments)]
fn evaluate<F>(
    prior: &Belief,
    root: &[Vec<f64>],
    z: &[f64],
    measurements: &IndependentMeasurements,
    model: &mut F,
    cx: &Cx<'_>,
    progress: &mut WorkProgress<'_, '_>,
    hash: &mut DomainHasher,
    evaluations: &mut usize,
) -> Result<Option<Point>, AssimError>
where
    F: FnMut(&[f64], &Cx<'_>) -> Result<ModelEvaluation, AssimError>,
{
    progress.checkpoint("nonlinear model evaluation")?;
    let n = prior.dim();
    let mut state = Vec::with_capacity(n);
    for (row, &mean) in root.iter().zip(prior.mean()) {
        state.push(finite(
            mean + dot(row, z, progress)?,
            "nonlinear physical state",
        )?);
    }
    *evaluations += 1;
    let value = model(&state, cx)?;
    progress.checkpoint("nonlinear model completion")?;
    hash_bytes(hash, b"evaluation", progress)?;
    for &x in &state {
        hash_f64(hash, x, progress)?;
    }
    let ModelEvaluation::Value {
        predictions,
        jacobian,
    } = value
    else {
        hash_u64(hash, 0, progress)?;
        return Ok(None);
    };
    hash_u64(hash, 1, progress)?;
    if predictions.len() != measurements.len()
        || jacobian.len() != measurements.len()
        || jacobian.iter().any(|row| row.len() != n)
    {
        return Err(invalid("nonlinear prediction/Jacobian shape"));
    }
    let mut gradient = z.to_vec();
    let mut objective = 0.5 * dot(z, z, progress)?;
    let mut residual = Vec::with_capacity(predictions.len());
    let mut transformed = Vec::with_capacity(predictions.len());
    for (i, (&prediction, row)) in predictions.iter().zip(&jacobian).enumerate() {
        finite(prediction, "nonlinear model prediction")?;
        hash_f64(hash, prediction, progress)?;
        for (component, &entry) in row.iter().enumerate() {
            if !entry.is_finite() {
                return Err(AssimError::NonFiniteObservationJacobian { component });
            }
            hash_f64(hash, entry, progress)?;
        }
        let sigma = measurements.variances[i].sqrt();
        let r = finite(
            (measurements.readings[i] - prediction) / sigma,
            "nonlinear standardized residual",
        )?;
        objective = finite(objective + 0.5 * r * r, "nonlinear objective")?;
        let mut j = vec![0.0; n];
        for column in 0..n {
            for k in 0..=column {
                j[column] = finite(
                    j[column] + row[k] * root[k][column],
                    "nonlinear whitened Jacobian",
                )?;
                progress.scalar("nonlinear Jacobian transform", 2)?;
            }
            j[column] = finite(j[column] / sigma, "nonlinear scaled Jacobian")?;
            gradient[column] = finite(gradient[column] - j[column] * r, "nonlinear gradient")?;
            progress.scalar("nonlinear gradient", 3)?;
        }
        residual.push(r);
        transformed.push(j);
    }
    Ok(Some(Point {
        state,
        residual,
        jacobian: transformed,
        gradient,
        objective,
    }))
}

fn affine_posterior(
    standard: &FactorBelief,
    z: &[f64],
    point: &Point,
    measurements: &IndependentMeasurements,
    cx: &Cx<'_>,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<FactorBelief, AssimError> {
    let mut observations = Vec::with_capacity(measurements.len());
    for (i, row) in point.jacobian.iter().enumerate() {
        progress.checkpoint("nonlinear affine row")?;
        // A zero Jacobian adds a constant likelihood term, already in the
        // objective. It cannot update the prior or be invented as sensitivity.
        if row.iter().all(|v| *v == 0.0) {
            continue;
        }
        let value = finite(
            point.residual[i] + dot(row, z, progress)?,
            "nonlinear affine target",
        )?;
        observations.push(Observation::new(
            row.clone(),
            value,
            1.0,
            measurements.instruments[i].clone(),
        )?);
    }
    if observations.is_empty() {
        progress.scalar("nonlinear unchanged factor", (z.len() * z.len()) as u128)?;
        Ok(standard.clone())
    } else {
        progress.scalar(
            "nonlinear factor sequence allowance",
            sequence_work(z.len(), observations.len())?,
        )?;
        Ok(assimilate_sequence(standard, &observations, cx)?
            .into_parts()
            .0)
    }
}

fn sequence_work(dimension: usize, rows: usize) -> Result<u128, AssimError> {
    let update = scalar_factor_work_estimate(dimension)?.update_with_receipt;
    Ok(rows as u128 * (update + 80 * dimension as u128 + 281))
}

fn dot(a: &[f64], b: &[f64], progress: &mut WorkProgress<'_, '_>) -> Result<f64, AssimError> {
    let mut value = 0.0;
    for (&a, &b) in a.iter().zip(b) {
        value = finite(a.mul_add(b, value), "nonlinear dot product")?;
        progress.scalar("nonlinear dot product", 1)?;
    }
    Ok(value)
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |a, b| a.max(b.abs()))
}
fn record(point: &Point, step_fraction: f64) -> IteratedRecord {
    IteratedRecord {
        objective: point.objective,
        gradient_norm: norm(&point.gradient),
        step_fraction,
    }
}
fn finite(value: f64, stage: &'static str) -> Result<f64, AssimError> {
    if value.is_finite() {
        Ok(canonicalize_zero(value))
    } else {
        Err(AssimError::NonFiniteComputation { stage })
    }
}
fn invalid(parameter: &'static str) -> AssimError {
    AssimError::InvalidLinearizationParameter { parameter }
}

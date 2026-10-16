//! Kernel learning with fixed, per-observation noise variances.
//!
//! Maximizes the exact dense GP log marginal likelihood in a declared log
//! parameter box. Noise is data, not an optimization variable. The derivative
//! is 1/2 tr[(alpha alpha^T - C^-1) dK], C = K + diag(noise), as in
//! Rasmussen & Williams (2006), eq. 5.9. No finite-difference refits are needed.
//! A bounded projected-gradient search retains the best evaluated model; it
//! does not assert global optimality or turn a posterior into a certificate.

use crate::gp::{Gp, Kernel, Matern};

/// Explicit training limits. Bounds are in physical units, NOT log units.
#[derive(Debug, Clone)]
pub struct HeteroFitConfig {
    /// One positive, finite (lower, upper) lengthscale interval per input.
    /// Equal endpoints freeze that parameter.
    pub lengthscale_bounds: Vec<(f64, f64)>,
    /// Positive signal-VARIANCE interval, in squared objective units.
    pub signal_bounds: (f64, f64),
    /// Starts INCLUDING the supplied initial kernel, followed by seeded starts.
    pub starts: usize,
    /// Maximum projected-gradient updates per start.
    pub max_iterations: usize,
    /// Total likelihood evaluations, INCLUDING initial, rejected and restart
    /// probes. Every evaluation makes at most one training factorization.
    pub max_evaluations: usize,
    /// Stop a start when the infinity norm of the unit projected log-gradient
    /// mapping is at most this value. This is not a global optimality test.
    pub gradient_tolerance: f64,
    /// Seed for deterministic restart points.
    pub seed: u64,
}

/// Admission, cancellation and initial-model failures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeteroFitError {
    /// Empty, ragged, non-finite or mismatched data; negative noise variance.
    InvalidData,
    /// Invalid kernel, bounds, initial point or work limits.
    InvalidConfig,
    /// The supplied initial covariance or likelihood/gradient is inadmissible.
    /// No noise floor is invented to repair it.
    InvalidInitialModel,
    /// The caller's continuation check refused further work.
    Cancelled,
}

impl std::fmt::Display for HeteroFitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "heteroscedastic GP training: {}", match self {
            Self::InvalidData => "invalid training data or observation variances",
            Self::InvalidConfig => "invalid kernel, parameter bounds or work limits",
            Self::InvalidInitialModel => "initial model is singular or numerically inadmissible",
            Self::Cancelled => "cancelled by caller",
        })
    }
}
impl std::error::Error for HeteroFitError {}

/// Best evaluated model plus actual work and local stationarity information.
/// No posterior fit is performed after consuming the evaluation allowance.
pub struct HeteroFitReport {
    /// Best finite evaluated GP; its LML cannot be lower than `initial_lml`.
    pub model: Gp,
    /// LML of the supplied initial kernel on THIS training set.
    pub initial_lml: f64,
    /// Actual likelihood probes, including failed trial factorizations.
    pub evaluations: usize,
    /// Starts actually evaluated, including the initial kernel.
    pub starts_attempted: usize,
    /// Accepted local updates, summed over starts.
    pub accepted_steps: usize,
    /// Projected log-gradient mapping at the RETURNED kernel, not at the
    /// last restart. A large value discloses an unfinished local fit.
    pub projected_gradient_inf: f64,
    /// Whether the exact likelihood-evaluation allowance was consumed.
    pub evaluation_limit_reached: bool,
}

impl HeteroFitConfig {
    /// Validate policy and the initial kernel without reading observations.
    pub fn validate(&self, initial: &Kernel) -> Result<(), HeteroFitError> {
        let valid_bound = |(lo, hi): (f64, f64)| {
            lo.is_finite() && hi.is_finite() && lo > 0.0 && lo <= hi
        };
        if initial.lengthscales.is_empty()
            || self.lengthscale_bounds.len() != initial.lengthscales.len()
            || self.starts == 0 || u32::try_from(self.starts).is_err()
            || self.max_evaluations == 0
            || !self.gradient_tolerance.is_finite() || self.gradient_tolerance < 0.0
            || !valid_bound(self.signal_bounds)
            || !initial.signal.is_finite()
            || !(self.signal_bounds.0..=self.signal_bounds.1).contains(&initial.signal)
            || !self.lengthscale_bounds.iter().copied().all(valid_bound)
            || initial.lengthscales.iter().zip(&self.lengthscale_bounds)
                .any(|(v, (lo, hi))| !v.is_finite() || *v < *lo || *v > *hi)
        {
            return Err(HeteroFitError::InvalidConfig);
        }
        Ok(())
    }

    fn log_bounds(&self) -> Vec<(f64, f64)> {
        self.lengthscale_bounds.iter().copied().chain([self.signal_bounds])
            .map(|(lo, hi)| (fs_math::det::ln(lo), fs_math::det::ln(hi))).collect()
    }

    fn kernel(&self, family: Matern, p: &[f64]) -> Kernel {
        let physical: Vec<f64> = p.iter().zip(
            self.lengthscale_bounds.iter().copied().chain([self.signal_bounds])
        ).map(|(v, (lo, hi))| fs_math::det::exp(*v).clamp(lo, hi)).collect();
        Kernel { family, signal: physical[physical.len() - 1],
            lengthscales: physical[..physical.len() - 1].to_vec() }
    }
}

fn validate_data(x: &[Vec<f64>], y: &[f64], noise: &[f64], dim: usize)
    -> Result<(), HeteroFitError>
{
    if x.is_empty() || x.len() != y.len() || x.len() != noise.len()
        || x.len().checked_mul(x.len()).is_none()
        || x.iter().any(|row| row.len() != dim || row.iter().any(|v| !v.is_finite()))
        || y.iter().any(|v| !v.is_finite())
        || noise.iter().any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(HeteroFitError::InvalidData);
    }
    Ok(())
}

fn checkpoint(keep_going: &mut dyn FnMut() -> bool) -> Result<(), HeteroFitError> {
    if keep_going() { Ok(()) } else { Err(HeteroFitError::Cancelled) }
}

// An exact duplicate with two noiseless observations makes C singular for
// EVERY kernel. Do not let a rounded positive Cholesky pivot admit it.
pub(crate) fn zero_noise_duplicates(x: &[Vec<f64>], noise: &[f64]) -> bool {
    noise.iter().enumerate().any(|(i, v)| *v == 0.0
        && (0..i).any(|j| noise[j] == 0.0 && x[i] == x[j]))
}

struct Evaluation {
    model: Gp,
    gradient: Vec<f64>,
}

// Differentiate log lengthscales and log SIGNAL VARIANCE (not amplitude).
// At coincident points every lengthscale derivative is exactly zero, including
// Matern-1/2: changing a lengthscale does not move the input coordinates.
fn kernel_gradient(kernel: &Kernel, x: &[f64], y: &[f64], out: &mut [f64]) {
    out.fill(0.0);
    let d = kernel.lengthscales.len();
    out[d] = kernel.eval(x, y);
    let mut r2 = 0.0f64;
    for ((a, b), ell) in x.iter().zip(y).zip(&kernel.lengthscales) {
        let z = (a - b) / ell;
        r2 = z.mul_add(z, r2);
    }
    let r = fs_math::det::sqrt(r2);
    // Match Kernel::eval's far-distance numerical branch exactly.
    if r == 0.0 || !r.is_finite() || r > 1e8 { return; }
    let radial = match kernel.family {
        Matern::Half => fs_math::det::exp(-r),
        Matern::ThreeHalves => {
            let a = fs_math::det::sqrt(3.0) * r;
            3.0 * r * fs_math::det::exp(-a)
        }
        Matern::FiveHalves => {
            let a = fs_math::det::sqrt(5.0) * r;
            (5.0 / 3.0) * r * (1.0 + a) * fs_math::det::exp(-a)
        }
    };
    for j in 0..d {
        let z = (x[j] - y[j]) / kernel.lengthscales[j];
        // z*(z/r) avoids division by a tiny r before multiplication by z^2.
        out[j] = kernel.signal * (radial * (z * (z / r)));
    }
}

fn evaluate(x: &[Vec<f64>], y: &[f64], noise: &[f64], kernel: Kernel,
    keep_going: &mut dyn FnMut() -> bool) -> Result<Option<Evaluation>, HeteroFitError>
{
    checkpoint(keep_going)?;
    let Some(model) = Gp::try_fit_diag(x, y, kernel, noise) else { return Ok(None); };
    if !model.lml.is_finite() { return Ok(None); }
    let (_, lower, alpha) = model.training_view();
    let n = x.len();
    let mut gradient = vec![0.0; model.kernel.lengthscales.len() + 1];
    let mut dk = gradient.clone();
    let mut column = vec![0.0f64; n];
    // Solve C * column = e_j using the ALREADY computed training Cholesky.
    // Accumulate the lower triangle of the trace; off-diagonal terms occur
    // twice in the symmetric trace and cancel its factor 1/2.
    for j in 0..n {
        checkpoint(keep_going)?;
        column.fill(0.0);
        column[j] = 1.0;
        for i in 0..n {
            let mut v = column[i];
            for k in 0..i { v = (-lower[i * n + k]).mul_add(column[k], v); }
            column[i] = v / lower[i * n + i];
        }
        for i in (0..n).rev() {
            let mut v = column[i];
            for k in i + 1..n { v = (-lower[k * n + i]).mul_add(column[k], v); }
            column[i] = v / lower[i * n + i];
        }
        for i in j..n {
            let weight = (alpha[i] * alpha[j] - column[i]) * if i == j { 0.5 } else { 1.0 };
            kernel_gradient(&model.kernel, &x[i], &x[j], &mut dk);
            for (g, derivative) in gradient.iter_mut().zip(&dk) {
                *g = weight.mul_add(*derivative, *g);
            }
        }
    }
    if gradient.iter().any(|v| !v.is_finite()) { return Ok(None); }
    Ok(Some(Evaluation { model, gradient }))
}

fn parameters(kernel: &Kernel) -> Vec<f64> {
    kernel.lengthscales.iter().copied().chain([kernel.signal]).map(fs_math::det::ln).collect()
}

fn projected_norm(p: &[f64], g: &[f64], bounds: &[(f64, f64)]) -> f64 {
    p.iter().zip(g).zip(bounds).map(|((v, grad), (lo, hi))|
        ((v + grad).clamp(*lo, *hi) - v).abs()).fold(0.0f64, f64::max)
}

/// Fit a fixed-noise GP from centered observations. The mean is NOT learned.
/// Uses the same policy as [`fit_heteroscedastic_controlled`] without a stop hook.
pub fn fit_heteroscedastic(x: &[Vec<f64>], y: &[f64], noise: &[f64], initial: &Kernel,
    config: &HeteroFitConfig) -> Result<HeteroFitReport, HeteroFitError>
{
    fit_heteroscedastic_controlled(x, y, noise, initial, config, &mut || true)
}

/// Bounded kernel learning with a caller-supplied continuation check.
///
/// Invalid data/policy refuses before numerical work. The initial kernel must
/// be admissible; trial singularities/non-finite likelihoods are rejected and
/// counted. Log-space projected Armijo steps and seeded restarts can only
/// replace the retained model with a higher measured LML. Frozen coordinates
/// remain frozen. Exhausting either cap returns the best model, not a claim
/// of convergence. Original noise variances are never rescaled or floored.
///
/// Checks occur before each fit, between inverse-column solves, between search
/// steps, and before returning. One dense Cholesky or triangular-column solve
/// is not preemptible. A caller may adapt its Cx checkpoint into this hook;
/// cancellation returns an error, never an apparent successful fit.
pub fn fit_heteroscedastic_controlled(x: &[Vec<f64>], y: &[f64], noise: &[f64], initial: &Kernel,
    config: &HeteroFitConfig, keep_going: &mut dyn FnMut() -> bool)
    -> Result<HeteroFitReport, HeteroFitError>
{
    config.validate(initial)?;
    validate_data(x, y, noise, initial.lengthscales.len())?;
    if zero_noise_duplicates(x, noise) { return Err(HeteroFitError::InvalidInitialModel); }
    let mut best = evaluate(x, y, noise, initial.clone(), keep_going)?
        .ok_or(HeteroFitError::InvalidInitialModel)?;
    let initial_lml = best.model.lml;
    let bounds = config.log_bounds();
    let mut evaluations = 1;
    let mut starts_attempted = 1;
    let mut accepted_steps = 0;
    let mut p = parameters(initial);
    let mut value = best.model.lml;
    let mut gradient = best.gradient.clone();
    for start in 0..config.starts {
        checkpoint(keep_going)?;
        if start > 0 {
            if evaluations == config.max_evaluations { break; }
            // Philox works for arbitrary ARD width; no fabricated QMC dims.
            let mut rng = fs_rand::StreamKey { seed: config.seed, kernel: 0x4845_5446,
                tile: u32::try_from(start).expect("admitted start count") }.stream();
            p = bounds.iter().map(|(lo, hi)| (hi - lo).mul_add(rng.next_f64(), *lo)).collect();
            evaluations += 1;
            starts_attempted += 1;
            let Some(scored) = evaluate(x, y, noise, config.kernel(initial.family, &p), keep_going)?
                else { continue; };
            value = scored.model.lml;
            gradient = scored.gradient.clone();
            p = parameters(&scored.model.kernel);
            if value > best.model.lml { best = scored; }
        }
        for _ in 0..config.max_iterations {
            checkpoint(keep_going)?;
            if evaluations == config.max_evaluations
                || projected_norm(&p, &gradient, &bounds) <= config.gradient_tolerance { break; }
            let active: Vec<f64> = p.iter().zip(&gradient).zip(&bounds)
                .map(|((v, g), (lo, hi))| {
                    if (*v <= *lo && *g < 0.0) || (*v >= *hi && *g > 0.0) {
                        0.0
                    } else { *g }
                }).collect();
            let scale = active.iter().map(|v| v.abs()).fold(1.0f64, f64::max);
            let direction: Vec<f64> = p.iter().zip(&active).zip(&bounds)
                .map(|((v, g), (lo, hi))| (v + g / scale).clamp(*lo, *hi) - v).collect();
            let mut step = 1.0f64;
            let mut accepted = false;
            for _ in 0..32 {
                if evaluations == config.max_evaluations { break; }
                let trial: Vec<f64> = p.iter().zip(&direction).zip(&bounds)
                    .map(|((v, d), (lo, hi))| step.mul_add(*d, *v).clamp(*lo, *hi)).collect();
                let slope: f64 = gradient.iter().zip(trial.iter().zip(&p))
                    .map(|(g, (a, b))| g * (a - b)).sum();
                if !slope.is_finite() || slope <= 0.0 { break; }
                evaluations += 1;
                if let Some(scored) = evaluate(x, y, noise, config.kernel(initial.family, &trial), keep_going)? {
                    let trial_value = scored.model.lml;
                    let sufficient = trial_value > value && trial_value - value >= 1e-4 * slope;
                    if sufficient {
                        p = parameters(&scored.model.kernel);
                        value = trial_value;
                        gradient = scored.gradient.clone();
                        accepted_steps += 1;
                        accepted = true;
                    }
                    // Even a finite rejected line-search probe is an evaluated
                    // model. Preserve a strict LML improvement before moving on.
                    if trial_value > best.model.lml { best = scored; }
                }
                if accepted { break; }
                step *= 0.5;
            }
            if !accepted { break; }
        }
    }
    checkpoint(keep_going)?;
    let projected_gradient_inf = projected_norm(&parameters(&best.model.kernel), &best.gradient, &bounds);
    Ok(HeteroFitReport { model: best.model, initial_lml, evaluations, starts_attempted,
        accepted_steps, projected_gradient_inf,
        evaluation_limit_reached: evaluations == config.max_evaluations })
}

#[cfg(test)]
mod tests;

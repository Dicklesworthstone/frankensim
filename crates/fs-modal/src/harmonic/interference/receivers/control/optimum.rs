//! Regularized optimal source gains for one physical interference window.
//!
//! Assemble the small real gain-space quadratic from original response bases,
//! retaining real/proper factor laws. The optimum solves (Q + tau W)x = -b,
//! tau = authored effort weight + nonnegative effort-limit multiplier.

use super::{
    ControlBank, HarmonicFactorStatistics, HarmonicReceiverControlEvaluation,
    HarmonicResponseError, MAX_RECEIVER_CONTROL_SOURCES, PreparedHarmonicReceiverControl, cap,
    checked, extent, invalid, observation, poll, rounded,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_la::factor::{Cholesky, FactorRunError, cholesky_cancellable};
use fs_math::{c64::C64, dd::Dd};

const MAX_VARIABLES: usize = 2 * MAX_RECEIVER_CONTROL_SOURCES;

/// Additional preparation/storage admission, excluding the borrowed controller.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiverOptimumLimits {
    /// Conservative retained arrays plus preparation/query factorization scratch.
    /// Allocator overhead, caller inputs and stack arrays are excluded.
    pub maximum_payload_bytes: usize,
    /// Maximum analytic scalar products in the gain quadratic:
    /// receivers * (d+1)*(d+2)/2 * sum_b(tones_b^2), d=2*sources.
    pub maximum_product_terms: usize,
}

/// A physical effort tradeoff and optional total effort cap.
///
/// Effort is E=x^T W x in the explicitly authored metric; x orders each gain's
/// real then imaginary component. This is a command/actuator effort model,
// not inferred mechanical input power or heat. The coefficient converts effort
/// units into the receiver objective's units. It is strictly positive and is
/// never silently inserted as a numerical ridge.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiverOptimumRequest {
    /// Positive coefficient lambda in J_receiver + lambda*E.
    pub effort_weight: f64,
    /// Optional strictly positive bound E <= budget, in the metric's units.
    pub maximum_effort: Option<f64>,
    /// Maximum gain-system factorizations; excludes the preparation's metric check.
    pub maximum_factorizations: usize,
    /// For a binding cap, maximum relative unused effort (budget-E)/budget.
    /// The published candidate always remains on the feasible side.
    pub relative_effort_tolerance: f64,
    /// Maximum original-field componentwise relative KKT stationarity residual.
    pub maximum_relative_stationarity_residual: f64,
}

/// A completed, independently checked source-gain optimum of the supplied model.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicReceiverOptimumSolution {
    /// Complex gains in the original source order and normalization.
    pub gains: Vec<C64>,
    /// Original receiver owner evaluated after field recomposition at these gains.
    pub receiver: HarmonicReceiverControlEvaluation,
    /// E=x^T W x, evaluated in the original authored metric.
    pub effort: f64,
    /// Original receiver objective plus the authored effort penalty lambda*E.
    pub regularized_objective: f64,
    /// Nonnegative multiplier mu of the total effort constraint; zero if unneeded.
    pub effort_multiplier: f64,
    /// Whether a positive multiplier was required to satisfy the effort cap.
    pub effort_limited: bool,
    /// Recomputed from the original physical receiver gradient, not Qx+b alone.
    /// This is backward/KKT evidence, not a forward-error or interval certificate.
    pub relative_stationarity_residual: f64,
    /// Gain-system factorizations used by this query.
    pub factorizations: usize,
}

/// Typed preparation, solve, accuracy or resource refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum HarmonicReceiverOptimumError {
    /// Existing source/window/cancellation/payload admission refused.
    Response(HarmonicResponseError),
    /// The authored effort metric failed the original equilibrated Cholesky.
    EffortNotPositiveDefinite,
    /// The regularized gain system could not be factored without alteration.
    RegularizedNotPositiveDefinite {
        /// Actual lambda+mu used by the failed factorization.
        effective_effort_weight: f64,
    },
    /// The explicit factorization budget ended before an acceptable optimum.
    FactorizationBudget {
        /// Completed factorizations.
        completed: usize,
        /// Caller-authored maximum.
        maximum: usize,
    },
    /// Binary64 multiplier resolution could not meet the authored effort tolerance.
    EffortResolution {
        /// Relative unused effort of the feasible endpoint.
        relative_deficit: f64,
    },
    /// The original physical-gradient KKT check exceeded its authored tolerance.
    Stationarity {
        /// Measured componentwise relative residual.
        measured: f64,
        /// Authored maximum.
        limit: f64,
    },
}

impl From<HarmonicResponseError> for HarmonicReceiverOptimumError {
    fn from(value: HarmonicResponseError) -> Self {
        Self::Response(value)
    }
}
impl core::fmt::Display for HarmonicReceiverOptimumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-RECEIVER-OPTIMUM: {self:?}")
    }
}
impl std::error::Error for HarmonicReceiverOptimumError {}

/// Cached physical gain-space quadratic and full effort metric for one window.
///
/// The mutable borrow keeps the original controller available for independent
/// field/gradient evaluation and prevents mismatched source replacement during
/// a query family. Different effort weights/budgets reuse this preparation.
///
/// This specialized convex harmonic solve uses existing cancellable Cholesky
/// factors and a bounded monotone multiplier bracket. It does not introduce a
/// generic optimization kernel, per-actuator saturation, nonlinear feedback or
/// source/actuator causality beyond the original harmonic model.
#[derive(Debug)]
pub struct PreparedHarmonicReceiverOptimum<'a> {
    controller: &'a mut PreparedHarmonicReceiverControl,
    window: HarmonicWindow,
    dimension: usize,
    // [bias, Re(g0), Im(g0), ...]. The constant bias/bias term is unused.
    gram: Vec<Dd>,
    effort: Vec<f64>,
    matrix: Vec<f64>,
    scales: Vec<f64>,
    rhs: Vec<f64>,
    dual_norm: f64,
}

impl PreparedHarmonicReceiverControl {
    /// Prepare optimal cancellation over this window and a full SPD effort law.
    ///
    /// effort_metric is a row-major (2*sources)^2 real, exactly symmetric
    /// matrix in [Re(g0), Im(g0), Re(g1), Im(g1), ...] order. It may retain
    /// cross-actuator and quadrature costs. The caller supplies its physical
    /// calibration and units; receiver response alone cannot determine effort.
    ///
    /// Source phases, shared factors, receiver weights and the mechanical pencil
    /// remain unchanged. The gain quadratic uses DD physical contractions, not
    /// differences of large objective values or finite-differenced gradients.
    pub fn prepare_optimum(
        &mut self,
        window: HarmonicWindow,
        effort_metric: &[f64],
        limits: HarmonicReceiverOptimumLimits,
        cx: &Cx<'_>,
    ) -> Result<PreparedHarmonicReceiverOptimum<'_>, HarmonicReceiverOptimumError> {
        poll(cx)?;
        let n = 2 * self.source_count;
        let width = n + 1;
        if effort_metric.len() != n * n || effort_metric.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite square source-effort metric required").into());
        }
        for i in 0..n {
            poll(cx)?;
            for j in 0..i {
                if effort_metric[i * n + j] != effort_metric[j * n + i] {
                    return Err(invalid("source-effort metric must be exactly symmetric").into());
                }
            }
            if effort_metric[i * n + i] <= 0.0 {
                return Err(HarmonicReceiverOptimumError::EffortNotPositiveDefinite);
            }
        }
        let mut pairs = 0usize;
        for bank in &self.banks {
            pairs = extent(
                bank.frequencies
                    .len()
                    .checked_mul(bank.frequencies.len())
                    .and_then(|x| pairs.checked_add(x)),
            )?;
        }
        let products = extent(
            pairs
                .checked_mul(self.weights.len())
                .and_then(|x| x.checked_mul(width * (width + 1) / 2)),
        )?;
        // n<=64. The original blocked factor's fixed GEMM packing arrays fit
        // within 2 MiB when a trailing panel exists (n>32). Six dense f64 arrays
        // cover retained metric/matrix, factor and panel/update temporaries.
        let packing = if n > 32 { 2 * 1024 * 1024 } else { 0 };
        let payload = width * width * core::mem::size_of::<Dd>()
            + (6 * n * n + 20 * n) * core::mem::size_of::<f64>()
            + packing
            + core::mem::size_of::<PreparedHarmonicReceiverOptimum<'_>>();
        cap(
            "receiver optimum products",
            products,
            limits.maximum_product_terms,
        )?;
        cap(
            "receiver optimum payload bytes",
            payload,
            limits.maximum_payload_bytes,
        )?;
        let mut scales = vec![0.0; n];
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            scales[i] = effort_metric[i * n + i].sqrt();
        }
        for i in 0..n {
            poll(cx)?;
            for j in 0..=i {
                let value = rounded(checked(
                    Dd::from_f64(effort_metric[i * n + j])
                        / Dd::from_f64(scales[i])
                        / Dd::from_f64(scales[j]),
                )?)?;
                matrix[i * n + j] = value;
                matrix[j * n + i] = value;
            }
        }
        let metric_factor = cholesky_cancellable(&matrix, n, || cx.checkpoint().is_err())
            .map_err(|error| factor_error(error, None))?;
        let mut gram = vec![Dd::ZERO; width * width];
        for bank in &self.banks {
            for (a, &omega_a) in bank.frequencies.iter().enumerate() {
                for (b, &omega_b) in bank.frequencies.iter().enumerate().skip(a) {
                    poll(cx)?;
                    match bank.statistics {
                        HarmonicFactorStatistics::Real => {
                            let kernel = window
                                .product_kernel(omega_a, omega_b)
                                .map_err(observation)?;
                            accumulate_gram(
                                bank,
                                a,
                                b,
                                &self.weights,
                                n,
                                &mut gram,
                                |x, y| kernel.mean_real_product_parts_dd(x, y).map_err(observation),
                                cx,
                            )?;
                        }
                        HarmonicFactorStatistics::ProperComplex => {
                            let kernel = window
                                .proper_product_kernel(omega_a, omega_b)
                                .map_err(observation)?;
                            accumulate_gram(
                                bank,
                                a,
                                b,
                                &self.weights,
                                n,
                                &mut gram,
                                |x, y| kernel.mean_real_product_parts_dd(x, y).map_err(observation),
                                cx,
                            )?;
                        }
                    }
                }
            }
        }
        for i in 0..width {
            for j in i + 1..width {
                gram[j * width + i] = gram[i * width + j];
            }
        }
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            rhs[i] = rounded(checked(gram[i + 1] / Dd::from_f64(scales[i]))?)?;
        }
        solve_factor(&metric_factor, &mut rhs, cx)?;
        let mut dual_square = Dd::ZERO;
        for i in 0..n {
            dual_square = checked(
                dual_square
                    + checked(gram[i + 1] / Dd::from_f64(scales[i]))? * Dd::from_f64(rhs[i]),
            )?;
        }
        if dual_square.lt(Dd::ZERO) {
            return Err(invalid("source-effort dual norm is not representable").into());
        }
        let dual_norm = rounded(checked(dual_square.sqrt())?)?;
        poll(cx)?;
        Ok(PreparedHarmonicReceiverOptimum {
            controller: self,
            window,
            dimension: n,
            gram,
            effort: effort_metric.to_vec(),
            matrix,
            scales,
            rhs,
            dual_norm,
        })
    }
}

#[derive(Clone, Copy)]
struct Trial {
    gains: [C64; MAX_RECEIVER_CONTROL_SOURCES],
    effort: Dd,
}

impl PreparedHarmonicReceiverOptimum<'_> {
    /// Solve the fixed physical problem with an explicit positive effort cost.
    ///
    /// Without a binding cap, one factorization gives the unique regularized
    /// optimum. A binding total cap uses the monotone effort of
    /// (Q+(lambda+mu)W)x=-b. The bracket always retains a feasible endpoint;
    /// exhausting work or effort resolution returns no approximate solution.
    ///
    /// Final vibration and its gradient come from the original field owner.
    /// Effort and KKT stationarity use the original authored metric. Preparation
    /// and query use ordinary floating linear algebra with DD contractions;
    /// their residuals do not constitute an interval optimality certificate.
    pub fn solve(
        &mut self,
        request: HarmonicReceiverOptimumRequest,
        cx: &Cx<'_>,
    ) -> Result<HarmonicReceiverOptimumSolution, HarmonicReceiverOptimumError> {
        poll(cx)?;
        if !request.effort_weight.is_finite()
            || request.effort_weight <= 0.0
            || request
                .maximum_effort
                .is_some_and(|x| !x.is_finite() || x <= 0.0)
            || request.maximum_factorizations == 0
            || !request.relative_effort_tolerance.is_finite()
            || request.relative_effort_tolerance <= 0.0
            || request.relative_effort_tolerance >= 1.0
            || !request.maximum_relative_stationarity_residual.is_finite()
            || request.maximum_relative_stationarity_residual <= 0.0
            || request.maximum_relative_stationarity_residual >= 1.0
        {
            return Err(invalid(
                "positive authored effort cost, cap, work and residual tolerances required",
            )
            .into());
        }
        let mut factorizations = 0;
        let mut effective = request.effort_weight;
        let mut best = self.solve_at(
            effective,
            &mut factorizations,
            request.maximum_factorizations,
            cx,
        )?;
        let mut limited = false;
        if let Some(budget) = request.maximum_effort {
            let bound = Dd::from_f64(budget);
            if bound.lt(best.effort) {
                limited = true;
                let mut lower = effective;
                let mut upper = (self.dual_norm / budget.sqrt()).max(effective);
                if upper <= lower {
                    upper = lower * 2.0;
                }
                loop {
                    if !upper.is_finite() {
                        return Err(
                            invalid("effort multiplier bracket is not representable").into()
                        );
                    }
                    best = self.solve_at(
                        upper,
                        &mut factorizations,
                        request.maximum_factorizations,
                        cx,
                    )?;
                    if !bound.lt(best.effort) {
                        break;
                    }
                    lower = upper;
                    upper *= 2.0;
                }
                loop {
                    let deficit = rounded(checked((bound - best.effort) / bound)?)?;
                    if deficit <= request.relative_effort_tolerance {
                        break;
                    }
                    let middle = lower + (upper - lower) * 0.5;
                    if middle <= lower || middle >= upper {
                        return Err(HarmonicReceiverOptimumError::EffortResolution {
                            relative_deficit: deficit,
                        });
                    }
                    let trial = self.solve_at(
                        middle,
                        &mut factorizations,
                        request.maximum_factorizations,
                        cx,
                    )?;
                    if bound.lt(trial.effort) {
                        lower = middle;
                    } else {
                        upper = middle;
                        best = trial;
                    }
                }
                effective = upper;
            }
        }
        let count = self.controller.source_count;
        let mut gradient = [C64::ZERO; MAX_RECEIVER_CONTROL_SOURCES];
        let receiver = self.controller.evaluate_into(
            &best.gains[..count],
            self.window,
            &mut gradient[..count],
            cx,
        )?;
        let residual = self.stationarity(&best.gains, &gradient, effective, cx)?;
        if residual > request.maximum_relative_stationarity_residual {
            return Err(HarmonicReceiverOptimumError::Stationarity {
                measured: residual,
                limit: request.maximum_relative_stationarity_residual,
            });
        }
        let effort = rounded(best.effort)?;
        let objective = rounded(checked(
            Dd::from_f64(receiver.mean_square) + Dd::from_f64(request.effort_weight) * best.effort,
        )?)?;
        poll(cx)?;
        Ok(HarmonicReceiverOptimumSolution {
            gains: best.gains[..count].to_vec(),
            receiver,
            effort,
            regularized_objective: objective,
            effort_multiplier: effective - request.effort_weight,
            effort_limited: limited,
            relative_stationarity_residual: residual,
            factorizations,
        })
    }

    fn solve_at(
        &mut self,
        effective: f64,
        completed: &mut usize,
        maximum: usize,
        cx: &Cx<'_>,
    ) -> Result<Trial, HarmonicReceiverOptimumError> {
        poll(cx)?;
        if *completed == maximum {
            return Err(HarmonicReceiverOptimumError::FactorizationBudget {
                completed: *completed,
                maximum,
            });
        }
        let n = self.dimension;
        let width = n + 1;
        let tau = Dd::from_f64(effective);
        for i in 0..n {
            let diagonal = checked(
                self.gram[(i + 1) * width + i + 1] + tau * Dd::from_f64(self.effort[i * n + i]),
            )?;
            if !Dd::ZERO.lt(diagonal) {
                return Err(
                    HarmonicReceiverOptimumError::RegularizedNotPositiveDefinite {
                        effective_effort_weight: effective,
                    },
                );
            }
            self.scales[i] = rounded(checked(diagonal.sqrt())?)?;
        }
        for i in 0..n {
            poll(cx)?;
            for j in 0..=i {
                let entry = checked(
                    self.gram[(i + 1) * width + j + 1] + tau * Dd::from_f64(self.effort[i * n + j]),
                )?;
                let value = rounded(checked(
                    entry / Dd::from_f64(self.scales[i]) / Dd::from_f64(self.scales[j]),
                )?)?;
                self.matrix[i * n + j] = value;
                self.matrix[j * n + i] = value;
            }
            self.rhs[i] = rounded(checked(-self.gram[i + 1] / Dd::from_f64(self.scales[i]))?)?;
        }
        let factor = cholesky_cancellable(&self.matrix, n, || cx.checkpoint().is_err())
            .map_err(|error| factor_error(error, Some(effective)))?;
        *completed += 1;
        solve_factor(&factor, &mut self.rhs, cx)?;
        let mut gains = [C64::ZERO; MAX_RECEIVER_CONTROL_SOURCES];
        for i in 0..n {
            let value = rounded(checked(
                Dd::from_f64(self.rhs[i]) / Dd::from_f64(self.scales[i]),
            )?)?;
            if i % 2 == 0 {
                gains[i / 2].re = value;
            } else {
                gains[i / 2].im = value;
            }
        }
        let mut effort = Dd::ZERO;
        for i in 0..n {
            poll(cx)?;
            for j in i..n {
                let term = checked(
                    Dd::from_f64(self.effort[i * n + j])
                        * Dd::from_f64(component(&gains, i))
                        * Dd::from_f64(component(&gains, j)),
                )?;
                effort = checked(effort + term * Dd::from_f64(if i == j { 1.0 } else { 2.0 }))?;
            }
        }
        if effort.lt(Dd::ZERO) {
            return Err(
                invalid("original actuator effort is negative after numerical solve").into(),
            );
        }
        Ok(Trial { gains, effort })
    }

    fn stationarity(
        &self,
        gains: &[C64],
        gradient: &[C64],
        effective: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, HarmonicResponseError> {
        let n = self.dimension;
        let width = n + 1;
        let tau = Dd::from_f64(effective);
        let mut maximum: f64 = 0.0;
        for i in 0..n {
            poll(cx)?;
            let mut wx = Dd::ZERO;
            let mut scale = self.gram[i + 1].abs();
            for j in 0..n {
                let x = Dd::from_f64(component(gains, j));
                let w = Dd::from_f64(self.effort[i * n + j]);
                wx = checked(wx + w * x)?;
                scale = checked(
                    scale
                        + checked(self.gram[(i + 1) * width + j + 1].abs() + tau * w.abs())?
                            * x.abs(),
                )?;
            }
            let residual =
                checked(Dd::from_f64(component(gradient, i)) + Dd::from_f64(2.0) * tau * wx)?.abs();
            let denominator = checked(Dd::from_f64(2.0) * scale)?;
            let ratio = if denominator == Dd::ZERO {
                if residual == Dd::ZERO {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                rounded(checked(residual / denominator)?)?
            };
            maximum = maximum.max(ratio);
        }
        Ok(maximum)
    }
}

#[allow(clippy::too_many_arguments)] // One physical tone pair and its gain-space contraction.
fn accumulate_gram(
    bank: &ControlBank,
    a: usize,
    b: usize,
    weights: &[f64],
    n: usize,
    gram: &mut [Dd],
    product: impl Fn([Dd; 2], [Dd; 2]) -> Result<Dd, HarmonicResponseError>,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    let width = n + 1;
    let stride = weights.len() * (n / 2 + 1);
    let mut left = [[Dd::ZERO; 2]; MAX_VARIABLES + 1];
    let mut right = [[Dd::ZERO; 2]; MAX_VARIABLES + 1];
    for (row, &weight) in weights.iter().enumerate() {
        poll(cx)?;
        if weight == 0.0 {
            continue;
        }
        left[0] = bank.responses[a * stride + row];
        right[0] = bank.responses[b * stride + row];
        for source in 0..n / 2 {
            let l = bank.responses[a * stride + (source + 1) * weights.len() + row];
            let r = bank.responses[b * stride + (source + 1) * weights.len() + row];
            left[2 * source + 1] = l;
            left[2 * source + 2] = [-l[1], l[0]];
            right[2 * source + 1] = r;
            right[2 * source + 2] = [-r[1], r[0]];
        }
        let weight = Dd::from_f64(weight);
        for i in 0..width {
            poll(cx)?;
            for j in i.max(1)..width {
                let mut value = product(left[i], right[j])?;
                // Swapping tone order also swaps the kernel's frequency slots.
                // Keep each response in its own slot when using the (a,b) kernel.
                if a != b {
                    value = checked(value + product(left[j], right[i])?)?;
                }
                gram[i * width + j] = checked(gram[i * width + j] + weight * value)?;
            }
        }
    }
    Ok(())
}

fn component(values: &[C64], i: usize) -> f64 {
    if i % 2 == 0 {
        values[i / 2].re
    } else {
        values[i / 2].im
    }
}

fn factor_error(error: FactorRunError, effective: Option<f64>) -> HarmonicReceiverOptimumError {
    match error {
        FactorRunError::Cancelled => HarmonicResponseError::Cancelled.into(),
        FactorRunError::Factor(_) => match effective {
            None => HarmonicReceiverOptimumError::EffortNotPositiveDefinite,
            Some(effective_effort_weight) => {
                HarmonicReceiverOptimumError::RegularizedNotPositiveDefinite {
                    effective_effort_weight,
                }
            }
        },
    }
}

fn solve_factor(
    factor: &Cholesky,
    rhs: &mut [f64],
    cx: &Cx<'_>,
) -> Result<(), HarmonicReceiverOptimumError> {
    factor
        .solve_cancellable(rhs, || cx.checkpoint().is_err())
        .map_err(|error| factor_error(error, None))
}

#[cfg(test)]
mod tests;

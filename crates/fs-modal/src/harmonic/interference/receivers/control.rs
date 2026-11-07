//! Exact source-gain derivatives of finite-window physical vibration.
//!
//! All response bases use the same declared random-factor coordinates.
//! Recompose their complex amplitudes before products. Differentiating that
//! residual field avoids subtracting large independently squared source powers.

use super::super::{
    HarmonicFactorStatistics, HarmonicResponseError, MAX_INTERFERENCE_FACTORS, cap, checked,
    finite, invalid, observation, parts, poll,
};
use super::PreparedHarmonicReceivers;
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd};

/// Hard ceiling on independently controlled complex source gains.
pub const MAX_RECEIVER_CONTROL_SOURCES: usize = 32;

/// Admission for a fixed physical receiver objective and its response bases.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiverControlLimits {
    /// Maximum controllable sources, bounded by MAX_RECEIVER_CONTROL_SOURCES.
    pub maximum_sources: usize,
    /// Maximum retained bank-frequency pairs after exact frequency coalescing.
    pub maximum_tones: usize,
    /// Owned response/mixed arrays, weights, gradients, and vector headers.
    /// Excludes input owners, allocator overhead, and stack variables.
    pub maximum_payload_bytes: usize,
    /// Conservative scalar products per evaluation:
    /// receivers * (1 + 2*sources) * sum_b(distinct_frequencies_b^2).
    pub maximum_product_terms: usize,
}

/// Expected weighted mean-square motion over the actual observation window.
///
/// The objective is <E[sum_r w_r y_r(t)^2]>. It is not the mean square of the
/// window-averaged reading. Weights retain caller-authored physical units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicReceiverControlEvaluation {
    /// Contribution of the controlled deterministic mean trajectory.
    pub coherent_mean_square: f64,
    /// Contribution of the controlled real/proper random factors.
    pub stochastic_mean_square: f64,
    /// Sum of coherent and stochastic contributions, rounded only at the end.
    /// Signed floating roundoff is retained, not clipped.
    pub mean_square: f64,
}

#[derive(Clone, Debug)]
struct ControlBank {
    statistics: HarmonicFactorStatistics,
    frequencies: Vec<f64>,
    // Frequency-major, then bias/source, then physical receiver.
    responses: Vec<[Dd; 2]>,
    // Frequency-major, then receiver. Rebuilt completely for every evaluation.
    mixed: Vec<[Dd; 2]>,
}

/// Prepared gain/phase control of already solved physical receiver responses.
///
/// The motion is the bias plus sum_j gain_j * source_j in complex phasors.
/// Each source may contain multiple exact frequencies and every declared
/// factor. A separate source basis can represent an independently controlled
/// actuator/frequency. Gains do not retune the mechanical pencil or its loads.
///
/// Inputs are projections of existing harmonic responses. Preparation performs
/// no solve, and evaluation uses no allocation or original-coordinate matrix
/// operation. DD field recomposition precedes analytic window products, even
/// when several source bases share a frequency. All phase and beat terms remain.
///
/// This is a linear harmonic source-control model on the supplied common
/// physical coordinates/clock. A complex gain on a finite tone bank does not
/// certify a causal broadband actuator/filter, a startup transient, actuator
/// constraints, source statistics, or relative accuracy at a null.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicReceiverControl {
    source_count: usize,
    weights: Vec<f64>,
    banks: Vec<ControlBank>,
    gradient: Vec<[Dd; 2]>,
    maximum_relative_response_residual: f64,
}

impl PreparedHarmonicReceiverControl {
    /// Bind a bias and controllable source bases on explicitly SHARED factors.
    ///
    /// Receiver order/units and original physical pencil/clock must agree.
    /// Every input must have the same number of banks and the same real/proper
    /// law at each index. Bank zero is the deterministic mean; random bank k
    /// means the SAME random variable in every input. Calling this constructor
    /// declares that joint model; index coincidence alone cannot establish it.
    ///
    /// Independent sources require distinct shared factor slots, with explicit
    /// empty factors in the other source preparations. Mean-only inputs need
    /// no random slots. Factor phases and square-root statistical weights remain
    /// in their original coefficients. No Gaussian or fourth-moment law is used.
    ///
    /// Receiver weights must be finite and nonnegative, with at least one
    /// positive entry. A strain/energy objective uses the original physical
    /// receiver map and its quadrature/material weights, not modal magnitudes.
    pub fn try_with_shared_factors(
        bias: &PreparedHarmonicReceivers,
        sources: &[&PreparedHarmonicReceivers],
        receiver_weights: &[f64],
        limits: HarmonicReceiverControlLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        if sources.is_empty()
            || limits.maximum_sources == 0
            || limits.maximum_sources > MAX_RECEIVER_CONTROL_SOURCES
        {
            return Err(invalid(
                "nonempty source controls and bounded source cap required",
            ));
        }
        cap(
            "receiver control sources",
            sources.len(),
            limits.maximum_sources,
        )?;
        let receivers = bias.receiver_count;
        let bank_count = bias.banks.len();
        if bank_count == 0
            || bank_count > MAX_INTERFERENCE_FACTORS + 1
            || receiver_weights.len() != receivers
            || receiver_weights.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !receiver_weights.iter().any(|x| *x > 0.0)
        {
            return Err(invalid(
                "shared receiver rows and a finite positive observation metric required",
            ));
        }
        let mut inputs = [bias; MAX_RECEIVER_CONTROL_SOURCES + 1];
        let mut maximum_relative_response_residual = bias.maximum_relative_response_residual;
        for (i, &source) in sources.iter().enumerate() {
            poll(cx)?;
            if source.receiver_count != receivers
                || source.banks.len() != bank_count
                || source
                    .banks
                    .iter()
                    .zip(&bias.banks)
                    .any(|(a, b)| a.statistics != b.statistics)
            {
                return Err(invalid(
                    "source bases must use the same receiver order and explicit factor laws",
                ));
            }
            maximum_relative_response_residual =
                maximum_relative_response_residual.max(source.maximum_relative_response_residual);
            inputs[i + 1] = source;
        }
        let inputs = &inputs[..sources.len() + 1];
        let mut counts = [0usize; MAX_INTERFERENCE_FACTORS + 1];
        let mut total_tones = 0usize;
        let mut pair_count = 0usize;
        for (bank, count) in counts.iter_mut().enumerate().take(bank_count) {
            for_frequency(inputs, bank, cx, |_, _| {
                *count = extent(count.checked_add(1))?;
                Ok(())
            })?;
            total_tones = extent(total_tones.checked_add(*count))?;
            pair_count = extent(
                count
                    .checked_mul(*count)
                    .and_then(|x| pair_count.checked_add(x)),
            )?;
        }
        cap("receiver control tones", total_tones, limits.maximum_tones)?;
        let product_terms = extent(
            sources
                .len()
                .checked_mul(2)
                .and_then(|x| x.checked_add(1))
                .and_then(|x| x.checked_mul(receivers))
                .and_then(|x| x.checked_mul(pair_count)),
        )?;
        let payload = extent(
            sources
                .len()
                .checked_add(2)
                .and_then(|x| x.checked_mul(receivers))
                .and_then(|x| x.checked_mul(core::mem::size_of::<[Dd; 2]>()))
                .and_then(|x| x.checked_add(core::mem::size_of::<f64>()))
                .and_then(|x| x.checked_mul(total_tones))
                .and_then(|x| {
                    x.checked_add(bank_count.checked_mul(core::mem::size_of::<ControlBank>())?)
                })
                .and_then(|x| x.checked_add(receivers.checked_mul(core::mem::size_of::<f64>())?))
                .and_then(|x| {
                    x.checked_add(sources.len().checked_mul(core::mem::size_of::<[Dd; 2]>())?)
                })
                .and_then(|x| x.checked_add(core::mem::size_of::<Self>())),
        )?;
        cap(
            "receiver control payload bytes",
            payload,
            limits.maximum_payload_bytes,
        )?;
        cap(
            "receiver control product terms",
            product_terms,
            limits.maximum_product_terms,
        )?;

        let stride = receivers * inputs.len(); // Extents were admitted above.
        let mut banks = Vec::with_capacity(bank_count);
        for (bank, &count) in counts.iter().enumerate().take(bank_count) {
            let mut frequencies = Vec::with_capacity(count);
            let mut responses = Vec::with_capacity(count * stride);
            for_frequency(inputs, bank, cx, |omega, positions| {
                frequencies.push(omega);
                for (input, &position) in inputs.iter().zip(positions) {
                    let tone = input.banks[bank]
                        .tones
                        .get(position)
                        .filter(|t| t.omega == omega);
                    if let Some(tone) = tone {
                        responses.extend_from_slice(&tone.amplitudes);
                    } else {
                        responses.resize(responses.len() + receivers, [Dd::ZERO; 2]);
                    }
                }
                Ok(())
            })?;
            banks.push(ControlBank {
                statistics: bias.banks[bank].statistics,
                frequencies,
                responses,
                mixed: vec![[Dd::ZERO; 2]; count * receivers],
            });
        }
        poll(cx)?;
        Ok(Self {
            source_count: sources.len(),
            weights: receiver_weights.to_vec(),
            banks,
            gradient: vec![[Dd::ZERO; 2]; sources.len()],
            maximum_relative_response_residual,
        })
    }

    /// Number of complex source gains and required gradient entries.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }

    /// Largest original response backward residual among all input owners.
    /// This is not a forward/relative bound on a nearly cancelled reading.
    #[must_use]
    pub const fn maximum_relative_response_residual(&self) -> f64 {
        self.maximum_relative_response_residual
    }

    /// Evaluate finite-window vibration and its exact Cartesian gain gradient.
    ///
    /// gradient[j].re = dJ/dRe(gain_j), gradient[j].im = dJ/dIm(gain_j).
    /// This pair represents two REAL derivatives, not a holomorphic derivative.
    /// For gain=A*exp(i*phi), dJ/dphi =
    /// -gain.im*gradient.re + gain.re*gradient.im. That chain rule retains the
    /// original exp(-i omega t) convention; no phase is averaged away.
    ///
    /// The derivative contracts twice the controlled residual field with the
    /// corresponding source response, under each bank's declared real/proper
    /// law. It does not differentiate a solver iteration or finite difference.
    /// Source phases, frequencies, receiver weights and the pencil are fixed.
    ///
    /// A refusal leaves the caller's gradient unchanged. Private scratch is
    /// completely rebuilt on the next call, including after cancellation.
    pub fn evaluate_into(
        &mut self,
        gains: &[C64],
        window: HarmonicWindow,
        gradient: &mut [C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicReceiverControlEvaluation, HarmonicResponseError> {
        poll(cx)?;
        if gains.len() != self.source_count
            || gradient.len() != self.source_count
            || gains.iter().any(|&gain| !finite(gain))
        {
            return Err(invalid(
                "one finite complex gain and gradient slot per source required",
            ));
        }
        let mut gains_dd = [[Dd::ZERO; 2]; MAX_RECEIVER_CONTROL_SOURCES];
        for (target, &gain) in gains_dd.iter_mut().zip(gains) {
            *target = parts(gain);
        }
        let receivers = self.weights.len();
        let stride = receivers * (self.source_count + 1);
        for bank in &mut self.banks {
            for (tone, current) in bank.mixed.chunks_exact_mut(receivers).enumerate() {
                poll(cx)?;
                let responses = &bank.responses[tone * stride..(tone + 1) * stride];
                for (row, value) in current.iter_mut().enumerate() {
                    if row % 32 == 0 {
                        poll(cx)?;
                    }
                    let mut sum = responses[row];
                    for (source, gain) in gains_dd.iter().enumerate().take(self.source_count) {
                        let amplitude = responses[(source + 1) * receivers + row];
                        sum[0] = checked(
                            sum[0] + checked(gain[0] * amplitude[0] - gain[1] * amplitude[1])?,
                        )?;
                        sum[1] = checked(
                            sum[1] + checked(gain[0] * amplitude[1] + gain[1] * amplitude[0])?,
                        )?;
                    }
                    *value = sum;
                }
            }
        }
        self.gradient.fill([Dd::ZERO; 2]);
        let mut coherent = Dd::ZERO;
        let mut stochastic = Dd::ZERO;
        for (index, bank) in self.banks.iter().enumerate() {
            let mut moment = Dd::ZERO;
            for (a, &omega_a) in bank.frequencies.iter().enumerate() {
                for (b, &omega_b) in bank.frequencies.iter().enumerate().skip(a) {
                    poll(cx)?;
                    match bank.statistics {
                        HarmonicFactorStatistics::Real => {
                            let kernel = window
                                .product_kernel(omega_a, omega_b)
                                .map_err(observation)?;
                            accumulate_pair(
                                bank,
                                a,
                                b,
                                &self.weights,
                                &mut self.gradient,
                                &mut moment,
                                |x, y| kernel.mean_real_product_parts_dd(x, y).map_err(observation),
                                cx,
                            )?;
                        }
                        HarmonicFactorStatistics::ProperComplex => {
                            let kernel = window
                                .proper_product_kernel(omega_a, omega_b)
                                .map_err(observation)?;
                            accumulate_pair(
                                bank,
                                a,
                                b,
                                &self.weights,
                                &mut self.gradient,
                                &mut moment,
                                |x, y| kernel.mean_real_product_parts_dd(x, y).map_err(observation),
                                cx,
                            )?;
                        }
                    }
                }
            }
            if index == 0 {
                coherent = moment;
            } else {
                stochastic = checked(stochastic + moment)?;
            }
        }
        let result = HarmonicReceiverControlEvaluation {
            coherent_mean_square: rounded(coherent)?,
            stochastic_mean_square: rounded(stochastic)?,
            mean_square: rounded(checked(coherent + stochastic)?)?,
        };
        let mut candidate = [C64::ZERO; MAX_RECEIVER_CONTROL_SOURCES];
        for (target, &derivative) in candidate.iter_mut().zip(&self.gradient) {
            *target = C64::new(rounded(derivative[0])?, rounded(derivative[1])?);
        }
        poll(cx)?;
        gradient.copy_from_slice(&candidate[..self.source_count]);
        Ok(result)
    }
}

// A merge over already sorted exact frequencies. Its bounded cursor array
// avoids allocating ordering records or imposing a nearest-frequency tolerance.
fn for_frequency(
    inputs: &[&PreparedHarmonicReceivers],
    bank: usize,
    cx: &Cx<'_>,
    mut visit: impl FnMut(f64, &[usize]) -> Result<(), HarmonicResponseError>,
) -> Result<(), HarmonicResponseError> {
    let mut positions = [0usize; MAX_RECEIVER_CONTROL_SOURCES + 1];
    loop {
        poll(cx)?;
        let mut next: Option<f64> = None;
        for (input, &position) in inputs.iter().zip(&positions) {
            if let Some(tone) = input.banks[bank].tones.get(position) {
                if next.is_none_or(|omega| tone.omega < omega) {
                    next = Some(tone.omega);
                }
            }
        }
        let Some(omega) = next else {
            return Ok(());
        };
        visit(omega, &positions[..inputs.len()])?;
        for (input, position) in inputs.iter().zip(&mut positions) {
            if input.banks[bank]
                .tones
                .get(*position)
                .is_some_and(|t| t.omega == omega)
            {
                *position += 1;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)] // One physical tone pair and its gradient accumulator.
fn accumulate_pair(
    bank: &ControlBank,
    a: usize,
    b: usize,
    weights: &[f64],
    gradient: &mut [[Dd; 2]],
    moment: &mut Dd,
    product: impl Fn([Dd; 2], [Dd; 2]) -> Result<Dd, HarmonicResponseError>,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    let receivers = weights.len();
    let stride = receivers * (gradient.len() + 1);
    let multiplicity = Dd::from_f64(if a == b { 1.0 } else { 2.0 });
    for (row, &weight) in weights.iter().enumerate() {
        if row % 32 == 0 {
            poll(cx)?;
        }
        if weight == 0.0 {
            continue;
        }
        let weight = Dd::from_f64(weight);
        let left = bank.mixed[a * receivers + row];
        let right = bank.mixed[b * receivers + row];
        *moment = checked(*moment + checked(weight * product(left, right)?)? * multiplicity)?;
        for (source, derivative) in gradient.iter_mut().enumerate() {
            let right_source = bank.responses[b * stride + (source + 1) * receivers + row];
            let mut real = product(left, right_source)?;
            let mut imaginary = product(left, [-right_source[1], right_source[0]])?;
            if a != b {
                let left_source = bank.responses[a * stride + (source + 1) * receivers + row];
                real = checked(real + product(left_source, right)?)?;
                imaginary =
                    checked(imaginary + product([-left_source[1], left_source[0]], right)?)?;
            }
            derivative[0] = checked(derivative[0] + checked(weight * real)? * Dd::from_f64(2.0))?;
            derivative[1] =
                checked(derivative[1] + checked(weight * imaginary)? * Dd::from_f64(2.0))?;
        }
    }
    Ok(())
}

fn extent(value: Option<usize>) -> Result<usize, HarmonicResponseError> {
    value.ok_or(invalid("receiver control extent overflow"))
}

fn rounded(value: Dd) -> Result<f64, HarmonicResponseError> {
    let value = value.to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("controlled receiver result is not representable"))
    }
}

#[cfg(test)]
mod tests;

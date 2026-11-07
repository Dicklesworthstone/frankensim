//! Physical receiver rows applied before finite-window mechanical products.
//!
//! Project the retained response, preserving its factor, tone, and phase.
//! The window covariance of the averaged reading is distinct from the window
//! average of the instantaneous covariance; both are evaluated directly.

use super::{
    HarmonicFactorStatistics, HarmonicResponseError, PreparedHarmonicInterference, Tone, cap,
    checked, invalid, observation, parts, poll,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::dd::Dd;

mod control;
pub use control::{
    HarmonicReceiverOptimumError, HarmonicReceiverOptimumLimits,
    HarmonicReceiverOptimumRequest, HarmonicReceiverOptimumSolution,
    PreparedHarmonicReceiverOptimum,
    HarmonicReceiverControlEvaluation, HarmonicReceiverControlLimits,
    MAX_RECEIVER_CONTROL_SOURCES, PreparedHarmonicReceiverControl,
};

/// Hard ceiling on simultaneously prepared physical receiver rows.
pub const MAX_INTERFERENCE_RECEIVERS: usize = 128;

/// Which motion the real receiver row observes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicReceiverQuantity {
    /// Original generalized displacement, including static displacement at DC.
    Displacement,
    /// Original generalized velocity under exp(-i omega t).
    Velocity,
    /// Original generalized acceleration, differentiated from the velocity.
    Acceleration,
}

/// One real, frequency-independent linear physical observation.
///
/// The row acts on the original pencil coordinates, not on modal magnitudes.
/// It may describe a point, an area average, a differential pickup, or a linear
/// strain/stress recovery. Its coefficients and selected motion determine
/// output units; no normalization, modal-shape interpolation, or material
/// recovery is inferred. Frequency-dependent receiver dynamics require their
/// own physical transfer owner.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiver<'a> {
    /// One finite signed coefficient per original pencil coordinate.
    pub weights: &'a [f64],
    /// Motion to which the row is applied.
    pub quantity: HarmonicReceiverQuantity,
}

/// Additional storage and work admitted when projecting an existing response.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReceiverControls {
    /// Maximum receiver rows, bounded by MAX_INTERFERENCE_RECEIVERS.
    pub maximum_receivers: usize,
    /// Retained projected amplitudes and headers; excludes the source owner,
    /// caller rows, allocator overhead, and stack variables.
    pub maximum_payload_bytes: usize,
    /// Maximum complex weighted-coordinate terms: n * receivers * total tones.
    pub maximum_projection_terms: usize,
    /// Conservative scalar products per pair observation:
    /// sum_b(L_b^2) + 4*sum_b(L_b) + 2*banks.
    pub maximum_product_terms: usize,
}

/// Joint moments of two real readings under one normalized physical window.
///
/// Brackets denote the supplied time weighting; E denotes the explicitly
/// authored ensemble. Using the same receiver twice gives mean-square motion
/// and variances. Different receivers give signed cross moments, including
/// displacement/velocity or opposite-polarity pickups. No negative roundoff is
/// clipped. Output product units are the product of the two receiver units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicReceiverPairMoments {
    /// <E[left(t)]>. Stochastic zero-mean factors do not shift this reading.
    pub mean_left: f64,
    /// <E[right(t)]>.
    pub mean_right: f64,
    /// <E[left(t)] * E[right(t)]>; preserves deterministic beats.
    pub mean_coherent_product: f64,
    /// <Cov(left(t), right(t))>, accumulated from random factors directly.
    pub mean_stochastic_product: f64,
    /// <E[left(t) * right(t)]>, including coherent and stochastic parts.
    pub mean_product: f64,
    /// Cov(<left(t)>, <right(t)>). This is the covariance of the window-averaged
    /// readings, not the window average of their instantaneous covariance.
    /// Each factor's complex window amplitude is summed before multiplication.
    pub averaged_covariance: f64,
}

#[derive(Clone, Debug)]
struct ReceiverTone {
    omega: f64,
    amplitudes: Vec<[Dd; 2]>,
}

#[derive(Clone, Debug)]
struct ReceiverBank {
    statistics: HarmonicFactorStatistics,
    tones: Vec<ReceiverTone>,
}

/// Immutable physical receiver projections of a coherent or statistical owner.
///
/// Preparation contracts each row in double-double arithmetic before any
/// temporal product. Distinct random factors remain separate. Acceleration is
/// obtained by multiplying the projected velocity by -i omega, avoiding an
/// unnecessary omega^2 intermediate. Original dynamic-solve residuals still
/// apply; projection cannot certify a relative error near a physical null.
///
/// Pair observation is allocation-free, with no solves, coordinate loops, or
/// matrix applications. It reuses the shared analytic uniform/triangular window
/// kernels. Only the admitted pair's two readings are contracted, so requesting
/// a small number of correlations does not build a dense receiver matrix.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicReceivers {
    receiver_count: usize,
    banks: Vec<ReceiverBank>,
    maximum_relative_response_residual: f64,
}

impl PreparedHarmonicInterference {
    /// Prepare real physical receiver rows on the complete coherent response.
    /// The returned owner is independent of the source and input row lifetimes.
    pub fn prepare_receivers(
        &self,
        receivers: &[HarmonicReceiver<'_>],
        controls: HarmonicReceiverControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedHarmonicReceivers, HarmonicResponseError> {
        PreparedHarmonicReceivers::prepare(
            self.n,
            receivers,
            core::iter::once((HarmonicFactorStatistics::Real, self.tones.as_slice())),
            self.maximum_relative_response_residual,
            controls,
            cx,
        )
    }
}

impl PreparedHarmonicReceivers {
    pub(super) fn prepare<'a>(
        n: usize,
        receivers: &[HarmonicReceiver<'_>],
        banks: impl Iterator<Item = (HarmonicFactorStatistics, &'a [Tone])> + Clone,
        maximum_relative_response_residual: f64,
        controls: HarmonicReceiverControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        if receivers.is_empty()
            || controls.maximum_receivers == 0
            || controls.maximum_receivers > MAX_INTERFERENCE_RECEIVERS
        {
            return Err(invalid(
                "nonempty receiver rows and a bounded receiver cap required",
            ));
        }
        cap("receiver rows", receivers.len(), controls.maximum_receivers)?;
        for receiver in receivers {
            poll(cx)?;
            if receiver.weights.len() != n || receiver.weights.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "finite receiver row in original coordinates required",
                ));
            }
        }
        let mut total_tones = 0usize;
        let mut product_terms = 0usize;
        let mut bank_count = 0usize;
        for (_, tones) in banks.clone() {
            poll(cx)?;
            bank_count = extent(bank_count.checked_add(1))?;
            total_tones = extent(total_tones.checked_add(tones.len()))?;
            product_terms = extent(
                tones
                    .len()
                    .checked_mul(tones.len())
                    .and_then(|x| product_terms.checked_add(x)),
            )?;
        }
        product_terms = extent(
            total_tones
                .checked_mul(4)
                .and_then(|x| x.checked_add(bank_count.checked_mul(2)?))
                .and_then(|x| x.checked_add(product_terms)),
        )?;
        let projection_terms = extent(
            n.checked_mul(receivers.len())
                .and_then(|x| x.checked_mul(total_tones)),
        )?;
        let payload = extent(
            total_tones
                .checked_mul(
                    receivers
                        .len()
                        .checked_mul(core::mem::size_of::<[Dd; 2]>())
                        .and_then(|x| x.checked_add(core::mem::size_of::<ReceiverTone>()))
                        .ok_or(invalid("receiver payload extent overflow"))?,
                )
                .and_then(|x| {
                    x.checked_add(bank_count.checked_mul(core::mem::size_of::<ReceiverBank>())?)
                })
                .and_then(|x| x.checked_add(core::mem::size_of::<Self>())),
        )?;
        cap(
            "receiver payload bytes",
            payload,
            controls.maximum_payload_bytes,
        )?;
        cap(
            "receiver projection terms",
            projection_terms,
            controls.maximum_projection_terms,
        )?;
        cap(
            "receiver product terms",
            product_terms,
            controls.maximum_product_terms,
        )?;

        let mut projected = Vec::with_capacity(bank_count);
        for (statistics, tones) in banks {
            poll(cx)?;
            let mut outputs = Vec::with_capacity(tones.len());
            for tone in tones {
                poll(cx)?;
                let mut amplitudes = Vec::with_capacity(receivers.len());
                for receiver in receivers {
                    let motion = match receiver.quantity {
                        HarmonicReceiverQuantity::Displacement => &tone.response.displacement,
                        HarmonicReceiverQuantity::Velocity
                        | HarmonicReceiverQuantity::Acceleration => &tone.response.velocity,
                    };
                    let mut amplitude = [Dd::ZERO; 2];
                    for (i, (&weight, &value)) in receiver.weights.iter().zip(motion).enumerate() {
                        if i % 64 == 0 {
                            poll(cx)?;
                        }
                        let weight = Dd::from_f64(weight);
                        let value = parts(value);
                        amplitude[0] = checked(amplitude[0] + weight * value[0])?;
                        amplitude[1] = checked(amplitude[1] + weight * value[1])?;
                    }
                    if receiver.quantity == HarmonicReceiverQuantity::Acceleration {
                        let frequency = Dd::from_f64(tone.omega);
                        amplitude = [
                            checked(amplitude[1] * frequency)?,
                            checked(-amplitude[0] * frequency)?,
                        ];
                    }
                    amplitudes.push(amplitude);
                }
                outputs.push(ReceiverTone {
                    omega: tone.omega,
                    amplitudes,
                });
            }
            projected.push(ReceiverBank {
                statistics,
                tones: outputs,
            });
        }
        poll(cx)?;
        Ok(Self {
            receiver_count: receivers.len(),
            banks: projected,
            maximum_relative_response_residual,
        })
    }

    /// Number of receiver rows, retaining their authored order.
    #[must_use]
    pub const fn receiver_count(&self) -> usize {
        self.receiver_count
    }

    /// The original owner's maximum dynamic-equation backward residual.
    /// This does not bound relative receiver error at destructive interference.
    #[must_use]
    pub const fn maximum_relative_response_residual(&self) -> f64 {
        self.maximum_relative_response_residual
    }

    /// Observe one pair of physical readings, or a mean square when left=right.
    ///
    /// The deterministic mean is bank zero. Every random bank contributes its
    /// declared real/proper second moment; no cross-factor coherence is inferred.
    /// Finite windows retain both real sum-frequency terms and slow beats.
    /// Cancellation/refusal publishes no partial result and mutates no state.
    pub fn observe_pair(
        &self,
        left: usize,
        right: usize,
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<HarmonicReceiverPairMoments, HarmonicResponseError> {
        poll(cx)?;
        if left >= self.receiver_count || right >= self.receiver_count {
            return Err(invalid("receiver index outside prepared rows"));
        }
        let mut means = [Dd::ZERO; 2];
        let mut coherent = Dd::ZERO;
        let mut stochastic = Dd::ZERO;
        let mut averaged_covariance = Dd::ZERO;
        for (index, bank) in self.banks.iter().enumerate() {
            poll(cx)?;
            let mut averaged = [[Dd::ZERO; 2]; 2];
            for tone in &bank.tones {
                poll(cx)?;
                // Re[z] and Re[-i*z] give both components without rounding the
                // analytic window transform to binary64 before cancellation.
                let kernel = window
                    .product_kernel(tone.omega, 0.0)
                    .map_err(observation)?;
                for (sum, receiver) in averaged.iter_mut().zip([left, right]) {
                    let amplitude = tone.amplitudes[receiver];
                    sum[0] = checked(
                        sum[0]
                            + kernel
                                .mean_real_product_parts_dd(amplitude, [Dd::ONE, Dd::ZERO])
                                .map_err(observation)?,
                    )?;
                    if index != 0 && bank.statistics == HarmonicFactorStatistics::ProperComplex {
                        sum[1] = checked(
                            sum[1]
                                + kernel
                                    .mean_real_product_parts_dd(
                                        [amplitude[1], -amplitude[0]],
                                        [Dd::ONE, Dd::ZERO],
                                    )
                                    .map_err(observation)?,
                        )?;
                    }
                }
            }
            if index == 0 {
                means = [averaged[0][0], averaged[1][0]];
            } else {
                let covariance = match bank.statistics {
                    HarmonicFactorStatistics::Real => checked(averaged[0][0] * averaged[1][0])?,
                    HarmonicFactorStatistics::ProperComplex => checked(
                        (averaged[0][0] * averaged[1][0] + averaged[0][1] * averaged[1][1])
                            * Dd::from_f64(0.5),
                    )?,
                };
                averaged_covariance = checked(averaged_covariance + covariance)?;
            }
            let mut moment = Dd::ZERO;
            for (i, a) in bank.tones.iter().enumerate() {
                for (j, b) in bank.tones.iter().enumerate() {
                    if left == right && j < i {
                        continue;
                    }
                    poll(cx)?;
                    let value = match bank.statistics {
                        HarmonicFactorStatistics::Real => window
                            .product_kernel(a.omega, b.omega)
                            .map_err(observation)?
                            .mean_real_product_parts_dd(a.amplitudes[left], b.amplitudes[right])
                            .map_err(observation)?,
                        HarmonicFactorStatistics::ProperComplex => window
                            .proper_product_kernel(a.omega, b.omega)
                            .map_err(observation)?
                            .mean_real_product_parts_dd(a.amplitudes[left], b.amplitudes[right])
                            .map_err(observation)?,
                    };
                    let multiplicity = if left == right && i != j { 2.0 } else { 1.0 };
                    moment = checked(moment + value * Dd::from_f64(multiplicity))?;
                }
            }
            if index == 0 {
                coherent = moment;
            } else {
                stochastic = checked(stochastic + moment)?;
            }
        }
        let total = checked(coherent + stochastic)?;
        poll(cx)?;
        Ok(HarmonicReceiverPairMoments {
            mean_left: means[0].to_f64(),
            mean_right: means[1].to_f64(),
            mean_coherent_product: coherent.to_f64(),
            mean_stochastic_product: stochastic.to_f64(),
            mean_product: total.to_f64(),
            averaged_covariance: averaged_covariance.to_f64(),
        })
    }
}

fn extent(value: Option<usize>) -> Result<usize, HarmonicResponseError> {
    value.ok_or(invalid("receiver extent overflow"))
}

#[cfg(test)]
mod tests;

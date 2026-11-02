//! Coherent delay-and-sum pressure channels from the original physical array.
//!
//! Every channel is y(t) = sum_r w_r p_r(t-delay_r), with authored real
//! dimensionless weights and nonnegative delays. Under exp(-i omega t), its
//! complex transfer is w_r exp(+i omega delay_r). Neither receiver powers nor
//! cross spectra replace this field sum. A channel is a processed signal,
//! not a new spatial receiver, pressure-to-intensity model or acoustic load.

use super::*;

/// Maximum processed pressure channels in one prepared bank.
pub const MAX_PRESSURE_ARRAY_CHANNELS: usize = 64;
/// Maximum authored taps over all channels, including silent/repeated taps.
pub const MAX_PRESSURE_ARRAY_TAPS: usize = 4096;
const MAX_CHANNEL_ID_BYTES: usize = 128;

/// One causal tap on an ORIGINAL physical receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPressureTap {
    pub receiver: usize,
    /// Real dimensionless gain. Signed and zero gains remain explicit.
    pub weight: f64,
    /// Actual delay, not an inferred propagation time or steering direction.
    pub delay_s: f64,
}

/// Authored pressure channel. Tap order is retained through compensated sums.
/// An empty tap list declares a dark channel. Gains are not normalized.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicPressureChannel {
    pub id: String,
    pub taps: Vec<HarmonicPressureTap>,
}

/// Whole preparation and repeated-window admission bounds.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPressureArrayLimits {
    pub maximum_channels: usize,
    pub maximum_taps: usize,
    /// Tap phasors, tap/field contractions and frequency/bank lookups. A lookup
    /// performs at most thirteen comparisons under the original 4096-line cap.
    pub maximum_prepare_products: usize,
    /// Linear real contractions, channel/frequency pairs and averaged-channel
    /// covariance terms. Conservative counts include both matrix triangles.
    pub maximum_query_products: usize,
    /// Logical retained payload plus preparation temporaries, query accumulators
    /// and one returned observation. Borrowed original fields and allocator
    /// overhead are excluded. Preparation and query temporaries do not coexist;
    /// counting both makes this a conservative peak bound.
    pub maximum_workspace_bytes: usize,
}

/// First moment and uncentered mean square of a processed pressure signal.
/// Mean square is E[average(y(t)^2)], not the variance of average(y).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicPressureChannelMoments {
    pub mean_pressure_pa: f64,
    pub mean_square_pressure_pa2: f64,
}

struct PressureLine {
    omega: f64,
    amplitudes: Vec<[Dd; 2]>,
}
struct PressureBank {
    lines: Vec<PressureLine>,
}

/// Prepared channel phasors, retaining original physical receiver provenance.
///
/// The optional deterministic mean and proper-complex factors come from the
/// original exposure owner. Each factor spans all receivers AND frequencies.
/// No independent-microphone, stationary or Gaussian approximation is made.
/// Queries retain carrier sum terms for the deterministic mean and the
/// difference-frequency coherence of proper factors.
///
/// Delay phases and all receiver/tap contractions are prepared once. Repeated
/// observations perform no BEM solve, delay trigonometry or spatial field
/// reconstruction. Original solver/source uncertainty is unchanged.
pub struct PreparedHarmonicPressureArray<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    channels: Vec<HarmonicPressureChannel>,
    // First bank is the deterministic mean, including when it is dark.
    banks: Vec<PressureBank>,
    workspace_bytes: usize,
    query_products: usize,
}

/// Processed pressure signals and their joint averaged-readout uncertainty.
pub struct WindowedHarmonicPressureArray<'a> {
    window: HarmonicWindow,
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    channels: &'a [HarmonicPressureChannel],
    moments: Vec<HarmonicPressureChannelMoments>,
    // Row-major, symmetric centered covariance of average(y_c).
    averaged_covariance_pa2: Vec<f64>,
}

impl WindowedHarmonicPressureArray<'_> {
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn physical_receivers(&self) -> &[ReceiverSpec] {
        &self.reference.receivers
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.reference.sources
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.reference.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        ProducerUncertainty::Unresolved {
            source_hash: self.reference.producer_hash,
        }
    }
    pub fn channels(&self) -> &[HarmonicPressureChannel] {
        self.channels
    }
    pub fn moments(&self) -> &[HarmonicPressureChannelMoments] {
        &self.moments
    }
    /// Row-major centered Cov(average(y_a), average(y_b)) [Pa²]. This is
    /// generally different from average of the instantaneous covariance.
    /// Deterministic means are never subtracted from large second moments:
    /// only integrated fluctuation factors contribute here.
    pub fn averaged_covariance_pa2(&self) -> &[f64] {
        &self.averaged_covariance_pa2
    }
}

impl<'a> PreparedHarmonicPressureArray<'a> {
    pub fn new(
        exposure: &PreparedHarmonicExposure<'a>,
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(exposure, channels, limits, || false)
    }

    /// Prepare authored causal taps on the same receiver/clock owner. Delays
    /// refer to the field's full authored harmonic continuation; this does not
    /// invent a finite pulse onset or a sampled fractional-delay filter.
    pub fn new_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        let n = channels.len();
        admit(n, MAX_PRESSURE_ARRAY_CHANNELS)?;
        admit(n, limits.maximum_channels)?;
        let mut taps = 0usize;
        let mut id_bytes = 0usize;
        for (i, channel) in channels.iter().enumerate() {
            poll(&mut cancelled)?;
            if channel.id.len() > MAX_CHANNEL_ID_BYTES
                || channel.id.trim().is_empty()
                || channels[..i].iter().any(|old| old.id == channel.id)
            {
                return Err(bad(
                    "pressure array channel identities must be unique and bounded",
                ));
            }
            taps = taps.checked_add(channel.taps.len()).ok_or_else(overflow)?;
            id_bytes = id_bytes
                .checked_add(channel.id.len())
                .ok_or_else(overflow)?;
            // Admit the cumulative count before scanning an unbounded tap list.
            admit(taps, MAX_PRESSURE_ARRAY_TAPS)?;
            admit(taps, limits.maximum_taps)?;
            for tap in &channel.taps {
                poll(&mut cancelled)?;
                if tap.receiver >= exposure.reference.receivers.len()
                    || !tap.weight.is_finite()
                    || !tap.delay_s.is_finite()
                    || tap.delay_s < 0.0
                {
                    return Err(bad(
                        "pressure array tap needs an original receiver, finite weight and nonnegative delay",
                    ));
                }
            }
        }
        let bank_count = exposure.factors.len() + 1;
        let mut lines = 0usize;
        let mut frequency_pairs = 0usize;
        for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
            lines = lines.checked_add(bank.lines.len()).ok_or_else(overflow)?;
            frequency_pairs = bank
                .lines
                .len()
                .checked_mul(bank.lines.len())
                .and_then(|p| frequency_pairs.checked_add(p))
                .ok_or_else(overflow)?;
        }
        // The number of unique frequencies is at most the retained line count.
        // Use this conservative admission before sorting or allocation.
        let prepare_products = lines
            .checked_mul(taps)
            .and_then(|p| p.checked_mul(2))
            .and_then(|p| p.checked_add(lines.checked_mul(bank_count)?))
            .ok_or_else(overflow)?;
        let query_products = lines
            .checked_mul(n)
            .and_then(|p| p.checked_mul(2))
            .and_then(|p| p.checked_add(frequency_pairs.checked_mul(n)?))
            .and_then(|p| p.checked_add(n.checked_mul(n)?.checked_mul(exposure.factors.len())?))
            .ok_or_else(overflow)?;
        admit(prepare_products, MAX_PRODUCTS)?;
        admit(prepare_products, limits.maximum_prepare_products)?;
        admit(query_products, MAX_PRODUCTS)?;
        admit(query_products, limits.maximum_query_products)?;

        let line_bytes = n
            .checked_mul(core::mem::size_of::<[Dd; 2]>())
            .and_then(|v| v.checked_add(core::mem::size_of::<PressureLine>()))
            .ok_or_else(overflow)?;
        let workspace_bytes = lines
            .checked_mul(line_bytes)
            .and_then(|v| {
                v.checked_add(bank_count.checked_mul(core::mem::size_of::<PressureBank>())?)
            })
            .and_then(|v| {
                v.checked_add(n.checked_mul(core::mem::size_of::<HarmonicPressureChannel>())?)
            })
            .and_then(|v| {
                v.checked_add(taps.checked_mul(core::mem::size_of::<HarmonicPressureTap>())?)
            })
            .and_then(|v| v.checked_add(id_bytes))
            .and_then(|v| v.checked_add(lines.checked_mul(core::mem::size_of::<f64>())?))
            .and_then(|v| v.checked_add(taps.checked_mul(core::mem::size_of::<[Dd; 2]>())?))
            .and_then(|v| {
                v.checked_add(n.checked_mul(
                    2 * core::mem::size_of::<[Dd; 2]>()
                        + core::mem::size_of::<HarmonicPressureChannelMoments>(),
                )?)
            })
            .and_then(|v| {
                v.checked_add(
                    n.checked_mul(n)?
                        .checked_mul(core::mem::size_of::<Dd>() + core::mem::size_of::<f64>())?,
                )
            })
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        let mut owned_channels = storage(n)?;
        for channel in channels {
            poll(&mut cancelled)?;
            let mut id = String::new();
            id.try_reserve_exact(channel.id.len())
                .map_err(|_| bad("pressure array identity allocation failed"))?;
            id.push_str(&channel.id);
            let mut owned_taps = storage(channel.taps.len())?;
            owned_taps.extend_from_slice(&channel.taps);
            owned_channels.push(HarmonicPressureChannel {
                id,
                taps: owned_taps,
            });
        }
        let mut frequencies = storage(lines)?;
        let mut banks = storage(bank_count)?;
        for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
            poll(&mut cancelled)?;
            let mut prepared = storage(bank.lines.len())?;
            for line in &bank.lines {
                poll(&mut cancelled)?;
                frequencies.push(line.omega);
                let mut amplitudes = storage(n)?;
                amplitudes.resize(n, [Dd::ZERO; 2]);
                prepared.push(PressureLine {
                    omega: line.omega,
                    amplitudes,
                });
            }
            banks.push(PressureBank { lines: prepared });
        }
        frequencies.sort_unstable_by(f64::total_cmp);
        frequencies.dedup_by(|a, b| a.to_bits() == b.to_bits());
        let mut weights = storage(taps)?;
        weights.resize(taps, [Dd::ZERO; 2]);
        for omega in frequencies {
            poll(&mut cancelled)?;
            let mut index = 0;
            for channel in channels {
                for &tap in &channel.taps {
                    poll(&mut cancelled)?;
                    weights[index] = delayed_weight(omega, tap.weight, tap.delay_s)?;
                    index += 1;
                }
            }
            for (bank_index, bank) in core::iter::once(&exposure.mean)
                .chain(exposure.factors.iter())
                .enumerate()
            {
                poll(&mut cancelled)?;
                let Ok(line_index) = bank
                    .lines
                    .binary_search_by(|line| line.omega.total_cmp(&omega))
                else {
                    continue;
                };
                let input = &bank.lines[line_index];
                let output = &mut banks[bank_index].lines[line_index];
                let mut offset = 0;
                for (channel_index, channel) in channels.iter().enumerate() {
                    poll(&mut cancelled)?;
                    let mut sum = [Dd::ZERO; 2];
                    for (j, tap) in channel.taps.iter().enumerate() {
                        if j % 64 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let a = component(&input.fields[tap.receiver], 0);
                        let [wr, wi] = weights[offset + j];
                        sum[0] = sum[0] + a[0] * wr - a[1] * wi;
                        sum[1] = sum[1] + a[0] * wi + a[1] * wr;
                    }
                    if sum.iter().any(|&v| !finite(v)) {
                        return Err(bad("pressure array coherent channel is unrepresentable"));
                    }
                    output.amplitudes[channel_index] = sum;
                    offset += channel.taps.len();
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            reference: exposure.reference,
            source_clock: exposure.source_clock,
            channels: owned_channels,
            banks,
            workspace_bytes,
            query_products,
        })
    }

    pub fn channels(&self) -> &[HarmonicPressureChannel] {
        &self.channels
    }
    pub fn physical_receivers(&self) -> &[ReceiverSpec] {
        &self.reference.receivers
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn query_products(&self) -> usize {
        self.query_products
    }

    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedHarmonicPressureArray<'_>, TransferError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Analytical moments on the output signal clock. Uniform and shaped
    /// detector windows preserve the original source/frequency coherence.
    /// No partial observation is returned after refusal or cancellation.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicPressureArray<'_>, TransferError> {
        poll(&mut cancelled)?;
        let n = self.channels.len();
        if n == 0 {
            return Ok(WindowedHarmonicPressureArray {
                window,
                reference: self.reference,
                source_clock: self.source_clock,
                channels: &self.channels,
                moments: Vec::new(),
                averaged_covariance_pa2: Vec::new(),
            });
        }
        let mut sums = storage(n)?;
        sums.resize(n, [Dd::ZERO; 2]); // Mean and mean square.
        let mut integrated = storage(n)?;
        integrated.resize(n, [Dd::ZERO; 2]);
        let mut covariance = storage(n * n)?;
        covariance.resize(n * n, Dd::ZERO);
        let zero = window
            .proper_product_kernel(0.0, 0.0)
            .map_err(window_error)?;
        for (bank_index, bank) in self.banks.iter().enumerate() {
            poll(&mut cancelled)?;
            integrated.fill([Dd::ZERO; 2]);
            for line in &bank.lines {
                poll(&mut cancelled)?;
                let dc = window
                    .product_kernel(line.omega, 0.0)
                    .map_err(window_error)?;
                for (sum, &a) in integrated.iter_mut().zip(&line.amplitudes) {
                    sum[0] = sum[0]
                        + dc.mean_real_product_parts_dd(a, [Dd::ONE, Dd::ZERO])
                            .map_err(window_error)?;
                    sum[1] = sum[1]
                        + dc.mean_real_product_parts_dd([a[1], -a[0]], [Dd::ONE, Dd::ZERO])
                            .map_err(window_error)?;
                }
            }
            if bank_index == 0 {
                for (sum, value) in sums.iter_mut().zip(&integrated) {
                    sum[0] = value[0];
                }
            } else {
                for row in 0..n {
                    poll(&mut cancelled)?;
                    for col in row..n {
                        covariance[row * n + col] = covariance[row * n + col]
                            + zero
                                .mean_real_product_parts_dd(integrated[row], integrated[col])
                                .map_err(window_error)?;
                    }
                }
            }
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled)?;
                    let kernel = if bank_index == 0 {
                        Kernel::Deterministic(
                            window
                                .product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        )
                    } else {
                        Kernel::Proper(
                            window
                                .proper_product_kernel(a.omega, b.omega)
                                .map_err(window_error)?,
                        )
                    };
                    let multiplicity = Dd::from_f64(if i == j { 1.0 } else { 2.0 });
                    for ((sum, &a), &b) in sums.iter_mut().zip(&a.amplitudes).zip(&b.amplitudes) {
                        sum[1] = sum[1] + multiplicity * kernel.product(a, b)?;
                    }
                }
            }
        }
        let mut moments = storage(n)?;
        for sum in sums {
            poll(&mut cancelled)?;
            let mean = value(sum[0])?;
            let square = value(sum[1])?;
            if square < 0.0 {
                return Err(bad("pressure array mean square became negative"));
            }
            moments.push(HarmonicPressureChannelMoments {
                mean_pressure_pa: mean,
                mean_square_pressure_pa2: square,
            });
        }
        let mut averaged_covariance_pa2 = storage(n * n)?;
        averaged_covariance_pa2.resize(n * n, 0.0);
        for row in 0..n {
            poll(&mut cancelled)?;
            for col in row..n {
                let result = value(covariance[row * n + col])?;
                if row == col && result < 0.0 {
                    return Err(bad("pressure array averaged variance became negative"));
                }
                averaged_covariance_pa2[row * n + col] = result;
                averaged_covariance_pa2[col * n + row] = result;
            }
        }
        poll(&mut cancelled)?;
        Ok(WindowedHarmonicPressureArray {
            window,
            reference: self.reference,
            source_clock: self.source_clock,
            channels: &self.channels,
            moments,
            averaged_covariance_pa2,
        })
    }
}

// Keep omega*delay and the near-unit cosine decrement through the entire tap
// contraction. Existing deterministic f64 trig supplies the leaves; DD does
// not certify relative error at an arbitrary beamforming null.
fn delayed_weight(omega: f64, weight: f64, delay: f64) -> Result<[Dd; 2], TransferError> {
    fn phase(x: f64) -> [Dd; 2] {
        let (s, c) = fs_math::det::sin_cos(x);
        let cosine = if c > 0.5 {
            let q = Dd::from_f64(fs_math::det::sin(0.5 * x));
            Dd::ONE - Dd::from_f64(2.0) * q * q
        } else if c < -0.5 {
            let q = Dd::from_f64(fs_math::det::cos(0.5 * x));
            -Dd::ONE + Dd::from_f64(2.0) * q * q
        } else {
            Dd::from_f64(c)
        };
        [cosine, Dd::from_f64(s)]
    }
    let angle = Dd::from_f64(omega) * Dd::from_f64(delay);
    if !finite(angle) {
        return Err(bad("pressure array delay phase is unrepresentable"));
    }
    let [ar, ai] = phase(angle.hi);
    let [br, bi] = phase(angle.lo);
    let weight = Dd::from_f64(weight);
    let result = [(ar * br - ai * bi) * weight, (ar * bi + ai * br) * weight];
    if result.iter().any(|&x| !finite(x)) {
        Err(bad("pressure array delayed gain is unrepresentable"))
    } else {
        Ok(result)
    }
}

fn value(v: Dd) -> Result<f64, TransferError> {
    let out = v.to_f64();
    if finite(v) && out.is_finite() {
        Ok(out)
    } else {
        Err(bad("pressure array moment is unrepresentable"))
    }
}

#[cfg(test)]
mod tests;

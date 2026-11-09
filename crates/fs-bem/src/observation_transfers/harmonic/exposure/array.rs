//! Coherent delay-and-sum pressure channels from the original physical array.
//!
//! Every channel is y(t) = sum_r w_r m_r(t-delay_r), with authored real
//! dimensionless weights and nonnegative delays. The default m_r is pressure;
//! optional ideal first-order or exact-frequency calibrated microphones project
//! actual pressure AND particle velocity before any channel sum or product.
//! Under exp(-i omega t), its
//! complex transfer is w_r exp(+i omega delay_r). Neither receiver powers nor
//! cross spectra replace this field sum. A channel is a processed signal,
//! not a new spatial receiver, pressure-to-intensity model or acoustic load.

use super::*;
use crate::near_field::FirstOrder;

mod band;

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
    /// Explicit microphone patterns add eight coefficient products per original
    /// receiver and four field projections per receiver and retained line.
    /// Calibrated microphones count sixteen real DD products per projection
    /// and complete frequency/response validation and conversion work.
    pub maximum_prepare_products: usize,
    /// Linear real contractions, temporal field products and averaged-channel
    /// covariance terms. Each bank admits the larger of its original pair work
    /// and a 17-mode projection bound; counts include both matrix triangles.
    pub maximum_query_products: usize,
    /// Logical retained payload plus preparation temporaries, query accumulators
    /// and one returned observation. Borrowed original fields and allocator
    /// overhead are excluded. Preparation and query temporaries do not coexist;
    /// counting both makes this a conservative peak bound. Fixed band and
    /// profile-preparation scratch is included when a bank can use projection.
    pub maximum_workspace_bytes: usize,
}

/// First moment and uncentered mean square of a processed pressure signal.
/// Mean square is E[average(y(t)^2)], not the variance of average(y).
/// Explicit first-order microphones report their calibrated Pa-equivalent
/// signal, including pressure/velocity interference, in these same units.
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
    statistics: HarmonicExposureFactorStatistics,
}

// Borrowed original fields for either an ensemble or one already fixed draw.
// The same microphone projection, delay phases, admission and channel sums
// serve both paths; no second directional or propagation model is introduced.
pub(super) struct PressureInput<'a, 'b> {
    pub(super) reference: &'a PreparedHarmonicObservation,
    pub(super) source_clock: &'a str,
    pub(super) mean: &'b [Line],
    pub(super) factors: &'b [Bank],
}

struct PressureBankInput<'a> {
    lines: &'a [Line],
    statistics: HarmonicExposureFactorStatistics,
}

impl<'a, 'b> PressureInput<'a, 'b> {
    fn from_exposure(exposure: &'b PreparedHarmonicExposure<'a>) -> Self {
        Self {
            reference: exposure.reference,
            source_clock: exposure.source_clock,
            mean: &exposure.mean.lines,
            factors: &exposure.factors,
        }
    }

    fn banks(&self) -> impl Iterator<Item = PressureBankInput<'_>> {
        core::iter::once(PressureBankInput {
            lines: self.mean,
            statistics: HarmonicExposureFactorStatistics::Real,
        })
        .chain(self.factors.iter().map(|bank| PressureBankInput {
            lines: &bank.lines,
            statistics: bank.statistics,
        }))
    }
}

mod signal;
pub use signal::{HarmonicPressureSampleLimits, PreparedHarmonicPressureSignal};

mod microphone;
pub use microphone::{HarmonicMicrophoneCalibrationLine, HarmonicMicrophoneResponse};

/// Prepared channel phasors, retaining original physical receiver provenance.
///
/// The deterministic mean and authored real/proper factors come from the
/// original exposure owner. Each factor spans all receivers AND frequencies;
/// its statistics survive the delay-and-sum transformation. No independent-
/// microphone, stationary or Gaussian approximation is made. Queries retain
/// carrier sum terms for the mean and real factors, and difference-frequency
/// coherence for every factor.
///
/// Delay phases and all receiver/tap contractions are prepared once. Repeated
/// observations perform no BEM solve, delay trigonometry or spatial field
/// reconstruction. Original solver/source uncertainty is unchanged.
pub struct PreparedHarmonicPressureArray<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    channels: Vec<HarmonicPressureChannel>,
    microphone_patterns: Option<Vec<FirstOrder>>,
    microphone_calibration: Option<Vec<HarmonicMicrophoneCalibrationLine>>,
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
    microphone_patterns: Option<&'a [FirstOrder]>,
    microphone_calibration: Option<&'a [HarmonicMicrophoneCalibrationLine]>,
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
    /// Ideal responses, when selected. Calibrated arrays instead expose
    /// microphone_calibration(); neither response is inferred from the field.
    pub fn microphone_patterns(&self) -> Option<&[FirstOrder]> {
        self.microphone_patterns
    }
    /// Complete authored complex response, at exact native angular frequencies.
    pub fn microphone_calibration(&self) -> Option<&[HarmonicMicrophoneCalibrationLine]> {
        self.microphone_calibration
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        Self::prepare(PressureInput::from_exposure(exposure), channels, None, None, limits, cancelled)
    }

    /// Prepare causal channels of the existing ideal first-order microphones.
    /// Patterns are in original physical receiver order, with their front axes
    /// pointing toward the front source. Each receiver is projected as
    /// alpha*p - (1-alpha)*rho*c*(front dot u) [Pa-equivalent], using actual
    /// retained particle velocity, including reactive near fields and standing
    /// waves. No plane-wave substitution u=p/(rho*c) is made.
    ///
    /// Projection occurs once per receiver/retained line before the signed
    /// delayed sums. Source and cross-frequency coherence, mixed real/proper
    /// statistics and the common clock survive. The original propagation
    /// medium supplies the same reference rho*c calibration as FirstOrder;
    /// no diaphragm response, electronic noise or receiver loading is inferred.
    pub fn new_with_microphones(
        exposure: &PreparedHarmonicExposure<'a>,
        patterns: &[FirstOrder],
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
    ) -> Result<Self, TransferError> {
        Self::new_with_microphones_cancellable(exposure, patterns, channels, limits, || false)
    }

    /// Cancellable counterpart. Patterns are owned by the prepared array;
    /// neither input patterns nor the source exposure need remain alive.
    /// Failure publishes no partial array and does not mutate either input.
    pub fn new_with_microphones_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        patterns: &[FirstOrder],
        channels: &[HarmonicPressureChannel],
        limits: HarmonicPressureArrayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        Self::prepare(PressureInput::from_exposure(exposure), channels, Some(patterns), None, limits, cancelled)
    }

    fn prepare(
        source: PressureInput<'a, '_>,
        channels: &[HarmonicPressureChannel],
        patterns: Option<&[FirstOrder]>,
        calibration: Option<&[HarmonicMicrophoneCalibrationLine]>,
        limits: HarmonicPressureArrayLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        if patterns.is_some() && calibration.is_some() {
            return Err(bad("choose ideal or calibrated microphone response"));
        }
        let microphones = if let Some(calibration) = calibration {
            microphone::admit_shape(
                calibration, source.reference.receivers.len(), &mut cancelled,
            )?;
            source.reference.receivers.len()
        } else {
            match patterns {
                Some(patterns) if patterns.len() == source.reference.receivers.len() => {
                    patterns.len()
                }
                Some(_) => {
                    return Err(bad(
                        "microphone patterns differ from original receiver order",
                    ));
                }
                None => 0,
            }
        };
        let n = channels.len();
        admit(n, MAX_PRESSURE_ARRAY_CHANNELS)?;
        admit(n, limits.maximum_channels)?;
        let mut taps = 0usize;
        let mut id_bytes = 0usize;
        let mut tapped_receivers = [false; crate::near_field::MAX_RECEIVERS];
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
                if tap.receiver >= source.reference.receivers.len()
                    || !tap.weight.is_finite()
                    || !tap.delay_s.is_finite()
                    || tap.delay_s < 0.0
                {
                    return Err(bad(
                        "pressure array tap needs an original receiver, finite weight and nonnegative delay",
                    ));
                }
                tapped_receivers[tap.receiver] = true;
            }
        }
        let bank_count = source.factors.len() + 1;
        let mut lines = 0usize;
        let mut temporal_products = 0usize;
        let mut band_scratch = 0usize;
        for bank in source.banks() {
            lines = lines.checked_add(bank.lines.len()).ok_or_else(overflow)?;
            let (products, scratch) = band::admission(bank.lines.len(), n, bank.statistics)?;
            temporal_products = temporal_products
                .checked_add(products)
                .ok_or_else(overflow)?;
            band_scratch = band_scratch.max(scratch);
        }
        // The number of unique frequencies is at most the retained line count.
        // Use this conservative admission before sorting or allocation.
        let microphone_products = if let Some(calibration) = calibration {
            microphone::extra_work(calibration.len(), microphones, lines)?
        } else {
            microphones.checked_mul(
                lines.checked_mul(4).and_then(|v| v.checked_add(8)).ok_or_else(overflow)?
            ).ok_or_else(overflow)?
        };
        let calibration_bytes = if let Some(calibration) = calibration {
            microphone::extra_bytes(calibration.len(), microphones)?
        } else {
            0
        };
        let prepare_products = lines
            .checked_mul(taps)
            .and_then(|p| p.checked_mul(2))
            .and_then(|p| p.checked_add(lines.checked_mul(bank_count)?))
            .and_then(|p| p.checked_add(microphone_products))
            .ok_or_else(overflow)?;
        let query_products = lines
            .checked_mul(n)
            .and_then(|p| p.checked_mul(2))
            .and_then(|p| p.checked_add(temporal_products))
            .and_then(|p| p.checked_add(n.checked_mul(n)?.checked_mul(source.factors.len())?))
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
            .and_then(|v| {
                v.checked_add(microphones.checked_mul(core::mem::size_of::<[Dd; 2]>())?)
            })
            .and_then(|v| {
                if patterns.is_some() {
                    v.checked_add(microphones.checked_mul(
                        core::mem::size_of::<FirstOrder>() + core::mem::size_of::<[Dd; 4]>(),
                    )?)
                } else {
                    Some(v)
                }
            })
            .and_then(|v| v.checked_add(calibration_bytes))
            .and_then(|v| v.checked_add(band_scratch))
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        let microphone_patterns = if let Some(patterns) = patterns {
            let mut owned = storage(microphones)?;
            owned.extend_from_slice(patterns);
            Some(owned)
        } else {
            None
        };
        let microphone_calibration = if let Some(calibration) = calibration {
            Some(microphone::copy_calibration(calibration, &mut cancelled)?)
        } else {
            None
        };
        let mut projections = storage(if patterns.is_some() { microphones } else { 0 })?;
        let mut calibrated_projections = storage(if calibration.is_some() { microphones } else { 0 })?;
        if calibration.is_some() {
            calibrated_projections.resize(microphones, [Dd::ZERO; 8]);
        }
        let mut projected = storage(microphones)?;
        projected.resize(microphones, [Dd::ZERO; 2]);
        if let Some(patterns) = patterns {
            let medium = source.reference.medium;
            let impedance = Dd::from_f64(medium.density) * Dd::from_f64(medium.sound_speed);
            for pattern in patterns {
                poll(&mut cancelled)?;
                let fraction = Dd::from_f64(pattern.pressure_fraction());
                let velocity_gain = -(Dd::ONE - fraction) * impedance;
                let front = pattern.front_axis();
                projections.push([
                    fraction,
                    velocity_gain * Dd::from_f64(front[0]),
                    velocity_gain * Dd::from_f64(front[1]),
                    velocity_gain * Dd::from_f64(front[2]),
                ]);
            }
        }
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
        for bank in source.banks() {
            poll(&mut cancelled)?;
            let mut prepared = storage(bank.lines.len())?;
            for line in bank.lines {
                poll(&mut cancelled)?;
                frequencies.push(line.omega);
                let mut amplitudes = storage(n)?;
                amplitudes.resize(n, [Dd::ZERO; 2]);
                prepared.push(PressureLine {
                    omega: line.omega,
                    amplitudes,
                });
            }
            banks.push(PressureBank {
                lines: prepared,
                statistics: bank.statistics,
            });
        }
        frequencies.sort_unstable_by(f64::total_cmp);
        frequencies.dedup_by(|a, b| a.to_bits() == b.to_bits());
        if let Some(calibration) = calibration {
            microphone::admit_frequencies(calibration, &frequencies, &mut cancelled)?;
        }
        let mut weights = storage(taps)?;
        weights.resize(taps, [Dd::ZERO; 2]);
        for (frequency, omega) in frequencies.into_iter().enumerate() {
            poll(&mut cancelled)?;
            if let Some(calibration) = calibration {
                microphone::coefficients(
                    &calibration[frequency], &mut calibrated_projections, &mut cancelled,
                )?;
            }
            let mut index = 0;
            for channel in channels {
                for &tap in &channel.taps {
                    poll(&mut cancelled)?;
                    weights[index] = delayed_weight(omega, tap.weight, tap.delay_s)?;
                    index += 1;
                }
            }
            for (bank_index, bank) in source.banks().enumerate()
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
                for (receiver, ((field, weights), signal)) in input
                    .fields
                    .iter()
                    .zip(&projections)
                    .zip(&mut projected)
                    .enumerate()
                {
                    poll(&mut cancelled)?;
                    if !tapped_receivers[receiver] {
                        continue;
                    }
                    *signal = [Dd::ZERO; 2];
                    for (c, &weight) in weights.iter().enumerate() {
                        let amplitude = component(field, c);
                        signal[0] = signal[0] + amplitude[0] * weight;
                        signal[1] = signal[1] + amplitude[1] * weight;
                    }
                    if signal.iter().any(|&v| !finite(v)) {
                        return Err(bad(
                            "pressure array microphone projection is unrepresentable",
                        ));
                    }
                }
                if calibration.is_some() {
                    microphone::project(
                        input, &calibrated_projections, &tapped_receivers,
                        &mut projected, &mut cancelled,
                    )?;
                }
                let mut offset = 0;
                for (channel_index, channel) in channels.iter().enumerate() {
                    poll(&mut cancelled)?;
                    let mut sum = [Dd::ZERO; 2];
                    for (j, tap) in channel.taps.iter().enumerate() {
                        if j % 64 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let a = if microphones == 0 {
                            component(&input.fields[tap.receiver], 0)
                        } else {
                            projected[tap.receiver]
                        };
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
            reference: source.reference,
            source_clock: source.source_clock,
            channels: owned_channels,
            microphone_patterns,
            microphone_calibration,
            banks,
            workspace_bytes,
            query_products,
        })
    }

    pub fn channels(&self) -> &[HarmonicPressureChannel] {
        &self.channels
    }
    /// Explicit ideal responses in original receiver order, when authored.
    /// Calibrated arrays expose microphone_calibration() instead.
    pub fn microphone_patterns(&self) -> Option<&[FirstOrder]> {
        self.microphone_patterns.as_deref()
    }
    /// Complete owned response data. None and no ideal patterns mean pressure.
    pub fn microphone_calibration(&self) -> Option<&[HarmonicMicrophoneCalibrationLine]> {
        self.microphone_calibration.as_deref()
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
    /// Original factor laws; the prepared channel bank owns these copies.
    pub fn factor_statistics(
        &self,
    ) -> impl ExactSizeIterator<Item = HarmonicExposureFactorStatistics> + '_ {
        self.banks.iter().skip(1).map(|bank| bank.statistics)
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
                microphone_patterns: self.microphone_patterns.as_deref(),
                microphone_calibration: self.microphone_calibration.as_deref(),
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
        for (bank_index, bank) in self.banks.iter().enumerate() {
            poll(&mut cancelled)?;
            integrated.fill([Dd::ZERO; 2]);
            let projected = band::accumulate(
                bank, window, &mut sums, &mut integrated, &mut cancelled,
            )?;
            // A proper temporal factor removes a common carrier only for
            // second moments. Its mode zero cannot replace the original-clock
            // field average used in the readout covariance.
            if !projected || bank.statistics == HarmonicExposureFactorStatistics::ProperComplex {
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
            }
            if bank_index == 0 {
                for (sum, value) in sums.iter_mut().zip(&integrated) {
                    sum[0] = value[0];
                }
            } else {
                let zero = bank.statistics.kernel(window, 0.0, 0.0)?;
                for row in 0..n {
                    poll(&mut cancelled)?;
                    for col in row..n {
                        covariance[row * n + col] = covariance[row * n + col]
                            + zero.product(integrated[row], integrated[col])?;
                    }
                }
            }
            if projected {
                continue;
            }
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled)?;
                    let kernel = bank.statistics.kernel(window, a.omega, b.omega)?;
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
            microphone_patterns: self.microphone_patterns.as_deref(),
            microphone_calibration: self.microphone_calibration.as_deref(),
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

mod correlation;
pub use correlation::{
    HarmonicPressureCorrelation, HarmonicPressureCorrelationLimits,
    HarmonicPressureCorrelationPair, MAX_PRESSURE_CORRELATION_PAIRS,
    PreparedHarmonicPressureCorrelations, WindowedHarmonicPressureCorrelations,
};

#[cfg(test)]
mod tests;

#[cfg(test)]
mod band_tests;

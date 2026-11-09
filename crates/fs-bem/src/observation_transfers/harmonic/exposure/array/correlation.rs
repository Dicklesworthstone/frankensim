//! Signed-lag, finite-window correlations of the retained physical array.
//!
//! R_ab(tau) = E[<y_a(t) y_b(t + tau)>]. Every requested pair retains
//! cross-frequency and cross-channel interference. Positive lag advances the
//! second signal; under exp(-i omega t) its multiplier is exp(-i omega tau).

use super::*;

mod band;

/// Ceiling on authored channel/lag requests, including repeats.
pub const MAX_PRESSURE_CORRELATION_PAIRS: usize = 65_536;

/// Ordered output channels and the signed lag of the second channel.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicPressureCorrelationPair {
    pub first: usize,
    pub second: usize,
    /// Evaluate the second channel at t + lag_s on the original source clock.
    /// This samples the authored harmonic continuation, not a pulse onset.
    pub lag_s: f64,
}

/// Whole preparation and repeated-window bounds. The borrowed array's
/// existing payload is excluded; no new propagation solve is performed.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPressureCorrelationLimits {
    pub maximum_channel_pairs: usize,
    /// One delay phasor and one complex contraction per pair/retained line.
    pub maximum_prepare_products: usize,
    /// Four linear real contractions per pair/line, one readout term per
    /// factor, and the larger of the original ordered-pair contraction and
    /// bounded field projection. Requested first channels share projected
    /// modes; signed-lag second channels retain their own field projections.
    pub maximum_query_products: usize,
    /// Logical retained pair/lag payload, bank headers, query accumulators
    /// and one returned observation. Eligible banks also reserve the fixed
    /// temporal factor, profile preparation and first-channel mode cache.
    /// Stack scalars and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Correlation and distinct source fluctuation statistics [Pa-equivalent].
/// Negative cross correlations are physical and retain their signs.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicPressureCorrelation {
    pub channels: HarmonicPressureCorrelationPair,
    pub mean_first_pressure_pa: f64,
    pub mean_lagged_second_pressure_pa: f64,
    /// E[<y_first(t) y_second(t + lag)>], including the deterministic field.
    pub mean_product_pa2: f64,
    /// <Cov(y_first(t), y_second(t + lag))>, from fluctuation factors only.
    /// This is ensemble centering at each time, not temporal detrending.
    pub mean_covariance_pa2: f64,
    /// Cov(<y_first(t)>, <y_second(t + lag)>). Averaging occurs BEFORE
    /// the product; this can vanish while mean_covariance_pa2 remains large.
    pub averaged_covariance_pa2: f64,
}

struct CorrelationBank {
    // Original line-major, requested pair-minor second-channel amplitudes.
    shifted: Vec<[Dd; 2]>,
}

/// Reusable signed-lag observation of an existing directional/delayed array.
///
/// All lag phases and complex channel shifts are prepared once. Repeated
/// windows project eligible coherent fields before their signed products;
/// wider bands use the original analytic real/proper pair kernels. Neither
/// path requires delay trigonometry, BEM solves, a source covariance matrix
/// or a time grid. Temporal factors and requested first-channel projections
/// are shared across pairs, without constructing a dense channel matrix.
pub struct PreparedHarmonicPressureCorrelations<'array, 'field> {
    array: &'array PreparedHarmonicPressureArray<'field>,
    pairs: Vec<HarmonicPressureCorrelationPair>,
    banks: Vec<CorrelationBank>,
    workspace_bytes: usize,
    query_products: usize,
}

/// Complete result in authored request order, preserving repeated pairs/lags.
pub struct WindowedHarmonicPressureCorrelations<'a> {
    window: HarmonicWindow,
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    channels: &'a [HarmonicPressureChannel],
    microphone_patterns: Option<&'a [FirstOrder]>,
    microphone_calibration: Option<&'a [HarmonicMicrophoneCalibrationLine]>,
    correlations: Vec<HarmonicPressureCorrelation>,
}

impl WindowedHarmonicPressureCorrelations<'_> {
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
    pub fn microphone_patterns(&self) -> Option<&[FirstOrder]> {
        self.microphone_patterns
    }
    /// Exact native-frequency response data retained by the original array.
    pub fn microphone_calibration(&self) -> Option<&[HarmonicMicrophoneCalibrationLine]> {
        self.microphone_calibration
    }
    pub fn correlations(&self) -> &[HarmonicPressureCorrelation] {
        &self.correlations
    }
}

impl<'array, 'field> PreparedHarmonicPressureCorrelations<'array, 'field> {
    pub fn new(
        array: &'array PreparedHarmonicPressureArray<'field>,
        pairs: &[HarmonicPressureCorrelationPair],
        limits: HarmonicPressureCorrelationLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(array, pairs, limits, || false)
    }

    /// Transactional preparation. Out-of-range channels, nonfinite lags,
    /// exceeded bounds and cancellation publish no partial observer.
    pub fn new_cancellable(
        array: &'array PreparedHarmonicPressureArray<'field>,
        pairs: &[HarmonicPressureCorrelationPair],
        limits: HarmonicPressureCorrelationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        admit(pairs.len(), MAX_PRESSURE_CORRELATION_PAIRS)?;
        admit(pairs.len(), limits.maximum_channel_pairs)?;
        let mut first_channels = [false; MAX_PRESSURE_ARRAY_CHANNELS];
        for pair in pairs {
            poll(&mut cancelled)?;
            if pair.first >= array.channels.len()
                || pair.second >= array.channels.len()
                || !pair.lag_s.is_finite()
            {
                return Err(bad(
                    "pressure correlation needs original channels and a finite signed lag",
                ));
            }
            first_channels[pair.first] = true;
        }
        let first_count = first_channels.iter().filter(|&&used| used).count();
        let mut lines = 0usize;
        let mut temporal_products = 0usize;
        let mut band_scratch = 0usize;
        for bank in &array.banks {
            lines = lines.checked_add(bank.lines.len()).ok_or_else(overflow)?;
            let (products, scratch) =
                band::admission(bank.lines.len(), pairs.len(), first_count, bank.statistics)?;
            temporal_products = temporal_products
                .checked_add(products)
                .ok_or_else(overflow)?;
            band_scratch = band_scratch.max(scratch);
        }
        let shifted_count = pairs.len().checked_mul(lines).ok_or_else(overflow)?;
        let prepare_products = shifted_count.checked_mul(2).ok_or_else(overflow)?;
        let query_products = lines
            .checked_mul(4)
            .and_then(|p| p.checked_add(array.banks.len() - 1))
            .and_then(|p| p.checked_mul(pairs.len()))
            .and_then(|p| p.checked_add(temporal_products))
            .ok_or_else(overflow)?;
        admit(prepare_products, MAX_PRODUCTS)?;
        admit(prepare_products, limits.maximum_prepare_products)?;
        admit(query_products, MAX_PRODUCTS)?;
        admit(query_products, limits.maximum_query_products)?;
        let workspace_bytes = shifted_count
            .checked_mul(core::mem::size_of::<[Dd; 2]>())
            .and_then(|v| {
                v.checked_add(
                    array
                        .banks
                        .len()
                        .checked_mul(core::mem::size_of::<CorrelationBank>())?,
                )
            })
            .and_then(|v| {
                v.checked_add(pairs.len().checked_mul(
                    core::mem::size_of::<HarmonicPressureCorrelationPair>()
                        + core::mem::size_of::<HarmonicPressureCorrelation>()
                        + core::mem::size_of::<[Dd; 5]>()
                        + core::mem::size_of::<[[Dd; 2]; 2]>(),
                )?)
            })
            .and_then(|v| {
                v.checked_add(
                    band_scratch.max(core::mem::size_of::<[bool; MAX_PRESSURE_ARRAY_CHANNELS]>()),
                )
            })
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        let mut owned_pairs = storage(pairs.len())?;
        owned_pairs.extend_from_slice(pairs);
        let mut banks = storage(array.banks.len())?;
        for bank in &array.banks {
            poll(&mut cancelled)?;
            let count = pairs
                .len()
                .checked_mul(bank.lines.len())
                .ok_or_else(overflow)?;
            let mut shifted = storage(count)?;
            for line in &bank.lines {
                poll(&mut cancelled)?;
                for (i, pair) in pairs.iter().enumerate() {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let [pr, pi] = delayed_weight(line.omega, 1.0, -pair.lag_s)?;
                    let [ar, ai] = line.amplitudes[pair.second];
                    let signal = [ar * pr - ai * pi, ar * pi + ai * pr];
                    if signal.iter().any(|&v| !finite(v)) {
                        return Err(bad("pressure correlation lagged signal is unrepresentable"));
                    }
                    shifted.push(signal);
                }
            }
            banks.push(CorrelationBank { shifted });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            array,
            pairs: owned_pairs,
            banks,
            workspace_bytes,
            query_products,
        })
    }

    pub fn pairs(&self) -> &[HarmonicPressureCorrelationPair] {
        &self.pairs
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
    ) -> Result<WindowedHarmonicPressureCorrelations<'_>, TransferError> {
        self.observe_window_cancellable(window, || false)
    }

    /// Integrate on the first channel's window. In a nonstationary field,
    /// swapping channels and negating the lag ALSO shifts the window by lag:
    /// R_ab(window, lag) = R_ba(window + lag, -lag). No stationary,
    /// independent-microphone, Gaussian or electronic-noise model is inferred.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicPressureCorrelations<'_>, TransferError> {
        poll(&mut cancelled)?;
        let n = self.pairs.len();
        if n == 0 {
            return Ok(WindowedHarmonicPressureCorrelations {
                window,
                reference: self.array.reference,
                source_clock: self.array.source_clock,
                channels: &self.array.channels,
                microphone_patterns: self.array.microphone_patterns.as_deref(),
                microphone_calibration: self.array.microphone_calibration.as_deref(),
                correlations: Vec::new(),
            });
        }
        let mut sums = storage(n)?;
        sums.resize(n, [Dd::ZERO; 5]); // two means, product, covariance, averaged covariance.
        let mut integrated = storage(n)?;
        integrated.resize(n, [[Dd::ZERO; 2]; 2]);
        for (bank_index, (bank, lagged)) in self.array.banks.iter().zip(&self.banks).enumerate() {
            poll(&mut cancelled)?;
            integrated.fill([[Dd::ZERO; 2]; 2]);
            let projected = band::accumulate(
                bank,
                lagged,
                &self.pairs,
                window,
                bank_index != 0,
                &mut sums,
                &mut integrated,
                &mut cancelled,
            )?;
            // A proper field factor removes only a common carrier from its
            // second moments. Its mode zero is not the absolute-clock mean
            // used for covariance of the two averaged readouts.
            if !projected || bank.statistics == HarmonicExposureFactorStatistics::ProperComplex {
                for (line_index, line) in bank.lines.iter().enumerate() {
                    poll(&mut cancelled)?;
                    let dc = window
                        .product_kernel(line.omega, 0.0)
                        .map_err(window_error)?;
                    for (i, (pair, sums)) in self.pairs.iter().zip(&mut integrated).enumerate() {
                        if i % 64 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let b = lagged.shifted[line_index * n + i];
                        for (sum, a) in sums.iter_mut().zip([line.amplitudes[pair.first], b]) {
                            sum[0] = sum[0]
                                + dc.mean_real_product_parts_dd(a, [Dd::ONE, Dd::ZERO])
                                    .map_err(window_error)?;
                            sum[1] = sum[1]
                                + dc.mean_real_product_parts_dd([a[1], -a[0]], [Dd::ONE, Dd::ZERO])
                                    .map_err(window_error)?;
                        }
                    }
                }
            }
            if bank_index == 0 {
                for (sum, a) in sums.iter_mut().zip(&integrated) {
                    sum[0] = a[0][0];
                    sum[1] = a[1][0];
                }
            } else {
                let zero = bank.statistics.kernel(window, 0.0, 0.0)?;
                for (i, (sum, a)) in sums.iter_mut().zip(&integrated).enumerate() {
                    if i % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    sum[4] = sum[4] + zero.product(a[0], a[1])?;
                }
            }
            if projected {
                continue;
            }
            for (i, a) in bank.lines.iter().enumerate() {
                for (j, b) in bank.lines[..=i].iter().enumerate() {
                    poll(&mut cancelled)?;
                    // Share each window kernel across ALL requested pairs.
                    let kernel = bank.statistics.kernel(window, a.omega, b.omega)?;
                    for (p, (pair, sum)) in self.pairs.iter().zip(&mut sums).enumerate() {
                        if p % 64 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let mut product =
                            kernel.product(a.amplitudes[pair.first], lagged.shifted[j * n + p])?;
                        if i != j {
                            // Off-diagonal contractions differ for ordered
                            // channels/lags. Doubling just one is incorrect.
                            product = product
                                + kernel
                                    .product(lagged.shifted[i * n + p], b.amplitudes[pair.first])?;
                        }
                        sum[2] = sum[2] + product;
                        if bank_index != 0 {
                            sum[3] = sum[3] + product;
                        }
                    }
                }
            }
        }
        let mut correlations = storage(n)?;
        for (&pair, sum) in self.pairs.iter().zip(sums) {
            poll(&mut cancelled)?;
            let result = HarmonicPressureCorrelation {
                channels: pair,
                mean_first_pressure_pa: value(sum[0])?,
                mean_lagged_second_pressure_pa: value(sum[1])?,
                mean_product_pa2: value(sum[2])?,
                mean_covariance_pa2: value(sum[3])?,
                averaged_covariance_pa2: value(sum[4])?,
            };
            if pair.first == pair.second
                && pair.lag_s == 0.0
                && (result.mean_product_pa2 < 0.0
                    || result.mean_covariance_pa2 < 0.0
                    || result.averaged_covariance_pa2 < 0.0)
            {
                return Err(bad(
                    "pressure correlation zero-lag variance became negative",
                ));
            }
            correlations.push(result);
        }
        poll(&mut cancelled)?;
        Ok(WindowedHarmonicPressureCorrelations {
            window,
            reference: self.array.reference,
            source_clock: self.array.source_clock,
            channels: &self.array.channels,
            microphone_patterns: self.array.microphone_patterns.as_deref(),
            microphone_calibration: self.array.microphone_calibration.as_deref(),
            correlations,
        })
    }
}

#[cfg(test)]
mod tests;

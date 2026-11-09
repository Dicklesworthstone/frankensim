//! Reusable local irreversible heat waveforms from the original film fields.
//!
//! Settled source harmonics retain their physical source clock. Ohmic and
//! Debye/Lorentz friction channels are sums of squared REAL reservoir fields; their
//! sum/difference tones are prepared by the shared analytic spectrum owner.
//! Reversible electromagnetic/polarization storage is deliberately separate.
//! Admitted uniform and triangular shutters retain the original positive reservoir
//! factors instead of recovering a weak heat residual from bright signed tones.

use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
    HarmonicProductTone, MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum,
    ProductSpectrumError, phase_rad_s_dd_from_clock_dd, phase_rad_s_from_clock_dd,
};

/// Original irreversible reservoir(s) contributing to a local heat load.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngularSpectrumExposureHeatChannel {
    /// sigma |E_real|² from the original conductor.
    Ohmic,
    /// One original Debye pole, using its original material ordering.
    DebyePole { index: usize },
    /// One original Lorentz pole, using its original material ordering.
    LorentzPole { index: usize },
    /// Sum of all original polarization friction reservoirs.
    Polarization,
    /// Ohmic plus all polarization friction heat, with no reversible storage rate.
    Total,
}

/// Additional admission for one grid-point heat spectrum.
///
/// Products count three Cartesian contractions for every ordered frequency
/// pair, selected reservoir and source bank. Generated tones count twice that
/// for deterministic/real banks and once for proper-complex banks.
/// Workspace conservatively includes the prepared layer's original complete
/// workspace, one retained point's reservoir fields and the unsimplified
/// shared spectrum. Borrowed source owners and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureHeatSpectrumLimits {
    pub maximum_products: usize,
    pub maximum_generated_terms: usize,
    pub maximum_workspace_bytes: usize,
}

/// Expected irreversible heat density [W/m³] at one original film grid point.
///
/// The original prepared layer, selected reservoir, transverse sample and
/// depth remain attached. This is a local density, not absorbed film watts.
/// A physical volume/quadrature is required before supplying a lumped thermal
/// power. Do not also debit the corresponding integrated film absorption.
///
/// Positive-time EM phasors are conjugated into the shared negative-time
/// spectrum convention. Angular frequencies retain their original clock
/// without a round trip through f64 Hz. Proper factors remove sum tones only
/// by their authored second moments; real factors retain phase-locked noise.
/// No switch-on state, material-temperature feedback, source calibration,
/// stationary spectrum or relative numerical certificate at a null is added.
/// The prepared owner retains its original DD reservoir fields and source-bank
/// boundaries within the already admitted preparation workspace. Instantaneous
/// queries sum each coherent reservoir field before its real/proper square;
/// admitted uniform and triangular observations reuse its positive temporal factor.
/// Neither query allocates, reconstructs a field FFT or solves a material.
/// The exported spectrum retains the original signed sum/difference tones.
pub struct AngularSpectrumExposureLayerHeatSpectrum<'p, 'e, 'b, 's, 'm> {
    layer: &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>,
    depth_m: f64,
    point_index: usize,
    channel: AngularSpectrumExposureHeatChannel,
    spectrum: PreparedHarmonicProductSpectrum,
    heat_banks: Vec<HeatBank>,
    heat_lines: Vec<HeatLine>,
    reservoirs_per_source: usize,
    workspace_bytes: usize,
}

impl<'p, 'e, 'b, 's, 'm> AngularSpectrumExposureLayerHeatSpectrum<'p, 'e, 'b, 's, 'm> {
    pub fn layer(&self) -> &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
        self.layer
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.layer.anchor.grid()
    }
    /// Original row-major sample, x varying fastest.
    pub fn point_index(&self) -> usize {
        self.point_index
    }
    /// Depth from the selected film's first face, in metres.
    pub fn depth_m(&self) -> f64 {
        self.depth_m
    }
    pub fn channel(&self) -> AngularSpectrumExposureHeatChannel {
        self.channel
    }
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    /// Conservative peak logical storage admitted during preparation.
    pub fn preparation_workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Actual or ensemble-expected local irreversible heating [W/m³].
    pub fn evaluate(&self, time_s: f64) -> Result<f64, AngularSpectrumError> {
        self.evaluate_cancellable(time_s, || false)
    }
    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        positive_heat_at_time(&self.heat_banks, &self.heat_lines, time_s, &mut cancelled)
    }
    /// Detector-weighted expected heat density [W/m³] for either profile.
    /// Each eligible source bank resolves coherent nulls before squaring its
    /// reservoirs. Other banks retain their original real/proper DD pair law;
    /// if none admit, the cached flattened spectrum remains the fast path.
    pub fn observe_window(&self, window: HarmonicWindow) -> Result<f64, AngularSpectrumError> {
        self.observe_window_cancellable(window, || false)
    }
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        if let Some(heat) = positive_heat_mean(
            &self.heat_banks,
            &self.heat_lines,
            self.reservoirs_per_source,
            window,
            &mut cancelled,
        )? {
            return Ok(heat);
        }
        nonnegative_heat(
            self.spectrum
                .mean_cancellable(window, cancelled)
                .map_err(spectrum_error)?,
        )
    }
}

#[derive(Clone, Copy)]
struct PointElectric {
    omega: f64,
    cache: usize,
    electric: [C64; 3],
}

#[derive(Clone, Copy)]
struct HeatLine {
    omega: f64,
    field: [[Dd; 2]; 3],
}

struct HeatBank {
    start: usize,
    end: usize,
    law: HarmonicProductLaw,
    left_scale: Dd,
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    /// Prepare one reusable local heat spectrum without a new film/material
    /// solve. Each original frequency plane is reconstructed once and only
    /// this point is retained. Reservoir filters are applied once per line,
    /// then lazy scalar products feed the shared preparation engine.
    ///
    /// For G distinct bank-frequency groups, selected R reservoirs and K
    /// generated terms, additional retained scratch is O(G*R + K). Work is
    /// original field reconstruction plus O(G*R + K log K). Original field
    /// work limits remain enforced by this immutable prepared owner.
    pub fn heat_spectrum<'p>(
        &'p self,
        depth_m: f64,
        point_index: usize,
        channel: AngularSpectrumExposureHeatChannel,
        limits: AngularSpectrumExposureHeatSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerHeatSpectrum<'p, 'e, 'b, 's, 'm>, AngularSpectrumError>
    {
        poll(&mut cancelled)?;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.thickness_m {
            return Err(bad("heat spectrum depth is outside the selected layer"));
        }
        if point_index >= self.anchor.modes.len() {
            return Err(bad("heat spectrum point is outside the original grid"));
        }
        let (first_reservoir, reservoirs) = reservoir_range(
            channel,
            self.law.conductivity,
            self.law.debye_poles.len(),
            self.law.lorentz_poles.len(),
        )?;
        let mut products = 0usize;
        let mut generated = 0usize;
        let mut line_count = 0usize;
        let mut bank_count = 0usize;
        let mut maximum_frequencies = 0usize;
        for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
            poll(&mut cancelled)?;
            let f = bank.frequency_count();
            if f == 0 || reservoirs == 0 {
                continue;
            }
            maximum_frequencies = maximum_frequencies.max(f);
            line_count = f
                .checked_mul(reservoirs)
                .and_then(|v| line_count.checked_add(v))
                .ok_or_else(|| bad("heat spectrum line count overflows"))?;
            bank_count = bank_count
                .checked_add(reservoirs)
                .ok_or_else(|| bad("heat spectrum bank count overflows"))?;
            let pairs = f
                .checked_mul(f)
                .and_then(|v| v.checked_mul(reservoirs))
                .and_then(|v| v.checked_mul(3))
                .ok_or_else(|| bad("heat spectrum product count overflows"))?;
            products = products
                .checked_add(pairs)
                .ok_or_else(|| bad("heat spectrum product count overflows"))?;
            generated = pairs
                .checked_mul(if proper { 1 } else { 2 })
                .and_then(|v| generated.checked_add(v))
                .ok_or_else(|| bad("heat spectrum tone count overflows"))?;
        }
        let tone_bytes = generated
            .checked_mul(core::mem::size_of::<HarmonicProductTone>())
            .ok_or_else(|| bad("heat spectrum storage overflows"))?;
        let workspace_bytes = line_count
            .checked_mul(core::mem::size_of::<HeatLine>())
            .and_then(|v| v.checked_add(bank_count.checked_mul(core::mem::size_of::<HeatBank>())?))
            .and_then(|v| {
                v.checked_add(
                    maximum_frequencies.checked_mul(core::mem::size_of::<PointElectric>())?,
                )
            })
            .and_then(|v| v.checked_add(tone_bytes))
            .and_then(|v| v.checked_add(self.workspace_bytes))
            .ok_or_else(|| bad("heat spectrum storage overflows"))?;
        if products > limits.maximum_products
            || generated > limits.maximum_generated_terms
            || generated > MAX_PRODUCT_SPECTRUM_TERMS
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad("heat spectrum exceeds its work or workspace budget"));
        }
        let mut point_lines = Vec::new();
        point_lines
            .try_reserve_exact(maximum_frequencies)
            .map_err(|_| bad("heat spectrum point allocation failed"))?;
        let mut heat_lines = Vec::new();
        heat_lines
            .try_reserve_exact(line_count)
            .map_err(|_| bad("heat spectrum reservoir allocation failed"))?;
        let mut heat_banks = Vec::new();
        heat_banks
            .try_reserve_exact(bank_count)
            .map_err(|_| bad("heat spectrum bank allocation failed"))?;
        if reservoirs != 0 {
            for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
                poll(&mut cancelled)?;
                if bank.frequency_count() == 0 {
                    continue;
                }
                point_lines.clear();
                for group in 0..bank.frequency_count() {
                    poll(&mut cancelled)?;
                    let cache = self
                        .frequencies
                        .binary_search_by(|value| value.omega.total_cmp(&bank.frequencies[group]))
                        .map_err(|_| bad("original film frequency cache is missing"))?;
                    let layer = &self.frequencies[cache].layer;
                    let fields =
                        bank.frequency_field_with(group, &mut cancelled, |line, cancelled| {
                            line.modal_layer_field(
                                layer.as_ref(),
                                self.layer_index,
                                depth_m,
                                cancelled,
                            )
                        })?;
                    let field = fields[point_index];
                    point_lines.push(PointElectric {
                        omega: bank.frequencies[group],
                        cache,
                        electric: [field[0], field[1], field[2]],
                    });
                }
                for reservoir in first_reservoir..first_reservoir + reservoirs {
                    poll(&mut cancelled)?;
                    let start = heat_lines.len();
                    let root_sigma = if reservoir == 0 {
                        det::sqrt(self.law.conductivity)
                    } else {
                        1.0
                    };
                    let root = Dd::from_f64(root_sigma);
                    // Balanced square-root scaling avoids amplifying one
                    // operand alone; this correction retains authored sigma
                    // instead of replacing it by a rounded sqrt(sigma)^2.
                    let left_scale = if reservoir == 0 {
                        Dd::from_f64(self.law.conductivity) / (root * root)
                    } else {
                        Dd::ONE
                    };
                    for (index, point) in point_lines.iter().enumerate() {
                        if index % 64 == 0 {
                            poll(&mut cancelled)?;
                        }
                        let transfer = if reservoir == 0 {
                            C64::from_re(root_sigma)
                        } else {
                            self.frequencies[point.cache].poles[reservoir - 1].heat
                        };
                        heat_lines.push(HeatLine {
                            omega: point.omega,
                            field: point
                                .electric
                                .map(|e| filtered_negative_time_parts(e, transfer)),
                        });
                    }
                    heat_banks.push(HeatBank {
                        start,
                        end: heat_lines.len(),
                        law: if proper {
                            HarmonicProductLaw::ProperComplex
                        } else {
                            HarmonicProductLaw::Real
                        },
                        left_scale,
                    });
                }
            }
        }
        let spectrum = prepare_products(
            &heat_banks,
            &heat_lines,
            HarmonicProductSpectrumLimits {
                maximum_products: products,
                maximum_generated_terms: generated,
                maximum_workspace_bytes: tone_bytes,
            },
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureLayerHeatSpectrum {
            layer: self,
            depth_m,
            point_index,
            channel,
            spectrum,
            heat_banks,
            heat_lines,
            reservoirs_per_source: reservoirs,
            workspace_bytes,
        })
    }
}

// Reservoir zero is the original conductor; subsequent indices are its
// original Debye poles. A zero-conductivity Ohmic channel is exactly dark.
pub(super) fn reservoir_range(
    channel: AngularSpectrumExposureHeatChannel,
    sigma: f64,
    debye_poles: usize,
    lorentz_poles: usize,
) -> Result<(usize, usize), AngularSpectrumError> {
    use AngularSpectrumExposureHeatChannel as Channel;
    let poles = debye_poles + lorentz_poles;
    match channel {
        Channel::Ohmic => Ok((0, usize::from(sigma > 0.0))),
        Channel::DebyePole { index } if index < debye_poles => Ok((index + 1, 1)),
        Channel::DebyePole { .. } => Err(bad(
            "heat spectrum Debye pole is outside the original material",
        )),
        Channel::LorentzPole { index } if index < lorentz_poles => Ok((debye_poles + index + 1, 1)),
        Channel::LorentzPole { .. } => Err(bad(
            "heat spectrum Lorentz pole is outside the original material",
        )),
        Channel::Polarization => Ok((1, poles)),
        Channel::Total => Ok((usize::from(sigma == 0.0), poles + usize::from(sigma > 0.0))),
    }
}

// Filter the complete original electric field first. Retain DD products and
// conjugate the positive-time EM response into the shared negative-time law.
fn filtered_negative_time_parts(e: C64, transfer: C64) -> [Dd; 2] {
    let (a, b) = (Dd::from_f64(e.re), Dd::from_f64(e.im));
    let (c, d) = (Dd::from_f64(transfer.re), Dd::from_f64(transfer.im));
    [a * c - b * d, -(a * d + b * c)]
}

// An instantaneous reservoir is a coherent vector field, not a cancellation
// of pre-squared tones. Factor out one original source clock before summation;
// a proper factor needs only relative phases, while a real field restores its
// absolute carrier once. Each independent source/reservoir is squared separately.
// Work is linear in the already retained reservoir lines, with fixed scratch.
fn positive_heat_at_time(
    banks: &[HeatBank],
    fields: &[HeatLine],
    time_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AngularSpectrumError> {
    poll(cancelled)?;
    if !time_s.is_finite() {
        return Err(bad("retained heat time is nonfinite"));
    }
    let time = Dd::from_f64(time_s);
    let mut heat = Dd::ZERO;
    for bank in banks {
        poll(cancelled)?;
        let lines = &fields[bank.start..bank.end];
        let Some(first) = lines.first() else {
            continue;
        };
        let reference = Dd::from_f64(first.omega);
        let mut field = [[Dd::ZERO; 2]; 3];
        for line in lines {
            poll(cancelled)?;
            let phase = phase_rad_s_dd_from_clock_dd(Dd::from_f64(line.omega) - reference, time)
                .map_err(|_| bad("retained heat relative phase is unrepresentable"))?;
            for (sum, [real, imaginary]) in field.iter_mut().zip(line.field) {
                sum[0] = sum[0] + real * phase[0] - imaginary * phase[1];
                sum[1] = sum[1] + real * phase[1] + imaginary * phase[0];
            }
        }
        match bank.law {
            HarmonicProductLaw::Real => {
                let phase = phase_rad_s_from_clock_dd(first.omega, time)
                    .map_err(|_| bad("retained heat carrier phase is unrepresentable"))?;
                for [real, imaginary] in field {
                    let value = real * phase[0] - imaginary * phase[1];
                    heat = heat + (value * bank.left_scale) * value;
                }
            }
            HarmonicProductLaw::ProperComplex => {
                // Apply the peak-amplitude and authored material normalization
                // before the square: the physical half-power can be finite
                // even when the unscaled analytic norm overflows.
                let scale = Dd::from_f64(0.5) * bank.left_scale;
                for [real, imaginary] in field {
                    heat = heat + (real * scale) * real + (imaginary * scale) * imaginary;
                }
            }
        }
    }
    let heat = heat.to_f64();
    if !heat.is_finite() {
        return Err(bad("retained instantaneous heat is unrepresentable"));
    }
    poll(cancelled)?;
    nonnegative_heat(heat)
}

// Preparation already admits and constructs these buffers. Retaining them
// preserves the positive reservoir factorization without another source/FFT
// query, a larger preparation peak, or an observation allocation. Consecutive
// reservoir groups share one original source bank and one temporal factor.
fn positive_heat_mean(
    banks: &[HeatBank],
    fields: &[HeatLine],
    reservoirs_per_source: usize,
    window: HarmonicWindow,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, AngularSpectrumError> {
    poll(cancelled)?;
    if banks.is_empty() {
        return Ok(Some(0.0));
    }
    if reservoirs_per_source == 0 || banks.len() % reservoirs_per_source != 0 {
        return Err(bad("retained heat reservoir groups are inconsistent"));
    }
    let mut frequencies = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    // Keep the cached spectrum fast path when no source benefits. An unrelated
    // narrow, wide-carrier or even silent bank must not disable an eligible
    // source's positive contraction.
    let mut any_positive = false;
    for source in banks.chunks_exact(reservoirs_per_source) {
        poll(cancelled)?;
        let first = &source[0];
        let lines = &fields[first.start..first.end];
        let proper = first.law == HarmonicProductLaw::ProperComplex;
        if positive_band_tones(lines.len(), proper).is_none() {
            continue;
        }
        for (frequency, line) in frequencies.iter_mut().zip(lines) {
            *frequency = line.omega;
        }
        let admitted = if proper {
            window.band_mode_count_rad_s(&frequencies[..lines.len()])
        } else {
            window.real_band_mode_count_rad_s(&frequencies[..lines.len()])
        }
        .map_err(|_| bad("retained heat temporal factor is unrepresentable"))?;
        if admitted.is_some() {
            any_positive = true;
            break;
        }
    }
    if !any_positive {
        return Ok(None);
    }
    let mut workspace = HarmonicBandWorkspace::new();
    let mut heat = Dd::ZERO;
    for source in banks.chunks_exact(reservoirs_per_source) {
        poll(cancelled)?;
        let first = &source[0];
        let lines = &fields[first.start..first.end];
        let proper = first.law == HarmonicProductLaw::ProperComplex;
        let positive = if positive_band_tones(lines.len(), proper).is_some() {
            for (frequency, line) in frequencies.iter_mut().zip(lines) {
                *frequency = line.omega;
            }
            PositiveLayerMoments::prepare(
                window,
                &frequencies[..lines.len()],
                proper,
                &mut workspace,
                cancelled,
            )?
        } else {
            None
        };
        let Some(positive) = positive else {
            heat = heat
                + if proper {
                    retained_pair_heat_mean(
                        source,
                        fields,
                        |a, b| window.proper_product_kernel(a, b),
                        ProperHarmonicProduct::mean_real_product_parts_dd,
                        cancelled,
                    )?
                } else {
                    retained_pair_heat_mean(
                        source,
                        fields,
                        |a, b| window.product_kernel(a, b),
                        HarmonicProduct::mean_real_product_parts_dd,
                        cancelled,
                    )?
                };
            continue;
        };
        for bank in source {
            let lines = &fields[bank.start..bank.end];
            for axis in 0..3 {
                heat = heat
                    + positive
                        .square_scaled(
                            |i| {
                                // The retained filter already uses negative
                                // time. Restore the layer helper's native EM
                                // input convention without rounding a field.
                                let [real, imaginary] = lines[i].field[axis];
                                [real, -imaginary]
                            },
                            bank.left_scale,
                            cancelled,
                        )
                        .map_err(moment_error)?;
            }
        }
    }
    let heat = heat.to_f64();
    if !heat.is_finite() {
        return Err(bad("retained heat moment is unrepresentable"));
    }
    poll(cancelled)?;
    nonnegative_heat(heat).map(Some)
}

// A mixed observation still honors unsupported banks' exact source statistics
// and shutter. Each kernel serves all reservoirs of that source; work is at
// most its already admitted ordered products, with no temporary spectrum.
fn retained_pair_heat_mean<P: Copy>(
    banks: &[HeatBank],
    fields: &[HeatLine],
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    product: impl Fn(P, [Dd; 2], [Dd; 2]) -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, AngularSpectrumError> {
    let first = &banks[0];
    let lines = &fields[first.start..first.end];
    let mut heat = Dd::ZERO;
    for (i, left) in lines.iter().enumerate() {
        for (j, right) in lines.iter().enumerate() {
            poll(cancelled)?;
            let kernel = kernel(left.omega, right.omega)
                .map_err(|_| bad("retained heat pair kernel is unrepresentable"))?;
            for bank in banks {
                for axis in 0..3 {
                    let a = fields[bank.start + i].field[axis].map(|value| value * bank.left_scale);
                    let b = fields[bank.start + j].field[axis];
                    heat = heat
                        + product(kernel, a, b)
                            .map_err(|_| bad("retained heat pair product is unrepresentable"))?;
                }
            }
        }
    }
    if !heat.to_f64().is_finite() {
        return Err(bad("retained heat pair moment is unrepresentable"));
    }
    Ok(heat)
}

fn prepare_products(
    banks: &[HeatBank],
    fields: &[HeatLine],
    limits: HarmonicProductSpectrumLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedHarmonicProductSpectrum, AngularSpectrumError> {
    let inputs = banks.iter().flat_map(|bank| {
        let lines = &fields[bank.start..bank.end];
        lines.iter().flat_map(move |a| {
            lines.iter().flat_map(move |b| {
                (0..3).map(move |axis| AngularHarmonicProductInput {
                    left_angular_frequency_rad_s: a.omega,
                    right_angular_frequency_rad_s: b.omega,
                    left_amplitude: a.field[axis].map(|value| value * bank.left_scale),
                    right_amplitude: b.field[axis],
                    law: bank.law,
                })
            })
        })
    });
    PreparedHarmonicProductSpectrum::from_angular_frequencies_cancellable(inputs, limits, cancelled)
        .map_err(spectrum_error)
}

pub(super) fn nonnegative_heat(value: f64) -> Result<f64, AngularSpectrumError> {
    if value < 0.0 {
        Err(bad(
            "irreversible heat spectrum lost numerical nonnegativity",
        ))
    } else {
        Ok(value)
    }
}

pub(super) fn spectrum_error(error: ProductSpectrumError) -> AngularSpectrumError {
    match error {
        ProductSpectrumError::Cancelled => AngularSpectrumError::Cancelled,
        _ => bad("heat spectrum is outside its work/storage or numerical domain"),
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod positive_window_tests;

#[cfg(test)]
mod instantaneous_tests;

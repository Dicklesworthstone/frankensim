//! Reusable local irreversible heat waveforms from the original film fields.
//!
//! Settled source harmonics retain their physical source clock. Ohmic and
//! Debye friction channels are sums of squared REAL reservoir fields; their
//! sum/difference tones are prepared by the shared analytic spectrum owner.
//! Reversible electromagnetic/polarization storage is deliberately separate.

use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
    HarmonicProductTone, MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum,
    ProductSpectrumError,
};

/// Original irreversible reservoir(s) contributing to a local heat load.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngularSpectrumExposureHeatChannel {
    /// sigma |E_real|² from the original conductor.
    Ohmic,
    /// One original Debye pole, using its original material ordering.
    DebyePole { index: usize },
    /// Sum of all original Debye friction reservoirs.
    Polarization,
    /// Ohmic plus all Debye friction heat, with no reversible storage rate.
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
pub struct AngularSpectrumExposureLayerHeatSpectrum<'p, 'e, 'b, 's, 'm> {
    layer: &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>,
    depth_m: f64,
    point_index: usize,
    channel: AngularSpectrumExposureHeatChannel,
    spectrum: PreparedHarmonicProductSpectrum,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        nonnegative_heat(
            self.spectrum
                .evaluate_cancellable(time_s, cancelled)
                .map_err(spectrum_error)?,
        )
    }
    /// Detector-weighted expected heat density [W/m³] for either profile.
    pub fn observe_window(&self, window: HarmonicWindow) -> Result<f64, AngularSpectrumError> {
        self.observe_window_cancellable(window, || false)
    }
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
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
        let (first_reservoir, reservoirs) =
            reservoir_range(channel, self.law.conductivity, self.law.poles.len())?;
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
            workspace_bytes,
        })
    }
}

// Reservoir zero is the original conductor; subsequent indices are its
// original Debye poles. A zero-conductivity Ohmic channel is exactly dark.
pub(super) fn reservoir_range(
    channel: AngularSpectrumExposureHeatChannel,
    sigma: f64,
    poles: usize,
) -> Result<(usize, usize), AngularSpectrumError> {
    use AngularSpectrumExposureHeatChannel as Channel;
    match channel {
        Channel::Ohmic => Ok((0, usize::from(sigma > 0.0))),
        Channel::DebyePole { index } if index < poles => Ok((index + 1, 1)),
        Channel::DebyePole { .. } => Err(bad(
            "heat spectrum Debye pole is outside the original material",
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

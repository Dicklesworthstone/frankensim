//! Positive layer integrals of the combined pressure/Couette shear stress.
//!
//! For a channel, tau(z) = tau_C - (delta_p*h/L)*(z-z0). Integrating
//! tau^2/mu over each constant-viscosity layer gives a squared mean stress
//! plus its positive within-layer variance. The signed cross term survives
//! within a layer even though it sums to zero over the complete gap.
//!
//! The cold reduction normalizes these two basis functions by the already
//! admitted pressure and shear work. Evaluating their local amplitudes avoids
//! large dimensional stress intermediates and subtraction of separate powers.

use core::ops::Range;

use super::super::stratified::MAX_LAYERS;
use super::{FilmError, GapViscosityProfile, MAX_CHANNELS, ShearChannelReport, State, Sum, fail};

/// Maximum flattened layer powers in one original film evaluation.
///
/// Channels remain in insertion order; each retains its lower-to-upper layer
/// order. A channel without an explicit stratified profile contributes one
/// implicit uniform layer. Scratch storage is bounded by the native geometry.
pub const MAX_LAYER_REPORTS: usize = MAX_CHANNELS * MAX_LAYERS;

#[derive(Debug, Clone, Copy)]
struct Layer {
    pressure_fraction: f64,
    shear_fraction: f64,
    pressure_mean: f64,
    pressure_spread: f64,
    shear_amplitude: f64,
}

#[derive(Debug, Clone)]
pub(super) struct LayerWork {
    ranges: Vec<Range<usize>>,
    layers: Vec<Layer>,
}

impl LayerWork {
    pub(super) fn compile(
        channels: usize,
        profiles: Option<&[GapViscosityProfile]>,
    ) -> Result<Self, FilmError> {
        let count = profiles.map_or(channels, |profiles| {
            profiles.iter().map(|profile| profile.layers().len()).sum()
        });
        if channels > MAX_CHANNELS
            || count > MAX_LAYER_REPORTS
            || profiles.is_some_and(|profiles| profiles.len() != channels)
        {
            return fail("sliding layer layout exceeds the original bounded channel profiles");
        }
        let mut result = Self {
            ranges: Vec::with_capacity(channels),
            layers: Vec::with_capacity(count),
        };
        for channel in 0..channels {
            let start = result.layers.len();
            match profiles {
                Some(profiles) => result.compile_profile(&profiles[channel])?,
                None => result.layers.push(Layer {
                    pressure_fraction: 1.0,
                    shear_fraction: 1.0,
                    pressure_mean: 0.0,
                    pressure_spread: 1.0,
                    shear_amplitude: 1.0,
                }),
            }
            result.ranges.push(start..result.layers.len());
        }
        Ok(result)
    }

    pub(super) fn len(&self) -> usize {
        self.layers.len()
    }

    pub(super) fn range(&self, channel: usize) -> Option<Range<usize>> {
        self.ranges.get(channel).cloned()
    }

    fn compile_profile(&mut self, profile: &GapViscosityProfile) -> Result<(), FilmError> {
        let layers = profile.layers();
        let minimum = layers
            .iter()
            .map(|layer| layer.dynamic_viscosity_pa_s)
            .fold(f64::INFINITY, f64::min);
        let mut lower = [0.0; MAX_LAYERS];
        let mut thickness = [0.0; MAX_LAYERS];
        let mut weights = [0.0; MAX_LAYERS];
        let mut previous = 0.0;
        let mut largest = 0;
        let mut total = Sum::default();
        for (i, layer) in layers.iter().enumerate() {
            lower[i] = previous;
            thickness[i] = layer.upper_fraction - previous;
            weights[i] = thickness[i] * (minimum / layer.dynamic_viscosity_pa_s);
            if weights[i] > weights[largest] {
                largest = i;
            }
            total.add(weights[i]);
            previous = layer.upper_fraction;
        }
        let total = total.value();
        let anchor = lower[largest];
        let mut centers = [0.0; MAX_LAYERS];
        let mut moment = Sum::default();
        for i in 0..layers.len() {
            // Retain local offsets before adding an almost-one global position,
            // exactly as the native stratified pressure mobility does.
            centers[i] = (lower[i] - anchor) + 0.5 * thickness[i];
            moment.add(weights[i] * centers[i]);
        }
        let uniform = layers
            .iter()
            .all(|layer| layer.dynamic_viscosity_pa_s == layers[0].dynamic_viscosity_pa_s);
        let center = if uniform {
            0.5 - anchor
        } else {
            moment.value() / total
        };
        let root_twelve = 12.0_f64.sqrt();
        for i in 0..layers.len() {
            let distance = centers[i] - center;
            let spread = thickness[i] / root_twelve;
            let scale = distance.abs().max(spread);
            if !scale.is_finite() || scale <= 0.0 {
                return fail("sliding layer central moment is not representable");
            }
            // Normalize before squaring: thin layers need no squared absolute
            // gap coordinate and no cancellation of raw second moments.
            let a = distance / scale;
            let b = spread / scale;
            let norm = a.hypot(b);
            let pressure_fraction = profile.layer_pressure_work_fractions()[i];
            let shear_fraction = weights[i] / total;
            let pressure_amplitude = pressure_fraction.sqrt();
            let layer = Layer {
                pressure_fraction,
                shear_fraction,
                pressure_mean: pressure_amplitude * (a / norm),
                pressure_spread: pressure_amplitude * (b / norm),
                shear_amplitude: shear_fraction.sqrt(),
            };
            if !shear_fraction.is_finite()
                || shear_fraction <= 0.0
                || [
                    layer.pressure_mean,
                    layer.pressure_spread,
                    layer.shear_amplitude,
                ]
                .iter()
                .any(|value| !value.is_finite())
            {
                return fail("sliding layer normalized work is not representable");
            }
            self.layers.push(layer);
        }
        Ok(())
    }

    /// Validate the complete report before changing caller storage. Two bounded
    /// arithmetic passes over the solved channels avoid a maximum-sized layer
    /// buffer or heap allocation on small films; the pressure solve is not repeated.
    pub(super) fn evaluate_into(
        &self,
        state: &State,
        channel_work: &[ShearChannelReport],
        pressure_drops: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        if out.len() != self.len()
            || channel_work.len() != self.ranges.len()
            || pressure_drops.len() != self.ranges.len()
        {
            return fail("invalid sliding layer work output shape");
        }
        let mut basis = [Amplitudes::default(); MAX_CHANNELS];
        for (channel, range) in self.ranges.iter().enumerate() {
            let report = channel_work[channel];
            if range.len() == 1 {
                continue;
            }
            let pressure = report.pressure_dissipated_power_w;
            let shear = report.shear_dissipated_power_w;
            let total = report.dissipated_power_w();
            basis[channel] = Amplitudes {
                pressure: pressure.sqrt().copysign(pressure_drops[channel]),
                shear: shear
                    .sqrt()
                    .copysign(state.channels[channel].relative_speed),
            };
            let mut sum = Sum::default();
            for index in range.clone() {
                let power = self.layers[index].power(pressure, shear, basis[channel]);
                if !power.is_finite() || power < 0.0 || (total > 0.0 && power == 0.0) {
                    return fail("positive combined layer work is not representable");
                }
                sum.add(power);
            }
            let sum = sum.value();
            let tolerance = 256.0 * f64::EPSILON * range.len() as f64 * total;
            if !sum.is_finite() || (sum - total).abs() > tolerance {
                return fail(
                    "sliding layer work does not close against native channel dissipation",
                );
            }
        }
        // All fallible work has finished. Reevaluate the same pure local
        // quadratic with identical immutable inputs, then publish its values.
        for (channel, range) in self.ranges.iter().enumerate() {
            let report = channel_work[channel];
            if range.len() == 1 {
                // A full-gap layer exactly recovers the native channel account.
                out[range.start] = report.dissipated_power_w();
            } else {
                for index in range.clone() {
                    out[index] = self.layers[index].power(
                        report.pressure_dissipated_power_w,
                        report.shear_dissipated_power_w,
                        basis[channel],
                    );
                }
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Amplitudes {
    pressure: f64,
    shear: f64,
}

impl Layer {
    fn power(self, pressure: f64, shear: f64, basis: Amplitudes) -> f64 {
        if shear == 0.0 {
            pressure * self.pressure_fraction
        } else if pressure == 0.0 {
            shear * self.shear_fraction
        } else {
            let mean = product_difference(
                basis.shear,
                self.shear_amplitude,
                basis.pressure,
                self.pressure_mean,
            );
            let spread = basis.pressure * self.pressure_spread;
            mean.mul_add(mean, spread * spread)
        }
    }
}

/// Compensate the two products and their subtraction. Near cancellation of the
/// mean stress leaves the positive within-layer variance, not a power deficit.
fn product_difference(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let left = a * b;
    let right = c * d;
    let left_error = a.mul_add(b, -left);
    let right_error = c.mul_add(d, -right);
    let result = left - right;
    let right_virtual = result - left;
    let sum_error = (left - (result - right_virtual)) + (-right - right_virtual);
    result + ((left_error - right_error) + sum_error)
}

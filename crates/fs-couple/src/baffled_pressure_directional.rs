//! Causal first-order microphones over the original baffled surface history.
//!
//! A volume-velocity element Q at range R produces, at t-R/c,
//! p = rho Q'/(2 pi R), u = Rhat (Q/R² + Q'/(c R))/(2 pi).
//! Thus alpha*p-(1-alpha)*rho*c*(front.u) has one acceleration term and
//! one velocity term. The latter is the reactive near field, not an EQ curve.
//! Both terms use the existing anti-alias, derivative and travel-time owners.

use super::path::{HALF, TAPS, taps, window};
use super::{BaffledPressure, CircularOutletReceiver, MAX_HISTORY, RayleighMedium, SurfaceSample};
use fs_bem::near_field::FirstOrder;
use fs_math::det;
use std::collections::BTreeMap;
use std::f64::consts::TAU;

impl BaffledPressure {
    /// Ideal first-order Pa-equivalent output from the full retarded Rayleigh
    /// field, including particle velocity's reactive 1/R² contribution.
    /// `microphone` is the same pressure/velocity projection as the BEM and
    /// harmonic receivers; its front axis points toward the front source.
    /// Uses bandlimited differentiation and fractional travel delays, with the
    /// same reported numerical latency as `new_differentiating`.
    ///
    /// This is a stationary ideal capsule in an inviscid homogeneous half-space,
    /// without measured capsule response, molecular absorption or backreaction.
    /// Spatial refinement and the source bandwidth remain caller obligations.
    ///
    /// # Errors
    /// The existing surface, medium, clock and causal-history admissions apply
    /// to both terms; a failure returns no partially prepared receiver.
    pub fn new_first_order(
        surface: &[SurfaceSample],
        inputs: usize,
        rate: u32,
        position_m: [f64; 3],
        medium: RayleighMedium,
        microphone: FirstOrder,
    ) -> Result<Self, String> {
        Self::new_multirate_first_order(surface, inputs, rate, rate, position_m, medium, microphone)
    }

    /// First-order observation of EVERY accepted mechanics-frame velocity.
    /// Both physical terms are anti-aliased before propagation at output rate.
    /// The velocity term receives the derivative filter's extra numerical
    /// delay, including fractional output-frame delay at non-dyadic ratios;
    /// projecting fields from different physical times would corrupt proximity
    /// response and directional cancellation. No future sample is requested.
    ///
    /// Alpha=1 delegates bit for bit to `new_multirate_differentiating` and
    /// allocates no second bank. Other patterns retain two bounded banks whose
    /// previews succeed together before either accepted history advances.
    /// Stepping and modal-output reconstruction allocate nothing.
    ///
    /// # Errors
    /// Invalid geometry, clocks, medium, unrepresentable weights, or history
    /// exhausted by physical travel plus the explicit numerical alignment.
    pub fn new_multirate_first_order(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        position_m: [f64; 3],
        medium: RayleighMedium,
        microphone: FirstOrder,
    ) -> Result<Self, String> {
        let alpha = microphone.pressure_fraction();
        if alpha == 1.0 {
            return Self::new_multirate_differentiating(
                surface,
                inputs,
                mechanics_rate,
                output_rate,
                position_m,
                medium,
            );
        }

        // Reuse admission and clock construction without preparing and then
        // discarding two complete bandlimited surface banks.
        let mut receiver = Self::new_multirate(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            medium,
        )?;
        let ratio = receiver.decimator.input_frames();
        receiver.differentiate_velocity = true;
        receiver.set_input_ratio(ratio)?;
        receiver.interpolation_delay = HALF;
        let mut reactive = Self::new_multirate(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            medium,
        )?;
        // Consume velocity directly, without finite differencing or integration
        // of an acceleration signal. The initial velocity history is at rest.
        reactive.differentiate_velocity = true;
        reactive.interpolation_delay = HALF;
        let alignment =
            receiver.decimator.delay_output_frames() - reactive.decimator.delay_output_frames();
        let front = microphone.front_axis();
        let window = window();
        let mut acceleration_bins = BTreeMap::<usize, Vec<f64>>::new();
        let mut velocity_bins = BTreeMap::<usize, Vec<f64>>::new();
        for sample in surface {
            let delta: [f64; 3] = std::array::from_fn(|i| position_m[i] - sample.position_m[i]);
            let radius = det::sqrt(delta.iter().map(|x| x * x).sum());
            let projection = (0..3).map(|i| front[i] * delta[i] / radius).sum::<f64>();
            let gain = medium.density * sample.area_m2 / (TAU * radius);
            let acceleration_gain = gain * (alpha - (1.0 - alpha) * projection);
            let velocity_gain = -gain * (1.0 - alpha) * medium.sound_speed * projection / radius;
            let travel = radius * f64::from(output_rate) / medium.sound_speed;
            for (bins, delay, weight) in [
                (&mut acceleration_bins, travel, acceleration_gain),
                (&mut velocity_bins, travel + alignment, velocity_gain),
            ] {
                if !delay.is_finite() || delay >= (MAX_HISTORY - TAPS) as f64 {
                    return Err(
                        "first-order propagation/alignment exceeds the causal history budget"
                            .into(),
                    );
                }
                let first = delay.floor() as usize;
                let interpolation = taps(delay - first as f64, &window);
                for (tap, interpolation_weight) in interpolation.iter().enumerate() {
                    let row = bins.entry(first + tap).or_insert_with(|| vec![0.0; inputs]);
                    for (value, shape) in row.iter_mut().zip(&sample.mode_shape) {
                        *value += weight * interpolation_weight * shape;
                    }
                }
            }
        }
        install(&mut receiver, acceleration_bins)?;
        install(&mut reactive, velocity_bins)?;
        receiver.delay_samples = (
            receiver.delay_samples.0.min(reactive.delay_samples.0),
            receiver.delay_samples.1.max(reactive.delay_samples.1),
        );
        receiver.directional_velocity = Some(Box::new(reactive));
        Ok(receiver)
    }

    /// A uniform circular outlet through the same first-order receiver. Input
    /// remains outward volume flow [m³/s]; the original positive-area disk
    /// quadrature, geometry and spatial phase-resolution checks are retained.
    /// The declared acoustic band must be within 0.45 of the output rate.
    ///
    /// # Errors
    /// Original disk/clock/history admissions or an unresolved output band.
    pub fn circular_outlet_first_order(
        radius_m: f64,
        mechanics_rate: u32,
        output_rate: u32,
        receiver: CircularOutletReceiver,
        medium: RayleighMedium,
        microphone: FirstOrder,
    ) -> Result<Self, String> {
        if receiver.maximum_frequency_hz > 0.45 * f64::from(output_rate) {
            return Err(
                "first-order circular receiver band exceeds its output reconstruction".into(),
            );
        }
        let surface = Self::circular_surface(radius_m, receiver, medium)?;
        Self::new_multirate_first_order(
            &surface,
            1,
            mechanics_rate,
            output_rate,
            receiver.position_m,
            medium,
            microphone,
        )
    }
}

fn install(receiver: &mut BaffledPressure, bins: BTreeMap<usize, Vec<f64>>) -> Result<(), String> {
    if bins.values().flatten().any(|v| !v.is_finite()) {
        return Err("first-order Rayleigh kernel overflow".into());
    }
    let first = *bins
        .first_key_value()
        .ok_or("empty first-order Rayleigh kernel")?
        .0;
    let last = *bins
        .last_key_value()
        .ok_or("empty first-order Rayleigh kernel")?
        .0;
    let frames = last + 1;
    let samples = frames
        .checked_mul(receiver.modes)
        .ok_or("first-order history size overflow")?;
    receiver.kernels = bins.into_iter().collect();
    receiver.history = vec![0.0; samples];
    receiver.history_frames = frames;
    receiver.delay_samples = (first, last);
    receiver.prepare_taps();
    Ok(())
}

#[cfg(test)]
#[path = "baffled_pressure_directional_tests.rs"]
mod tests;

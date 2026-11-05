//! Causal first-order microphones over the original baffled surface history.
//!
//! A volume-velocity element Q at range R produces, at t-R/c,
//! p = rho Q'/(2 pi R), u = Rhat (Q/R² + Q'/(c R))/(2 pi).
//! Thus alpha*p-(1-alpha)*rho*c*(front.u) has one acceleration term and
//! one velocity term. The latter is the reactive near field, not an EQ curve.
//! Both terms use the existing anti-alias, derivative and travel-time owners.

use super::path::{HALF, TAPS, KernelBins, taps, window};
use super::{BaffledPressure, CircularOutletReceiver, MAX_HISTORY, RayleighMedium, SurfaceSample};
use fs_bem::near_field::FirstOrder;
use fs_math::{dd::Dd, det};

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
        let mut acceleration_bins = KernelBins::default();
        let mut velocity_bins = KernelBins::default();
        add_paths(surface, output_rate, position_m, medium, microphone, alignment,
            Dd::ONE, &window(), &mut acceleration_bins, &mut velocity_bins, &mut || false)?;
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

/// Shared point/aperture path preparation. Only the receiver quadrature weight
/// varies: both images retain the same reactive term and causal alignment.
#[allow(clippy::too_many_arguments)]
pub(super) fn add_paths(
    surface: &[SurfaceSample],
    output_rate: u32,
    position_m: [f64; 3],
    medium: RayleighMedium,
    microphone: FirstOrder,
    alignment: f64,
    aperture_weight: Dd,
    kernel_window: &[f64; TAPS],
    acceleration_bins: &mut KernelBins,
    velocity_bins: &mut KernelBins,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), String> {
    let alpha = microphone.pressure_fraction();
    let front = microphone.front_axis();
    for (index, sample) in surface.iter().enumerate() {
        if index % 64 == 0 && cancelled() {
            return Err("Rayleigh aperture preparation cancelled".into());
        }
        let delta: [f64; 3] = std::array::from_fn(|i| position_m[i] - sample.position_m[i]);
        let radius = det::sqrt(delta.iter().map(|x| x * x).sum());
        let projection = (0..3).map(|i| front[i] * delta[i] / radius).sum::<f64>();
        let gain = medium.density * sample.area_m2 / (TAU * radius);
        let acceleration_gain = gain * (alpha - (1.0 - alpha) * projection);
        let velocity_gain = -gain * (1.0 - alpha) * medium.sound_speed * projection / radius;
        let travel = radius * f64::from(output_rate) / medium.sound_speed;
        for (bank, delay, weight) in [
            (0, travel, acceleration_gain),
            (1, travel + alignment, velocity_gain),
        ] {
            if bank == 1 && alpha == 1.0 { continue; }
            if !delay.is_finite() || delay < 1.0 || delay >= (MAX_HISTORY - TAPS) as f64 {
                return Err("first-order propagation/alignment exceeds the causal history budget".into());
            }
            let first = delay.floor() as usize;
            let interpolation = taps(delay - first as f64, kernel_window);
            let bins = if bank == 0 { &mut *acceleration_bins } else { &mut *velocity_bins };
            for (tap, interpolation_weight) in interpolation.iter().enumerate() {
                let scale = Dd::from_f64(weight) * Dd::from_f64(*interpolation_weight) * aperture_weight;
                bins.add_scaled(first + tap, scale, &sample.mode_shape);
            }
        }
    }
    Ok(())
}

pub(super) fn install(receiver: &mut BaffledPressure, bins: KernelBins) -> Result<(), String> {
    let kernels = bins.finish()?;
    let first = kernels
        .first()
        .ok_or("empty first-order Rayleigh kernel")?
        .0;
    let last = kernels
        .last()
        .ok_or("empty first-order Rayleigh kernel")?
        .0;
    let frames = last + 1;
    let samples = frames
        .checked_mul(receiver.modes)
        .ok_or("first-order history size overflow")?;
    receiver.kernels = kernels;
    receiver.history = vec![0.0; samples];
    receiver.history_frames = frames;
    receiver.delay_samples = (first, last);
    receiver.prepare_taps();
    Ok(())
}

#[cfg(test)]
#[path = "baffled_pressure_directional_tests.rs"]
mod tests;

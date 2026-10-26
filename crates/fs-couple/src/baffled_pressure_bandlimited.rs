//! Bandlimited fractional travel times for the existing one-way Rayleigh owner.
//! Windowed-sinc interpolation: Julius O. Smith, Physical Audio Signal Processing,
//! https://www.dsprelated.com/freebooks/pasp/Windowed_Sinc_Interpolation.html
//! This changes the delay approximation, never the source geometry or SI gain.
use super::{BaffledPressure, RayleighMedium, SurfaceSample, MAX_HISTORY};
use fs_math::det;
use std::f64::consts::TAU;
use super::path::{HALF, TAPS, taps, window};

impl BaffledPressure {
    /// Rayleigh pressure with a 65-tap Blackman-windowed sinc fractional delay.
    /// The original `new` remains the two-tap linear comparison image.
    ///
    /// The added interpolation latency is 32 output samples, separate from
    /// physical travel and included in `filter_delay_samples`. No kernel tap
    /// precedes floor(travel_time * rate); no future sample is requested.
    /// DC gain is normalized to retain physical pressure units. This is a
    /// finite approximation, not an exact bandlimited reconstruction: focused
    /// frequency-grid tests bound complex error through 0.45 of the output rate.
    /// It is not an anti-alias filter, structural bandwidth extension, room
    /// response, acoustic backreaction or a perceptual calibration.
    ///
    /// # Errors
    /// Same admission as `new`, with the wider kernel in the same history budget.
    pub fn new_bandlimited(surface: &[SurfaceSample], inputs: usize, rate: u32,
        position_m: [f64; 3], medium: RayleighMedium) -> Result<Self, String> {
        Self::new_multirate_bandlimited(surface, inputs, rate, rate, position_m, medium)
    }

    /// Observe all mechanics substeps, anti-alias as before, then propagate with
    /// bandlimited fractional delays on the output clock. Construction is cold;
    /// stepping and modal pressure reconstruction use the existing prepared loop.
    /// At 4x oversampling, filter latency is 44 + 32 = 76 output samples.
    /// The original `new_multirate` keeps its exact linear-delay arithmetic.
    ///
    /// # Errors
    /// Invalid geometry, clocks, medium or the bounded causal history allocation.
    pub fn new_multirate_bandlimited(surface: &[SurfaceSample], inputs: usize,
        mechanics_rate: u32, output_rate: u32, position_m: [f64; 3],
        medium: RayleighMedium) -> Result<Self, String> {
        // Reuse complete existing admission, clock setup and buffer ownership.
        // Nothing is installed into a running receiver or published on failure.
        let mut receiver = Self::new_multirate(surface, inputs, mechanics_rate,
            output_rate, position_m, medium)?;
        let window = window();
        let mut bins = std::collections::BTreeMap::<usize, Vec<f64>>::new();
        for p in surface {
            let radius = det::sqrt((0..3).map(|i| (position_m[i] - p.position_m[i]).powi(2)).sum());
            let delay = radius * f64::from(output_rate) / medium.sound_speed;
            if delay >= (MAX_HISTORY - TAPS) as f64 {
                return Err("bandlimited propagation kernel exceeds the causal history budget".into());
            }
            let first = delay.floor() as usize;
            let h = taps(delay - first as f64, &window);
            let gain = medium.density * p.area_m2 / (TAU * radius);
            for (tap, weight) in h.iter().enumerate() {
                let row = bins.entry(first + tap).or_insert_with(|| vec![0.0; inputs]);
                for (i, shape) in p.mode_shape.iter().enumerate() {
                    row[i] += gain * weight * shape;
                }
            }
        }
        if bins.values().flatten().any(|v| !v.is_finite()) {
            return Err("bandlimited radiation kernel overflow".into());
        }
        let first = *bins.first_key_value().ok_or("empty radiation kernel")?.0;
        let last = *bins.last_key_value().ok_or("empty radiation kernel")?.0;
        let frames = last + 1;
        let samples = frames.checked_mul(inputs).ok_or("microphone history size overflow")?;
        receiver.kernels = bins.into_iter().collect();
        receiver.history = vec![0.0; samples];
        receiver.history_frames = frames;
        receiver.delay_samples = (first, last);
        receiver.interpolation_delay = HALF;
        Ok(receiver)
    }

    /// Additional fixed fractional-delay latency, not distance or anti-alias lag.
    pub fn interpolation_delay_samples(&self) -> usize { self.interpolation_delay }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn medium() -> RayleighMedium { RayleighMedium { density: 1.2, sound_speed: 320.0 } }
    fn patch(shape: Vec<f64>) -> SurfaceSample {
        SurfaceSample { position_m: [0.0; 3], area_m2: 0.1, mode_shape: shape }
    }
    #[test]
    fn complex_response_retains_upper_partials_and_fractional_phase() {
        let w = window();
        for j in 0..=64 {
            let fraction = j as f64 / 65.0;
            let h = taps(fraction, &w);
            assert!((h.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            assert_eq!((h[0], h[TAPS - 1]), (0.0, 0.0));
            for bin in 0..=180 {
                let omega = TAU * 0.45 * bin as f64 / 180.0;
                let (mut re, mut im) = (0.0, 0.0);
                for (n, value) in h.iter().enumerate() {
                    let phase = omega * (n as f64 - HALF as f64 - fraction);
                    re += value * phase.cos(); im -= value * phase.sin();
                }
                assert!((re - 1.0).hypot(im) < 4e-4, "fraction={fraction} omega={omega}");
            }
        }
        for fraction in [f64::EPSILON, 1.0 - f64::EPSILON] {
            let h = taps(fraction, &w);
            assert!(h.iter().all(|v| v.is_finite()));
            assert!((h.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            let index = if fraction < 0.5 { HALF } else { HALF + 1 };
            assert!((h[index] - 1.0).abs() < 1e-14);
        }
    }
    #[test]
    fn actual_receiver_matches_independent_delayed_discrete_acceleration() {
        let rate = 48_000;
        let delay = 75.5;
        let position = [0.0, 0.0, delay * 320.0 / f64::from(rate)];
        let gain = 1.2 * 0.1 / (TAU * position[2]);
        for frequency in [4800.0, 10000.0, 15000.0, 21600.0] {
            let mut mic = BaffledPressure::new_bandlimited(&[patch(vec![1.0])], 1, rate,
                position, medium()).unwrap();
            let omega = TAU * frequency / f64::from(rate);
            let amplitude = gain * 2.0 * f64::from(rate) * (0.5 * omega).sin();
            let (mut error, mut signal) = (0.0, 0.0);
            for n in 0..2400 {
                let p = mic.step(&[(omega * n as f64).sin()]).unwrap();
                if n > 400 {
                    let exact = amplitude * (omega * (n as f64 - 0.5 - delay - HALF as f64)).cos();
                    error += (p - exact).powi(2); signal += exact.powi(2);
                }
            }
            assert!((error / signal).sqrt() < 4e-4, "frequency={frequency}");
        }
    }
    #[test]
    fn integer_delays_are_the_original_waveform_shifted_by_the_declared_latency() {
        let surface = [patch(vec![1.0])];
        let mut old = BaffledPressure::new(&surface, 1, 32_000, [0.0, 0.0, 1.0], medium()).unwrap();
        let mut new = BaffledPressure::new_bandlimited(&surface, 1, 32_000,
            [0.0, 0.0, 1.0], medium()).unwrap();
        assert_eq!(old.filter_delay_samples(), 0.0);
        assert_eq!(new.filter_delay_samples(), HALF as f64);
        let mut previous = Vec::new();
        for n in 0..600 {
            let v = (n as f64 * 0.19).sin();
            previous.push(old.step(&[v]).unwrap());
            let p = new.step(&[v]).unwrap();
            if n < HALF { assert_eq!(p, 0.0); }
            else { assert_eq!(p, previous[n - HALF]); }
        }
    }
    #[test]
    fn causal_support_and_constant_acceleration_retain_physical_units() {
        let surface = [patch(vec![1.0])];
        let z = 1.005;
        let mut mic = BaffledPressure::new_bandlimited(&surface, 1, 32_000,
            [0.0, 0.0, z], medium()).unwrap();
        let gain = 1.2 * 0.1 / (TAU * z);
        for n in 0..600 {
            let p = mic.step(&[(n + 1) as f64 / 32_000.0]).unwrap();
            if n < mic.delay_samples.0 { assert_eq!(p, 0.0); }
            if n > mic.delay_samples.1 { assert!((p / gain - 1.0).abs() < 1e-11); }
        }
    }
    #[test]
    fn multirate_modal_sum_and_rejected_steps_preserve_one_output_clock() {
        let surface = [patch(vec![1.0, -0.5])];
        let mut a = BaffledPressure::new_multirate_bandlimited(&surface, 2, 192_000, 48_000,
            [0.0, 0.0, 1.005], medium()).unwrap();
        let mut b = BaffledPressure::new_multirate_bandlimited(&surface, 2, 192_000, 48_000,
            [0.0, 0.0, 1.005], medium()).unwrap();
        assert_eq!(a.filter_delay_samples(), 76.0);
        assert_eq!(a.interpolation_delay_samples(), 32);
        let pointers = (a.history.as_ptr(), a.acceleration.as_ptr(), a.previous_velocity.as_ptr());
        let mut contributions = [0.0; 2];
        let mut heard = false;
        for frame in 0..700 {
            assert!(a.step_trace(&[f64::NAN; 8]).is_err());
            assert!(a.step_trace(&[f64::MAX; 8]).is_err());
            let mut trace = [0.0; 8];
            for s in 0..4 {
                trace[s * 2] = ((frame * 4 + s) as f64 * 0.03).sin();
                trace[s * 2 + 1] = ((frame * 4 + s) as f64 * 0.071).sin();
            }
            let p = a.step_trace(&trace).unwrap();
            assert_eq!(p, b.step_trace(&trace).unwrap());
            a.mode_pressures(&mut contributions).unwrap();
            assert!((contributions.iter().sum::<f64>() - p).abs() < 1e-10);
            heard |= p.abs() > 1e-5;
        }
        assert!(heard);
        assert_eq!(pointers, (a.history.as_ptr(), a.acceleration.as_ptr(), a.previous_velocity.as_ptr()));
    }
    #[test]
    fn wider_kernels_are_included_in_history_admission() {
        let surface = [patch(vec![1.0])];
        let position = [0.0, 0.0, (MAX_HISTORY - 10) as f64 * 320.0 / 32_000.0];
        assert!(BaffledPressure::new(&surface, 1, 32_000, position, medium()).is_ok());
        assert!(BaffledPressure::new_bandlimited(&surface, 1, 32_000, position, medium()).is_err());
        assert!(BaffledPressure::new_multirate_bandlimited(&surface, 1, 48_001, 48_000,
            [0.0, 0.0, 1.0], medium()).is_err());
    }
}

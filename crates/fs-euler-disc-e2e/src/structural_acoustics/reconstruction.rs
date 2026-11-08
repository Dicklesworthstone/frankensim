//! Bounded fractional-delay reconstruction on the original source clock.

use fs_math::det;

/// Finite window applied to the retarded source's normalized sinc stencil.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum RetardedReconstructionWindow {
    /// Sinc window with a comparatively narrow transition near source Nyquist.
    Lanczos = 0,
    /// Three-term Blackman window: lower passband ripple, but a wider transition.
    /// The caller must account for this transition in its source guard band.
    Blackman = 1,
}

/// One reusable allocation-free stencil; all source channels share its weights.
pub(super) struct PreparedWindowedSinc {
    radius: usize,
    window: RetardedReconstructionWindow,
    step_sine: f64,
    step_cosine: f64,
    weights: [f64; 128],
}

impl PreparedWindowedSinc {
    pub(super) fn new(radius: u8, window: RetardedReconstructionWindow) -> Option<Self> {
        if !matches!(radius, 8 | 16 | 32 | 64) {
            return None;
        }
        let (step_sine, step_cosine) = det::sin_cos(core::f64::consts::PI / f64::from(radius));
        Some(Self {
            radius: usize::from(radius),
            window,
            step_sine,
            step_cosine,
            weights: [0.0; 128],
        })
    }

    pub(super) fn radius(&self) -> usize {
        self.radius
    }

    pub(super) fn tap_count(&self) -> usize {
        2 * self.radius
    }

    /// Weights for offsets `1-radius ..= radius` from `floor(source_frame)`.
    /// Five deterministic scalar trig evaluations suffice for any admitted
    /// radius (one sine and two shared sine/cosine reductions). Recurrences
    /// run outward from the nearest samples to avoid cancellation near them.
    pub(super) fn weights(&mut self, phase: f64) -> Option<&[f64]> {
        if !(phase.is_finite() && (0.0..1.0).contains(&phase)) {
            return None;
        }
        let count = self.tap_count();
        if phase == 0.0 {
            self.weights[..count].fill(0.0);
            self.weights[self.radius - 1] = 1.0;
            return Some(&self.weights[..count]);
        }
        let pi = core::f64::consts::PI;
        // Reflection keeps sin(pi * phase) accurate near either integer.
        let mut carrier = det::sin(pi * phase.min(1.0 - phase));
        let step = pi / self.radius as f64;
        let (mut left_sine, mut left_cosine) = det::sin_cos(step * phase);
        let (mut right_sine, mut right_cosine) = det::sin_cos(step * (1.0 - phase));
        for offset in 0..self.radius {
            self.weights[self.radius - 1 - offset] =
                self.weight(phase + offset as f64, carrier, left_sine, left_cosine);
            self.weights[self.radius + offset] = self.weight(
                (1.0 - phase) + offset as f64,
                carrier,
                right_sine,
                right_cosine,
            );
            (left_sine, left_cosine) = self.advance(left_sine, left_cosine);
            (right_sine, right_cosine) = self.advance(right_sine, right_cosine);
            carrier = -carrier;
        }
        // A finite sinc window is not exactly a partition of unity. Normalize
        // once, before applying these same weights to every radiation channel.
        let mut sum = 0.0_f64;
        let mut compensation = 0.0_f64;
        for &weight in &self.weights[..count] {
            let next = sum + weight;
            compensation += if sum.abs() >= weight.abs() {
                (sum - next) + weight
            } else {
                (weight - next) + sum
            };
            sum = next;
        }
        let normalization = sum + compensation;
        if !(normalization > 0.0 && normalization.is_finite()) {
            return None;
        }
        for weight in &mut self.weights[..count] {
            *weight /= normalization;
        }
        Some(&self.weights[..count])
    }

    fn advance(&self, sine: f64, cosine: f64) -> (f64, f64) {
        (
            sine.mul_add(self.step_cosine, cosine * self.step_sine),
            cosine.mul_add(self.step_cosine, -sine * self.step_sine),
        )
    }

    fn weight(&self, distance: f64, carrier: f64, sine: f64, cosine: f64) -> f64 {
        if distance >= self.radius as f64 {
            return 0.0;
        }
        let pix = core::f64::consts::PI * distance;
        let sinc = if distance < 1.0e-3 {
            sinc_near_zero(pix)
        } else {
            carrier / pix
        };
        let window = match self.window {
            RetardedReconstructionWindow::Lanczos => {
                let angle = pix / self.radius as f64;
                if distance < 1.0e-3 {
                    sinc_near_zero(angle)
                } else {
                    sine / angle
                }
            }
            RetardedReconstructionWindow::Blackman => {
                // 0.42 + 0.5 cos(angle) + 0.08 cos(2 angle), factored to
                // retain its zero endpoint without a second trig recurrence.
                (1.0 + cosine).max(0.0) * 0.16_f64.mul_add(cosine, 0.34)
            }
        };
        sinc * window
    }
}

fn sinc_near_zero(angle: f64) -> f64 {
    let square = angle * angle;
    square.mul_add(
        square.mul_add(square.mul_add(-1.0 / 5_040.0, 1.0 / 120.0), -1.0 / 6.0),
        1.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const WINDOWS: [RetardedReconstructionWindow; 2] = [
        RetardedReconstructionWindow::Lanczos,
        RetardedReconstructionWindow::Blackman,
    ];

    // Independent direct formula: no recurrence, reflection, polynomial sinc,
    // factored Blackman window, or production normalization accumulator.
    fn direct_weights(radius: usize, window: RetardedReconstructionWindow, phase: f64) -> Vec<f64> {
        let mut weights: Vec<_> = (1 - radius as isize..=radius as isize)
            .map(|tap| {
                let x = phase - tap as f64;
                if x == 0.0 {
                    return 1.0;
                }
                if x.abs() >= radius as f64 {
                    return 0.0;
                }
                let pix = core::f64::consts::PI * x;
                let angle = pix / radius as f64;
                let window = match window {
                    RetardedReconstructionWindow::Lanczos => det::sin(angle) / angle,
                    RetardedReconstructionWindow::Blackman => {
                        0.42 + 0.5 * det::cos(angle) + 0.08 * det::cos(2.0 * angle)
                    }
                };
                (det::sin(pix) / pix) * window
            })
            .collect();
        let sum: f64 = weights.iter().sum();
        for weight in &mut weights {
            *weight /= sum;
        }
        weights
    }

    fn complex_tone_error(weights: &[f64], phase: f64, frequency: f64) -> f64 {
        let radius = weights.len() / 2;
        let (mut real, mut imaginary) = (0.0_f64, 0.0_f64);
        for (index, &weight) in weights.iter().enumerate() {
            let x = index as f64 + 1.0 - radius as f64 - phase;
            let (sine, cosine) = det::sin_cos(core::f64::consts::TAU * frequency * x);
            real = weight.mul_add(cosine, real);
            imaginary = weight.mul_add(sine, imaginary);
        }
        det::sqrt((real - 1.0).powi(2) + imaginary.powi(2))
    }

    #[test]
    fn admitted_radii_and_phases_are_bounded() {
        for window in WINDOWS {
            for radius in [0, 1, 7, 9, 15, 17, 31, 33, 63, 65, 128, 255] {
                assert!(PreparedWindowedSinc::new(radius, window).is_none());
            }
            for radius in [8, 16, 32, 64] {
                let mut kernel = PreparedWindowedSinc::new(radius, window).unwrap();
                assert_eq!(kernel.radius(), usize::from(radius));
                assert_eq!(kernel.tap_count(), 2 * usize::from(radius));
                for phase in [-1.0, -f64::EPSILON, 1.0, f64::INFINITY, f64::NAN] {
                    assert!(kernel.weights(phase).is_none());
                }
                assert!(kernel.weights(0.5).is_some());
            }
        }
    }

    #[test]
    fn prepared_weights_match_independent_direct_sinc_for_all_radii_and_windows() {
        for window in WINDOWS {
            for radius in [8, 16, 32, 64] {
                let mut kernel = PreparedWindowedSinc::new(radius, window).unwrap();
                let phases = (0..64).map(|phase| f64::from(phase) / 64.0).chain([
                    1.0e-12,
                    1.0e-6,
                    1.0e-3,
                    1.0 - 1.0e-3,
                    1.0 - 1.0e-6,
                    1.0 - 1.0e-12,
                ]);
                for phase in phases {
                    let reference = direct_weights(usize::from(radius), window, phase);
                    let weights = kernel.weights(phase).unwrap();
                    for (&actual, expected) in weights.iter().zip(reference) {
                        assert!(
                            (actual - expected).abs() < 5.0e-13,
                            "{window:?} radius={radius} phase={phase}: {actual} != {expected}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn normalization_preserves_dc_and_integer_phases_are_exact() {
        for window in WINDOWS {
            for radius in [8, 16, 32, 64] {
                let mut kernel = PreparedWindowedSinc::new(radius, window).unwrap();
                for phase in 0..64 {
                    let weights = kernel.weights(f64::from(phase) / 64.0).unwrap();
                    assert!((weights.iter().sum::<f64>() - 1.0).abs() < 2.0e-15);
                    if phase == 0 {
                        for (tap, &weight) in weights.iter().enumerate() {
                            let expected = if tap == usize::from(radius) - 1 {
                                1.0
                            } else {
                                0.0
                            };
                            assert_eq!(weight, expected);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn complex_tone_grid_measures_passband_and_transition_tradeoffs() {
        // Discrete regression grid, not a continuum response certificate.
        let bands = [0.375, 0.45, 11.0 / 24.0, 23.0 / 48.0, 0.49];
        let mut measured = [[[0.0_f64; 5]; 4]; 2];
        for (window_index, window) in WINDOWS.into_iter().enumerate() {
            for (radius_index, radius) in [8, 16, 32, 64].into_iter().enumerate() {
                let mut kernel = PreparedWindowedSinc::new(radius, window).unwrap();
                for phase_index in 0..16 {
                    let phase = f64::from(phase_index) / 16.0;
                    let weights = kernel.weights(phase).unwrap();
                    for (band_index, band) in bands.into_iter().enumerate() {
                        for frequency_index in 0..=256 {
                            let error = complex_tone_error(
                                weights,
                                phase,
                                band * f64::from(frequency_index) / 256.0,
                            );
                            let maximum = &mut measured[window_index][radius_index][band_index];
                            *maximum = maximum.max(error);
                        }
                    }
                }
                println!(
                    "{window:?} radius={radius}: {:?}",
                    measured[window_index][radius_index]
                );
            }
        }
        // Larger admitted stencils improve each sampled band for both windows.
        for window in measured {
            for radii in window.windows(2) {
                for band in 0..bands.len() {
                    assert!(radii[1][band] < radii[0][band]);
                }
            }
        }
        let lanczos64 = measured[0][3];
        let blackman64 = measured[1][3];
        for (actual, ceiling) in lanczos64
            .into_iter()
            .zip([0.00024, 0.0015, 0.0022, 0.0092, 0.022])
        {
            assert!(actual < ceiling);
        }
        for (actual, ceiling) in blackman64
            .into_iter()
            .zip([0.000009, 0.000095, 0.00015, 0.0007, 0.165])
        {
            assert!(actual < ceiling);
        }
        for band in 0..4 {
            assert!(blackman64[band] < 0.1 * lanczos64[band]);
        }
        // Blackman's lower ripple trades away more of the last Nyquist margin.
        assert!(blackman64[4] > 5.0 * lanczos64[4]);
    }

    #[test]
    fn time_warped_high_frequency_tones_retain_source_clock_phase() {
        let source_frequency = 22_000.0 / 48_000.0;
        let source: Vec<_> = (0..1_024)
            .map(|frame| det::sin_cos(core::f64::consts::TAU * source_frequency * f64::from(frame)))
            .collect();
        for warp in [1.0, 1.08 / 0.92, 0.94 / 1.08] {
            for window in WINDOWS {
                let mut kernel = PreparedWindowedSinc::new(64, window).unwrap();
                let mut maximum = 0.0_f64;
                for frame in 0..512 {
                    // Source 48 kHz, receiver 96 kHz, with a smoothly varying
                    // source-clock rate as well as approaching/receding cases.
                    let t = f64::from(frame);
                    let coordinate = 128.137 + 0.5 * warp * t + 0.00001 * t * t;
                    let base = coordinate.floor() as usize;
                    let weights = kernel.weights(coordinate - base as f64).unwrap();
                    let (mut sine, mut cosine) = (0.0_f64, 0.0_f64);
                    for (&weight, &(sample_sine, sample_cosine)) in
                        weights.iter().zip(&source[base + 1 - 64..base + 1 + 64])
                    {
                        sine = weight.mul_add(sample_sine, sine);
                        cosine = weight.mul_add(sample_cosine, cosine);
                    }
                    let (expected_sine, expected_cosine) =
                        det::sin_cos(core::f64::consts::TAU * source_frequency * coordinate);
                    let error = det::sqrt(
                        (sine - expected_sine).powi(2) + (cosine - expected_cosine).powi(2),
                    );
                    maximum = maximum.max(error);
                }
                let ceiling = match window {
                    RetardedReconstructionWindow::Lanczos => 0.0017,
                    RetardedReconstructionWindow::Blackman => 0.00013,
                };
                println!("{window:?} time-warp={warp}: maximum complex error={maximum}");
                assert!(maximum < ceiling);
            }
        }
    }

    #[test]
    fn exact_nyquist_remains_phase_ambiguous_for_every_finite_window() {
        for window in WINDOWS {
            for radius in [8, 16, 32, 64] {
                let mut kernel = PreparedWindowedSinc::new(radius, window).unwrap();
                let weights = kernel.weights(0.5).unwrap();
                assert!(complex_tone_error(weights, 0.5, 0.5) > 0.99);
            }
        }
    }
}

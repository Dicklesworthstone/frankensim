//! Positive temporal Jones coherency at a deep, resolved optical beat null.
use super::*;
use fs_math::dd::Dd;
use fs_render::camera::PolarizedDiffractionCoherency;

const C: f64 = 299_792_458.0;
const F0: f64 = 562_949_953_421_312.0;
const DISTANCE: f64 = 0.03;

// On-axis pupil normalization varies with wavelength. Author each coefficient
// as two original waves carrying the exact DD product of the OTHER pupils.
// All three detector amplitudes then share h0*h1*h2, without a rounded inverse
// pupil gain manufacturing a much larger residual than the short-window beat.
fn compensated_waves(camera: &PhysicalCamera) -> ([[CoherentJonesPupilWave; 2]; 3], f64) {
    let h: [f64; 3] = std::array::from_fn(|i| {
        let wavelength =
            (Dd::from_f64(C) / Dd::from_f64(F0 + i as f64) * Dd::from_f64(1e9)).to_f64();
        let value = camera
            .in_focus_diffraction_amplitude([0.0; 2], wavelength)
            .unwrap();
        assert_eq!(value.im, 0.0);
        value.re
    });
    let waves = std::array::from_fn(|i| {
        let a = Dd::from_f64(h[(i + 1) % 3])
            * Dd::from_f64(h[(i + 2) % 3])
            * Dd::from_f64([1.0, -2.0, 1.0][i]);
        [a.hi, a.lo].map(|part| jones_wave([0.0; 2], [C64::from_re(part), C64::new(0.0, part)]))
    });
    let amplitude = (Dd::from_f64(h[0]) * Dd::from_f64(h[1]) * Dd::from_f64(h[2])).to_f64();
    (waves, amplitude)
}

fn shutter(duration: f64) -> CoherentDiffractionShutter {
    let start =
        Dd::from_f64(DISTANCE) / Dd::from_f64(C) - Dd::from_f64(duration) * Dd::from_f64(0.5);
    CoherentDiffractionShutter::try_new(start.to_f64(), duration).unwrap()
}

// Independent real-time envelope: |1-2 exp(-i2pi t)+exp(-i4pi t)|^2
// = 16 sin^4(pi t). Integrate its t^4 limit analytically, retaining the actual
// binary64 shutter's small offset from the retarded detector clock. The omitted
// relative correction is below 1e-16 for these 1 Hz, <=1 ns intervals.
fn expected_axis_power(amplitude: f64, shutter: CoherentDiffractionShutter) -> f64 {
    let duration = shutter.duration_s();
    let centre = (Dd::from_f64(shutter.start_s()) - Dd::from_f64(DISTANCE) / Dd::from_f64(C)
        + Dd::from_f64(duration) * Dd::from_f64(0.5))
    .to_f64();
    let t2 = duration * duration;
    let c2 = centre * centre;
    16.0 * core::f64::consts::PI.powi(4)
        * amplitude
        * amplitude
        * (c2 * c2 + 0.5 * c2 * t2 + t2 * t2 / 80.0)
}

fn circular_matrix(actual: PolarizedDiffractionCoherency, power: f64, tolerance: f64) {
    assert!(power > 0.0);
    let [x, y] = actual.axis_powers();
    let xy = actual.cross();
    for (value, expected) in [x, y, xy.re, xy.im].into_iter().zip([1.0, 1.0, 0.0, -1.0]) {
        assert!(
            (value / power - expected).abs() < tolerance,
            "{actual:?}, expected {power}"
        );
    }
    assert!((actual.stokes().degree() - 1.0).abs() < 2e-14);
}

#[test]
fn g1_g3_short_shutter_keeps_positive_jones_matrix_at_three_tone_null() {
    let camera = focused();
    let (waves, amplitude) = compensated_waves(&camera);
    let matched = CoherentJonesDetector::analyzer([C64::ONE, C64::new(0.0, 1.0)]).unwrap();
    let rejected = CoherentJonesDetector::analyzer([C64::ONE, C64::new(0.0, -1.0)]).unwrap();
    for duration in [1e-9, 1e-12] {
        let shutter = shutter(duration);
        let expected = expected_axis_power(amplitude, shutter);
        let tones: [_; 3] = std::array::from_fn(|i| jones_tone(F0 + i as f64, &waves[i], &[]));
        for tones in [tones, [tones[2], tones[1], tones[0]]] {
            let query = camera
                .prepare_coherent_jones_diffraction_exposure(
                    &tones,
                    DISTANCE,
                    shutter,
                    EXPOSURE_LIMITS,
                    || false,
                )
                .unwrap();
            let value = query.observe_coherency([0.0; 2], || false).unwrap();
            // The 1 ps signal is ~1e-48 of an isolated carrier's power.
            // Its ~1e-24 field residue remains resolved by DD source products.
            circular_matrix(value.coherency_density_w, expected, 2e-6);
            assert_eq!(value.fluctuation_coherency_density_w.stokes().i, 0.0);
            for detector in [JONES_TOTAL, matched] {
                let power = query
                    .observe([0.0; 2], detector, || false)
                    .unwrap()
                    .power_density_w;
                assert!((power / (2.0 * expected) - 1.0).abs() < 2e-6);
            }
            assert_eq!(
                query
                    .observe([0.0; 2], rejected, || false)
                    .unwrap()
                    .power_density_w,
                0.0
            );
        }
    }
}

#[test]
fn g3_short_jones_shutter_preserves_source_covariance_columns() {
    let camera = focused();
    let (waves, amplitude) = compensated_waves(&camera);
    let shutter = shutter(1e-12);
    let shared: [_; 3] = std::array::from_fn(|i| [JonesPupilCoherenceFactor { waves: &waves[i] }]);
    let tones: [_; 3] = std::array::from_fn(|i| jones_tone(F0 + i as f64, &[], &shared[i]));
    let correlated = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            DISTANCE,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe_coherency([0.0; 2], || false)
        .unwrap();
    assert_eq!(correlated.mean_coherency_density_w.stokes().i, 0.0);
    circular_matrix(
        correlated.fluctuation_coherency_density_w,
        expected_axis_power(amplitude, shutter),
        2e-6,
    );
    let independent: [[JonesPupilCoherenceFactor<'_>; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|bank| JonesPupilCoherenceFactor {
            waves: if i == bank { &waves[i] } else { &[] },
        })
    });
    let tones: [_; 3] = std::array::from_fn(|i| jones_tone(F0 + i as f64, &[], &independent[i]));
    let independent = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            DISTANCE,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe_coherency([0.0; 2], || false)
        .unwrap();
    circular_matrix(
        independent.fluctuation_coherency_density_w,
        6.0 * amplitude * amplitude,
        2e-12,
    );
    assert_eq!(independent.mean_coherency_density_w.stokes().i, 0.0);
}

#[test]
fn g1_weighted_subnormal_jones_power_survives_negligible_temporal_modes() {
    let camera = focused();
    let h: [f64; 2] = std::array::from_fn(|i| {
        let wavelength =
            (Dd::from_f64(C) / Dd::from_f64(F0 + i as f64) * Dd::from_f64(1e9)).to_f64();
        camera
            .in_focus_diffraction_amplitude([0.0; 2], wavelength)
            .unwrap()
            .re
    });
    let scale = Dd::from_f64(4e-160) / (Dd::from_f64(h[0]) * Dd::from_f64(h[1]));
    let waves: [[CoherentJonesPupilWave; 2]; 2] = std::array::from_fn(|i| {
        let amplitude = Dd::from_f64(h[1 - i]) * scale * Dd::from_f64([1.0, -1.0][i]);
        [amplitude.hi, amplitude.lo]
            .map(|part| jones_wave([0.0; 2], [C64::from_re(part), C64::new(0.0, part)]))
    });
    let tones = [
        jones_tone(F0, &waves[0], &[]),
        jones_tone(F0 + 1.0, &waves[1], &[]),
        // An explicitly empty third line selects the multi-line band without
        // changing the physical two-source waveform used by the oracle.
        jones_tone(F0 + 2.0, &[], &[]),
    ];
    let value = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            DISTANCE,
            shutter(0.01),
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe_coherency([0.0; 2], || false)
        .unwrap()
        .coherency_density_w;
    // Exact two-line mean = 2*A^2*(1-sinc(pi*T)). Its stable Taylor
    // polynomial is independent of the temporal projection and positive form.
    let x2 = (core::f64::consts::PI * 0.01).powi(2);
    let expected = 2.0 * 4e-160_f64.powi(2) * (x2 / 6.0 - x2 * x2 / 120.0 + x2.powi(3) / 5040.0);
    assert!(expected > 0.0 && expected < f64::MIN_POSITIVE);
    for actual in value.axis_powers().into_iter().chain([-value.cross().im]) {
        assert!(actual > 0.0);
        assert!(
            (actual - expected).abs() <= 3.0 * f64::from_bits(1),
            "{actual} != {expected}"
        );
    }
    assert_eq!(value.cross().re, 0.0);
}

#[test]
fn g0_g4_short_jones_coherency_counts_projection_work_and_cancels_transactionally() {
    let camera = focused();
    let waves = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1.0)])];
    let tones: [_; 3] = std::array::from_fn(|i| jones_tone(F0 + i as f64, &waves, &[]));
    let shutter = shutter(1e-12);
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            DISTANCE,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let footprint = pixel([-0.1; 2], [0.1; 2]);
    let baseline = query
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let per_bank = baseline.pair_contractions / baseline.average_power.quadrature_samples;
    assert!(per_bank > 0);
    let limited = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            DISTANCE,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: per_bank - 1,
                ..EXPOSURE_LIMITS
            },
            || false,
        )
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        limited
            .observe_coherency([0.0; 2], || {
                polls += 1;
                false
            })
            .unwrap_err(),
        CoherentDiffractionError::Limit {
            quantity: "Jones coherency pair contractions",
            requested: per_bank,
            maximum: per_bank - 1,
        }
    );
    assert_eq!(polls, 1);
    let exact = CoherentDiffractionExposurePixelLimits {
        maximum_pair_contractions: baseline.pair_contractions,
        ..EXPOSURE_PIXEL_LIMITS
    };
    assert_eq!(
        query
            .integrate_pixel_coherency(footprint, exact, || false)
            .unwrap(),
        baseline
    );
    let mut polls = 0;
    assert_eq!(
        query
            .integrate_pixel_coherency(
                footprint,
                CoherentDiffractionExposurePixelLimits {
                    maximum_pair_contractions: baseline.pair_contractions - 1,
                    ..exact
                },
                || {
                    polls += 1;
                    false
                }
            )
            .unwrap_err(),
        CoherentDiffractionError::Limit {
            quantity: "pixel exposure pair contractions",
            requested: baseline.pair_contractions,
            maximum: baseline.pair_contractions - 1,
        }
    );
    assert_eq!(polls, 1);
    let mut polls = 0;
    let original = query
        .observe_coherency([0.0; 2], || {
            polls += 1;
            false
        })
        .unwrap();
    let late = polls - 1;
    assert!(late > 1);
    let mut polls = 0;
    assert_eq!(
        query.observe_coherency([0.0; 2], || {
            polls += 1;
            polls == late
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe_coherency([0.0; 2], || false).unwrap(),
        original
    );
}

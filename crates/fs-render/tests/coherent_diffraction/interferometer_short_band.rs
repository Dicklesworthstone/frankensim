//! Multiline shutter nulls retain physical two-arm moments and honest work.
use super::*;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits, CoherentJonesArm,
};

const C: f64 = 299_792_458.0;
const DISTANCE: f64 = 0.03;
const FREQUENCIES: [f64; 3] = [5.0e14, 5.0e14 + 1.0, 5.0e14 + 2.0];

fn receiver(fractions: [f64; 2], phase: f64) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], DISTANCE, fractions, 0.0, phase)
        .unwrap()
}

fn window(duration: f64) -> CoherentDiffractionShutter {
    // Source phases are at the pupil, so center the shutter after propagation.
    CoherentDiffractionShutter::try_new(DISTANCE / C - 0.5 * duration, duration).unwrap()
}

fn work(shutter: CoherentDiffractionShutter) -> [usize; 3] {
    let band = HarmonicWindow::new(shutter.start_s(), shutter.duration_s())
        .unwrap()
        .prepare_uniform_band_hz(&FREQUENCIES)
        .unwrap()
        .unwrap();
    [band.power_work(), band.cross_work(), band.coherency_work()]
}

fn limits(contractions: usize) -> CoherentInterferometerExposureLimits {
    CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 12,
        maximum_pair_contractions: contractions,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        expected != 0.0,
        "use an exact assertion for a physical zero"
    );
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn coefficients() -> [f64; 3] {
    // The on-axis pupil amplitude is proportional to frequency. Compensate
    // it so the arriving envelope approaches 1 - 2 exp(-i t) + exp(-2i t).
    [
        1.0,
        -2.0 * FREQUENCIES[0] / FREQUENCIES[1],
        FREQUENCIES[0] / FREQUENCIES[2],
    ]
}

#[test]
fn g3_short_three_line_scalar_interferometers_preserve_nulls_and_passive_energy() {
    let camera = focused();
    let amplitudes = coefficients();
    let mean = amplitudes.map(|a| [wave([0.0; 2], C64::new(a, 0.0))]);
    // One shared covariance column is coherent across frequencies, and is
    // independent of the deterministic mean. Its half amplitude adds 1/4 power.
    let random = amplitudes.map(|a| [wave([0.0; 2], C64::new(0.0, 0.5 * a))]);
    let factors = std::array::from_fn::<_, 3, _>(|i| [PupilCoherenceFactor { waves: &random[i] }]);
    let tones = std::array::from_fn::<_, 3, _>(|i| tone(FREQUENCIES[i], &mean[i], &factors[i]));
    let single_line = core::f64::consts::PI * RADIUS * RADIUS * (FREQUENCIES[0] / C).powi(2);

    for duration in [1e-7, 1e-9, 1e-12] {
        let shutter = window(duration);
        let [scalar, cross, _] = work(shutter);
        let both_work = 2 * (4 * scalar + cross);
        let query = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                DISTANCE,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let total = query.observe([0.0; 2], || false).unwrap().power_density_w;
        // Below binary64 bright-power subtraction accuracy, including at
        // durations whose ideal T^4 null is below authored f64 coefficient error.
        assert!(total > 0.0 && total < 1e-24 * single_line, "{total:.16e}");
        if duration == 1e-7 {
            // Integral |1-2e^-iwt+e^-2iwt|^2 = (wT)^4/80 + O((wT)^6).
            // Pupil wavelength and authored coefficients are binary64; the
            // shorter windows intentionally test physical invariants instead.
            let leading = 1.25 * single_line * (core::f64::consts::TAU * duration).powi(4) / 80.0;
            relative(total, leading, 0.02);
        }
        let matched = receiver([0.5; 2], 0.0);
        let sample = query
            .observe_interferometer(matched, limits(both_work), || false)
            .unwrap();
        relative(sample.output_power_density_w[0], total, 3e-11);
        assert_eq!(sample.output_power_density_w[1], 0.0);
        assert_eq!(sample.output_energy_density_j[1], 0.0);
        for arm in 0..2 {
            relative(sample.arm_power_density_w[arm], 0.5 * total, 3e-11);
            relative(
                sample.arm_energy_density_j[arm],
                0.5 * total * duration,
                3e-11,
            );
        }
        relative(sample.mutual_coherence_density_w.re, 0.5 * total, 3e-11);
        assert_eq!(sample.mutual_coherence_density_w.im, 0.0);
        relative(sample.fringe_visibility, 1.0, 3e-11);
        relative(sample.output_energy_density_j[0], total * duration, 3e-11);
        assert_eq!(sample.pair_contractions, both_work);
        assert_eq!(sample.pupil_evaluations, 12);
        assert!(matches!(
            query.observe_interferometer(matched, limits(both_work - 1), || false),
            Err(CoherentDiffractionError::Limit { .. })
        ));

        let one = query
            .observe_interferometer(receiver([0.5, 0.0], 0.0), limits(2 * scalar), || false)
            .unwrap();
        assert_eq!(one.pair_contractions, 2 * scalar);
        assert_eq!(one.pupil_evaluations, 6);
        assert_eq!(one.arm_power_density_w[1], 0.0);
        assert_eq!(one.fringe_visibility, 0.0);
        for power in one.output_power_density_w {
            relative(power, 0.25 * total, 3e-11);
        }

        if duration == 1e-9 {
            let phase = 0.37_f64;
            let shifted = query
                .observe_interferometer(receiver([0.5; 2], phase), limits(both_work), || false)
                .unwrap();
            relative(
                shifted.mutual_coherence_density_w.re,
                0.5 * total * phase.cos(),
                3e-11,
            );
            relative(
                shifted.mutual_coherence_density_w.im,
                -0.5 * total * phase.sin(),
                3e-11,
            );
            relative(shifted.fringe_visibility, 1.0, 3e-11);
            relative(shifted.output_power_density_w.iter().sum(), total, 3e-11);

            let offsets = pixel([-0.02, -0.015], [0.02, 0.015]);
            let collected = query
                .integrate_interferometer_pixel(matched, offsets, EXPOSURE_PIXEL_LIMITS, || false)
                .unwrap();
            let power = collected.average_power;
            assert!(power.output_power_w[0] > 0.0 && power.output_power_w[0].is_finite());
            assert_eq!(power.output_power_w[1], 0.0);
            relative(
                power.output_power_w[0],
                power.arm_power_w.iter().sum(),
                3e-11,
            );
            relative(
                power.mutual_coherence_w.re,
                0.5 * power.output_power_w[0],
                3e-11,
            );
            relative(power.fringe_visibility, 1.0, 3e-11);
            relative(
                collected.output_energy_j[0],
                power.output_power_w[0] * duration,
                3e-11,
            );
            assert_eq!(
                collected.pair_contractions,
                both_work * power.quadrature_samples
            );
            let exact = CoherentDiffractionExposurePixelLimits {
                maximum_pair_contractions: collected.pair_contractions,
                ..EXPOSURE_PIXEL_LIMITS
            };
            assert_eq!(
                query
                    .integrate_interferometer_pixel(matched, offsets, exact, || false)
                    .unwrap(),
                collected
            );
            assert!(matches!(
                query.integrate_interferometer_pixel(
                    matched,
                    offsets,
                    CoherentDiffractionExposurePixelLimits {
                        maximum_pair_contractions: collected.pair_contractions - 1,
                        ..exact
                    },
                    || false,
                ),
                Err(CoherentDiffractionError::Limit { .. })
            ));
        }
    }
}

#[test]
fn g3_short_three_line_jones_interferometers_keep_full_coherency_physical() {
    let camera = focused();
    let amplitudes = coefficients();
    let mean = amplitudes.map(|a| [jones_wave([0.0; 2], [C64::new(a, 0.0), C64::new(0.0, a)])]);
    // Opposite circular polarization in the independent covariance column
    // gives an explicitly mixed coherency: equal populations, V/I = -3/5.
    let random = amplitudes.map(|a| {
        [jones_wave(
            [0.0; 2],
            [C64::new(0.5 * a, 0.0), C64::new(0.0, -0.5 * a)],
        )]
    });
    let factors =
        std::array::from_fn::<_, 3, _>(|i| [JonesPupilCoherenceFactor { waves: &random[i] }]);
    let tones =
        std::array::from_fn::<_, 3, _>(|i| jones_tone(FREQUENCIES[i], &mean[i], &factors[i]));
    let clear = [CoherentJonesArm::IDENTITY; 2];

    for duration in [1e-7, 1e-9, 1e-12] {
        let shutter = window(duration);
        let [scalar, cross, coherency] = work(shutter);
        let both_work = 2 * (4 * scalar + 2 * coherency + 2 * cross);
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &tones,
                DISTANCE,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let reference = query
            .observe_coherency([0.0; 2], || false)
            .unwrap()
            .coherency_density_w;
        let total = reference.stokes().i;
        assert!(total > 0.0 && total.is_finite());
        let matched = receiver([0.5; 2], 0.0);
        let sample = query
            .observe_interferometer(matched, clear, limits(both_work), || false)
            .unwrap();
        let bright = sample.output_coherency_density_w[0];
        let dark = sample.output_coherency_density_w[1];
        relative(bright.stokes().i, total, 3e-11);
        assert_eq!(dark.axis_powers(), [0.0; 2]);
        assert_eq!(dark.cross(), C64::ZERO);
        for power in bright.axis_powers() {
            relative(power, 0.5 * total, 3e-11);
        }
        assert!(bright.cross().re.abs() <= 3e-11 * total);
        relative(bright.cross().im, -0.3 * total, 3e-11);
        for arm in 0..2 {
            relative(sample.arm_power_density_w[arm], 0.5 * total, 3e-11);
            relative(
                sample.arm_energy_density_j[arm],
                0.5 * total * duration,
                3e-11,
            );
        }
        relative(sample.mutual_coherence_density_w.re, 0.5 * total, 3e-11);
        assert_eq!(sample.mutual_coherence_density_w.im, 0.0);
        relative(sample.fringe_visibility, 1.0, 3e-11);
        relative(
            sample.output_coherency_energy_density_j[0].stokes().i,
            total * duration,
            3e-11,
        );
        assert_eq!(sample.pair_contractions, both_work);
        assert_eq!(sample.pupil_evaluations, 12);
        assert!(matches!(
            query.observe_interferometer(matched, clear, limits(both_work - 1), || false),
            Err(CoherentDiffractionError::Limit { .. })
        ));

        let one = query
            .observe_interferometer(
                receiver([0.5, 0.0], 0.0),
                clear,
                limits(2 * coherency),
                || false,
            )
            .unwrap();
        assert_eq!(one.pair_contractions, 2 * coherency);
        assert_eq!(one.pupil_evaluations, 6);
        assert_eq!(one.fringe_visibility, 0.0);
        for port in one.output_coherency_density_w {
            relative(port.stokes().i, 0.25 * total, 3e-11);
            relative(port.cross().im, -0.075 * total, 3e-11);
        }

        if duration == 1e-9 {
            let offsets = pixel([-0.02, -0.015], [0.02, 0.015]);
            let collected = query
                .integrate_interferometer_pixel(
                    matched,
                    clear,
                    offsets,
                    EXPOSURE_PIXEL_LIMITS,
                    || false,
                )
                .unwrap();
            let bright = collected.output_coherency_w[0];
            let dark = collected.output_coherency_w[1];
            assert!(bright.stokes().i > 0.0 && bright.stokes().i.is_finite());
            assert_eq!(dark.axis_powers(), [0.0; 2]);
            assert_eq!(dark.cross(), C64::ZERO);
            relative(bright.stokes().i, collected.arm_power_w.iter().sum(), 3e-11);
            relative(bright.cross().im, -0.3 * bright.stokes().i, 3e-11);
            relative(
                collected.mutual_coherence_w.re,
                0.5 * bright.stokes().i,
                3e-11,
            );
            relative(collected.fringe_visibility, 1.0, 3e-11);
            relative(
                collected.output_coherency_j[0].stokes().i,
                bright.stokes().i * duration,
                3e-11,
            );
            assert_eq!(
                collected.pair_contractions,
                both_work * collected.quadrature_samples
            );
            let exact = CoherentDiffractionExposurePixelLimits {
                maximum_pair_contractions: collected.pair_contractions,
                ..EXPOSURE_PIXEL_LIMITS
            };
            assert_eq!(
                query
                    .integrate_interferometer_pixel(matched, clear, offsets, exact, || false)
                    .unwrap(),
                collected
            );
            assert!(matches!(
                query.integrate_interferometer_pixel(
                    matched,
                    clear,
                    offsets,
                    CoherentDiffractionExposurePixelLimits {
                        maximum_pair_contractions: collected.pair_contractions - 1,
                        ..exact
                    },
                    || false,
                ),
                Err(CoherentDiffractionError::Limit { .. })
            ));
        }
    }
}

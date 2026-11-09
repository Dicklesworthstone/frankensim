//! G1 physical shutter throughput and independent polarization-beat integration;
//! G3 coherent nulls and retained source clocks; G0/G4 bounded refusal/reuse.
use super::*;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_math::dd::Dd;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits, CoherentJonesArm,
};

const C: f64 = 299_792_458.0;
const SENSOR: f64 = 0.03;
const FREQUENCIES: [f64; 3] = [5e14, 5e14 + 1.0, 5e14 + 2.0];
const POINT_LIMITS: CoherentInterferometerExposureLimits = CoherentInterferometerExposureLimits {
    maximum_pupil_evaluations: 1024,
    maximum_pair_contractions: 100_000,
};

fn triangular(start: Dd, duration: f64) -> CoherentDiffractionShutter {
    CoherentDiffractionShutter::try_with_profile_from_clock_dd(
        start,
        duration,
        HarmonicWindowProfile::Triangular,
    )
    .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual > 0.0 && (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_unit_peak_shutter_mass_reaches_every_scalar_jones_and_two_port_energy() {
    let camera = focused();
    let duration = 0.125;
    let shutter = triangular(Dd::from_f64(-0.2), duration);
    assert_eq!(shutter.profile(), HarmonicWindowProfile::Triangular);
    assert_eq!(shutter.duration_s(), duration);
    assert_eq!(shutter.effective_duration_s(), duration * 0.5);
    let footprint = pixel([-0.2; 2], [0.2; 2]);
    let receiver =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], SENSOR, [0.4, 0.6], 0.0, 0.37)
            .unwrap();
    let mean = [wave([0.0; 2], C64::ONE)];
    let noise = [wave([0.0; 2], C64::new(0.0, 0.5))];
    let factors = [PupilCoherenceFactor { waves: &noise }];
    let tones = [tone(FREQUENCIES[0], &mean, &factors)];
    let scalar = camera
        .prepare_coherent_diffraction_exposure(&tones, SENSOR, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let collected = scalar
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    for (power, energy) in [
        (
            collected.average_power.mean_power_w,
            collected.mean_energy_j,
        ),
        (
            collected.average_power.fluctuation_power_w,
            collected.fluctuation_energy_j,
        ),
        (collected.average_power.power_w, collected.energy_j),
    ] {
        relative(energy, power * duration * 0.5, 3e-14);
    }
    let point = scalar
        .observe_interferometer(receiver, POINT_LIMITS, || false)
        .unwrap();
    let area = scalar
        .integrate_interferometer_pixel(receiver, footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    for i in 0..2 {
        relative(
            point.arm_energy_density_j[i],
            point.arm_power_density_w[i] * duration * 0.5,
            3e-14,
        );
        relative(
            point.output_energy_density_j[i],
            point.output_power_density_w[i] * duration * 0.5,
            3e-14,
        );
        relative(
            area.arm_energy_j[i],
            area.average_power.arm_power_w[i] * duration * 0.5,
            3e-14,
        );
        relative(
            area.output_energy_j[i],
            area.average_power.output_power_w[i] * duration * 0.5,
            3e-14,
        );
    }
    relative(
        area.output_energy_j.iter().sum(),
        area.arm_energy_j.iter().sum(),
        3e-14,
    );

    let jones = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 0.5)])];
    let tones = [jones_tone(FREQUENCIES[0], &jones, &[])];
    let vector = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            SENSOR,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let analyzed = vector
        .integrate_pixel(footprint, JONES_TOTAL, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    relative(
        analyzed.energy_j,
        analyzed.average_power.power_w * duration * 0.5,
        3e-14,
    );
    let coherency = vector
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let scale = coherency.average_power.coherency_w.stokes().i;
    let p = coherency.average_power.coherency_w.stokes();
    let e = coherency.coherency_j.stokes();
    for (p, e) in [p.i, p.q, p.u, p.v].into_iter().zip([e.i, e.q, e.u, e.v]) {
        assert!((e - p * duration * 0.5).abs() < 3e-14 * scale * duration);
    }
    let point = vector
        .observe_interferometer(
            receiver,
            [CoherentJonesArm::IDENTITY; 2],
            POINT_LIMITS,
            || false,
        )
        .unwrap();
    let area = vector
        .integrate_interferometer_pixel(
            receiver,
            [CoherentJonesArm::IDENTITY; 2],
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    for i in 0..2 {
        relative(
            point.arm_energy_density_j[i],
            point.arm_power_density_w[i] * duration * 0.5,
            3e-14,
        );
        relative(
            point.output_coherency_energy_density_j[i].stokes().i,
            point.output_coherency_density_w[i].stokes().i * duration * 0.5,
            3e-14,
        );
        relative(
            area.arm_energy_j[i],
            area.arm_power_w[i] * duration * 0.5,
            3e-14,
        );
        relative(
            area.output_coherency_j[i].stokes().i,
            area.output_coherency_w[i].stokes().i * duration * 0.5,
            3e-14,
        );
    }
}

// Two authored pupil waves retain the DD inverse chromatic on-axis transfer.
// This prevents rounded source normalization from masking the physical null.
fn normalized(camera: &PhysicalCamera, f: f64, coefficient: f64) -> [CoherentPupilWave; 2] {
    let wavelength = (Dd::from_f64(C) / Dd::from_f64(f) * Dd::from_f64(1e9)).to_f64();
    let unit = [wave([0.0; 2], C64::ONE)];
    let transfer = camera
        .prepare_coherent_diffraction(
            wavelength,
            CoherentPupilIllumination {
                mean: &unit,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap()
        .mean_field;
    assert_eq!(transfer.im, 0.0);
    let amplitude = Dd::from_f64(coefficient) / Dd::from_f64(transfer.re);
    [
        wave([0.0; 2], C64::new(amplitude.hi, 0.0)),
        wave([0.0; 2], C64::new(amplitude.lo, 0.0)),
    ]
}

#[test]
fn g1_triangular_polarization_beats_match_independent_real_time_quadrature() {
    let camera = focused();
    let frequencies = [FREQUENCIES[0], FREQUENCIES[0] + 1e9];
    let delta = frequencies[1] - frequencies[0];
    let start = SENSOR / C - 0.13 / delta;
    let duration = 0.71 / delta;
    let shutter = CoherentDiffractionShutter::try_with_profile(
        start,
        duration,
        HarmonicWindowProfile::Triangular,
    )
    .unwrap();
    let scalar: [_; 2] = std::array::from_fn(|i| normalized(&camera, frequencies[i], 1.0));
    let waves: [[CoherentJonesPupilWave; 2]; 2] = std::array::from_fn(|i| {
        scalar[i].map(|w| {
            let mut amplitude = [C64::ZERO; 2];
            amplitude[i] = w.amplitude_sqrt_w;
            CoherentJonesPupilWave {
                tangent_direction: w.tangent_direction,
                amplitude_sqrt_w: amplitude,
            }
        })
    });
    let tones = [
        jones_tone(frequencies[0], &waves[0], &[]),
        jones_tone(frequencies[1], &waves[1], &[]),
    ];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            SENSOR,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let actual = query
        .observe_coherency([0.0; 2], || false)
        .unwrap()
        .coherency_density_w;
    let mut cross = C64::ZERO;
    let n = 4096;
    for i in 0..=n {
        let u = i as f64 / n as f64;
        let weight = (if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }) * 2.0
            * (1.0 - (2.0 * u - 1.0).abs())
            / (3.0 * n as f64);
        let phase = core::f64::consts::TAU * delta * (start + duration * u - SENSOR / C);
        cross = cross + C64::new(phase.cos(), phase.sin()).scale(weight);
    }
    close(actual.axis_powers()[0], 1.0, 3e-13);
    close(actual.axis_powers()[1], 1.0, 3e-13);
    close(actual.cross().re, cross.re, 3e-12);
    close(actual.cross().im, cross.im, 3e-12);
    let plus = CoherentJonesDetector::analyzer([C64::ONE, C64::ONE]).unwrap();
    close(
        query
            .observe([0.0; 2], plus, || false)
            .unwrap()
            .power_density_w,
        1.0 + cross.re,
        3e-12,
    );
    let expected_visibility =
        sinc_reference(core::f64::consts::PI * delta * duration * 0.5).powi(2);
    close(actual.cross().abs(), expected_visibility, 3e-12);
}

#[test]
fn g1_g3_positive_triangle_modes_retain_fourth_order_scalar_and_jones_nulls() {
    let camera = focused();
    let waves: [_; 3] =
        std::array::from_fn(|i| normalized(&camera, FREQUENCIES[i], [1.0, -2.0, 1.0][i]));
    let tones: [_; 3] = std::array::from_fn(|i| tone(FREQUENCIES[i], &waves[i], &[]));
    let vector = waves.map(|bank| {
        bank.map(|w| CoherentJonesPupilWave {
            tangent_direction: w.tangent_direction,
            amplitude_sqrt_w: [w.amplitude_sqrt_w, w.amplitude_sqrt_w * C64::new(0.0, 0.5)],
        })
    });
    let jones: [_; 3] = std::array::from_fn(|i| jones_tone(FREQUENCIES[i], &vector[i], &[]));
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let start =
            Dd::from_f64(SENSOR) / Dd::from_f64(C) - Dd::from_f64(duration) * Dd::from_f64(0.5);
        let shutter = triangular(start, duration);
        let query = camera
            .prepare_coherent_diffraction_exposure(&tones, SENSOR, shutter, EXPOSURE_LIMITS, || {
                false
            })
            .unwrap();
        // |1-2 exp(-i 2 pi t)+exp(-i 4 pi t)|² = 16 sin(pi t)^4.
        // The normalized symmetric triangle has E[t^4]=T^4/240.
        let expected = core::f64::consts::PI.powi(4) * duration.powi(4) / 15.0;
        relative(
            query.observe([0.0; 2], || false).unwrap().power_density_w,
            expected,
            3e-6,
        );
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &jones,
                SENSOR,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let actual = query
            .observe_coherency([0.0; 2], || false)
            .unwrap()
            .coherency_density_w;
        relative(actual.stokes().i, 1.25 * expected, 3e-6);
        relative(-actual.stokes().v, expected, 3e-6);
    }
}

#[test]
fn g3_triangular_shutter_retains_the_sub_ulp_entrance_pupil_epoch() {
    let camera = focused();
    let waves: [_; 3] =
        std::array::from_fn(|i| normalized(&camera, FREQUENCIES[i], [1.0, -2.0, 1.0][i]));
    let tones: [_; 3] = std::array::from_fn(|i| tone(FREQUENCIES[i], &waves[i], &[]));
    let epoch = Dd::from_pair(2.0_f64.powi(49), 0.017);
    let observe = |start| {
        let shutter = triangular(start, 0.02);
        assert_eq!(shutter.start_clock_dd(), start);
        camera
            .prepare_coherent_diffraction_exposure(&tones, SENSOR, shutter, EXPOSURE_LIMITS, || {
                false
            })
            .unwrap()
            .observe([0.0; 2], || false)
            .unwrap()
            .power_density_w
    };
    let precise = observe(epoch);
    relative(precise, observe(Dd::from_f64(0.017)), 3e-13);
    let rounded = observe(Dd::from_f64(epoch.to_f64()));
    assert!((rounded - precise).abs() > 0.1 * precise);
}

#[test]
fn g0_g4_profile_mass_work_limits_and_cancelled_reuse_are_explicit() {
    assert_eq!(
        CoherentDiffractionShutter::try_new(0.0, 0.25)
            .unwrap()
            .effective_duration_s(),
        0.25
    );
    assert!(
        CoherentDiffractionShutter::try_with_profile(
            0.0,
            f64::from_bits(1),
            HarmonicWindowProfile::Triangular
        )
        .is_err()
    );
    assert!(
        CoherentDiffractionShutter::try_with_profile(0.0, 0.0, HarmonicWindowProfile::Triangular)
            .is_err()
    );
    let camera = focused();
    let waves: [_; 3] = std::array::from_fn(|i| normalized(&camera, FREQUENCIES[i], 1.0));
    let tones: [_; 3] = std::array::from_fn(|i| tone(FREQUENCIES[i], &waves[i], &[]));
    let shutter = triangular(Dd::ZERO, 1e-9);
    let modes = HarmonicWindow::with_profile(0.0, 1e-9, HarmonicWindowProfile::Triangular)
        .unwrap()
        .band_mode_count_hz(&FREQUENCIES)
        .unwrap()
        .unwrap();
    let work = modes * 4;
    assert!(matches!(
        camera.prepare_coherent_diffraction_exposure(
            &tones,
            SENSOR,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: work - 1,
                ..EXPOSURE_LIMITS
            },
            || false
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "exposure pair contractions",
            ..
        })
    ));
    let query = camera
        .prepare_coherent_diffraction_exposure(
            &tones,
            SENSOR,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: work,
                ..EXPOSURE_LIMITS
            },
            || false,
        )
        .unwrap();
    let baseline = query.observe([0.0; 2], || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        query.observe([0.0; 2], || {
            polls += 1;
            polls == 8
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(query.observe([0.0; 2], || false).unwrap(), baseline);
    let sample = query
        .integrate_pixel(pixel([-0.1; 2], [0.1; 2]), EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    assert_eq!(
        sample.pair_contractions,
        work * sample.average_power.quadrature_samples
    );
    relative(
        sample.energy_j,
        sample.average_power.power_w * 0.5e-9,
        3e-14,
    );
}

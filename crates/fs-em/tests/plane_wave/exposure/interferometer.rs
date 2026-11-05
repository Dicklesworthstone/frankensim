//! Finite-shutter multi-arm interference through the original Maxwell owners.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureInterferenceArm as Arm, AngularSpectrumFactorStatistics as Statistics,
    AngularSpectrumIlluminationLimits, AngularSpectrumIlluminationSide as Face,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn phase(value: f64) -> C64 {
    let (s, c) = value.sin_cos();
    C64::new(c, s)
}

fn arm(x: f64) -> Arm {
    Arm {
        side: Face::OutgoingRight,
        distance_from_face_m: 0.0,
        position_m: [x, 0.0],
        analyzer: [C64::ZERO, C64::ONE, C64::ZERO],
    }
}

// Independent elementary sum/difference-frequency integral. The production
// kernel's DD centered Gram calculation is not used by this oracle.
fn mean_phase(omega: f64, window: HarmonicWindow) -> C64 {
    let x = omega * window.duration_s() * 0.5;
    let sinc = |x: f64| if x == 0.0 { 1.0 } else { x.sin() / x };
    let weight = match window.profile() {
        HarmonicWindowProfile::Uniform => sinc(x),
        HarmonicWindowProfile::Triangular => sinc(x * 0.5).powi(2),
    };
    phase(omega * (window.start_s() + 0.5 * window.duration_s())).scale(weight)
}

fn reference(omega: &[f64], amplitudes: &[C64], window: HarmonicWindow, proper: bool) -> f64 {
    let mut result = 0.0;
    for (i, &a) in amplitudes.iter().enumerate() {
        for (j, &b) in amplitudes.iter().enumerate() {
            result += 0.5 * (a * b.conj() * mean_phase(omega[i] - omega[j], window)).re;
            if !proper {
                result += 0.5 * (a * b * mean_phase(omega[i] + omega[j], window)).re;
            }
        }
    }
    result
}

fn delayed(omega: &[f64], source: &[C64], gains: &[C64], delays: &[f64]) -> Vec<C64> {
    omega
        .iter()
        .zip(source)
        .map(|(&w, &amplitude)| {
            gains
                .iter()
                .zip(delays)
                .fold(C64::ZERO, |sum, (&gain, &delay)| {
                    sum + gain * amplitude * phase(-w * delay)
                })
        })
        .collect()
}

#[test]
fn g1_broadband_arm_delays_retain_real_beats_and_both_shutter_profiles() {
    let medium = material(1.0, 1.0);
    let omega = [2.0, 3.0];
    let stacks = omega.map(|w| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), w, &[])
            .unwrap()
    });
    let amplitudes = [phase(0.2), phase(-0.4).scale(0.7)];
    let lines = [
        constant(&stacks[0], [C64::ZERO, amplitudes[0]]),
        constant(&stacks[1], [C64::ZERO, amplitudes[1]]),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let arms = [arm(0.0), arm(1.0), arm(-2.0)];
    let gains = [C64::from_re(0.8), C64::new(0.0, -0.4), C64::ONE];
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.7, profile).unwrap();
        let mut prepared = bank
            .prepare_interferometer(
                &arms,
                window,
                bank.interferometer_requirements(arms.len()).unwrap(),
                || false,
            )
            .unwrap();
        assert_eq!(prepared.arm_count(), 3);
        assert_eq!(prepared.frequency_group_count(), 2);
        assert_eq!(prepared.window(), window);
        for step in 0..13 {
            let delays = [0.0, (step as f64 - 4.0) * 0.125, 0.37];
            let expected = reference(
                &omega,
                &delayed(&omega, &amplitudes, &gains, &delays),
                window,
                false,
            );
            let actual = prepared
                .electric_mean_square_v2_m2(&gains, &delays, || false)
                .unwrap();
            near(actual, expected, 3e-12);
            assert_eq!(
                actual.to_bits(),
                prepared
                    .electric_mean_square_v2_m2(&gains, &delays, || false)
                    .unwrap()
                    .to_bits()
            );
        }
    }
}

#[test]
fn g1_broadband_mixed_factors_keep_their_authored_joint_second_moments() {
    let medium = material(1.0, 1.0);
    let omega = [2.0, 3.0];
    let stacks = omega.map(|w| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), w, &[])
            .unwrap()
    });
    let mean_amplitudes = [C64::from_re(0.3), C64::ZERO];
    let real_amplitudes = [C64::ONE, phase(0.4).scale(0.5)];
    let proper_amplitudes = [C64::new(0.0, 0.2), C64::from_re(-0.4)];
    let mean_lines = [
        constant(&stacks[0], [C64::ZERO, mean_amplitudes[0]]),
        constant(&stacks[1], [C64::ZERO, mean_amplitudes[1]]),
    ];
    let real_lines = [
        constant(&stacks[0], [C64::ZERO, real_amplitudes[0]]),
        constant(&stacks[1], [C64::ZERO, real_amplitudes[1]]),
    ];
    let proper_lines = [
        constant(&stacks[0], [C64::ZERO, proper_amplitudes[0]]),
        constant(&stacks[1], [C64::ZERO, proper_amplitudes[1]]),
    ];
    let mean =
        AngularSpectrumExposure::try_new(&mean_lines[0], &mean_lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new(&mean_lines[0], &real_lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new(&mean_lines[0], &proper_lines, limits(), || false)
            .unwrap(),
    ];
    let statistics = [Statistics::Real, Statistics::ProperComplex];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &mean_lines[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window = HarmonicWindow::new(-0.07, 0.2).unwrap();
    let mut prepared = ensemble
        .prepare_interferometer(
            &[arm(0.0); 2],
            window,
            ensemble.interferometer_requirements(2).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(prepared.frequency_group_count(), 6);
    let gains = [C64::ONE, phase(0.3).scale(0.6)];
    for delays in [[0.0, 0.0], [0.2, 0.7], [-0.4, 1.0]] {
        let expected = [
            (mean_amplitudes, false),
            (real_amplitudes, false),
            (proper_amplitudes, true),
        ]
        .into_iter()
        .map(|(a, proper)| {
            reference(
                &omega,
                &delayed(&omega, &a, &gains, &delays),
                window,
                proper,
            )
        })
        .sum();
        near(
            prepared
                .electric_mean_square_v2_m2(&gains, &delays, || false)
                .unwrap(),
            expected,
            2e-12,
        );
    }
}

#[test]
fn g1_broadband_dark_port_preserves_original_sub_ulp_modal_amplitudes() {
    let medium = material(1.0, 1.0);
    let omega = 2.0;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let epsilon = 2.0_f64.powi(-80);
    let mut modes = [[C64::ZERO; 2]; 4];
    modes[0][1] = C64::ONE;
    modes[1][1] = C64::from_re(epsilon);
    let lines = [PlanarAngularSpectrum::prepare_modes(
        &stack,
        grid(),
        &modes,
        AngularSpectrumLimits {
            maximum_grid_points: 128,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let arms = [arm(0.0), arm(16.0), arm(0.0)];
    let gains = [C64::ONE, C64::ONE, C64::from_re(-2.0)];
    let duration = 0.1;
    let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
    let mut prepared = bank
        .prepare_interferometer(
            &arms,
            window,
            bank.interferometer_requirements(3).unwrap(),
            || false,
        )
        .unwrap();
    let expected =
        4.0 * epsilon * epsilon * (0.5 + 0.5 * (omega * duration).sin() / (omega * duration));
    near(
        prepared
            .electric_mean_square_v2_m2(&gains, &[0.0; 3], || false)
            .unwrap()
            / expected,
        1.0,
        2e-13,
    );
    // A proper source has no carrier sum tone, even for the same short shutter.
    let factors = [bank];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let mut prepared = ensemble
        .prepare_interferometer(
            &arms,
            window,
            ensemble.interferometer_requirements(3).unwrap(),
            || false,
        )
        .unwrap();
    near(
        prepared
            .electric_mean_square_v2_m2(&gains, &[0.0; 3], || false)
            .unwrap()
            / (2.0 * epsilon * epsilon),
        1.0,
        2e-13,
    );
}

#[test]
fn g1_broadband_interferometer_retains_both_native_illumination_faces() {
    let medium = material(1.0, 1.0);
    let omega = [2.0, 3.0];
    let stacks = omega.map(|w| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), w, &[])
            .unwrap()
    });
    let left = [C64::ONE, C64::new(0.0, 0.4)];
    let right = [C64::new(0.0, 0.3), C64::from_re(0.2)];
    let beams = [
        constant(&stacks[0], [C64::ZERO, left[0]]),
        constant(&stacks[1], [C64::ZERO, left[1]]),
    ];
    let illumination: [_; 2] = std::array::from_fn(|i| {
        beams[i]
            .prepare_illumination(
                &[[C64::ZERO, right[i]]; 4],
                AngularSpectrumIlluminationLimits {
                    maximum_layer_mode_pairs: 4096,
                    maximum_workspace_bytes: 8 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
    });
    let bank =
        AngularSpectrumExposure::try_new_illumination(&beams[0], &illumination, limits(), || false)
            .unwrap();
    let arms = [
        Arm {
            side: Face::IncomingLeft,
            ..arm(0.0)
        },
        Arm {
            side: Face::IncomingRight,
            ..arm(0.0)
        },
    ];
    let window = HarmonicWindow::new(0.1, 0.4).unwrap();
    let mut prepared = bank
        .prepare_interferometer(
            &arms,
            window,
            bank.interferometer_requirements(2).unwrap(),
            || false,
        )
        .unwrap();
    let gains = [phase(0.2), phase(-0.3).scale(0.7)];
    let delays = [0.0, 0.37];
    let amplitudes: Vec<_> = (0..2)
        .map(|i| left[i] * gains[0] + right[i] * gains[1] * phase(-omega[i] * delays[1]))
        .collect();
    near(
        prepared
            .electric_mean_square_v2_m2(&gains, &delays, || false)
            .unwrap(),
        reference(&omega, &amplitudes, window, false),
        2e-12,
    );
}

#[test]
fn g0_g4_broadband_interferometer_admits_budgets_and_keeps_dark_refusals() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ZERO, C64::ONE])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let arms = [arm(0.0); 2];
    let exact = bank.interferometer_requirements(2).unwrap();
    let mut short = [exact; 3];
    short[0].maximum_mode_evaluations -= 1;
    short[1].maximum_frequency_pairs -= 1;
    short[2].maximum_workspace_bytes -= 1;
    for limit in short {
        assert!(
            bank.prepare_interferometer(&arms, window, limit, || false)
                .is_err()
        );
    }
    assert!(bank.interferometer_requirements(0).is_err());
    assert!(bank.interferometer_requirements(4097).is_err());
    assert!(matches!(
        bank.prepare_interferometer(&arms, window, exact, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let mut prepared = bank
        .prepare_interferometer(&arms, window, exact, || false)
        .unwrap();
    assert!(matches!(
        prepared.electric_mean_square_v2_m2(&[C64::ONE; 2], &[0.0; 2], || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(
        prepared
            .electric_mean_square_v2_m2(&[C64::ONE; 2], &[1e308; 2], || false)
            .is_err()
    );
    assert_eq!(
        prepared
            .electric_mean_square_v2_m2(&[C64::ONE, -C64::ONE], &[0.0; 2], || false)
            .unwrap(),
        0.0
    );
    let dark = AngularSpectrumExposure::try_new(&lines[0], &[], limits(), || false).unwrap();
    let exact = dark.interferometer_requirements(2).unwrap();
    let mut prepared = dark
        .prepare_interferometer(&arms, window, exact, || false)
        .unwrap();
    assert_eq!(
        prepared
            .electric_mean_square_v2_m2(&[C64::ONE; 2], &[0.0; 2], || false)
            .unwrap(),
        0.0
    );
    assert!(
        prepared
            .electric_mean_square_v2_m2(&[C64::ONE], &[0.0; 2], || false)
            .is_err()
    );
    assert!(
        prepared
            .electric_mean_square_v2_m2(&[C64::ONE; 2], &[f64::NAN; 2], || false)
            .is_err()
    );
    assert!(
        prepared
            .electric_mean_square_v2_m2(&[C64::new(f64::NAN, 0.0); 2], &[0.0; 2], || false)
            .is_err()
    );
    let invalid = [
        Arm {
            position_m: [f64::INFINITY, 0.0],
            ..arms[0]
        },
        arms[1],
    ];
    assert!(
        dark.prepare_interferometer(&invalid, window, exact, || false)
            .is_err()
    );
}

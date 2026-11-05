//! Finite-aperture E/H interference through the original Maxwell owners.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposureApertureArm as Arm,
    AngularSpectrumExposureApertureMoments as Moments,
    AngularSpectrumFactorStatistics as Statistics, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn cell() -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 2,
        height: 2,
        spacing_m: [8.0, 9.0],
    }
}
fn phase(x: f64) -> C64 {
    let (s, c) = x.sin_cos();
    C64::new(c, s)
}
fn arm(depth: f64, shift: [f64; 2]) -> Arm {
    Arm {
        distance_from_face_m: depth,
        translation_m: shift,
    }
}
fn amplitudes(frequency: usize) -> [[C64; 2]; 4] {
    core::array::from_fn(|m| {
        [
            C64::new(0.08 * (m + 1) as f64, -0.03 * (frequency + 1) as f64),
            phase(0.19 * (m + frequency) as f64).scale(0.21 + 0.13 * m as f64),
        ]
    })
}
fn modes<'s, 'm>(
    stack: &'s PlanarDielectricLayers<'m>,
    values: &[[C64; 2]],
) -> PlanarAngularSpectrum<'s, 'm> {
    PlanarAngularSpectrum::prepare_modes(
        stack,
        cell(),
        values,
        AngularSpectrumLimits {
            maximum_grid_points: 16,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 16,
        },
        || false,
    )
    .unwrap()
}
fn wavevector(omega: f64, mode: usize) -> [f64; 3] {
    let kx = if mode % 2 == 0 {
        0.0
    } else {
        -core::f64::consts::TAU / 16.0
    };
    let ky = if mode / 2 == 0 {
        0.0
    } else {
        -core::f64::consts::TAU / 18.0
    };
    [kx, ky, (omega * omega - kx * kx - ky * ky).sqrt()]
}
// Independent homogeneous-medium Maxwell plane wave, epsilon=mu=1.
fn field(omega: f64, mode: usize, tangential: [C64; 2]) -> [C64; 6] {
    let [x, y, z] = wavevector(omega, mode);
    let [ex, ey] = tangential;
    let ez = (ex.scale(-x) + ey.scale(-y)).scale(1.0 / z);
    [
        ex,
        ey,
        ez,
        (ez.scale(y) - ey.scale(z)).scale(1.0 / omega),
        (ex.scale(z) - ez.scale(x)).scale(1.0 / omega),
        (ey.scale(x) - ex.scale(y)).scale(1.0 / omega),
    ]
}
fn mean_phase(frequency: f64, origin: f64, extent: f64, triangular: bool) -> C64 {
    let sinc = |x: f64| if x == 0.0 { 1.0 } else { x.sin() / x };
    let x = frequency * extent * 0.5;
    let weight = if triangular {
        sinc(x * 0.5).powi(2)
    } else {
        sinc(x)
    };
    phase(frequency * (origin + 0.5 * extent)).scale(weight)
}
#[allow(clippy::too_many_arguments)]
fn reference(
    frequencies: &[f64],
    arms: &[Arm],
    gains: &[C64],
    delays: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    proper: bool,
) -> [[f64; 6]; 6] {
    let mut waves = Vec::new();
    for (f, &omega) in frequencies.iter().enumerate() {
        for (mode, tangential) in amplitudes(f).into_iter().enumerate() {
            let [kx, ky, kz] = wavevector(omega, mode);
            let gain =
                arms.iter()
                    .zip(gains)
                    .zip(delays)
                    .fold(C64::ZERO, |sum, ((arm, gain), delay)| {
                        sum + *gain
                            * phase(
                                -kz * arm.distance_from_face_m
                                    - kx * arm.translation_m[0]
                                    - ky * arm.translation_m[1]
                                    - omega * delay,
                            )
                    });
            waves.push((
                [kx, ky, omega],
                field(omega, mode, tangential).map(|v| v * gain),
            ));
        }
    }
    let [x, y] = aperture.origin_m();
    let [width, height] = aperture.extent_m();
    let mut result = [[0.0; 6]; 6];
    for (a, left) in &waves {
        for (b, right) in &waves {
            for sign in [-1.0, 1.0] {
                if proper && sign > 0.0 {
                    continue;
                }
                let kernel = mean_phase(-(a[0] + sign * b[0]), x, width, false)
                    * mean_phase(-(a[1] + sign * b[1]), y, height, false)
                    * mean_phase(
                        a[2] + sign * b[2],
                        window.start_s(),
                        window.duration_s(),
                        window.profile() == HarmonicWindowProfile::Triangular,
                    );
                for i in 0..6 {
                    for j in 0..6 {
                        let other = if sign < 0.0 {
                            right[j].conj()
                        } else {
                            right[j]
                        };
                        result[i][j] += 0.5 * (left[i] * other * kernel).re;
                    }
                }
            }
        }
    }
    result
}
fn compare(actual: Moments, expected: [[f64; 6]; 6], aperture: AngularSpectrumAperture) {
    for (block, row, column) in [
        (actual.electric_second_moment_v2_m2(), 0, 0),
        (actual.magnetic_second_moment_a2_m2(), 3, 3),
        (actual.electric_magnetic_second_moment_w_m2(), 0, 3),
    ] {
        for i in 0..3 {
            for j in 0..3 {
                near(block[i][j], expected[row + i][column + j], 3e-11);
            }
        }
    }
    near(
        actual.reading().normal_power_w,
        (expected[0][4] - expected[1][3]) * aperture.area_m2(),
        8e-10,
    );
}

#[test]
fn g1_delayed_2d_aperture_preserves_all_eh_moments_and_independent_factor_laws() {
    let medium = material(1.0, 1.0);
    let frequencies = [2.0, 3.0];
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let lines = [
        modes(&stacks[0], &amplitudes(0)),
        modes(&stacks[1], &amplitudes(1)),
    ];
    let mean = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new(&lines[0], &lines[..1], limits(), || false).unwrap(),
    ];
    let statistics = [Statistics::Real, Statistics::ProperComplex];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let aperture = AngularSpectrumAperture::try_new([0.2, -0.4], [7.3, 5.1]).unwrap();
    let arms = [
        arm(0.11, [0.0, 0.0]),
        arm(0.31, [1.7, -0.6]),
        arm(0.11, [-0.4, 0.9]),
    ];
    let gains = [C64::ONE, C64::new(0.2, -0.7), C64::from_re(-0.3)];
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.51, profile).unwrap();
        let mut coherent = mean
            .prepare_aperture_interferometer(
                Face::OutgoingRight,
                &arms,
                aperture,
                window,
                mean.aperture_interferometer_requirements(3).unwrap(),
                || false,
            )
            .unwrap();
        let mut mixed = ensemble
            .prepare_aperture_interferometer(
                Face::OutgoingRight,
                &arms,
                aperture,
                window,
                ensemble.aperture_interferometer_requirements(3).unwrap(),
                || false,
            )
            .unwrap();
        for delays in [[0.0, 0.0, 0.0], [0.1, 0.37, -0.2], [-0.8, 0.03, 0.4]] {
            let real = reference(
                &frequencies,
                &arms,
                &gains,
                &delays,
                aperture,
                window,
                false,
            );
            compare(
                coherent.measure_moments(&gains, &delays, || false).unwrap(),
                real,
                aperture,
            );
            let proper = reference(
                &frequencies[..1],
                &arms,
                &gains,
                &delays,
                aperture,
                window,
                true,
            );
            let expected =
                core::array::from_fn(|i| core::array::from_fn(|j| 2.0 * real[i][j] + proper[i][j]));
            let reading = mixed.measure_moments(&gains, &delays, || false).unwrap();
            compare(reading, expected, aperture);
            assert_eq!(
                reading,
                mixed.measure_moments(&gains, &delays, || false).unwrap()
            );
        }
    }
}

#[test]
fn g1_aperture_watts_retain_backward_flux_and_counterpropagating_field_cancellation() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let left = constant(&stack, [C64::ZERO, C64::ONE]);
    let illumination = [left
        .prepare_illumination(
            &[[C64::ZERO, C64::ONE]; 4],
            AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 128,
                maximum_workspace_bytes: 1 << 20,
                maximum_fft_work: 1 << 16,
            },
            || false,
        )
        .unwrap()];
    let bank =
        AngularSpectrumExposure::try_new_illumination(&left, &illumination, limits(), || false)
            .unwrap();
    let aperture = AngularSpectrumAperture::try_new([-0.2, 0.3], [2.0, 1.5]).unwrap();
    let window = HarmonicWindow::new(-0.1, 0.2).unwrap();
    let limits = bank.aperture_interferometer_requirements(2).unwrap();
    for face in [Face::IncomingRight, Face::TotalLeft] {
        let mut detector = bank
            .prepare_aperture_interferometer(
                face,
                &[arm(0.0, [0.0; 2]); 2],
                aperture,
                window,
                limits,
                || false,
            )
            .unwrap();
        let reading = detector
            .measure_moments(&[C64::from_re(0.6), C64::from_re(0.4)], &[0.0; 2], || false)
            .unwrap()
            .reading();
        assert!(reading.mean.electric_mean_squared_v2_m2[1] > 0.0);
        if face == Face::IncomingRight {
            near(
                reading.normal_power_w
                    / (-aperture.area_m2() * reading.mean.electric_mean_squared_v2_m2[1]),
                1.0,
                2e-13,
            );
        } else {
            near(reading.mean.magnetic_mean_squared_a2_m2[0], 0.0, 1e-28);
            near(reading.normal_power_w, 0.0, 1e-28);
        }
    }
}

#[test]
fn g1_aperture_dark_port_retains_sub_ulp_spatial_field_and_signed_power() {
    let medium = material(1.0, 1.0);
    let omega = 2.0;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let epsilon = 2.0_f64.powi(-80);
    let mut values = [[C64::ZERO; 2]; 4];
    values[0][1] = C64::ONE;
    values[1][1] = C64::from_re(epsilon);
    let lines = [modes(&stack, &values)];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let arms = [arm(0.0, [0.0; 2]), arm(0.0, [8.0, 0.0]), arm(0.0, [0.0; 2])];
    let aperture = AngularSpectrumAperture::try_new([0.3, -0.1], [3.7, 2.1]).unwrap();
    let mut detector = ensemble
        .prepare_aperture_interferometer(
            Face::OutgoingRight,
            &arms,
            aperture,
            HarmonicWindow::new(0.17, 0.11).unwrap(),
            ensemble.aperture_interferometer_requirements(3).unwrap(),
            || false,
        )
        .unwrap();
    let reading = detector
        .measure_moments(&[C64::ONE, C64::ONE, C64::from_re(-2.0)], &[0.0; 3], || {
            false
        })
        .unwrap()
        .reading();
    let expected = 2.0 * epsilon * epsilon;
    near(
        reading.mean.electric_mean_squared_v2_m2[1] / expected,
        1.0,
        3e-13,
    );
    near(
        reading.normal_power_w / (expected * wavevector(omega, 1)[2] / omega * aperture.area_m2()),
        1.0,
        3e-13,
    );
}

#[test]
fn g3_aperture_interference_enforces_budgets_dark_controls_and_cancellable_reuse() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [modes(&stack, &amplitudes(0))];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [2.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.1, 0.3).unwrap();
    let arms = [arm(0.0, [0.0; 2]); 2];
    let needed = bank.aperture_interferometer_requirements(2).unwrap();
    for which in 0..4 {
        let mut cap = needed;
        match which {
            0 => cap.maximum_mode_evaluations -= 1,
            1 => cap.maximum_source_field_entries -= 1,
            2 => cap.maximum_gram_products -= 1,
            _ => cap.maximum_workspace_bytes -= 1,
        }
        assert!(
            bank.prepare_aperture_interferometer(
                Face::OutgoingRight,
                &arms,
                aperture,
                window,
                cap,
                || false
            )
            .is_err()
        );
    }
    let mut prepared = bank
        .prepare_aperture_interferometer(
            Face::OutgoingRight,
            &arms,
            aperture,
            window,
            needed,
            || false,
        )
        .unwrap();
    let expected = prepared
        .measure_moments(&[C64::ONE; 2], &[0.0; 2], || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.measure_moments(&[C64::ONE; 2], &[0.0; 2], || {
            polls += 1;
            polls == 15
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        expected,
        prepared
            .measure_moments(&[C64::ONE; 2], &[0.0; 2], || false)
            .unwrap()
    );
    let dark = AngularSpectrumExposure::try_new(&lines[0], &[], limits(), || false).unwrap();
    let needed = dark.aperture_interferometer_requirements(2).unwrap();
    let mut prepared = dark
        .prepare_aperture_interferometer(
            Face::OutgoingRight,
            &arms,
            aperture,
            window,
            needed,
            || false,
        )
        .unwrap();
    assert_eq!(
        prepared
            .measure_moments(&[C64::ONE; 2], &[0.0; 2], || false)
            .unwrap()
            .reading()
            .normal_power_w,
        0.0
    );
    assert!(
        prepared
            .measure_moments(&[C64::ONE; 1], &[0.0; 2], || false)
            .is_err()
    );
    assert!(
        prepared
            .measure_moments(&[C64::ONE, C64::from_re(f64::NAN)], &[0.0; 2], || false)
            .is_err()
    );
    assert!(
        prepared
            .measure_moments(&[C64::ONE; 2], &[0.0, f64::INFINITY], || false)
            .is_err()
    );
    for invalid in [arm(-1.0, [0.0; 2]), arm(0.0, [f64::INFINITY, 0.0])] {
        assert!(
            dark.prepare_aperture_interferometer(
                Face::OutgoingRight,
                &[invalid; 2],
                aperture,
                window,
                needed,
                || false
            )
            .is_err()
        );
    }
}

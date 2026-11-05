//! Arbitrary-angle, finite-shutter fields against independent real Maxwell waves.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposureApertureLimits,
    AngularSpectrumFactorStatistics, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

#[derive(Clone, Copy)]
struct Wave {
    omega: f64,
    k: [f64; 3],
    field: [C64; 6],
}

fn rotation(angle: f64) -> C64 {
    let (s, c) = angle.sin_cos();
    C64::new(c, s)
}

// Independently impose k.E=0 and H=k cross E/(omega*mu), with mu=1.
fn wave(omega: f64, k: [f64; 3], transverse: [C64; 2]) -> Wave {
    let [ex, ey] = transverse;
    let ez = -(ex.scale(k[0]) + ey.scale(k[1])).scale(k[2].recip());
    Wave {
        omega,
        k,
        field: [
            ex,
            ey,
            ez,
            (ez.scale(k[1]) - ey.scale(k[2])).scale(omega.recip()),
            (ex.scale(k[2]) - ez.scale(k[0])).scale(omega.recip()),
            (ey.scale(k[0]) - ex.scale(k[1])).scale(omega.recip()),
        ],
    }
}

// Sixteen-point Gauss-Legendre quadrature. It integrates the actual real
// signals; no production harmonic kernels or modal products enter this oracle.
fn average(start: f64, extent: f64, mut f: impl FnMut(f64) -> f64) -> f64 {
    let mut sum = 0.0;
    for (node, weight) in [
        (0.09501250983763744, 0.1894506104550685),
        (0.2816035507792589, 0.1826034150449236),
        (0.4580167776572274, 0.16915651939500254),
        (0.6178762444026438, 0.14959598881657673),
        (0.755404408355003, 0.12462897125553388),
        (0.8656312023878318, 0.09515851168249279),
        (0.9445750230732326, 0.06225352393864789),
        (0.9894009349916499, 0.027152459411754095),
    ] {
        for sign in [-1.0, 1.0] {
            sum += 0.5 * weight * f(start + 0.5 * extent * (1.0 + sign * node));
        }
    }
    sum
}

fn shutter(window: HarmonicWindow, mut f: impl FnMut(f64) -> f64) -> f64 {
    match window.profile() {
        HarmonicWindowProfile::Uniform => average(window.start_s(), window.duration_s(), f),
        HarmonicWindowProfile::Triangular => {
            // Split at the triangle's cusp before smooth quadrature.
            average(0.0, 0.5, |u| {
                2.0 * u * f(window.start_s() + window.duration_s() * u)
            }) + average(0.5, 0.5, |u| {
                2.0 * (1.0 - u) * f(window.start_s() + window.duration_s() * u)
            })
        }
    }
}

fn moment(waves: &[Wave], position: [f64; 3], time: f64, channel: usize, proper: bool) -> f64 {
    let mut field = [C64::ZERO; 6];
    for wave in waves {
        let phase = rotation(
            wave.omega * time
                - wave.k[0] * position[0]
                - wave.k[1] * position[1]
                - wave.k[2] * position[2],
        );
        for (sum, amplitude) in field.iter_mut().zip(wave.field) {
            *sum = *sum + amplitude * phase;
        }
    }
    let product = |a: usize, b: usize| {
        if proper {
            0.5 * (field[a].re * field[b].re + field[a].im * field[b].im)
        } else {
            field[a].re * field[b].re
        }
    };
    if channel < 6 {
        product(channel, channel)
    } else {
        let component = channel - 6;
        let (a, b) = ((component + 1) % 3, (component + 2) % 3);
        product(a, b + 3) - product(b, a + 3)
    }
}

fn reference_power(
    waves: &[Wave],
    depth: f64,
    area: AngularSpectrumAperture,
    window: HarmonicWindow,
    proper: bool,
) -> f64 {
    let [x, y] = area.origin_m();
    let [width, height] = area.extent_m();
    area.area_m2()
        * average(x, width, |x| {
            average(y, height, |y| {
                shutter(window, |t| moment(waves, [x, y, depth], t, 8, proper))
            })
        })
}

fn values(point: &fs_em::plane_wave::angular_spectrum::AngularSpectrumExposurePoint) -> [f64; 9] {
    let e = point.electric_mean_squared_v2_m2;
    let h = point.magnetic_mean_squared_a2_m2;
    let p = point.poynting_w_m2;
    [e[0], e[1], e[2], h[0], h[1], h[2], p[0], p[1], p[2]]
}

fn beam_limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 32,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}

fn aperture_limits() -> AngularSpectrumExposureApertureLimits {
    AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 4096,
        maximum_workspace_bytes: 8 << 20,
    }
}

#[test]
fn g1_modulated_carrier_refraction_matches_real_area_and_shutter_quadrature() {
    let air = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stacks = [1.6, 2.1].map(|omega| {
        PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&glass), omega, &[])
            .unwrap()
    });
    let cell = AngularSpectrumGrid {
        width: 4,
        height: 2,
        spacing_m: [4.0, 5.0],
    };
    let carrier = [0.21, -0.13];
    let offsets = [0.0, core::f64::consts::TAU / 16.0];
    let amplitudes = [
        [
            [C64::new(0.6, -0.2), C64::new(0.3, 0.5)],
            [C64::new(-0.4, 0.7), C64::new(0.2, -0.3)],
        ],
        [
            [C64::new(-0.1, 0.3), C64::new(0.7, -0.2)],
            [C64::new(0.5, 0.1), C64::new(-0.3, 0.4)],
        ],
    ];
    let mut waves = Vec::new();
    let mut lines = Vec::new();
    for (frequency, stack) in stacks.iter().enumerate() {
        let omega = [1.6, 2.1][frequency];
        let mut samples = vec![[C64::ZERO; 2]; 8];
        for (mode, offset) in offsets.iter().enumerate() {
            let amplitude = amplitudes[frequency][mode];
            for (point, sample) in samples.iter_mut().enumerate() {
                let phase = rotation(-offset * (point % cell.width) as f64 * cell.spacing_m[0]);
                for component in 0..2 {
                    sample[component] = sample[component] + amplitude[component] * phase;
                }
            }
            // Independent TE/TM half-space Fresnel factors, preserving the
            // same Cartesian tangential electric basis at both frequencies.
            let [kx, ky] = [carrier[0] + offset, carrier[1]];
            let kt2 = kx * kx + ky * ky;
            let qi = (omega * omega - kt2).sqrt();
            let qt = (2.25 * omega * omega - kt2).sqrt();
            let p = [kx / kt2.sqrt(), ky / kt2.sqrt()];
            let s = [-p[1], p[0]];
            let te =
                (amplitude[0].scale(s[0]) + amplitude[1].scale(s[1])).scale(2.0 * qi / (qi + qt));
            let tm = (amplitude[0].scale(p[0]) + amplitude[1].scale(p[1]))
                .scale(2.0 * qt / (qt + 2.25 * qi));
            waves.push(wave(
                omega,
                [kx, ky, qt],
                [
                    te.scale(s[0]) + tm.scale(p[0]),
                    te.scale(s[1]) + tm.scale(p[1]),
                ],
            ));
        }
        lines.push(
            PlanarAngularSpectrum::prepare_with_carrier(
                stack,
                cell,
                carrier,
                &samples,
                beam_limits(),
                || false,
            )
            .unwrap(),
        );
    }
    let factors = [
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap(),
    ];
    let bank = &factors[0];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors,
        &[
            AngularSpectrumFactorStatistics::Real,
            AngularSpectrumFactorStatistics::ProperComplex,
        ],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let full = AngularSpectrumAperture::try_new([0.0; 2], [16.0, 10.0]).unwrap();
    let patch = AngularSpectrumAperture::try_new([-0.7, 0.4], [3.5, 1.8]).unwrap();
    let depth = 0.37;
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let field = bank
            .field(Side::Transmitted, depth, window, || false)
            .unwrap();
        let statistical = ensemble
            .field(Side::Transmitted, depth, window, || false)
            .unwrap();
        for (point, actual) in field.points().iter().enumerate() {
            let position = [(point % 4) as f64 * 4.0, (point / 4) as f64 * 5.0, depth];
            for (channel, actual) in values(actual).iter().enumerate() {
                let real = shutter(window, |t| moment(&waves, position, t, channel, false));
                let proper = shutter(window, |t| moment(&waves, position, t, channel, true));
                near(*actual, real, 3e-12);
                near(
                    values(&statistical.points()[point])[channel],
                    real + proper,
                    5e-12,
                );
            }
        }
        let real = reference_power(&waves, depth, full, window, false);
        let proper = reference_power(&waves, depth, full, window, true);
        near(field.normal_power_w(), real, 2e-9);
        near(statistical.normal_power_w(), real + proper, 3e-9);
        for area in [full, patch] {
            let measured = bank
                .measure_aperture(
                    Side::Transmitted,
                    depth,
                    area,
                    window,
                    aperture_limits(),
                    || false,
                )
                .unwrap();
            let expected = reference_power(&waves, depth, area, window, false);
            near(measured.normal_power_w, expected, 2e-9);
        }
    }
    assert_eq!(bank.point_frequency_pairs(), 8 * 2 * 2 + 2 * 8 * (4 + 2));
    for bounded in [
        AngularSpectrumExposureLimits {
            maximum_point_frequency_pairs: bank.point_frequency_pairs() - 1,
            ..limits()
        },
        AngularSpectrumExposureLimits {
            maximum_workspace_bytes: bank.workspace_bytes() - 1,
            ..limits()
        },
    ] {
        assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, bounded, || false).is_err());
    }
}

#[test]
fn g1_two_sided_carrier_exposure_preserves_film_phase_and_signed_flux() {
    let medium = material(1.0, 1.0);
    let thickness = 0.41;
    let stacks = [1.2, 1.7].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&medium),
                thickness_m: thickness,
            }],
        )
        .unwrap()
    });
    let carrier = [0.3_f64, 0.2];
    let left = [
        [C64::new(0.7, 0.2), C64::new(-0.3, 0.4)],
        [C64::new(-0.2, 0.6), C64::new(0.5, -0.1)],
    ];
    let right = [
        [C64::new(-0.4, 0.5), C64::new(0.1, -0.2)],
        [C64::new(0.3, -0.1), C64::new(-0.6, 0.2)],
    ];
    let lines: Vec<_> = stacks
        .iter()
        .enumerate()
        .map(|(i, stack)| {
            PlanarAngularSpectrum::prepare_with_carrier(
                stack,
                grid(),
                carrier,
                &[left[i]; 4],
                beam_limits(),
                || false,
            )
            .unwrap()
        })
        .collect();
    let illumination: Vec<_> = lines
        .iter()
        .enumerate()
        .map(|(i, beam)| {
            beam.prepare_illumination(
                &[right[i]; 4],
                AngularSpectrumIlluminationLimits {
                    maximum_layer_mode_pairs: 4096,
                    maximum_workspace_bytes: 8 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
        })
        .collect();
    let bank =
        AngularSpectrumExposure::try_new_illumination(&lines[0], &illumination, limits(), || false)
            .unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let mut waves = Vec::new();
    for i in 0..2 {
        let omega = [1.2, 1.7][i];
        let q = (omega * omega - carrier[0] * carrier[0] - carrier[1] * carrier[1]).sqrt();
        waves.push(wave(omega, [carrier[0], carrier[1], q], left[i]));
        waves.push(wave(
            omega,
            [carrier[0], carrier[1], -q],
            right[i].map(|e| e * rotation(-q * thickness)),
        ));
    }
    let full = AngularSpectrumAperture::try_new([0.0; 2], [32.0, 3.0]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.37, profile).unwrap();
        for (face, distance, z) in [
            (Face::TotalLeft, 0.17, -0.17),
            (Face::TotalRight, 0.21, thickness + 0.21),
        ] {
            let field = bank
                .illumination_field(face, distance, window, || false)
                .unwrap();
            near(
                field.normal_power_w(),
                reference_power(&waves, z, full, window, false),
                3e-9,
            );
            let aperture = bank
                .measure_illumination_aperture(
                    face,
                    distance,
                    full,
                    window,
                    aperture_limits(),
                    || false,
                )
                .unwrap();
            near(aperture.normal_power_w, field.normal_power_w(), 3e-12);
        }
        for depth in [0.0, 0.19, thickness] {
            let field = layer.field(depth, window, || false).unwrap();
            near(
                field.field().normal_power_w(),
                reference_power(&waves, depth, full, window, false),
                3e-9,
            );
            for (channel, actual) in values(&field.field().points()[0]).iter().enumerate() {
                near(
                    *actual,
                    shutter(window, |t| {
                        moment(&waves, [0.0, 0.0, depth], t, channel, false)
                    }),
                    3e-12,
                );
            }
        }
    }
}

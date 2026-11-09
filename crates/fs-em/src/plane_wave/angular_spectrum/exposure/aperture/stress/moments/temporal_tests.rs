use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn wave(amplitude: C64) -> [C64; 6] {
    [
        amplitude,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        amplitude.scale(-1.0),
        C64::ZERO,
    ]
}

#[allow(clippy::too_many_arguments)]
fn observe<K: MomentKernel>(
    frequencies: &[f64],
    kx: &[f64],
    ky: &[f64],
    fields: &[Vec<[C64; 6]>],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    positive: bool,
) -> AngularSpectrumExposureApertureWrench {
    let grid = AngularSpectrumGrid {
        width: kx.len(),
        height: ky.len(),
        spacing_m: [1.0; 2],
    };
    let spatial = MomentSpatial::<K>::from_axes(
        grid,
        aperture,
        |i| Dd::from_f64(kx[i]),
        |i| Dd::from_f64(ky[i]),
        &mut || false,
    )
    .unwrap();
    let mut parts: Vec<Vec<FieldParts>> = fields
        .iter()
        .map(|field| {
            field
                .iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect()
        })
        .collect();
    if K::CENTER_FIELDS {
        for (f, field) in parts.iter_mut().enumerate() {
            center_fields(
                field,
                frequencies[f],
                |i| kx[i % kx.len()],
                |i| ky[i / kx.len()],
                aperture,
                window,
                &mut || false,
            )
            .unwrap();
        }
    }
    let mut sum = Wrench::new(Dd::ONE, Dd::ONE);
    if positive {
        let modes = if K::CENTER_FIELDS {
            window.real_band_mode_count_rad_s(frequencies)
        } else {
            window.band_mode_count_rad_s(frequencies)
        }
        .unwrap()
        .unwrap();
        let fields = temporal::project::<K>(
            window,
            frequencies,
            &parts,
            grid.width * grid.height,
            modes,
            &mut || false,
        )
        .unwrap();
        integrate_projected_moments(&spatial, &fields, modes, &mut sum, &mut || false).unwrap();
    } else {
        let time = pair_kernels(frequencies.len(), window, |i| -frequencies[i], &mut || {
            false
        })
        .unwrap();
        integrate_moments(&spatial, &time, &parts, &mut sum, &mut || false).unwrap();
    }
    sum.finish(aperture).unwrap()
}

#[test]
fn short_real_two_tone_window_retains_surface_traction() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [0.7, 0.9]).unwrap();
    let fields = [vec![wave(C64::ONE)], vec![wave(C64::from_re(-1.0))]];
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 9.0 / 320.0),
        (HarmonicWindowProfile::Triangular, 9.0 / 960.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let result =
            observe::<RealKernel>(&[1.0, 2.0], &[0.0], &[0.0], &fields, area, window, true);
        let expected = -area.area_m2() * coefficient * duration.powi(4);
        close(result.integrated_traction_n()[2] / expected, 1.0, 2e-12);
        assert_eq!(result.centered_traction_moment_nm(), [0.0; 3]);
    }
}

#[test]
fn dark_proper_interference_retains_signed_fringe_torque_before_power() {
    let duration: f64 = 1e-12;
    let alpha = 0.6;
    let area = AngularSpectrumAperture::try_new([-0.4, -0.3], [0.8, 0.6]).unwrap();
    // S(x)=1+alpha*sin(x) is represented by the original +/- spatial modes.
    let shape = [
        C64::ONE,
        C64::new(0.0, alpha / 2.0),
        C64::new(0.0, -alpha / 2.0),
    ];
    let fields: Vec<_> = [1.0, -2.0, 1.0]
        .map(|a| shape.map(|v| wave(v.scale(a))).to_vec())
        .into();
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 1.0 / 160.0),
        (HarmonicWindowProfile::Triangular, 1.0 / 480.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let result = observe::<ProperKernel>(
            &[1e15, 1e15 + 1.0, 1e15 + 2.0],
            &[0.0, 1.0, -1.0],
            &[0.0],
            &fields,
            area,
            window,
            true,
        );
        let temporal = coefficient * duration.powi(4);
        let width = area.extent_m()[0];
        let z = width / 2.0;
        let mean_square = 1.0 + alpha * alpha * (1.0 - width.sin() / width) / 2.0;
        let first = 2.0 * alpha * (z.sin() - z * z.cos()) / z;
        close(
            result.integrated_traction_n()[2] / (-area.area_m2() * mean_square * temporal),
            1.0,
            3e-12,
        );
        close(
            result.centered_traction_moment_nm()[1] / (area.area_m2() * first * temporal),
            1.0,
            3e-12,
        );
        assert!(result.centered_traction_moment_nm()[1] > 0.0);
        let offset = [0.2, -0.1, 0.3];
        let shifted = result.traction_moment_about_nm(offset).unwrap();
        close(
            shifted[1]
                / (result.centered_traction_moment_nm()[1]
                    + offset[0] * result.integrated_traction_n()[2]),
            1.0,
            3e-12,
        );
    }
}

#[test]
fn projected_real_and_proper_first_moments_match_resolved_original_kernel() {
    let (frequencies, kx, ky) = ([1.3, 2.1], [0.2, -0.7], [0.3, -0.5]);
    let fields: Vec<_> = (0..2)
        .map(|f| {
            (0..4)
                .map(|m| {
                    core::array::from_fn(|c| {
                        C64::new(
                            ((f + 1) * (c + 1)) as f64 / (7.0 * (m + 1) as f64),
                            (c as f64 - 2.0) * (m as f64 - 1.5) / (11.0 * (f + 1) as f64),
                        )
                    })
                })
                .collect()
        })
        .collect();
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for proper in [false, true] {
            let eval = |positive| {
                if proper {
                    observe::<ProperKernel>(&frequencies, &kx, &ky, &fields, area, window, positive)
                } else {
                    observe::<RealKernel>(&frequencies, &kx, &ky, &fields, area, window, positive)
                }
            };
            let (actual, original) = (eval(true), eval(false));
            for c in 0..3 {
                close(
                    actual.integrated_traction_n()[c],
                    original.integrated_traction_n()[c],
                    3e-15,
                );
                close(
                    actual.centered_traction_moment_nm()[c],
                    original.centered_traction_moment_nm()[c],
                    3e-15,
                );
            }
        }
    }
}

#[test]
fn projected_first_moments_poll_cancellation() {
    let area = AngularSpectrumAperture::try_new([-0.4, -0.3], [0.8, 0.6]).unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let spatial = MomentSpatial::<RealKernel>::from_axes(
        grid,
        area,
        |i| Dd::from_f64(i as f64),
        |_| Dd::ZERO,
        &mut || false,
    )
    .unwrap();
    let fields = vec![[[Dd::ONE, Dd::ZERO]; 6]; 8];
    let mut polls = 0;
    assert!(matches!(
        integrate_projected_moments(
            &spatial,
            &fields,
            4,
            &mut Wrench::new(Dd::ONE, Dd::ONE),
            &mut || {
                polls += 1;
                polls == 4
            }
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}

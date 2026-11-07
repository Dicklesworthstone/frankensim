use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
#[allow(clippy::too_many_arguments)]
fn accumulate<K: MomentKernel>(
    values: &[Vec<[C64; 6]>],
    omega: &[f64],
    kx: &[f64],
    ky: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    sum: &mut Wrench,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
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
        cancelled,
    )?;
    let temporal = pair_kernels(omega.len(), window, |i| -omega[i], cancelled)?;
    let mut fields: Vec<Vec<FieldParts>> = values
        .iter()
        .map(|line| {
            line.iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect()
        })
        .collect();
    if K::CENTER_FIELDS {
        for (a, line) in fields.iter_mut().enumerate() {
            center_fields(
                line,
                omega[a],
                |i| kx[i % kx.len()],
                |i| ky[i / kx.len()],
                aperture,
                window,
                cancelled,
            )?;
        }
    }
    integrate_moments(&spatial, &temporal, &fields, sum, cancelled)
}
fn observe<K: MomentKernel>(
    values: &[Vec<[C64; 6]>],
    omega: &[f64],
    kx: &[f64],
    ky: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> AngularSpectrumExposureApertureWrench {
    let mut sum = Wrench::new(Dd::from_f64(2.0), Dd::from_f64(8.0));
    accumulate::<K>(
        values,
        omega,
        kx,
        ky,
        aperture,
        window,
        &mut sum,
        &mut || false,
    )
    .unwrap();
    sum.finish(aperture).unwrap()
}

#[test]
fn coherent_fringe_has_analytic_moment_that_separate_mode_powers_lose() {
    let a: f64 = 0.7;
    let alpha = (1.0 - (a / 4.0).powi(2)).sqrt();
    let first = [0.0, 1.0, 0.0, -alpha / 2.0, 0.0, a / 8.0].map(C64::from_re);
    let second = [0.0, 1.0, 0.0, -alpha / 2.0, 0.0, -a / 8.0].map(|v| C64::new(0.0, v));
    let values = vec![vec![first, second]];
    let aperture = AngularSpectrumAperture::try_new([-0.6, -0.4], [1.2, 0.8]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let proper = observe::<ProperKernel>(&values, &[1.0], &[a, -a], &[0.0], aperture, window);
    let real = observe::<RealKernel>(&values, &[1.0], &[a, -a], &[0.0], aperture, window);
    // Direct integration of Tzz=-2*alpha^2 + 2*sin(2*a*x).
    let z = a * aperture.extent_m()[0];
    let moment_x_sine = 0.5 * aperture.extent_m()[0] * (z.sin() - z * z.cos()) / (z * z);
    let expected = -2.0 * aperture.area_m2() * moment_x_sine;
    close(proper.centered_traction_moment_nm()[1], expected, 2e-14);
    assert!(expected.abs() > 0.1);
    close(
        proper.integrated_traction_n()[2],
        -2.0 * alpha * alpha * aperture.area_m2(),
        2e-14,
    );
    for c in 0..3 {
        close(
            real.integrated_traction_n()[c],
            proper.integrated_traction_n()[c],
            2e-14,
        );
        close(
            real.centered_traction_moment_nm()[c],
            proper.centered_traction_moment_nm()[c],
            2e-14,
        );
    }
    for (field, k) in [(first, a), (second, -a)] {
        let alone = observe::<ProperKernel>(&[vec![field]], &[1.0], &[k], &[0.0], aperture, window);
        assert_eq!(alone.centered_traction_moment_nm(), [0.0; 3]);
    }
}

fn scene(omega: &[f64], kx: &[f64], ky: &[f64]) -> Vec<Vec<[C64; 6]>> {
    omega
        .iter()
        .enumerate()
        .map(|(line, &w)| {
            ky.iter()
                .enumerate()
                .flat_map(|(iy, &y)| {
                    kx.iter().enumerate().map(move |(ix, &x)| {
                        let k = 4.0 * w;
                        let transverse = x.hypot(y);
                        let sign = if (ix + iy + line) % 2 == 0 { 1.0 } else { -1.0 };
                        let n = [x / k, y / k, sign * (1.0 - (transverse / k).powi(2)).sqrt()];
                        let te = [-y / transverse, x / transverse, 0.0];
                        let tm = [
                            n[1] * te[2] - n[2] * te[1],
                            n[2] * te[0] - n[0] * te[2],
                            n[0] * te[1] - n[1] * te[0],
                        ];
                        let s = C64::new(0.3 + 0.09 * (line + ix) as f64, -0.13 * (iy + 1) as f64);
                        let p = C64::new(-0.17 * (ix + 1) as f64, 0.11 * (line + iy + 1) as f64);
                        let e: [C64; 3] = core::array::from_fn(|c| {
                            s * C64::from_re(te[c]) + p * C64::from_re(tm[c])
                        });
                        let h: [C64; 3] = core::array::from_fn(|c| {
                            let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                            (e[v] * C64::from_re(n[u]) - e[u] * C64::from_re(n[v]))
                                * C64::from_re(0.5)
                        });
                        [e[0], e[1], e[2], h[0], h[1], h[2]]
                    })
                })
                .collect()
        })
        .collect()
}

// Independent tensor-product Simpson rule on real Maxwell fields. It does not
// use production aperture, moment, stress or harmonic-product kernels.
fn quadrature(
    values: &[Vec<[C64; 6]>],
    omega: &[f64],
    kx: &[f64],
    ky: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> ([f64; 3], [f64; 3]) {
    let n = 96usize;
    let weight = |i| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let origin = aperture.origin_m();
    let extent = aperture.extent_m();
    let (mut force, mut moment) = ([0.0; 3], [0.0; 3]);
    for ix in 0..=n {
        for iy in 0..=n {
            for it in 0..=n {
                let u = [
                    ix as f64 / n as f64,
                    iy as f64 / n as f64,
                    it as f64 / n as f64,
                ];
                let (x, y, t) = (
                    origin[0] + extent[0] * u[0],
                    origin[1] + extent[1] * u[1],
                    window.start_s() + window.duration_s() * u[2],
                );
                let mut field = [0.0; 6];
                for (a, line) in values.iter().enumerate() {
                    for (mode, value) in line.iter().enumerate() {
                        let phase =
                            omega[a] * t - kx[mode % kx.len()] * x - ky[mode / kx.len()] * y;
                        let (s, c) = phase.sin_cos();
                        for j in 0..6 {
                            field[j] += value[j].re * c - value[j].im * s;
                        }
                    }
                }
                let energy = (0..3)
                    .map(|c| field[c] * field[c] + 4.0 * field[c + 3] * field[c + 3])
                    .sum::<f64>();
                let traction: [f64; 3] = core::array::from_fn(|c| {
                    2.0 * field[c] * field[2] + 8.0 * field[c + 3] * field[5]
                        - if c == 2 { energy } else { 0.0 }
                });
                let shutter = match window.profile() {
                    HarmonicWindowProfile::Uniform => 1.0,
                    HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u[2] - 1.0).abs()),
                };
                let w = weight(ix) * weight(iy) * weight(it) * shutter * aperture.area_m2()
                    / (27 * n * n * n) as f64;
                let lever = [extent[0] * (u[0] - 0.5), extent[1] * (u[1] - 0.5), 0.0];
                for c in 0..3 {
                    force[c] += w * traction[c];
                    let (j, k) = ((c + 1) % 3, (c + 2) % 3);
                    moment[c] += w * (lever[j] * traction[k] - lever[k] * traction[j]);
                }
            }
        }
    }
    (force, moment)
}

#[test]
fn mixed_polarization_two_dimensional_fringe_moments_match_real_field_quadrature() {
    let (omega, kx, ky) = ([3.2, 4.7], [0.7, -0.8], [0.4, -0.9]);
    let values = scene(&omega, &kx, &ky);
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.23], [0.63, 0.49]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.31, 0.42, profile).unwrap();
        let actual = observe::<RealKernel>(&values, &omega, &kx, &ky, aperture, window);
        let (force, moment) = quadrature(&values, &omega, &kx, &ky, aperture, window);
        for c in 0..3 {
            close(actual.integrated_traction_n()[c], force[c], 2e-8);
            close(actual.centered_traction_moment_nm()[c], moment[c], 2e-8);
        }
    }
}

#[test]
fn proper_first_moments_equal_independent_real_quadratures() {
    let (omega, kx, ky) = ([3.2, 4.7], [0.7, -0.8], [0.4, -0.9]);
    let values = scene(&omega, &kx, &ky);
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.23], [0.63, 0.49]).unwrap();
    let window =
        HarmonicWindow::with_profile(0.31, 0.42, HarmonicWindowProfile::Triangular).unwrap();
    let proper = observe::<ProperKernel>(&values, &omega, &kx, &ky, aperture, window);
    let mut sum = Wrench::new(Dd::from_f64(2.0), Dd::from_f64(8.0));
    for scale in [
        C64::from_re(core::f64::consts::FRAC_1_SQRT_2),
        C64::new(0.0, core::f64::consts::FRAC_1_SQRT_2),
    ] {
        let rotated: Vec<Vec<_>> = values
            .iter()
            .map(|line| line.iter().map(|v| v.map(|z| z * scale)).collect())
            .collect();
        accumulate::<RealKernel>(
            &rotated,
            &omega,
            &kx,
            &ky,
            aperture,
            window,
            &mut sum,
            &mut || false,
        )
        .unwrap();
    }
    let quadratures = sum.finish(aperture).unwrap();
    for c in 0..3 {
        close(
            proper.integrated_traction_n()[c],
            quadratures.integrated_traction_n()[c],
            3e-14,
        );
        close(
            proper.centered_traction_moment_nm()[c],
            quadratures.centered_traction_moment_nm()[c],
            3e-14,
        );
    }
}

#[test]
fn partitioned_apertures_add_after_translating_moment_origins() {
    let (omega, kx, ky) = ([3.2, 4.7], [0.7, -0.8], [0.4, -0.9]);
    let values = scene(&omega, &kx, &ky);
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.23], [0.63, 0.49]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.42).unwrap();
    let whole = observe::<RealKernel>(&values, &omega, &kx, &ky, aperture, window);
    let (mut force, mut moment) = ([0.0; 3], [0.0; 3]);
    for (origin, width) in [(0.37, 0.29), (0.66, 0.34)] {
        let patch = AngularSpectrumAperture::try_new([origin, -0.23], [width, 0.49]).unwrap();
        let part = observe::<RealKernel>(&values, &omega, &kx, &ky, patch, window);
        let shifted = part
            .traction_moment_about_nm([0.37 + 0.63 / 2.0 - origin - width / 2.0, 0.0, 0.0])
            .unwrap();
        for c in 0..3 {
            force[c] += part.integrated_traction_n()[c];
            moment[c] += shifted[c];
        }
    }
    for c in 0..3 {
        close(force[c], whole.integrated_traction_n()[c], 3e-14);
        close(moment[c], whole.centered_traction_moment_nm()[c], 3e-14);
    }
    let offset = [0.41, -0.2, 0.7];
    let translated = whole.traction_moment_about_nm(offset).unwrap();
    let f = whole.integrated_traction_n();
    let m = whole.centered_traction_moment_nm();
    for c in 0..3 {
        let (u, v) = ((c + 1) % 3, (c + 2) % 3);
        close(
            translated[c],
            m[c] - offset[u] * f[v] + offset[v] * f[u],
            3e-14,
        );
    }
    assert!(
        whole
            .traction_moment_about_nm([f64::INFINITY, 0.0, 0.0])
            .is_err()
    );
}

#[test]
fn narrow_interval_first_moments_keep_the_resolved_derivative() {
    let k = 1e-10;
    let width = 1e-5;
    let actual = centered_sine_moment(Dd::from_f64(k), width)
        .unwrap()
        .to_f64();
    let expected = k * width * width / 12.0;
    assert!((actual / expected - 1.0).abs() < 1e-15);
    assert_eq!(centered_sine_moment(Dd::ZERO, width).unwrap(), Dd::ZERO);
    let large = centered_sine_moment(Dd::from_f64(1e160), 1.0).unwrap();
    assert!(finite_dd(large) && large.to_f64().is_finite());
    assert!(centered_sine_moment(Dd::from_f64(1e-300), 1e-300).is_err());
    assert!(centered_sine_moment(Dd::from_f64(f64::INFINITY), width).is_err());

    // Moving the coordinate origin into an authored field phase must leave
    // the physical first moment unchanged, even for a very distant detector.
    let q = Dd::from_f64(1.1);
    let origin = 1e30;
    let width = 1e-7;
    let distant =
        ProperKernel::prepare_first(HarmonicWindow::new(origin, width).unwrap(), q, Dd::ZERO)
            .unwrap()
            .product_first([Dd::ONE, Dd::ZERO], [Dd::ONE, Dd::ZERO]);
    let (s, c) = fs_math::det::sin_cos_dd(q * Dd::from_f64(origin));
    let local = ProperKernel::prepare_first(HarmonicWindow::new(0.0, width).unwrap(), q, Dd::ZERO)
        .unwrap()
        .product_first([c, -s], [Dd::ONE, Dd::ZERO]);
    let scale = centered_sine_moment(q, width).unwrap().to_f64();
    assert!(((distant - local).to_f64() / scale).abs() < 1e-14);
}

#[test]
fn moment_work_storage_cancellation_and_derived_units_are_bounded() {
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 2,
        spacing_m: [1.0; 2],
    };
    let extra = moment_cache_bytes(grid, true, true).unwrap();
    let limits = AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 64,
        maximum_workspace_bytes: 1 << 20,
    };
    assert!(admit::<Wrench>(grid, 4, extra, true, true, limits).is_ok());
    assert!(
        admit::<Wrench>(
            grid,
            4,
            extra,
            true,
            true,
            AngularSpectrumExposureApertureLimits {
                maximum_mode_frequency_pairs: 63,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        admit::<Wrench>(
            grid,
            4,
            extra,
            true,
            true,
            AngularSpectrumExposureApertureLimits {
                maximum_workspace_bytes: extra,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        moment_cache_bytes(
            AngularSpectrumGrid {
                width: usize::MAX,
                ..grid
            },
            true,
            true
        )
        .is_err()
    );
    let (omega, kx, ky) = ([3.2, 4.7], [0.7, -0.8], [0.4, -0.9]);
    let values = scene(&omega, &kx, &ky);
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.23], [0.63, 0.49]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.42).unwrap();
    let mut sum = Wrench::new(Dd::from_f64(2.0), Dd::from_f64(8.0));
    let mut polls = 0;
    assert!(
        accumulate::<RealKernel>(
            &values,
            &omega,
            &kx,
            &ky,
            aperture,
            window,
            &mut sum,
            &mut || {
                polls += 1;
                polls == 12
            }
        )
        .is_err()
    );
    let a = observe::<RealKernel>(&values, &omega, &kx, &ky, aperture, window);
    assert_eq!(
        a,
        observe::<RealKernel>(&values, &omega, &kx, &ky, aperture, window)
    );
    let mut huge = Wrench::new(Dd::ONE, Dd::ONE);
    huge.force[2] = Dd::from_f64(25.0);
    assert!(
        huge.finish(AngularSpectrumAperture::try_new([0.0; 2], [1e154; 2]).unwrap())
            .is_err()
    );
    let large = AngularSpectrumExposureApertureWrench {
        aperture,
        force: [Dd::from_f64(1e308); 3],
        moment: [Dd::ZERO; 3],
    };
    assert!(large.traction_moment_about_nm([0.0, 10.0, 0.0]).is_err());
}

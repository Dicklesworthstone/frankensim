use super::*;

fn fields(input: &[Vec<[C64; 6]>]) -> Vec<Vec<FieldParts>> {
    input
        .iter()
        .map(|line| {
            line.iter()
                .map(|value| value.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect()
        })
        .collect()
}

fn accumulate<K: Kernel>(
    values: &[Vec<[C64; 6]>],
    omega: &[f64],
    kx: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    sum: &mut Stress,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let grid = AngularSpectrumGrid {
        width: kx.len(),
        height: 1,
        spacing_m: [1.0; 2],
    };
    let spatial = Spatial::<K>::from_axes(grid, aperture, |i| kx[i], |_| 0.0, cancelled)?;
    let temporal = pair_kernels(omega.len(), window, |i| -omega[i], cancelled)?;
    let mut values = fields(values);
    if K::CENTER_FIELDS {
        for (i, line) in values.iter_mut().enumerate() {
            center_fields(
                line,
                omega[i],
                |i| kx[i],
                |_| 0.0,
                aperture,
                window,
                cancelled,
            )?;
        }
    }
    integrate(&spatial, &temporal, &values, sum, cancelled)
}

fn observe<K: Kernel>(
    values: &[Vec<[C64; 6]>],
    omega: &[f64],
    kx: &[f64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> AngularSpectrumExposureApertureStress {
    let mut sum = Stress::new(2.0, 8.0);
    accumulate::<K>(values, omega, kx, aperture, window, &mut sum, &mut || false).unwrap();
    sum.finish(aperture).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

#[test]
fn plane_wave_energy_flux_tensor_and_signed_surface_traction() {
    let aperture = AngularSpectrumAperture::try_new([0.3, -0.8], [0.4, 0.7]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let values = vec![vec![[3.0, 4.0, 0.0, -2.0, 1.5, 0.0].map(C64::from_re)]];
    let result = observe::<RealKernel>(&values, &[1.0], &[0.0], aperture, window);
    // epsilon=2, mu=8: Z=2, c=1/4, peak |E|=5.
    close(result.mean_field_energy_density_j_m3(), 25.0, 2e-14);
    close(result.mean_poynting_w_m2()[2], 6.25, 2e-14);
    let tensor = result.mean_maxwell_stress_pa();
    for i in 0..3 {
        for j in 0..3 {
            close(
                tensor[i][j],
                if i == 2 && j == 2 { -25.0 } else { 0.0 },
                2e-14,
            );
        }
    }
    let forward = result.integrated_traction_n([0.0, 0.0, 1.0]).unwrap();
    let reverse = result.integrated_traction_n([0.0, 0.0, -1.0]).unwrap();
    for i in 0..3 {
        assert_eq!(forward[i], -reverse[i]);
    }
    close(forward[2], -25.0 * aperture.area_m2(), 2e-14);
    close(
        result.normal_power_w().unwrap(),
        6.25 * aperture.area_m2(),
        2e-14,
    );
    close(
        result.integrated_stress_n().unwrap()[2][2],
        forward[2],
        2e-14,
    );
    assert!(result.mean_traction_pa([1.0, 1.0, 0.0]).is_err());
    assert!(result.mean_traction_pa([f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn standing_wave_node_retains_magnetic_energy_and_anisotropic_stress() {
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    // The original counterpropagating fields have cancelled E at this plane.
    let values = vec![vec![[0.0, 0.0, 0.0, 0.0, 2.0, 0.0].map(C64::from_re)]];
    let result = observe::<RealKernel>(&values, &[1.0], &[0.0], aperture, window);
    assert_eq!(result.mean_poynting_w_m2(), [0.0; 3]);
    close(result.mean_field_energy_density_j_m3(), 8.0, 2e-14);
    let tensor = result.mean_maxwell_stress_pa();
    for i in 0..3 {
        close(tensor[i][i], if i == 1 { 8.0 } else { -8.0 }, 2e-14);
    }
}

fn scene() -> Vec<Vec<[C64; 6]>> {
    (0..2)
        .map(|line| {
            (0..2)
                .map(|mode| {
                    core::array::from_fn(|c| {
                        C64::new(
                            (1 + c + 2 * line) as f64 * (if mode == 0 { 0.17 } else { -0.11 }),
                            (1 + mode + line + c % 3) as f64
                                * (if c % 2 == 0 { 0.09 } else { -0.13 }),
                        )
                    })
                })
                .collect()
        })
        .collect()
}

#[test]
fn finite_patch_carrier_beats_and_longitudinal_products_match_independent_quadrature() {
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.2], [0.63, 0.4]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.42).unwrap();
    let values = scene();
    let omega = [3.2, 4.7];
    let kx = [0.0, 1.7];
    let result = observe::<RealKernel>(&values, &omega, &kx, aperture, window);
    // Independent tensor construction on a refined tensor-product Simpson rule.
    // No production harmonic-product or stress helper is used by the oracle.
    let n = 96;
    let weight = |i| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let (mut stress, mut flux, mut energy) = ([[0.0; 3]; 3], [0.0; 3], 0.0);
    for ix in 0..=n {
        for it in 0..=n {
            let x = 0.37 + 0.63 * ix as f64 / n as f64;
            let t = 0.31 + 0.42 * it as f64 / n as f64;
            let mut field = [0.0; 6];
            for a in 0..2 {
                for b in 0..2 {
                    let (sin, cos) = (omega[a] * t - kx[b] * x).sin_cos();
                    for c in 0..6 {
                        field[c] += values[a][b][c].re * cos - values[a][b][c].im * sin;
                    }
                }
            }
            let u = (0..3)
                .map(|i| field[i] * field[i] + 4.0 * field[i + 3] * field[i + 3])
                .sum::<f64>();
            let w = weight(ix) * weight(it) / (9 * n * n) as f64;
            energy += w * u;
            for i in 0..3 {
                let (a, b) = ((i + 1) % 3, (i + 2) % 3);
                flux[i] += w * (field[a] * field[b + 3] - field[b] * field[a + 3]);
                for j in 0..3 {
                    stress[i][j] += w
                        * (2.0 * field[i] * field[j] + 8.0 * field[i + 3] * field[j + 3]
                            - if i == j { u } else { 0.0 });
                }
            }
        }
    }
    close(result.mean_field_energy_density_j_m3(), energy, 2e-8);
    for i in 0..3 {
        close(result.mean_poynting_w_m2()[i], flux[i], 2e-8);
        for j in 0..3 {
            close(result.mean_maxwell_stress_pa()[i][j], stress[i][j], 2e-8);
        }
    }
    let shifted = observe::<RealKernel>(
        &values,
        &omega,
        &kx,
        AngularSpectrumAperture::try_new([0.81, -0.2], [0.63, 0.4]).unwrap(),
        window,
    );
    assert!((shifted.mean_field_energy_density_j_m3() - energy).abs() > 0.01);
}

#[test]
fn proper_factors_equal_two_real_quadratures_and_independent_banks_do_not_cancel() {
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.2], [0.63, 0.4]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.42).unwrap();
    let values = scene();
    let omega = [3.2, 4.7];
    let kx = [0.0, 1.7];
    let proper = observe::<ProperKernel>(&values, &omega, &kx, aperture, window);
    let mut sum = Stress::new(2.0, 8.0);
    for multiplier in [
        C64::from_re(core::f64::consts::FRAC_1_SQRT_2),
        C64::new(0.0, core::f64::consts::FRAC_1_SQRT_2),
    ] {
        let rotated: Vec<_> = values
            .iter()
            .map(|v| v.iter().map(|x| x.map(|z| z * multiplier)).collect())
            .collect();
        accumulate::<RealKernel>(
            &rotated,
            &omega,
            &kx,
            aperture,
            window,
            &mut sum,
            &mut || false,
        )
        .unwrap();
    }
    let split = sum.finish(aperture).unwrap();
    close(
        split.mean_field_energy_density_j_m3(),
        proper.mean_field_energy_density_j_m3(),
        2e-14,
    );
    for i in 0..3 {
        for j in 0..3 {
            close(
                split.mean_maxwell_stress_pa()[i][j],
                proper.mean_maxwell_stress_pa()[i][j],
                2e-14,
            );
        }
    }
    let mut independent = Stress::new(2.0, 8.0);
    for sign in [1.0, -1.0] {
        let v: Vec<_> = values
            .iter()
            .map(|v| {
                v.iter()
                    .map(|x| x.map(|z| z * C64::from_re(sign)))
                    .collect()
            })
            .collect();
        accumulate::<ProperKernel>(
            &v,
            &omega,
            &kx,
            aperture,
            window,
            &mut independent,
            &mut || false,
        )
        .unwrap();
    }
    close(
        independent
            .finish(aperture)
            .unwrap()
            .mean_field_energy_density_j_m3(),
        2.0 * proper.mean_field_energy_density_j_m3(),
        2e-14,
    );
    let dark = vec![vec![[C64::ZERO; 6]; 2]; 2];
    assert_eq!(
        observe::<ProperKernel>(&dark, &omega, &kx, aperture, window)
            .mean_field_energy_density_j_m3(),
        0.0
    );
}

#[test]
fn a_narrow_dark_fringe_retains_positive_finite_area_stress() {
    let width = 1e-8;
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [width, 1.0]).unwrap();
    let window = HarmonicWindow::new(0.0, 0.2).unwrap();
    let a = [1.0, 0.0, 0.0, 0.0, 0.5, 0.0].map(C64::from_re);
    let values = vec![vec![a, a.map(|z| -z)]];
    let result = observe::<ProperKernel>(&values, &[3.0], &[0.0, 1.0], aperture, window);
    // epsilon*<|1-exp(-ix)|²>/2 = width²/3 + O(width^4).
    let expected = width * width / 3.0;
    assert!((result.mean_field_energy_density_j_m3() / expected - 1.0).abs() < 1e-12);
    assert!((result.mean_maxwell_stress_pa()[2][2] / (-expected) - 1.0).abs() < 1e-12);
}

#[test]
fn work_limits_cancellation_replay_and_derived_unit_overflow_refuse() {
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let limits = AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 16,
        maximum_workspace_bytes: 1 << 20,
    };
    assert!(admit::<Stress>(grid, 4, 0, true, true, limits).is_ok());
    assert!(
        admit::<Stress>(
            grid,
            4,
            0,
            true,
            true,
            AngularSpectrumExposureApertureLimits {
                maximum_mode_frequency_pairs: 15,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        admit::<Stress>(
            grid,
            4,
            0,
            true,
            true,
            AngularSpectrumExposureApertureLimits {
                maximum_workspace_bytes: 0,
                ..limits
            }
        )
        .is_err()
    );
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [0.5; 2]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.42).unwrap();
    let mut polls = 0;
    let mut scratch = Stress::new(2.0, 8.0);
    assert!(
        accumulate::<RealKernel>(
            &scene(),
            &[3.2, 4.7],
            &[0.0, 1.7],
            aperture,
            window,
            &mut scratch,
            &mut || {
                polls += 1;
                polls == 7
            }
        )
        .is_err()
    );
    let a = observe::<RealKernel>(&scene(), &[3.2, 4.7], &[0.0, 1.7], aperture, window);
    let b = observe::<RealKernel>(&scene(), &[3.2, 4.7], &[0.0, 1.7], aperture, window);
    assert_eq!(a, b);
    let huge = AngularSpectrumAperture::try_new([0.0; 2], [1e154; 2]).unwrap();
    let mut sum = Stress::new(2.0, 8.0);
    sum.dyad = [25.0, 25.0, 0.0, 0.0, 0.0, 0.0].map(Dd::from_f64);
    sum.poynting[2] = Dd::from_f64(6.25);
    let reading = sum.finish(huge).unwrap();
    assert!(reading.integrated_stress_n().is_err());
    assert!(reading.integrated_traction_n([0.0, 0.0, 1.0]).is_err());
    assert!(reading.normal_power_w().is_err());
}

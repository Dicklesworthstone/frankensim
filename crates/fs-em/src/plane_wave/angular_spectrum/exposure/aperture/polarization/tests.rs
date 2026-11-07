use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn observe<K: Kernel, A: Accumulator>(
    kx: &[f64],
    ky: &[f64],
    omega: &[f64],
    fields: &[Vec<FieldParts>],
    area: AngularSpectrumAperture,
    window: HarmonicWindow,
    mut sum: A,
) -> A::Output {
    let grid = AngularSpectrumGrid {
        width: kx.len(),
        height: ky.len(),
        spacing_m: [1.0; 2],
    };
    let spatial = Spatial::<K>::from_axes(grid, area, |i| kx[i], |i| ky[i], &mut || false).unwrap();
    let temporal = pair_kernels::<K>(omega.len(), window, |i| -omega[i], &mut || false).unwrap();
    let mut fields = fields.to_vec();
    for (i, field) in fields.iter_mut().enumerate() {
        sum.prepare_field(field, &mut || false).unwrap();
        if K::CENTER_FIELDS {
            center_fields(
                field,
                omega[i],
                |m| kx[m % kx.len()],
                |m| ky[m / kx.len()],
                area,
                window,
                &mut || false,
            )
            .unwrap();
        }
    }
    integrate(&spatial, &temporal, &fields, &mut sum, &mut || false).unwrap();
    sum.finish(area).unwrap()
}

fn parts(fields: &[Vec<[C64; 6]>]) -> Vec<Vec<FieldParts>> {
    fields
        .iter()
        .map(|f| {
            f.iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect()
        })
        .collect()
}

fn reading_values(reading: AngularSpectrumExposureApertureReading) -> [f64; 10] {
    let p = reading.mean;
    [
        p.electric_mean_squared_v2_m2[0],
        p.electric_mean_squared_v2_m2[1],
        p.electric_mean_squared_v2_m2[2],
        p.magnetic_mean_squared_a2_m2[0],
        p.magnetic_mean_squared_a2_m2[1],
        p.magnetic_mean_squared_a2_m2[2],
        p.poynting_w_m2[0],
        p.poynting_w_m2[1],
        p.poynting_w_m2[2],
        reading.normal_power_w,
    ]
}

#[test]
fn all_real_correlations_match_independent_space_time_and_phase_quadrature() {
    let kx = [0.2, -0.7];
    let ky = [0.3, -0.5];
    let omega = [1.3, 2.1];
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
    let fields_dd = parts(&fields);
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    let projection =
        AngularSpectrumExposureApertureProjection::try_new([0.6, -0.8, 0.0], [0.3, 0.4, -0.5])
            .unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for proper in [false, true] {
            let (full, direct, basic) = if proper {
                (
                    observe::<ProperKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        Full([Dd::ZERO; 21]),
                    ),
                    observe::<ProperKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        Projected::new(projection),
                    ),
                    observe::<ProperKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        [Dd::ZERO; 9],
                    ),
                )
            } else {
                (
                    observe::<RealKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        Full([Dd::ZERO; 21]),
                    ),
                    observe::<RealKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        Projected::new(projection),
                    ),
                    observe::<RealKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields_dd,
                        area,
                        window,
                        [Dd::ZERO; 9],
                    ),
                )
            };
            let expected = super::super::tests::reference_matrix(
                &kx, &ky, &omega, &fields, area, window, proper,
            );
            for (matrix, row, column) in [
                (full.electric_second_moment_v2_m2(), 0, 0),
                (full.magnetic_second_moment_a2_m2(), 3, 3),
                (full.electric_magnetic_second_moment_w_m2(), 0, 3),
            ] {
                for i in 0..3 {
                    for j in 0..3 {
                        close(
                            matrix[i][j],
                            expected[row + i][column + j],
                            3e-7 * matrix[i][j].abs().max(1.0),
                        );
                        if row == column {
                            assert_eq!(matrix[i][j], matrix[j][i]);
                        }
                    }
                }
            }
            for (a, b) in reading_values(full.reading())
                .into_iter()
                .zip(reading_values(basic))
            {
                close(a, b, 2e-14);
            }
            let contracted = full.project(projection).unwrap();
            close(
                direct.electric_mean_squared_v2_m2,
                contracted.electric_mean_squared_v2_m2,
                2e-14,
            );
            close(
                direct.magnetic_mean_squared_a2_m2,
                contracted.magnetic_mean_squared_a2_m2,
                2e-14,
            );
            close(
                direct.electric_magnetic_mean_product_w_m2,
                contracted.electric_magnetic_mean_product_w_m2,
                2e-14,
            );
        }
    }
}

#[test]
fn real_linear_analyzers_obey_malus_scaling_and_proper_circular_isotropy() {
    let phi: f64 = 0.37;
    let (s, c) = phi.sin_cos();
    let linear = parts(&[vec![[
        C64::from_re(c),
        C64::from_re(s),
        C64::ZERO,
        C64::from_re(-s),
        C64::from_re(c),
        C64::ZERO,
    ]]]);
    let circular = parts(&[vec![[
        C64::ONE,
        C64::new(0.0, 1.0),
        C64::ZERO,
        C64::new(0.0, -1.0),
        C64::ONE,
        C64::ZERO,
    ]]]);
    let area = AngularSpectrumAperture::try_new([-0.3, 0.2], [0.7, 0.4]).unwrap();
    let cycle = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let subcycle =
        HarmonicWindow::with_profile(0.31, 0.07, HarmonicWindowProfile::Triangular).unwrap();
    for theta in [0.0_f64, 0.53, 1.19] {
        let (s, c) = theta.sin_cos();
        let projection = AngularSpectrumExposureApertureProjection::try_new(
            [2.0 * c, 2.0 * s, 0.0],
            [s, -c, 0.0],
        )
        .unwrap();
        let actual = observe::<RealKernel, _>(
            &[0.0],
            &[0.0],
            &[1.0],
            &linear,
            area,
            cycle,
            Projected::new(projection),
        );
        let malus = 0.5 * (theta - phi).cos().powi(2);
        close(actual.electric_mean_squared_v2_m2, 4.0 * malus, 2e-15);
        close(actual.magnetic_mean_squared_a2_m2, malus, 2e-15);
        close(
            actual.electric_magnetic_mean_product_w_m2,
            -2.0 * malus,
            2e-15,
        );
        let actual = observe::<ProperKernel, _>(
            &[0.0],
            &[0.0],
            &[1.0],
            &circular,
            area,
            subcycle,
            Projected::new(projection),
        );
        close(actual.electric_mean_squared_v2_m2, 2.0, 2e-15);
        close(actual.magnetic_mean_squared_a2_m2, 0.5, 2e-15);
        close(actual.electric_magnetic_mean_product_w_m2, -1.0, 2e-15);
    }
}

#[test]
fn full_field_tensors_transform_with_a_rotated_cartesian_basis() {
    let fields = [vec![[
        C64::new(0.7, -0.2),
        C64::new(-0.3, 0.5),
        C64::new(0.2, 0.4),
        C64::new(-0.6, 0.1),
        C64::new(0.4, 0.3),
        C64::new(0.5, -0.7),
    ]]];
    let (s, c) = 0.41_f64.sin_cos();
    let (u, v) = 0.73_f64.sin_cos();
    let rotation = [[c, -s * v, s * u], [s, c * v, -c * u], [0.0, u, v]];
    let rotated = [vec![core::array::from_fn(|i| {
        let offset = 3 * (i / 3);
        (0..3).fold(C64::ZERO, |a, j| {
            a + fields[0][0][offset + j].scale(rotation[i % 3][j])
        })
    })]];
    let area = AngularSpectrumAperture::try_new([-0.2, 0.1], [0.7, 0.3]).unwrap();
    let window = HarmonicWindow::new(0.17, 0.29).unwrap();
    let a = observe::<RealKernel, _>(
        &[0.4],
        &[-0.3],
        &[1.7],
        &parts(&fields),
        area,
        window,
        Full([Dd::ZERO; 21]),
    );
    let b = observe::<RealKernel, _>(
        &[0.4],
        &[-0.3],
        &[1.7],
        &parts(&rotated),
        area,
        window,
        Full([Dd::ZERO; 21]),
    );
    for (old, new) in [
        (
            a.electric_second_moment_v2_m2(),
            b.electric_second_moment_v2_m2(),
        ),
        (
            a.magnetic_second_moment_a2_m2(),
            b.magnetic_second_moment_a2_m2(),
        ),
        (
            a.electric_magnetic_second_moment_w_m2(),
            b.electric_magnetic_second_moment_w_m2(),
        ),
    ] {
        for i in 0..3 {
            for j in 0..3 {
                let mut expected = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        expected += rotation[i][k] * old[k][l] * rotation[j][l];
                    }
                }
                close(new[i][j], expected, 4e-15);
            }
        }
    }
    let flux = a.reading().mean.poynting_w_m2;
    for i in 0..3 {
        close(
            b.reading().mean.poynting_w_m2[i],
            (0..3).map(|j| rotation[i][j] * flux[j]).sum(),
            4e-15,
        );
    }
}

#[test]
fn direct_projection_preserves_grouped_weak_residual_before_quadratic_cancellation() {
    let tiny = 1e-20;
    // Original equal-frequency line summation can retain these low parts:
    // a strong [1,1] field plus a weak [tiny,0] field. A crossed real analyzer
    // must project the grouped amplitudes before forming their second moments.
    let mut field = [[Dd::ZERO; 2]; 6];
    field[0][0] = Dd::from_pair(1.0, tiny);
    field[1][0] = Dd::ONE;
    field[3][0] = Dd::from_pair(1.0, 2.0 * tiny);
    field[4][0] = Dd::ONE;
    let projection =
        AngularSpectrumExposureApertureProjection::try_new([1.0, -1.0, 0.0], [1.0, -1.0, 0.0])
            .unwrap();
    let area = AngularSpectrumAperture::try_new([1.0, -0.3], [1e-16, 0.7]).unwrap();
    let time = HarmonicWindow::new(0.27, core::f64::consts::TAU).unwrap();
    for proper in [false, true] {
        let value = if proper {
            observe::<ProperKernel, _>(
                &[0.0],
                &[0.0],
                &[1.0],
                &[vec![field]],
                area,
                time,
                Projected::new(projection),
            )
        } else {
            observe::<RealKernel, _>(
                &[0.0],
                &[0.0],
                &[1.0],
                &[vec![field]],
                area,
                time,
                Projected::new(projection),
            )
        };
        close(
            value.electric_mean_squared_v2_m2 / (0.5 * tiny * tiny),
            1.0,
            3e-14,
        );
        close(
            value.magnetic_mean_squared_a2_m2 / (2.0 * tiny * tiny),
            1.0,
            3e-14,
        );
        close(
            value.electric_magnetic_mean_product_w_m2 / (tiny * tiny),
            1.0,
            3e-14,
        );
    }
}

#[test]
fn correlation_and_projection_limits_and_cancellation_are_explicit() {
    assert!(
        AngularSpectrumExposureApertureProjection::try_new([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err()
    );
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 2,
        spacing_m: [1.0; 2],
    };
    let limits = AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 64,
        maximum_workspace_bytes: 0,
    };
    assert!(admit::<Full>(grid, 4, 128, true, true, limits).is_err());
    assert!(admit::<Projected>(grid, 4, 128, true, true, limits).is_err());
    let projection =
        AngularSpectrumExposureApertureProjection::try_new([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let projected = Projected::new(projection);
    assert!(matches!(
        projected.prepare_field(&mut [[[Dd::ZERO; 2]; 6]], &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let overflow = Projected::new(
        AngularSpectrumExposureApertureProjection::try_new([2.0, 0.0, 0.0], [0.0; 3]).unwrap(),
    );
    let mut huge = [[[Dd::ZERO; 2]; 6]];
    huge[0][0][0] = Dd::from_f64(f64::MAX);
    assert!(overflow.prepare_field(&mut huge, &mut || false).is_err());
}

#[test]
fn jones_analyzers_resolve_circular_handedness_without_cycle_averaging() {
    let field = parts(&[vec![[
        C64::ONE,
        C64::new(0.0, 1.0),
        C64::ZERO,
        C64::new(0.0, -1.0),
        C64::ONE,
        C64::ZERO,
    ]]]);
    let area = AngularSpectrumAperture::try_new([0.0; 2], [0.7, 0.4]).unwrap();
    let window = HarmonicWindow::new(0.31, 0.07).unwrap();
    let r = core::f64::consts::FRAC_1_SQRT_2;
    for handedness in [-1.0, 1.0] {
        let weights = [C64::from_re(r), C64::new(0.0, handedness * r), C64::ZERO];
        let analyzer =
            AngularSpectrumExposureApertureJonesProjection::try_new(weights, weights).unwrap();
        let actual = observe::<RealKernel, _>(
            &[0.0],
            &[0.0],
            &[1.0],
            &field,
            area,
            window,
            Projected::new(analyzer),
        );
        let proper = observe::<ProperKernel, _>(
            &[0.0],
            &[0.0],
            &[1.0],
            &field,
            area,
            window,
            Projected::new(analyzer),
        );
        let (electric, magnetic, mixed, cycle) = if handedness > 0.0 {
            // a†E=sqrt(2), a†H=-i sqrt(2). The REAL signals are
            // sqrt(2) cos(t), sqrt(2) sin(t); retain the subcycle carrier.
            let start = window.start_s();
            let end = start + window.duration_s();
            let cosine = ((2.0 * end).sin() - (2.0 * start).sin()) / (2.0 * window.duration_s());
            let sine = ((2.0 * start).cos() - (2.0 * end).cos()) / (2.0 * window.duration_s());
            (1.0 + cosine, 1.0 - cosine, sine, 1.0)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };
        close(actual.electric_mean_squared_v2_m2, electric, 3e-15);
        close(actual.magnetic_mean_squared_a2_m2, magnetic, 3e-15);
        close(actual.electric_magnetic_mean_product_w_m2, mixed, 3e-15);
        close(proper.electric_mean_squared_v2_m2, cycle, 3e-15);
        close(proper.magnetic_mean_squared_a2_m2, cycle, 3e-15);
        close(proper.electric_magnetic_mean_product_w_m2, 0.0, 3e-15);
    }
}

#[test]
fn complex_elliptical_analyzers_match_space_time_phase_quadrature() {
    let kx = [0.2, -0.7];
    let ky = [0.3];
    let omega = [1.3, 2.1];
    let fields: Vec<Vec<[C64; 6]>> = (0..2)
        .map(|f| {
            (0..2)
                .map(|m| {
                    core::array::from_fn(|c| {
                        C64::new(
                            (1 + f + c) as f64 / (7.0 * (m + 1) as f64),
                            (c as f64 - 2.0) * (m as f64 - 1.5) / (11.0 * (f + 1) as f64),
                        )
                    })
                })
                .collect()
        })
        .collect();
    let electric = [C64::new(0.3, -0.7), C64::new(-0.4, 0.2), C64::new(0.1, 0.3)];
    let magnetic = [C64::new(-0.2, 0.1), C64::new(0.6, -0.4), C64::new(0.5, 0.2)];
    let analyzer =
        AngularSpectrumExposureApertureJonesProjection::try_new(electric, magnetic).unwrap();
    assert_eq!(analyzer.electric_weights(), electric);
    assert_eq!(analyzer.magnetic_weights(), magnetic);
    let projected: Vec<Vec<[C64; 6]>> = fields
        .iter()
        .map(|line| {
            line.iter()
                .map(|v| {
                    // Independent ordinary-complex projection; the oracle integrates the
                    // actual real signals over space, time and four proper-factor phases.
                    let mut out = [C64::ZERO; 6];
                    for c in 0..3 {
                        out[0] = out[0] + electric[c].conj() * v[c];
                        out[3] = out[3] + magnetic[c].conj() * v[c + 3];
                    }
                    out
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
            let actual = if proper {
                observe::<ProperKernel, _>(
                    &kx,
                    &ky,
                    &omega,
                    &parts(&fields),
                    area,
                    window,
                    Projected::new(analyzer),
                )
            } else {
                observe::<RealKernel, _>(
                    &kx,
                    &ky,
                    &omega,
                    &parts(&fields),
                    area,
                    window,
                    Projected::new(analyzer),
                )
            };
            let expected = super::super::tests::reference_matrix(
                &kx, &ky, &omega, &projected, area, window, proper,
            );
            for (value, expected) in [
                (actual.electric_mean_squared_v2_m2, expected[0][0]),
                (actual.magnetic_mean_squared_a2_m2, expected[3][3]),
                (actual.electric_magnetic_mean_product_w_m2, expected[0][3]),
            ] {
                // Match the scale-aware Simpson tolerance of the original
                // full-moment oracle; these authored weights need not be unit.
                close(value, expected, 3e-7 * expected.abs().max(1.0));
            }
        }
    }
}

#[test]
fn crossed_circular_analyzer_preserves_dd_residuals_and_refuses_invalid_weights() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            AngularSpectrumExposureApertureJonesProjection::try_new(
                [C64::new(0.0, value), C64::ZERO, C64::ZERO],
                [C64::ZERO; 3],
            )
            .is_err()
        );
    }
    let analyzer = AngularSpectrumExposureApertureJonesProjection::try_new(
        [C64::ONE, C64::new(0.0, -1.0), C64::ZERO],
        [C64::ZERO; 3],
    )
    .unwrap();
    let tiny = 1e-20;
    let mut field = [[Dd::ZERO; 2]; 6];
    field[0][0] = Dd::from_pair(1.0, tiny);
    field[1][1] = Dd::ONE;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let actual = observe::<RealKernel, _>(
        &[0.0],
        &[0.0],
        &[1.0],
        &[vec![field]],
        area,
        window,
        Projected::new(analyzer),
    );
    close(
        actual.electric_mean_squared_v2_m2 / (0.5 * tiny * tiny),
        1.0,
        3e-14,
    );
    assert_eq!(actual.magnetic_mean_squared_a2_m2, 0.0);
    assert!(matches!(
        Projected::new(analyzer).prepare_field(&mut [field], &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
}

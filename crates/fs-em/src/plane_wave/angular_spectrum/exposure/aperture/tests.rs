use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn grid(width: usize, height: usize) -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width,
        height,
        spacing_m: [1.0; 2],
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn evaluate<K: Kernel>(
    kx: &[f64],
    ky: &[f64],
    omega: &[f64],
    fields: &[Vec<[C64; 6]>],
    area: AngularSpectrumAperture,
    time: HarmonicWindow,
) -> [f64; 9] {
    let spatial = Spatial::<K>::from_axes(
        grid(kx.len(), ky.len()),
        area,
        |i| kx[i],
        |i| ky[i],
        &mut || false,
    )
    .unwrap();
    let temporal = pair_kernels::<K>(omega.len(), time, |i| -omega[i], &mut || false).unwrap();
    let mut fields: Vec<Vec<FieldParts>> = fields
        .iter()
        .map(|f| {
            f.iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect()
        })
        .collect();
    if K::CENTER_FIELDS {
        for (f, field) in fields.iter_mut().enumerate() {
            center_fields(
                field,
                omega[f],
                |i| kx[i % kx.len()],
                |i| ky[i / kx.len()],
                area,
                time,
                &mut || false,
            )
            .unwrap();
        }
    }
    let mut sum = [Dd::ZERO; 9];
    integrate(&spatial, &temporal, &fields, &mut sum, &mut || false).unwrap();
    sum.map(Dd::to_f64)
}

fn wave(amplitude: C64, backwards: bool) -> [C64; 6] {
    [
        amplitude,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        amplitude.scale(if backwards { -1.0 } else { 1.0 }),
        C64::ZERO,
    ]
}

// Independent three-dimensional Simpson integration of the actual real field.
// Proper factors are averaged over z={1,-1,i,-i}, whose joint moments are
// exactly E[z]=E[z²]=0 and E[|z|²]=1. No analytic kernel is used by this oracle.
pub(super) fn reference_matrix(
    kx: &[f64],
    ky: &[f64],
    omega: &[f64],
    fields: &[Vec<[C64; 6]>],
    area: AngularSpectrumAperture,
    time: HarmonicWindow,
    proper: bool,
) -> [[f64; 6]; 6] {
    let steps = 24;
    let phases = if proper {
        vec![
            C64::ONE,
            C64::from_re(-1.0),
            C64::new(0.0, 1.0),
            C64::new(0.0, -1.0),
        ]
    } else {
        vec![C64::ONE]
    };
    let weight = |i| {
        if i == 0 || i == steps {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let mut sum = [[0.0; 6]; 6];
    for ix in 0..=steps {
        let x = area.origin_m()[0] + area.extent_m()[0] * ix as f64 / steps as f64;
        for iy in 0..=steps {
            let y = area.origin_m()[1] + area.extent_m()[1] * iy as f64 / steps as f64;
            for it in 0..=steps {
                let u = it as f64 / steps as f64;
                let t = time.start_s() + time.duration_s() * u;
                let sensitivity = match time.profile() {
                    HarmonicWindowProfile::Uniform => 1.0,
                    HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
                };
                let w = weight(ix) * weight(iy) * weight(it) * sensitivity;
                let mut phasor = [C64::ZERO; 6];
                for (f, field) in fields.iter().enumerate() {
                    for (m, values) in field.iter().enumerate() {
                        let angle = omega[f] * t - kx[m % kx.len()] * x - ky[m / kx.len()] * y;
                        let (s, c) = angle.sin_cos();
                        for component in 0..6 {
                            phasor[component] =
                                phasor[component] + values[component] * C64::new(c, s);
                        }
                    }
                }
                for z in &phases {
                    let real = phasor.map(|v| (v * *z).re);
                    for i in 0..6 {
                        for j in 0..6 {
                            sum[i][j] += w * real[i] * real[j];
                        }
                    }
                }
            }
        }
    }
    let normalization = (3.0 * steps as f64).powi(3) * phases.len() as f64;
    sum.map(|row| row.map(|v| v / normalization))
}

fn reference(
    kx: &[f64],
    ky: &[f64],
    omega: &[f64],
    fields: &[Vec<[C64; 6]>],
    area: AngularSpectrumAperture,
    time: HarmonicWindow,
    proper: bool,
) -> [f64; 9] {
    let matrix = reference_matrix(kx, ky, omega, fields, area, time, proper);
    let mut result = [0.0; 9];
    for c in 0..3 {
        result[c] = matrix[c][c];
        result[c + 3] = matrix[c + 3][c + 3];
        let (a, b) = ((c + 1) % 3, (c + 2) % 3);
        result[c + 6] = matrix[a][b + 3] - matrix[b][a + 3];
    }
    result
}

#[test]
fn vector_moments_match_independent_space_time_and_phase_quadrature() {
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
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let time = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for proper in [false, true] {
            let actual = if proper {
                evaluate::<ProperKernel>(&kx, &ky, &omega, &fields, area, time)
            } else {
                evaluate::<RealKernel>(&kx, &ky, &omega, &fields, area, time)
            };
            let expected = reference(&kx, &ky, &omega, &fields, area, time, proper);
            for c in 0..9 {
                close(actual[c], expected[c], 2e-7 * expected[c].abs().max(1.0));
            }
        }
    }
}

#[test]
fn real_gram_preserves_joint_subcycle_and_subfringe_nulls() {
    let fields = [vec![wave(C64::new(0.0, 1.0), true)]];
    let tiny = 1e-20;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [tiny, 1.0]).unwrap();
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 1.0 / 6.0),
        (HarmonicWindowProfile::Triangular, 1.0 / 8.0),
    ] {
        let time = HarmonicWindow::with_profile(0.0, tiny, profile).unwrap();
        // E=sin(2*x-2*t). Its squared linear term has zero mean phase and
        // variance 4*(Var(x)+Var(t)); the next term is O(tiny^4).
        let actual = evaluate::<RealKernel>(&[2.0], &[0.0], &[2.0], &fields, area, time);
        let expected = 4.0 * tiny * tiny * coefficient;
        close(actual[0] / expected, 1.0, 2e-14);
        close(actual[4] / expected, 1.0, 2e-14);
        close(actual[8] / expected, -1.0, 2e-14);
        let proper = evaluate::<ProperKernel>(&[2.0], &[0.0], &[2.0], &fields, area, time);
        close(proper[0], 0.5, 1e-15);
        close(proper[8], -0.5, 1e-15);
    }
}

fn translated_null(origin: f64, width: f64) {
    let fields = [vec![wave(C64::new(0.0, 1.0), true)]];
    let area = AngularSpectrumAperture::try_new([origin, 0.0], [width, 1.0]).unwrap();
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 2.0 / 3.0),
        (HarmonicWindowProfile::Triangular, 0.5),
    ] {
        let time = HarmonicWindow::with_profile(origin, width, profile).unwrap();
        // The common translation leaves E=sin(2*x-2*t) unchanged. Its leading
        // mean square is 4*(Var(x)+Var(t)), without a center-phase contribution.
        let expected = coefficient * width * width;
        let values = evaluate::<RealKernel>(&[2.0], &[0.0], &[2.0], &fields, area, time);
        close(values[0] / expected, 1.0, 3e-14);
        close(values[4] / expected, 1.0, 3e-14);
        close(values[8] / expected, -1.0, 3e-14);
    }
}

#[test]
fn translated_null_does_not_round_positive_power_to_darkness() {
    // Uncentered Gram composition returned exactly zero for this observation.
    translated_null(1.0, 1e-16);
}

#[test]
fn translated_null_does_not_create_negative_mean_square_or_refuse() {
    // Uncentered Gram composition returned -6.16e-33 instead of +6.67e-41.
    translated_null(0.25, 1e-20);
    // Unequal very short intervals also retain their distinct midpoint phases
    // even when neither half-width fits beside the common start in DD.
    let tiny = 1e-50;
    let fields = [vec![wave(C64::new(0.0, 1.0), false)]];
    let area = AngularSpectrumAperture::try_new([1.0, 0.0], [2.0 * tiny, 1.0]).unwrap();
    let time = HarmonicWindow::new(1.0, tiny).unwrap();
    let actual = evaluate::<RealKernel>(&[2.0], &[0.0], &[2.0], &fields, area, time);
    close(actual[0] / ((8.0 / 3.0) * tiny * tiny), 1.0, 3e-14);
}

#[test]
fn proper_fringe_and_beat_washout_are_joint_not_separate_power_averages() {
    let fields = [
        vec![wave(C64::ONE, false), [C64::ZERO; 6]],
        vec![[C64::ZERO; 6], wave(C64::ONE, false)],
    ];
    let sinc = |x: f64| x.sin() / x;
    for (width, duration) in [
        (1.2, 0.9),
        (core::f64::consts::TAU, 0.9),
        (1.2, core::f64::consts::TAU),
    ] {
        let area = AngularSpectrumAperture::try_new([-0.4, 0.2], [width, 0.3]).unwrap();
        let time = HarmonicWindow::new(0.7, duration).unwrap();
        let expected = 1.0
            + sinc(width / 2.0)
                * sinc(duration / 2.0)
                * ((0.7 + duration / 2.0) - (-0.4 + width / 2.0)).cos();
        let actual =
            evaluate::<ProperKernel>(&[0.0, 1.0], &[0.0], &[2.0, 3.0], &fields, area, time);
        close(actual[0], expected, 2e-15);
        close(actual[8], expected, 2e-15);
        let rotated: Vec<_> = fields
            .iter()
            .map(|f| {
                f.iter()
                    .map(|v| v.map(|z| z * C64::new(0.0, 1.0)))
                    .collect()
            })
            .collect();
        let phase_invariant =
            evaluate::<ProperKernel>(&[0.0, 1.0], &[0.0], &[2.0, 3.0], &rotated, area, time);
        assert_eq!(actual, phase_invariant);
    }
}

#[test]
fn aperture_partition_preserves_all_moments_and_full_cell_nyquist_mean() {
    let kx = [0.0, -core::f64::consts::PI];
    let fields = [vec![
        wave(C64::new(0.2, 0.3), false),
        wave(C64::new(-0.7, 0.4), true),
    ]];
    let time = HarmonicWindow::new(0.13, 0.29).unwrap();
    let whole = AngularSpectrumAperture::try_new([0.0; 2], [2.0, 1.0]).unwrap();
    let left = AngularSpectrumAperture::try_new([0.0; 2], [0.3, 1.0]).unwrap();
    let right = AngularSpectrumAperture::try_new([0.3, 0.0], [1.7, 1.0]).unwrap();
    for proper in [false, true] {
        let observe = |area| {
            if proper {
                evaluate::<ProperKernel>(&kx, &[0.0], &[2.0], &fields, area, time)
            } else {
                evaluate::<RealKernel>(&kx, &[0.0], &[2.0], &fields, area, time)
            }
        };
        let (all, a, b) = (observe(whole), observe(left), observe(right));
        for c in 0..9 {
            close(all[c], (0.3 * a[c] + 1.7 * b[c]) / 2.0, 3e-15);
        }
        let dc = if proper {
            0.5 * (0.2_f64.powi(2) + 0.3_f64.powi(2))
        } else {
            time.product_kernel(-2.0, -2.0)
                .unwrap()
                .mean_real_product(C64::new(0.2, 0.3), C64::new(0.2, 0.3))
                .unwrap()
        };
        close(
            all[0],
            dc + 0.5 * (0.7_f64.powi(2) + 0.4_f64.powi(2)),
            3e-15,
        );
        close(
            all[8],
            dc - 0.5 * (0.7_f64.powi(2) + 0.4_f64.powi(2)),
            3e-15,
        );
    }
}

#[test]
fn query_budgets_overflow_and_cancellation_are_bounded() {
    let limits = AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 64,
        maximum_workspace_bytes: 1 << 20,
    };
    assert!(admit::<[Dd; 9]>(grid(2, 2), 4, 128, true, true, limits).is_ok());
    assert!(
        admit::<[Dd; 9]>(
            grid(2, 2),
            4,
            128,
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
        admit::<[Dd; 9]>(
            grid(2, 2),
            4,
            128,
            true,
            true,
            AngularSpectrumExposureApertureLimits {
                maximum_workspace_bytes: 0,
                ..limits
            }
        )
        .is_err()
    );
    assert!(admit::<[Dd; 9]>(grid(usize::MAX, 2), 4, 128, true, true, limits).is_err());
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    assert!(matches!(
        Spatial::<RealKernel>::from_axes(grid(2, 1), area, |i| i as f64, |_| 0.0, &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let spatial =
        Spatial::<RealKernel>::from_axes(grid(2, 1), area, |i| i as f64, |_| 0.0, &mut || false)
            .unwrap();
    let temporal = pair_kernels::<RealKernel>(
        1,
        HarmonicWindow::new(0.0, 1.0).unwrap(),
        |_| 1.0,
        &mut || false,
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        integrate(
            &spatial,
            &temporal,
            &[vec![
                wave(C64::ONE, false)
                    .map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]);
                2
            ]],
            &mut [Dd::ZERO; 9],
            &mut || {
                polls += 1;
                polls == 3
            }
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}

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

fn observe<K: Kernel>(
    kx: &[f64],
    ky: &[f64],
    frequencies: &[f64],
    fields: &[Vec<[C64; 6]>],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> [f64; 9] {
    let spatial = Spatial::<K>::from_axes(
        grid(kx.len(), ky.len()),
        aperture,
        |i| kx[i],
        |i| ky[i],
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
    let mut headroom = usize::MAX;
    let modes = admit(
        window,
        frequencies,
        kx.len() * ky.len(),
        !K::CENTER_FIELDS,
        &mut headroom,
        usize::MAX,
    )
    .unwrap()
    .unwrap();
    let mut sum = [Dd::ZERO; 9];
    integrate(
        &spatial,
        window,
        frequencies,
        &parts,
        modes,
        &mut sum,
        &mut || false,
    )
    .unwrap();
    sum.map(Dd::to_f64)
}

#[test]
fn real_two_tone_dark_residual_is_formed_before_power() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let fields = [vec![wave(C64::ONE)], vec![wave(C64::from_re(-1.0))]];
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 9.0 / 320.0),
        (HarmonicWindowProfile::Triangular, 9.0 / 960.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let actual = observe::<RealKernel>(&[0.0], &[0.0], &[1.0, 2.0], &fields, area, window);
        // cos(t)-cos(2t) = 3*t²/2 + O(t⁴), squared before the exact
        // uniform/triangular centered fourth moment is applied.
        let expected = coefficient * duration.powi(4);
        close(actual[0] / expected, 1.0, 2e-12);
        close(actual[4] / expected, 1.0, 2e-12);
        close(actual[8] / expected, -1.0, 2e-12);
    }
}

#[test]
fn proper_three_tone_null_retains_optical_carrier_invariance_and_half_power() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let fields: Vec<_> = [1.0, -2.0, 1.0].map(|a| vec![wave(C64::from_re(a))]).into();
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 1.0 / 160.0),
        (HarmonicWindowProfile::Triangular, 1.0 / 480.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let low = observe::<ProperKernel>(&[0.0], &[0.0], &[1.0, 2.0, 3.0], &fields, area, window);
        let optical = observe::<ProperKernel>(
            &[0.0],
            &[0.0],
            &[1e15, 1e15 + 1.0, 1e15 + 2.0],
            &fields,
            area,
            window,
        );
        assert_eq!(low, optical);
        let expected = coefficient * duration.powi(4);
        close(optical[0] / expected, 1.0, 2e-12);
        close(optical[8] / expected, -1.0, 2e-12);
    }
}

#[test]
fn continuous_spatial_and_temporal_parities_match_independent_field_quadrature() {
    let (kx, ky, frequencies) = ([0.2, -0.7], [0.3, -0.5], [1.3, 2.1]);
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
            let actual = if proper {
                observe::<ProperKernel>(&kx, &ky, &frequencies, &fields, area, window)
            } else {
                observe::<RealKernel>(&kx, &ky, &frequencies, &fields, area, window)
            };
            let matrix = super::super::tests::reference_matrix(
                &kx,
                &ky,
                &frequencies,
                &fields,
                area,
                window,
                proper,
            );
            for c in 0..3 {
                close(actual[c], matrix[c][c], 2e-7 * matrix[c][c].abs().max(1.0));
                close(
                    actual[c + 3],
                    matrix[c + 3][c + 3],
                    2e-7 * matrix[c + 3][c + 3].abs().max(1.0),
                );
                let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                close(actual[c + 6], matrix[u][v + 3] - matrix[v][u + 3], 2e-7);
            }
        }
    }
}

#[test]
fn unpaired_nyquist_and_rectangle_partition_remain_continuous() {
    let kx = [0.0, -core::f64::consts::PI];
    let fields = [
        vec![wave(C64::new(0.2, 0.3)), wave(C64::new(-0.7, 0.4))],
        vec![wave(C64::new(-0.1, 0.1)), wave(C64::new(0.2, -0.3))],
    ];
    let window = HarmonicWindow::new(0.13, 0.29).unwrap();
    let whole = AngularSpectrumAperture::try_new([0.0; 2], [2.0, 1.0]).unwrap();
    let left = AngularSpectrumAperture::try_new([0.0; 2], [0.3, 1.0]).unwrap();
    let right = AngularSpectrumAperture::try_new([0.3, 0.0], [1.7, 1.0]).unwrap();
    for proper in [false, true] {
        let eval = |area| {
            if proper {
                observe::<ProperKernel>(&kx, &[0.0], &[1.0, 2.0], &fields, area, window)
            } else {
                observe::<RealKernel>(&kx, &[0.0], &[1.0, 2.0], &fields, area, window)
            }
        };
        let (all, a, b) = (eval(whole), eval(left), eval(right));
        for c in 0..9 {
            close(all[c], (0.3 * a[c] + 1.7 * b[c]) / 2.0, 3e-15);
        }
        assert_eq!(all, eval(whole));
    }
}

#[test]
fn projection_work_and_workspace_boundaries_fall_back_without_spending() {
    let window = HarmonicWindow::new(0.0, 1e-12).unwrap();
    let frequencies = [1.0, 2.0];
    for proper in [false, true] {
        let modes = if proper {
            window.band_mode_count_rad_s(&frequencies)
        } else {
            window.real_band_mode_count_rad_s(&frequencies)
        }
        .unwrap()
        .unwrap();
        let bytes = FIXED_BYTES + modes * core::mem::size_of::<FieldParts>();
        let extra = 3 * modes - 4;
        let mut work = extra;
        assert_eq!(
            admit(window, &frequencies, 1, proper, &mut work, bytes).unwrap(),
            Some(modes)
        );
        assert_eq!(work, 0);
        let mut work = extra;
        assert_eq!(
            admit(window, &frequencies, 1, proper, &mut work, bytes - 1).unwrap(),
            None
        );
        assert_eq!(work, extra);
        let mut work = extra - 1;
        assert_eq!(
            admit(window, &frequencies, 1, proper, &mut work, bytes).unwrap(),
            None
        );
        assert_eq!(work, extra - 1);
    }
    let mut headroom = usize::MAX;
    assert_eq!(
        admit(window, &[1.0], 1, false, &mut headroom, usize::MAX).unwrap(),
        None
    );
    assert_eq!(
        admit(
            HarmonicWindow::new(0.0, 2.0).unwrap(),
            &frequencies,
            1,
            false,
            &mut headroom,
            usize::MAX
        )
        .unwrap(),
        None
    );
    assert_eq!(
        admit(
            window,
            &frequencies,
            usize::MAX,
            true,
            &mut headroom,
            usize::MAX
        )
        .unwrap(),
        None
    );
}

#[test]
fn temporal_projection_polls_cancellation_before_and_during_work() {
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let spatial =
        Spatial::<RealKernel>::from_axes(grid(2, 1), area, |i| i as f64, |_| 0.0, &mut || false)
            .unwrap();
    let fields = vec![vec![[[Dd::ONE, Dd::ZERO]; 6]; 2]; 2];
    let modes = window
        .real_band_mode_count_rad_s(&[1.0, 2.0])
        .unwrap()
        .unwrap();
    for stop in [1, 4, 10] {
        let mut polls = 0;
        assert!(matches!(
            integrate(
                &spatial,
                window,
                &[1.0, 2.0],
                &fields,
                modes,
                &mut [Dd::ZERO; 9],
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(AngularSpectrumError::Cancelled)
        ));
    }
}

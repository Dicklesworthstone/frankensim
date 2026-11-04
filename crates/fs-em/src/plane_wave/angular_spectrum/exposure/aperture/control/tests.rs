use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

#[derive(Clone)]
struct Block {
    control: Option<usize>,
    frequency: usize,
    fields: Vec<[C64; 6]>,
}

fn action(
    controls: usize,
    omega: &[f64],
    kx: &[f64],
    blocks: &[Block],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> EnergyAction {
    let grid = AngularSpectrumGrid {
        width: kx.len(),
        height: 1,
        spacing_m: [1.0; 2],
    };
    let spatial = (!omega.is_empty()).then(|| {
        Spatial::<RealKernel>::from_axes(grid, aperture, |i| kx[i], |_| 0.0, &mut || false).unwrap()
    });
    let temporal = pair_kernels(omega.len(), window, |i| -omega[i], &mut || false).unwrap();
    let groups = blocks
        .iter()
        .map(|block| {
            let mut fields: Vec<_> = block
                .fields
                .iter()
                .map(|field| field.map(|v| [Dd::from_f64(v.re), Dd::from_f64(v.im)]))
                .collect();
            center_fields(
                &mut fields,
                omega[block.frequency],
                |i| kx[i],
                |_| 0.0,
                aperture,
                window,
                &mut || false,
            )
            .unwrap();
            SourceGroup {
                control: block.control,
                frequency: block.frequency,
                fields,
            }
        })
        .collect();
    EnergyAction::new(
        controls,
        kx.len(),
        omega.to_vec(),
        spatial,
        temporal,
        groups,
        2.0,
        8.0,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

fn scene() -> Vec<Block> {
    [
        (None, 0),
        (None, 1),
        (Some(0), 0),
        (Some(0), 1),
        (Some(1), 1),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (control, frequency))| Block {
        control,
        frequency,
        fields: (0..2)
            .map(|mode| {
                core::array::from_fn(|component| {
                    C64::new(
                        (index + component + 1) as f64 * if mode == 0 { 0.07 } else { -0.11 },
                        (1 + mode + frequency + component % 3) as f64
                            * if index % 2 == 0 { 0.03 } else { -0.05 },
                    )
                })
            })
            .collect(),
    })
    .collect()
}

/// Independent composite Simpson integral of actual real fields. This oracle
/// uses no production harmonic kernel, center rotation, gradient or action.
fn quadrature(
    omega: &[f64],
    kx: &[f64],
    blocks: &[Block],
    coefficients: &[C64],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
) -> f64 {
    let n = 192;
    let weight = |i| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let mut value = 0.0;
    for ix in 0..=n {
        let x = aperture.origin_m()[0] + aperture.extent_m()[0] * ix as f64 / n as f64;
        for it in 0..=n {
            let u = it as f64 / n as f64;
            let t = window.start_s() + window.duration_s() * u;
            let mut field = [0.0; 6];
            for block in blocks {
                let a = block.control.map_or(C64::ONE, |index| coefficients[index]);
                for (mode, values) in block.fields.iter().enumerate() {
                    let (s, c) = (omega[block.frequency] * t - kx[mode] * x).sin_cos();
                    for component in 0..6 {
                        let z = a * values[component];
                        field[component] += z.re * c - z.im * s;
                    }
                }
            }
            let energy = (0..3)
                .map(|i| field[i] * field[i] + 4.0 * field[i + 3] * field[i + 3])
                .sum::<f64>();
            let profile = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            value += weight(ix) * weight(it) * profile * energy / (9 * n * n) as f64;
        }
    }
    value
}

#[test]
fn finite_patch_energy_and_complex_source_derivatives_match_independent_integration() {
    let omega = [3.2, 4.7];
    let kx = [0.0, 1.7];
    let blocks = scene();
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.2], [0.63, 0.4]).unwrap();
    let coefficients = [C64::new(0.3, -0.7), C64::new(-0.6, 0.2)];
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.31, 0.42, profile).unwrap();
        let mut owner = action(2, &omega, &kx, &blocks, aperture, window);
        let actual = owner.evaluate(&coefficients, true, &mut || false).unwrap();
        let expected = quadrature(&omega, &kx, &blocks, &coefficients, aperture, window);
        close(actual, expected, 3e-9);
        let gradient = owner.gradient.clone();
        for source in 0..2 {
            for imaginary in [false, true] {
                let (mut plus, mut minus) = (coefficients, coefficients);
                let h = 1e-4;
                if imaginary {
                    plus[source].im += h;
                    minus[source].im -= h;
                } else {
                    plus[source].re += h;
                    minus[source].re -= h;
                }
                let derivative = (quadrature(&omega, &kx, &blocks, &plus, aperture, window)
                    - quadrature(&omega, &kx, &blocks, &minus, aperture, window))
                    / (2.0 * h);
                close(
                    if imaginary {
                        gradient[source].im
                    } else {
                        gradient[source].re
                    },
                    derivative,
                    4e-8,
                );
            }
        }
    }
}

#[test]
fn hessian_action_is_the_symmetric_positive_real_cartesian_derivative() {
    let blocks = scene();
    let aperture = AngularSpectrumAperture::try_new([-0.3, 0.0], [0.8, 1.0]).unwrap();
    let window = HarmonicWindow::new(0.17, 0.21).unwrap();
    let mut owner = action(2, &[3.2, 4.7], &[0.0, 1.7], &blocks, aperture, window);
    let x = [C64::new(0.8, -0.2), C64::new(-0.3, 0.7)];
    let v = [C64::new(0.2, 0.9), C64::new(-0.6, 0.3)];
    let w = [C64::new(-0.7, 0.1), C64::new(0.4, -0.2)];
    owner.evaluate(&v, false, &mut || false).unwrap();
    let hv = owner.gradient.clone();
    owner.evaluate(&w, false, &mut || false).unwrap();
    let hw = owner.gradient.clone();
    let dot = |a: &[C64], b: &[C64]| {
        a.iter()
            .zip(b)
            .map(|(x, y)| x.re * y.re + x.im * y.im)
            .sum::<f64>()
    };
    close(dot(&w, &hv), dot(&v, &hw), 3e-14);
    assert!(dot(&v, &hv) > 0.0);
    let h = 1e-4;
    let plus = core::array::from_fn::<_, 2, _>(|i| x[i] + v[i].scale(h));
    let minus = core::array::from_fn::<_, 2, _>(|i| x[i] - v[i].scale(h));
    owner.evaluate(&plus, true, &mut || false).unwrap();
    let gp = owner.gradient.clone();
    owner.evaluate(&minus, true, &mut || false).unwrap();
    for i in 0..2 {
        close(
            (gp[i].re - owner.gradient[i].re) / (2.0 * h),
            hv[i].re,
            3e-11,
        );
        close(
            (gp[i].im - owner.gradient[i].im) / (2.0 * h),
            hv[i].im,
            3e-11,
        );
    }
    let iv = v.map(|z| C64::new(-z.im, z.re));
    owner.evaluate(&iv, false, &mut || false).unwrap();
    assert!(
        (owner.gradient[0].re + hv[0].im).abs() > 1e-3,
        "finite-shutter Hessian must retain carrier-sum quadrature coupling"
    );
}

#[test]
fn fixed_illumination_cancels_before_energy_and_empty_controls_keep_their_order() {
    let field = [
        C64::from_re(3.0),
        C64::from_re(4.0),
        C64::ZERO,
        C64::from_re(-2.0),
        C64::from_re(1.5),
        C64::ZERO,
    ];
    let blocks = vec![
        Block {
            control: None,
            frequency: 0,
            fields: vec![field],
        },
        Block {
            control: Some(1),
            frequency: 0,
            fields: vec![field],
        },
    ];
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let mut owner = action(3, &[1.0], &[0.0], &blocks, aperture, window);
    close(
        owner
            .evaluate(&[C64::ZERO; 3], true, &mut || false)
            .unwrap(),
        25.0,
        2e-14,
    );
    close(owner.gradient[1].re, 50.0, 2e-14);
    let cancel = [C64::new(9.0, -3.0), C64::from_re(-1.0), C64::new(-2.0, 8.0)];
    assert_eq!(owner.evaluate(&cancel, true, &mut || false).unwrap(), 0.0);
    assert_eq!(owner.gradient, vec![C64::ZERO; 3]);
    // Cancellation of the fixed bank does not erase source curvature.
    owner.evaluate(&cancel, false, &mut || false).unwrap();
    close(owner.gradient[1].re, -50.0, 2e-14);
    assert_eq!(owner.gradient[0], C64::ZERO);
    assert_eq!(owner.gradient[2], C64::ZERO);
}

#[test]
fn short_shutter_quadratures_retain_the_small_positive_energy_at_a_field_node() {
    let mut field = [C64::ZERO; 6];
    field[0] = C64::ONE;
    field[4] = C64::from_re(0.5);
    let blocks = vec![Block {
        control: Some(0),
        frequency: 0,
        fields: vec![field],
    }];
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let width = 1e-8;
    let window = HarmonicWindow::new(-width * 0.5, width).unwrap();
    let mut owner = action(1, &[1.0], &[0.0], &blocks, aperture, window);
    let value = owner
        .evaluate(&[C64::new(0.0, 1.0)], true, &mut || false)
        .unwrap();
    let expected = width * width / 6.0;
    assert!(value > 0.0);
    assert!((value / expected - 1.0).abs() < 3e-14);
    assert!((owner.gradient[0].im / (2.0 * expected) - 1.0).abs() < 3e-14);
    close(owner.gradient[0].re, 0.0, 1e-25);
}

#[test]
fn cancellation_bad_coefficients_and_overflow_leave_the_action_reusable() {
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, 0.2).unwrap();
    let mut owner = action(2, &[3.2, 4.7], &[0.0, 1.7], &scene(), aperture, window);
    let a = [C64::new(0.2, -0.4), C64::new(0.3, 0.7)];
    let expected = owner.evaluate(&a, true, &mut || false).unwrap();
    let expected_gradient = owner.gradient.clone();
    let mut polls = 0;
    assert!(
        owner
            .evaluate(&a, true, &mut || {
                polls += 1;
                polls == 18
            })
            .is_err()
    );
    assert!(polls >= 18);
    assert!(owner.evaluate(&[], true, &mut || false).is_err());
    assert!(
        owner
            .evaluate(&[C64::new(f64::NAN, 0.0); 2], true, &mut || false)
            .is_err()
    );
    assert!(
        owner
            .evaluate(&[C64::new(f64::MAX, f64::MAX); 2], true, &mut || false)
            .is_err()
    );
    assert_eq!(owner.evaluate(&a, true, &mut || false).unwrap(), expected);
    assert_eq!(owner.gradient, expected_gradient);
    let mut dark = action(2, &[], &[0.0], &[], aperture, window);
    assert_eq!(dark.evaluate(&a, true, &mut || false).unwrap(), 0.0);
    assert_eq!(dark.gradient, vec![C64::ZERO; 2]);
}

#[test]
fn separable_xy_time_action_matches_the_original_full_pair_integrator() {
    let (nx, ny) = (2, 3);
    let grid = AngularSpectrumGrid {
        width: nx,
        height: ny,
        spacing_m: [1.0; 2],
    };
    let kx = [-0.9, 0.6];
    let ky = [-0.4, 0.2, 0.8];
    let omega = [2.7, 4.1];
    let aperture = AngularSpectrumAperture::try_new([0.31, -0.7], [0.6, 0.9]).unwrap();
    let window =
        HarmonicWindow::with_profile(-0.2, 0.8, HarmonicWindowProfile::Triangular).unwrap();
    let blocks: Vec<_> = [(None, 0), (None, 1), (Some(0), 1), (Some(1), 0)]
        .into_iter()
        .enumerate()
        .map(|(b, (control, frequency))| Block {
            control,
            frequency,
            fields: (0..nx * ny)
                .map(|mode| {
                    core::array::from_fn(|c| {
                        C64::new(
                            (1 + b + mode + c) as f64 * if mode % 2 == 0 { 0.03 } else { -0.05 },
                            (1 + 2 * b + mode + c % 3) as f64
                                * if c % 2 == 0 { -0.04 } else { 0.02 },
                        )
                    })
                })
                .collect(),
        })
        .collect();
    let spatial = || {
        Spatial::<RealKernel>::from_axes(grid, aperture, |i| kx[i], |i| ky[i], &mut || false)
            .unwrap()
    };
    let temporal = || pair_kernels(2, window, |i| -omega[i], &mut || false).unwrap();
    let groups = blocks
        .iter()
        .map(|block| {
            let mut fields: Vec<_> = block
                .fields
                .iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect();
            center_fields(
                &mut fields,
                omega[block.frequency],
                |i| kx[i % nx],
                |i| ky[i / nx],
                aperture,
                window,
                &mut || false,
            )
            .unwrap();
            SourceGroup {
                control: block.control,
                frequency: block.frequency,
                fields,
            }
        })
        .collect();
    let mut owner = EnergyAction::new(
        2,
        nx * ny,
        omega.to_vec(),
        Some(spatial()),
        temporal(),
        groups,
        2.0,
        8.0,
    )
    .unwrap();
    struct Direct(Dd);
    impl Accumulator for Direct {
        type Output = f64;
        fn add<K: Kernel>(&mut self, k: K, a: &FieldParts, b: &FieldParts) {
            for c in 0..3 {
                self.0 = self.0
                    + k.product(a[c], b[c])
                    + Dd::from_f64(4.0) * k.product(a[c + 3], b[c + 3]);
            }
        }
        fn finish(self, _: AngularSpectrumAperture) -> Result<f64, AngularSpectrumError> {
            Ok(self.0.to_f64())
        }
    }
    let direct = |coefficients: &[C64]| {
        let mut fields = vec![vec![[C64::ZERO; 6]; nx * ny]; 2];
        for block in &blocks {
            let a = block.control.map_or(C64::ONE, |i| coefficients[i]);
            for (mode, field) in block.fields.iter().enumerate() {
                for c in 0..6 {
                    fields[block.frequency][mode][c] =
                        fields[block.frequency][mode][c] + a * field[c];
                }
            }
        }
        let mut parts: Vec<Vec<FieldParts>> = fields
            .iter()
            .map(|line| {
                line.iter()
                    .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                    .collect()
            })
            .collect();
        for f in 0..2 {
            center_fields(
                &mut parts[f],
                omega[f],
                |i| kx[i % nx],
                |i| ky[i / nx],
                aperture,
                window,
                &mut || false,
            )
            .unwrap();
        }
        let mut sum = Direct(Dd::ZERO);
        integrate(&spatial(), &temporal(), &parts, &mut sum, &mut || false).unwrap();
        sum.finish(aperture).unwrap()
    };
    let a = [C64::new(-0.2, 0.7), C64::new(0.3, -0.6)];
    let value = owner.evaluate(&a, true, &mut || false).unwrap();
    close(value, direct(&a), 5e-14);
    let direction = [C64::new(0.4, -0.1), C64::new(-0.5, 0.8)];
    let derivative = owner
        .gradient
        .iter()
        .zip(direction)
        .map(|(g, v)| g.re * v.re + g.im * v.im)
        .sum::<f64>();
    let h = 1e-4;
    let plus = core::array::from_fn::<_, 2, _>(|i| a[i] + direction[i].scale(h));
    let minus = core::array::from_fn::<_, 2, _>(|i| a[i] - direction[i].scale(h));
    close(
        derivative,
        (direct(&plus) - direct(&minus)) / (2.0 * h),
        3e-11,
    );
}

use super::*;

fn config(averaging: DynamicAveraging) -> DynamicSgsConfig {
    DynamicSgsConfig {
        averaging,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn field(n: usize) -> (Vec<f64>, Vec<[f64; 3]>) {
    cartesian_field([n; 3])
}

fn cartesian_field(shape: [usize; 3]) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut rho = Vec::new();
    let mut u = Vec::new();
    for z in 0..shape[2] {
        for y in 0..shape[1] {
            for x in 0..shape[0] {
                let [x, y, z] = core::array::from_fn::<_, 3, _>(|axis| {
                    core::f64::consts::TAU * ([x, y, z][axis] as f64 + 0.5) / shape[axis] as f64
                });
                rho.push(1.2 + 0.2 * x.sin() + 0.1 * y.cos());
                u.push([
                    0.3 * x.sin() + 0.2 * (y + z).sin() + 0.1 * (2.0 * x - y).cos(),
                    0.2 * y.cos() + 0.3 * (z + x).sin() + 0.1 * (x - 2.0 * z).sin(),
                    0.2 * z.sin() + 0.1 * (x + y).cos() + 0.3 * (2.0 * y - z).sin(),
                ]);
            }
        }
    }
    (rho, u)
}

// Independent full 3x3 physical-unit reference. It uses raw first/second
// moments, explicit modulo indexing and no packed tensors or normalization.
fn reference(
    shape: [usize; 3],
    spacing: [f64; 3],
    rho: &[f64],
    u: &[[f64; 3]],
    averaging: DynamicAveraging,
) -> (Vec<f64>, Vec<f64>) {
    let delta_squared = spacing.iter().product::<f64>().cbrt().powi(2);
    let at = |x: usize, y: usize, z: usize| {
        ((z % shape[2]) * shape[1] + y % shape[1]) * shape[0] + x % shape[0]
    };
    let xyz = |c: usize| {
        [
            c % shape[0],
            c / shape[0] % shape[1],
            c / (shape[0] * shape[1]),
        ]
    };
    let neighbors = |c: usize| -> Vec<usize> {
        let [x, y, z] = xyz(c);
        (0..3)
            .flat_map(move |k| {
                (0..3).flat_map(move |j| {
                    (0..3).map(move |i| {
                        at(
                            x + shape[0] + i - 1,
                            y + shape[1] + j - 1,
                            z + shape[2] + k - 1,
                        )
                    })
                })
            })
            .collect()
    };
    let strain = |v: &[[f64; 3]], c: usize| -> (f64, [[f64; 3]; 3]) {
        let p = xyz(c);
        let mut g = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut a = p;
            let mut b = p;
            a[j] += shape[j] - 1;
            b[j] += 1;
            for i in 0..3 {
                g[i][j] =
                    (v[at(b[0], b[1], b[2])][i] - v[at(a[0], a[1], a[2])][i]) / (2.0 * spacing[j]);
            }
        }
        let mut s = core::array::from_fn::<_, 3, _>(|i| {
            core::array::from_fn::<_, 3, _>(|j| (g[i][j] + g[j][i]) / 2.0)
        });
        let norm = (2.0 * s.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt();
        let trace = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
        for i in 0..3 {
            s[i][i] -= trace;
        }
        (norm, s)
    };
    let mut filtered = vec![[0.0; 3]; u.len()];
    let mut rh = vec![0.0; u.len()];
    let mut base = vec![[[0.0; 3]; 3]; u.len()];
    let mut rate = vec![0.0; u.len()];
    for c in 0..u.len() {
        for d in neighbors(c) {
            rh[c] += rho[d] / 27.0;
            for i in 0..3 {
                filtered[c][i] += rho[d] * u[d][i] / 27.0;
            }
        }
        for i in 0..3 {
            filtered[c][i] /= rh[c];
        }
        let (r, s) = strain(u, c);
        rate[c] = r;
        for i in 0..3 {
            for j in 0..3 {
                base[c][i][j] = rho[c] * delta_squared * r * s[i][j];
            }
        }
    }
    let mut pairs = vec![[0.0; 2]; u.len()];
    for c in 0..u.len() {
        let (r, s) = strain(&filtered, c);
        let mut l = [[0.0; 3]; 3];
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut product = 0.0;
                let mut a = 0.0;
                for d in neighbors(c) {
                    product += rho[d] * u[d][i] * u[d][j] / 27.0;
                    a += base[d][i][j] / 27.0;
                }
                l[i][j] = product - rh[c] * filtered[c][i] * filtered[c][j];
                m[i][j] = 2.0 * (a - 9.0 * delta_squared * rh[c] * r * s[i][j]);
            }
        }
        let trace = (l[0][0] + l[1][1] + l[2][2]) / 3.0;
        for i in 0..3 {
            l[i][i] -= trace;
        }
        for i in 0..3 {
            for j in 0..3 {
                pairs[c][0] += l[i][j] * m[i][j];
                pairs[c][1] += m[i][j] * m[i][j];
            }
        }
    }
    let domain = [
        pairs.iter().map(|p| p[0]).sum::<f64>(),
        pairs.iter().map(|p| p[1]).sum::<f64>(),
    ];
    let coefficients = (0..u.len())
        .map(|c| {
            let p = match averaging {
                DynamicAveraging::Domain => domain,
                DynamicAveraging::LocalBox => neighbors(c).iter().fold([0.0; 2], |mut a, &d| {
                    a[0] += pairs[d][0] / 27.0;
                    a[1] += pairs[d][1] / 27.0;
                    a
                }),
            };
            if p[1] == 0.0 { 0.0 } else { p[0] / p[1] }
        })
        .collect::<Vec<_>>();
    let nu = coefficients
        .iter()
        .zip(rate)
        .map(|(c, r)| c.max(0.0) * delta_squared * r)
        .collect();
    (coefficients, nu)
}

#[test]
fn dynamic_favre_filter_matches_independent_full_tensor_reference() {
    let n = 7;
    let (rho, u) = field(n);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        for density in [&rho, &vec![1.0; u.len()]] {
            let mut work = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
            let report = work
                .prepare(config(averaging), density, &u, || false)
                .unwrap();
            let (coefficient, expected) = reference([n; 3], [1.0; 3], density, &u, averaging);
            for (a, b) in work.viscosities().iter().zip(expected) {
                assert!((a - b).abs() < 2e-14, "{averaging:?}: {a} != {b}");
            }
            assert!(
                (report.minimum_raw_coefficient_squared
                    - coefficient.iter().copied().fold(f64::INFINITY, f64::min))
                .abs()
                    < 2e-14
            );
            if averaging == DynamicAveraging::LocalBox {
                assert!(report.maximum_eddy_viscosity > 1e-5);
                assert!(report.backscatter_limited_cells > 0);
            }
        }
    }
}

#[test]
fn cartesian_favre_momentum_matches_independent_physical_tensor_reference() {
    let shape = [5, 6, 7];
    let spacing = [0.3, 0.7, 1.1];
    let (rho, u) = cartesian_field(shape);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        for density in [&rho, &vec![1.0; u.len()]] {
            let mut work = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
            let report = work
                .prepare(config(averaging), density, &u, || false)
                .unwrap();
            let (coefficient, expected) = reference(shape, spacing, density, &u, averaging);
            for (a, b) in work.viscosities().iter().zip(expected) {
                assert!((a - b).abs() < 2e-13, "{averaging:?}: {a} != {b}");
            }
            let minimum = coefficient.iter().copied().fold(f64::INFINITY, f64::min);
            assert!((report.minimum_raw_coefficient_squared - minimum).abs() < 2e-13);
            if averaging == DynamicAveraging::LocalBox {
                assert!(report.maximum_eddy_viscosity > 1e-5);
                assert!(report.backscatter_limited_cells > 0);
            }
        }
    }
}

#[test]
fn cartesian_geometry_preserves_bytes_and_refuses_unrepresentable_metrics() {
    let shape = [5, 6, 7];
    let bytes = PeriodicDynamicSgs::required_bytes(shape.iter().product()).unwrap();
    assert_eq!(bytes, 136 * shape.iter().product::<usize>());
    for scale in [1.0, 1e-250, 1e250] {
        let spacing = [scale, 2.0 * scale, 4.0 * scale];
        let work = PeriodicDynamicSgs::new_cartesian(shape, spacing, bytes).unwrap();
        assert_eq!(work.shape(), shape);
        assert_eq!(work.spacing(), spacing);
        assert!((work.filter_width() / scale - 2.0).abs() < 3e-15);
        assert_eq!(work.workspace_bytes(), bytes);
    }
    let tiny = f64::from_bits(1);
    let work =
        PeriodicDynamicSgs::new_cartesian(shape, [tiny, 2.0 * tiny, 4.0 * tiny], bytes).unwrap();
    assert_eq!(work.filter_width(), 2.0 * tiny);
    let work =
        PeriodicDynamicSgs::new_cartesian(shape, [f64::MAX / 4.0, f64::MAX / 2.0, f64::MAX], bytes)
            .unwrap();
    assert!((work.filter_width() / (f64::MAX / 2.0) - 1.0).abs() < 3e-15);
    for ulps in [1, 2, 4, 16, 1024] {
        let spacing = [
            f64::MAX,
            f64::from_bits(f64::MAX.to_bits() - ulps),
            f64::from_bits(f64::MAX.to_bits() - 2 * ulps),
        ];
        let work = PeriodicDynamicSgs::new_cartesian(shape, spacing, bytes).unwrap();
        assert!(work.filter_width() >= spacing[2] && work.filter_width() <= spacing[0]);
    }
    for offsets in [[0, 0, 1], [0, 1, 1]] {
        let spacing = offsets.map(|ulps| f64::from_bits(f64::MAX.to_bits() - ulps));
        let work = PeriodicDynamicSgs::new_cartesian(shape, spacing, bytes).unwrap();
        assert!(work.filter_width() >= spacing[2] && work.filter_width() <= spacing[0]);
    }
    for spacing in [
        [0.0, 1.0, 2.0],
        [-1.0, 1.0, 2.0],
        [f64::NAN, 1.0, 2.0],
        [1.0, f64::INFINITY, 2.0],
        [tiny, 1.0, f64::MAX],
    ] {
        assert!(PeriodicDynamicSgs::new_cartesian(shape, spacing, bytes).is_err());
    }
    assert!(matches!(
        PeriodicDynamicSgs::new_cartesian(shape, [0.3, 0.7, 1.1], bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    assert!(PeriodicDynamicSgs::new_cartesian([2, 6, 7], [0.3, 0.7, 1.1], bytes).is_err());
    let mut extreme =
        PeriodicDynamicSgs::new_cartesian(shape, [1e-150, 1.0, 1e150], bytes).unwrap();
    let (rho, u) = cartesian_field(shape);
    assert!(matches!(
        extreme.prepare(config(DynamicAveraging::LocalBox), &rho, &u, || false),
        Err(DynamicSgsError::Cell { .. })
    ));
}

#[test]
fn dynamic_filter_has_the_declared_box_fourier_response() {
    let n = 9;
    let len = n * n * n;
    let mut u = vec![[0.0; 3]; len];
    let k = core::f64::consts::TAU / n as f64;
    for c in 0..len {
        u[c][0] = (k * (c % n) as f64).sin();
    }
    let mut work = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    work.prepare(
        config(DynamicAveraging::LocalBox),
        &vec![1.0; len],
        &u,
        || false,
    )
    .unwrap();
    let scale = u.iter().map(|v| v[0].abs()).fold(0.0, f64::max);
    let response = (1.0 + 2.0 * k.cos()) / 3.0;
    for c in 0..len {
        assert!((work.filtered[c][0] * scale - response * u[c][0]).abs() < 1e-15);
    }
}

#[test]
fn dynamic_translation_and_laminar_periodic_shear_vanish() {
    let n = 7;
    let (rho, mut u) = field(n);
    u.fill([3.0, -4.0, 5.0]);
    let mut work = PeriodicDynamicSgs::new([n; 3], 0.25, usize::MAX).unwrap();
    let allocation = work.viscosities().as_ptr();
    let report = work
        .prepare(config(DynamicAveraging::LocalBox), &rho, &u, || false)
        .unwrap();
    assert_eq!(report.undetermined_cells, u.len());
    assert!(work.viscosities().iter().all(|&v| v == 0.0));
    for (c, v) in u.iter_mut().enumerate() {
        *v = [
            (core::f64::consts::TAU * (c / n % n) as f64 / n as f64).sin(),
            0.0,
            0.0,
        ];
    }
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        work.prepare(config(averaging), &vec![1.0; u.len()], &u, || false)
            .unwrap();
        assert!(work.viscosities().iter().all(|&v| v == 0.0));
        assert_eq!(allocation, work.viscosities().as_ptr());
    }
}

#[test]
fn dynamic_coefficients_are_boost_invariant_and_scale_safe() {
    let n = 7;
    let (rho, u) = field(n);
    let cfg = config(DynamicAveraging::LocalBox);
    let mut work = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    work.prepare(cfg, &rho, &u, || false).unwrap();
    let expected = work.viscosities().to_vec();
    let boosted = u
        .iter()
        .map(|v| [v[0] + 1e5, v[1] - 3e5, v[2] + 2e5])
        .collect::<Vec<_>>();
    work.prepare(cfg, &rho, &boosted, || false).unwrap();
    for (a, b) in work.viscosities().iter().zip(&expected) {
        assert!((a - b).abs() < 2e-11);
    }
    for (u_scale, rho_scale, h) in [(1e150, 1e-150, 0.25), (1e-150, 1e150, 2.0)] {
        let mut scaled = PeriodicDynamicSgs::new([n; 3], h, usize::MAX).unwrap();
        let rho = rho.iter().map(|v| v * rho_scale).collect::<Vec<_>>();
        let u = u.iter().map(|v| v.map(|x| x * u_scale)).collect::<Vec<_>>();
        scaled.prepare(cfg, &rho, &u, || false).unwrap();
        for (a, b) in scaled.viscosities().iter().zip(&expected) {
            assert!((a / (u_scale * h) - b).abs() < 2e-14);
        }
    }
}

#[test]
fn dynamic_negative_coefficient_is_reported_and_positive_cap_refuses() {
    let n = 7;
    let (rho, u) = field(n);
    let mut work = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let cfg = config(DynamicAveraging::LocalBox);
    let report = work.prepare(cfg, &rho, &u, || false).unwrap();
    assert!(report.backscatter_limited_cells > 0 && report.maximum_coefficient_squared > 0.0);
    assert!(work.viscosities().iter().all(|&v| v >= 0.0));
    assert!(matches!(
        work.prepare(
            DynamicSgsConfig {
                maximum_coefficient_squared: 0.0,
                ..cfg
            },
            &rho,
            &u,
            || false
        ),
        Err(DynamicSgsError::Coefficient { .. })
    ));
}

#[test]
fn dynamic_admission_and_cancellation_are_bounded() {
    let n = 5;
    let (mut rho, u) = field(n);
    let bytes = PeriodicDynamicSgs::required_bytes(u.len()).unwrap();
    assert!(matches!(
        PeriodicDynamicSgs::new([n; 3], 1.0, bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    assert!(PeriodicDynamicSgs::new([2, n, n], 1.0, usize::MAX).is_err());
    let mut work = PeriodicDynamicSgs::new([n; 3], 1.0, bytes).unwrap();
    let cfg = config(DynamicAveraging::Domain);
    assert_eq!(
        work.prepare(cfg, &rho, &u, || true),
        Err(DynamicSgsError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        work.prepare(cfg, &rho, &u, || {
            polls += 1;
            polls == 5
        }),
        Err(DynamicSgsError::Cancelled)
    );
    assert_eq!(polls, 5);
    rho[17] = f64::NAN;
    assert!(matches!(
        work.prepare(cfg, &rho, &u, || false),
        Err(DynamicSgsError::Cell { cell: 17, .. })
    ));
    assert!(
        DynamicSgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..cfg
        }
        .validate()
        .is_err()
    );
}

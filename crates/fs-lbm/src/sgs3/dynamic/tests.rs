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
    let mut rho = Vec::new();
    let mut u = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let [x, y, z] =
                    [x, y, z].map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / n as f64);
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
    n: usize,
    rho: &[f64],
    u: &[[f64; 3]],
    averaging: DynamicAveraging,
) -> (Vec<f64>, Vec<f64>) {
    let at = |x: usize, y: usize, z: usize| ((z % n) * n + y % n) * n + x % n;
    let neighbors = |c: usize| -> Vec<usize> {
        let [x, y, z] = [c % n, c / n % n, c / (n * n)];
        (0..3)
            .flat_map(move |k| {
                (0..3).flat_map(move |j| {
                    (0..3).map(move |i| at(x + n + i - 1, y + n + j - 1, z + n + k - 1))
                })
            })
            .collect()
    };
    let strain = |v: &[[f64; 3]], c: usize| -> (f64, [[f64; 3]; 3]) {
        let p = [c % n, c / n % n, c / (n * n)];
        let mut g = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut a = p;
            let mut b = p;
            a[j] += n - 1;
            b[j] += 1;
            for i in 0..3 {
                g[i][j] = (v[at(b[0], b[1], b[2])][i] - v[at(a[0], a[1], a[2])][i]) / 2.0;
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
                base[c][i][j] = rho[c] * r * s[i][j];
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
                m[i][j] = 2.0 * (a - 9.0 * rh[c] * r * s[i][j]);
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
        .map(|(c, r)| c.max(0.0) * r)
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
            let (coefficient, expected) = reference(n, density, &u, averaging);
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

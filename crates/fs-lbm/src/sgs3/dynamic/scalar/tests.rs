use super::*;

fn settings(averaging: DynamicAveraging) -> (DynamicSgsConfig, DynamicScalarConfig) {
    (
        DynamicSgsConfig {
            averaging,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        DynamicScalarConfig {
            averaging,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    )
}

fn field(n: usize) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut phi = Vec::new();
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
                phi.push(
                    0.7 * (x + y).sin() + 0.9 * (2.0 * y - z).cos() + 0.4 * (x - 2.0 * z).sin(),
                );
            }
        }
    }
    (rho, u, phi)
}

// Independent full physical-unit reference: raw first/second moments, explicit
// modulo indexing, all nine strain entries, no production filter helpers or
// normalization. The h^2 and test-width factor are present in the model tensor.
fn reference(
    n: usize,
    h: f64,
    rho: &[f64],
    u: &[[f64; 3]],
    phi: &[f64],
    averaging: DynamicAveraging,
) -> (Vec<f64>, Vec<f64>) {
    let len = u.len();
    let at = |p: [usize; 3]| ((p[2] % n) * n + p[1] % n) * n + p[0] % n;
    let xyz = |c: usize| [c % n, c / n % n, c / (n * n)];
    let box_cells = |c: usize| {
        let p = xyz(c);
        let mut result = Vec::new();
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    result.push(at([p[0] + n + x - 1, p[1] + n + y - 1, p[2] + n + z - 1]));
                }
            }
        }
        result
    };
    let derivative = |v: &[f64], c: usize, axis: usize| {
        let mut a = xyz(c);
        let mut b = a;
        a[axis] += n - 1;
        b[axis] += 1;
        (v[at(b)] - v[at(a)]) / (2.0 * h)
    };
    let rates = |velocity: &[[f64; 3]]| {
        let components: Vec<Vec<f64>> = (0..3)
            .map(|i| velocity.iter().map(|v| v[i]).collect())
            .collect();
        (0..len)
            .map(|c| {
                let mut norm = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let s = (derivative(&components[i], c, j)
                            + derivative(&components[j], c, i))
                            / 2.0;
                        norm += 2.0 * s * s;
                    }
                }
                norm.sqrt()
            })
            .collect::<Vec<_>>()
    };
    let rate = rates(u);
    let mut mean_rho = vec![0.0; len];
    let mut mean_u = vec![[0.0; 3]; len];
    let mut mean_phi = vec![0.0; len];
    for c in 0..len {
        for j in box_cells(c) {
            mean_rho[c] += rho[j] / 27.0;
            mean_phi[c] += rho[j] * phi[j] / 27.0;
            for i in 0..3 {
                mean_u[c][i] += rho[j] * u[j][i] / 27.0;
            }
        }
        mean_phi[c] /= mean_rho[c];
        for i in 0..3 {
            mean_u[c][i] /= mean_rho[c];
        }
    }
    let mean_rate = rates(&mean_u);
    let mut pairs = vec![[0.0; 2]; len];
    for c in 0..len {
        for i in 0..3 {
            let mut moment = 0.0;
            let mut a = 0.0;
            for j in box_cells(c) {
                moment += rho[j] * u[j][i] * phi[j] / 27.0;
                a += rho[j] * h * h * rate[j] * derivative(phi, j, i) / 27.0;
            }
            let l = moment - mean_rho[c] * mean_u[c][i] * mean_phi[c];
            let m = a - 9.0 * h * h * mean_rho[c] * mean_rate[c] * derivative(&mean_phi, c, i);
            pairs[c][0] += l * m;
            pairs[c][1] += m * m;
        }
    }
    let domain = [
        pairs.iter().map(|p| p[0]).sum(),
        pairs.iter().map(|p| p[1]).sum(),
    ];
    let coefficient = (0..len)
        .map(|c| {
            let p = match averaging {
                DynamicAveraging::Domain => domain,
                DynamicAveraging::LocalBox => {
                    box_cells(c).into_iter().fold([0.0; 2], |mut a, j| {
                        a[0] += pairs[j][0] / 27.0;
                        a[1] += pairs[j][1] / 27.0;
                        a
                    })
                }
            };
            if p[1] == 0.0 { 0.0 } else { p[0] / p[1] }
        })
        .collect::<Vec<_>>();
    let alpha = coefficient
        .iter()
        .zip(rate)
        .map(|(c, s)| c.max(0.0) * h * h * s)
        .collect();
    (coefficient, alpha)
}

#[test]
fn g1_scalar_germano_matches_independent_physical_unit_reference() {
    let n = 7;
    let h = 0.37;
    let (rho, u, phi) = field(n);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let (mc, sc) = settings(averaging);
        for density in [&rho, &vec![1.0; u.len()]] {
            let mut flow = PeriodicDynamicSgs::new([n; 3], h, usize::MAX).unwrap();
            let mut scalar = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
            let (_, report) = flow
                .prepare_with_scalar(mc, sc, density, &u, &phi, &mut scalar, || false)
                .unwrap();
            let (coefficient, expected) = reference(n, h, density, &u, &phi, averaging);
            for (a, b) in scalar.diffusivities().iter().zip(expected) {
                assert!((a - b).abs() < 2e-13, "{averaging:?}: {a} != {b}");
            }
            let minimum = coefficient.into_iter().fold(f64::INFINITY, f64::min);
            assert!((report.minimum_raw_coefficient - minimum).abs() < 2e-13);
            if averaging == DynamicAveraging::LocalBox {
                assert!(report.maximum_diffusivity > 1e-4);
                assert!(report.backscatter_limited_cells > 0);
            }
        }
    }
}

#[test]
fn g3_scalar_fit_preserves_momentum_bits_and_affine_unit_invariance() {
    let n = 7;
    let (rho, u, phi) = field(n);
    let (mc, sc) = settings(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
    let ordinary = flow.prepare(mc, &rho, &u, || false).unwrap();
    let momentum = flow.viscosities().to_vec();
    let (paired, _) = flow
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
        .unwrap();
    assert_eq!(paired, ordinary);
    assert_eq!(flow.viscosities(), momentum);
    let expected = scalar.diffusivities().to_vec();
    let allocation = scalar.diffusivities().as_ptr();
    let boosted = u
        .iter()
        .map(|v| [v[0] + 1e4, v[1] - 3e4, v[2] + 2e4])
        .collect::<Vec<_>>();
    let changed = phi.iter().map(|v| -2.5 * v + 273.15).collect::<Vec<_>>();
    flow.prepare_with_scalar(mc, sc, &rho, &boosted, &changed, &mut scalar, || false)
        .unwrap();
    for (a, b) in scalar.diffusivities().iter().zip(&expected) {
        assert!((a - b).abs() < 3e-11);
    }
    assert_eq!(scalar.diffusivities().as_ptr(), allocation);
    for (us, rs, ps, h) in [(1e150, 1e-150, 1e-150, 0.25), (1e-150, 1e150, 1e150, 2.0)] {
        let mut scaled = PeriodicDynamicSgs::new([n; 3], h, usize::MAX).unwrap();
        let rho = rho.iter().map(|v| v * rs).collect::<Vec<_>>();
        let u = u.iter().map(|v| v.map(|q| q * us)).collect::<Vec<_>>();
        let phi = phi.iter().map(|v| v * ps).collect::<Vec<_>>();
        scaled
            .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
            .unwrap();
        for (a, b) in scalar.diffusivities().iter().zip(&expected) {
            assert!((a / (us * h) - b).abs() < 2e-13);
        }
    }
}

#[test]
fn g0_scalar_diffusion_can_be_active_with_exactly_zero_momentum_fit() {
    let n = 9;
    let len = n * n * n;
    let rho = vec![1.0; len];
    let mut u = vec![[0.0; 3]; len];
    let mut phi = vec![0.0; len];
    for c in 0..len {
        let x = core::f64::consts::TAU * (c % n) as f64 / n as f64;
        let y = core::f64::consts::TAU * (c / n % n) as f64 / n as f64;
        u[c][0] = 0.3 * y.sin();
        phi[c] = (x + y).sin();
    }
    let (mc, sc) = settings(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
    let (momentum, heat) = flow
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
        .unwrap();
    assert_eq!(momentum.maximum_eddy_viscosity, 0.0);
    assert!(flow.viscosities().iter().all(|&v| v == 0.0));
    assert!(heat.maximum_diffusivity > 1e-4);
    assert!(heat.backscatter_limited_cells > 0);
}

#[test]
fn g0_uniform_fields_clear_prior_diffusivity_and_still_validate_inputs() {
    let n = 5;
    let (rho, u, mut phi) = field(n);
    let len = u.len();
    let (mc, sc) = settings(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
    for constant_velocity in [true, false] {
        flow.prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
            .unwrap();
        assert!(scalar.diffusivities().iter().any(|&v| v > 0.0));
        let v = if constant_velocity {
            vec![[1.0, -2.0, 3.0]; len]
        } else {
            u.clone()
        };
        let s = if constant_velocity {
            phi.clone()
        } else {
            vec![300.0; len]
        };
        let (_, report) = flow
            .prepare_with_scalar(mc, sc, &rho, &v, &s, &mut scalar, || false)
            .unwrap();
        assert_eq!(report.undetermined_cells, len);
        assert!(scalar.diffusivities().iter().all(|&v| v == 0.0));
    }
    phi[11] = f64::NAN;
    assert!(matches!(
        flow.prepare_with_scalar(
            mc,
            sc,
            &rho,
            &vec![[0.0; 3]; len],
            &phi,
            &mut scalar,
            || false
        ),
        Err(DynamicSgsError::Cell { cell: 11, .. })
    ));
    assert!(
        flow.prepare_with_scalar(mc, sc, &rho, &u, &phi[..len - 1], &mut scalar, || false)
            .is_err()
    );
}

#[test]
fn g0_g4_scalar_caps_budgets_and_cancellation_refuse_and_retry_exactly() {
    let n = 7;
    let (rho, u, phi) = field(n);
    let len = u.len();
    let bytes = PeriodicDynamicScalar::required_bytes(len).unwrap();
    assert_eq!(bytes, 64 * len);
    assert!(matches!(
        PeriodicDynamicScalar::new([n; 3], bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    assert!(PeriodicDynamicScalar::new([2, n, n], usize::MAX).is_err());
    let (mc, sc) = settings(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new([n; 3], 1.0, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new([n; 3], bytes).unwrap();
    assert_eq!(scalar.workspace_bytes(), bytes);
    let mut count = 0;
    let reference = flow
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || {
            count += 1;
            false
        })
        .unwrap();
    let expected = scalar.diffusivities().to_vec();
    assert!(reference.1.minimum_raw_coefficient < 0.0);
    for cutoff in [1, count / 2, count - 2, count] {
        let mut polls = 0;
        assert_eq!(
            flow.prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || {
                polls += 1;
                polls == cutoff
            }),
            Err(DynamicSgsError::Cancelled)
        );
        assert_eq!(polls, cutoff);
        assert_eq!(
            flow.prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
                .unwrap(),
            reference
        );
        assert_eq!(scalar.diffusivities(), expected);
    }
    assert!(matches!(
        flow.prepare_with_scalar(
            mc,
            DynamicScalarConfig {
                maximum_coefficient: 0.0,
                ..sc
            },
            &rho,
            &u,
            &phi,
            &mut scalar,
            || false
        ),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    for bad in [
        DynamicScalarConfig {
            usage: SgsUse::AcceptanceRender,
            ..sc
        },
        DynamicScalarConfig {
            maximum_coefficient: f64::NAN,
            ..sc
        },
        DynamicScalarConfig {
            maximum_diffusivity_to_molecular: -1.0,
            ..sc
        },
    ] {
        assert!(bad.validate().is_err());
    }
    let mut wrong = PeriodicDynamicScalar::new([n + 1, n, n], usize::MAX).unwrap();
    assert!(
        flow.prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut wrong, || false)
            .is_err()
    );
}

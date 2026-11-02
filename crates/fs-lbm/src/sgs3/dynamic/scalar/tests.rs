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
    cartesian_field([n; 3])
}

fn cartesian_field(shape: [usize; 3]) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut phi = Vec::new();
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
// normalization. Delta^2 and the test-width factor are present in the model.
fn reference(
    shape: [usize; 3],
    spacing: [f64; 3],
    rho: &[f64],
    u: &[[f64; 3]],
    phi: &[f64],
    averaging: DynamicAveraging,
) -> (Vec<f64>, Vec<f64>) {
    let len = u.len();
    let delta_squared = spacing.iter().product::<f64>().cbrt().powi(2);
    let at = |p: [usize; 3]| {
        ((p[2] % shape[2]) * shape[1] + p[1] % shape[1]) * shape[0] + p[0] % shape[0]
    };
    let xyz = |c: usize| {
        [
            c % shape[0],
            c / shape[0] % shape[1],
            c / (shape[0] * shape[1]),
        ]
    };
    let box_cells = |c: usize| {
        let p = xyz(c);
        let mut result = Vec::new();
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    result.push(at([
                        p[0] + shape[0] + x - 1,
                        p[1] + shape[1] + y - 1,
                        p[2] + shape[2] + z - 1,
                    ]));
                }
            }
        }
        result
    };
    let derivative = |v: &[f64], c: usize, axis: usize| {
        let mut a = xyz(c);
        let mut b = a;
        a[axis] += shape[axis] - 1;
        b[axis] += 1;
        (v[at(b)] - v[at(a)]) / (2.0 * spacing[axis])
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
                a += rho[j] * delta_squared * rate[j] * derivative(phi, j, i) / 27.0;
            }
            let l = moment - mean_rho[c] * mean_u[c][i] * mean_phi[c];
            let m =
                a - 9.0 * delta_squared * mean_rho[c] * mean_rate[c] * derivative(&mean_phi, c, i);
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
        .map(|(c, s)| c.max(0.0) * delta_squared * s)
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
            let (coefficient, expected) = reference([n; 3], [h; 3], density, &u, &phi, averaging);
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
fn cartesian_scalar_matches_independent_physical_unit_reference() {
    let shape = [5, 6, 7];
    let spacing = [0.3, 0.7, 1.1];
    let (rho, u, phi) = cartesian_field(shape);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let (mc, sc) = settings(averaging);
        for density in [&rho, &vec![1.0; u.len()]] {
            let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
            let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
            let (_, report) = flow
                .prepare_with_scalar(mc, sc, density, &u, &phi, &mut scalar, || false)
                .unwrap();
            let (coefficient, expected) = reference(shape, spacing, density, &u, &phi, averaging);
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
fn cartesian_equal_spacing_preserves_both_fits_bit_for_bit() {
    let shape = [5, 6, 7];
    let (rho, u, phi) = cartesian_field(shape);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let (mc, sc) = settings(averaging);
        for h in [0.37, 1.0, 1e-150, 1e150] {
            let mut old = PeriodicDynamicSgs::new(shape, h, usize::MAX).unwrap();
            let mut new = PeriodicDynamicSgs::new_cartesian(shape, [h; 3], usize::MAX).unwrap();
            let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
            let expected = old
                .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
                .unwrap();
            let alpha = scalar
                .diffusivities()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>();
            let observed = new
                .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
                .unwrap();
            assert_eq!(new.filter_width().to_bits(), h.to_bits());
            assert_eq!(observed, expected);
            assert!(
                new.viscosities()
                    .iter()
                    .zip(old.viscosities())
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
            assert!(
                scalar
                    .diffusivities()
                    .iter()
                    .zip(alpha)
                    .all(|(a, b)| a.to_bits() == b)
            );
        }
    }
}

#[test]
fn cartesian_fits_covary_under_axis_permutation() {
    let shape = [5, 6, 7];
    let spacing = [0.2, 0.35, 0.8];
    let (rho, u, phi) = cartesian_field(shape);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let (mc, sc) = settings(averaging);
        let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
        let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
        let expected_report = flow
            .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
            .unwrap();
        for permutation in [[1, 0, 2], [2, 0, 1]] {
            let new_shape = permutation.map(|axis| shape[axis]);
            let new_spacing = permutation.map(|axis| spacing[axis]);
            let mut source = Vec::new();
            for z in 0..new_shape[2] {
                for y in 0..new_shape[1] {
                    for x in 0..new_shape[0] {
                        let mut p = [0; 3];
                        for axis in 0..3 {
                            p[permutation[axis]] = [x, y, z][axis];
                        }
                        source.push((p[2] * shape[1] + p[1]) * shape[0] + p[0]);
                    }
                }
            }
            let rotated_rho = source.iter().map(|&c| rho[c]).collect::<Vec<_>>();
            let rotated_u = source
                .iter()
                .map(|&c| permutation.map(|axis| u[c][axis]))
                .collect::<Vec<_>>();
            let rotated_phi = source.iter().map(|&c| phi[c]).collect::<Vec<_>>();
            let mut rotated =
                PeriodicDynamicSgs::new_cartesian(new_shape, new_spacing, usize::MAX).unwrap();
            let mut rotated_scalar = PeriodicDynamicScalar::new(new_shape, usize::MAX).unwrap();
            let report = rotated
                .prepare_with_scalar(
                    mc,
                    sc,
                    &rotated_rho,
                    &rotated_u,
                    &rotated_phi,
                    &mut rotated_scalar,
                    || false,
                )
                .unwrap();
            assert_eq!(
                rotated.filter_width().to_bits(),
                flow.filter_width().to_bits()
            );
            for (c, &original) in source.iter().enumerate() {
                assert!((rotated.viscosities()[c] - flow.viscosities()[original]).abs() < 2e-13);
                assert!(
                    (rotated_scalar.diffusivities()[c] - scalar.diffusivities()[original]).abs()
                        < 2e-13
                );
            }
            assert!(
                (report.0.minimum_raw_coefficient_squared
                    - expected_report.0.minimum_raw_coefficient_squared)
                    .abs()
                    < 2e-13
            );
            assert!(
                (report.1.minimum_raw_coefficient - expected_report.1.minimum_raw_coefficient)
                    .abs()
                    < 2e-13
            );
        }
    }
}

#[test]
fn cartesian_extreme_length_units_reuse_current_scalar_geometry() {
    let shape = [5, 6, 7];
    let spacing = [0.125, 0.25, 0.5];
    let (rho, u, phi) = cartesian_field(shape);
    let (mc, sc) = settings(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
    let bytes = scalar.workspace_bytes();
    let allocation = scalar.diffusivities().as_ptr();
    let reference = flow
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
        .unwrap();
    let alpha = scalar.diffusivities().to_vec();
    for (hs, us, rs, ps) in [
        (1e-250, 1e250, 1e-150, 1e150),
        (1e250, 1e-250, 1e150, -1e-150),
        (4.0, 2.0, 3.0, -5.0),
    ] {
        let mut scaled =
            PeriodicDynamicSgs::new_cartesian(shape, spacing.map(|h| h * hs), usize::MAX).unwrap();
        let r = rho.iter().map(|v| v * rs).collect::<Vec<_>>();
        let v = u.iter().map(|v| v.map(|q| q * us)).collect::<Vec<_>>();
        let p = phi.iter().map(|v| v * ps).collect::<Vec<_>>();
        let report = scaled
            .prepare_with_scalar(mc, sc, &r, &v, &p, &mut scalar, || false)
            .unwrap();
        for (a, b) in scaled.viscosities().iter().zip(flow.viscosities()) {
            assert!((a / (hs * us) - b).abs() < 2e-13);
        }
        for (a, b) in scalar.diffusivities().iter().zip(&alpha) {
            assert!((a / (hs * us) - b).abs() < 2e-13);
        }
        assert!(
            (report.0.maximum_coefficient_squared - reference.0.maximum_coefficient_squared).abs()
                < 2e-13
        );
        assert!((report.1.maximum_coefficient - reference.1.maximum_coefficient).abs() < 2e-13);
        assert_eq!(scalar.diffusivities().as_ptr(), allocation);
        assert_eq!(scalar.workspace_bytes(), bytes);
    }
    // The scalar owner carries no metric: reusing it with a different aspect
    // ratio must agree with a fresh owner on the same current flow geometry.
    let mut changed =
        PeriodicDynamicSgs::new_cartesian(shape, [0.2, 0.7, 0.4], usize::MAX).unwrap();
    let mut fresh = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
    let expected = changed
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut fresh, || false)
        .unwrap();
    let observed = changed
        .prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false)
        .unwrap();
    assert_eq!(observed, expected);
    assert_eq!(scalar.diffusivities(), fresh.diffusivities());
}

#[test]
fn cartesian_underflow_is_refused_instead_of_becoming_an_undetermined_fit() {
    let shape = [7, 8, 9];
    let len = shape.iter().product();
    let rho = vec![1.0; len];
    let mut u = vec![[0.0; 3]; len];
    let mut phi = vec![0.0; len];
    for c in 0..len {
        let x = core::f64::consts::TAU * (c % shape[0]) as f64 / shape[0] as f64;
        let y = core::f64::consts::TAU * (c / shape[0] % shape[1]) as f64 / shape[1] as f64;
        u[c][0] = y.sin();
        phi[c] = (x + y).sin();
    }
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let (mc, sc) = settings(averaging);
        for spacing in [[1e100, 1e100, 1e-200], [1e-200, 1e200, 1.0]] {
            let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
            let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
            assert!(matches!(
                flow.prepare_with_scalar(mc, sc, &rho, &u, &phi, &mut scalar, || false),
                Err(DynamicSgsError::Cell { .. })
            ));
            // A mathematically zero velocity strain remains a valid zero fit
            // on the very same geometry; the refusal is numerical, not an
            // arbitrary aspect-ratio limit.
            let (momentum, heat) = flow
                .prepare_with_scalar(
                    mc,
                    sc,
                    &rho,
                    &vec![[1.0; 3]; len],
                    &phi,
                    &mut scalar,
                    || false,
                )
                .unwrap();
            assert_eq!(momentum.undetermined_cells, len);
            assert_eq!(heat.undetermined_cells, len);
            assert!(flow.viscosities().iter().all(|v| *v == 0.0));
            assert!(scalar.diffusivities().iter().all(|v| *v == 0.0));
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

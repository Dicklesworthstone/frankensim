use super::*;

fn settings() -> (ChannelDynamicSgsConfig, ChannelDynamicScalarConfig) {
    (
        ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        ChannelDynamicScalarConfig {
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    )
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut h = ((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn fields(shape: [usize; 3], normal: usize) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let cells = shape.iter().product();
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut phi = Vec::new();
    for c in 0..cells {
        let p = [
            c % shape[0],
            c / shape[0] % shape[1],
            c / (shape[0] * shape[1]),
        ];
        let eta = (p[normal] as f64 + 0.5) / shape[normal] as f64;
        let envelope = 4.0 * eta * (1.0 - eta);
        rho.push(1.0 + 0.15 * noise(c, 911));
        u.push(core::array::from_fn(|j| {
            0.02 * envelope * noise(c, 31 * j as u64 + 7)
        }));
        phi.push(1.0 - eta + 0.05 * envelope * noise(c, 101));
    }
    (rho, u, phi)
}

fn xyz(shape: [usize; 3], cell: usize) -> [usize; 3] {
    [
        cell % shape[0],
        cell / shape[0] % shape[1],
        cell / (shape[0] * shape[1]),
    ]
}

fn at(shape: [usize; 3], p: [usize; 3]) -> usize {
    ((p[2] % shape[2]) * shape[1] + p[1] % shape[1]) * shape[0] + p[0] % shape[0]
}

// Physical-unit reference uses full tensors, raw moments, ordinary cell sums,
// and explicit stencils. It shares neither production normalization nor filters.
#[allow(clippy::too_many_arguments)]
fn derivative(
    shape: [usize; 3],
    spacing: [f64; 3],
    normal: usize,
    field: &[f64],
    walls: [ScalarWall; 2],
    filtered: bool,
    cell: usize,
    axis: usize,
) -> f64 {
    let p = xyz(shape, cell);
    if axis == normal && (p[axis] == 0 || p[axis] == shape[axis] - 1) {
        let low = p[axis] == 0;
        let side = if low { 0 } else { 1 };
        let direction = if low { 1.0 } else { -1.0 };
        let mut next = p;
        next[axis] = if low { 1 } else { shape[axis] - 2 };
        let f = field[cell];
        let g = field[at(shape, next)];
        return direction / spacing[axis]
            * match walls[side] {
                ScalarWall::Value(w) => (-8.0 * w + 6.0 * f + 2.0 * g) / 6.0,
                ScalarWall::ZeroGradient if !filtered => (g - f) / 2.0,
                ScalarWall::ZeroGradient => {
                    let mut far = p;
                    far[axis] = if low { 2 } else { shape[axis] - 3 };
                    (-3.0 * f + 4.0 * g - field[at(shape, far)]) / 2.0
                }
            };
    }
    let mut low = p;
    let mut high = p;
    low[axis] += shape[axis] - 1;
    high[axis] += 1;
    (field[at(shape, high)] - field[at(shape, low)]) / (2.0 * spacing[axis])
}

#[allow(clippy::too_many_arguments)]
fn reference(
    shape: [usize; 3],
    spacing: [f64; 3],
    normal: usize,
    rho: &[f64],
    u: &[[f64; 3]],
    wall_u: [[f64; 3]; 2],
    phi: &[f64],
    walls: [ScalarWall; 2],
) -> (Vec<f64>, Vec<f64>) {
    let n = u.len();
    let delta2 = spacing.iter().product::<f64>().cbrt().powi(2);
    let test_delta2 = (9.0 * spacing.iter().product::<f64>()).cbrt().powi(2);
    let tangent = (0..3).filter(|j| *j != normal).collect::<Vec<_>>();
    let samples = |cell| {
        let base = xyz(shape, cell);
        let mut result = Vec::new();
        for q in 0..3 {
            for p in 0..3 {
                let mut point = base;
                point[tangent[0]] += shape[tangent[0]] + p - 1;
                point[tangent[1]] += shape[tangent[1]] + q - 1;
                result.push(at(shape, point));
            }
        }
        result
    };
    let strain = |velocity: &[[f64; 3]]| {
        let components: [Vec<f64>; 3] =
            core::array::from_fn(|j| velocity.iter().map(|u| u[j]).collect());
        (0..n)
            .map(|c| {
                let g: [[f64; 3]; 3] = core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        derivative(
                            shape,
                            spacing,
                            normal,
                            &components[i],
                            wall_u.map(|w| ScalarWall::Value(w[i])),
                            false,
                            c,
                            j,
                        )
                    })
                });
                let mut s: [[f64; 3]; 3] =
                    core::array::from_fn(|i| core::array::from_fn(|j| (g[i][j] + g[j][i]) / 2.0));
                let rate = (2.0 * s.iter().flatten().map(|x| x * x).sum::<f64>()).sqrt();
                let third = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
                for i in 0..3 {
                    s[i][i] -= third;
                }
                (rate, s)
            })
            .collect::<Vec<_>>()
    };
    let base = strain(u);
    let mut mean_rho = vec![0.0; n];
    let mut mean_u = vec![[0.0; 3]; n];
    let mut mean_phi = vec![0.0; n];
    for c in 0..n {
        for j in samples(c) {
            mean_rho[c] += rho[j] / 9.0;
            mean_phi[c] += rho[j] * phi[j] / 9.0;
            for axis in 0..3 {
                mean_u[c][axis] += rho[j] * u[j][axis] / 9.0;
            }
        }
        mean_phi[c] /= mean_rho[c];
        for axis in 0..3 {
            mean_u[c][axis] /= mean_rho[c];
        }
    }
    let coarse = strain(&mean_u);
    let mut momentum = vec![[0.0; 2]; n];
    let mut scalar = vec![[0.0; 2]; n];
    for c in 0..n {
        let mut l = [[0.0; 3]; 3];
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                for j in samples(c) {
                    l[i][k] += rho[j] * u[j][i] * u[j][k] / 9.0;
                    m[i][k] += 2.0 * rho[j] * delta2 * base[j].0 * base[j].1[i][k] / 9.0;
                }
                l[i][k] -= mean_rho[c] * mean_u[c][i] * mean_u[c][k];
                m[i][k] -= 2.0 * mean_rho[c] * test_delta2 * coarse[c].0 * coarse[c].1[i][k];
            }
            let mut ls = -mean_rho[c] * mean_u[c][i] * mean_phi[c];
            let mut ms = -mean_rho[c]
                * test_delta2
                * coarse[c].0
                * derivative(shape, spacing, normal, &mean_phi, walls, true, c, i);
            for j in samples(c) {
                ls += rho[j] * u[j][i] * phi[j] / 9.0;
                ms += rho[j]
                    * delta2
                    * base[j].0
                    * derivative(shape, spacing, normal, phi, walls, false, j, i)
                    / 9.0;
            }
            scalar[c][0] += ls * ms;
            scalar[c][1] += ms * ms;
        }
        let third = (l[0][0] + l[1][1] + l[2][2]) / 3.0;
        for i in 0..3 {
            l[i][i] -= third;
        }
        for i in 0..3 {
            for j in 0..3 {
                momentum[c][0] += l[i][j] * m[i][j];
                momentum[c][1] += m[i][j] * m[i][j];
            }
        }
    }
    let fit = |pairs: Vec<[f64; 2]>| {
        let mut planes = vec![[0.0; 2]; shape[normal]];
        for (c, pair) in pairs.iter().enumerate() {
            let layer = xyz(shape, c)[normal];
            planes[layer][0] += pair[0];
            planes[layer][1] += pair[1];
        }
        (0..n)
            .map(|c| {
                let p = planes[xyz(shape, c)[normal]];
                let coefficient = if p[1] == 0.0 {
                    0.0
                } else {
                    (p[0] / p[1]).max(0.0)
                };
                coefficient * delta2 * base[c].0
            })
            .collect()
    };
    (fit(momentum), fit(scalar))
}

#[test]
fn channel_favre_momentum_and_heat_match_independent_physical_reference() {
    let shape = [5, 4, 5];
    let spacing = [0.3, 0.7, 1.1];
    let (mc, sc) = settings();
    for normal in 0..3 {
        let (rho, u, phi) = fields(shape, normal);
        let mut wall_u = [[0.03, -0.01, 0.02], [-0.02, 0.015, -0.01]];
        wall_u[0][normal] = 0.0;
        wall_u[1][normal] = 0.0;
        for walls in [
            [ScalarWall::Value(1.0), ScalarWall::Value(0.0)],
            [ScalarWall::ZeroGradient; 2],
            [ScalarWall::Value(1.0), ScalarWall::ZeroGradient],
        ] {
            let mut flow =
                ChannelDynamicSgs::new_cartesian(shape, spacing, normal, usize::MAX).unwrap();
            let mut scalar = ChannelDynamicScalar::new(shape, normal, usize::MAX).unwrap();
            flow.prepare_with_scalar(mc, sc, &rho, &u, wall_u, &phi, walls, &mut scalar, || false)
                .unwrap();
            let expected = reference(shape, spacing, normal, &rho, &u, wall_u, &phi, walls);
            for (actual, expected) in flow
                .viscosities()
                .iter()
                .zip(expected.0)
                .chain(scalar.diffusivities().iter().zip(expected.1))
            {
                assert!(
                    (actual - expected).abs() < 3e-13,
                    "axis {normal}: {actual} != {expected}"
                );
            }
        }
    }
}

#[test]
fn channel_fit_has_independent_heat_and_respects_boosts_and_scalar_offsets() {
    let shape = [5, 4, 5];
    let (mut rho, u, phi) = fields(shape, 2);
    rho.fill(1.0);
    let (mc, sc) = settings();
    let walls = [ScalarWall::Value(1.0), ScalarWall::Value(0.0)];
    let mut flow = ChannelDynamicSgs::new(shape, 1.0, 2, usize::MAX).unwrap();
    let mut scalar = ChannelDynamicScalar::new(shape, 2, usize::MAX).unwrap();
    let report = flow
        .prepare_with_scalar(
            mc,
            sc,
            &rho,
            &u,
            [[0.0; 3]; 2],
            &phi,
            walls,
            &mut scalar,
            || false,
        )
        .unwrap();
    assert!(report.0.maximum_eddy_viscosity > 0.0 && report.1.maximum_diffusivity > 0.0);
    assert!(
        flow.viscosities()
            .iter()
            .zip(scalar.diffusivities())
            .any(|(nu, alpha)| *nu == 0.0 && *alpha > 0.0)
    );
    let expected_u = flow.viscosities().to_vec();
    let expected_t = scalar.diffusivities().to_vec();
    let boost = [4.0, -2.0, 0.0];
    let shifted: Vec<[f64; 3]> = u
        .iter()
        .map(|v| core::array::from_fn(|j| v[j] + boost[j]))
        .collect();
    let shifted_phi: Vec<_> = phi.iter().map(|t| t + 8.0).collect();
    flow.prepare_with_scalar(
        mc,
        sc,
        &rho,
        &shifted,
        [boost; 2],
        &shifted_phi,
        [ScalarWall::Value(9.0), ScalarWall::Value(8.0)],
        &mut scalar,
        || false,
    )
    .unwrap();
    for (a, b) in flow
        .viscosities()
        .iter()
        .zip(expected_u)
        .chain(scalar.diffusivities().iter().zip(expected_t))
    {
        assert!((a - b).abs() < 3e-13);
    }
}

#[test]
fn channel_walls_participate_in_uniform_field_normalization() {
    let shape = [3, 3, 2];
    let cells = 18;
    let (mc, sc) = settings();
    let mut flow = ChannelDynamicSgs::new(shape, 1.0, 2, usize::MAX).unwrap();
    let mut scalar = ChannelDynamicScalar::new(shape, 2, usize::MAX).unwrap();
    let reports = flow
        .prepare_with_scalar(
            mc,
            sc,
            &vec![1.0; cells],
            &vec![[0.0; 3]; cells],
            [[1.0, 0.0, 0.0], [0.0; 3]],
            &vec![0.0; cells],
            [ScalarWall::Value(1.0), ScalarWall::Value(0.0)],
            &mut scalar,
            || false,
        )
        .unwrap();
    assert_eq!(flow.inner.strain_rate[0], 4.0 / 3.0);
    assert!(reports.0.undetermined_cells < cells && reports.1.undetermined_cells < cells);
    assert!(
        flow.viscosities()
            .iter()
            .chain(scalar.diffusivities())
            .all(|v| *v == 0.0)
    );
}

#[test]
fn channel_test_box_and_quadratic_wall_derivatives_preserve_geometry() {
    let shape = [4, 5, 6];
    let h = 0.7;
    for normal in 0..3 {
        let geometry = FilterGeometry::Channel(normal);
        let length = shape[normal] as f64 * h;
        let polynomial = |y: f64| 3.0 - 2.0 * y + 0.4 * y * y;
        let values: Vec<_> = (0..shape.iter().product())
            .map(|c| polynomial((xyz(shape, c)[normal] as f64 + 0.5) * h))
            .collect();
        for c in 0..values.len() {
            let mut count = 0;
            geometry.each_box(shape, c, |j| {
                count += 1;
                assert_eq!(xyz(shape, j)[normal], xyz(shape, c)[normal]);
            });
            assert_eq!(count, 9);
            let y = (xyz(shape, c)[normal] as f64 + 0.5) * h;
            let d = geometry
                .stencil(shape, [1.0 / h; 3], c, xyz(shape, c), normal)
                .derivative(
                    |j| values[j],
                    [
                        ScalarWall::Value(polynomial(0.0)),
                        ScalarWall::Value(polynomial(length)),
                    ],
                )
                .unwrap();
            assert!((d - (-2.0 + 0.8 * y)).abs() < 3e-14);
        }
        assert!((geometry.ratio_squared() - 9.0_f64.powf(2.0 / 3.0)).abs() < 2e-14);
    }
}

#[test]
fn filtered_neumann_favre_gradient_is_not_assigned_zero_wall_derivative() {
    // rho=1+a(x)y and phi=b(x), with mean(a)=mean(b)=0 and
    // mean(a*b)=1/6. Physical dphi/dy=0, yet Favre(phi)=y/6.
    // This is an exact manufactured covariance contribution at the wall.
    let shape = [3, 3, 5];
    let h = 0.2;
    let filtered: Vec<_> = (0..45)
        .map(|c| {
            let y = (xyz(shape, c)[2] as f64 + 0.5) * h;
            let mut moment = 0.0;
            let mut density = 0.0;
            for a in [-0.5, 0.0, 0.5] {
                moment += (1.0 + a * y) * a / 3.0;
                density += (1.0 + a * y) / 3.0;
            }
            moment / density
        })
        .collect();
    let g = FilterGeometry::Channel(2)
        .stencil(shape, [1.0 / h; 3], 0, xyz(shape, 0), 2)
        .filtered_scalar_derivative(&filtered, [ScalarWall::ZeroGradient; 2])
        .unwrap();
    assert!((g - 1.0 / 6.0).abs() < 2e-15);
    let errors: Vec<_> = [0.2, 0.1, 0.05]
        .into_iter()
        .map(|h| {
            let field: Vec<_> = (0..45)
                .map(|c| {
                    let y = (xyz(shape, c)[2] as f64 + 0.5) * h;
                    y / 6.0 + y.powi(3)
                })
                .collect();
            let g = FilterGeometry::Channel(2)
                .stencil(shape, [1.0 / h; 3], 0, xyz(shape, 0), 2)
                .filtered_scalar_derivative(&field, [ScalarWall::ZeroGradient; 2])
                .unwrap();
            (g - (1.0 / 6.0 + 3.0 * (0.5 * h).powi(2))).abs()
        })
        .collect();
    assert!((errors[0] / errors[1] - 4.0).abs() < 1e-12);
    assert!((errors[1] / errors[2] - 4.0).abs() < 1e-12);
}

#[test]
fn channel_workspace_caps_cancellation_and_boundary_admission() {
    let shape = [3, 3, 2];
    let cells = 18;
    let (mc, sc) = settings();
    assert!(ChannelDynamicSgs::new(shape, 1.0, 3, usize::MAX).is_err());
    assert!(ChannelDynamicSgs::new([2, 3, 2], 1.0, 2, usize::MAX).is_err());
    let bytes = ChannelDynamicSgs::required_bytes(cells).unwrap();
    assert_eq!(bytes, cells * 136);
    assert!(ChannelDynamicSgs::new(shape, 1.0, 2, bytes - 1).is_err());
    let mut flow = ChannelDynamicSgs::new(shape, 1.0, 2, bytes).unwrap();
    let mut scalar = ChannelDynamicScalar::new(shape, 2, cells * 64).unwrap();
    assert_eq!(
        flow.workspace_bytes() + scalar.workspace_bytes(),
        cells * 200
    );
    let rho = vec![1.0; cells];
    let u = vec![[0.0; 3]; cells];
    let phi = vec![0.0; cells];
    assert!(matches!(
        flow.prepare(mc, &rho, &u, [[0.0; 3]; 2], || true),
        Err(DynamicSgsError::Cancelled)
    ));
    assert!(
        flow.prepare(mc, &rho, &u, [[0.0, 0.0, 0.01]; 2], || false)
            .is_err()
    );
    assert!(
        flow.prepare_with_scalar(
            mc,
            sc,
            &rho,
            &u,
            [[0.0; 3]; 2],
            &phi,
            [ScalarWall::ZeroGradient; 2],
            &mut scalar,
            || false
        )
        .is_err()
    );
    assert!(
        flow.prepare_with_scalar(
            mc,
            sc,
            &rho,
            &u,
            [[0.0; 3]; 2],
            &phi,
            [ScalarWall::Value(f64::NAN), ScalarWall::Value(0.0)],
            &mut scalar,
            || false
        )
        .is_err()
    );
}

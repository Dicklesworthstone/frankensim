use super::*;
use crate::sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicSgsConfig};
use crate::sgs3::dynamic::scalar::{DynamicScalarConfig, PeriodicDynamicScalar};

fn policy() -> DynamicWaleConfig {
    DynamicWaleConfig {
        sensor_threshold: 0.09,
        near_wall_coefficient: 0.5,
    }
}

fn controls(averaging: DynamicAveraging) -> DynamicSgsConfig {
    DynamicSgsConfig {
        averaging,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-12 * expected.abs().max(1e-5),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

// Independent full 3x3 tensor evaluation, in physical units. It deliberately
// does not call the production WALE kernel, packed contractions, or sensor.
fn tensor(g: [[f64; 3]; 3]) -> (f64, f64, [[f64; 3]; 3]) {
    let mut s = [[0.0; 3]; 3];
    let mut square = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            s[i][j] = (g[i][j] + g[j][i]) / 2.0;
            for k in 0..3 {
                square[i][j] += g[i][k] * g[k][j];
            }
        }
    }
    let trace = (square[0][0] + square[1][1] + square[2][2]) / 3.0;
    let mut a = 0.0;
    let mut b = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let d = (square[i][j] + square[j][i]) / 2.0 - if i == j { trace } else { 0.0 };
            a += s[i][j] * s[i][j];
            b += d * d;
        }
    }
    let w = if b == 0.0 {
        0.0
    } else {
        b.powf(1.5) / (a.powf(2.5) + b.powf(1.25))
    };
    let svs = if a == 0.0 && b == 0.0 {
        0.0
    } else {
        b.powf(1.5) / (b.powf(1.5) + a.powi(3))
    };
    let trace_s = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
    for i in 0..3 {
        s[i][i] -= trace_s;
    }
    (w, svs, s)
}

fn at(shape: [usize; 3], i: usize) -> [usize; 3] {
    [
        i % shape[0],
        i / shape[0] % shape[1],
        i / (shape[0] * shape[1]),
    ]
}

fn flat(shape: [usize; 3], p: [usize; 3]) -> usize {
    p[0] + shape[0] * (p[1] + shape[1] * p[2])
}

fn donors(shape: [usize; 3], cell: usize, normal: Option<usize>) -> Vec<usize> {
    let p = at(shape, cell);
    let mut result = Vec::new();
    for dz in -1_isize..=1 {
        for dy in -1_isize..=1 {
            for dx in -1_isize..=1 {
                let offsets = [dx, dy, dz];
                if normal.is_some_and(|axis| offsets[axis] != 0) {
                    continue;
                }
                let q = core::array::from_fn(|axis| {
                    (p[axis] as isize + offsets[axis]).rem_euclid(shape[axis] as isize) as usize
                });
                result.push(flat(shape, q));
            }
        }
    }
    result
}

fn gradient(
    shape: [usize; 3],
    spacing: [f64; 3],
    normal: Option<usize>,
    walls: [[f64; 3]; 2],
    u: &[[f64; 3]],
    cell: usize,
) -> [[f64; 3]; 3] {
    let p = at(shape, cell);
    core::array::from_fn(|component| {
        core::array::from_fn(|axis| {
            let mut low = p;
            let mut high = p;
            if normal == Some(axis) && (p[axis] == 0 || p[axis] + 1 == shape[axis]) {
                let side = usize::from(p[axis] != 0);
                high[axis] = if side == 0 { 1 } else { shape[axis] - 2 };
                let sign = if side == 0 { 1.0 } else { -1.0 };
                sign * (u[cell][component] + u[flat(shape, high)][component] / 3.0
                    - 4.0 * walls[side][component] / 3.0)
                    / spacing[axis]
            } else {
                low[axis] = (p[axis] + shape[axis] - 1) % shape[axis];
                high[axis] = (p[axis] + 1) % shape[axis];
                (u[flat(shape, high)][component] - u[flat(shape, low)][component])
                    / (2.0 * spacing[axis])
            }
        })
    })
}

fn field(shape: [usize; 3]) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let cells = shape.iter().product();
    let rho = (0..cells)
        .map(|i| 1.0 + ((i * 7 + 3) % 11) as f64 / 32.0)
        .collect();
    let u = (0..cells)
        .map(|i| {
            core::array::from_fn(|a| {
                (((i * (17 + 8 * a) + i * i * (a + 1) + 13 * a) % 127) as f64 - 63.0) / 512.0
            })
        })
        .collect();
    let theta = (0..cells)
        .map(|i| 3.0 + ((i * i * 3 + 7 * i) % 41) as f64 / 64.0)
        .collect();
    (rho, u, theta)
}

// G1: independent physical-space Favre filter, full tensors, and physical
// Delta squared in M. No production derivative/filter/fitting helper is used.
fn reference(
    shape: [usize; 3],
    spacing: [f64; 3],
    normal: Option<usize>,
    walls: [[f64; 3]; 2],
    averaging: DynamicAveraging,
    rho: &[f64],
    u: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let n = rho.len();
    let delta2 = (spacing.iter().product::<f64>()).cbrt().powi(2);
    let ratio2 = if normal.is_some() {
        9.0_f64.powf(2.0 / 3.0)
    } else {
        9.0
    };
    let rates: Vec<_> = (0..n)
        .map(|i| tensor(gradient(shape, spacing, normal, walls, u, i)))
        .collect();
    let samples: Vec<_> = (0..n).map(|i| donors(shape, i, normal)).collect();
    let filtered_rho: Vec<_> = samples
        .iter()
        .map(|cells| cells.iter().map(|&j| rho[j]).sum::<f64>() / cells.len() as f64)
        .collect();
    let filtered_u: Vec<_> = samples
        .iter()
        .enumerate()
        .map(|(i, cells)| {
            core::array::from_fn(|a| {
                cells.iter().map(|&j| rho[j] * u[j][a]).sum::<f64>()
                    / (cells.len() as f64 * filtered_rho[i])
            })
        })
        .collect();
    let pairs: Vec<[f64; 2]> = (0..n)
        .map(|cell| {
            let test = tensor(gradient(shape, spacing, normal, walls, &filtered_u, cell));
            let mut l = [[0.0; 3]; 3];
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for &donor in &samples[cell] {
                        let weight = rho[donor] / samples[cell].len() as f64;
                        l[i][j] += weight
                            * (u[donor][i] - filtered_u[cell][i])
                            * (u[donor][j] - filtered_u[cell][j]);
                        m[i][j] += weight * rates[donor].0 * rates[donor].2[i][j];
                    }
                    m[i][j] = 2.0
                        * delta2
                        * (m[i][j] - ratio2 * filtered_rho[cell] * test.0 * test.2[i][j]);
                }
            }
            let tr = (l[0][0] + l[1][1] + l[2][2]) / 3.0;
            let mut pair = [0.0; 2];
            for i in 0..3 {
                for j in 0..3 {
                    pair[0] += (l[i][j] - if i == j { tr } else { 0.0 }) * m[i][j];
                    pair[1] += m[i][j] * m[i][j];
                }
            }
            pair
        })
        .collect();
    (0..n)
        .map(|cell| {
            let selection: Vec<_> = if let Some(axis) = normal {
                (0..n)
                    .filter(|&j| at(shape, j)[axis] == at(shape, cell)[axis])
                    .collect()
            } else if averaging == DynamicAveraging::Domain {
                (0..n).collect()
            } else {
                samples[cell].clone()
            };
            let pair = selection
                .iter()
                .fold([0.0; 2], |a, &j| [a[0] + pairs[j][0], a[1] + pairs[j][1]]);
            let raw = if pair[1] == 0.0 {
                0.0
            } else {
                pair[0] / pair[1]
            };
            let coefficient = if rates[cell].1 <= policy().sensor_threshold {
                policy().near_wall_coefficient.powi(2)
            } else {
                raw.max(0.0)
            };
            [
                coefficient,
                coefficient * delta2 * rates[cell].0,
                rates[cell].1,
            ]
        })
        .collect()
}

#[test]
fn g0_tensor_sensor_limits_and_g1_near_wall_cubic_rate() {
    for g in [
        [[0.0; 3]; 3],
        [[0.0, 2.0, 0.0], [0.0; 3], [0.0; 3]],
        [[0.0, -2.0, 0.0], [2.0, 0.0, 0.0], [0.0; 3]],
        [[0.3, 0.7, -0.2], [-0.1, -0.4, 0.6], [0.8, -0.5, 0.1]],
    ] {
        let expected = tensor(g);
        close(rate(g, 0).unwrap(), expected.0);
        close(sensor(g), expected.1);
        for scale in [1e-100, 1e100] {
            let scaled = g.map(|row| row.map(|v| v * scale));
            close(sensor(scaled), expected.1);
            close(rate(scaled, 0).unwrap() / scale, expected.0);
        }
    }
    assert_eq!(sensor([[0.0, 2.0, 0.0], [0.0; 3], [0.0; 3]]), 0.0);
    assert_eq!(sensor([[0.0, -2.0, 0.0], [2.0, 0.0, 0.0], [0.0; 3]]), 1.0);
    let values = [0.01, 0.005, 0.0025].map(|y| {
        let g = [
            [0.2 * y, 1.0, 0.0],
            [0.0, -0.5 * y, 0.0],
            [0.0, 0.4, 0.3 * y],
        ];
        assert!(sensor(g) < policy().sensor_threshold);
        policy().near_wall_coefficient.powi(2) * rate(g, 0).unwrap()
    });
    for p in values.windows(2) {
        assert!(((p[0] / p[1]).log2() - 3.0).abs() < 0.01);
    }
}

#[test]
fn g1_periodic_and_all_channel_axes_match_independent_favre_wale_tensors() {
    let shape = [5, 6, 7];
    let spacing = [0.5, 0.75, 1.25];
    let (rho, u, _) = field(shape);
    for normal in [None, Some(0), Some(1), Some(2)] {
        let mut walls = [[0.03, -0.02, 0.01], [-0.01, 0.025, -0.015]];
        if let Some(axis) = normal {
            walls[0][axis] = 0.0;
            walls[1][axis] = 0.0;
        }
        for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
            if normal.is_some() && averaging == DynamicAveraging::LocalBox {
                continue;
            }
            let expected = reference(shape, spacing, normal, walls, averaging, &rho, &u);
            let (coefficients, viscosities) = if let Some(axis) = normal {
                let mut flow = ChannelDynamicSgs::new_cartesian(shape, spacing, axis, usize::MAX)
                    .unwrap()
                    .with_wale(policy(), usize::MAX)
                    .unwrap();
                flow.prepare(
                    ChannelDynamicSgsConfig {
                        maximum_coefficient_squared: 100.0,
                        maximum_eddy_to_molecular: 1000.0,
                        usage: SgsUse::Exploration,
                    },
                    &rho,
                    &u,
                    walls,
                    || false,
                )
                .unwrap();
                (
                    flow.coefficients_squared().collect::<Vec<_>>(),
                    flow.viscosities().to_vec(),
                )
            } else {
                let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX)
                    .unwrap()
                    .with_wale(policy(), usize::MAX)
                    .unwrap();
                flow.prepare(controls(averaging), &rho, &u, || false)
                    .unwrap();
                (
                    flow.coefficients_squared().collect::<Vec<_>>(),
                    flow.viscosities().to_vec(),
                )
            };
            assert!(expected.iter().any(|v| v[2] <= policy().sensor_threshold));
            assert!(expected.iter().any(|v| v[2] > policy().sensor_threshold));
            assert!(viscosities.iter().any(|v| *v > 0.0));
            for i in 0..rho.len() {
                close(coefficients[i], expected[i][0]);
                close(viscosities[i], expected[i][1]);
            }
        }
    }
}

#[test]
fn g3_independent_heat_fit_keeps_its_original_bits_under_wale_momentum() {
    let shape = [5, 6, 7];
    let (rho, u, theta) = field(shape);
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let mut smag =
            PeriodicDynamicSgs::new_cartesian(shape, [0.5, 0.75, 1.25], usize::MAX).unwrap();
        let mut wale = PeriodicDynamicSgs::new_cartesian(shape, [0.5, 0.75, 1.25], usize::MAX)
            .unwrap()
            .with_wale(policy(), usize::MAX)
            .unwrap();
        let mut a = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
        let mut b = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
        let heat = DynamicScalarConfig {
            averaging,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        };
        let (_, expected) = smag
            .prepare_with_scalar(controls(averaging), heat, &rho, &u, &theta, &mut a, || {
                false
            })
            .unwrap();
        let (_, observed) = wale
            .prepare_with_scalar(controls(averaging), heat, &rho, &u, &theta, &mut b, || {
                false
            })
            .unwrap();
        assert_eq!(observed, expected);
        assert_eq!(a.diffusivities(), b.diffusivities());
        assert_eq!(
            a.coefficients().collect::<Vec<_>>(),
            b.coefficients().collect::<Vec<_>>()
        );
        assert_ne!(smag.viscosities(), wale.viscosities());
    }
}

#[test]
fn g3_dynamic_wale_preserves_boost_and_length_velocity_density_units() {
    let shape = [5, 6, 7];
    let h = [0.5, 0.75, 1.25];
    let (rho, u, _) = field(shape);
    let cfg = controls(DynamicAveraging::LocalBox);
    let mut base = PeriodicDynamicSgs::new_cartesian(shape, h, usize::MAX)
        .unwrap()
        .with_wale(policy(), usize::MAX)
        .unwrap();
    base.prepare(cfg, &rho, &u, || false).unwrap();
    for (length, velocity, density) in
        [(1.0, 1.0, 1.0), (1e-80, 1e60, 1e150), (1e80, 1e-60, 1e-150)]
    {
        let r: Vec<_> = rho.iter().map(|r| r * density).collect();
        let v: Vec<_> = u
            .iter()
            .map(|row| row.map(|value| (value + 8.0) * velocity))
            .collect();
        let mut scaled =
            PeriodicDynamicSgs::new_cartesian(shape, h.map(|v| v * length), usize::MAX)
                .unwrap()
                .with_wale(policy(), usize::MAX)
                .unwrap();
        scaled.prepare(cfg, &r, &v, || false).unwrap();
        for i in 0..rho.len() {
            close(
                scaled.coefficient_squared(i).unwrap(),
                base.coefficient_squared(i).unwrap(),
            );
            close(
                scaled.viscosities()[i] / (length * velocity),
                base.viscosities()[i],
            );
        }
    }
}

#[test]
fn g0_uniform_and_laminar_shear_clear_stale_wale_output() {
    let shape = [5, 6, 7];
    let (rho, u, _) = field(shape);
    let cfg = controls(DynamicAveraging::LocalBox);
    let mut flow = PeriodicDynamicSgs::new(shape, 1.0, usize::MAX)
        .unwrap()
        .with_wale(policy(), usize::MAX)
        .unwrap();
    flow.prepare(cfg, &rho, &u, || false).unwrap();
    assert!(flow.viscosities().iter().any(|v| *v > 0.0));
    for velocity in [
        vec![[0.25, 0.5, -0.125]; rho.len()],
        (0..rho.len())
            .map(|i| [at(shape, i)[1] as f64 / 16.0, 0.0, 0.0])
            .collect(),
    ] {
        flow.prepare(cfg, &vec![1.0; rho.len()], &velocity, || false)
            .unwrap();
        assert!(flow.viscosities().iter().all(|v| *v == 0.0));
        assert!(flow.coefficients_squared().all(|c| c == 0.25));
    }
}

#[test]
fn g0_wale_admission_and_g4_reuse_after_every_cancellation_poll() {
    let shape = [4; 3];
    let n = 64;
    let (rho, u, _) = field(shape);
    assert_eq!(PeriodicDynamicSgs::required_bytes(n).unwrap(), 136 * n);
    assert_eq!(PeriodicDynamicSgs::required_wale_bytes(n).unwrap(), 152 * n);
    assert!(matches!(
        PeriodicDynamicSgs::new(shape, 1.0, usize::MAX)
            .unwrap()
            .with_wale(policy(), 152 * n - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    for bad in [
        DynamicWaleConfig {
            sensor_threshold: 0.0,
            ..policy()
        },
        DynamicWaleConfig {
            sensor_threshold: 1.0,
            ..policy()
        },
        DynamicWaleConfig {
            sensor_threshold: f64::NAN,
            ..policy()
        },
        DynamicWaleConfig {
            near_wall_coefficient: -0.5,
            ..policy()
        },
        DynamicWaleConfig {
            near_wall_coefficient: 1e-200,
            ..policy()
        },
        DynamicWaleConfig {
            near_wall_coefficient: 1e200,
            ..policy()
        },
    ] {
        assert!(bad.validate().is_err());
    }
    let mut flow = PeriodicDynamicSgs::new(shape, 1.0, usize::MAX)
        .unwrap()
        .with_wale(policy(), 152 * n)
        .unwrap();
    let cfg = controls(DynamicAveraging::LocalBox);
    assert!(matches!(
        flow.prepare(
            DynamicSgsConfig {
                maximum_coefficient_squared: 0.24,
                ..cfg
            },
            &rho,
            &u,
            || false
        ),
        Err(DynamicSgsError::Invalid(_))
    ));
    let mut polls = 0;
    let expected = flow
        .prepare(cfg, &rho, &u, || {
            polls += 1;
            false
        })
        .unwrap();
    let nu = flow.viscosities().to_vec();
    let coefficients: Vec<_> = flow.coefficients_squared().collect();
    let pointers = (
        flow.viscosity.as_ptr(),
        flow.basis.as_ptr(),
        flow.wale.as_ref().unwrap().response.as_ptr(),
    );
    for stop in 1..=polls {
        let mut seen = 0;
        assert_eq!(
            flow.prepare(cfg, &rho, &u, || {
                seen += 1;
                seen == stop
            }),
            Err(DynamicSgsError::Cancelled)
        );
        assert_eq!(flow.prepare(cfg, &rho, &u, || false).unwrap(), expected);
        assert_eq!(flow.viscosities(), nu);
        assert_eq!(
            flow.coefficients_squared().collect::<Vec<_>>(),
            coefficients
        );
        assert_eq!(
            pointers,
            (
                flow.viscosity.as_ptr(),
                flow.basis.as_ptr(),
                flow.wale.as_ref().unwrap().response.as_ptr()
            )
        );
    }
    assert_eq!(flow.workspace_bytes(), 152 * n);
}

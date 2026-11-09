use super::*;

fn controls(averaging: DynamicAveraging) -> (DynamicSgsConfig, DynamicScalarConfig) {
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

fn donors(shape: [usize; 3], normal: usize, cell: usize) -> Vec<usize> {
    let p = at(shape, cell);
    let mut result = Vec::new();
    for z in -1_isize..=1 {
        for y in -1_isize..=1 {
            for x in -1_isize..=1 {
                let d = [x, y, z];
                if d[normal] != 0 {
                    continue;
                }
                result.push(flat(
                    shape,
                    core::array::from_fn(|a| {
                        (p[a] as isize + d[a]).rem_euclid(shape[a] as isize) as usize
                    }),
                ));
            }
        }
    }
    result
}

fn derivative(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    values: &[f64],
    ghosts: [f64; 2],
    cell: usize,
    axis: usize,
) -> f64 {
    let p = at(shape, cell);
    let mut low = p;
    let mut high = p;
    low[axis] = (p[axis] + shape[axis] - 1) % shape[axis];
    high[axis] = (p[axis] + 1) % shape[axis];
    let a = if axis == normal && p[axis] == 0 {
        ghosts[0]
    } else {
        values[flat(shape, low)]
    };
    let b = if axis == normal && p[axis] + 1 == shape[axis] {
        ghosts[1]
    } else {
        values[flat(shape, high)]
    };
    (b - a) / (2.0 * h[axis])
}

fn tensors(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    velocity: &[[f64; 3]],
    ghosts: [[f64; 3]; 2],
) -> Vec<(f64, [[f64; 3]; 3])> {
    let components: [Vec<_>; 3] = core::array::from_fn(|a| velocity.iter().map(|u| u[a]).collect());
    (0..velocity.len())
        .map(|cell| {
            let g: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    derivative(
                        shape,
                        h,
                        normal,
                        &components[i],
                        ghosts.map(|v| v[i]),
                        cell,
                        j,
                    )
                })
            });
            let mut s: [[f64; 3]; 3] =
                core::array::from_fn(|i| core::array::from_fn(|j| 0.5 * (g[i][j] + g[j][i])));
            let rate = (2.0 * s.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt();
            let tr = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
            for i in 0..3 {
                s[i][i] -= tr;
            }
            (rate, s)
        })
        .collect()
}

fn field(shape: [usize; 3]) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let n = shape.iter().product();
    (
        (0..n)
            .map(|i| 1.0 + ((i * 7 + 3) % 11) as f64 / 32.0)
            .collect(),
        (0..n)
            .map(|i| {
                core::array::from_fn(|a| {
                    (((i * (17 + 8 * a) + i * i * (a + 1) + 13 * a) % 127) as f64 - 63.0) / 512.0
                })
            })
            .collect(),
        (0..n)
            .map(|i| 3.0 + ((i * i * 3 + 7 * i) % 41) as f64 / 64.0)
            .collect(),
    )
}

// Independent physical-unit full-tensor Germano and vector scalar identities.
// No production gradient, packed tensor, filter or fitting helper is called.
fn reference(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    averaging: DynamicAveraging,
    rho: &[f64],
    u: &[[f64; 3]],
    ghosts: [[f64; 3]; 2],
    phi: &[f64],
    scalar_ghosts: [f64; 2],
) -> Vec<[f64; 4]> {
    let n = rho.len();
    let delta2 = h.iter().product::<f64>().cbrt().powi(2);
    let ratio = 9.0_f64.powf(2.0 / 3.0);
    let samples: Vec<_> = (0..n).map(|i| donors(shape, normal, i)).collect();
    let r: Vec<_> = samples
        .iter()
        .map(|d| d.iter().map(|&j| rho[j]).sum::<f64>() / 9.0)
        .collect();
    let v: Vec<_> = samples
        .iter()
        .enumerate()
        .map(|(i, d)| {
            core::array::from_fn(|a| {
                d.iter().map(|&j| rho[j] * u[j][a]).sum::<f64>() / (9.0 * r[i])
            })
        })
        .collect();
    let p: Vec<_> = samples
        .iter()
        .enumerate()
        .map(|(i, d)| d.iter().map(|&j| rho[j] * phi[j]).sum::<f64>() / (9.0 * r[i]))
        .collect();
    let s = tensors(shape, h, normal, u, ghosts);
    let sh = tensors(shape, h, normal, &v, ghosts);
    let grad: Vec<[f64; 3]> = (0..n)
        .map(|i| core::array::from_fn(|a| derivative(shape, h, normal, phi, scalar_ghosts, i, a)))
        .collect();
    let pairs: Vec<[f64; 4]> = (0..n)
        .map(|i| {
            let mut l = [[0.0; 3]; 3];
            let mut m = [[0.0; 3]; 3];
            let mut q = [0.0; 3];
            let mut b = [0.0; 3];
            for &j in &samples[i] {
                let weight = rho[j] / 9.0;
                for a in 0..3 {
                    q[a] += weight * (u[j][a] - v[i][a]) * (phi[j] - p[i]);
                    b[a] += weight * s[j].0 * grad[j][a];
                    for c in 0..3 {
                        l[a][c] += weight * (u[j][a] - v[i][a]) * (u[j][c] - v[i][c]);
                        m[a][c] += weight * s[j].0 * s[j].1[a][c];
                    }
                }
            }
            let tr = (l[0][0] + l[1][1] + l[2][2]) / 3.0;
            let mut pair = [0.0; 4];
            for a in 0..3 {
                let m_phi = delta2
                    * (b[a]
                        - ratio
                            * r[i]
                            * sh[i].0
                            * derivative(shape, h, normal, &p, scalar_ghosts, i, a));
                pair[2] += q[a] * m_phi;
                pair[3] += m_phi * m_phi;
                for c in 0..3 {
                    let stress = 2.0 * delta2 * (m[a][c] - ratio * r[i] * sh[i].0 * sh[i].1[a][c]);
                    pair[0] += (l[a][c] - if a == c { tr } else { 0.0 }) * stress;
                    pair[1] += stress * stress;
                }
            }
            pair
        })
        .collect();
    (0..n)
        .map(|i| {
            let selection = if averaging == DynamicAveraging::Domain {
                (0..n)
                    .filter(|&j| at(shape, j)[normal] == at(shape, i)[normal])
                    .collect()
            } else {
                samples[i].clone()
            };
            let pair = selection.iter().fold([0.0; 4], |mut sum, &j| {
                for k in 0..4 {
                    sum[k] += pairs[j][k];
                }
                sum
            });
            let c = if pair[1] == 0.0 {
                0.0
            } else {
                (pair[0] / pair[1]).max(0.0)
            };
            let cp = if pair[3] == 0.0 {
                0.0
            } else {
                (pair[2] / pair[3]).max(0.0)
            };
            [c, cp, c * delta2 * s[i].0, cp * delta2 * s[i].0]
        })
        .collect()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 5e-12 * b.abs().max(1e-5),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g1_reservoir_momentum_and_heat_match_independent_physical_tensors_on_every_axis() {
    let shape = [5, 6, 7];
    let h = [0.5, 0.75, 1.25];
    let (rho, u, phi) = field(shape);
    let ghosts = [[0.09, -0.04, 0.02], [-0.02, 0.03, -0.05]];
    let temperature = [2.7, 3.8];
    let mut active_momentum = false;
    let mut active_heat = false;
    for normal in 0..3 {
        for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
            let (mc, hc) = controls(averaging);
            let expected = reference(
                shape,
                h,
                normal,
                averaging,
                &rho,
                &u,
                ghosts,
                &phi,
                temperature,
            );
            let mut flow = ReservoirDynamicSgs::new(shape, h, normal, usize::MAX).unwrap();
            let mut scalar = ReservoirDynamicScalar::new(shape, normal, usize::MAX).unwrap();
            flow.prepare_with_scalar(
                mc,
                hc,
                &rho,
                &u,
                ghosts,
                &phi,
                temperature,
                &mut scalar,
                || false,
            )
            .unwrap();
            for i in 0..rho.len() {
                close(flow.coefficient_squared(i).unwrap(), expected[i][0]);
                close(scalar.coefficient(i).unwrap(), expected[i][1]);
                close(flow.viscosities()[i], expected[i][2]);
                close(scalar.diffusivities()[i], expected[i][3]);
            }
            active_momentum |= flow.viscosities().iter().any(|v| *v > 0.0);
            active_heat |= scalar.diffusivities().iter().any(|v| *v > 0.0);
        }
    }
    assert!(active_momentum && active_heat);
}

#[test]
fn g1_full_spacing_ghosts_resolve_quadratic_shear_and_exact_laminar_zero_fit() {
    let shape = [5; 3];
    let h = [0.5, 0.75, 1.25];
    let rho = vec![1.0; 125];
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        let value = |position: f64| {
            let mut u = [0.0; 3];
            u[normal] = 0.04;
            u[tangent] = 0.02 * position + 0.01 * position * position;
            u
        };
        let u: Vec<_> = (0..125)
            .map(|i| value((at(shape, i)[normal] as f64 + 0.5) * h[normal]))
            .collect();
        let ghosts = [
            value(-0.5 * h[normal]),
            value((shape[normal] as f64 + 0.5) * h[normal]),
        ];
        let phi: Vec<_> = (0..125)
            .map(|i| 3.0 + 0.01 * (at(shape, i)[normal] as f64 + 0.5))
            .collect();
        let mut flow = ReservoirDynamicSgs::new(shape, h, normal, usize::MAX).unwrap();
        let mut scalar = ReservoirDynamicScalar::new(shape, normal, usize::MAX).unwrap();
        let (mc, hc) = controls(DynamicAveraging::LocalBox);
        flow.prepare_with_scalar(
            mc,
            hc,
            &rho,
            &u,
            ghosts,
            &phi,
            [2.995, 3.055],
            &mut scalar,
            || false,
        )
        .unwrap();
        assert!(flow.viscosities().iter().all(|v| *v == 0.0));
        assert!(scalar.diffusivities().iter().all(|v| *v == 0.0));
        let scale = u
            .iter()
            .chain(ghosts.iter())
            .flat_map(|row| (0..3).map(|a| (row[a] - u[0][a]).abs()))
            .fold(0.0_f64, f64::max);
        for i in 0..125 {
            let physical = 0.02 + 0.02 * (at(shape, i)[normal] as f64 + 0.5) * h[normal];
            close(
                flow.inner.strain_rate[i] * scale / flow.filter_width(),
                physical,
            );
        }
    }
}

#[test]
fn g1_ghost_updates_refresh_current_heat_without_wrapping_or_overwriting_momentum() {
    let shape = [7, 5, 5];
    let (rho, u, phi) = field(shape);
    let ghosts = [[0.09, -0.04, 0.02], [-0.02, 0.03, -0.05]];
    let (mc, hc) = controls(DynamicAveraging::LocalBox);
    let mut flow = ReservoirDynamicSgs::new(shape, [1.0; 3], 0, usize::MAX).unwrap();
    let mut scalar = ReservoirDynamicScalar::new(shape, 0, usize::MAX).unwrap();
    flow.prepare_with_scalar(
        mc,
        hc,
        &rho,
        &u,
        ghosts,
        &phi,
        [2.7, 3.8],
        &mut scalar,
        || false,
    )
    .unwrap();
    let nu = flow.viscosities().to_vec();
    let previous = scalar.diffusivities().to_vec();
    flow.prepare_with_scalar(
        mc,
        hc,
        &rho,
        &u,
        ghosts,
        &phi,
        [4.1, 3.8],
        &mut scalar,
        || false,
    )
    .unwrap();
    assert_eq!(nu, flow.viscosities());
    assert!(
        (0..rho.len())
            .filter(|&i| at(shape, i)[0] == 0)
            .any(|i| (scalar.diffusivities()[i] - previous[i]).abs() > 1e-9)
    );
    for i in 0..rho.len() {
        if at(shape, i)[0] > 0 {
            close(scalar.diffusivities()[i], previous[i]);
        }
    }
    let mut changed = ghosts;
    changed[1] = [0.2, -0.15, 0.11];
    flow.prepare_with_scalar(
        mc,
        hc,
        &rho,
        &u,
        changed,
        &phi,
        [2.7, 3.8],
        &mut scalar,
        || false,
    )
    .unwrap();
    for i in 0..rho.len() {
        if at(shape, i)[0] + 1 < shape[0] {
            close(flow.viscosities()[i], nu[i]);
        }
    }
}

#[test]
fn g3_reservoir_covariance_preserves_boosts_scalar_origins_and_physical_units() {
    let shape = [5; 3];
    let h = [0.5, 0.75, 1.25];
    let (rho, u, phi) = field(shape);
    let ghosts = [[0.125, -0.25, 0.0625], [-0.125, 0.125, -0.0625]];
    let gs = [2.5, 3.75];
    let (mc, hc) = controls(DynamicAveraging::LocalBox);
    let mut base = ReservoirDynamicSgs::new(shape, h, 1, usize::MAX).unwrap();
    let mut heat = ReservoirDynamicScalar::new(shape, 1, usize::MAX).unwrap();
    base.prepare_with_scalar(mc, hc, &rho, &u, ghosts, &phi, gs, &mut heat, || false)
        .unwrap();
    for (length, speed, density) in [(1.0, 1.0, 1.0), (1e-60, 1e40, 1e100), (1e60, 1e-40, 1e-100)] {
        let r: Vec<_> = rho.iter().map(|v| v * density).collect();
        let v: Vec<_> = u.iter().map(|row| row.map(|v| (v + 8.0) * speed)).collect();
        let gv = ghosts.map(|row| row.map(|v| (v + 8.0) * speed));
        let p: Vec<_> = phi.iter().map(|v| 4.0 * v - 16.0).collect();
        let mut tested =
            ReservoirDynamicSgs::new(shape, h.map(|v| v * length), 1, usize::MAX).unwrap();
        let mut scalar = ReservoirDynamicScalar::new(shape, 1, usize::MAX).unwrap();
        tested
            .prepare_with_scalar(
                mc,
                hc,
                &r,
                &v,
                gv,
                &p,
                gs.map(|v| 4.0 * v - 16.0),
                &mut scalar,
                || false,
            )
            .unwrap();
        for i in 0..rho.len() {
            close(
                tested.viscosities()[i] / (length * speed),
                base.viscosities()[i],
            );
            close(
                scalar.diffusivities()[i] / (length * speed),
                heat.diffusivities()[i],
            );
        }
    }
}

#[test]
fn g0_g4_reservoir_controls_bytes_and_every_cancellation_poll_refuse_and_retry() {
    let shape = [3; 3];
    let (rho, u, phi) = field(shape);
    let saved = (rho.clone(), u.clone(), phi.clone());
    let (mc, hc) = controls(DynamicAveraging::LocalBox);
    let ghosts = [[0.03; 3], [-0.02; 3]];
    let mut flow = ReservoirDynamicSgs::new(shape, [1.0; 3], 2, 136 * 27).unwrap();
    let mut scalar = ReservoirDynamicScalar::new(shape, 2, 64 * 27).unwrap();
    assert_eq!(flow.workspace_bytes() + scalar.workspace_bytes(), 200 * 27);
    assert!(ReservoirDynamicSgs::new(shape, [1.0; 3], 2, 136 * 27 - 1).is_err());
    assert!(ReservoirDynamicSgs::new(shape, [1.0; 3], 3, usize::MAX).is_err());
    let mut polls = 0;
    let report = flow
        .prepare_with_scalar(
            mc,
            hc,
            &rho,
            &u,
            ghosts,
            &phi,
            [2.7, 3.8],
            &mut scalar,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let nu = flow.viscosities().to_vec();
    let alpha = scalar.diffusivities().to_vec();
    let pointers = (
        flow.inner.velocity.as_ptr(),
        flow.inner.viscosity.as_ptr(),
        scalar.inner.diffusivity.as_ptr(),
    );
    for stop in 1..=polls {
        let mut seen = 0;
        assert_eq!(
            flow.prepare_with_scalar(
                mc,
                hc,
                &rho,
                &u,
                ghosts,
                &phi,
                [2.7, 3.8],
                &mut scalar,
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(DynamicSgsError::Cancelled)
        );
        assert_eq!(
            flow.prepare_with_scalar(
                mc,
                hc,
                &rho,
                &u,
                ghosts,
                &phi,
                [2.7, 3.8],
                &mut scalar,
                || false
            )
            .unwrap(),
            report
        );
        assert_eq!(flow.viscosities(), nu);
        assert_eq!(scalar.diffusivities(), alpha);
        assert_eq!(
            pointers,
            (
                flow.inner.velocity.as_ptr(),
                flow.inner.viscosity.as_ptr(),
                scalar.inner.diffusivity.as_ptr()
            )
        );
    }
    assert_eq!((rho, u, phi), saved);
}

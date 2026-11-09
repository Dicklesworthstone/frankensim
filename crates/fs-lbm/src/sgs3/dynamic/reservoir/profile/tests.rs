use super::*;

fn controls(mode: DynamicAveraging) -> (DynamicSgsConfig, DynamicScalarConfig) {
    (
        DynamicSgsConfig {
            averaging: mode,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        DynamicScalarConfig {
            averaging: mode,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    )
}

fn close(a: f64, b: f64) {
    // The independent unscaled tensor contraction cancels near zero, whereas
    // production contracts normalized differences. Allow absolute roundoff
    // there while retaining a relative check on resolved coefficients.
    assert!(
        (a - b).abs() < 3e-15 + 8e-12 * b.abs(),
        "{a:.17e} != {b:.17e}"
    );
}

fn donor_box(shape: [usize; 3], normal: usize, cell: usize) -> Vec<usize> {
    let p = coordinates(shape, cell);
    let mut result = Vec::new();
    for z in -1..=1 {
        for y in -1..=1 {
            for x in -1..=1 {
                let shift = [x, y, z];
                if shift[normal] == 0 {
                    result.push(index(
                        shape,
                        core::array::from_fn(|a| {
                            (p[a] as isize + shift[a]).rem_euclid(shape[a] as isize) as usize
                        }),
                    ));
                }
            }
        }
    }
    result
}

fn fields(shape: [usize; 3]) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
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

fn exterior(shape: [usize; 3], normal: usize, side: usize, cell: usize) -> ReservoirPrimitive {
    let mut p = coordinates(shape, cell);
    p[normal] = side * (shape[normal] - 1);
    let i = index(shape, p);
    ReservoirPrimitive {
        density: 0.5 + ((i * 5 + side * 7) % 13) as f64 / 8.0,
        velocity: core::array::from_fn(|a| {
            (((i * (3 + 5 * a) + side * 19 + 7 * a) % 43) as f64 - 21.0) / 128.0
        }),
        scalar: 2.5 + ((i * 11 + side * 13) % 29) as f64 / 16.0,
    }
}

fn derivative(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    field: &[f64],
    ghosts: &[Vec<f64>; 2],
    cell: usize,
    axis: usize,
) -> f64 {
    let p = coordinates(shape, cell);
    let mut low = p;
    let mut high = p;
    low[axis] = (p[axis] + shape[axis] - 1) % shape[axis];
    high[axis] = (p[axis] + 1) % shape[axis];
    let a = if axis == normal && p[axis] == 0 {
        ghosts[0][cell]
    } else {
        field[index(shape, low)]
    };
    let b = if axis == normal && p[axis] + 1 == shape[axis] {
        ghosts[1][cell]
    } else {
        field[index(shape, high)]
    };
    (b - a) / (2.0 * h[axis])
}

fn tensors(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    u: &[[f64; 3]],
    ghosts: &[Vec<ReservoirPrimitive>; 2],
) -> Vec<(f64, [[f64; 3]; 3])> {
    let components: [Vec<f64>; 3] = core::array::from_fn(|a| u.iter().map(|v| v[a]).collect());
    let exterior: [[Vec<f64>; 2]; 3] = core::array::from_fn(|a| {
        core::array::from_fn(|s| ghosts[s].iter().map(|v| v.velocity[a]).collect())
    });
    (0..u.len())
        .map(|cell| {
            let g: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    derivative(shape, h, normal, &components[i], &exterior[i], cell, j)
                })
            });
            let mut s: [[f64; 3]; 3] =
                core::array::from_fn(|i| core::array::from_fn(|j| 0.5 * (g[i][j] + g[j][i])));
            let rate = (2.0 * s.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt();
            let trace = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
            for i in 0..3 {
                s[i][i] -= trace;
            }
            (rate, s)
        })
        .collect()
}

// Independent physical-unit, unpacked tensor/vector identities. Each exterior
// test sample is filtered independently of the interior field. The deliberately
// wrong branch supplies an unfiltered exterior at the test level.
fn reference(
    shape: [usize; 3],
    h: [f64; 3],
    normal: usize,
    mode: DynamicAveraging,
    rho: &[f64],
    u: &[[f64; 3]],
    phi: &[f64],
    filter_ghosts: bool,
) -> Vec<[f64; 4]> {
    let n = rho.len();
    let delta2 = h.iter().product::<f64>().cbrt().powi(2);
    let ratio = 9.0_f64.powf(2.0 / 3.0);
    let boxes: Vec<_> = (0..n).map(|i| donor_box(shape, normal, i)).collect();
    let r: Vec<f64> = boxes
        .iter()
        .map(|d| d.iter().map(|&j| rho[j]).sum::<f64>() / 9.0)
        .collect();
    let v: Vec<[f64; 3]> = boxes
        .iter()
        .enumerate()
        .map(|(i, d)| {
            core::array::from_fn(|a| {
                d.iter().map(|&j| rho[j] * u[j][a]).sum::<f64>() / (9.0 * r[i])
            })
        })
        .collect();
    let p: Vec<f64> = boxes
        .iter()
        .enumerate()
        .map(|(i, d)| d.iter().map(|&j| rho[j] * phi[j]).sum::<f64>() / (9.0 * r[i]))
        .collect();
    let raw: [Vec<_>; 2] =
        core::array::from_fn(|s| (0..n).map(|i| exterior(shape, normal, s, i)).collect());
    let filtered: [Vec<_>; 2] = core::array::from_fn(|s| {
        (0..n)
            .map(|i| {
                if !filter_ghosts {
                    return raw[s][i];
                }
                let mass = boxes[i].iter().map(|&j| raw[s][j].density).sum::<f64>();
                ReservoirPrimitive {
                    density: mass / 9.0,
                    velocity: core::array::from_fn(|a| {
                        boxes[i]
                            .iter()
                            .map(|&j| raw[s][j].density * raw[s][j].velocity[a])
                            .sum::<f64>()
                            / mass
                    }),
                    scalar: boxes[i]
                        .iter()
                        .map(|&j| raw[s][j].density * raw[s][j].scalar)
                        .sum::<f64>()
                        / mass,
                }
            })
            .collect()
    });
    let scalar_raw: [Vec<_>; 2] =
        core::array::from_fn(|s| raw[s].iter().map(|g| g.scalar).collect());
    let scalar_filtered: [Vec<_>; 2] =
        core::array::from_fn(|s| filtered[s].iter().map(|g| g.scalar).collect());
    let strain = tensors(shape, h, normal, u, &raw);
    let test_strain = tensors(shape, h, normal, &v, &filtered);
    let grad: Vec<[f64; 3]> = (0..n)
        .map(|i| core::array::from_fn(|a| derivative(shape, h, normal, phi, &scalar_raw, i, a)))
        .collect();
    let contractions: Vec<[f64; 4]> = (0..n)
        .map(|i| {
            let mut l = [[0.0; 3]; 3];
            let mut b = [[0.0; 3]; 3];
            let mut q = [0.0; 3];
            let mut heat = [0.0; 3];
            for &j in &boxes[i] {
                let weight = rho[j] / 9.0;
                for a in 0..3 {
                    q[a] += weight * (u[j][a] - v[i][a]) * (phi[j] - p[i]);
                    heat[a] += weight * strain[j].0 * grad[j][a];
                    for c in 0..3 {
                        l[a][c] += weight * (u[j][a] - v[i][a]) * (u[j][c] - v[i][c]);
                        b[a][c] += weight * strain[j].0 * strain[j].1[a][c];
                    }
                }
            }
            let trace = (l[0][0] + l[1][1] + l[2][2]) / 3.0;
            let mut pair = [0.0; 4];
            for a in 0..3 {
                let m = delta2
                    * (heat[a]
                        - ratio
                            * r[i]
                            * test_strain[i].0
                            * derivative(shape, h, normal, &p, &scalar_filtered, i, a));
                pair[2] += q[a] * m;
                pair[3] += m * m;
                for c in 0..3 {
                    let m = 2.0
                        * delta2
                        * (b[a][c] - ratio * r[i] * test_strain[i].0 * test_strain[i].1[a][c]);
                    pair[0] += (l[a][c] - if a == c { trace } else { 0.0 }) * m;
                    pair[1] += m * m;
                }
            }
            pair
        })
        .collect();
    (0..n)
        .map(|i| {
            let selection: Vec<usize> = if mode == DynamicAveraging::Domain {
                (0..n)
                    .filter(|&j| coordinates(shape, j)[normal] == coordinates(shape, i)[normal])
                    .collect()
            } else {
                boxes[i].clone()
            };
            let pair = selection.iter().fold([0.0; 4], |mut sum, &j| {
                for a in 0..4 {
                    sum[a] += contractions[j][a];
                }
                sum
            });
            let c = if pair[1] == 0.0 {
                0.0
            } else {
                (pair[0] / pair[1]).max(0.0)
            };
            let heat = if pair[3] == 0.0 {
                0.0
            } else {
                (pair[2] / pair[3]).max(0.0)
            };
            [
                c,
                heat,
                c * delta2 * strain[i].0,
                heat * delta2 * strain[i].0,
            ]
        })
        .collect()
}

#[test]
fn g1_spatial_ghosts_match_independent_favre_tensors_and_reject_unfiltered_test_ghosts() {
    let shape = [5, 6, 7];
    let h = [0.5, 0.75, 1.25];
    let (rho, u, phi) = fields(shape);
    let mut observable = false;
    for normal in 0..3 {
        for mode in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
            let (mc, hc) = controls(mode);
            let expected = reference(shape, h, normal, mode, &rho, &u, &phi, true);
            let wrong = reference(shape, h, normal, mode, &rho, &u, &phi, false);
            let mut flow = ReservoirDynamicSgs::new(shape, h, normal, usize::MAX).unwrap();
            let mut heat = ReservoirDynamicScalar::new(shape, normal, usize::MAX).unwrap();
            flow.prepare_with_profile(
                mc,
                hc,
                &rho,
                &u,
                &phi,
                |s, i| exterior(shape, normal, s, i),
                &mut heat,
                || false,
            )
            .unwrap();
            for i in 0..rho.len() {
                close(flow.coefficient_squared(i).unwrap(), expected[i][0]);
                close(heat.coefficient(i).unwrap(), expected[i][1]);
                close(flow.viscosities()[i], expected[i][2]);
                close(heat.diffusivities()[i], expected[i][3]);
                observable |= (wrong[i][2] - expected[i][2]).abs() > 1e-6
                    || (wrong[i][3] - expected[i][3]).abs() > 1e-6;
            }
        }
    }
    assert!(
        observable,
        "the fixture must distinguish physical from test-filtered ghosts"
    );
}

#[test]
fn g3_profiles_preserve_uniform_limits_and_unit_covariance() {
    let shape = [5; 3];
    let h = [0.5, 0.75, 1.25];
    let (rho, u, phi) = fields(shape);
    let (mc, hc) = controls(DynamicAveraging::LocalBox);
    let mut base = ReservoirDynamicSgs::new(shape, h, 1, usize::MAX).unwrap();
    let mut base_heat = ReservoirDynamicScalar::new(shape, 1, usize::MAX).unwrap();
    base.prepare_with_profile(
        mc,
        hc,
        &rho,
        &u,
        &phi,
        |s, i| exterior(shape, 1, s, i),
        &mut base_heat,
        || false,
    )
    .unwrap();
    for (length, speed, mass) in [(1.0, 1.0, 1.0), (1e-60, 1e40, 1e100), (1e60, 1e-40, 1e-100)] {
        let r: Vec<_> = rho.iter().map(|r| r * mass).collect();
        let v: Vec<_> = u.iter().map(|v| v.map(|v| (v + 8.0) * speed)).collect();
        let p: Vec<_> = phi.iter().map(|p| 4.0 * p - 16.0).collect();
        let mut flow =
            ReservoirDynamicSgs::new(shape, h.map(|h| h * length), 1, usize::MAX).unwrap();
        let mut heat = ReservoirDynamicScalar::new(shape, 1, usize::MAX).unwrap();
        flow.prepare_with_profile(
            mc,
            hc,
            &r,
            &v,
            &p,
            |s, i| {
                let g = exterior(shape, 1, s, i);
                ReservoirPrimitive {
                    density: g.density * mass,
                    velocity: g.velocity.map(|v| (v + 8.0) * speed),
                    scalar: 4.0 * g.scalar - 16.0,
                }
            },
            &mut heat,
            || false,
        )
        .unwrap();
        for i in 0..rho.len() {
            close(
                flow.viscosities()[i] / (length * speed),
                base.viscosities()[i],
            );
            close(
                heat.diffusivities()[i] / (length * speed),
                base_heat.diffusivities()[i],
            );
        }
    }
    let ghosts = [[0.125, -0.25, 0.0625], [-0.125, 0.125, -0.0625]];
    let scalar = [2.5, 3.75];
    base.prepare_with_scalar(
        mc,
        hc,
        &rho,
        &u,
        ghosts,
        &phi,
        scalar,
        &mut base_heat,
        || false,
    )
    .unwrap();
    let mut profile = ReservoirDynamicSgs::new(shape, h, 1, usize::MAX).unwrap();
    let mut heat = ReservoirDynamicScalar::new(shape, 1, usize::MAX).unwrap();
    profile
        .prepare_with_profile(
            mc,
            hc,
            &rho,
            &u,
            &phi,
            |s, i| ReservoirPrimitive {
                density: exterior(shape, 1, s, i).density,
                velocity: ghosts[s],
                scalar: scalar[s],
            },
            &mut heat,
            || false,
        )
        .unwrap();
    for i in 0..rho.len() {
        close(profile.viscosities()[i], base.viscosities()[i]);
        close(heat.diffusivities()[i], base_heat.diffusivities()[i]);
    }
}

#[test]
fn g0_g4_profile_admission_and_every_cancellation_poll_preserve_sources_and_retry() {
    let shape = [3; 3];
    let (rho, u, phi) = fields(shape);
    let saved = (rho.clone(), u.clone(), phi.clone());
    let (mc, hc) = controls(DynamicAveraging::LocalBox);
    let mut flow = ReservoirDynamicSgs::new(shape, [1.0; 3], 2, 136 * 27).unwrap();
    let mut heat = ReservoirDynamicScalar::new(shape, 2, 64 * 27).unwrap();
    let mut polls = 0;
    let expected = flow
        .prepare_with_profile(
            mc,
            hc,
            &rho,
            &u,
            &phi,
            |s, i| exterior(shape, 2, s, i),
            &mut heat,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let nu = flow.viscosities().to_vec();
    let alpha = heat.diffusivities().to_vec();
    let address = (
        flow.inner.velocity.as_ptr(),
        heat.inner.diffusivity.as_ptr(),
    );
    for stop in 1..=polls {
        let mut seen = 0;
        assert!(matches!(
            flow.prepare_with_profile(
                mc,
                hc,
                &rho,
                &u,
                &phi,
                |s, i| exterior(shape, 2, s, i),
                &mut heat,
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(DynamicSgsError::Cancelled)
        ));
        assert_eq!(
            flow.prepare_with_profile(
                mc,
                hc,
                &rho,
                &u,
                &phi,
                |s, i| exterior(shape, 2, s, i),
                &mut heat,
                || false
            )
            .unwrap(),
            expected
        );
        assert_eq!(flow.viscosities(), nu);
        assert_eq!(heat.diffusivities(), alpha);
    }
    assert_eq!(
        address,
        (
            flow.inner.velocity.as_ptr(),
            heat.inner.diffusivity.as_ptr()
        )
    );
    assert_eq!((rho, u, phi), saved);
    let zero = vec![[0.0; 3]; 27];
    let constant = vec![3.0; 27];
    let density = vec![1.0; 27];
    for invalid in 0..3 {
        assert!(
            flow.prepare_with_profile(
                mc,
                hc,
                &density,
                &zero,
                &constant,
                |s, i| {
                    let mut g = ReservoirPrimitive {
                        density: 1.0,
                        velocity: [0.0; 3],
                        scalar: 3.0,
                    };
                    if s == 1 && i == 26 {
                        match invalid {
                            0 => g.density = 0.0,
                            1 => g.velocity[1] = f64::NAN,
                            _ => g.scalar = f64::INFINITY,
                        }
                    }
                    g
                },
                &mut heat,
                || false
            )
            .is_err()
        );
    }
}

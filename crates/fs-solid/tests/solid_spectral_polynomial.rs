//! Exact cubic force histories, independently integrated moments and work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{CompleteSpectralSolid, CompleteSpectralSolidControls};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 39,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 128,
        maximum_preparation_bytes: 16 * 1024 * 1024,
        maximum_preparation_multiply_terms: 30_000_000,
        maximum_step_multiply_terms: 2_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 1000,
    }
}
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    let (mut stiffness, mut mass) = (Coo::new(n, n), Coo::new(n, n));
    for i in 0..n {
        for j in 0..n {
            stiffness.push(i, j, k[i * n + j]);
            mass.push(i, j, m[i * n + j]);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t,
        "actual={a:e}, expected={b:e}, allowance={t:e}"
    );
}
fn bits(s: &CompleteSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([
            s.time_s(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
            s.cumulative_work_j(),
        ])
        .map(f64::to_bits)
        .collect()
}

// Independent continuous solution in s=t/h. For well-separated frequencies,
// solve a cubic particular solution plus sine/cosine and integrate its terms.
// Small/zero frequencies instead use the original ODE coefficient recurrence.
fn oracle(lambda: f64, h: f64, q0: f64, v0: f64, g: [f64; 4]) -> (f64, f64, [f64; 3], f64, f64) {
    let z = lambda.sqrt() * h;
    if z < 0.5 {
        let mut a = [0.0; 48];
        a[0] = q0;
        a[1] = h * v0;
        for k in 0..46 {
            let force = if k < 4 { h * h * g[k] } else { 0.0 };
            a[k + 2] = (force - z * z * a[k]) / ((k + 1) * (k + 2)) as f64;
        }
        let q = a.iter().sum();
        let v = (1..48).map(|k| k as f64 * a[k] / h).sum();
        let moments = core::array::from_fn(|r| (0..48).map(|k| a[k] / (k + r + 1) as f64).sum());
        let work = (0..4)
            .map(|r| {
                g[r] * (1..48)
                    .map(|k| k as f64 * a[k] / (r + k) as f64)
                    .sum::<f64>()
            })
            .sum();
        let acceleration = g.iter().sum::<f64>() - lambda * q;
        return (q, v, moments, work, acceleration);
    }
    let z2 = z * z;
    let p = [
        g[0] / lambda - 2.0 * g[2] / (lambda * z2),
        g[1] / lambda - 6.0 * g[3] / (lambda * z2),
        g[2] / lambda,
        g[3] / lambda,
    ];
    let a = q0 - p[0];
    let b = (h * v0 - p[1]) / z;
    let (s, c) = z.sin_cos();
    let q = p.iter().sum::<f64>() + a * c + b * s;
    let v = (p[1] + 2.0 * p[2] + 3.0 * p[3] - z * a * s + z * b * c) / h;
    // Integrals of s^r cos(zs), s^r sin(zs), r=0..3, by parts.
    let mut ci = [0.0; 4];
    let mut si = [0.0; 4];
    ci[0] = s / z;
    si[0] = (1.0 - c) / z;
    for r in 1..4 {
        ci[r] = s / z - r as f64 * si[r - 1] / z;
        si[r] = -c / z + r as f64 * ci[r - 1] / z;
    }
    let moments = core::array::from_fn(|k| {
        (0..4).map(|r| p[r] / (r + k + 1) as f64).sum::<f64>() + a * ci[k] + b * si[k]
    });
    let work = (0..4)
        .map(|r| {
            g[r] * ((1..4)
                .map(|j| j as f64 * p[j] / (r + j) as f64)
                .sum::<f64>()
                - z * a * si[r]
                + z * b * ci[r])
        })
        .sum();
    (q, v, moments, work, g.iter().sum::<f64>() - lambda * q)
}

#[test]
fn g1_cubic_moments_and_work_match_continuous_motion_across_phase_regimes() {
    with_cx(false, |cx| {
        let h = 0.17;
        let g = [0.7, -1.3, 2.1, -0.9];
        let coefficients = g.map(|v| [v]);
        let source = core::array::from_fn(|i| coefficients[i].as_slice());
        let mut largest_work_error = 0.0_f64;
        for z in [
            0.0_f64, 1e-10, 0.499, 0.501, 3.999999, 4.0, 4.000001, 12.3, 37.8,
        ] {
            let lambda = (z / h).powi(2);
            let pencil = assembly(&[lambda], &[1.0], 1);
            let mut solid = CompleteSpectralSolid::prepare_cubic_forces_with_basis(
                &pencil,
                &[lambda],
                &[1.0],
                h,
                &[0.2],
                &[-0.13],
                controls(),
                cx,
            )
            .unwrap();
            let result = solid.advance_cubic(source, cx).unwrap();
            let (q, v, moments, work, acceleration) = oracle(lambda, h, 0.2, -0.13, g);
            near(solid.displacement_m()[0], q, 3e-14);
            near(solid.velocity_mps()[0], v, 3e-12);
            near(result.acceleration_mps2[0], acceleration, 5e-10);
            for k in 0..3 {
                near(result.displacement_moments_m[k][0], moments[k], 5e-14);
            }
            largest_work_error =
                largest_work_error.max((result.motion.energy.supplied_work_j - work).abs());
            near(result.motion.energy.supplied_work_j, work, 2e-13);
            near(result.motion.energy.balance_error_j, 0.0, 2e-10);
        }
        println!("cubic exact phase ladder: max independent work error={largest_work_error:e} J");
    });
}

#[test]
fn g1_complete_coupled_coordinates_retain_cubic_forcing_and_original_work() {
    with_cx(false, |cx| {
        let n = 96;
        let s = [[1.0, 0.5], [0.25, 1.0]];
        let inv = [[8.0 / 7.0, -4.0 / 7.0], [-2.0 / 7.0, 8.0 / 7.0]];
        let lambda: Vec<_> = (0..n).map(|i| (7.0 + 0.4 * i as f64).powi(2)).collect();
        let (mut k, mut m, mut basis) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]);
        for first in (0..n).step_by(2) {
            for a in 0..2 {
                for b in 0..2 {
                    basis[(first + a) * n + first + b] = s[a][b];
                    for j in 0..2 {
                        k[(first + a) * n + first + b] += inv[j][a] * lambda[first + j] * inv[j][b];
                        m[(first + a) * n + first + b] += inv[j][a] * inv[j][b];
                    }
                }
            }
        }
        let q0: Vec<_> = (0..n).map(|j| 0.001 + 0.0001 * (j % 5) as f64).collect();
        let v0: Vec<_> = (0..n).map(|j| -0.01 + 0.001 * (j % 7) as f64).collect();
        let g: [Vec<f64>; 4] = core::array::from_fn(|r| {
            (0..n)
                .map(|j| 0.01 * (1 + r) as f64 - 0.001 * (j % 11) as f64)
                .collect()
        });
        let transform = |x: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| (0..n).map(|j| basis[i * n + j] * x[j]).sum())
                .collect()
        };
        let f: [Vec<f64>; 4] = core::array::from_fn(|r| {
            let mut f = vec![0.0; n];
            for first in (0..n).step_by(2) {
                for a in 0..2 {
                    f[first + a] = (0..2).map(|j| inv[j][a] * g[r][first + j]).sum();
                }
            }
            f
        });
        let pencil = assembly(&k, &m, n);
        let h = 0.11;
        let mut solid = CompleteSpectralSolid::prepare_cubic_forces_with_basis(
            &pencil,
            &lambda,
            &basis,
            h,
            &transform(&q0),
            &transform(&v0),
            controls(),
            cx,
        )
        .unwrap();
        let result = solid
            .advance_cubic(core::array::from_fn(|r| f[r].as_slice()), cx)
            .unwrap();
        let mut q = vec![0.0; n];
        let mut v = vec![0.0; n];
        let mut work = 0.0;
        let mut moments = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
        for i in 0..n {
            let point = oracle(
                lambda[i],
                h,
                q0[i],
                v0[i],
                core::array::from_fn(|r| g[r][i]),
            );
            q[i] = point.0;
            v[i] = point.1;
            work += point.3;
            for j in 0..3 {
                moments[j][i] = point.2[j];
            }
        }
        for (a, b) in solid
            .displacement_m()
            .iter()
            .zip(transform(&q))
            .chain(solid.velocity_mps().iter().zip(transform(&v)))
        {
            near(*a, b, 1e-13);
        }
        for r in 0..3 {
            for (a, b) in result.displacement_moments_m[r]
                .iter()
                .zip(transform(&moments[r]))
            {
                near(*a, b, 2e-14);
            }
        }
        near(result.motion.energy.supplied_work_j, work, 2e-14);
        near(result.motion.energy.balance_error_j, 0.0, 2e-13);
        assert_eq!(solid.free_dofs().len(), n);
    });
}

#[test]
fn g1_cubic_partition_preserves_all_three_continuous_moments() {
    with_cx(false, |cx| {
        let (lambda, h) = (49.0, 0.3);
        let pencil = assembly(&[lambda], &[1.0], 1);
        let prepare = |dt| {
            CompleteSpectralSolid::prepare_cubic_forces_with_basis(
                &pencil,
                &[lambda],
                &[1.0],
                dt,
                &[0.17],
                &[-0.2],
                controls(),
                cx,
            )
            .unwrap()
        };
        let (mut coarse, mut fine) = (prepare(h), prepare(h / 8.0));
        let f = [0.4, -0.7, 0.9, -0.3];
        let arrays = f.map(|v| [v]);
        let whole = coarse
            .advance_cubic(core::array::from_fn(|r| arrays[r].as_slice()), cx)
            .unwrap();
        let choose = [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [1.0, 2.0, 1.0, 0.0],
            [1.0, 3.0, 3.0, 1.0],
        ];
        let (mut moments, mut work) = ([0.0; 3], 0.0);
        let b = 1.0 / 8.0_f64;
        for part in 0..8 {
            let a = part as f64 * b;
            let local: [f64; 4] = core::array::from_fn(|k| {
                (k..4)
                    .map(|r| choose[r][k] * f[r] * a.powi((r - k) as i32) * b.powi(k as i32))
                    .sum()
            });
            let arrays = local.map(|v| [v]);
            let step = fine
                .advance_cubic(core::array::from_fn(|r| arrays[r].as_slice()), cx)
                .unwrap();
            work += step.motion.energy.supplied_work_j;
            for k in 0..3 {
                moments[k] += b
                    * (0..=k)
                        .map(|j| {
                            choose[k][j]
                                * a.powi((k - j) as i32)
                                * b.powi(j as i32)
                                * step.displacement_moments_m[j][0]
                        })
                        .sum::<f64>();
            }
        }
        near(coarse.displacement_m()[0], fine.displacement_m()[0], 2e-14);
        near(coarse.velocity_mps()[0], fine.velocity_mps()[0], 2e-13);
        for k in 0..3 {
            near(whole.displacement_moments_m[k][0], moments[k], 2e-14);
        }
        near(whole.motion.energy.supplied_work_j, work, 2e-14);
        // Degree-one sources remain equivalent to the existing affine path.
        let (mut cubic, mut affine) = (prepare(h), prepare(h));
        let a = cubic
            .advance_cubic([&[0.4], &[-0.7], &[0.0], &[0.0]], cx)
            .unwrap();
        let b = affine.advance_affine(&[0.4], &[-0.3], cx).unwrap();
        near(cubic.displacement_m()[0], affine.displacement_m()[0], 2e-14);
        near(cubic.velocity_mps()[0], affine.velocity_mps()[0], 2e-13);
        near(
            a.motion.energy.supplied_work_j,
            b.energy.supplied_work_j,
            2e-14,
        );
    });
}

#[test]
fn g4_cubic_capacity_source_energy_and_cancellation_refusals_preserve_state() {
    with_cx(false, |cx| {
        let pencil = assembly(&[1.0], &[1.0], 1);
        let mut cap = controls();
        cap.maximum_free_dofs = 1;
        cap.maximum_preparation_bytes = 2112;
        cap.maximum_preparation_multiply_terms = 4132;
        cap.maximum_step_multiply_terms = 1048;
        let make = |c| {
            CompleteSpectralSolid::prepare_cubic_forces_with_basis(
                &pencil,
                &[1.0],
                &[1.0],
                0.1,
                &[0.2],
                &[-0.1],
                c,
                cx,
            )
        };
        let mut solid = make(cap).unwrap();
        for field in 0..4 {
            let mut bad = cap;
            match field {
                0 => bad.maximum_free_dofs -= 1,
                1 => bad.maximum_preparation_bytes -= 1,
                2 => bad.maximum_preparation_multiply_terms -= 1,
                _ => bad.maximum_step_multiply_terms -= 1,
            }
            assert!(make(bad).is_err());
        }
        let accepted = bits(&solid);
        assert!(
            solid
                .advance_cubic([&[], &[0.0], &[0.0], &[0.0]], cx)
                .is_err()
        );
        assert_eq!(bits(&solid), accepted);
        assert!(
            solid
                .advance_cubic([&[f64::NAN], &[0.0], &[0.0], &[0.0]], cx)
                .is_err()
        );
        assert_eq!(bits(&solid), accepted);
        assert!(
            solid
                .advance_cubic([&[f64::MAX], &[f64::MAX], &[f64::MAX], &[f64::MAX]], cx)
                .is_err()
        );
        assert_eq!(bits(&solid), accepted);
        with_cx(true, |cancelled| {
            assert!(
                solid
                    .advance_cubic([&[0.1], &[0.2], &[0.3], &[0.4]], cancelled)
                    .is_err()
            )
        });
        assert_eq!(bits(&solid), accepted);
        let mut retry = make(cap).unwrap();
        solid
            .advance_cubic([&[0.1], &[0.2], &[0.3], &[0.4]], cx)
            .unwrap();
        retry
            .advance_cubic([&[0.1], &[0.2], &[0.3], &[0.4]], cx)
            .unwrap();
        assert_eq!(bits(&solid), bits(&retry));
        let mut standard = CompleteSpectralSolid::prepare_with_basis(
            &pencil,
            &[1.0],
            &[1.0],
            0.1,
            &[0.2],
            &[-0.1],
            cap,
            cx,
        )
        .unwrap();
        assert!(standard.advance_cubic([&[0.0]; 4], cx).is_err());
        let mut loose = controls();
        loose.maximum_relative_operator_error = 0.01;
        loose.maximum_work_error_j = 1e-14;
        let mut inaccurate = CompleteSpectralSolid::prepare_cubic_forces_with_basis(
            &pencil,
            &[1.001],
            &[1.0],
            0.1,
            &[0.2],
            &[-0.1],
            loose,
            cx,
        )
        .unwrap();
        let before = bits(&inaccurate);
        assert!(inaccurate.advance_cubic([&[0.0]; 4], cx).is_err());
        assert_eq!(bits(&inaccurate), before);
    });
}

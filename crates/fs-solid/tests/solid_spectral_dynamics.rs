//! Complete spectral propagation beyond the augmented-exponential dimension cap.
//! Independent continuous solutions check all coordinates and applied work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolid, CompleteSpectralSolidControls, CompleteSpectralSolidError,
};
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
                kernel_id: 35,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 128,
        maximum_preparation_bytes: 8 * 1024 * 1024,
        maximum_preparation_multiply_terms: 20_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-10,
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
        free_dofs: (0..n).map(|i| 2 * i + 1).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "actual={a:e}, expected={b:e}, allowance={allowance:e}"
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
fn multiply(a: &[f64], q: &[f64]) -> Vec<f64> {
    let n = q.len();
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * q[j]).sum())
        .collect()
}

// Independent continuous solution q''+w²q=g0+rate*t, including integral q dt.
fn ramp(w: f64, q0: f64, v0: f64, g0: f64, rate: f64, t: f64) -> (f64, f64, f64) {
    let (sin, cos) = (w * t).sin_cos();
    let (a, b) = (g0 / (w * w), rate / (w * w));
    let (c, d) = (q0 - a, (v0 - b) / w);
    (
        a + b * t + c * cos + d * sin,
        b - w * c * sin + w * d * cos,
        a * t + 0.5 * b * t * t + c * sin / w + d * (1.0 - cos) / w,
    )
}

#[test]
fn g1_all_128_coupled_coordinates_match_continuous_ramp_and_work() {
    with_cx(false, |cx| {
        let n = 128;
        // Each two-coordinate block is physically coupled in BOTH M and K.
        // Every one of the 128 modal coordinates is excited and retained.
        let s = [[1.0, 0.5], [0.25, 1.0]];
        let inv = [[8.0 / 7.0, -4.0 / 7.0], [-2.0 / 7.0, 8.0 / 7.0]];
        let omega: Vec<_> = (0..n).map(|i| 10.0 + 0.3 * i as f64).collect();
        let lambda: Vec<_> = omega.iter().map(|w| w * w).collect();
        let (mut m, mut k, mut basis) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]);
        for first in (0..n).step_by(2) {
            for a in 0..2 {
                for b in 0..2 {
                    basis[(first + a) * n + first + b] = s[a][b];
                    for j in 0..2 {
                        m[(first + a) * n + first + b] += inv[j][a] * inv[j][b];
                        k[(first + a) * n + first + b] += inv[j][a] * lambda[first + j] * inv[j][b];
                    }
                }
            }
        }
        let pencil = assembly(&k, &m, n);
        let q0: Vec<_> = (0..n).map(|i| 1e-3 * (1.0 + (i % 7) as f64)).collect();
        let v0: Vec<_> = (0..n).map(|i| -0.01 + 0.001 * (i % 11) as f64).collect();
        let g0: Vec<_> = (0..n).map(|i| 0.03 - 0.001 * (i % 13) as f64).collect();
        let rate: Vec<_> = (0..n).map(|i| 0.001 * ((i % 9) as f64 - 4.0)).collect();
        let force = |t: f64| -> Vec<f64> {
            let mut f = vec![0.0; n];
            for first in (0..n).step_by(2) {
                for a in 0..2 {
                    f[first + a] = (0..2)
                        .map(|j| inv[j][a] * (g0[first + j] + rate[first + j] * t))
                        .sum();
                }
            }
            f
        };
        let dt = 0.017;
        let mut solid = CompleteSpectralSolid::prepare_with_basis(
            &pencil,
            &lambda,
            &basis,
            dt,
            &multiply(&basis, &q0),
            &multiply(&basis, &v0),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(solid.free_dofs(), pencil.free_dofs);
        let mut maximum_error = 0.0_f64;
        let mut maximum_work_error = 0.0_f64;
        for interval in 0..12 {
            let (a, b) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let step = solid.advance_affine(&force(a), &force(b), cx).unwrap();
            let (mut q, mut v, mut work) = (vec![0.0; n], vec![0.0; n], 0.0);
            for j in 0..n {
                let (qa, _, ia) = ramp(omega[j], q0[j], v0[j], g0[j], rate[j], a);
                let (qb, vb, ib) = ramp(omega[j], q0[j], v0[j], g0[j], rate[j], b);
                q[j] = qb;
                v[j] = vb;
                work +=
                    (g0[j] + rate[j] * b) * qb - (g0[j] + rate[j] * a) * qa - rate[j] * (ib - ia);
            }
            for (actual, expected) in solid
                .displacement_m()
                .iter()
                .zip(multiply(&basis, &q))
                .chain(solid.velocity_mps().iter().zip(multiply(&basis, &v)))
            {
                maximum_error = maximum_error.max((actual - expected).abs());
                near(*actual, expected, 3e-13);
            }
            maximum_work_error = maximum_work_error.max((step.energy.supplied_work_j - work).abs());
            near(step.energy.supplied_work_j, work, 2e-14);
            near(step.energy.balance_error_j, 0.0, 3e-13);
            assert_eq!(step.energy.dissipated_j, 0.0);
            assert_eq!(step.energy.cumulative_dissipated_j, 0.0);
        }
        println!(
            "complete spectral n={n}: max trajectory error={maximum_error:e}, max independent work error={maximum_work_error:e}"
        );
    });
}

// Independent power-series solution from the original second-order ODE.
fn series(lambda: f64, q0: f64, v0: f64, g0: f64, rate: f64, t: f64) -> (f64, f64, f64) {
    let mut a = [0.0; 20];
    a[0] = q0;
    a[1] = v0;
    for k in 0..18 {
        let forcing = if k == 0 {
            g0
        } else if k == 1 {
            rate
        } else {
            0.0
        };
        a[k + 2] = (forcing - lambda * a[k]) / ((k + 1) * (k + 2)) as f64;
    }
    let q = (0..20).map(|k| a[k] * t.powi(k as i32)).sum();
    let v = (1..20)
        .map(|k| k as f64 * a[k] * t.powi(k as i32 - 1))
        .sum();
    let integral = (0..20)
        .map(|k| a[k] * t.powi(k as i32 + 1) / (k + 1) as f64)
        .sum();
    (q, v, integral)
}

#[test]
fn g1_rigid_soft_modes_partitioned_ramps_and_force_jumps_are_physical() {
    with_cx(false, |cx| {
        let lambda = [0.0, 1e-20, 1e-4];
        let masses: [f64; 3] = [2.0, 3.0, 5.0];
        let (mut m, mut k, mut basis) = ([0.0; 9], [0.0; 9], [0.0; 9]);
        for i in 0..3 {
            m[4 * i] = masses[i];
            k[4 * i] = masses[i] * lambda[i];
            basis[4 * i] = 1.0 / masses[i].sqrt();
        }
        let pencil = assembly(&k, &m, 3);
        let (q0, v0, g0, rate) = (
            [0.02, 0.03, -0.01],
            [0.1, -0.2, 0.3],
            [0.5, -0.3, 0.4],
            [0.2, 0.4, -0.1],
        );
        let force = |t: f64| -> Vec<f64> {
            (0..3)
                .map(|i| masses[i].sqrt() * (g0[i] + rate[i] * t))
                .collect()
        };
        let prepare = |h: f64| {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &lambda,
                &basis,
                h,
                &multiply(&basis, &q0),
                &multiply(&basis, &v0),
                controls(),
                cx,
            )
            .unwrap()
        };
        let mut whole = prepare(0.5);
        let mut pieces = prepare(0.0625);
        let step = whole.advance_affine(&force(0.0), &force(0.5), cx).unwrap();
        for i in 0..8 {
            pieces
                .advance_affine(
                    &force(i as f64 * 0.0625),
                    &force((i + 1) as f64 * 0.0625),
                    cx,
                )
                .unwrap();
        }
        let mut work = 0.0;
        for i in 0..3 {
            let (q, v, integral) = series(lambda[i], q0[i], v0[i], g0[i], rate[i], 0.5);
            near(whole.displacement_m()[i], q / masses[i].sqrt(), 3e-15);
            near(whole.velocity_mps()[i], v / masses[i].sqrt(), 3e-15);
            near(whole.displacement_m()[i], pieces.displacement_m()[i], 3e-15);
            near(whole.velocity_mps()[i], pieces.velocity_mps()[i], 3e-15);
            work += (g0[i] + 0.5 * rate[i]) * q - g0[i] * q0[i] - rate[i] * integral;
        }
        near(step.energy.supplied_work_j, work, 3e-15);
        near(whole.cumulative_work_j(), pieces.cumulative_work_j(), 3e-15);
        let (rigid_u, rigid_v) = (whole.displacement_m()[0], whole.velocity_mps()[0]);
        let release = whole.advance(&[0.0; 3], cx).unwrap();
        near(whole.velocity_mps()[0], rigid_v, 0.0);
        near(whole.displacement_m()[0], rigid_u + 0.5 * rigid_v, 3e-16);
        assert_eq!(release.energy.supplied_work_j, 0.0);
    });
}

#[test]
fn g0_incomplete_wrong_and_singular_bases_refuse_without_a_unit_scale_floor() {
    with_cx(false, |cx| {
        let pencil = assembly(&[1.0, 0.0, 0.0, 2.0], &[1.0, 0.0, 0.0, 1.0], 2);
        let prepare = |values: &[f64], basis: &[f64]| {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                values,
                basis,
                0.1,
                &[0.0; 2],
                &[0.0; 2],
                controls(),
                cx,
            )
        };
        assert!(prepare(&[1.0], &[1.0, 0.0]).is_err());
        assert!(prepare(&[1.0, 1.0], &[1.0, 1.0, 0.0, 0.0]).is_err());
        assert!(prepare(&[1.0, 3.0], &[1.0, 0.0, 0.0, 1.0]).is_err());
        assert!(prepare(&[1.0, -2.0], &[1.0, 0.0, 0.0, 1.0]).is_err());
        let tiny = assembly(&[1e-30], &[1.0], 1);
        assert!(matches!(
            CompleteSpectralSolid::prepare_with_basis(
                &tiny,
                &[2e-30],
                &[1.0],
                0.1,
                &[0.0],
                &[0.0],
                controls(),
                cx
            ),
            Err(CompleteSpectralSolidError::Operator { .. })
        ));
        let singular = assembly(&[0.0], &[0.0], 1);
        assert!(
            CompleteSpectralSolid::prepare_with_basis(
                &singular,
                &[0.0],
                &[1.0],
                0.1,
                &[0.0],
                &[0.0],
                controls(),
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn g4_capacity_energy_and_cancellation_refusals_preserve_accepted_state() {
    with_cx(false, |cx| {
        let pencil = assembly(&[1.0], &[1.0], 1);
        let prepare = |c| {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[1.0],
                &[1.0],
                0.1,
                &[0.1],
                &[0.2],
                c,
                cx,
            )
        };
        let mut cap = controls();
        cap.maximum_free_dofs = 1;
        cap.maximum_preparation_bytes = 1088;
        cap.maximum_preparation_multiply_terms = 548;
        cap.maximum_step_multiply_terms = 528;
        assert!(prepare(cap).is_ok());
        for which in 0..4 {
            let mut too_small = cap;
            match which {
                0 => too_small.maximum_free_dofs -= 1,
                1 => too_small.maximum_preparation_bytes -= 1,
                2 => too_small.maximum_preparation_multiply_terms -= 1,
                _ => too_small.maximum_step_multiply_terms -= 1,
            }
            assert!(matches!(
                prepare(too_small),
                Err(CompleteSpectralSolidError::Budget { .. })
            ));
        }
        let mut solid = prepare(controls()).unwrap();
        let before = bits(&solid);
        assert!(solid.advance_affine(&[-f64::MAX], &[f64::MAX], cx).is_err());
        assert!(solid.advance(&[], cx).is_err());
        assert!(solid.advance(&[f64::NAN], cx).is_err());
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.advance(&[0.0], cancelled),
                Err(CompleteSpectralSolidError::Cancelled)
            ));
            assert!(matches!(
                CompleteSpectralSolid::prepare_with_basis(
                    &pencil,
                    &[1.0],
                    &[1.0],
                    0.1,
                    &[0.1],
                    &[0.2],
                    controls(),
                    cancelled
                ),
                Err(CompleteSpectralSolidError::Cancelled)
            ));
        });
        assert_eq!(bits(&solid), before);
        let mut replay = prepare(controls()).unwrap();
        solid.advance_affine(&[0.2], &[-0.1], cx).unwrap();
        replay.advance_affine(&[0.2], &[-0.1], cx).unwrap();
        assert_eq!(bits(&solid), bits(&replay));
        // A deliberately loose basis residual cannot replace original K's
        // physical energy with the supplied diagonal eigenvalues.
        let mut loose = controls();
        loose.maximum_relative_operator_error = 0.02;
        loose.maximum_work_error_j = 1e-14;
        let mut mismatch = CompleteSpectralSolid::prepare_with_basis(
            &pencil,
            &[1.01],
            &[1.0],
            0.3,
            &[0.1],
            &[0.2],
            loose,
            cx,
        )
        .unwrap();
        let before = bits(&mismatch);
        assert!(matches!(
            mismatch.advance(&[0.0], cx),
            Err(CompleteSpectralSolidError::Energy { .. })
        ));
        assert_eq!(bits(&mismatch), before);
        let mut one = controls();
        one.maximum_steps = 1;
        let mut limited = prepare(one).unwrap();
        limited.advance(&[0.0], cx).unwrap();
        let before = bits(&limited);
        assert!(limited.advance(&[0.0], cx).is_err());
        assert_eq!(bits(&limited), before);
    });
}

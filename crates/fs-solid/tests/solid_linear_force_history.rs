//! Exact affine-in-time force histories in the complete solid coordinates.
//! Work oracles integrate force times velocity independently of total energy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls, ExactLinearSolidError};
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
                kernel_id: 30,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> ExactLinearSolidControls {
    ExactLinearSolidControls {
        maximum_free_dofs: 85,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_exponential_multiply_terms: 10_000_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_work_error_j: 1e-9,
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
fn near(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "actual={a:e}, expected={b:e}, allowance={allowance:e}"
    );
}
fn bits(s: &ExactLinearSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([s.time_s(), s.stored_energy_j(), s.kinetic_energy_j()])
        .map(f64::to_bits)
        .collect()
}

// q''+omega^2*q=f0+rate*t, unit modal mass. The integral is from zero to t.
fn ramp(omega: f64, q0: f64, v0: f64, f0: f64, rate: f64, t: f64) -> (f64, f64, f64) {
    let (sin, cos) = (omega * t).sin_cos();
    let (a, b) = (f0 / (omega * omega), rate / (omega * omega));
    let (c, d) = (q0 - a, (v0 - b) / omega);
    (
        a + b * t + c * cos + d * sin,
        b - omega * c * sin + omega * d * cos,
        a * t + 0.5 * b * t * t + c * sin / omega + d * (1.0 - cos) / omega,
    )
}

#[test]
fn g1_coupled_ramps_match_analytic_motion_and_independent_work() {
    with_cx(false, |cx| {
        // u=S*q; the original non-diagonal matrices retain both modes.
        let s = [[1.0, 0.4], [0.2, 0.8]];
        let inverse = [[0.8 / 0.72, -0.4 / 0.72], [-0.2 / 0.72, 1.0 / 0.72]];
        let omega = [20.0, 47.0];
        let (mut m, mut k) = ([0.0; 4], [0.0; 4]);
        for a in 0..2 {
            for b in 0..2 {
                for i in 0..2 {
                    m[2 * a + b] += inverse[i][a] * inverse[i][b];
                    k[2 * a + b] += inverse[i][a] * omega[i] * omega[i] * inverse[i][b];
                }
            }
        }
        let mut pencil = assembly(&k, &m, 2);
        pencil.free_dofs = vec![1, 5];
        let project =
            |q: [f64; 2]| -> [f64; 2] { core::array::from_fn(|i| s[i][0] * q[0] + s[i][1] * q[1]) };
        let q0 = [0.04, -0.025];
        let v0 = [0.1, -0.2];
        let f0 = [0.7, -0.2];
        let rate = [2.0, -1.3];
        let dt = 0.03;
        let load: [f64; 2] = core::array::from_fn(|i| s[0][i] * f0[0] + s[1][i] * f0[1]);
        let slope: [f64; 2] = core::array::from_fn(|i| s[0][i] * rate[0] + s[1][i] * rate[1]);
        let mut solid = ExactLinearSolid::prepare_affine_forces(
            &pencil,
            dt,
            &project(q0),
            &project(v0),
            1.0,
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(solid.free_dofs(), &[1, 5]);
        let mut trapezoid_discrepancy = 0.0_f64;
        for interval in 0..12 {
            let (a, b) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let force_a: [f64; 2] = core::array::from_fn(|i| f0[i] + rate[i] * a);
            let force_b: [f64; 2] = core::array::from_fn(|i| f0[i] + rate[i] * b);
            let before = solid.displacement_m().to_vec();
            let step = solid.advance_affine(&force_a, &force_b, cx).unwrap();
            let (mut q, mut v, mut work) = ([0.0; 2], [0.0; 2], 0.0);
            for i in 0..2 {
                let (qa, _, ia) = ramp(omega[i], q0[i], v0[i], load[i], slope[i], a);
                let (qb, vb, ib) = ramp(omega[i], q0[i], v0[i], load[i], slope[i], b);
                q[i] = qb;
                v[i] = vb;
                // Integral f*q' = [f*q] - rate*integral q, independently
                // evaluated from the continuous analytic displacement.
                work += (load[i] + slope[i] * b) * qb
                    - (load[i] + slope[i] * a) * qa
                    - slope[i] * (ib - ia);
            }
            for i in 0..2 {
                near(solid.displacement_m()[i], project(q)[i], 2e-12);
                near(solid.velocity_mps()[i], project(v)[i], 2e-10);
                near(
                    step.mean_velocity_mps[i] * dt,
                    solid.displacement_m()[i] - before[i],
                    1e-15,
                );
            }
            near(step.energy.supplied_work_j, work, 2e-12);
            near(step.energy.balance_error_j, 0.0, 2e-12);
            let trapezoid: f64 = (0..2)
                .map(|i| 0.5 * (force_a[i] + force_b[i]) * (solid.displacement_m()[i] - before[i]))
                .sum();
            trapezoid_discrepancy = trapezoid_discrepancy.max((trapezoid - work).abs());
        }
        assert!(
            trapezoid_discrepancy > 1e-5,
            "the fixture must distinguish exact ramp work from mean force times displacement"
        );
    });
}

#[test]
fn g1_zero_stiffness_has_exact_cubic_motion_and_polynomial_work() {
    with_cx(false, |cx| {
        let m = [2.0, 0.25, 0.25, 1.0];
        let inverse = |f: [f64; 2]| {
            [
                (f[0] - 0.25 * f[1]) / 1.9375,
                (2.0 * f[1] - 0.25 * f[0]) / 1.9375,
            ]
        };
        let f0 = [0.7, -0.3];
        let rate = [0.9, 0.4];
        let a0 = inverse(f0);
        let jerk = inverse(rate);
        let u0 = [0.4, -0.2];
        let v0 = [0.8, 0.3];
        let t = 0.4_f64;
        let f1 = core::array::from_fn::<_, 2, _>(|i| f0[i] + rate[i] * t);
        let mut solid = ExactLinearSolid::prepare_affine_forces(
            &assembly(&[0.0; 4], &m, 2),
            t,
            &u0,
            &v0,
            1.0,
            controls(),
            cx,
        )
        .unwrap();
        let step = solid.advance_affine(&f0, &f1, cx).unwrap();
        let mut work = 0.0;
        for i in 0..2 {
            near(
                solid.displacement_m()[i],
                u0[i] + v0[i] * t + 0.5 * a0[i] * t * t + jerk[i] * t.powi(3) / 6.0,
                2e-13,
            );
            near(
                solid.velocity_mps()[i],
                v0[i] + a0[i] * t + 0.5 * jerk[i] * t * t,
                2e-13,
            );
            work += f0[i] * v0[i] * t
                + 0.5 * (f0[i] * a0[i] + rate[i] * v0[i]) * t * t
                + (0.5 * f0[i] * jerk[i] + rate[i] * a0[i]) * t.powi(3) / 3.0
                + rate[i] * jerk[i] * t.powi(4) / 8.0;
        }
        assert_eq!(step.energy.stored_after_j, 0.0);
        near(step.energy.supplied_work_j, work, 2e-13);
        near(step.energy.balance_error_j, 0.0, 2e-13);
    });
}

#[test]
fn g1_affine_partition_invariance_includes_force_jumps_without_impulses() {
    with_cx(false, |cx| {
        let pencil = assembly(&[18.0], &[2.0], 1);
        let (dt, u0, v0) = (0.12, 0.04, -0.1);
        let make = |h| {
            ExactLinearSolid::prepare_affine_forces(&pencil, h, &[u0], &[v0], 1.0, controls(), cx)
                .unwrap()
        };
        let (mut coarse, mut fine) = (make(dt), make(dt / 4.0));
        let (mut q, mut v, mut coarse_work, mut fine_work) = (u0, v0, 0.0, 0.0);
        // The second start force deliberately jumps from1.1N to3N. A finite
        // force jump changes acceleration but introduces no velocity impulse.
        for (start, end) in [(-0.4, 1.1), (3.0, -2.3)] {
            let slope = (end - start) / dt;
            let (next_q, next_v, integral) = ramp(3.0, q, v, start / 2.0, slope / 2.0, dt);
            let work = end * next_q - start * q - slope * integral;
            let step = coarse.advance_affine(&[start], &[end], cx).unwrap();
            near(coarse.displacement_m()[0], next_q, 2e-13);
            near(coarse.velocity_mps()[0], next_v, 2e-13);
            near(step.energy.supplied_work_j, work, 2e-13);
            coarse_work += step.energy.supplied_work_j;
            for j in 0..4 {
                fine_work += fine
                    .advance_affine(
                        &[start + (end - start) * j as f64 / 4.0],
                        &[start + (end - start) * (j + 1) as f64 / 4.0],
                        cx,
                    )
                    .unwrap()
                    .energy
                    .supplied_work_j;
            }
            near(fine.displacement_m()[0], coarse.displacement_m()[0], 2e-13);
            near(fine.velocity_mps()[0], coarse.velocity_mps()[0], 2e-13);
            near(fine_work, coarse_work, 2e-13);
            (q, v) = (next_q, next_v);
        }
    });
}

#[test]
fn g0_affine_controls_and_failed_histories_preserve_accepted_state() {
    with_cx(false, |cx| {
        let pencil = assembly(&[4.0], &[2.0], 1);
        let make = |scale| {
            ExactLinearSolid::prepare_affine_forces(
                &pencil,
                0.03,
                &[0.1],
                &[0.2],
                scale,
                controls(),
                cx,
            )
        };
        for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(make(scale).is_err());
        }
        // Retaining force as source coordinates changes the complete bound
        // from 3*N+1 to 4*N+1; N=64 must refuse before dense preparation.
        let mut identity = vec![0.0; 64 * 64];
        for i in 0..64 {
            identity[64 * i + i] = 1.0;
        }
        let capped = ExactLinearSolid::prepare_affine_forces(
            &assembly(&identity, &identity, 64),
            0.03,
            &[0.0; 64],
            &[0.0; 64],
            1.0,
            controls(),
            cx,
        );
        assert!(matches!(
            capped,
            Err(ExactLinearSolidError::Budget {
                required: 257,
                maximum: 256,
            })
        ));
        let mut solid = make(1.0).unwrap();
        let before = bits(&solid);
        for (a, b) in [
            (vec![], vec![0.0]),
            (vec![0.0], vec![]),
            (vec![f64::NAN], vec![0.0]),
            (vec![0.0], vec![f64::INFINITY]),
            (vec![f64::MAX], vec![-f64::MAX]),
        ] {
            assert!(solid.advance_affine(&a, &b, cx).is_err());
            assert_eq!(bits(&solid), before);
        }
        assert!(solid.advance(&[0.3], cx).is_err());
        with_cx(true, |cancelled| {
            assert!(solid.advance_affine(&[0.0], &[0.5], cancelled).is_err());
        });
        assert_eq!(bits(&solid), before);
        let mut held =
            ExactLinearSolid::prepare(&pencil, 0.03, &[0.1], &[0.2], controls(), cx).unwrap();
        let held_before = bits(&held);
        assert!(held.advance_affine(&[0.0], &[0.5], cx).is_err());
        assert_eq!(bits(&held), held_before);
        let mut twin = make(1.0).unwrap();
        for (a, b) in [(0.0, 0.5), (1.2, -0.3)] {
            solid.advance_affine(&[a], &[b], cx).unwrap();
            twin.advance_affine(&[a], &[b], cx).unwrap();
        }
        assert_eq!(bits(&solid), bits(&twin));
    });
}

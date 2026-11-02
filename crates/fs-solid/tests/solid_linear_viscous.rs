//! Exact memoryless viscosity in the original coupled mechanical coordinates.
//! Independent continuous trajectories supply force work and viscosity loss.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls, ExactLinearSolidError};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

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
                kernel_id: 31,
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
        maximum_work_error_j: 1e-10,
        maximum_steps: 1000,
    }
}
fn matrix(a: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            coo.push(i, j, a[n * i + j]);
        }
    }
    coo.assemble()
}
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    TetElasticAssembly {
        stiffness: matrix(k, n),
        mass: matrix(m, n),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}
fn bits(solid: &ExactLinearSolid) -> Vec<u64> {
    solid
        .displacement_m()
        .iter()
        .chain(solid.velocity_mps())
        .copied()
        .chain([
            solid.time_s(),
            solid.stored_energy_j(),
            solid.kinetic_energy_j(),
            solid.cumulative_dissipated_energy_j(),
        ])
        .map(f64::to_bits)
        .collect()
}
fn scalar(m: f64, c: f64, k: f64, q0: f64, v0: f64, t: f64) -> (f64, f64) {
    let gamma = c / (2.0 * m);
    let delta = k / m - gamma * gamma;
    if delta > 0.0 {
        let omega = delta.sqrt();
        let (sin, cos) = (omega * t).sin_cos();
        let b = (v0 + gamma * q0) / omega;
        let q = q0 * cos + b * sin;
        let v = -omega * q0 * sin + omega * b * cos - gamma * q;
        let decay = (-gamma * t).exp();
        (decay * q, decay * v)
    } else if delta == 0.0 {
        let b = v0 + gamma * q0;
        let q = q0 + b * t;
        let decay = (-gamma * t).exp();
        (decay * q, decay * (b - gamma * q))
    } else {
        let root = (-delta).sqrt();
        let (r1, r2) = (-gamma + root, -gamma - root);
        let a = (v0 - r2 * q0) / (r1 - r2);
        let b = q0 - a;
        (
            a * (r1 * t).exp() + b * (r2 * t).exp(),
            r1 * a * (r1 * t).exp() + r2 * b * (r2 * t).exp(),
        )
    }
}
fn integral(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    // Independent composite Simpson integration of a continuous closed-form
    // trajectory, not a production endpoint or its energy difference.
    let panels = 1024;
    let h = (b - a) / panels as f64;
    let mut sum = f(a) + f(b);
    for i in 1..panels {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(a + i as f64 * h);
    }
    h * sum / 3.0
}

#[test]
fn g1_under_critical_and_overdamped_motion_matches_continuous_loss() {
    with_cx(false, |cx| {
        let (m, k, q0, v0, dt) = (2.0, 18.0, 0.2, -0.1, 0.037);
        let pencil = assembly(&[k], &[m], 1);
        for c in [2.0, 12.0, 20.0] {
            let mut solid = ExactLinearSolid::prepare_viscous(
                &pencil,
                &matrix(&[c], 1),
                dt,
                &[q0],
                &[v0],
                controls(),
                cx,
            )
            .unwrap();
            let mut total_loss = 0.0;
            for interval in 0..40 {
                let (a, b) = (interval as f64 * dt, (interval + 1) as f64 * dt);
                let (q, v) = scalar(m, c, k, q0, v0, b);
                let loss = integral(a, b, |t| c * scalar(m, c, k, q0, v0, t).1.powi(2));
                let step = solid.advance(&[0.0], cx).unwrap();
                total_loss += loss;
                near(solid.displacement_m()[0], q, 2e-13);
                near(solid.velocity_mps()[0], v, 2e-13);
                near(step.energy.dissipated_j, loss, 2e-13);
                near(step.energy.cumulative_dissipated_j, total_loss, 2e-12);
                near(step.energy.balance_error_j, 0.0, 2e-13);
                assert!(step.energy.dissipated_j >= 0.0);
            }
            near(solid.cumulative_dissipated_energy_j(), total_loss, 2e-12);
        }
    });
}

#[test]
fn g1_rigid_drag_retains_terminal_velocity_and_original_force_work() {
    with_cx(false, |cx| {
        let (m, c, f, q0, v0, dt) = (2.0, 3.0, 1.5, 0.2, -0.4, 0.15);
        let pencil = assembly(&[0.0], &[m], 1);
        let mut solid = ExactLinearSolid::prepare_viscous(
            &pencil,
            &matrix(&[c], 1),
            dt,
            &[q0],
            &[v0],
            controls(),
            cx,
        )
        .unwrap();
        let (lambda, terminal, departure) = (c / m, f / c, v0 - f / c);
        let position =
            |t: f64| q0 + terminal * t + departure * (1.0 - (-lambda * t).exp()) / lambda;
        for interval in 0..20 {
            let (a, b) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let exp1 = ((-lambda * a).exp() - (-lambda * b).exp()) / lambda;
            let exp2 = ((-2.0 * lambda * a).exp() - (-2.0 * lambda * b).exp()) / (2.0 * lambda);
            let loss = c
                * (terminal * terminal * (b - a)
                    + 2.0 * terminal * departure * exp1
                    + departure * departure * exp2);
            let step = solid.advance(&[f], cx).unwrap();
            near(solid.displacement_m()[0], position(b), 2e-13);
            near(
                solid.velocity_mps()[0],
                terminal + departure * (-lambda * b).exp(),
                2e-13,
            );
            near(
                step.energy.supplied_work_j,
                f * (position(b) - position(a)),
                2e-13,
            );
            near(step.energy.dissipated_j, loss, 2e-13);
            near(step.energy.balance_error_j, 0.0, 2e-13);
        }
    });
}

#[test]
fn g1_nonproportional_viscosity_preserves_the_coupled_consistent_pencil() {
    with_cx(false, |cx| {
        // Two genuine eigenvectors of the quadratic pencil, not of elastic K:
        // (r²I+rC+K)[1,0]=0 at r=-1 and (...) [1,1]=0 at r=-2.
        // C and K do not commute, so a modal-diagonal damping shortcut is wrong.
        let k = [[1.5, 0.5], [0.5, 2.5]];
        let c = [[2.5, 0.5], [0.5, 3.0]];
        let commutator =
            k[0][0] * c[0][1] + k[0][1] * c[1][1] - c[0][0] * k[0][1] - c[0][1] * k[1][1];
        assert!(f64::abs(commutator) > 0.1);
        let s = [[1.0, 0.4], [0.2, 0.8]];
        let inv = [[0.8 / 0.72, -0.4 / 0.72], [-0.2 / 0.72, 1.0 / 0.72]];
        let congruence = |a: [[f64; 2]; 2]| -> [f64; 4] {
            core::array::from_fn(|i| {
                let (row, col) = (i / 2, i % 2);
                (0..2)
                    .map(|j| {
                        (0..2)
                            .map(|l| inv[j][row] * a[j][l] * inv[l][col])
                            .sum::<f64>()
                    })
                    .sum()
            })
        };
        let transform =
            |q: [f64; 2]| -> [f64; 2] { core::array::from_fn(|i| s[i][0] * q[0] + s[i][1] * q[1]) };
        let (a, b, dt) = (0.1, -0.06, 0.043);
        let q = |t: f64| [a * (-t).exp() + b * (-2.0 * t).exp(), b * (-2.0 * t).exp()];
        let v = |t: f64| {
            [
                -a * (-t).exp() - 2.0 * b * (-2.0 * t).exp(),
                -2.0 * b * (-2.0 * t).exp(),
            ]
        };
        let pencil = assembly(&congruence(k), &congruence([[1.0, 0.0], [0.0, 1.0]]), 2);
        let mut solid = ExactLinearSolid::prepare_viscous(
            &pencil,
            &matrix(&congruence(c), 2),
            dt,
            &transform(q(0.0)),
            &transform(v(0.0)),
            controls(),
            cx,
        )
        .unwrap();
        for interval in 0..30 {
            let (start, end) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let exp_integral = |rate: f64| ((-rate * start).exp() - (-rate * end).exp()) / rate;
            let loss = 2.5 * a * a * exp_integral(2.0)
                + 12.0 * a * b * exp_integral(3.0)
                + 26.0 * b * b * exp_integral(4.0);
            let step = solid.advance(&[0.0, 0.0], cx).unwrap();
            for i in 0..2 {
                near(solid.displacement_m()[i], transform(q(end))[i], 2e-13);
                near(solid.velocity_mps()[i], transform(v(end))[i], 2e-13);
            }
            near(step.energy.dissipated_j, loss, 2e-13);
            near(step.energy.cumulative_balance_error_j, 0.0, 2e-13);
        }
    });
}

#[test]
fn g1_affine_force_work_and_viscous_loss_match_independent_integrals_and_partitions() {
    with_cx(false, |cx| {
        let (m, c, k, q0, v0, f0, rate, dt) = (2.0, 3.0, 18.0, 0.12, -0.07, 0.7, -0.9, 0.23);
        let pencil = assembly(&[k], &[m], 1);
        let damping = matrix(&[c], 1);
        let make = |h| {
            ExactLinearSolid::prepare_viscous_affine_forces(
                &pencil,
                &damping,
                h,
                &[q0],
                &[v0],
                1.0,
                controls(),
                cx,
            )
            .unwrap()
        };
        let (mut coarse, mut fine) = (make(dt), make(dt / 4.0));
        let slope = rate / k;
        let offset = (f0 - c * slope) / k;
        let solution = |t: f64| {
            let (q, v) = scalar(m, c, k, q0 - offset, v0 - slope, t);
            (offset + slope * t + q, slope + v)
        };
        let step = coarse.advance_affine(&[f0], &[f0 + rate * dt], cx).unwrap();
        let work = integral(0.0, dt, |t| (f0 + rate * t) * solution(t).1);
        let loss = integral(0.0, dt, |t| c * solution(t).1.powi(2));
        near(coarse.displacement_m()[0], solution(dt).0, 2e-13);
        near(coarse.velocity_mps()[0], solution(dt).1, 2e-13);
        near(step.energy.supplied_work_j, work, 2e-13);
        near(step.energy.dissipated_j, loss, 2e-13);
        near(step.energy.balance_error_j, 0.0, 2e-13);
        let mut fine_work = 0.0;
        for i in 0..4 {
            let (a, b) = (i as f64 * dt / 4.0, (i + 1) as f64 * dt / 4.0);
            fine_work += fine
                .advance_affine(&[f0 + rate * a], &[f0 + rate * b], cx)
                .unwrap()
                .energy
                .supplied_work_j;
        }
        near(fine_work, work, 2e-13);
        near(fine.cumulative_dissipated_energy_j(), loss, 2e-13);
        near(fine.displacement_m()[0], coarse.displacement_m()[0], 2e-13);
        near(fine.velocity_mps()[0], coarse.velocity_mps()[0], 2e-13);
        let endpoint_shortcut = dt * c * ((v0 + coarse.velocity_mps()[0]) * 0.5).powi(2);
        assert!(
            (loss - endpoint_shortcut).abs() > 1e-5,
            "fixture must distinguish original continuous viscosity from midpoint damping"
        );
    });
}

#[test]
fn g0_viscosity_admission_budgets_and_failed_steps_preserve_accepted_history() {
    with_cx(false, |cx| {
        let pencil = assembly(&[18.0], &[2.0], 1);
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                ExactLinearSolid::prepare_viscous(
                    &pencil,
                    &matrix(&[bad], 1),
                    0.1,
                    &[0.1],
                    &[0.2],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        let two = assembly(&[2.0, 0.0, 0.0, 3.0], &[1.0, 0.0, 0.0, 1.0], 2);
        for bad in [[1.0, 2.0, 2.0, 1.0], [1.0, 0.1, 0.2, 1.0]] {
            assert!(
                ExactLinearSolid::prepare_viscous(
                    &two,
                    &matrix(&bad, 2),
                    0.1,
                    &[0.1, 0.2],
                    &[0.0, 0.0],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        assert!(
            ExactLinearSolid::prepare_viscous(
                &pencil,
                &matrix(&[1.0, 0.0, 0.0, 1.0], 2),
                0.1,
                &[0.1],
                &[0.2],
                controls(),
                cx
            )
            .is_err()
        );
        let damping = matrix(&[3.0], 1);
        let mut limited = controls();
        limited.maximum_preparation_bytes = 52 * 4 * 4 * 8;
        assert!(
            ExactLinearSolid::prepare_viscous(&pencil, &damping, 0.1, &[0.1], &[0.2], limited, cx)
                .is_ok()
        );
        limited.maximum_preparation_bytes -= 1;
        assert!(matches!(
            ExactLinearSolid::prepare_viscous(&pencil, &damping, 0.1, &[0.1], &[0.2], limited, cx),
            Err(ExactLinearSolidError::Budget { .. })
        ));
        limited = controls();
        limited.maximum_exponential_multiply_terms = 100_000;
        assert!(
            ExactLinearSolid::prepare_viscous_affine_forces(
                &pencil,
                &damping,
                0.1,
                &[0.1],
                &[0.2],
                1.0,
                limited,
                cx
            )
            .is_err()
        );
        let make = || {
            ExactLinearSolid::prepare_viscous(
                &pencil,
                &damping,
                0.1,
                &[0.1],
                &[0.2],
                controls(),
                cx,
            )
            .unwrap()
        };
        let (mut solid, mut replay) = (make(), make());
        solid.advance(&[0.7], cx).unwrap();
        replay.advance(&[0.7], cx).unwrap();
        assert_eq!(bits(&solid), bits(&replay));
        let before = bits(&solid);
        for force in [vec![], vec![f64::NAN], vec![f64::MAX]] {
            assert!(solid.advance(&force, cx).is_err());
            assert_eq!(bits(&solid), before);
        }
        assert!(solid.advance_affine(&[0.0], &[0.0], cx).is_err());
        assert_eq!(bits(&solid), before);
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.advance(&[0.0], cancelled),
                Err(ExactLinearSolidError::Cancelled)
            ));
            assert!(matches!(
                ExactLinearSolid::prepare_viscous(
                    &pencil,
                    &damping,
                    0.1,
                    &[0.1],
                    &[0.2],
                    controls(),
                    cancelled
                ),
                Err(ExactLinearSolidError::Cancelled)
            ));
        });
        assert_eq!(bits(&solid), before);
        let zero = matrix(&[0.0], 1);
        let mut unmodified =
            ExactLinearSolid::prepare(&pencil, 0.1, &[0.1], &[0.2], controls(), cx).unwrap();
        let mut explicit_zero =
            ExactLinearSolid::prepare_viscous(&pencil, &zero, 0.1, &[0.1], &[0.2], controls(), cx)
                .unwrap();
        unmodified.advance(&[0.7], cx).unwrap();
        let step = explicit_zero.advance(&[0.7], cx).unwrap();
        assert_eq!(bits(&unmodified), bits(&explicit_zero));
        assert_eq!(step.energy.dissipated_j, 0.0);
    });
}

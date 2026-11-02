//! Complete-coordinate linear solid motion through the original exact owner.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls};
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
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    });
}
fn controls() -> ExactLinearSolidControls {
    ExactLinearSolidControls {
        maximum_free_dofs: 85,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_exponential_multiply_terms: 10_000_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 1_000_001,
    }
}
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    let mut stiffness = Coo::new(n, n);
    let mut mass = Coo::new(n, n);
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
fn bits(s: &ExactLinearSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([s.time_s(), s.stored_energy_j(), s.kinetic_energy_j()])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn g1_five_khz_exact_phase_survives_a_million_complete_state_steps() {
    with_cx(false, |cx| {
        let omega = core::f64::consts::TAU * 5000.0;
        let dt = 1.0 / 48000.0;
        let pencil = assembly(&[omega * omega], &[1.0], 1);
        let mut solid =
            ExactLinearSolid::prepare(&pencil, dt, &[1e-6], &[0.0], controls(), cx).unwrap();
        let energy0 = solid.stored_energy_j();
        let (mut angle, mut correction, mut drift) = (0.0, 0.0, 0.0_f64);
        let (mut old_q, mut old_p) = (omega * solid.displacement_m()[0], -solid.velocity_mps()[0]);
        for _ in 0..1_000_000 {
            solid.advance(&[0.0], cx).unwrap();
            let (q, p) = (omega * solid.displacement_m()[0], -solid.velocity_mps()[0]);
            let delta = (old_q * p - old_p * q).atan2(old_q * q + old_p * p);
            let add = delta - correction;
            let next = angle + add;
            correction = (next - angle) - add;
            angle = next;
            drift = drift.max(
                (solid.stored_energy_j() + solid.kinetic_energy_j() - energy0).abs() / energy0,
            );
            old_q = q;
            old_p = p;
        }
        let cents = 1200.0 * (angle / (1_000_000.0 * omega * dt)).log2();
        near(cents, 0.0, 1e-6);
        assert!(drift < 1e-8, "relative energy drift={drift:e}");
        // Independent Cayley recurrence is a falsifier for accidental midpoint.
        let x = 0.5 * omega * dt;
        let c = (1.0 - x * x) / (1.0 + x * x);
        let s = 2.0 * x / (1.0 + x * x);
        let measured = 1200.0 * (s.atan2(c) / (omega * dt)).log2();
        let warped = 1200.0 * (x.atan() / x).log2();
        near(measured, warped, 1e-10);
        assert!(measured.abs() > 50.0);
        println!(
            "exact 5kHz/48kHz million steps: pitch={cents:e} cents, energy_relative={drift:e}; midpoint={measured:e} cents"
        );
    });
}

#[test]
fn g1_coupled_mass_and_stiffness_retain_both_frequencies_and_held_force_work() {
    with_cx(false, |cx| {
        // u=S*q, M=S^-T*S^-1, K=S^-T*diag(w²)*S^-1. The two analytic
        // coordinates are independent oracles; the runtime sees only full K,M.
        let s = [[1.0, 0.4], [0.2, 0.8]];
        let det = 0.72;
        let inv = [[0.8 / det, -0.4 / det], [-0.2 / det, 1.0 / det]];
        let omega = [
            core::f64::consts::TAU * 5000.0,
            core::f64::consts::TAU * 11000.0,
        ];
        let mut m = [0.0; 4];
        let mut k = [0.0; 4];
        for a in 0..2 {
            for b in 0..2 {
                for i in 0..2 {
                    m[2 * a + b] += inv[i][a] * inv[i][b];
                    k[2 * a + b] += inv[i][a] * omega[i] * omega[i] * inv[i][b];
                }
            }
        }
        let mut pencil = assembly(&k, &m, 2);
        pencil.free_dofs = vec![1, 5];
        let q0 = [1e-5, -2e-5];
        let v0 = [0.03, -0.08];
        let force = [100.0, -60.0];
        let project =
            |q: [f64; 2]| -> [f64; 2] { core::array::from_fn(|i| s[i][0] * q[0] + s[i][1] * q[1]) };
        let load: [f64; 2] = core::array::from_fn(|i| s[0][i] * force[0] + s[1][i] * force[1]);
        let mut solid = ExactLinearSolid::prepare(
            &pencil,
            1.0 / 48000.0,
            &project(q0),
            &project(v0),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(solid.free_dofs(), &[1, 5]);
        for _ in 0..192 {
            let before = solid.displacement_m().to_vec();
            let result = solid.advance(&force, cx).unwrap();
            let t = result.time_s;
            let mut q = [0.0; 2];
            let mut v = [0.0; 2];
            for i in 0..2 {
                let (sin, cos) = (omega[i] * t).sin_cos();
                let equilibrium = load[i] / (omega[i] * omega[i]);
                q[i] = equilibrium + (q0[i] - equilibrium) * cos + v0[i] / omega[i] * sin;
                v[i] = -omega[i] * (q0[i] - equilibrium) * sin + v0[i] * cos;
            }
            for i in 0..2 {
                near(solid.displacement_m()[i], project(q)[i], 2e-16);
                near(solid.velocity_mps()[i], project(v)[i], 2e-11);
                near(
                    result.mean_velocity_mps[i] * solid.duration_s(),
                    solid.displacement_m()[i] - before[i],
                    1e-20,
                );
            }
            let work: f64 = (0..2)
                .map(|i| force[i] * (solid.displacement_m()[i] - before[i]))
                .sum();
            near(result.energy.supplied_work_j, work, 1e-15);
            near(result.energy.balance_error_j, 0.0, 2e-12);
        }
    });
}

#[test]
fn g1_zero_stiffness_keeps_rigid_motion_and_constant_acceleration() {
    with_cx(false, |cx| {
        let mass = [2.0, 0.25, 0.25, 1.0];
        let force = [0.7, -0.3];
        let determinant = 1.9375;
        let acceleration = [
            (force[0] - 0.25 * force[1]) / determinant,
            (2.0 * force[1] - 0.25 * force[0]) / determinant,
        ];
        let u0 = [0.4, -0.2];
        let v0 = [0.8, 0.3];
        let pencil = assembly(&[0.0; 4], &mass, 2);
        let mut solid = ExactLinearSolid::prepare(&pencil, 0.02, &u0, &v0, controls(), cx).unwrap();
        for _ in 0..8 {
            let step = solid.advance(&force, cx).unwrap();
            let t = step.time_s;
            for i in 0..2 {
                near(
                    solid.displacement_m()[i],
                    u0[i] + t * v0[i] + 0.5 * t * t * acceleration[i],
                    1e-13,
                );
                near(solid.velocity_mps()[i], v0[i] + t * acceleration[i], 1e-13);
            }
            assert_eq!(step.energy.stored_after_j, 0.0);
            near(
                step.energy.supplied_work_j,
                step.energy.kinetic_after_j - step.energy.kinetic_before_j,
                1e-14,
            );
        }
    });
}

#[test]
fn g1_reference_p1_and_signed_p2_consistent_mass_pencils_are_not_lumped() {
    with_cx(false, |cx| {
        let edges = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];
        for nodes in [4, 10] {
            let n = 3 * nodes;
            let mut mass = vec![0.0; n * n];
            for a in 0..nodes {
                for b in 0..nodes {
                    let fraction = if nodes == 4 {
                        if a == b { 0.1 } else { 0.05 }
                    } else if a < 4 && b < 4 {
                        if a == b { 6.0 / 420.0 } else { 1.0 / 420.0 }
                    } else if a < 4 {
                        if edges[b - 4].contains(&a) {
                            -4.0 / 420.0
                        } else {
                            -6.0 / 420.0
                        }
                    } else if b < 4 {
                        if edges[a - 4].contains(&b) {
                            -4.0 / 420.0
                        } else {
                            -6.0 / 420.0
                        }
                    } else if a == b {
                        32.0 / 420.0
                    } else if edges[a - 4].iter().any(|i| edges[b - 4].contains(i)) {
                        16.0 / 420.0
                    } else {
                        8.0 / 420.0
                    };
                    for c in 0..3 {
                        mass[(3 * a + c) * n + 3 * b + c] = fraction;
                    }
                }
            }
            // A distributed elastic foundation K=w²M has a known full-space
            // trajectory for every shape, including quadratic nonaffine fields.
            let omega = 1300.0;
            let k: Vec<_> = mass.iter().map(|v| omega * omega * v).collect();
            let pencil = assembly(&k, &mass, n);
            let u: Vec<_> = (0..n).map(|i| ((i * 7 % 17) as f64 - 8.0) * 1e-5).collect();
            let mut solid =
                ExactLinearSolid::prepare(&pencil, 1e-4, &u, &vec![0.0; n], controls(), cx)
                    .unwrap();
            let initial = solid.stored_energy_j();
            for _ in 0..8 {
                solid.advance(&vec![0.0; n], cx).unwrap();
                let (sin, cos) = (omega * solid.time_s()).sin_cos();
                for i in 0..n {
                    near(solid.displacement_m()[i], u[i] * cos, 3e-16);
                    near(solid.velocity_mps()[i], -omega * u[i] * sin, 3e-12);
                }
                near(
                    solid.stored_energy_j() + solid.kinetic_energy_j(),
                    initial,
                    1e-13,
                );
            }
        }
    });
}

#[test]
fn g0_invalid_operators_controls_and_failed_candidates_do_not_publish_state() {
    with_cx(false, |cx| {
        let good = assembly(&[4.0, 0.0, 0.0, 9.0], &[1.0, 0.0, 0.0, 2.0], 2);
        let make =
            |config| ExactLinearSolid::prepare(&good, 0.03, &[0.1, 0.2], &[0.0, 0.0], config, cx);
        for c in [
            ExactLinearSolidControls {
                maximum_free_dofs: 1,
                ..controls()
            },
            ExactLinearSolidControls {
                maximum_preparation_bytes: 1,
                ..controls()
            },
            ExactLinearSolidControls {
                maximum_exponential_multiply_terms: 1,
                ..controls()
            },
            ExactLinearSolidControls {
                maximum_work_error_j: f64::NAN,
                ..controls()
            },
        ] {
            assert!(make(c).is_err());
        }
        for (k, m) in [
            ([4.0, 1.0, 0.0, 9.0], [1.0, 0.0, 0.0, 2.0]),
            ([-4.0, 0.0, 0.0, 9.0], [1.0, 0.0, 0.0, 2.0]),
            ([4.0, 0.0, 0.0, 9.0], [1.0, 0.0, 0.0, -2.0]),
        ] {
            assert!(
                ExactLinearSolid::prepare(
                    &assembly(&k, &m, 2),
                    0.03,
                    &[0.0; 2],
                    &[0.0; 2],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        let mut solid = make(controls()).unwrap();
        let before = bits(&solid);
        assert!(solid.advance(&[f64::NAN, 0.0], cx).is_err());
        assert_eq!(bits(&solid), before);
        assert!(solid.advance(&[0.0], cx).is_err());
        assert_eq!(bits(&solid), before);
        with_cx(true, |cancelled| {
            assert!(solid.advance(&[0.0; 2], cancelled).is_err());
        });
        assert_eq!(bits(&solid), before);
        let mut twin = make(controls()).unwrap();
        for _ in 0..8 {
            solid.advance(&[0.3, -0.1], cx).unwrap();
            twin.advance(&[0.3, -0.1], cx).unwrap();
        }
        assert_eq!(bits(&solid), bits(&twin));
        let mut one = make(ExactLinearSolidControls {
            maximum_steps: 1,
            ..controls()
        })
        .unwrap();
        one.advance(&[0.0; 2], cx).unwrap();
        let before = bits(&one);
        assert!(one.advance(&[0.0; 2], cx).is_err());
        assert_eq!(bits(&one), before);
    });
}

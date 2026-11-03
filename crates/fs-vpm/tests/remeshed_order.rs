//! G0/G1/G3/G4: signed inviscid stages inside one scheduled-remesh transaction.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::{
    CompositionOrder,
    ambient::{self, AffineFlow2, AmbientMomentChange},
    remeshed::{self, RemeshSchedule, RemeshedVpmBudget, RemeshedVpmError},
};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::remesh::{CartesianGrid2, RemeshBudget, RemeshError, remesh};
use fs_vpm::{VortexParticle as Particle, VpmBudget, VpmError};
use std::sync::atomic::{AtomicU64, Ordering};

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let result = pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 53,
                kernel_id: 13,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        ))
    });
    assert!(pool.stats().quiescent());
    result
}

fn caps() -> RemeshedVpmBudget {
    RemeshedVpmBudget {
        max_input_particles: 1000,
        max_grid_particles: 1000,
        max_steps: 100,
        max_pair_evaluations: 2_000_000,
        max_work_units: 3_000_000,
        max_live_bytes: 64_000,
    }
}

fn remap_caps() -> RemeshBudget {
    RemeshBudget {
        max_input_particles: 1000,
        max_grid_particles: 1000,
        max_work_units: 3_000_000,
        max_live_bytes: 64_000,
    }
}

fn grids() -> [CartesianGrid2; 3] {
    [
        CartesianGrid2::new([-0.3; 2], [0.1; 2], [7; 2]).unwrap(),
        CartesianGrid2::new([-0.5, -0.45], [0.08, 0.09], [14, 13]).unwrap(),
        CartesianGrid2::new([-0.8; 2], [0.1, 0.11], [18, 17]).unwrap(),
    ]
}

fn particles() -> [Particle; 3] {
    [
        Particle::new([-0.08, -0.03], 0.04),
        Particle::new([0.06, 0.05], -0.03),
        Particle::new([0.04, 0.07], 0.01),
    ]
}

fn flows() -> [AffineFlow2; 2] {
    [
        AffineFlow2::new([0.01, -0.03], [0.02, -0.01], [0.0, 0.4, 0.0]).unwrap(),
        AffineFlow2::new([-0.02, 0.04], [-0.01, 0.02], [0.05, -0.35, 0.25]).unwrap(),
    ]
}

fn moments(particles: &[Particle]) -> [f64; 4] {
    let mut result = [0.0; 4];
    for p in particles {
        let [x, y] = p.pos;
        for (out, monomial) in result.iter_mut().zip([1.0, x, y, x * x + y * y]) {
            *out += p.circulation * monomial;
        }
    }
    result
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-13, "{a:e} != {b:e}");
}

#[test]
fn g0_changing_grids_and_held_fields_match_existing_order_owners_and_counts() {
    let input = particles();
    let grids = grids();
    let flows = flows();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            let combined = remeshed::simulate_with_affine_flows_order_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                order,
                schedule,
                &flows,
                caps(),
            )
            .unwrap();
            let mut explicit = remesh(&input, grids[0], remap_caps(), || false).unwrap();
            let mut source = [0.0; 3];
            for (k, steps) in [2, 1].into_iter().enumerate() {
                let moved = ambient::simulate_with_order_cx(
                    cx,
                    &explicit.particles,
                    0.001,
                    steps,
                    0.1,
                    flows[k],
                    order,
                    VpmBudget::new(1000, 10, 1_000_000, 64_000),
                )
                .unwrap();
                for (out, delta) in source.iter_mut().zip([
                    moved.ambient_moment_change.first_moment_m3_s[0],
                    moved.ambient_moment_change.first_moment_m3_s[1],
                    moved.ambient_moment_change.quadratic_moment_m4_s,
                ]) {
                    *out += delta;
                }
                explicit =
                    remesh(&moved.run.particles, grids[k + 1], remap_caps(), || false).unwrap();
            }
            let run = &combined.run;
            assert_eq!(run.particles, explicit.particles);
            assert_eq!(run.areas_m2, explicit.areas_m2);
            assert_eq!(run.grid, grids[2]);
            assert_eq!(run.steps_completed, 3);
            assert_eq!(run.remeshes_completed, 3);
            let (sweeps, passes) = match order {
                CompositionOrder::Second => (1, 2),
                CompositionOrder::Fourth => (3, 4),
            };
            assert_eq!(run.pair_evaluations, sweeps * (2 * 49 * 48 + 182 * 181));
            assert_eq!(combined.source_particle_updates, passes * (2 * 49 + 182));
            assert_eq!(run.stencil_deposits, 16 * (3 + 49 + 182));
            let remap_work = 17 * (3 + 49 + 182) + 2 * (49 + 182 + 306);
            assert_eq!(
                run.work_units,
                remap_work + run.pair_evaluations + 3 + combined.source_particle_updates
            );
            assert_eq!(run.budget.planned_cost, run.work_units);
            assert_eq!(run.budget.cost_charged, run.work_units);
            assert_eq!(run.peak_live_bytes, 2 * 306 * (24 + 8));
            for (actual, expected) in [
                combined.ambient_moment_change.first_moment_m3_s[0],
                combined.ambient_moment_change.first_moment_m3_s[1],
                combined.ambient_moment_change.quadratic_moment_m4_s,
            ]
            .into_iter()
            .zip(source)
            {
                near(actual, expected);
            }
            if order == CompositionOrder::Second {
                let original = remeshed::simulate_with_affine_flows_cx(
                    cx,
                    &input,
                    0.001,
                    3,
                    0.1,
                    PseDiffusion::new(0.0, 0.16).unwrap(),
                    schedule,
                    &flows,
                    caps(),
                )
                .unwrap();
                assert_eq!(combined, original);
            }
        }
    });
}

#[test]
fn g0_piecewise_analytic_moment_law_and_signed_impulse_survive_changing_resolution() {
    let dipole = [
        Particle::new([-0.07, 0.04], 0.08),
        Particle::new([0.06, -0.03], -0.08),
    ];
    let nonzero = particles();
    let grids = grids();
    let shear = AffineFlow2::new([0.02, -0.01], [0.2, -0.1], [0.0, 0.7, 0.0]).unwrap();
    let rotation = AffineFlow2::new([-0.03, 0.02], [0.0; 2], [0.0, -0.8, 0.8]).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for input in [dipole.as_slice(), nonzero.as_slice()] {
            let before = moments(input);
            let schedule = RemeshSchedule {
                every_steps: 2,
                grids: &grids,
            };
            let unforced = remeshed::simulate_with_order_cx(
                cx,
                input,
                0.015,
                3,
                0.1,
                CompositionOrder::Fourth,
                schedule,
                caps(),
            )
            .unwrap();
            for (actual, expected) in moments(&unforced.particles).into_iter().zip(before) {
                near(actual, expected);
            }
            let run = remeshed::simulate_with_affine_flows_order_cx(
                cx,
                input,
                0.015,
                3,
                0.1,
                CompositionOrder::Fourth,
                schedule,
                &[shear, rotation],
                caps(),
            )
            .unwrap();
            // Integrate M'=A(M-Gamma*origin)+Gamma*U independently of pair
            // splitting, source exponentials and per-particle moment receipts.
            let gamma = before[0];
            let t = 0.03;
            let mx = before[1]
                + 0.7 * t * (before[2] + gamma * 0.01)
                + gamma * (0.2 * t - 0.5 * 0.7 * 0.1 * t * t);
            let my = before[2] - gamma * 0.1 * t;
            let (sin, cos) = (0.8_f64 * 0.015).sin_cos();
            let relative = [mx + gamma * 0.03, my - gamma * 0.02];
            let expected = [
                -gamma * 0.03 + cos * relative[0] - sin * relative[1],
                gamma * 0.02 + sin * relative[0] + cos * relative[1],
            ];
            let after = moments(&run.run.particles);
            near(after[0], gamma);
            for axis in 0..2 {
                near(after[axis + 1], expected[axis]);
                near(
                    run.ambient_moment_change.first_moment_m3_s[axis],
                    expected[axis] - before[axis + 1],
                );
            }
            near(
                after[3] - before[3],
                run.ambient_moment_change.quadratic_moment_m4_s,
            );
        }
    });
}

// Independent ODE and RK4 reference: no production pair/source time map.
fn derivative(p: &[Particle], forced: bool) -> Vec<[f64; 2]> {
    p.iter()
        .enumerate()
        .map(|(i, target)| {
            let [x, y] = [target.pos[0] - 0.02, target.pos[1] + 0.03];
            let mut velocity = if forced {
                [0.12 + 0.25 * x + 0.8 * y, -0.08 - 0.35 * x - 0.25 * y]
            } else {
                [0.0; 2]
            };
            for (j, source) in p.iter().enumerate() {
                if i == j {
                    continue;
                }
                let [x, y] = [target.pos[0] - source.pos[0], target.pos[1] - source.pos[1]];
                let f =
                    source.circulation / (std::f64::consts::TAU * (x * x + y * y + 0.15 * 0.15));
                velocity[0] -= f * y;
                velocity[1] += f * x;
            }
            velocity
        })
        .collect()
}

fn reference(initial: &[Particle], time: f64, steps: usize, forced: bool) -> Vec<Particle> {
    let mut p = initial.to_vec();
    let dt = time / steps as f64;
    for _ in 0..steps {
        let shift = |k: &[[f64; 2]], fraction: f64| -> Vec<_> {
            p.iter()
                .zip(k)
                .map(|(p, k)| {
                    Particle::new(
                        [
                            p.pos[0] + fraction * dt * k[0],
                            p.pos[1] + fraction * dt * k[1],
                        ],
                        p.circulation,
                    )
                })
                .collect()
        };
        let k1 = derivative(&p, forced);
        let k2 = derivative(&shift(&k1, 0.5), forced);
        let k3 = derivative(&shift(&k2, 0.5), forced);
        let k4 = derivative(&shift(&k3, 1.0), forced);
        for (i, p) in p.iter_mut().enumerate() {
            for a in 0..2 {
                p.pos[a] += dt * (k1[i][a] + 2.0 * k2[i][a] + 2.0 * k3[i][a] + k4[i][a]) / 6.0;
            }
        }
    }
    p
}

fn strength_error(a: &[Particle], b: &[Particle]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(a, b)| (a.circulation - b.circulation).abs())
        .fold(0.0, f64::max)
}

#[test]
fn g1_fixed_endpoint_remaps_retain_fourth_order_accuracy_at_lower_equal_work_error() {
    let input = [
        Particle::new([-0.24, -0.11], 0.8),
        Particle::new([0.28, -0.26], -0.3),
        Particle::new([0.18, 0.27], 0.25),
        Particle::new([-0.22, 0.22], -0.75),
    ];
    let grids = [
        CartesianGrid2::new([-0.75; 2], [0.25; 2], [7; 2]).unwrap(),
        CartesianGrid2::new([-1.8; 2], [0.3; 2], [13; 2]).unwrap(),
    ];
    let initial = remesh(&input, grids[0], remap_caps(), || false).unwrap();
    let time = 0.5;
    for forced in [false, true] {
        let flow = if forced {
            AffineFlow2::new([0.02, -0.03], [0.12, -0.08], [0.25, 0.8, -0.35]).unwrap()
        } else {
            AffineFlow2::default()
        };
        let truth = remesh(
            &reference(&initial.particles, time, 4096, forced),
            grids[1],
            remap_caps(),
            || false,
        )
        .unwrap();
        let coarser = remesh(
            &reference(&initial.particles, time, 2048, forced),
            grids[1],
            remap_caps(),
            || false,
        )
        .unwrap();
        assert!(strength_error(&truth.particles, &coarser.particles) < 5e-13);
        let evolve = |steps, order| {
            with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
                remeshed::simulate_with_affine_flows_order_cx(
                    cx,
                    &input,
                    time / steps as f64,
                    steps,
                    0.15,
                    order,
                    RemeshSchedule {
                        every_steps: steps,
                        grids: &grids,
                    },
                    &[flow],
                    caps(),
                )
            })
            .unwrap()
        };
        let mut errors = [0.0; 3];
        for (i, steps) in [4, 8, 16].into_iter().enumerate() {
            let run = evolve(steps, CompositionOrder::Fourth);
            errors[i] = strength_error(&run.run.particles, &truth.particles);
            assert_eq!(run.run.remeshes_completed, 2);
        }
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!(
                (12.0..20.0).contains(&ratio),
                "forced={forced}, errors={errors:?}, ratio={ratio}"
            );
        }
        let second = evolve(48, CompositionOrder::Second);
        let fourth = evolve(16, CompositionOrder::Fourth);
        assert_eq!(second.run.pair_evaluations, fourth.run.pair_evaluations);
        assert!(fourth.run.work_units <= second.run.work_units);
        assert!(
            errors[2] < strength_error(&second.run.particles, &truth.particles) / 4.0,
            "forced={forced}, fourth={errors:?}, second={}",
            strength_error(&second.run.particles, &truth.particles)
        );
        eprintln!(
            "scheduled refinement forced={forced}: {errors:?}; equal-pair-work second={:e}",
            strength_error(&second.run.particles, &truth.particles)
        );
    }
}

#[test]
fn g0_zero_sources_time_and_steps_keep_original_receipts_and_empty_quadrature() {
    let input = particles();
    let grids = grids();
    let law = PseDiffusion::new(0.0, 0.16).unwrap();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let huge = AffineFlow2::new([0.0; 2], [0.0; 2], [1e308, 0.0, 0.0]).unwrap();
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            for (dt, flows) in [(0.001, [AffineFlow2::default(); 2]), (0.0, [huge; 2])] {
                let plain = remeshed::simulate_with_order_cx(
                    cx,
                    &input,
                    dt,
                    3,
                    0.1,
                    order,
                    schedule,
                    caps(),
                )
                .unwrap();
                let source = remeshed::simulate_with_affine_flows_order_cx(
                    cx,
                    &input,
                    dt,
                    3,
                    0.1,
                    order,
                    schedule,
                    &flows,
                    caps(),
                )
                .unwrap();
                assert_eq!(source.run, plain);
                assert_eq!(source.source_particle_updates, 0);
                assert_eq!(source.ambient_moment_change, AmbientMomentChange::default());
                if order == CompositionOrder::Second || dt == 0.0 {
                    assert_eq!(
                        plain,
                        remeshed::simulate_with_cx(cx, &input, dt, 3, 0.1, law, schedule, caps())
                            .unwrap()
                    );
                }
            }
            let initial_only = RemeshSchedule {
                every_steps: 2,
                grids: &grids[..1],
            };
            let initial = remeshed::simulate_with_affine_flows_order_cx(
                cx,
                &input,
                0.001,
                0,
                0.1,
                order,
                initial_only,
                &[],
                caps(),
            )
            .unwrap();
            assert_eq!(
                initial.run,
                remeshed::simulate_with_cx(cx, &input, 0.001, 0, 0.1, law, initial_only, caps())
                    .unwrap()
            );
            assert_eq!(initial.run.peak_live_bytes, 49 * 32);
        }
        let empty = remeshed::simulate_with_affine_flows_order_cx(
            cx,
            &[],
            0.001,
            3,
            0.1,
            CompositionOrder::Fourth,
            schedule,
            &flows(),
            caps(),
        )
        .unwrap();
        assert!(empty.run.particles.iter().all(|p| p.circulation == 0.0));
        assert_eq!(empty.source_particle_updates, 4 * (2 * 49 + 182));
        assert_eq!(empty.ambient_moment_change, AmbientMomentChange::default());
        let mixed = remeshed::simulate_with_affine_flows_order_cx(
            cx,
            &input,
            0.001,
            3,
            0.1,
            CompositionOrder::Fourth,
            schedule,
            &[flows()[0], AffineFlow2::default()],
            caps(),
        )
        .unwrap();
        assert_eq!(mixed.source_particle_updates, 4 * 2 * 49);
    });
}

#[test]
fn g3_whole_signed_stage_work_and_original_two_banks_are_admitted_up_front() {
    let input = particles();
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let flow = flows();
    let run = |cx: &Cx<'_>, local| {
        remeshed::simulate_with_affine_flows_order_cx(
            cx,
            &input,
            0.001,
            3,
            0.1,
            CompositionOrder::Fourth,
            schedule,
            &flow,
            local,
        )
    };
    let reference = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        run(cx, caps())
    })
    .unwrap();
    let exact = RemeshedVpmBudget {
        max_input_particles: input.len(),
        max_grid_particles: 306,
        max_steps: 3,
        max_pair_evaluations: reference.run.pair_evaluations,
        max_work_units: reference.run.work_units,
        max_live_bytes: reference.run.peak_live_bytes,
    };
    let admitted = with_cx(
        &CancelGate::new_clock_free(),
        Budget::INFINITE.with_cost_quota(exact.max_work_units),
        |cx| run(cx, exact),
    )
    .unwrap();
    assert_eq!(admitted.run.particles, reference.run.particles);
    assert_eq!(admitted.run.budget.cost_charged, exact.max_work_units);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for local in [
            RemeshedVpmBudget {
                max_pair_evaluations: exact.max_pair_evaluations - 1,
                ..exact
            },
            RemeshedVpmBudget {
                max_work_units: exact.max_work_units - 1,
                ..exact
            },
            RemeshedVpmBudget {
                max_live_bytes: exact.max_live_bytes - 1,
                ..exact
            },
        ] {
            assert!(matches!(
                run(cx, local),
                Err(RemeshedVpmError::Evolution(VpmError::CapExceeded { .. }))
            ));
        }
        assert!(
            remeshed::simulate_with_order_cx(
                cx,
                &input,
                -0.001,
                3,
                0.1,
                CompositionOrder::Fourth,
                schedule,
                caps()
            )
            .is_err()
        );
        assert_eq!(
            remeshed::simulate_with_affine_flows_order_cx(
                cx,
                &input,
                0.0,
                3,
                0.1,
                CompositionOrder::Fourth,
                schedule,
                &flow[..1],
                caps(),
            ),
            Err(RemeshedVpmError::Evolution(VpmError::LengthMismatch {
                field: "interval affine flows",
                expected: 2,
                actual: 1
            }))
        );
    });
    for quota in [
        exact.max_work_units - 1,
        exact.max_work_units - reference.source_particle_updates,
    ] {
        let refused = with_cx(
            &CancelGate::new_clock_free(),
            Budget::INFINITE.with_cost_quota(quota),
            |cx| run(cx, exact),
        );
        assert!(matches!(
            refused,
            Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
                BudgetRefusal::CostPlanExceedsQuota { .. }
            )))
        ));
    }
    assert_eq!(input, particles());
}

#[test]
fn g4_late_signed_map_and_support_failure_return_no_partial_state_and_retry_exactly() {
    let input = particles();
    let grids = grids();
    let flow = flows();
    let huge = AffineFlow2::new([0.0; 2], [0.0; 2], [1e308, 0.0, 0.0]).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let run = |grids: &[CartesianGrid2], flows: &[AffineFlow2]| {
            remeshed::simulate_with_affine_flows_order_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                CompositionOrder::Fourth,
                RemeshSchedule {
                    every_steps: 2,
                    grids,
                },
                flows,
                caps(),
            )
        };
        let reference = run(&grids, &flow).unwrap();
        assert_eq!(
            run(&grids, &[flow[0], huge]),
            Err(RemeshedVpmError::Evolution(VpmError::NonFiniteResult {
                phase: "fs-vpm.ambient-map",
                particle: 0,
            }))
        );
        assert!(matches!(
            run(&[grids[0], grids[1], grids[0]], &flow),
            Err(RemeshedVpmError::Remesh(
                RemeshError::SupportOutsideGrid { .. }
            ))
        ));
        assert_eq!(run(&grids, &flow).unwrap(), reference);
    });
    assert_eq!(input, particles());
}

struct TripClock<'a> {
    gate: &'a CancelGate,
    calls: AtomicU64,
    cancel_on: u64,
}
impl TimeSource for TripClock<'_> {
    fn now(&self) -> Time {
        if self.calls.fetch_add(1, Ordering::Relaxed) + 1 == self.cancel_on {
            self.gate.request();
        }
        Time::ZERO
    }
}

#[test]
fn g4_cancellation_inside_bridge_negative_pair_stage_and_final_remap_is_atomic() {
    let input = particles();
    // 144 nodes make the first bridge start at work 144 + 144*143 = 81*256.
    let grids = [
        CartesianGrid2::new([-0.55; 2], [0.1; 2], [12; 2]).unwrap(),
        CartesianGrid2::new([-0.8; 2], [0.1; 2], [17; 2]).unwrap(),
    ];
    let flow = [flows()[0]];
    let budget = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let run = |gate: &CancelGate, clock: &TripClock<'_>, steps| {
        with_cx(gate, budget, |cx| {
            let timed = Cx::new(gate, cx.arena(), cx.stream_key(), budget, cx.mode())
                .with_time_source(clock);
            remeshed::simulate_with_affine_flows_order_cx(
                &timed,
                &input,
                0.001,
                steps,
                0.1,
                CompositionOrder::Fourth,
                RemeshSchedule {
                    every_steps: 1,
                    grids: &grids[..steps + 1],
                },
                &flow[..steps],
                caps(),
            )
        })
    };
    let active = CancelGate::new_clock_free();
    let counter = TripClock {
        gate: &active,
        calls: AtomicU64::new(0),
        cancel_on: u64::MAX,
    };
    run(&active, &counter, 0).unwrap();
    let first_source = counter.calls.swap(0, Ordering::Relaxed);
    let reference = run(&active, &counter, 1).unwrap();
    let polls = counter.calls.load(Ordering::Relaxed);
    // Clock requests at one checkpoint are observed at the next checkpoint.
    // The second probe therefore stops inside the negative-time middle sweep.
    for (cancel_on, phase) in [
        (first_source + 80, "fs-vpm.ambient-source"),
        (first_source + 81, "fs-vpm.conservative-forward"),
        (polls - 4, "fs-vpm.remeshed-remap"),
    ] {
        let gate = CancelGate::new_clock_free();
        let clock = TripClock {
            gate: &gate,
            calls: AtomicU64::new(0),
            cancel_on,
        };
        assert_eq!(
            run(&gate, &clock, 1),
            Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
                BudgetRefusal::Cancelled { phase },
            )))
        );
        assert_eq!(input, particles());
    }
    assert_eq!(run(&active, &counter, 1).unwrap(), reference);
}

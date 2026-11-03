//! G0/G3/G4 checks of the complete native remap–evolution trajectory.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::{
    remeshed::{RemeshSchedule, RemeshedVpmBudget, RemeshedVpmError, simulate_with_cx},
    simulate_viscous_with_cx,
};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::remesh::{CartesianGrid2, RemeshBudget, RemeshError, remesh};
use fs_vpm::{VortexParticle, VpmBudget, VpmError};
use std::sync::atomic::{AtomicU64, Ordering};

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let result = pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 9,
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
        max_steps: 10,
        max_pair_evaluations: 1_000_000,
        max_work_units: 2_000_000,
        max_live_bytes: 64_000,
    }
}

fn remap_caps() -> RemeshBudget {
    RemeshBudget {
        max_input_particles: 1000,
        max_grid_particles: 1000,
        max_work_units: 2_000_000,
        max_live_bytes: 64_000,
    }
}

fn grids() -> [CartesianGrid2; 3] {
    [
        CartesianGrid2::new([-0.3; 2], [0.1; 2], [7; 2]).unwrap(),
        CartesianGrid2::new([-0.5; 2], [0.1; 2], [11; 2]).unwrap(),
        CartesianGrid2::new([-0.7; 2], [0.1; 2], [15; 2]).unwrap(),
    ]
}

fn sources() -> [VortexParticle; 3] {
    [
        VortexParticle::new([-0.08, -0.03], 0.04),
        VortexParticle::new([0.06, 0.05], -0.03),
        VortexParticle::new([0.04, 0.07], 0.01),
    ]
}

fn moments(particles: &[VortexParticle]) -> [f64; 6] {
    let mut result = [0.0; 6];
    for p in particles {
        let [x, y] = p.pos;
        for (out, monomial) in result.iter_mut().zip([1.0, x, y, x * x, x * y, y * y]) {
            *out += p.circulation * monomial;
        }
    }
    result
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-13, "{a:e} != {b:e}");
}

#[test]
fn g0_scheduled_run_matches_original_owners_with_one_complete_work_account() {
    let input = sources();
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let run = simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, caps()).unwrap();
        let mut explicit = remesh(&input, grids[0], remap_caps(), || false).unwrap();
        for (target, steps) in [(grids[1], 2), (grids[2], 1)] {
            let moved = simulate_viscous_with_cx(
                cx,
                &explicit.particles,
                &explicit.areas_m2,
                0.001,
                steps,
                0.1,
                law,
                VpmBudget::new(1000, 10, 1_000_000, 64_000),
            )
            .unwrap();
            explicit = remesh(&moved.particles, target, remap_caps(), || false).unwrap();
        }
        assert_eq!(run.particles, explicit.particles);
        assert_eq!(run.areas_m2, explicit.areas_m2);
        assert_eq!(run.grid, grids[2]);
        assert_eq!(run.steps_completed, 3);
        assert_eq!(run.remeshes_completed, 3);
        assert_eq!(run.pair_evaluations, 3 * (2 * 49 * 48 + 121 * 120));
        assert_eq!(run.stencil_deposits, 16 * (3 + 49 + 121));
        let remap_work = 17 * (3 + 49 + 121) + 2 * (49 + 121 + 225);
        assert_eq!(run.work_units, run.pair_evaluations + 3 + remap_work);
        assert_eq!(run.budget.planned_cost, run.work_units);
        assert_eq!(run.budget.cost_charged, run.work_units);
        assert_eq!(run.peak_live_bytes, 2 * 225 * (24 + 8));
        near(moments(&run.particles)[0], moments(&input)[0]);

        let exact = RemeshedVpmBudget {
            max_live_bytes: run.peak_live_bytes,
            max_work_units: run.work_units,
            max_pair_evaluations: run.pair_evaluations,
            ..caps()
        };
        let bounded = with_cx(
            &CancelGate::new_clock_free(),
            Budget::INFINITE.with_cost_quota(run.work_units),
            |cx| simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, exact),
        )
        .unwrap();
        assert_eq!(bounded.particles, run.particles);
        assert_eq!(bounded.budget.cost_charged, run.work_units);
    });
}

#[test]
fn g0_zero_net_inviscid_cloud_retains_linear_and_angular_impulse_through_remaps() {
    let input = [
        VortexParticle::new([-0.07, 0.04], 0.08),
        VortexParticle::new([0.06, -0.03], -0.08),
    ];
    let before = moments(&input);
    let grids = grids();
    let run = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        simulate_with_cx(
            cx,
            &input,
            0.004,
            4,
            0.1,
            PseDiffusion::new(0.0, 0.16).unwrap(),
            RemeshSchedule {
                every_steps: 2,
                grids: &grids,
            },
            caps(),
        )
    })
    .unwrap();
    let after = moments(&run.particles);
    for axis in 0..3 {
        near(after[axis], before[axis]);
    }
    near(after[3] + after[5], before[3] + before[5]);
    assert_eq!(run.pair_evaluations, 2 * (49 * 48 + 121 * 120));
    assert!(run.particles.iter().any(|p| p.circulation < 0.0));
    assert!(run.particles.iter().any(|p| p.circulation > 0.0));
}

#[test]
fn g3_zero_time_resolution_changes_keep_moments_and_largest_bank_admission() {
    let input = sources();
    let grids = [
        CartesianGrid2::new([-0.5; 2], [0.1; 2], [11; 2]).unwrap(),
        CartesianGrid2::new([-0.8; 2], [0.2; 2], [9; 2]).unwrap(),
        CartesianGrid2::new([-0.8; 2], [0.2; 2], [9; 2]).unwrap(),
    ];
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let run = simulate_with_cx(
            cx,
            &input,
            0.0,
            3,
            0.1,
            law,
            RemeshSchedule {
                every_steps: 2,
                grids: &grids,
            },
            caps(),
        )
        .unwrap();
        for (actual, expected) in moments(&run.particles).into_iter().zip(moments(&input)) {
            near(actual, expected);
        }
        assert_eq!(run.pair_evaluations, 0);
        assert_eq!(run.areas_m2, vec![0.2 * 0.2; 81]);
        assert_eq!(run.peak_live_bytes, 2 * 121 * 32);
        assert!(run.particles.capacity() >= 121 && run.areas_m2.capacity() >= 121);
        let initial_only = simulate_with_cx(
            cx,
            &input,
            0.1,
            0,
            0.1,
            law,
            RemeshSchedule {
                every_steps: 2,
                grids: &grids[..1],
            },
            caps(),
        )
        .unwrap();
        let original = remesh(&input, grids[0], remap_caps(), || false).unwrap();
        assert_eq!(initial_only.particles, original.particles);
        assert_eq!(initial_only.areas_m2, original.areas_m2);
        assert_eq!(initial_only.peak_live_bytes, 121 * 32);
        assert_eq!(initial_only.work_units, original.work_units);
    });
}

#[test]
fn g4_combined_caps_and_ambient_cost_admit_the_whole_schedule() {
    let input = sources();
    let saved = input;
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let reference = simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, caps()).unwrap();
        let mut bad = [caps(); 6];
        bad[0].max_input_particles = input.len() - 1;
        bad[1].max_grid_particles = 224;
        bad[2].max_steps = 2;
        bad[3].max_pair_evaluations = reference.pair_evaluations - 1;
        bad[4].max_work_units = reference.work_units - 1;
        bad[5].max_live_bytes = reference.peak_live_bytes - 1;
        for cap in bad {
            assert!(matches!(
                simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, cap),
                Err(RemeshedVpmError::Evolution(VpmError::CapExceeded { .. }))
            ));
        }
        let refused = with_cx(
            &CancelGate::new_clock_free(),
            Budget::INFINITE.with_cost_quota(reference.work_units - 1),
            |cx| simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, caps()),
        );
        assert!(matches!(
            refused,
            Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
                BudgetRefusal::CostPlanExceedsQuota { .. }
            )))
        ));
    });
    assert_eq!(input, saved);
}

#[test]
fn g4_invalid_schedule_and_late_escaping_support_publish_no_cloud() {
    let input = sources();
    let saved = input;
    let grids = grids();
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for schedule in [
            RemeshSchedule {
                every_steps: 0,
                grids: &grids,
            },
            RemeshSchedule {
                every_steps: 2,
                grids: &grids[..2],
            },
        ] {
            assert!(simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, caps()).is_err());
        }
        let schedule = RemeshSchedule {
            every_steps: 2,
            grids: &grids,
        };
        assert!(simulate_with_cx(cx, &input, -0.001, 3, 0.1, law, schedule, caps()).is_err());
        let insufficient = [grids[0], grids[1], grids[0]];
        assert!(matches!(
            simulate_with_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                law,
                RemeshSchedule {
                    every_steps: 2,
                    grids: &insufficient
                },
                caps(),
            ),
            Err(RemeshedVpmError::Remesh(
                RemeshError::SupportOutsideGrid { .. }
            ))
        ));
    });
    assert_eq!(input, saved);
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
fn g4_cancellation_inside_the_final_remap_discards_all_prior_stages() {
    let input = sources();
    let saved = input;
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    let budget = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let gate = CancelGate::new_clock_free();
    let clock = TripClock {
        gate: &gate,
        calls: AtomicU64::new(0),
        cancel_on: u64::MAX,
    };
    let reference = with_cx(&gate, budget, |cx| {
        let timed =
            Cx::new(&gate, cx.arena(), cx.stream_key(), budget, cx.mode()).with_time_source(&clock);
        simulate_with_cx(&timed, &input, 0.001, 3, 0.1, law, schedule, caps())
    })
    .unwrap();
    let observations = clock.calls.load(Ordering::Relaxed);
    assert!(observations > 20);
    let stop_gate = CancelGate::new_clock_free();
    let stop_clock = TripClock {
        gate: &stop_gate,
        calls: AtomicU64::new(0),
        cancel_on: observations - 4,
    };
    let stopped = with_cx(&stop_gate, budget, |cx| {
        let timed = Cx::new(&stop_gate, cx.arena(), cx.stream_key(), budget, cx.mode())
            .with_time_source(&stop_clock);
        simulate_with_cx(&timed, &input, 0.001, 3, 0.1, law, schedule, caps())
    });
    assert_eq!(
        stopped,
        Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
            BudgetRefusal::Cancelled {
                phase: "fs-vpm.remeshed-remap"
            }
        )))
    );
    assert_eq!(input, saved);
    let retry = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        simulate_with_cx(cx, &input, 0.001, 3, 0.1, law, schedule, caps())
    })
    .unwrap();
    assert_eq!(retry.particles, reference.particles);
    assert_eq!(retry.areas_m2, reference.areas_m2);
}

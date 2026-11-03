//! G0/G3/G4 checks of one source–PSE–advection–remap trajectory.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::{
    ambient::{self, AffineFlow2, AmbientMomentChange},
    remeshed::{
        RemeshSchedule, RemeshedVpmBudget, RemeshedVpmError, simulate_with_affine_flows_cx,
        simulate_with_cx,
    },
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

// Each target changes the window, spacing and particle count. Padding admits
// the previous grid's PSE-filled boundary without clipping its signed tails.
fn grids() -> [CartesianGrid2; 3] {
    [
        CartesianGrid2::new([-0.3; 2], [0.1; 2], [7; 2]).unwrap(),
        CartesianGrid2::new([-0.5, -0.45], [0.08, 0.09], [14, 13]).unwrap(),
        CartesianGrid2::new([-0.8; 2], [0.1, 0.11], [18, 17]).unwrap(),
    ]
}

fn particles() -> [VortexParticle; 3] {
    [
        VortexParticle::new([-0.08, -0.03], 0.04),
        VortexParticle::new([0.06, 0.05], -0.03),
        VortexParticle::new([0.04, 0.07], 0.01),
    ]
}

fn flows() -> [AffineFlow2; 2] {
    [
        AffineFlow2::new([0.01, -0.03], [0.02, -0.01], [0.0, 0.4, 0.0]).unwrap(),
        AffineFlow2::new([-0.02, 0.04], [-0.01, 0.02], [0.05, -0.35, 0.25]).unwrap(),
    ]
}

fn moments(particles: &[VortexParticle]) -> [f64; 4] {
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
fn g0_changing_sources_and_grids_match_original_owners_and_complete_counts() {
    let input = particles();
    let grids = grids();
    let flows = flows();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for viscosity in [0.0, 0.004] {
            let law = PseDiffusion::new(viscosity, 0.16).unwrap();
            let combined = simulate_with_affine_flows_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                law,
                schedule,
                &flows,
                caps(),
            )
            .unwrap();
            let mut explicit = remesh(&input, grids[0], remap_caps(), || false).unwrap();
            let mut source = [0.0; 3];
            for (k, steps) in [2, 1].into_iter().enumerate() {
                let moved = ambient::simulate_viscous_with_cx(
                    cx,
                    &explicit.particles,
                    &explicit.areas_m2,
                    0.001,
                    steps,
                    0.1,
                    flows[k],
                    law,
                    VpmBudget::new(1000, 10, 1_000_000, 64_000),
                )
                .unwrap();
                for (out, change) in source.iter_mut().zip([
                    moved.ambient_moment_change.first_moment_m3_s[0],
                    moved.ambient_moment_change.first_moment_m3_s[1],
                    moved.ambient_moment_change.quadratic_moment_m4_s,
                ]) {
                    *out += change;
                }
                explicit =
                    remesh(&moved.run.particles, grids[k + 1], remap_caps(), || false).unwrap();
            }
            let run = &combined.run;
            assert_eq!(run.particles, explicit.particles);
            assert_eq!(run.areas_m2, explicit.areas_m2);
            assert_eq!(run.grid, grids[2]);
            assert!(run.areas_m2.iter().all(|&area| area == 0.1 * 0.11));
            assert_eq!(run.steps_completed, 3);
            assert_eq!(run.remeshes_completed, 3);
            let sweeps = if viscosity > 0.0 { 3 } else { 1 };
            assert_eq!(run.pair_evaluations, sweeps * (2 * 49 * 48 + 182 * 181));
            assert_eq!(run.stencil_deposits, 16 * (3 + 49 + 182));
            let remap_work = 17 * (3 + 49 + 182) + 2 * (49 + 182 + 306);
            assert_eq!(combined.source_particle_updates, 2 * (2 * 49 + 182));
            assert_eq!(
                run.work_units,
                run.pair_evaluations + 3 + remap_work + combined.source_particle_updates
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
            near(moments(&run.particles)[0], moments(&input)[0]);
            if viscosity > 0.0 {
                // The source receipt must not absorb the distinct PSE defect.
                let total = moments(&run.particles)[3] - moments(&input)[3];
                assert!((total - source[2]).abs() > 1e-9);
            }
        }
    });
}

#[test]
fn g0_independent_piecewise_moment_law_includes_signed_zero_net_clouds() {
    let dipole = [
        VortexParticle::new([-0.07, 0.04], 0.08),
        VortexParticle::new([0.06, -0.03], -0.08),
    ];
    let nonzero = particles();
    let grids = grids();
    let shear = AffineFlow2::new([0.02, -0.01], [0.2, -0.1], [0.0, 0.7, 0.0]).unwrap();
    let rotation = AffineFlow2::new([-0.03, 0.02], [0.0; 2], [0.0, -0.8, 0.8]).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for input in [dipole.as_slice(), nonzero.as_slice()] {
            let before = moments(input);
            let run = simulate_with_affine_flows_cx(
                cx,
                input,
                0.015,
                3,
                0.1,
                PseDiffusion::new(0.0, 0.16).unwrap(),
                RemeshSchedule {
                    every_steps: 2,
                    grids: &grids,
                },
                &[shear, rotation],
                caps(),
            )
            .unwrap();
            // Integrate M'=A(M-Gamma*origin)+Gamma*U analytically for each
            // held interval, independently of the per-particle source code.
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

#[test]
fn g0_inactive_sources_preserve_original_receipts_and_grid_only_paths() {
    let input = particles();
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let huge = AffineFlow2::new([0.0; 2], [0.0; 2], [1e308, 0.0, 0.0]).unwrap();
        for (dt, sources) in [(0.001, [AffineFlow2::default(); 2]), (0.0, [huge; 2])] {
            let original = simulate_with_cx(cx, &input, dt, 3, 0.1, law, schedule, caps()).unwrap();
            let extended = simulate_with_affine_flows_cx(
                cx,
                &input,
                dt,
                3,
                0.1,
                law,
                schedule,
                &sources,
                caps(),
            )
            .unwrap();
            assert_eq!(extended.run, original);
            assert_eq!(extended.source_particle_updates, 0);
            assert_eq!(
                extended.ambient_moment_change,
                AmbientMomentChange::default()
            );
        }
        let initial_only = RemeshSchedule {
            every_steps: 2,
            grids: &grids[..1],
        };
        let original =
            simulate_with_cx(cx, &input, 0.001, 0, 0.1, law, initial_only, caps()).unwrap();
        let extended = simulate_with_affine_flows_cx(
            cx,
            &input,
            0.001,
            0,
            0.1,
            law,
            initial_only,
            &[],
            caps(),
        )
        .unwrap();
        assert_eq!(extended.run, original);
        assert_eq!(extended.run.peak_live_bytes, 49 * (24 + 8));
        assert_eq!(extended.source_particle_updates, 0);

        // Empty vorticity still has authored quadrature nodes to transport.
        let empty =
            simulate_with_affine_flows_cx(cx, &[], 0.001, 3, 0.1, law, schedule, &flows(), caps())
                .unwrap();
        assert!(empty.run.particles.iter().all(|p| p.circulation == 0.0));
        assert_eq!(empty.source_particle_updates, 2 * (2 * 49 + 182));
        assert_eq!(empty.ambient_moment_change, AmbientMomentChange::default());

        let mixed = simulate_with_affine_flows_cx(
            cx,
            &input,
            0.001,
            3,
            0.1,
            law,
            schedule,
            &[flows()[0], AffineFlow2::default()],
            caps(),
        )
        .unwrap();
        assert_eq!(mixed.source_particle_updates, 2 * 2 * 49);
    });
}

#[test]
fn g3_whole_source_cost_and_two_bank_payload_are_admitted_before_evolution() {
    let input = particles();
    let saved = input;
    let grids = grids();
    let schedule = RemeshSchedule {
        every_steps: 2,
        grids: &grids,
    };
    let flows = flows();
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    let reference = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        simulate_with_affine_flows_cx(cx, &input, 0.001, 3, 0.1, law, schedule, &flows, caps())
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
        Budget::INFINITE.with_cost_quota(reference.run.work_units),
        |cx| simulate_with_affine_flows_cx(cx, &input, 0.001, 3, 0.1, law, schedule, &flows, exact),
    )
    .unwrap();
    assert_eq!(admitted.run.particles, reference.run.particles);
    assert_eq!(admitted.run.budget.cost_charged, reference.run.work_units);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for cap in [
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
                simulate_with_affine_flows_cx(
                    cx, &input, 0.001, 3, 0.1, law, schedule, &flows, cap,
                ),
                Err(RemeshedVpmError::Evolution(VpmError::CapExceeded { .. }))
            ));
        }
        for dt in [0.0, 0.001] {
            assert_eq!(
                simulate_with_affine_flows_cx(
                    cx,
                    &input,
                    dt,
                    3,
                    0.1,
                    law,
                    schedule,
                    &flows[..1],
                    caps(),
                ),
                Err(RemeshedVpmError::Evolution(VpmError::LengthMismatch {
                    field: "interval affine flows",
                    expected: 2,
                    actual: 1,
                }))
            );
        }
    });
    let missing_source_cost = reference.run.work_units - reference.source_particle_updates;
    let refused = with_cx(
        &CancelGate::new_clock_free(),
        Budget::INFINITE.with_cost_quota(missing_source_cost),
        |cx| {
            simulate_with_affine_flows_cx(cx, &input, 0.001, 3, 0.1, law, schedule, &flows, caps())
        },
    );
    assert!(matches!(
        refused,
        Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
            BudgetRefusal::CostPlanExceedsQuota { .. }
        )))
    ));
    assert_eq!(input, saved);
}

#[test]
fn g4_late_source_map_and_final_support_failure_publish_no_trajectory() {
    let input = particles();
    let saved = input;
    let grids = grids();
    let good = flows();
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    let huge = AffineFlow2::new([0.0; 2], [0.0; 2], [1e308, 0.0, 0.0]).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let schedule = RemeshSchedule {
            every_steps: 2,
            grids: &grids,
        };
        let reference =
            simulate_with_affine_flows_cx(cx, &input, 0.001, 3, 0.1, law, schedule, &good, caps())
                .unwrap();
        assert_eq!(
            simulate_with_affine_flows_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                law,
                schedule,
                &[good[0], huge],
                caps(),
            ),
            Err(RemeshedVpmError::Evolution(VpmError::NonFiniteResult {
                phase: "fs-vpm.ambient-map",
                particle: 0,
            }))
        );
        let insufficient = [grids[0], grids[1], grids[0]];
        assert!(matches!(
            simulate_with_affine_flows_cx(
                cx,
                &input,
                0.001,
                3,
                0.1,
                law,
                RemeshSchedule {
                    every_steps: 2,
                    grids: &insufficient,
                },
                &good,
                caps(),
            ),
            Err(RemeshedVpmError::Remesh(
                RemeshError::SupportOutsideGrid { .. }
            ))
        ));
        let retry =
            simulate_with_affine_flows_cx(cx, &input, 0.001, 3, 0.1, law, schedule, &good, caps())
                .unwrap();
        assert_eq!(retry, reference);
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
fn g4_cancellation_during_source_and_final_remap_is_atomic_and_retryable() {
    let input = particles();
    let saved = input;
    let grids = [
        CartesianGrid2::new([-0.8; 2], [0.1; 2], [17; 2]).unwrap(),
        CartesianGrid2::new([-1.0; 2], [0.1; 2], [21; 2]).unwrap(),
    ];
    let law = PseDiffusion::new(0.004, 0.16).unwrap();
    let flow = [flows()[0]];
    let schedule = RemeshSchedule {
        every_steps: 1,
        grids: &grids,
    };
    let budget = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let count_gate = CancelGate::new_clock_free();
    let count_clock = TripClock {
        gate: &count_gate,
        calls: AtomicU64::new(0),
        cancel_on: u64::MAX,
    };
    // The first source checkpoint replaces initial-only final publication.
    // Request cancellation there, then the next source tile must observe it.
    with_cx(&count_gate, budget, |cx| {
        let timed = Cx::new(&count_gate, cx.arena(), cx.stream_key(), budget, cx.mode())
            .with_time_source(&count_clock);
        simulate_with_affine_flows_cx(
            &timed,
            &input,
            0.001,
            0,
            0.1,
            law,
            RemeshSchedule {
                every_steps: 1,
                grids: &grids[..1],
            },
            &[],
            caps(),
        )
        .unwrap();
    });
    let source_checkpoint = count_clock.calls.swap(0, Ordering::Relaxed);
    let reference = with_cx(&count_gate, budget, |cx| {
        let timed = Cx::new(&count_gate, cx.arena(), cx.stream_key(), budget, cx.mode())
            .with_time_source(&count_clock);
        simulate_with_affine_flows_cx(&timed, &input, 0.001, 1, 0.1, law, schedule, &flow, caps())
    })
    .unwrap();
    let observations = count_clock.calls.load(Ordering::Relaxed);
    for (cancel_on, phase) in [
        (source_checkpoint, "fs-vpm.ambient-source"),
        (observations - 4, "fs-vpm.remeshed-remap"),
    ] {
        let gate = CancelGate::new_clock_free();
        let clock = TripClock {
            gate: &gate,
            calls: AtomicU64::new(0),
            cancel_on,
        };
        let result = with_cx(&gate, budget, |cx| {
            let timed = Cx::new(&gate, cx.arena(), cx.stream_key(), budget, cx.mode())
                .with_time_source(&clock);
            simulate_with_affine_flows_cx(
                &timed,
                &input,
                0.001,
                1,
                0.1,
                law,
                schedule,
                &flow,
                caps(),
            )
        });
        assert_eq!(
            result,
            Err(RemeshedVpmError::Evolution(VpmError::ExecutionBudget(
                BudgetRefusal::Cancelled { phase }
            )))
        );
        assert_eq!(input, saved);
    }
    let retry = with_cx(&count_gate, budget, |cx| {
        let timed = Cx::new(&count_gate, cx.arena(), cx.stream_key(), budget, cx.mode())
            .with_time_source(&count_clock);
        simulate_with_affine_flows_cx(&timed, &input, 0.001, 1, 0.1, law, schedule, &flow, caps())
    })
    .unwrap();
    assert_eq!(retry, reference);
}

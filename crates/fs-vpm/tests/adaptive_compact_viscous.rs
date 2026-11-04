//! Dense-path parity, moving support and bounded compact-neighbor work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::conservative::adaptive::{
    AdaptiveBudget, AdaptiveError, CompactAdaptiveRun, CompactNeighborBudget,
    ViscousAdaptiveOptions, simulate_viscous_compact_with_cx, simulate_viscous_with_cx,
};
use fs_vpm::diffusion::{PseDiffusion, adaptive::AdaptiveOptions};
use fs_vpm::{VortexParticle, VpmError, total_circulation};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5350_4152_5345,
                kernel_id: 0x5650_4d,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn options() -> ViscousAdaptiveOptions {
    ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: 1.0,
            min_dt_s: 1e-8,
            max_dt_s: 1.0,
            absolute_vorticity_tolerance_s_inv: 1e-8,
            relative_tolerance: 0.0,
        },
        absolute_position_tolerance_m: 1e-8,
    }
}

fn budget(n: usize) -> AdaptiveBudget {
    AdaptiveBudget {
        max_particles: n,
        max_attempts: 10_000,
        max_pair_evaluations: 10_000_000,
        max_work_units: 100_000_000,
        max_live_bytes: 8 * 1024 * 1024,
    }
}

fn neighbors(max_edges: usize) -> CompactNeighborBudget {
    CompactNeighborBudget {
        max_edges,
        max_candidate_visits: 10_000_000,
    }
}

fn run(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas: &[f64],
    duration: f64,
    law: PseDiffusion,
    caps: AdaptiveBudget,
    graph: CompactNeighborBudget,
) -> Result<CompactAdaptiveRun, AdaptiveError> {
    simulate_viscous_compact_with_cx(
        cx,
        particles,
        areas,
        duration,
        0.1,
        law,
        options(),
        caps,
        graph,
    )
}

fn assert_bits(a: &[VortexParticle], b: &[VortexParticle]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

fn crossing_cloud() -> [VortexParticle; 3] {
    [
        VortexParticle::new([-0.5, 0.0], 1.0),
        VortexParticle::new([0.5, 0.0], 1.0),
        VortexParticle::new([0.0, 1.1], 0.0),
    ]
}

#[test]
fn g1_separated_cloud_retains_dense_trajectory_with_far_fewer_diffusion_visits() {
    let initial: Vec<_> = (0..32)
        .flat_map(|i| {
            [
                VortexParticle::new([4.0 * f64::from(i), 0.05], 0.001),
                VortexParticle::new([4.0 * f64::from(i), -0.05], -0.001),
            ]
        })
        .collect();
    let n = initial.len();
    let areas = vec![0.01; n];
    let law = PseDiffusion::compact_wendland(0.001, 0.3).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let dense =
            simulate_viscous_with_cx(cx, &initial, &areas, 0.05, 0.1, law, options(), budget(n))
                .unwrap();
        let compact = run(cx, &initial, &areas, 0.05, law, budget(n), neighbors(n / 2)).unwrap();
        assert_bits(&compact.run.particles, &dense.particles);
        assert_eq!(compact.run.elapsed_s.to_bits(), dense.elapsed_s.to_bits());
        assert_eq!(compact.run.accepted_steps, dense.accepted_steps);
        assert_eq!(compact.run.rejected_steps, dense.rejected_steps);
        assert_eq!(
            compact.run.max_accepted_error_ratio.to_bits(),
            dense.max_accepted_error_ratio.to_bits()
        );
        let advection = 5 * (n * (n - 1)) as u64 * compact.run.attempts as u64;
        let diffusion = compact.run.pair_evaluations - advection;
        assert_eq!(diffusion, 2 * compact.run.fraction_evaluations);
        assert!(30 * diffusion < dense.pair_evaluations - advection);
        assert_eq!(compact.graph_rebuilds, 6 * compact.run.attempts);
        assert!(compact.candidate_visits > 0);
        // Graph probes/search/sorts are charged in addition to physical work.
        let physical_work = n as u64
            + (4 * n as u64 + 5) * compact.run.attempts as u64
            + compact.run.pair_evaluations
            + compact.run.fraction_evaluations;
        assert!(compact.run.work_units > physical_work);
        assert_eq!(compact.run.budget.cost_charged, compact.run.work_units);
        let replay = run(cx, &initial, &areas, 0.05, law, budget(n), neighbors(n / 2)).unwrap();
        assert_eq!(compact, replay);
    });
}

#[test]
fn g1_moved_particles_enter_compact_support_and_match_dense_rejected_trials() {
    let initial = crossing_cloud();
    let law = PseDiffusion::compact_wendland(1.0, 1.19).unwrap();
    // The third particle starts outside both supports. Its advected position
    // enters support during the window, so freezing the initial graph is wrong.
    assert!((initial[2].pos[0] - initial[1].pos[0]).hypot(initial[2].pos[1]) > 1.19);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let dense =
            simulate_viscous_with_cx(cx, &initial, &[1.0; 3], 1.0, 0.1, law, options(), budget(3))
                .unwrap();
        let compact = run(cx, &initial, &[1.0; 3], 1.0, law, budget(3), neighbors(3)).unwrap();
        assert_bits(&compact.run.particles, &dense.particles);
        assert!(compact.run.particles[2].circulation > 1e-6);
        assert!(compact.run.rejected_steps > 0);
        assert_eq!(compact.run.attempts, dense.attempts);
        assert_eq!(compact.graph_rebuilds, 6 * compact.run.attempts);
        assert!((total_circulation(&compact.run.particles) - 2.0).abs() < 2e-13);
    });
}

#[test]
fn g4_aggregate_caps_late_geometry_refusals_and_cancellation_leave_input_reusable() {
    let initial = crossing_cloud();
    let saved = initial;
    let law = PseDiffusion::compact_wendland(1.0, 1.19).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let good = run(cx, &initial, &[1.0; 3], 1.0, law, budget(3), neighbors(3)).unwrap();
        let exact = AdaptiveBudget {
            max_particles: 3,
            max_attempts: good.run.attempts,
            max_pair_evaluations: good.run.pair_evaluations,
            max_work_units: good.run.work_units,
            max_live_bytes: good.run.peak_live_bytes,
        };
        let graph = CompactNeighborBudget {
            max_edges: 3,
            max_candidate_visits: good.candidate_visits,
        };
        let bounded = run(cx, &initial, &[1.0; 3], 1.0, law, exact, graph).unwrap();
        assert_bits(&good.run.particles, &bounded.run.particles);
        for caps in [
            AdaptiveBudget {
                max_attempts: exact.max_attempts - 1,
                ..exact
            },
            AdaptiveBudget {
                max_pair_evaluations: exact.max_pair_evaluations - 1,
                ..exact
            },
            AdaptiveBudget {
                max_work_units: exact.max_work_units - 1,
                ..exact
            },
            AdaptiveBudget {
                max_live_bytes: exact.max_live_bytes - 1,
                ..exact
            },
        ] {
            assert!(matches!(
                run(cx, &initial, &[1.0; 3], 1.0, law, caps, graph),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
            ));
        }
        for graph in [
            CompactNeighborBudget {
                max_candidate_visits: graph.max_candidate_visits - 1,
                ..graph
            },
            // One edge fits initially; the newly supported moved edge must refuse.
            CompactNeighborBudget {
                max_edges: 1,
                ..graph
            },
        ] {
            assert!(matches!(
                run(cx, &initial, &[1.0; 3], 1.0, law, budget(3), graph),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
            ));
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, Budget::INFINITE, |cancelled_cx| {
            assert!(matches!(
                run(cancelled_cx, &initial, &[1.0; 3], 1.0, law, exact, graph),
                Err(AdaptiveError::Engine(VpmError::ExecutionBudget(
                    BudgetRefusal::Cancelled { .. }
                )))
            ));
        });
        with_cx(&gate, Budget::INFINITE.with_poll_quota(7), |limited_cx| {
            assert!(matches!(
                run(limited_cx, &initial, &[1.0; 3], 1.0, law, exact, graph),
                Err(AdaptiveError::Engine(VpmError::ExecutionBudget(
                    BudgetRefusal::PollsExhausted { .. }
                )))
            ));
        });
        assert_bits(&initial, &saved);
        let retry = run(cx, &initial, &[1.0; 3], 1.0, law, exact, graph).unwrap();
        assert_bits(&good.run.particles, &retry.run.particles);
    });
}

#[test]
fn g0_identity_and_zero_viscosity_do_not_allocate_or_build_neighbor_graphs() {
    let mut initial = crossing_cloud();
    initial[0].pos[0] = -0.0;
    let graph = CompactNeighborBudget {
        max_edges: 0,
        max_candidate_visits: 0,
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (n, nu, duration) in [(0, 1.0, 1.0), (1, 1.0, 1.0), (3, 1.0, 0.0), (3, 0.0, 1.0)] {
            let particles = &initial[..n];
            let law = PseDiffusion::compact_wendland(nu, 1.19).unwrap();
            let dense = simulate_viscous_with_cx(
                cx,
                particles,
                &vec![1.0; n],
                duration,
                0.1,
                law,
                options(),
                budget(n),
            )
            .unwrap();
            let compact = run(
                cx,
                particles,
                &vec![1.0; n],
                duration,
                law,
                budget(n),
                graph,
            )
            .unwrap();
            assert_eq!(compact.run, dense);
            assert_bits(&compact.run.particles, &dense.particles);
            assert_eq!(compact.graph_rebuilds, 0);
            assert_eq!(compact.candidate_visits, 0);
        }
        assert!(matches!(
            run(
                cx,
                &initial,
                &[1.0; 3],
                0.0,
                PseDiffusion::new(1.0, 1.0).unwrap(),
                budget(3),
                graph
            ),
            Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
        ));
        assert!(matches!(
            run(
                cx,
                &initial,
                &[1.0, 0.0, 1.0],
                0.0,
                PseDiffusion::compact_wendland(1.0, 1.19).unwrap(),
                budget(3),
                graph
            ),
            Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
        ));
    });
}

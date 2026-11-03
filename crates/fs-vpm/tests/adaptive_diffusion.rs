//! Independent finite-cloud eigenmodes for adaptive compact PSE.
//! These exercise temporal error control; no continuum/spatial-error bound is
//! inferred from the finite graph or the step-doubling acceptance estimate.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::adaptive::{AdaptiveBudget, AdaptiveError, AdaptiveOptions, AdaptiveRun};
use fs_vpm::diffusion::{PreparationBudget, PreparedPse, PseDiffusion};
use fs_vpm::{VortexParticle, VpmError, total_circulation};
use std::f64::consts::PI;
use std::mem::size_of;

const STREAM: StreamKey = StreamKey {
    seed: 0x4144_4150_5449_5645,
    kernel_id: 0x5053_45,
    tile: 0,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx =
            Cx::new(gate, arena, STREAM, budget, ExecMode::Deterministic).with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn options(tolerance: f64) -> AdaptiveOptions {
    AdaptiveOptions {
        initial_dt_s: 1.0,
        min_dt_s: 1e-8,
        max_dt_s: 1.0,
        absolute_vorticity_tolerance_s_inv: tolerance,
        relative_tolerance: 0.0,
    }
}

fn budget(n: usize) -> AdaptiveBudget {
    AdaptiveBudget {
        max_particles: n,
        max_attempts: 10_000,
        max_pair_evaluations: 1_000_000,
        max_work_units: 10_000_000,
        max_live_bytes: 1024 * 1024,
    }
}

fn prepare(cx: &Cx<'_>, particles: &[VortexParticle], areas: &[f64], nu: f64) -> PreparedPse {
    PreparedPse::prepare_with_cx(
        cx,
        particles,
        areas,
        PseDiffusion::compact_wendland(nu, 2.0).unwrap(),
        PreparationBudget::new(particles.len(), 100, 10_000, 100_000, 1024 * 1024),
    )
    .unwrap()
}

fn chain() -> [VortexParticle; 3] {
    [
        VortexParticle::new([0.0, 0.0], 2.0),
        VortexParticle::new([1.0, 0.0], -1.0),
        VortexParticle::new([2.0, 0.0], 0.5),
    ]
}

fn assert_bits(left: &[VortexParticle], right: &[VortexParticle]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_eq!(left.pos.map(f64::to_bits), right.pos.map(f64::to_bits));
        assert_eq!(left.circulation.to_bits(), right.circulation.to_bits());
    }
}

fn exact_chain(time: f64) -> [f64; 3] {
    // Nearest-neighbor q=1/2 gives shape 3/16, R^4=16. The path
    // Laplacian has independent eigenvectors (1,0,-1), (1,-2,1).
    let rate = 1008.0 / (5.0 * PI) / 16.0 * (3.0 / 16.0);
    let slow = 0.75 * (-rate * time).exp();
    let fast = 0.75 * (-3.0 * rate * time).exp();
    [0.5 + slow + fast, 0.5 - 2.0 * fast, 0.5 - slow + fast]
}

fn chain_error(run: &AdaptiveRun) -> f64 {
    run.particles
        .iter()
        .zip(exact_chain(run.elapsed_s))
        .map(|(actual, exact)| (actual.circulation - exact).abs())
        .fold(0.0, f64::max)
}

#[test]
fn g1_tolerance_refinement_resolves_two_independent_diffusion_eigenmodes() {
    let initial = chain();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let graph = prepare(cx, &initial, &[1.0; 3], 1.0);
        assert_eq!(graph.edges(), &[[0, 1], [1, 2]]);
        let loose = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, options(1e-3), budget(3))
            .unwrap();
        let tight = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, options(1e-6), budget(3))
            .unwrap();
        let coarse_error = chain_error(&loose);
        let fine_error = chain_error(&tight);
        assert!(fine_error < 2e-5, "fine error={fine_error}");
        assert!(
            coarse_error > 20.0 * fine_error,
            "{coarse_error} -> {fine_error}"
        );
        assert!(tight.accepted_steps > loose.accepted_steps);
        assert!(tight.rejected_steps > 0);
        for run in [&loose, &tight] {
            assert_eq!(run.elapsed_s.to_bits(), 2.5_f64.to_bits());
            assert_eq!(run.attempts, run.accepted_steps + run.rejected_steps);
            assert!((0.0..=1.0).contains(&run.max_accepted_error_ratio));
            assert!((total_circulation(&run.particles) - 1.5).abs() < 1e-13);
            let final_enstrophy: f64 = run.particles.iter().map(|p| p.circulation.powi(2)).sum();
            assert!(final_enstrophy <= 5.25 + 1e-13);
            assert!(
                run.particles
                    .iter()
                    .all(|p| (-1.0..=2.0).contains(&p.circulation))
            );
            assert_eq!(run.pair_evaluations, 12 * run.attempts as u64);
            assert_eq!(run.fraction_evaluations % 4, 0);
            assert!(run.fraction_evaluations <= 4 * run.attempts as u64);
            assert_eq!(
                run.work_units,
                3 + 19 * run.attempts as u64 + run.fraction_evaluations / 2,
            );
            assert_eq!(run.budget.cost_charged, run.work_units);
            assert_eq!(
                run.peak_live_bytes,
                9 * size_of::<VortexParticle>() + 4 * size_of::<f64>()
            );
        }
        // A stable accepted step size reuses both edge relaxation fractions.
        assert!(tight.fraction_evaluations < 4 * tight.attempts as u64);
        let replay = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, options(1e-6), budget(3))
            .unwrap();
        assert_bits(&tight.particles, &replay.particles);
        assert_eq!(tight.attempts, replay.attempts);
        assert_eq!(tight.work_units, replay.work_units);
        assert_eq!(tight.fraction_evaluations, replay.fraction_evaluations);
    });
}

#[test]
fn g1_unequal_areas_relax_exactly_and_short_terminal_step_reaches_target() {
    let initial = [
        VortexParticle::new([0.0, 0.0], 1.5),
        VortexParticle::new([1.0, 0.0], -0.25),
        VortexParticle::new([6.0, 0.0], -0.0),
    ];
    let areas = [0.25, 0.75, 0.3];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let graph = prepare(cx, &initial, &areas, 0.4);
        let settings = AdaptiveOptions {
            initial_dt_s: 0.9,
            min_dt_s: 0.2,
            max_dt_s: 0.9,
            ..options(1e-11)
        };
        let actual = graph
            .advance_adaptive_with_cx(cx, &initial, 1.0, settings, budget(3))
            .unwrap();
        let rate = 1008.0 / (5.0 * PI) * 0.4 / 16.0 * (3.0 / 16.0);
        let equilibrium = 0.25 * 1.25;
        let expected = equilibrium + (1.5 - equilibrium) * (-rate).exp();
        assert!((actual.particles[0].circulation - expected).abs() < 2e-14);
        assert!((actual.particles[1].circulation - (1.25 - expected)).abs() < 2e-14);
        assert_bits(&actual.particles[2..], &initial[2..]);
        assert_eq!(actual.elapsed_s.to_bits(), 1.0_f64.to_bits());
        assert_eq!(actual.accepted_steps, 2);
        assert_eq!(actual.rejected_steps, 0);
    });
}

#[test]
fn g0_inactive_diffusion_preserves_bits_without_trials() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        for (n, nu, duration) in [(0, 1.0, 2.0), (1, 1.0, 2.0), (3, 0.0, 2.0), (3, 1.0, 0.0)] {
            let mut initial = chain()[..n].to_vec();
            if let Some(first) = initial.first_mut() {
                first.circulation = -0.0;
            }
            let graph = prepare(cx, &initial, &vec![1.0; n], nu);
            let actual = graph
                .advance_adaptive_with_cx(cx, &initial, duration, options(1e-6), budget(n))
                .unwrap();
            assert_bits(&actual.particles, &initial);
            assert_eq!(actual.elapsed_s.to_bits(), duration.to_bits());
            assert_eq!(actual.attempts, 0);
            assert_eq!(actual.pair_evaluations, 0);
            assert_eq!(actual.fraction_evaluations, 0);
            assert_eq!(actual.work_units, n as u64);
            assert_eq!(actual.peak_live_bytes, n * size_of::<VortexParticle>());
        }
    });
}

#[test]
fn g4_failed_trials_caps_and_geometry_leave_inputs_and_graph_reusable() {
    let initial = chain();
    let original = initial;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let graph = prepare(cx, &initial, &[1.0; 3], 1.0);
        let settings = options(1e-6);
        let good = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, settings, budget(3))
            .unwrap();
        let exact = AdaptiveBudget {
            max_particles: 3,
            max_attempts: good.attempts,
            max_pair_evaluations: good.pair_evaluations,
            max_work_units: good.work_units,
            max_live_bytes: good.peak_live_bytes,
        };
        let bounded = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, settings, exact)
            .unwrap();
        assert_bits(&good.particles, &bounded.particles);
        for limited in [
            AdaptiveBudget {
                max_particles: 2,
                ..exact
            },
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
                graph.advance_adaptive_with_cx(cx, &initial, 2.5, settings, limited),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. })),
            ));
        }
        let impossible = AdaptiveOptions {
            min_dt_s: 1.0,
            ..options(1e-12)
        };
        assert!(matches!(
            graph.advance_adaptive_with_cx(cx, &initial, 2.5, impossible, budget(3)),
            Err(AdaptiveError::ToleranceUnreachable { .. }),
        ));
        let mut moved = initial;
        moved[0].pos[0] = -0.0;
        assert!(matches!(
            graph.advance_adaptive_with_cx(cx, &moved, 0.0, settings, budget(3)),
            Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. })),
        ));
        let mut invalid = initial;
        invalid[1].circulation = f64::NAN;
        assert!(matches!(
            graph.advance_adaptive_with_cx(cx, &invalid, 0.0, settings, budget(3)),
            Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. })),
        ));
        assert!(matches!(
            graph.advance_adaptive_with_cx(cx, &initial[..2], 0.0, settings, budget(3)),
            Err(AdaptiveError::Engine(VpmError::LengthMismatch { .. })),
        ));
        for invalid_options in [
            AdaptiveOptions {
                absolute_vorticity_tolerance_s_inv: 0.0,
                ..settings
            },
            AdaptiveOptions {
                relative_tolerance: f64::NAN,
                ..settings
            },
            AdaptiveOptions {
                min_dt_s: 2.0,
                ..settings
            },
            AdaptiveOptions {
                max_dt_s: 0.5,
                ..settings
            },
        ] {
            assert!(matches!(
                graph.advance_adaptive_with_cx(cx, &initial, 0.0, invalid_options, budget(3)),
                Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. })),
            ));
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, Budget::INFINITE, |cancelled_cx| {
            assert!(matches!(
                graph.advance_adaptive_with_cx(cancelled_cx, &initial, 2.5, settings, budget(3)),
                Err(AdaptiveError::Engine(VpmError::ExecutionBudget(
                    BudgetRefusal::Cancelled { .. }
                ))),
            ));
        });
        let retry = graph
            .advance_adaptive_with_cx(cx, &initial, 2.5, settings, budget(3))
            .unwrap();
        assert_bits(&good.particles, &retry.particles);
        assert_bits(&initial, &original);
        assert_eq!(graph.positions()[0][0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(graph.areas_m2(), &[1.0; 3]);
    });
}

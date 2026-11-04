//! Moving-cloud time refinement against an analytic decaying dipole and the
//! original fixed-step viscous owner. These are temporal/finite-cloud checks;
//! neither quadrature refinement nor a continuum error bound is implied.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::conservative::adaptive::{ViscousAdaptiveOptions, simulate_viscous_with_cx};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::diffusion::adaptive::{AdaptiveBudget, AdaptiveError, AdaptiveOptions};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, total_circulation};
use std::f64::consts::PI;
use std::mem::size_of;

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x4144_5645_4354,
                kernel_id: 0x5053_45,
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

fn options(position_tolerance: f64) -> ViscousAdaptiveOptions {
    ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: 1.0,
            min_dt_s: 1e-8,
            max_dt_s: 1.0,
            absolute_vorticity_tolerance_s_inv: 1e-7,
            relative_tolerance: 0.0,
        },
        absolute_position_tolerance_m: position_tolerance,
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

fn dipole() -> [VortexParticle; 2] {
    [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ]
}

fn assert_bits(left: &[VortexParticle], right: &[VortexParticle]) {
    assert_eq!(left.len(), right.len());
    for (a, b) in left.iter().zip(right) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

#[test]
fn g1_position_control_resolves_motion_when_dipole_strength_is_already_exact() {
    let initial = dipole();
    let time = 2.5;
    let nu = 0.3;
    // Equal areas A=1 and fixed separation d=1 give Gamma'= -2*A*kappa*Gamma.
    // Translation integrates Gamma(t)/(2*pi*(d*d+core*core)) in positive x.
    let laws = [
        (
            PseDiffusion::new(nu, 1.0).unwrap(),
            8.0 * nu / (PI * std::f64::consts::E),
        ),
        (
            PseDiffusion::compact_wendland(nu, 2.0).unwrap(),
            2.0 * 1008.0 * nu / (5.0 * PI * 16.0) * (3.0 / 16.0),
        ),
    ];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (law, rate) in laws {
            let loose = simulate_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                time,
                0.1,
                law,
                options(1e-3),
                budget(2),
            )
            .unwrap();
            let tight = simulate_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                time,
                0.1,
                law,
                options(1e-7),
                budget(2),
            )
            .unwrap();
            let exact_gamma = (-rate * time).exp();
            let exact_x = -(-rate * time).exp_m1() / (rate * 2.0 * PI * 1.01);
            let loose_error = (loose.particles[0].pos[0] - exact_x).abs();
            let tight_error = (tight.particles[0].pos[0] - exact_x).abs();
            assert!(tight_error < 3e-6, "position error={tight_error}");
            assert!(
                loose_error > 100.0 * tight_error,
                "{loose_error} -> {tight_error}"
            );
            assert!(tight.accepted_steps > loose.accepted_steps);
            assert!(tight.rejected_steps > 0);
            for run in [&loose, &tight] {
                for (i, p) in run.particles.iter().enumerate() {
                    assert!((p.circulation - initial[i].circulation * exact_gamma).abs() < 1e-13);
                    assert_eq!(p.pos[1], initial[i].pos[1]);
                }
                assert!(total_circulation(&run.particles).abs() < 1e-13);
                assert!((run.particles[0].pos[0] - run.particles[1].pos[0]).abs() < 1e-13);
                assert_eq!(run.elapsed_s.to_bits(), time.to_bits());
                assert_eq!(run.attempts, run.accepted_steps + run.rejected_steps);
                assert!((0.0..=1.0).contains(&run.max_accepted_error_ratio));
                assert_eq!(run.pair_evaluations, 24 * run.attempts as u64);
                assert_eq!(run.fraction_evaluations, 14 * run.attempts as u64);
                assert_eq!(run.work_units, 2 + 37 * run.attempts as u64);
                assert_eq!(run.budget.cost_charged, run.work_units);
                assert_eq!(run.peak_live_bytes, 6 * size_of::<VortexParticle>());
            }
            let replay = simulate_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                time,
                0.1,
                law,
                options(1e-7),
                budget(2),
            )
            .unwrap();
            assert_bits(&tight.particles, &replay.particles);
            assert_eq!(tight.attempts, replay.attempts);
            assert_eq!(tight.work_units, replay.work_units);
        }
    });
}

#[test]
fn g1_stiff_diffusion_must_resolve_transport_before_extinguishing_the_dipole() {
    let initial = dipole();
    let nu = 1000.0;
    let law = PseDiffusion::new(nu, 1.0).unwrap();
    let rate = 8.0 * nu / (PI * std::f64::consts::E);
    let exact_x = -(-rate).exp_m1() / (rate * 2.0 * PI * 1.01);
    // With h=1, both ordinary Strang candidates diffuse to zero before
    // appreciable advection; their difference misses the entire translation.
    // The independent reversed splitting must force initial time refinement.
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let results = [1e-5, 1e-8].map(|tolerance| {
            let actual = simulate_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                1.0,
                0.1,
                law,
                options(tolerance),
                budget(2),
            )
            .unwrap();
            assert_eq!(actual.elapsed_s.to_bits(), 1.0_f64.to_bits());
            assert!(actual.rejected_steps > 0);
            assert!(actual.particles[0].pos[0] > 0.9 * exact_x);
            for (p, original) in actual.particles.iter().zip(initial) {
                assert_eq!(p.pos[1], original.pos[1]);
                assert!(p.circulation.abs() < 1e-13);
            }
            assert!(total_circulation(&actual.particles).abs() < 1e-13);
            assert_eq!(actual.pair_evaluations, 24 * actual.attempts as u64);
            assert_eq!(actual.fraction_evaluations, 14 * actual.attempts as u64);
            assert_eq!(actual.work_units, 2 + 37 * actual.attempts as u64);
            assert_eq!(actual.budget.cost_charged, actual.work_units);
            assert_eq!(actual.peak_live_bytes, 6 * size_of::<VortexParticle>());
            (
                (actual.particles[0].pos[0] - exact_x).abs(),
                actual.accepted_steps,
            )
        });
        let [(loose_error, loose_steps), (tight_error, tight_steps)] = results;
        assert!(tight_error < 3e-8, "stiff translation error={tight_error}");
        assert!(
            loose_error > 50.0 * tight_error,
            "{loose_error} -> {tight_error}"
        );
        assert!(tight_steps > loose_steps);
    });
}

#[test]
fn g1_deforming_unequal_area_cloud_matches_refined_original_viscous_owner() {
    let initial = [
        VortexParticle::new([-0.7, -0.1], 1.0),
        VortexParticle::new([0.5, -0.3], -0.35),
        VortexParticle::new([0.2, 0.8], -0.2),
    ];
    let areas = [0.3, 0.7, 1.4];
    let law = PseDiffusion::new(0.4, 1.0).unwrap();
    let enstrophy = |particles: &[VortexParticle]| {
        particles
            .iter()
            .zip(areas)
            .map(|(p, a)| p.circulation.powi(2) / a)
            .sum::<f64>()
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let settings = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                absolute_vorticity_tolerance_s_inv: 1e-8,
                ..options(1e-8).time
            },
            ..options(1e-8)
        };
        let actual =
            simulate_viscous_with_cx(cx, &initial, &areas, 0.5, 0.1, law, settings, budget(3))
                .unwrap();
        let steps = 4096;
        let reference = fs_vpm::conservative::simulate_viscous_with_cx(
            cx,
            &initial,
            &areas,
            0.5 / steps as f64,
            steps,
            0.1,
            law,
            VpmBudget::new(3, steps, u64::MAX, 3 * size_of::<VortexParticle>()),
        )
        .unwrap();
        for (p, r) in actual.particles.iter().zip(&reference.particles) {
            assert!((p.pos[0] - r.pos[0]).hypot(p.pos[1] - r.pos[1]) < 1e-7);
            assert!((p.circulation - r.circulation).abs() < 1e-7);
        }
        assert!((actual.particles[1].pos[1] - initial[1].pos[1]).abs() > 0.01);
        assert!((total_circulation(&actual.particles) - total_circulation(&initial)).abs() < 2e-13);
        assert!(enstrophy(&actual.particles) < enstrophy(&initial));
        for (p, area) in actual.particles.iter().zip(areas) {
            assert!((-0.5 - 2e-13..=1.0 / 0.3 + 2e-13).contains(&(p.circulation / area)));
        }
        assert_eq!(areas, [0.3, 0.7, 1.4]);
    });
}

#[test]
fn g0_zero_viscosity_keeps_advection_and_short_final_steps_reach_target() {
    let initial = dipole();
    let law = PseDiffusion::new(0.0, 1.0).unwrap();
    let settings = ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: 0.9,
            min_dt_s: 0.2,
            max_dt_s: 0.9,
            ..options(1e-8).time
        },
        ..options(1e-8)
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let actual =
            simulate_viscous_with_cx(cx, &initial, &[1.0; 2], 1.0, 0.1, law, settings, budget(2))
                .unwrap();
        for (p, original) in actual.particles.iter().zip(initial) {
            assert!((p.pos[0] - 1.0 / (2.0 * PI * 1.01)).abs() < 2e-14);
            assert_eq!(p.pos[1], original.pos[1]);
            assert_eq!(p.circulation.to_bits(), original.circulation.to_bits());
        }
        assert_eq!(actual.elapsed_s.to_bits(), 1.0_f64.to_bits());
        assert_eq!(actual.accepted_steps, 2);
        assert_eq!(actual.rejected_steps, 0);
        assert_eq!(actual.pair_evaluations, 6 * actual.attempts as u64);
        assert_eq!(actual.fraction_evaluations, 0);
        assert_eq!(actual.work_units, 2 + 14 * actual.attempts as u64);
        assert_eq!(actual.peak_live_bytes, 6 * size_of::<VortexParticle>());
        for (n, duration) in [(0, 2.0), (1, 2.0), (2, 0.0)] {
            let mut initial = dipole()[..n].to_vec();
            if let Some(p) = initial.first_mut() {
                p.circulation = -0.0;
                p.pos[0] = -0.0;
            }
            let identity = simulate_viscous_with_cx(
                cx,
                &initial,
                &vec![1.0; n],
                duration,
                0.1,
                PseDiffusion::new(0.3, 1.0).unwrap(),
                options(1e-7),
                budget(n),
            )
            .unwrap();
            assert_bits(&identity.particles, &initial);
            assert_eq!(identity.elapsed_s.to_bits(), duration.to_bits());
            assert_eq!(identity.attempts, 0);
            assert_eq!(identity.pair_evaluations, 0);
            assert_eq!(identity.fraction_evaluations, 0);
            assert_eq!(identity.work_units, n as u64);
            assert_eq!(identity.peak_live_bytes, n * size_of::<VortexParticle>());
        }
    });
}

#[test]
fn g4_caps_position_floor_invalid_inputs_and_cancellation_publish_no_partial_run() {
    let initial = dipole();
    let saved = initial;
    let law = PseDiffusion::new(0.3, 1.0).unwrap();
    let settings = options(1e-7);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let run = |cx, settings, caps| {
            simulate_viscous_with_cx(cx, &initial, &[1.0; 2], 2.5, 0.1, law, settings, caps)
        };
        let good = run(cx, settings, budget(2)).unwrap();
        let exact = AdaptiveBudget {
            max_particles: 2,
            max_attempts: good.attempts,
            max_pair_evaluations: good.pair_evaluations,
            max_work_units: good.work_units,
            max_live_bytes: good.peak_live_bytes,
        };
        assert_bits(
            &good.particles,
            &run(cx, settings, exact).unwrap().particles,
        );
        for limited in [
            AdaptiveBudget {
                max_particles: 1,
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
                run(cx, settings, limited),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
            ));
        }
        let impossible = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                min_dt_s: 1.0,
                ..settings.time
            },
            absolute_position_tolerance_m: 1e-12,
        };
        assert!(matches!(
            run(cx, impossible, budget(2)),
            Err(AdaptiveError::ToleranceUnreachable { .. })
        ));
        for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let invalid = ViscousAdaptiveOptions {
                absolute_position_tolerance_m: tolerance,
                ..settings
            };
            assert!(matches!(
                run(cx, invalid, budget(2)),
                Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
            ));
        }
        for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                simulate_viscous_with_cx(
                    cx,
                    &initial,
                    &[1.0, area],
                    0.0,
                    0.1,
                    law,
                    settings,
                    budget(2)
                ),
                Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
            ));
        }
        assert!(matches!(
            simulate_viscous_with_cx(cx, &initial, &[1.0], 0.0, 0.1, law, settings, budget(2)),
            Err(AdaptiveError::Engine(VpmError::LengthMismatch { .. }))
        ));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, Budget::INFINITE, |cancelled_cx| {
            assert!(matches!(
                simulate_viscous_with_cx(
                    cancelled_cx,
                    &initial,
                    &[1.0; 2],
                    2.5,
                    0.1,
                    law,
                    settings,
                    budget(2)
                ),
                Err(AdaptiveError::Engine(VpmError::ExecutionBudget(
                    BudgetRefusal::Cancelled { .. }
                )))
            ));
        });
        with_cx(&gate, Budget::INFINITE.with_poll_quota(4), |limited_cx| {
            assert!(matches!(
                simulate_viscous_with_cx(
                    limited_cx,
                    &initial,
                    &[1.0; 2],
                    2.5,
                    0.1,
                    law,
                    settings,
                    budget(2)
                ),
                Err(AdaptiveError::Engine(VpmError::ExecutionBudget(
                    BudgetRefusal::PollsExhausted { .. }
                )))
            ));
        });
        assert_bits(&initial, &saved);
        assert_bits(
            &good.particles,
            &run(cx, settings, exact).unwrap().particles,
        );
    });
}

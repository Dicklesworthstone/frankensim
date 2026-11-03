//! G0/G1/G3/G4 checks of original PSE plus original conservative advection.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::conservative::{simulate_viscous_with_cx, simulate_with_cx};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::{VortexParticle, VpmBudget, VpmError, VpmRun};
use std::f64::consts::PI;
use std::mem::size_of;

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5053_455f_5041_4952,
                kernel_id: 0x0050_5345,
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

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, n * size_of::<VortexParticle>())
}

fn evolve(
    particles: &[VortexParticle],
    areas: &[f64],
    dt: f64,
    steps: usize,
    core: f64,
    law: PseDiffusion,
) -> VpmRun {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        simulate_viscous_with_cx(
            cx,
            particles,
            areas,
            dt,
            steps,
            core,
            law,
            caps(particles.len(), steps),
        )
    })
    .unwrap()
}

#[test]
fn g0_zero_viscosity_keeps_original_advection_bits_and_resources() {
    let initial = [
        VortexParticle::new([-0.7, -0.1], 1.0),
        VortexParticle::new([0.5, -0.3], -0.5),
        VortexParticle::new([0.2, 0.8], -0.5),
    ];
    let law = PseDiffusion::new(0.0, 0.7).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (dt, steps) in [(0.03, 20), (0.0, 3), (0.1, 0)] {
            let original = simulate_with_cx(cx, &initial, dt, steps, 0.1, caps(3, steps)).unwrap();
            let composed = simulate_viscous_with_cx(
                cx,
                &initial,
                &[0.1, 0.7, 2.0],
                dt,
                steps,
                0.1,
                law,
                caps(3, steps),
            )
            .unwrap();
            for (a, b) in original.particles.iter().zip(&composed.particles) {
                assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
                assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
            }
            assert_eq!(original, composed);
        }
    });
}

#[test]
fn g1_viscous_dipole_decay_and_travel_follow_independent_time_integrals() {
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let nu = 0.3;
    let law = PseDiffusion::new(nu, 1.0).unwrap();
    let rate = 8.0 * nu / (PI * std::f64::consts::E);
    let time = 2.0;
    let strength = (-rate * time).exp();
    let travel = -(-rate * time).exp_m1() / (rate * 2.0 * PI * 1.01);
    let errors = [10, 20, 40].map(|steps| {
        let result = evolve(&initial, &[1.0; 2], time / steps as f64, steps, 0.1, law);
        assert!((result.particles[0].circulation - strength).abs() < 5e-14);
        assert!((result.particles[1].circulation + strength).abs() < 5e-14);
        assert_eq!(result.particles[0].pos[1], 0.5);
        assert_eq!(result.particles[1].pos[1], -0.5);
        assert!((result.particles[0].pos[0] - result.particles[1].pos[0]).abs() < 2e-14);
        assert_eq!(result.pair_evaluations, 6 * steps as u64);
        assert_eq!(result.work_units, 2 + 7 * steps as u64);
        assert_eq!(result.budget.cost_charged, result.work_units);
        assert_eq!(result.peak_live_bytes, 2 * size_of::<VortexParticle>());
        (result.particles[0].pos[0] - travel).abs()
    });
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((3.8..4.2).contains(&ratio), "errors={errors:?}");
    }
    assert!(errors[2] < 3e-6);
}

#[test]
fn g0_signed_cloud_retains_pse_circulation_bounds_and_enstrophy() {
    let areas = [0.03, 0.4, 1.7, 0.1];
    let initial = [
        VortexParticle::new([-0.7, 0.1], -3.0 * areas[0]),
        VortexParticle::new([0.2, 0.6], 4.0 * areas[1]),
        VortexParticle::new([0.8, -0.4], -0.5 * areas[2]),
        VortexParticle::new([-0.3, -0.5], 2.0 * areas[3]),
    ];
    let circulation = |p: &[VortexParticle]| p.iter().map(|p| p.circulation).sum::<f64>();
    let enstrophy = |p: &[VortexParticle]| {
        p.iter()
            .zip(areas)
            .map(|(p, a)| p.circulation * p.circulation / a)
            .sum::<f64>()
    };
    let moment = |p: &[VortexParticle]| p.iter().map(|p| p.circulation * p.pos[0]).sum::<f64>();
    let law = PseDiffusion::new(0.8, 1.0).unwrap();
    let result = evolve(&initial, &areas, 0.05, 4, 0.1, law);
    assert!((circulation(&result.particles) - circulation(&initial)).abs() < 2e-13);
    assert!(enstrophy(&result.particles) < enstrophy(&initial));
    for (p, area) in result.particles.iter().zip(areas) {
        assert!((-3.0 - 2e-13..=4.0 + 2e-13).contains(&(p.circulation / area)));
    }
    // PSE on this finite, unequal cloud does change its first moment. The
    // composition must retain the actual physical/quadrature evolution rather
    // than projecting it back to a falsely inviscid constant-impulse result.
    assert!((moment(&result.particles) - moment(&initial)).abs() > 1e-5);
    assert_eq!(result, evolve(&initial, &areas, 0.05, 4, 0.1, law));
}

#[test]
fn g4_zero_motion_still_admits_areas_and_late_refusals_publish_nothing() {
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let law = PseDiffusion::new(0.3, 1.0).unwrap();
    let active = CancelGate::new_clock_free();
    with_cx(&active, Budget::INFINITE, |cx| {
        let limited = VpmBudget::new(2, 1, 5, 2 * size_of::<VortexParticle>());
        assert!(matches!(
            simulate_viscous_with_cx(cx, &initial, &[1.0; 2], 0.1, 1, 0.1, law, limited),
            Err(VpmError::CapExceeded { .. })
        ));
        for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                simulate_viscous_with_cx(
                    cx,
                    &initial,
                    &[1.0, area],
                    0.0,
                    0,
                    0.1,
                    PseDiffusion::new(0.0, 1.0).unwrap(),
                    caps(2, 0)
                ),
                Err(VpmError::InvalidScalar {
                    particle: Some(1),
                    ..
                })
            ));
        }
        assert!(matches!(
            simulate_viscous_with_cx(cx, &initial, &[1.0], 0.1, 1, 0.1, law, caps(2, 1)),
            Err(VpmError::LengthMismatch { .. })
        ));
        assert!(matches!(
            simulate_viscous_with_cx(cx, &initial, &[1.0; 2], -0.1, 1, 0.1, law, caps(2, 1)),
            Err(VpmError::InvalidScalar { .. })
        ));
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        with_cx(&cancelled, Budget::INFINITE, |cx| {
            simulate_viscous_with_cx(cx, &initial, &[1.0; 2], 0.1, 1, 0.1, law, caps(2, 1))
        }),
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. }))
    ));

    let cloud: Vec<_> = (0..17)
        .map(|i| VortexParticle::new([0.1 * f64::from(i), 0.0], f64::from(i) - 8.0))
        .collect();
    let saved = cloud.clone();
    let refused = with_cx(&active, Budget::INFINITE.with_poll_quota(4), |cx| {
        simulate_viscous_with_cx(cx, &cloud, &[0.1; 17], 0.1, 1, 0.1, law, caps(17, 1))
    });
    assert_eq!(
        refused,
        Err(VpmError::ExecutionBudget(BudgetRefusal::PollsExhausted {
            phase: "fs-vpm.diffusion-pair",
            quota: 4,
        }))
    );
    assert_eq!(cloud, saved);
    assert_eq!(initial[0].pos, [0.0, 0.5]);
}

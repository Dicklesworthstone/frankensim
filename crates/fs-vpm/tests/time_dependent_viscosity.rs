//! Physical exposure, chronology, source ownership and atomic native evolution.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_material::fluid_transport::FluidTransportState;
use fs_vpm::conservative::{
    TimeDependentVpmError, TimeDependentVpmRun, simulate_time_dependent_viscous_with_cx,
    simulate_viscous_with_cx,
};
use fs_vpm::diffusion::viscosity_history::{
    MaterialPseViscosityHistory, PseViscosityError, PseViscosityHistory, ViscosityHistoryError,
    ViscosityKnot,
};
use fs_vpm::diffusion::{PseDiffusion, PseKernel};
use fs_vpm::{VortexParticle, VpmBudget, VpmError};
use std::f64::consts::{E, PI};
use std::mem::size_of;

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x4e55_4849_5354,
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
fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, n * size_of::<VortexParticle>())
}
fn knot(time_s: f64, nu: f64) -> ViscosityKnot {
    ViscosityKnot {
        time_s,
        kinematic_viscosity_m2_s: nu,
    }
}
fn history(knots: &[ViscosityKnot]) -> PseViscosityHistory<'_> {
    PseViscosityHistory::piecewise_linear(knots, 1.0, PseKernel::Gaussian).unwrap()
}
#[allow(clippy::too_many_arguments)]
fn evolve(
    particles: &[VortexParticle],
    areas: &[f64],
    start: f64,
    end: f64,
    maximum_step: f64,
    core: f64,
    history: PseViscosityHistory<'_>,
) -> TimeDependentVpmRun {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        simulate_time_dependent_viscous_with_cx(
            cx,
            particles,
            areas,
            start,
            end,
            maximum_step,
            core,
            history,
            caps(particles.len(), 1024),
        )
    })
    .unwrap()
}

#[test]
fn autonomous_histories_preserve_original_native_bits_and_accounts() {
    let initial = [
        VortexParticle::new([-0.7, -0.1], 1.0),
        VortexParticle::new([0.5, -0.3], -0.5),
        VortexParticle::new([0.2, 0.8], -0.5),
    ];
    for nu in [0.0, 0.15] {
        for kernel in [PseKernel::Gaussian, PseKernel::CompactWendland] {
            let knots = [knot(0.0, nu), knot(1.0, nu)];
            let history = PseViscosityHistory::piecewise_linear(&knots, 1.5, kernel).unwrap();
            let law = match kernel {
                PseKernel::Gaussian => PseDiffusion::new(nu, 1.5).unwrap(),
                PseKernel::CompactWendland => PseDiffusion::compact_wendland(nu, 1.5).unwrap(),
            };
            with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
                let original = simulate_viscous_with_cx(
                    cx,
                    &initial,
                    &[0.1, 0.7, 2.0],
                    0.125,
                    8,
                    0.1,
                    law,
                    caps(3, 8),
                )
                .unwrap();
                let evolving = simulate_time_dependent_viscous_with_cx(
                    cx,
                    &initial,
                    &[0.1, 0.7, 2.0],
                    0.0,
                    1.0,
                    0.125,
                    0.1,
                    history,
                    caps(3, 8),
                )
                .unwrap();
                assert_eq!(evolving.time_s, [0.0, 1.0]);
                assert_eq!(evolving.run, original);
                for (a, b) in evolving.run.particles.iter().zip(&original.particles) {
                    assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
                    assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
                }
            });
        }
    }
}

#[test]
fn linear_viscosity_dipole_has_exact_exposure_and_second_order_physical_travel() {
    // An equal-area opposite-sign pair keeps unit separation. The independent
    // continuous equations are Gamma'= -c nu(t) Gamma and x'=Gamma/(2pi*1.01).
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let knots = [knot(0.0, 0.05), knot(2.0, 0.85)];
    let history = history(&knots);
    let c = 8.0 / (PI * E);
    let strength = |time: f64| (-c * (0.05 * time + 0.2 * time * time)).exp();
    let one = evolve(&initial, &[1.0; 2], 0.0, 2.0, 2.0, 0.1, history);
    // A full-step-averaged viscosity used on both sides would fail this check:
    // the center advection must see the actual FIRST half's diffusion exposure.
    let one_step_travel = 2.0 * strength(1.0) / (2.0 * PI * 1.01);
    assert!((one.run.particles[0].pos[0] - one_step_travel).abs() < 3e-14);
    assert!((one.run.particles[0].circulation - strength(2.0)).abs() < 3e-14);

    // Independent scalar Simpson quadrature, not another PSE splitting loop.
    let panels = 32768;
    let h = 2.0 / f64::from(panels);
    let integral = h / 3.0
        * (strength(0.0)
            + strength(2.0)
            + (1..panels)
                .map(|i| {
                    let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
                    weight * strength(f64::from(i) * h)
                })
                .sum::<f64>());
    let exact_travel = integral / (2.0 * PI * 1.01);
    let errors = [8, 16, 32].map(|steps| {
        let result = evolve(
            &initial,
            &[1.0; 2],
            0.0,
            2.0,
            2.0 / f64::from(steps),
            0.1,
            history,
        );
        assert_eq!(result.run.steps_completed, steps as usize);
        assert!((result.run.particles[0].circulation - strength(2.0)).abs() < 8e-14);
        assert!((result.run.particles[1].circulation + strength(2.0)).abs() < 8e-14);
        assert_eq!(result.run.particles[0].pos[1], 0.5);
        assert_eq!(result.run.particles[1].pos[1], -0.5);
        assert_eq!(result.run.pair_evaluations, 6 * steps as u64);
        (result.run.particles[0].pos[0] - exact_travel).abs()
    });
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((3.95..4.05).contains(&ratio), "errors={errors:?}");
    }
    assert!(errors[2] < 8e-6);
}

#[test]
fn exact_coefficient_events_match_segmented_calls_and_conserve_signed_inventory() {
    let areas = [0.2, 0.7, 1.3];
    let initial = [
        VortexParticle::new([-0.5, -0.2], -1.2 * areas[0]),
        VortexParticle::new([0.4, -0.1], 3.2 * areas[1]),
        VortexParticle::new([0.1, 0.7], -0.1 * areas[2]),
    ];
    let knots = [
        knot(0.0, 0.0),
        knot(0.25, 0.0),
        knot(0.75, 0.2),
        knot(1.0, 0.1),
    ];
    let history = history(&knots);
    let full = evolve(&initial, &areas, 0.0, 1.0, 0.125, 0.1, history);
    let mut state = initial.to_vec();
    let mut pairs = 0;
    for window in knots.windows(2) {
        let segment = evolve(
            &state,
            &areas,
            window[0].time_s,
            window[1].time_s,
            0.125,
            0.1,
            history,
        );
        assert_eq!(segment.time_s[1], window[1].time_s);
        pairs += segment.run.pair_evaluations;
        state = segment.run.particles;
    }
    assert_eq!(state, full.run.particles);
    assert_eq!(full.run.steps_completed, 8);
    assert_eq!(full.run.pair_evaluations, 120);
    assert_eq!(full.run.pair_evaluations, pairs);
    assert_eq!(full.run.work_units, 131);
    let circulation = |p: &[VortexParticle]| p.iter().map(|p| p.circulation).sum::<f64>();
    let enstrophy = |p: &[VortexParticle]| {
        p.iter()
            .zip(areas)
            .map(|(p, area)| p.circulation * p.circulation / area)
            .sum::<f64>()
    };
    assert!((circulation(&state) - circulation(&initial)).abs() < 5e-14);
    assert!(enstrophy(&state) <= enstrophy(&initial) * (1.0 + 1e-14));
    for (particle, area) in state.iter().zip(areas) {
        assert!((-1.2 - 1e-14..=3.2 + 1e-14).contains(&(particle.circulation / area)));
    }
}

#[test]
fn source_histories_retain_snapshots_and_require_fixed_incompressible_density() {
    let sources = [
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, 100.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, 500.0).unwrap(),
    ];
    let binding = MaterialPseViscosityHistory::piecewise_linear(
        &sources,
        &[0.0, 1.0],
        1.0,
        PseKernel::Gaussian,
    )
    .unwrap();
    assert!(core::ptr::eq(binding.sources().as_ptr(), sources.as_ptr()));
    assert_eq!(binding.transport_states(), &sources);
    assert!(
        (binding
            .history()
            .mean_law(0.0, 1.0)
            .unwrap()
            .kinematic_viscosity_m2_s()
            - 0.3)
            .abs()
            < 1e-15
    );
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let knots = [knot(0.0, 0.1), knot(1.0, 0.5)];
    assert_eq!(
        evolve(&initial, &[1.0; 2], 0.0, 1.0, 0.125, 0.1, binding.history()),
        evolve(&initial, &[1.0; 2], 0.0, 1.0, 0.125, 0.1, history(&knots)),
    );
    let changing_density = [
        sources[0],
        FluidTransportState::try_new(999.0, 4200.0, 0.6, 500.0).unwrap(),
    ];
    assert!(matches!(
        MaterialPseViscosityHistory::piecewise_linear(
            &changing_density,
            &[0.0, 1.0],
            1.0,
            PseKernel::Gaussian
        ),
        Err(PseViscosityError::DensityChanged { sample: 1 }),
    ));
}

#[test]
fn unsupported_windows_and_early_or_late_refusals_publish_no_cloud() {
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let saved = initial;
    let knots = [knot(0.0, 0.1), knot(1.0, 0.5)];
    let history = history(&knots);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        assert!(matches!(
            simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                0.0,
                2.0,
                0.125,
                0.1,
                history,
                caps(2, 32)
            ),
            Err(TimeDependentVpmError::Coefficient(
                PseViscosityError::History(ViscosityHistoryError::OutsideSupport)
            )),
        ));
        assert!(matches!(
            simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                0.0,
                1.0,
                0.125,
                0.1,
                history,
                caps(2, 3)
            ),
            Err(TimeDependentVpmError::Evolution(VpmError::CapExceeded {
                resource: "steps",
                ..
            })),
        ));
        let identity = simulate_time_dependent_viscous_with_cx(
            cx,
            &initial,
            &[1.0; 2],
            0.25,
            0.25,
            0.125,
            0.1,
            history,
            caps(2, 0),
        )
        .unwrap();
        assert_eq!(identity.run.particles, initial);
        assert_eq!(identity.run.steps_completed, 0);
        assert!(matches!(
            simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[1.0, -1.0],
                0.25,
                0.25,
                0.125,
                0.1,
                history,
                caps(2, 0)
            ),
            Err(TimeDependentVpmError::Evolution(VpmError::InvalidScalar {
                particle: Some(1),
                ..
            })),
        ));
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        with_cx(&cancelled, Budget::INFINITE, |cx| {
            simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                0.0,
                1.0,
                0.125,
                0.1,
                history,
                caps(2, 8),
            )
        }),
        Err(TimeDependentVpmError::Evolution(VpmError::ExecutionBudget(
            BudgetRefusal::Cancelled { .. }
        )))
    ));
    let cloud: Vec<_> = (0..17)
        .map(|i| VortexParticle::new([0.1 * f64::from(i), 0.0], f64::from(i) - 8.0))
        .collect();
    let saved_cloud = cloud.clone();
    let refused = with_cx(
        &CancelGate::new_clock_free(),
        Budget::INFINITE.with_poll_quota(4),
        |cx| {
            simulate_time_dependent_viscous_with_cx(
                cx,
                &cloud,
                &[0.1; 17],
                0.0,
                1.0,
                0.125,
                0.1,
                history,
                caps(17, 8),
            )
        },
    );
    assert!(matches!(
        refused,
        Err(TimeDependentVpmError::Evolution(VpmError::ExecutionBudget(
            BudgetRefusal::PollsExhausted { .. }
        )))
    ));
    assert_eq!(cloud, saved_cloud);
    assert_eq!(initial, saved);
}

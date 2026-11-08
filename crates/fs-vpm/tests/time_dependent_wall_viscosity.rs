//! Changing-viscosity wall loss and image-driven travel against independent ODEs.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::viscosity_history::{PseViscosityHistory, ViscosityKnot};
use fs_vpm::diffusion::{PseDiffusion, PseKernel};
use fs_vpm::planar_wall::{self, PlanarWall2, TimeDependentVpmError, TimeDependentVpmRun};
use fs_vpm::{VortexParticle, VpmBudget, VpmError};
use std::f64::consts::PI;
use std::mem::size_of;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5741_4c4c_4e55,
                kernel_id: 0x5053_45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(
        n,
        steps,
        u64::MAX,
        n * (5 * size_of::<VortexParticle>() + 4 * size_of::<[f64; 2]>()),
    )
}

fn wall() -> PlanarWall2 {
    PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap()
}

fn knot(time_s: f64, nu: f64) -> ViscosityKnot {
    ViscosityKnot {
        time_s,
        kinematic_viscosity_m2_s: nu,
    }
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
    with_cx(&CancelGate::new_clock_free(), |cx| {
        planar_wall::simulate_time_dependent_viscous_with_cx(
            cx,
            particles,
            areas,
            start,
            end,
            maximum_step,
            core,
            history,
            wall(),
            caps(particles.len(), 1024),
        )
    })
    .unwrap()
}

#[test]
fn constant_histories_retain_original_wall_particles_and_receipts() {
    let initial = [
        VortexParticle::new([-0.7, 1.1], 1.0),
        VortexParticle::new([0.5, 0.7], -0.5),
        VortexParticle::new([0.2, 1.8], -0.5),
    ];
    for nu in [0.0, 0.15] {
        for kernel in [PseKernel::Gaussian, PseKernel::CompactWendland] {
            let knots = [knot(0.0, nu), knot(1.0, nu)];
            let history = PseViscosityHistory::piecewise_linear(&knots, 1.5, kernel).unwrap();
            let law = match kernel {
                PseKernel::Gaussian => PseDiffusion::new(nu, 1.5).unwrap(),
                PseKernel::CompactWendland => PseDiffusion::compact_wendland(nu, 1.5).unwrap(),
            };
            with_cx(&CancelGate::new_clock_free(), |cx| {
                let original = planar_wall::simulate_viscous_with_cx(
                    cx,
                    &initial,
                    &[0.1, 0.7, 2.0],
                    0.125,
                    8,
                    0.1,
                    law,
                    wall(),
                    caps(3, 8),
                )
                .unwrap();
                let evolving = planar_wall::simulate_time_dependent_viscous_with_cx(
                    cx,
                    &initial,
                    &[0.1, 0.7, 2.0],
                    0.0,
                    1.0,
                    0.125,
                    0.1,
                    history,
                    wall(),
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
fn changing_viscosity_has_exact_wall_decay_and_second_order_image_travel() {
    let height: f64 = 0.3;
    let core: f64 = 0.08;
    let area: f64 = 0.8;
    let width: f64 = 0.9;
    let gamma: f64 = 0.8;
    let initial = [VortexParticle::new([-0.2, height], gamma)];
    let knots = [knot(0.0, 0.05), knot(2.0, 0.85)];
    let history =
        PseViscosityHistory::piecewise_linear(&knots, width, PseKernel::Gaussian).unwrap();
    // Independently, y is constant, Gamma'=-c*nu(t)*Gamma, and
    // x'=y*Gamma/[pi*(4*y^2+core^2)]. Here nu(t)=0.05+0.4*t.
    let c = 8.0 * area / (PI * width.powi(4)) * (-4.0 * height * height / (width * width)).exp();
    let velocity_per_gamma = height / (PI * (4.0 * height * height + core * core));
    let strength = |time: f64| gamma * (-c * (0.05 * time + 0.2 * time * time)).exp();
    let one = evolve(&initial, &[area], 0.0, 2.0, 2.0, core, history);
    // This detects substituting the full-window mean in both diffusion halves.
    let one_step_x = initial[0].pos[0] + 2.0 * velocity_per_gamma * strength(1.0);
    assert!((one.run.particles[0].pos[0] - one_step_x).abs() < 3e-14);
    assert!((one.run.particles[0].circulation - strength(2.0)).abs() < 3e-14);
    assert!(one.run.particles[0].circulation < gamma);

    // Scalar Simpson integration of the independent continuous ODE; this does
    // not call a PSE kernel, source interpolator or split-step owner.
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
    let exact_x = initial[0].pos[0] + velocity_per_gamma * integral;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let run = evolve(
            &initial,
            &[area],
            0.0,
            2.0,
            2.0 / f64::from(steps),
            core,
            history,
        )
        .run;
        assert_eq!(run.particles.len(), 1);
        assert!((run.particles[0].pos[1] - height).abs() < 1e-15);
        assert!((run.particles[0].circulation - strength(2.0)).abs() < 3e-14);
        assert_eq!(run.pair_evaluations, 12 * steps as u64);
        errors.push((run.particles[0].pos[0] - exact_x).abs());
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

#[test]
fn knot_cuts_and_refusals_preserve_the_original_wall_transaction() {
    let initial = [VortexParticle::new([0.0, 0.4], 0.7)];
    let saved = initial;
    let knots = [
        knot(0.0, 0.1),
        knot(0.25, 0.3),
        knot(0.75, 0.0),
        knot(1.0, 0.2),
    ];
    let history = PseViscosityHistory::piecewise_linear(&knots, 1.0, PseKernel::Gaussian).unwrap();
    let whole = evolve(&initial, &[0.6], 0.0, 1.0, 1.0, 0.1, history);
    assert_eq!(whole.run.steps_completed, 3);
    let mut segmented = initial.to_vec();
    for window in knots.windows(2) {
        segmented = evolve(
            &segmented,
            &[0.6],
            window[0].time_s,
            window[1].time_s,
            1.0,
            0.1,
            history,
        )
        .run
        .particles;
    }
    assert_eq!(whole.run.particles, segmented);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            planar_wall::simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[0.6],
                0.0,
                1.0,
                1.0,
                0.1,
                history,
                wall(),
                caps(1, 2),
            ),
            Err(TimeDependentVpmError::Evolution(VpmError::CapExceeded {
                resource: "steps",
                ..
            })),
        ));
        assert!(matches!(
            planar_wall::simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[0.6],
                0.0,
                2.0,
                1.0,
                0.1,
                history,
                wall(),
                caps(1, 8),
            ),
            Err(TimeDependentVpmError::Coefficient(_)),
        ));
        let empty = planar_wall::simulate_time_dependent_viscous_with_cx(
            cx,
            &initial,
            &[0.6],
            0.25,
            0.25,
            1.0,
            0.1,
            history,
            wall(),
            caps(1, 0),
        )
        .unwrap();
        assert_eq!(empty.run.particles, initial);
        assert_eq!(empty.run.peak_live_bytes, size_of::<VortexParticle>());
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        with_cx(&cancelled, |cx| {
            planar_wall::simulate_time_dependent_viscous_with_cx(
                cx,
                &initial,
                &[0.6],
                0.0,
                1.0,
                1.0,
                0.1,
                history,
                wall(),
                caps(1, 3),
            )
        }),
        Err(TimeDependentVpmError::Evolution(VpmError::ExecutionBudget(
            BudgetRefusal::Cancelled { .. }
        ))),
    ));
    assert_eq!(initial, saved);
}

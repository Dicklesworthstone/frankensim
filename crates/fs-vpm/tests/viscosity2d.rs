//! G0/G1/G3/G4 evidence for actual viscous particle trajectories.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::{
    VortexParticle, VpmBudget, VpmError, VpmRun,
    diffusion::{PseDiffusion, simulate_viscous_with_cx},
    simulate_with_cx, total_circulation,
};
use std::f64::consts::PI;

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 24,
                kernel_id: 7,
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

fn law(nu: f64, width: f64) -> PseDiffusion {
    PseDiffusion::new(nu, width).unwrap()
}

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

fn run(
    particles: &[VortexParticle],
    areas: &[f64],
    dt: f64,
    steps: usize,
    law: PseDiffusion,
) -> VpmRun {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        simulate_viscous_with_cx(
            cx,
            particles,
            areas,
            dt,
            steps,
            0.1,
            law,
            caps(particles.len(), steps),
        )
        .unwrap()
    })
}

fn enstrophy(particles: &[VortexParticle], areas: &[f64]) -> f64 {
    particles
        .iter()
        .zip(areas)
        .map(|(p, a)| 0.5 * p.circulation * p.circulation / a)
        .sum()
}

#[test]
fn unequal_area_exchange_matches_independent_two_particle_decay() {
    // Coincident quadrature samples isolate the exchange ODE from advection.
    // Omega_i relaxes toward the area-weighted mean, not the mean circulation.
    let particles = [
        VortexParticle::new([0.0; 2], 0.8),
        VortexParticle::new([0.0; 2], -0.3),
    ];
    let areas = [0.4, 0.6];
    let material = law(0.1, 0.7);
    let t = 0.8;
    let coefficient = 4.0 * 0.1 / (PI * 0.7_f64.powi(4));
    let mean = 0.5;
    let expected = areas[0] * (mean + areas[1] * 2.5 * (-coefficient * t).exp());
    let coarse = run(&particles, &areas, 0.1, 8, material);
    let fine = run(&particles, &areas, 0.05, 16, material);
    let coarse_error = (coarse.particles[0].circulation - expected).abs();
    let fine_error = (fine.particles[0].circulation - expected).abs();
    // The pair maps solve this isolated mode analytically at either time step.
    assert!(fine_error < 1e-14, "fine error {fine_error}");
    assert!(coarse_error < 1e-14, "coarse error {coarse_error}");
    assert!((total_circulation(&fine.particles) - 0.5).abs() < 2e-15);
    assert!(enstrophy(&fine.particles, &areas) < enstrophy(&particles, &areas));

    let uniform = [
        VortexParticle::new([0.0; 2], 0.5),
        VortexParticle::new([1.0, 0.0], 1.0),
    ];
    let stationary_strengths = run(&uniform, &[0.25, 0.5], 0.02, 4, material);
    for (initial, final_) in uniform.iter().zip(&stationary_strengths.particles) {
        assert_eq!(initial.circulation.to_bits(), final_.circulation.to_bits());
    }
}

#[test]
fn coupled_viscosity_conserves_circulation_dissipates_and_changes_the_actual_wake() {
    let areas = [0.25, 0.5, 0.75, 0.375];
    let particles = [
        VortexParticle::new([0.0, 0.0], 0.25),
        VortexParticle::new([0.5, 0.0], -0.25),
        VortexParticle::new([0.0, 0.5], 0.1875),
        VortexParticle::new([0.5, 0.5], 0.0),
    ];
    let viscous = run(&particles, &areas, 0.02, 20, law(0.2, 0.8));
    let inviscid = run(&particles, &areas, 0.02, 20, law(0.0, 0.8));
    assert!((total_circulation(&viscous.particles) - total_circulation(&particles)).abs() < 2e-15);
    assert!(enstrophy(&viscous.particles, &areas) < 0.7 * enstrophy(&particles, &areas));
    for (p, area) in viscous.particles.iter().zip(areas) {
        assert!((-0.5 - 1e-14..=1.0 + 1e-14).contains(&(p.circulation / area)));
    }
    assert!((viscous.particles[0].pos[0] - inviscid.particles[0].pos[0]).abs() > 1e-4);
}

#[test]
fn translation_units_and_permutation_refinement_preserve_the_viscous_physics() {
    let areas = [0.125, 0.5, 0.25];
    let particles = [
        VortexParticle::new([-0.5, 0.0], 0.025),
        VortexParticle::new([0.5, 0.0], -0.015),
        VortexParticle::new([0.0, 0.5], 0.005),
    ];
    let material = law(0.1, 0.8);
    let base = run(&particles, &areas, 0.025, 8, material);
    let translated =
        particles.map(|p| VortexParticle::new([p.pos[0] + 8.0, p.pos[1] - 4.0], p.circulation));
    let shifted = run(&translated, &areas, 0.025, 8, material);
    for (a, b) in base.particles.iter().zip(shifted.particles) {
        assert!((a.pos[0] - (b.pos[0] - 8.0)).abs() < 5e-15);
        assert!((a.pos[1] - (b.pos[1] + 4.0)).abs() < 5e-15);
        assert!((a.circulation - b.circulation).abs() < 2e-16);
    }
    // Symmetric pair splitting is second order, but its finite-step error
    // depends on pair order. Relabeling must converge to the same dynamics.
    let permutation_error = |reference: &VpmRun, dt, steps| {
        let reverse = run(
            &[particles[2], particles[1], particles[0]],
            &[areas[2], areas[1], areas[0]],
            dt,
            steps,
            material,
        );
        reference
            .particles
            .iter()
            .zip(reverse.particles.iter().rev())
            .map(|(a, b)| (a.circulation - b.circulation).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let coarse_error = permutation_error(&base, 0.025, 8);
    let fine = run(&particles, &areas, 0.0125, 16, material);
    let fine_error = permutation_error(&fine, 0.0125, 16);
    assert!(fine_error > 1e-15);
    assert!((3.6..4.4).contains(&(coarse_error / fine_error)));
    // Length unit scaling: x,eps,core -> 2x; A,Gamma,nu -> 4x.
    let scaled_particles =
        particles.map(|p| VortexParticle::new(p.pos.map(|x| 2.0 * x), 4.0 * p.circulation));
    let scaled_areas = areas.map(|a| 4.0 * a);
    let scaled = with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        simulate_viscous_with_cx(
            cx,
            &scaled_particles,
            &scaled_areas,
            0.025,
            8,
            0.2,
            law(0.4, 1.6),
            caps(3, 8),
        )
        .unwrap()
    });
    for (a, b) in base.particles.iter().zip(scaled.particles) {
        assert_eq!(a.pos.map(|x| (2.0 * x).to_bits()), b.pos.map(f64::to_bits));
        assert_eq!((4.0 * a.circulation).to_bits(), b.circulation.to_bits());
    }
}

fn gaussian_error(cells: usize, width: f64) -> f64 {
    // Tiny circulation suppresses quadrature-distortion errors from tangential
    // advection while the actual production flow/diffusion transaction runs.
    let h = 8.0 / (cells - 1) as f64;
    let total_gamma = 1e-6;
    let mut particles = Vec::new();
    for iy in 0..cells {
        for ix in 0..cells {
            let x = -4.0 + ix as f64 * h;
            let y = -4.0 + iy as f64 * h;
            let omega = total_gamma / PI * (-(x * x + y * y)).exp();
            particles.push(VortexParticle::new([x, y], omega * h * h));
        }
    }
    let areas = vec![h * h; particles.len()];
    let result = run(&particles, &areas, 0.05, 2, law(0.2, width));
    let a2 = 1.0 + 4.0 * 0.2 * 0.1;
    let mut difference = 0.0;
    let mut reference = 0.0;
    for p in &result.particles {
        let [x, y] = p.pos;
        let exact = total_gamma / (PI * a2) * (-(x * x + y * y) / a2).exp();
        difference += (p.circulation / (h * h) - exact).powi(2);
        reference += exact * exact;
    }
    assert!((total_circulation(&result.particles) - total_circulation(&particles)).abs() < 2e-20);
    assert!(enstrophy(&result.particles, &areas) < enstrophy(&particles, &areas));
    (difference / reference).sqrt()
}

#[test]
fn gaussian_vortex_diffusion_converges_to_lamb_oseen_width() {
    let coarse = gaussian_error(13, 1.2);
    let medium = gaussian_error(25, 0.7);
    let fine = gaussian_error(41, 0.5);
    eprintln!("Gaussian relative L2 errors: {coarse:.6e}, {medium:.6e}, {fine:.6e}");
    assert!(medium < 0.6 * coarse);
    assert!(fine < 0.65 * medium);
    assert!(fine < 0.015);
}

#[test]
fn zero_viscosity_keeps_legacy_bits_and_active_costs_match_the_actual_work() {
    let particles = [
        VortexParticle::new([-0.5, 0.0], 0.1),
        VortexParticle::new([0.5, 0.0], -0.1),
    ];
    let legacy = with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        simulate_with_cx(cx, &particles, 0.02, 3, 0.1, caps(2, 3)).unwrap()
    });
    let zero = run(&particles, &[0.5, 0.5], 0.02, 3, law(0.0, 0.8));
    for (a, b) in legacy.particles.iter().zip(&zero.particles) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
    assert_eq!(legacy.work_units, zero.work_units);
    assert_eq!(legacy.peak_live_bytes, zero.peak_live_bytes);
    assert_eq!(legacy.pair_evaluations, zero.pair_evaluations);
    let active = run(&particles, &[0.5, 0.5], 0.02, 3, law(0.1, 0.8));
    assert_eq!(active.pair_evaluations, 4 * 3 * 2 * 2 + 2 * 3 * 2);
    assert_eq!(active.work_units, 2 + 3 * (6 * 2 * 2 + 2 * 2 + 1));
    assert_eq!(active.work_units, active.budget.cost_charged);
    assert_eq!(active.peak_live_bytes, legacy.peak_live_bytes);
    let empty = run(&[], &[], 0.02, 3, law(0.1, 0.8));
    assert!(empty.particles.is_empty());
    let no_steps = run(&particles, &[0.5, 0.5], 0.02, 0, law(0.1, 0.8));
    assert_eq!(no_steps.particles, particles);
    // A half-weight average would round this odd subnormal circulation.
    let tiny = particles.map(|p| VortexParticle::new(p.pos, f64::from_bits(3)));
    let no_time = run(&tiny, &[0.5, 0.5], 0.0, 2, law(0.1, 0.8));
    assert_eq!(no_time.particles, tiny);
    assert_eq!(no_time.pair_evaluations, 4 * 2 * 2 * 2);
    let constant_tiny = run(&tiny, &[0.5, 0.5], 0.02, 2, law(0.1, 0.8));
    assert_eq!(constant_tiny.particles, tiny);
    assert_eq!(constant_tiny.pair_evaluations, 4 * 2 * 2 * 2 + 2 * 2 * 2);
}

#[test]
fn invalid_inputs_budget_and_cancellation_refuse_without_publishing_partial_particles() {
    let particles = [
        VortexParticle::new([0.0; 2], 0.1),
        VortexParticle::new([0.0; 2], -0.2),
    ];
    let before = particles;
    let areas = [0.5, 1.0];
    let material = law(0.1, 0.5);
    let attempt = |dt, area: &[f64], law, caps, gate: &CancelGate, budget| {
        with_cx(gate, budget, |cx| {
            simulate_viscous_with_cx(cx, &particles, area, dt, 2, 0.1, law, caps)
        })
    };
    let gate = CancelGate::new();
    // Exact pair relaxation remains dissipative above any explicit CFL limit.
    let stiff = attempt(2.0, &areas, material, caps(2, 2), &gate, Budget::INFINITE).unwrap();
    assert!(enstrophy(&stiff.particles, &areas) < enstrophy(&particles, &areas));
    assert!(matches!(
        attempt(-0.1, &areas, material, caps(2, 2), &gate, Budget::INFINITE),
        Err(VpmError::InvalidScalar { .. })
    ));
    assert!(matches!(
        attempt(
            0.01,
            &areas[..1],
            material,
            caps(2, 2),
            &gate,
            Budget::INFINITE
        ),
        Err(VpmError::LengthMismatch { .. })
    ));
    assert!(matches!(
        attempt(
            0.01,
            &[0.5, 0.0],
            material,
            caps(2, 2),
            &gate,
            Budget::INFINITE
        ),
        Err(VpmError::InvalidScalar { .. })
    ));
    assert!(matches!(
        PseDiffusion::new(-0.1, 0.5),
        Err(VpmError::InvalidScalar { .. })
    ));
    let limited = VpmBudget::new(2, 2, 4 * 2 * 2 * 2, usize::MAX); // Inviscid-only pair allowance.
    assert!(matches!(
        attempt(0.01, &areas, material, limited, &gate, Budget::INFINITE),
        Err(VpmError::CapExceeded {
            resource: "pair evaluations",
            ..
        })
    ));
    assert!(matches!(
        attempt(
            0.01,
            &areas,
            material,
            caps(2, 2),
            &gate,
            Budget::INFINITE.with_poll_quota(6)
        ),
        Err(VpmError::ExecutionBudget(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        attempt(
            0.01,
            &areas,
            material,
            caps(2, 2),
            &cancelled,
            Budget::INFINITE
        ),
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. }))
    ));
    assert_eq!(particles, before);
    assert_eq!(
        run(&particles, &areas, 0.01, 2, material).steps_completed,
        2
    );
}

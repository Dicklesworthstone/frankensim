//! Independent physical checks for 2-D particle-strength diffusion. The oracles
//! are the exact two-site relaxation, three-site Laplacian eigenmodes, the
//! Gaussian heat equation, and a translating viscous dipole. Particle areas
//! are quadrature weights, so bounds and enstrophy concern Gamma / area.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::{PseDiffusion, diffuse_with_cx, simulate_viscous_with_cx};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, VpmRun, simulate_with_cx, total_circulation};
use std::f64::consts::PI;
use std::mem::size_of;

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x5053_455F_4449_4646,
    kernel_id: 0x0050_5345,
    tile: 0,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, TEST_STREAM, budget, ExecMode::Deterministic)
            .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn roomy_budget(particles: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(particles, steps, u64::MAX, usize::MAX)
}

fn diffuse(
    particles: &[VortexParticle],
    areas: &[f64],
    dt: f64,
    steps: usize,
    law: PseDiffusion,
) -> VpmRun {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        diffuse_with_cx(
            cx,
            particles,
            areas,
            dt,
            steps,
            law,
            roomy_budget(particles.len(), steps),
        )
    })
    .expect("finite physical diffusion fixture must admit")
}

fn enstrophy(particles: &[VortexParticle], areas: &[f64]) -> f64 {
    particles
        .iter()
        .zip(areas)
        .map(|(p, area)| p.circulation * p.circulation / area)
        .sum()
}

fn assert_particles_bitwise_equal(left: &[VortexParticle], right: &[VortexParticle]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_eq!(left.pos.map(f64::to_bits), right.pos.map(f64::to_bits));
        assert_eq!(left.circulation.to_bits(), right.circulation.to_bits());
    }
}

#[test]
fn g0_unequal_area_pair_has_exact_relaxation_even_beyond_explicit_stability() {
    let initial = [
        VortexParticle::new([0.0, 0.0], 2.1),
        VortexParticle::new([0.4, -0.2], -0.6),
    ];
    let areas = [0.3, 1.7];
    let nu = 0.7;
    let width: f64 = 0.8;
    let law = PseDiffusion::new(nu, width).unwrap();
    let k = 4.0 * nu / (PI * width.powi(4)) * (-0.2 / width.powi(2)).exp();
    let total = total_circulation(&initial);
    let equilibrium = areas[0] * total / (areas[0] + areas[1]);
    for time in [0.05, 2.5] {
        let actual = diffuse(&initial, &areas, time, 1, law);
        let first = equilibrium
            + (initial[0].circulation - equilibrium) * (-k * (areas[0] + areas[1]) * time).exp();
        assert!((actual.particles[0].circulation - first).abs() < 2e-14);
        assert!((actual.particles[1].circulation - (total - first)).abs() < 2e-14);
        for (before, after) in initial.iter().zip(&actual.particles) {
            assert_eq!(before.pos.map(f64::to_bits), after.pos.map(f64::to_bits));
        }
        assert_eq!(actual.pair_evaluations, 2);
        assert_eq!(actual.work_units, 5);
        assert_eq!(actual.budget.cost_charged, actual.work_units);
        assert_eq!(actual.peak_live_bytes, 2 * size_of::<VortexParticle>());
    }
}

#[test]
fn g0_extreme_area_ratio_retains_representable_transfer_and_uniform_vorticity() {
    let areas = [1e-200, 1e200];
    let initial = [
        VortexParticle::new([0.0, 0.0], 0.0),
        VortexParticle::new([0.0, 0.0], 1e200),
    ];
    let law = PseDiffusion::new(1.0, 1.0).unwrap();
    let actual = diffuse(&initial, &areas, 1e-200, 1, law);
    // A0 / A1 underflows, but its product with Gamma1 is representable.
    // The exact rate-time product is (4/pi) * (A0+A1) * dt = 4/pi
    // to working precision; the small site's equilibrium circulation is A0.
    let expected = 1e-200 * -(-4.0 / PI).exp_m1();
    let received = actual.particles[0].circulation;
    assert!(received.is_finite() && received > 0.0);
    assert!((received / expected - 1.0).abs() < 1e-12);

    let uniform = [
        VortexParticle::new([0.0, 0.0], 1e-200),
        VortexParticle::new([0.0, 0.0], 1e200),
    ];
    let stationary = diffuse(&uniform, &areas, 1e-200, 1, law);
    for (before, after) in uniform.iter().zip(&stationary.particles) {
        assert!(after.circulation.is_finite() && after.circulation > 0.0);
        assert!((after.circulation / before.circulation - 1.0).abs() < 1e-12);
    }
}

#[test]
fn g0_signed_vorticity_conserves_circulation_and_contracts_bounds_and_enstrophy() {
    let positions = [
        [-0.7, 0.1],
        [0.2, 0.6],
        [0.8, -0.4],
        [-0.3, -0.5],
        [0.1, 0.0],
    ];
    let areas = [0.03, 0.4, 1.7, 0.1, 2.5];
    let vorticity = [-3.0, 4.0, -0.5, 2.0, -1.0];
    let mut particles: Vec<_> = positions
        .iter()
        .zip(areas.iter().zip(vorticity))
        .map(|(position, (area, omega))| VortexParticle::new(*position, area * omega))
        .collect();
    let law = PseDiffusion::new(0.4, 1.1).unwrap();
    let initial_total = total_circulation(&particles);
    let initial_enstrophy = enstrophy(&particles, &areas);
    // Deliberately large steps: stability cannot come from a hidden explicit
    // diffusion CFL or from clipping negative vorticity.
    for _ in 0..4 {
        let previous_enstrophy = enstrophy(&particles, &areas);
        particles = diffuse(&particles, &areas, 50.0, 1, law).particles;
        assert!((total_circulation(&particles) - initial_total).abs() < 2e-13);
        assert!(enstrophy(&particles, &areas) <= previous_enstrophy + 2e-13);
        for (p, area) in particles.iter().zip(areas) {
            let omega = p.circulation / area;
            assert!((-3.0 - 2e-14..=4.0 + 2e-14).contains(&omega));
        }
    }
    assert!(enstrophy(&particles, &areas) < 0.2 * initial_enstrophy);

    // Uniform vorticity is stationary even though the circulations differ by
    // almost two orders of magnitude on this unequal-area cloud.
    let uniform: Vec<_> = positions
        .iter()
        .zip(areas)
        .map(|(position, area)| VortexParticle::new(*position, -1.25 * area))
        .collect();
    let stationary = diffuse(&uniform, &areas, 50.0, 4, law);
    for (before, after) in uniform.iter().zip(&stationary.particles) {
        assert!((before.circulation - after.circulation).abs() < 2e-14);
    }
}

fn gaussian_error(half_count: i32) -> f64 {
    let h = 3.0 / f64::from(half_count);
    let area = h * h;
    let nu = 0.2;
    let time = 0.25;
    let initial_width2 = 1.0;
    let final_width2 = initial_width2 + 4.0 * nu * time;
    let gaussian = |position: [f64; 2], width2: f64| {
        (-(position[0] * position[0] + position[1] * position[1]) / width2).exp() / (PI * width2)
    };
    let mut particles = Vec::new();
    for j in -half_count..=half_count {
        for i in -half_count..=half_count {
            let position = [f64::from(i) * h, f64::from(j) * h];
            particles.push(VortexParticle::new(
                position,
                area * gaussian(position, initial_width2),
            ));
        }
    }
    let areas = vec![area; particles.len()];
    let law = PseDiffusion::new(nu, 1.5 * h).unwrap();
    let actual = diffuse(&particles, &areas, time / 8.0, 8, law);
    assert!((total_circulation(&actual.particles) - total_circulation(&particles)).abs() < 2e-12);
    let mut squared_error = 0.0;
    let mut squared_norm = 0.0;
    for particle in &actual.particles {
        let exact = gaussian(particle.pos, final_width2);
        let omega = particle.circulation / area;
        squared_error += area * (omega - exact).powi(2);
        squared_norm += area * exact * exact;
        assert!(omega >= 0.0);
    }
    (squared_error / squared_norm).sqrt()
}

#[test]
fn g1_gaussian_heat_solution_converges_under_cloud_and_kernel_refinement() {
    // The initial unit-circulation Gaussian is a Lamb-Oseen vorticity profile.
    // Its heat-equation width squared grows by 4 nu t. Padding to |x|,|y|=3
    // keeps finite-cloud boundary error below the resolved kernel error here.
    // eps/h=1.5 resolves Gaussian quadrature; both lengths halve together.
    let coarse = gaussian_error(6);
    let fine = gaussian_error(12);
    assert!(fine < 0.015, "relative L2 error {fine}");
    assert!(coarse > 2.8 * fine, "coarse={coarse}, fine={fine}");
    assert!(coarse < 0.05, "relative L2 error {coarse}");
}

#[test]
fn g1_three_site_analytic_modes_show_second_order_time_convergence() {
    let area = 0.7;
    let nu = 0.4;
    let width: f64 = 0.9;
    let distance: f64 = 0.5;
    let time = 0.8;
    let initial = [
        VortexParticle::new([-distance, 0.0], area),
        VortexParticle::new([0.0, 0.0], -0.4 * area),
        VortexParticle::new([distance, 0.0], 0.3 * area),
    ];
    let law = PseDiffusion::new(nu, width).unwrap();
    let factor = 4.0 * nu * area / (PI * width.powi(4));
    let neighbor = factor * (-distance.powi(2) / width.powi(2)).exp();
    let across = factor * (-4.0 * distance.powi(2) / width.powi(2)).exp();
    // The complete three-site operator has eigenvectors (1,1,1), (1,0,-1),
    // and (1,-2,1), with rates 0, neighbor+2*across, and 3*neighbor.
    // This oracle does not repeat any ordered pair sweep.
    let odd = 0.35 * (-(neighbor + 2.0 * across) * time).exp();
    let even = 0.35 * (-3.0 * neighbor * time).exp();
    let exact = [
        area * (0.3 + odd + even),
        area * (0.3 - 2.0 * even),
        area * (0.3 - odd + even),
    ];
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let result = diffuse(&initial, &[area; 3], time / steps as f64, steps, law);
        let error = result
            .particles
            .iter()
            .zip(exact)
            .map(|(particle, exact)| (particle.circulation - exact).powi(2))
            .sum::<f64>()
            .sqrt();
        errors.push(error);
    }
    assert!(
        errors[2] > 1e-9,
        "refinement must exceed roundoff: {errors:?}"
    );
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.7..4.3).contains(&ratio), "errors={errors:?}");
    }
}

#[test]
fn g3_viscous_dynamics_obey_joint_length_time_and_circulation_scaling() {
    let initial = [
        VortexParticle::new([-0.4, 0.1], 0.8),
        VortexParticle::new([0.3, 0.5], -0.2),
        VortexParticle::new([0.7, -0.4], 0.5),
    ];
    let areas = [0.2, 0.7, 0.4];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let reference = simulate_viscous_with_cx(
            cx,
            &initial,
            &areas,
            0.015,
            7,
            0.05,
            PseDiffusion::new(0.1, 0.7).unwrap(),
            roomy_budget(3, 7),
        )
        .unwrap();
        for (length, time) in [(1e3, 60.0), (1e-3, 0.01)] {
            let strength = length * length / time;
            let scaled = initial
                .map(|p| VortexParticle::new(p.pos.map(|x| x * length), p.circulation * strength));
            let result = simulate_viscous_with_cx(
                cx,
                &scaled,
                &areas.map(|a| a * length * length),
                0.015 * time,
                7,
                0.05 * length,
                PseDiffusion::new(0.1 * strength, 0.7 * length).unwrap(),
                roomy_budget(3, 7),
            )
            .unwrap();
            for (expected, actual) in reference.particles.iter().zip(&result.particles) {
                for axis in 0..2 {
                    assert!((actual.pos[axis] / length - expected.pos[axis]).abs() < 2e-13);
                }
                assert!((actual.circulation / strength - expected.circulation).abs() < 2e-13);
            }
        }
    });
}

#[test]
fn g2_diffusing_dipole_changes_the_integrated_trajectory_and_zero_nu_recovers_advection() {
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let areas = [1.0; 2];
    let dt = 0.025;
    let steps = 40;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let inviscid =
            simulate_with_cx(cx, &initial, dt, steps, 0.0, roomy_budget(2, steps)).unwrap();
        let zero = simulate_viscous_with_cx(
            cx,
            &initial,
            &areas,
            dt,
            steps,
            0.0,
            PseDiffusion::new(0.0, 1.0).unwrap(),
            roomy_budget(2, steps),
        )
        .unwrap();
        assert_particles_bitwise_equal(&zero.particles, &inviscid.particles);
        assert_eq!(zero.pair_evaluations, inviscid.pair_evaluations);
        assert_eq!(zero.work_units, inviscid.work_units);

        let nu = 0.6;
        let actual = simulate_viscous_with_cx(
            cx,
            &initial,
            &areas,
            dt,
            steps,
            0.0,
            PseDiffusion::new(nu, 1.0).unwrap(),
            roomy_budget(2, steps),
        )
        .unwrap();
        // The dipole stays one metre apart. Gamma(t)=exp(-lambda*t),
        // and its speed Gamma(t)/(2*pi) integrates analytically.
        let rate = 8.0 * nu / (PI * std::f64::consts::E);
        let time = dt * steps as f64;
        let strength = (-rate * time).exp();
        let travel = -(-rate * time).exp_m1() / (2.0 * PI * rate);
        assert!((actual.particles[0].circulation - strength).abs() < 2e-14);
        assert!((actual.particles[1].circulation + strength).abs() < 2e-14);
        for particle in &actual.particles {
            assert!((particle.pos[0] - travel).abs() < 2e-6);
        }
        assert_eq!(actual.particles[0].pos[1], 0.5);
        assert_eq!(actual.particles[1].pos[1], -0.5);
        assert!(inviscid.particles[0].pos[0] - actual.particles[0].pos[0] > 0.02);
        assert_eq!(
            actual.pair_evaluations - inviscid.pair_evaluations,
            4 * steps as u64
        );
    });
}

#[test]
fn g4_invalid_domains_caps_and_mid_sweep_refusal_publish_no_state() {
    for (nu, width) in [
        (-0.1, 1.0),
        (f64::NAN, 1.0),
        (0.1, 0.0),
        (0.1, f64::INFINITY),
    ] {
        assert!(PseDiffusion::new(nu, width).is_err());
    }
    let initial = [
        VortexParticle::new([-0.5, 0.0], 1.0),
        VortexParticle::new([0.5, 0.0], -0.2),
    ];
    let original = initial;
    let law = PseDiffusion::new(0.1, 1.0).unwrap();
    let active = CancelGate::new_clock_free();
    with_cx(&active, Budget::INFINITE, |cx| {
        assert_eq!(
            diffuse_with_cx(cx, &initial, &[1.0], 0.1, 1, law, roomy_budget(2, 1)),
            Err(VpmError::LengthMismatch {
                field: "particle areas",
                expected: 2,
                actual: 1
            })
        );
        for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let result = diffuse_with_cx(
                cx,
                &initial,
                &[1.0, area],
                0.0,
                0,
                PseDiffusion::new(0.0, 1.0).unwrap(),
                roomy_budget(2, 0),
            );
            assert_eq!(
                result,
                Err(VpmError::InvalidScalar {
                    field: "particle area",
                    bits: area.to_bits(),
                    particle: Some(1),
                })
            );
        }
        assert_eq!(
            diffuse_with_cx(cx, &initial, &[1.0; 2], -0.1, 1, law, roomy_budget(2, 1)),
            Err(VpmError::InvalidScalar {
                field: "diffusion dt",
                bits: (-0.1_f64).to_bits(),
                particle: None,
            })
        );
        assert_eq!(
            diffuse_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                0.1,
                1,
                law,
                VpmBudget::new(2, 1, 1, usize::MAX)
            ),
            Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: 2,
                maximum: 1,
            })
        );
        assert_eq!(
            simulate_viscous_with_cx(
                cx,
                &initial,
                &[1.0; 2],
                0.1,
                1,
                0.0,
                law,
                VpmBudget::new(2, 1, 19, usize::MAX)
            ),
            Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: 20,
                maximum: 19,
            })
        );
    });

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let refusal = with_cx(&cancelled, Budget::INFINITE, |cx| {
        diffuse_with_cx(cx, &initial, &[1.0; 2], 0.1, 1, law, roomy_budget(2, 1))
    });
    assert_eq!(
        refusal,
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled {
            phase: "fs-vpm.diffusion-admission",
        }))
    );

    // Admission, input start/end, and the first pair checkpoint are admitted.
    // The quota is exhausted at the next stride inside the 17*16 pair sweep,
    // after real private numerical updates have happened.
    let cloud: Vec<_> = (0..17)
        .map(|i| VortexParticle::new([0.1 * f64::from(i), 0.0], f64::from(i) - 8.0))
        .collect();
    let saved = cloud.clone();
    let refusal = with_cx(&active, Budget::INFINITE.with_poll_quota(4), |cx| {
        diffuse_with_cx(cx, &cloud, &[0.1; 17], 0.1, 1, law, roomy_budget(17, 1))
    });
    assert_eq!(
        refusal,
        Err(VpmError::ExecutionBudget(BudgetRefusal::PollsExhausted {
            phase: "fs-vpm.diffusion-pair",
            quota: 4,
        }))
    );
    assert_particles_bitwise_equal(&cloud, &saved);
    assert_particles_bitwise_equal(&initial, &original);
}

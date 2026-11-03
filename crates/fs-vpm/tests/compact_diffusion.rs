//! Independent physical regressions for compact particle-strength exchange.
//! Kernel moments are recovered from actual pair relaxation and integrated by
//! quadrature; Gaussian heat evolution and dipole travel have continuum oracles.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::{PseDiffusion, PseKernel, diffuse_with_cx, simulate_viscous_with_cx};
use fs_vpm::{VortexParticle, VpmBudget, VpmRun, simulate_with_cx, total_circulation};
use std::f64::consts::PI;

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x434F_4D50_4143_545F,
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
    .expect("finite compact diffusion fixture must admit")
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

fn gaussian_error(half_count: i32) -> f64 {
    let h = 3.0 / f64::from(half_count);
    let area = h * h;
    let nu = 0.2;
    let time = 0.25;
    let final_width2 = 1.0 + 4.0 * nu * time;
    let gaussian = |position: [f64; 2], width2: f64| {
        (-(position[0] * position[0] + position[1] * position[1]) / width2).exp()
            / (PI * width2)
    };
    let mut particles = Vec::new();
    for j in -half_count..=half_count {
        for i in -half_count..=half_count {
            let position = [f64::from(i) * h, f64::from(j) * h];
            particles.push(VortexParticle::new(position, area * gaussian(position, 1.0)));
        }
    }
    let areas = vec![area; particles.len()];
    let law = PseDiffusion::compact_wendland(nu, 4.0 * h).unwrap();
    let actual = diffuse(&particles, &areas, time / 8.0, 8, law);
    assert!((total_circulation(&actual.particles) - total_circulation(&particles)).abs() < 2e-12);
    let mut squared_error = 0.0;
    let mut squared_norm = 0.0;
    for particle in &actual.particles {
        let exact = gaussian(particle.pos, final_width2);
        let omega = particle.circulation / area;
        assert!(omega >= 0.0);
        squared_error += area * (omega - exact).powi(2);
        squared_norm += area * exact * exact;
    }
    (squared_error / squared_norm).sqrt()
}

#[test]
fn g1_compact_second_moment_and_heat_equation_refinement_have_independent_oracles() {
    let law = PseDiffusion::compact_wendland(1.0, 1.0).unwrap();
    assert_eq!(law.kernel(), PseKernel::CompactWendland);
    assert_eq!(law.support_radius_m(), Some(1.0));
    let gaussian = PseDiffusion::new(1.0, 1.0).unwrap();
    assert_eq!(gaussian.kernel(), PseKernel::Gaussian);
    assert_eq!(gaussian.support_radius_m(), None);

    let gate = CancelGate::new_clock_free();
    let second_moment = with_cx(&gate, Budget::INFINITE, |cx| {
        let intervals = 128;
        let dt = 1e-4;
        let mut radial_integral = 0.0;
        for i in 0..=intervals {
            let radius = f64::from(i) / f64::from(intervals);
            let pair = [
                VortexParticle::new([0.0, 0.0], 0.0),
                VortexParticle::new([radius, 0.0], 1.0),
            ];
            let actual = diffuse_with_cx(
                cx, &pair, &[1.0; 2], dt, 1, law, roomy_budget(2, 1),
            )
            .unwrap();
            // For equal unit areas, Gamma0=(1-exp(-2*k*dt))/2.
            // Recover the implemented radial rate without copying its kernel.
            let rate = -(-2.0 * actual.particles[0].circulation).ln_1p() / (2.0 * dt);
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            radial_integral += weight * radius.powi(3) * rate;
        }
        // Angular integration of x^2 contributes pi, so this is the full
        // Cartesian second moment, not the unweighted radial kernel mass.
        PI * radial_integral / (3.0 * f64::from(intervals))
    });
    assert!((second_moment - 2.0).abs() < 5e-8, "moment={second_moment}");

    // Unit-circulation Gaussian heat evolution broadens width^2 by 4*nu*t.
    // Both cloud spacing and support radius halve; four spacings per radius
    // keep the quadrature error small while exposing finite-kernel error.
    let coarse = gaussian_error(6);
    let fine = gaussian_error(12);
    assert!(coarse < 0.05, "coarse relative L2 error={coarse}");
    assert!(fine < 0.015, "fine relative L2 error={fine}");
    assert!(coarse > 2.5 * fine, "coarse={coarse}, fine={fine}");
}

#[test]
fn g0_compact_pairs_relax_exactly_and_disconnected_signed_clouds_preserve_invariants() {
    let initial = [
        VortexParticle::new([0.0, 0.0], 1.5),
        VortexParticle::new([0.4, 0.0], -0.25),
        VortexParticle::new([3.0, 0.0], -0.7),
    ];
    let areas = [0.25, 0.75, 0.2];
    let nu = 0.07;
    let radius: f64 = 0.8;
    let law = PseDiffusion::compact_wendland(nu, radius).unwrap();
    // At separation radius/2, the Wendland shape is exactly 3/16.
    let pair_rate = 1008.0 / (5.0 * PI) * nu / radius.powi(4) * (3.0 / 16.0);
    let pair_total = initial[0].circulation + initial[1].circulation;
    let equilibrium = areas[0] * pair_total / (areas[0] + areas[1]);
    for time in [0.03, 4.0] {
        let actual = diffuse(&initial, &areas, time, 1, law);
        let expected = equilibrium
            + (initial[0].circulation - equilibrium)
                * (-pair_rate * (areas[0] + areas[1]) * time).exp();
        assert!((actual.particles[0].circulation - expected).abs() < 3e-14);
        assert!((actual.particles[1].circulation - (pair_total - expected)).abs() < 3e-14);
        assert_particles_bitwise_equal(&actual.particles[2..], &initial[2..]);
        assert!((total_circulation(&actual.particles) - total_circulation(&initial)).abs() < 3e-14);
    }

    let positions = [[0.0, 0.0], [0.3, 0.1], [0.5, -0.1], [1.1, 0.0]];
    let areas = [0.03, 0.4, 1.7, 0.1];
    let omega = [-3.0, 4.0, -0.5, 2.0];
    let cloud: Vec<_> = positions
        .iter()
        .zip(areas.iter().zip(omega))
        .map(|(position, (area, value))| VortexParticle::new(*position, area * value))
        .collect();
    let actual = diffuse(&cloud, &areas, 50.0, 3, law);
    assert!((total_circulation(&actual.particles) - total_circulation(&cloud)).abs() < 2e-13);
    assert!(enstrophy(&actual.particles, &areas) < enstrophy(&cloud, &areas));
    for ((before, after), area) in cloud.iter().zip(&actual.particles).zip(areas) {
        assert_eq!(before.pos.map(f64::to_bits), after.pos.map(f64::to_bits));
        assert!((-3.0 - 2e-14..=4.0 + 2e-14).contains(&(after.circulation / area)));
    }
    let uniform: Vec<_> = positions
        .iter()
        .zip(areas)
        .map(|(position, area)| VortexParticle::new(*position, -1.25 * area))
        .collect();
    let stationary = diffuse(&uniform, &areas, 50.0, 3, law);
    for (before, after) in uniform.iter().zip(&stationary.particles) {
        assert!((after.circulation - before.circulation).abs() < 2e-14);
    }
}

#[test]
fn g2_g3_compact_viscosity_changes_dipole_travel_and_respects_physical_scaling() {
    let initial = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    let areas = [1.0; 2];
    let dt = 0.0125;
    let steps = 80;
    let nu = 0.6;
    let radius: f64 = 2.0;
    let law = PseDiffusion::compact_wendland(nu, radius).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let actual = simulate_viscous_with_cx(
            cx, &initial, &areas, dt, steps, 0.0, law, roomy_budget(2, steps),
        )
        .unwrap();
        let rate = 2.0 * 1008.0 / (5.0 * PI) * nu / radius.powi(4) * (3.0 / 16.0);
        let time = dt * steps as f64;
        let strength = (-rate * time).exp();
        let travel = -(-rate * time).exp_m1() / (2.0 * PI * rate);
        assert!((actual.particles[0].circulation - strength).abs() < 4e-14);
        assert!((actual.particles[1].circulation + strength).abs() < 4e-14);
        for (before, after) in initial.iter().zip(&actual.particles) {
            assert_eq!(before.pos[1], after.pos[1]);
            assert!((after.pos[0] - travel).abs() < 1e-6);
        }
        let inviscid = simulate_with_cx(
            cx, &initial, dt, steps, 0.0, roomy_budget(2, steps),
        )
        .unwrap();
        let zero = simulate_viscous_with_cx(
            cx,
            &initial,
            &areas,
            dt,
            steps,
            0.0,
            PseDiffusion::compact_wendland(0.0, radius).unwrap(),
            roomy_budget(2, steps),
        )
        .unwrap();
        assert_particles_bitwise_equal(&zero.particles, &inviscid.particles);
        assert!(inviscid.particles[0].pos[0] - actual.particles[0].pos[0] > 0.03);

        for (length, time_scale) in [(1e3, 60.0), (1e-3, 0.01)] {
            let circulation_scale = length * length / time_scale;
            let scaled = initial.map(|p| {
                VortexParticle::new(
                    p.pos.map(|x| x * length),
                    p.circulation * circulation_scale,
                )
            });
            let result = simulate_viscous_with_cx(
                cx,
                &scaled,
                &areas.map(|area| area * length * length),
                dt * time_scale,
                steps,
                0.0,
                PseDiffusion::compact_wendland(
                    nu * circulation_scale, radius * length,
                )
                .unwrap(),
                roomy_budget(2, steps),
            )
            .unwrap();
            for (expected, value) in actual.particles.iter().zip(&result.particles) {
                for axis in 0..2 {
                    assert!((value.pos[axis] / length - expected.pos[axis]).abs() < 3e-12);
                }
                assert!((value.circulation / circulation_scale - expected.circulation).abs() < 3e-12);
            }
        }
    });
}

#[test]
fn g0_support_edge_keeps_resolvable_weak_exchange_and_exactly_zero_exterior_flux() {
    let inside = f64::from_bits(1.0_f64.to_bits() - 1);
    let outside = f64::from_bits(1.0_f64.to_bits() + 1);
    let initial = [
        VortexParticle::new([0.0, 0.0], 0.0),
        VortexParticle::new([inside, 0.0], 1.0),
    ];
    let law = PseDiffusion::compact_wendland(1e60, 1.0).unwrap();
    let gap = 0.5 * f64::EPSILON;
    let rate = 1008.0 / (5.0 * PI) * 1e60 * gap.powi(4) * (5.0 - 4.0 * gap);
    let expected = 0.5 * -(-2.0 * rate).exp_m1();
    let actual = diffuse(&initial, &[1.0; 2], 1.0, 1, law);
    let received = actual.particles[0].circulation;
    assert!(received.is_finite() && received > 0.0);
    assert!((received / expected - 1.0).abs() < 2e-12);

    for separation in [1.0, outside, 2.0] {
        let disconnected = [
            initial[0],
            VortexParticle::new([separation, 0.0], 1.0),
        ];
        let unchanged = diffuse(&disconnected, &[1.0; 2], 1.0, 3, law);
        assert_particles_bitwise_equal(&unchanged.particles, &disconnected);
    }
}

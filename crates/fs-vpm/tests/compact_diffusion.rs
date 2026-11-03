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

fn preparation_budget(particles: usize) -> fs_vpm::diffusion::sparse::PreparationBudget {
    fs_vpm::diffusion::sparse::PreparationBudget {
        max_particles: particles,
        max_edges: usize::MAX,
        max_candidate_visits: u64::MAX,
        max_work_units: u64::MAX,
        max_live_bytes: usize::MAX,
    }
}

#[test]
fn g0_prepared_compact_graph_matches_independent_membership_and_dense_bits() {
    use fs_vpm::diffusion::sparse::PreparedPse;

    // Input order is deliberately unrelated to spatial bins. Include negative
    // coordinates, coincident sites, an exact support boundary, and an isolated
    // negative-zero circulation whose bits must survive omitted zero exchanges.
    let positions = [
        [1.1, 0.02],
        [-0.1, 0.0],
        [0.0, 0.0],
        [0.75, -0.3],
        [-0.75, 0.2],
        [2.5, 0.0],
        [0.1, 0.85],
        [0.2, 0.3],
        [1.0, 0.0],
        [0.0, 0.0],
    ];
    let strengths = [1.0, -0.3, 0.7, 0.2, -0.2, -0.0, 0.9, -1.1, 0.15, 0.65];
    let areas = [0.25, 0.5, 0.125, 0.2, 0.4, 1.0, 0.35, 0.3, 0.45, 0.15];
    let particles: Vec<_> = positions
        .iter()
        .zip(strengths)
        .map(|(position, strength)| VortexParticle::new(*position, strength))
        .collect();
    let law = PseDiffusion::compact_wendland(0.03, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prepared = PreparedPse::prepare_with_cx(
            cx, &particles, &areas, law, preparation_budget(particles.len()),
        )
        .unwrap();
        let mut expected_edges = Vec::new();
        for i in 0..particles.len() {
            for j in i + 1..particles.len() {
                let dx = particles[i].pos[0] - particles[j].pos[0];
                let dy = particles[i].pos[1] - particles[j].pos[1];
                if dx.hypot(dy) < 1.0 {
                    expected_edges.push([i, j]);
                }
            }
        }
        assert_eq!(prepared.edges(), expected_edges.as_slice());
        assert_eq!(prepared.areas_m2(), &areas);
        assert_eq!(prepared.law().kernel(), law.kernel());
        assert_eq!(prepared.law().support_radius_m(), law.support_radius_m());
        assert_eq!(
            prepared.law().kinematic_viscosity_m2_s(), law.kinematic_viscosity_m2_s(),
        );
        for (expected, actual) in positions.iter().zip(prepared.positions()) {
            assert_eq!(expected.map(f64::to_bits), actual.map(f64::to_bits));
        }
        assert_eq!(
            prepared.preparation().budget.cost_charged,
            prepared.preparation().work_units,
        );

        for (case, dt, steps) in [(0, 0.025, 7), (1, 0.4, 3)] {
            let mut state = particles.clone();
            if case == 1 {
                for (i, particle) in state.iter_mut().enumerate() {
                    particle.circulation = if i == 5 { -0.0 } else { 0.15 * i as f64 - 1.0 };
                }
            }
            let dense = diffuse_with_cx(
                cx, &state, &areas, dt, steps, law, roomy_budget(state.len(), steps),
            )
            .unwrap();
            let sparse = prepared
                .diffuse_with_cx(cx, &state, dt, steps, roomy_budget(state.len(), steps))
                .unwrap();
            assert_particles_bitwise_equal(&sparse.particles, &dense.particles);
            assert_eq!(sparse.particles[5].circulation.to_bits(), (-0.0_f64).to_bits());
            let edges = prepared.edges().len() as u64;
            assert_eq!(sparse.pair_evaluations, 2 * edges * steps as u64);
            assert!(sparse.pair_evaluations < dense.pair_evaluations);
            assert_eq!(
                sparse.work_units,
                state.len() as u64 + edges + steps as u64 * (2 * edges + 1),
            );
            assert_eq!(sparse.budget.cost_charged, sparse.work_units);
        }
    });
}

#[test]
fn g4_prepared_geometry_caps_and_cancellation_refuse_without_poisoning_reuse() {
    use fs_exec::BudgetRefusal;
    use fs_vpm::VpmError;
    use fs_vpm::diffusion::sparse::{PreparationBudget, PreparedPse};

    let particles = [
        VortexParticle::new([0.0, 0.0], 1.0),
        VortexParticle::new([0.4, 0.0], -0.3),
        VortexParticle::new([0.8, 0.0], 0.7),
        VortexParticle::new([2.0, 0.0], -0.2),
    ];
    let original = particles;
    let areas = [0.2, 0.4, 0.3, 0.1];
    let law = PseDiffusion::compact_wendland(0.03, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        assert!(matches!(
            PreparedPse::prepare_with_cx(
                cx, &particles, &areas,
                PseDiffusion::new(0.03, 1.0).unwrap(), preparation_budget(4),
            ),
            Err(VpmError::InvalidScalar {
                field: "prepared diffusion requires compact support", ..
            }),
        ));
        let prepared = PreparedPse::prepare_with_cx(
            cx, &particles, &areas, law, preparation_budget(4),
        )
        .unwrap();
        for (budget, resource) in [
            (
                PreparationBudget {
                    max_edges: prepared.edges().len() - 1,
                    ..preparation_budget(4)
                },
                "prepared edges",
            ),
            (
                PreparationBudget {
                    max_candidate_visits: prepared.preparation().candidate_visits - 1,
                    ..preparation_budget(4)
                },
                "preparation candidate visits",
            ),
            (
                PreparationBudget {
                    max_work_units: prepared.preparation().work_units - 1,
                    ..preparation_budget(4)
                },
                "preparation work",
            ),
        ] {
            assert!(matches!(
                PreparedPse::prepare_with_cx(cx, &particles, &areas, law, budget),
                Err(VpmError::CapExceeded { resource: observed, .. }) if observed == resource,
            ));
        }

        let mut moved = particles;
        moved[0].pos[0] = -0.0;
        assert!(matches!(
            prepared.diffuse_with_cx(cx, &moved, 0.1, 1, roomy_budget(4, 1)),
            Err(VpmError::InvalidScalar {
                field: "prepared position.x", particle: Some(0), ..
            }),
        ));
        assert!(matches!(
            prepared.diffuse_with_cx(cx, &particles[..3], 0.1, 1, roomy_budget(4, 1)),
            Err(VpmError::LengthMismatch {
                field: "prepared particles", expected: 4, actual: 3,
            }),
        ));
        let mut invalid = particles;
        invalid[1].circulation = f64::NAN;
        assert!(matches!(
            prepared.diffuse_with_cx(cx, &invalid, 0.1, 1, roomy_budget(4, 1)),
            Err(VpmError::InvalidScalar { field: "circulation", particle: Some(1), .. }),
        ));

        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, Budget::INFINITE, |cancelled_cx| {
            assert!(matches!(
                PreparedPse::prepare_with_cx(
                    cancelled_cx, &particles, &areas, law, preparation_budget(4),
                ),
                Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. })),
            ));
            assert!(matches!(
                prepared.diffuse_with_cx(
                    cancelled_cx, &particles, 0.1, 1, roomy_budget(4, 1),
                ),
                Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. })),
            ));
        });

        let retry = prepared
            .diffuse_with_cx(cx, &particles, 0.1, 2, roomy_budget(4, 2))
            .unwrap();
        let dense = diffuse_with_cx(
            cx, &particles, &areas, 0.1, 2, law, roomy_budget(4, 2),
        )
        .unwrap();
        assert_particles_bitwise_equal(&retry.particles, &dense.particles);
        assert_particles_bitwise_equal(&particles, &original);
        assert_eq!(prepared.positions()[0][0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(prepared.areas_m2(), &areas);
    });
}

#[test]
fn g3_sparse_exchange_work_and_storage_follow_neighbors_instead_of_all_pairs() {
    use fs_vpm::diffusion::sparse::PreparedPse;
    use std::mem::size_of;

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let mut previous_work = None;
        for n in [64, 256] {
            // An embedded chain isolates graph complexity: only adjacent
            // particles are inside the support radius as the domain grows.
            let particles: Vec<_> = (0..n)
                .map(|i| {
                    VortexParticle::new([i as f64, 0.0], (i % 7) as f64 - 2.75)
                })
                .collect();
            let areas = vec![1.0; n];
            let law = PseDiffusion::compact_wendland(0.03, 1.25).unwrap();
            let prepared = PreparedPse::prepare_with_cx(
                cx, &particles, &areas, law, preparation_budget(n),
            )
            .unwrap();
            assert_eq!(prepared.edges().len(), n - 1);
            assert!(prepared.preparation().candidate_visits <= 8 * n as u64);
            assert!(prepared.preparation().peak_live_bytes <= 80 * n);
            assert_eq!(
                prepared.preparation().budget.cost_charged,
                prepared.preparation().work_units,
            );
            let steps = 5;
            let actual = prepared
                .diffuse_with_cx(cx, &particles, 0.1, steps, roomy_budget(n, steps))
                .unwrap();
            assert_eq!(actual.pair_evaluations, 2 * (n - 1) as u64 * steps as u64);
            assert!(actual.pair_evaluations < (n * (n - 1) * steps / 16) as u64);
            assert_eq!(
                actual.peak_live_bytes,
                n * size_of::<VortexParticle>() + (n - 1) * size_of::<f64>(),
            );
            assert_eq!(actual.budget.cost_charged, actual.work_units);
            if let Some(previous) = previous_work {
                assert!(actual.work_units < 5 * previous);
            }
            previous_work = Some(actual.work_units);
            assert!(
                (total_circulation(&actual.particles) - total_circulation(&particles)).abs()
                    < 2e-12,
            );
        }
    });
}

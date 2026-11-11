use super::*;
use crate::d3q19::freesurface3::{PlanarWall3, PlanarWallFieldControls3, PlanarWallVelocityField3};
use crate::d3q19::{CollisionError3, OPP3};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn controls(n: usize) -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: n * n * n,
        maximum_wall_faces: n * n,
        maximum_relative_mass_defect: 1e-10,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn stats_bits(stats: RheologyStats) -> (u32, u32, u64, u64, u64) {
    (
        stats.floored,
        stats.capped,
        stats.max_shear.to_bits(),
        stats.tau_range.0.to_bits(),
        stats.tau_range.1.to_bits(),
    )
}

fn map_addresses(run: &FreeSurfaceRheology3) -> Vec<usize> {
    let mut addresses: Vec<_> = run
        .tau
        .values()
        .chain(run.scratch_tau.values())
        .map(|value| value as *const f64 as usize)
        .collect();
    addresses.sort_unstable();
    addresses
}

fn manufactured(
    rho: f64,
    velocity: [f64; 3],
    force: [f64; 3],
    tau: f64,
    strain: [[f64; 3]; 3],
) -> [f64; Q3] {
    let mut f = equilibrium3(rho, velocity);
    for (q, e) in E3.iter().enumerate() {
        let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
        let eu: f64 = (0..3).map(|a| e[a] * velocity[a]).sum();
        let mut contraction = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                let hermite = e[a] * e[b] - if a == b { crate::CS2 } else { 0.0 };
                contraction += hermite * (-2.0 * rho * crate::CS2 * tau * strain[a][b]);
            }
        }
        let source: f64 = (0..3)
            .map(|a| (3.0 * (e[a] - velocity[a]) + 9.0 * eu * e[a]) * force[a])
            .sum();
        f[q] += W3[q] * (contraction / (2.0 * crate::CS2 * crate::CS2) - 0.5 * source);
        // This odd third-order mode has zero density, momentum and second
        // moment. It makes a BGK substitution or a frozen TRT odd rate visible.
        f[q] += 0.012 * W3[q] * e[0] * (e[1] * e[1] - crate::CS2);
    }
    f
}

#[test]
fn trt_material_roots_drive_actual_even_and_odd_forced_collision() {
    let gate = CancelGate::new_clock_free();
    let force = [0.0002, -0.0001, 0.0003];
    let velocity = [0.03, -0.02, 0.01];
    let magic = 3.0 / 16.0;
    let original_tau = 1.7;
    let mut lattice = FreeSurface3::new(
        4,
        4,
        4,
        original_tau,
        force,
        0.0,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    lattice.set_trt_magic(magic).unwrap();
    let families = laws();
    let mut field = vec![families[0]; 64];
    let mut expected_tau = [0.0; 64];
    for lane in 0..64 {
        let x = lane % 4;
        let y = (lane / 4) % 4;
        let z = lane / 16;
        let strain = [
            [0.002, -0.001, 0.002],
            [-0.001, -0.003, 0.001],
            [0.002, 0.001, 0.002],
        ]
        .map(|row| row.map(|s| s * (1.0 + 0.02 * lane as f64)));
        // The material shear rate is sqrt(2 dev(D):dev(D)) and leaves out the trace, which is dilation (445f820cba,
        // 2026-11-05, separated material shear rheology from volumetric strain); this strain has trace 0.001 x (1 + 0.02 lane).
        let mean = (strain[0][0] + strain[1][1] + strain[2][2]) / 3.0;
        let mut rate_squared = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                let deviator = strain[a][b] - if a == b { mean } else { 0.0 };
                rate_squared += 2.0 * deviator * deviator;
            }
        }
        let rate = rate_squared.sqrt();
        field[lane] = families[(x + y + z) % families.len()];
        expected_tau[lane] = 0.5 + field[lane].viscosity(rate) / crate::CS2;
        let rho = 0.9 + 0.01 * x as f64 + 0.02 * z as f64;
        lattice.grid.set_populations(
            0,
            lane,
            manufactured(rho, velocity, force, expected_tau[lane], strain),
        );
    }
    let mut run = FreeSurfaceRheology3::with_field(
        lattice,
        RheologyField::new([4; 3], &field).unwrap(),
        &gate,
    )
    .unwrap();
    for (lane, expected) in expected_tau.iter().enumerate() {
        near(
            run.relaxation_at([lane % 4, (lane / 4) % 4, lane / 16])
                .unwrap(),
            *expected,
            4e-12,
        );
    }
    let mut expected = [0.0; Q3];
    let mut bgk_defect: f64 = 0.0;
    let mut frozen_odd_defect: f64 = 0.0;
    for (q, e) in E3.iter().enumerate() {
        let p = [2 - e.0, 2 - e.1, 2 - e.2];
        let (slot, lane) = run
            .lattice
            .grid
            .resolve_source(i64::from(p[0]), i64::from(p[1]), i64::from(p[2]))
            .unwrap();
        let rho = 0.9 + 0.01 * f64::from(p[0]) + 0.02 * f64::from(p[2]);
        let f = run.lattice.grid.populations(slot, lane);
        let eq = equilibrium3(rho, velocity);
        let opposite = OPP3[q];
        let even = 0.5 * (f[q] - eq[q] + f[opposite] - eq[opposite]);
        let odd = 0.5 * (f[q] - eq[q] - f[opposite] + eq[opposite]);
        let source = |direction: usize| {
            let e = E3[direction];
            let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
            let eu: f64 = (0..3).map(|a| e[a] * velocity[a]).sum();
            W3[direction]
                * (0..3)
                    .map(|a| (3.0 * (e[a] - velocity[a]) + 9.0 * eu * e[a]) * force[a])
                    .sum::<f64>()
        };
        let even_source = 0.5 * (source(q) + source(opposite));
        let odd_source = 0.5 * (source(q) - source(opposite));
        let tau = expected_tau[lane];
        let collide = |odd_tau: f64| {
            f[q] - even / tau - odd / odd_tau
                + (1.0 - 0.5 / tau) * even_source
                + (1.0 - 0.5 / odd_tau) * odd_source
        };
        expected[q] = collide(0.5 + magic / (tau - 0.5));
        bgk_defect = bgk_defect.max((expected[q] - collide(tau)).abs());
        frozen_odd_defect = frozen_odd_defect
            .max((expected[q] - collide(0.5 + magic / (original_tau - 0.5))).abs());
    }
    assert!(bgk_defect > 1e-6);
    assert!(frozen_odd_defect > 1e-6);
    let mass = run.lattice.ledger_mass();
    run.step_with_gate(&gate).unwrap();
    let (slot, lane) = run.lattice.grid.resolve_source(2, 2, 2).unwrap();
    for (actual, expected) in run
        .lattice
        .grid
        .populations(slot, lane)
        .into_iter()
        .zip(expected)
    {
        near(actual, expected, 3e-13);
    }
    near(run.lattice.ledger_mass(), mass, 3e-11);
    assert_eq!(
        run.lattice.collision_model(),
        CollisionModel3::Trt {
            tau: original_tau,
            magic
        }
    );
}

#[test]
fn newtonian_wall_fields_retain_native_bits_for_bgk_and_trt_on_every_face() {
    let gate = CancelGate::new_clock_free();
    let law = Rheology::Newtonian {
        nu: crate::CS2 * 0.5,
    };
    let velocity: Vec<_> = (0..16)
        .map(|i| {
            [
                0.002 + 0.0001 * (i % 4) as f64,
                -0.001 + 0.0001 * (i / 4) as f64,
            ]
        })
        .collect();
    for trt in [false, true] {
        for wall in [
            PlanarWall3::XMin,
            PlanarWall3::XMax,
            PlanarWall3::YMin,
            PlanarWall3::YMax,
            PlanarWall3::ZMin,
            PlanarWall3::ZMax,
        ] {
            let mut initial = FreeSurface3::new(
                4,
                4,
                4,
                1.0,
                [1e-6, 0.0, -2e-6],
                0.0,
                ContactModel3::Neutral,
                |_, _, _| true,
            )
            .unwrap();
            if trt {
                initial.set_trt_magic(0.19).unwrap();
            }
            let collision = initial.collision_model();
            let mut native = initial.clone();
            let mut run = FreeSurfaceRheology3::uniform(initial, law, &gate).unwrap();
            let motion = PlanarWallVelocityField3 {
                wall,
                moment_origin: [0.5, -0.75, 1.25],
                tangential_velocity: &velocity,
            };
            for _ in 0..3 {
                let expected = native
                    .step_with_planar_wall_velocity_field(motion, controls(4), &gate)
                    .unwrap();
                let (applied, report) = run
                    .step_with_planar_wall_velocity_field(motion, controls(4), &gate)
                    .unwrap();
                assert_eq!(applied.tau_range, (1.0, 1.0));
                assert_eq!(report, expected);
                assert_same_state(run.lattice(), &native);
                assert_eq!(run.lattice().collision_model(), collision);
            }
            // Repeated source refresh retains both sets of allocated wet keys.
            let addresses = map_addresses(&run);
            assert_eq!(addresses.len(), 128);
            for _ in 0..4 {
                run.set_uniform_law(law, &gate).unwrap();
            }
            assert_eq!(map_addresses(&run), addresses);
        }
    }
}

#[test]
fn actual_moving_wall_stress_refreshes_carreau_viscosity_and_exchange_together() {
    let gate = CancelGate::new_clock_free();
    let law = laws()[2];
    let velocities: Vec<_> = (0..64)
        .map(|i| {
            [
                0.02 + 0.0005 * (i % 8) as f64,
                0.005 - 0.0003 * (i / 8) as f64,
            ]
        })
        .collect();
    let motion = PlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin: [0.5, -0.75, 1.25],
        tangential_velocity: &velocities,
    };
    for trt in [false, true] {
        let mut lattice = FreeSurface3::new(
            8,
            8,
            8,
            1.7,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            |_, _, z| z < 6,
        )
        .unwrap()
        .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal);
        if trt {
            lattice.set_trt_magic(3.0 / 16.0).unwrap();
        }
        let mut run = FreeSurfaceRheology3::uniform(lattice, law, &gate).unwrap();
        let resting_tau = run.relaxation_at([3, 3, 0]).unwrap();
        let mass = run.lattice.ledger_mass();
        for _ in 0..4 {
            let applied = stats_bits(run.stats());
            let (native, expected_exchange) = run
                .lattice
                .planar_wall_field_candidate_with_relaxation(
                    motion,
                    controls(8),
                    Some(&run.tau),
                    &gate,
                )
                .unwrap();
            let (report, exchange) = run
                .step_with_planar_wall_velocity_field(motion, controls(8), &gate)
                .unwrap();
            assert_eq!(stats_bits(report), applied);
            assert_eq!(exchange, expected_exchange);
            assert_same_state(&run.lattice, &native);
            near(exchange.fluid_work, -exchange.wall_work, 2e-14);
            for face in &exchange.faces {
                let work: f64 = (0..3)
                    .map(|a| face.velocity[a] * face.fluid_impulse[a])
                    .sum();
                near(face.fluid_work, work, 2e-14);
            }
            // Independently recover every stress component at a moving-wall
            // liquid cell. Actual reflected populations supply its strain.
            let sample = run.lattice.sample_at(3, 3, 0).unwrap();
            let (slot, lane) = run.lattice.grid.resolve_source(3, 3, 0).unwrap();
            let f = run.lattice.grid.populations(slot, lane);
            let eq = equilibrium3(sample.density, sample.velocity);
            let tau = run.relaxation_at([3, 3, 0]).unwrap();
            let mut strain = [[0.0_f64; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    let stress: f64 = (0..Q3)
                        .map(|q| {
                            let e = [E3[q].0, E3[q].1, E3[q].2];
                            (f[q] - eq[q]) * f64::from(e[a] * e[b])
                        })
                        .sum();
                    strain[a][b] = -stress / (2.0 * sample.density * crate::CS2 * tau);
                }
            }
            // The material shear rate is sqrt(2 dev(D):dev(D)) without the trace (445f820cba, 2026-11-05).
            let mean = (strain[0][0] + strain[1][1] + strain[2][2]) / 3.0;
            let mut rate_squared = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let deviator = strain[a][b] - if a == b { mean } else { 0.0 };
                    rate_squared += 2.0 * deviator * deviator;
                }
            }
            let rate = rate_squared.sqrt();
            assert!(rate > 0.0);
            near(
                tau,
                (0.5 + law.viscosity(rate) / crate::CS2)
                    .clamp(crate::rheology::TAU_FLOOR, crate::rheology::TAU_CAP),
                4e-12,
            );
        }
        assert!(run.relaxation_at([3, 3, 0]).unwrap() < resting_tau - 1e-5);
        near(run.lattice.ledger_mass(), mass, 2e-9);
    }
}

#[test]
fn trt_derived_odd_refusal_preserves_accepted_material_and_flow() {
    let gate = CancelGate::new_clock_free();
    let mut lattice = FreeSurface3::new(
        4,
        4,
        4,
        0.65,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    lattice.set_trt_magic(1e-16).unwrap();
    let mut run =
        FreeSurfaceRheology3::uniform(lattice, Rheology::Newtonian { nu: 0.05 }, &gate).unwrap();
    let before = run.lattice.clone();
    let tau = run.tau.clone();
    let stats = stats_bits(run.stats());
    assert!(
        matches!(run.set_uniform_law(Rheology::Newtonian { nu: 1.0 }, &gate),
        Err(FreeSurfaceRheologyError3::Flow(FreeSurfaceError3::Collision {
            source: CollisionError3::InvalidRelaxationTime { tau }, ..
        })) if tau == 0.5)
    );
    assert_same_state(&run.lattice, &before);
    assert_eq!(run.tau, tau);
    assert_eq!(stats_bits(run.stats()), stats);
    run.set_uniform_law(Rheology::Newtonian { nu: 0.06 }, &gate)
        .unwrap();
    let expected =
        FreeSurfaceRheology3::uniform(before, Rheology::Newtonian { nu: 0.06 }, &gate).unwrap();
    assert_eq!(run.tau, expected.tau);
    assert_eq!(
        run.lattice.collision_model(),
        expected.lattice.collision_model()
    );
}

#[test]
fn moving_wall_endpoint_material_refusal_rolls_back_and_retries_exactly() {
    let gate = CancelGate::new_clock_free();
    let mut lattice = FreeSurface3::new(
        4,
        4,
        4,
        crate::rheology::TAU_FLOOR,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    lattice.set_trt_magic(1e-18).unwrap();
    let law = Rheology::PowerLaw { k: 0.5, n: 1.5 };
    let mut control = FreeSurfaceRheology3::uniform(lattice.clone(), law, &gate).unwrap();
    let mut run = FreeSurfaceRheology3::uniform(lattice.clone(), law, &gate).unwrap();
    let tau = run.tau.clone();
    let stats = stats_bits(run.stats());
    let velocities = [[0.02, 0.0]; 16];
    let motion = PlanarWallVelocityField3 {
        wall: PlanarWall3::ZMax,
        moment_origin: [0.0; 3],
        tangential_velocity: &velocities,
    };
    // The complete native flow/wall transaction succeeds. Only the new endpoint
    // material rate makes the derived TRT odd relaxation unrepresentable.
    let (_, exchange) = run
        .lattice
        .planar_wall_field_candidate_with_relaxation(motion, controls(4), Some(&tau), &gate)
        .unwrap();
    assert!(exchange.fluid_work > 0.0);
    assert!(
        matches!(run.step_with_planar_wall_velocity_field(motion, controls(4), &gate),
        Err(FreeSurfaceRheologyError3::Flow(FreeSurfaceError3::Collision {
            source: CollisionError3::InvalidRelaxationTime { tau }, ..
        })) if tau == 0.5)
    );
    assert_same_state(&run.lattice, &lattice);
    assert_eq!(run.tau, tau);
    assert_eq!(stats_bits(run.stats()), stats);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run.step_with_planar_wall_velocity_field(motion, controls(4), &cancelled),
        Err(FreeSurfaceRheologyError3::Flow(
            FreeSurfaceError3::Cancelled
        ))
    ));
    assert_same_state(&run.lattice, &lattice);
    let still = [[0.0; 2]; 16];
    let stationary = PlanarWallVelocityField3 {
        tangential_velocity: &still,
        ..motion
    };
    let (a, exchange_a) = run
        .step_with_planar_wall_velocity_field(stationary, controls(4), &gate)
        .unwrap();
    let (b, exchange_b) = control
        .step_with_planar_wall_velocity_field(stationary, controls(4), &gate)
        .unwrap();
    assert_eq!(stats_bits(a), stats_bits(b));
    assert_eq!(exchange_a, exchange_b);
    assert_same_state(&run.lattice, &control.lattice);
    assert_eq!(run.tau, control.tau);
}

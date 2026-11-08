//! Independent mechanical kicks around the native finite-piston wave field.

use super::*;
use fs_gas::body_force::{BodyForce, ForceControl};
use fs_gas::dynamics::nonideal::{RealGasWallContact, RealGasWallOperator};
use fs_gas::multicomponent::wall_exchange::{WallContact, WallControl, WallDrag, WallThermal};

fn force_controls() -> ForceControl {
    ForceControl { maximum_mach: 4.0 }
}

fn chamber(gas: RealGas, states: Vec<Conserved>) -> RealGasEuler1d {
    RealGasEuler1d::new(
        gas,
        Length::new(1.0 / 3.0),
        0.07,
        states,
        RealGasBoundaries::Ends {
            left: RealGasBoundary::Reflecting,
            right: RealGasBoundary::Reflecting,
        },
        || false,
    )
    .unwrap()
}

// An independent kinetic-energy difference, not the production midpoint-work
// kernel. The source cannot change density or internal energy at fixed volume.
fn kick(q: Conserved, impulse_density: f64) -> (Conserved, f64) {
    let momentum = q.momentum + impulse_density;
    let work = (momentum.powi(2) - q.momentum.powi(2)) / (2.0 * q.density);
    (
        Conserved {
            momentum,
            energy: q.energy + work,
            ..q
        },
        work,
    )
}

#[test]
fn distributed_force_drives_an_initially_pressure_balanced_piston_in_the_same_step() {
    let mut flow = field(gas(), false);
    let before = flow.totals(|| false).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 0.0,
    };
    let acceleration = 200_000.0;
    let mut piston = Piston::new(5.0, 0.0).unwrap();
    let report = flow
        .step_with_piston_and_body_force(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            load,
            mechanical_controls(),
            BodyForce::UniformAcceleration(acceleration),
            force_controls(),
            || false,
        )
        .unwrap();
    assert!(piston.velocity_m_s() > 0.0);
    assert!(report.displacement_m > 0.0);
    assert!(report.gas.source_change.work > 0.0);
    assert_eq!(report.gas.wall_heat, 0.0);
    near(
        report.gas.source_change.impulse,
        before.mass * acceleration * report.gas.elapsed.value(),
        1e-12,
    );
    near(
        report.gas.after.energy - before.energy
            + report.kinetic_energy_change_j
            + report.exterior_pressure_work_j
            - report.gas.source_change.work,
        0.0,
        3e-7,
    );
    near(report.gas.conservation_residual.energy, 0.0, 3e-7);
    near(report.combined_momentum_residual_n_s, 0.0, 2e-10);
}

#[test]
fn compression_matches_independent_half_kicks_with_moved_volume_and_signed_braking_work() {
    let dt = 2e-5;
    for force_density in [false, true] {
        let gas = gas();
        let q = state(&gas, 10.0, -6.0, 450.0);
        let mut flow = chamber(gas.clone(), vec![q; 3]);
        let initial_volume: f64 = flow.mesh().volumes().iter().sum();
        let acceleration = 200_000.0;
        let density_force = q.density * acceleration;
        let (first, first_work) = kick(q, 0.5 * dt * density_force);
        let mut manual = chamber(gas, vec![first; 3]);
        let load = PistonLoad {
            exterior_pressure_pa: flow.primitives()[0].pressure,
            external_force_n: -200.0,
        };
        let mut expected_piston = Piston::new(50.0, -10.0).unwrap();
        let pressure = manual
            .step_with_piston(
                &mut expected_piston,
                Time::new(dt),
                muscl(),
                load,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(pressure.gas.elapsed.value(), dt);
        let expected: Vec<_> = manual
            .states()
            .iter()
            .map(|q| {
                let density_force = if force_density {
                    density_force
                } else {
                    q.density * acceleration
                };
                kick(*q, 0.5 * dt * density_force)
            })
            .collect();
        let expected_work = first_work * initial_volume
            + expected
                .iter()
                .zip(manual.mesh().volumes())
                .map(|((_, work), volume)| work * volume)
                .sum::<f64>();
        let mut piston = Piston::new(50.0, -10.0).unwrap();
        let force = if force_density {
            BodyForce::UniformForceDensity(density_force)
        } else {
            BodyForce::UniformAcceleration(acceleration)
        };
        let actual = flow
            .step_with_piston_and_body_force(
                &mut piston,
                Time::new(dt),
                muscl(),
                load,
                mechanical_controls(),
                force,
                force_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(actual.gas.elapsed.value(), dt);
        assert!(flow.mesh().positions()[3] < 1.0);
        near(piston.velocity_m_s(), expected_piston.velocity_m_s(), 2e-8);
        near(
            actual.gas.boundary_change.energy,
            pressure.gas.boundary_change.energy,
            2e-6,
        );
        for (i, (expected, _)) in expected.iter().enumerate() {
            near(flow.states()[i].density, expected.density, 1e-10);
            near(flow.states()[i].momentum, expected.momentum, 1e-8);
            near(flow.states()[i].energy, expected.energy, 2e-5);
            near(
                flow.primitives()[i].temperature,
                manual.primitives()[i].temperature,
                2e-8,
            );
            near(flow.mesh().volumes()[i], manual.mesh().volumes()[i], 1e-12);
        }
        let final_volume: f64 = flow.mesh().volumes().iter().sum();
        let impulse = if force_density {
            0.5 * dt * density_force * (initial_volume + final_volume)
        } else {
            dt * acceleration * actual.gas.before.mass
        };
        near(actual.gas.source_change.impulse, impulse, 1e-10);
        near(actual.gas.source_change.work, expected_work, 2e-7);
        assert!(actual.gas.source_change.work < 0.0);
        near(actual.total_energy_residual_j, 0.0, 3e-7);
        near(actual.gas.conservation_residual.energy, 0.0, 3e-7);
    }
}

fn wall_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-7,
        relative_energy_tolerance: 1e-10,
        maximum_iterations: 128,
    }
}

fn sourced_controls() -> PistonControl {
    PistonControl {
        momentum_tolerance_n_s: 2e-9,
        energy_tolerance_j: 5e-6,
        ..mechanical_controls()
    }
}

fn thermal_case(drag: WallDrag) -> (RealGasEuler1d, RealGasWallOperator<'static>) {
    let gas = gas();
    let q = state(&gas, 10.0, 0.0, 450.0);
    let flow = chamber(gas, vec![q; 3]);
    let wall = RealGasWallContact {
        contact: WallContact {
            drag,
            thermal: WallThermal::Lumped {
                temperature_k: 500.0,
                heat_capacity_j_per_k: 80.0,
                conductance_w_per_k: 1000.0,
            },
        },
        heat_storage: None,
    };
    let walls = RealGasWallOperator::new(&flow, vec![wall; 3], || false).unwrap();
    (flow, walls)
}

#[test]
fn body_work_wall_storage_sliding_drive_and_suspension_close_one_physical_balance() {
    let wall_speed = 25.0;
    let (mut flow, mut walls) = thermal_case(WallDrag::MovingLinear {
        rate_per_s: 3000.0,
        wall_velocity_m_s: wall_speed,
    });
    let mut piston = Piston::new(3.0, -1.0).unwrap();
    let spring = PistonSuspension::new(4000.0, 30.0, 0.95).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 200.0,
    };
    let before = flow.totals(|| false).unwrap();
    let report = flow
        .step_with_piston_sources(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            load,
            Some(spring),
            sourced_controls(),
            BodyForce::UniformAcceleration(200_000.0),
            force_controls(),
            Some((&mut walls, wall_controls())),
            || false,
        )
        .unwrap();
    let thermal = walls.last_report().unwrap();
    let after = flow.totals(|| false).unwrap();
    let wall_storage: f64 = walls
        .contacts()
        .iter()
        .map(|contact| {
            let WallThermal::Lumped {
                temperature_k,
                heat_capacity_j_per_k,
                ..
            } = contact.thermal
            else {
                panic!("finite wall must retain its thermal owner")
            };
            heat_capacity_j_per_k * (temperature_k - 500.0)
        })
        .sum();
    assert!(thermal.wall_storage_change_j < 0.0);
    assert!(thermal.drag_heating_j > 0.0);
    assert!(report.damping_energy_j > 0.0);
    assert!(report.piston.gas.source_change.work > 0.0);
    assert_eq!(thermal.external_heat_to_gas_j, 0.0);
    near(thermal.wall_storage_change_j, wall_storage, 1e-10);
    near(
        thermal.mechanical_work_to_gas_j,
        -wall_speed * thermal.wall_impulse_ns,
        2e-9,
    );
    near(
        report.piston.gas.source_change.impulse,
        before.mass * 200_000.0 * report.piston.gas.elapsed.value(),
        1e-11,
    );
    near(
        after.energy - before.energy
            + 1.5 * (piston.velocity_m_s().powi(2) - 1.0)
            + spring.stored_energy_j(flow.mesh().positions()[3]).unwrap()
            - spring.stored_energy_j(1.0).unwrap()
            + report.damping_energy_j
            + report.piston.exterior_pressure_work_j
            - report.applied_force_work_j
            + wall_storage
            - thermal.mechanical_work_to_gas_j
            - report.piston.gas.source_change.work,
        0.0,
        5e-6,
    );
    near(
        after.momentum - before.momentum
            + 3.0 * (piston.velocity_m_s() + 1.0)
            + report.piston.fixed_wall_impulse_n_s
            + thermal.wall_impulse_ns
            - report.suspension_impulse_n_s
            - report.piston.gas.elapsed.value()
                * (load.external_force_n - 0.07 * load.exterior_pressure_pa)
            - report.piston.gas.source_change.impulse,
        0.0,
        2e-9,
    );
    near(thermal.gas_wall_storage_residual_j, 0.0, 5e-6);
    near(report.total_energy_residual_j, 0.0, 5e-6);
}

#[test]
fn a_late_mach_refusal_rolls_back_all_owners_then_retries_the_entire_interval() {
    let (mut flow, mut walls) = thermal_case(WallDrag::None);
    let (mut reference, mut reference_walls) = thermal_case(WallDrag::None);
    let states = flow.states().to_vec();
    let primitive = flow.primitives().to_vec();
    let positions = flow.mesh().positions().to_vec();
    let volumes = flow.mesh().volumes().to_vec();
    let contacts = walls.contacts().to_vec();
    let mut piston = Piston::new(5.0, 0.0).unwrap();
    let mut reference_piston = piston;
    let load = PistonLoad {
        exterior_pressure_pa: primitive[0].pressure,
        external_force_n: 0.0,
    };
    // The first half supplies 1 m/s and passes; the second reaches about 2 m/s.
    let force_controls = ForceControl {
        maximum_mach: 1.5 / primitive[0].sound_speed,
    };
    let error = flow
        .step_with_piston_sources(
            &mut piston,
            Time::new(2e-5),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            load,
            None,
            sourced_controls(),
            BodyForce::UniformAcceleration(100_000.0),
            force_controls,
            Some((&mut walls, wall_controls())),
            || false,
        )
        .unwrap_err();
    let PistonError::StepRetries { last, .. } = error else {
        panic!("late source refusal must exhaust the complete interval")
    };
    assert!(matches!(*last, PistonError::RealGas(error)
        if *error == fs_gas::dynamics::nonideal::RealGasFlowError::Invalid(
            "body-force state exceeds its authored Mach bound")));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.primitives(), primitive);
    assert_eq!(flow.mesh().positions(), positions);
    assert_eq!(flow.mesh().volumes(), volumes);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(piston, reference_piston);
    assert_eq!(walls.contacts(), contacts);
    assert_eq!(walls.last_report(), None);
    let retried = flow
        .step_with_piston_sources(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            load,
            None,
            sourced_controls(),
            BodyForce::UniformAcceleration(100_000.0),
            force_controls,
            Some((&mut walls, wall_controls())),
            || false,
        )
        .unwrap();
    assert_eq!(retried.piston.gas.halvings, 1);
    assert!(!retried.piston.gas.first_order_retry);
    let direct = reference
        .step_with_piston_sources(
            &mut reference_piston,
            Time::new(1e-5),
            muscl(),
            load,
            None,
            sourced_controls(),
            BodyForce::UniformAcceleration(100_000.0),
            force_controls,
            Some((&mut reference_walls, wall_controls())),
            || false,
        )
        .unwrap();
    assert_eq!(retried.piston.gas.elapsed, direct.piston.gas.elapsed);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.primitives(), reference.primitives());
    assert_eq!(flow.mesh().positions(), reference.mesh().positions());
    assert_eq!(piston, reference_piston);
    assert_eq!(walls.contacts(), reference_walls.contacts());
    assert_eq!(walls.last_report(), reference_walls.last_report());
}

#[test]
fn zero_force_preserves_the_existing_thermal_suspension_path_exactly() {
    let (mut flow, mut walls) = thermal_case(WallDrag::None);
    let (mut reference, mut reference_walls) = thermal_case(WallDrag::None);
    let mut piston = Piston::new(5.0, -1.0).unwrap();
    let mut reference_piston = piston;
    let spring = PistonSuspension::new(4000.0, 30.0, 0.95).unwrap();
    let expected = reference
        .step_with_suspended_piston_and_walls(
            &mut reference_piston,
            Time::new(2e-5),
            muscl(),
            PistonLoad::default(),
            spring,
            sourced_controls(),
            &mut reference_walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    let actual = flow
        .step_with_piston_sources(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            PistonLoad::default(),
            Some(spring),
            sourced_controls(),
            BodyForce::UniformAcceleration(0.0),
            force_controls(),
            Some((&mut walls, wall_controls())),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(piston, reference_piston);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.primitives(), reference.primitives());
    assert_eq!(flow.mesh().volumes(), reference.mesh().volumes());
    assert_eq!(walls.contacts(), reference_walls.contacts());
    assert_eq!(walls.last_report(), reference_walls.last_report());
}

//! Heat-driven work on the original real-gas field and finite wall inventories.

use super::*;
use fs_gas::dynamics::piston::{Piston, PistonControl, PistonError, PistonLoad, PistonSuspension};

fn piston_controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-20.0, 20.0],
        momentum_tolerance_n_s: 2e-9,
        energy_tolerance_j: 5e-6,
        maximum_iterations: 64,
    }
}

fn muscl() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::Muscl,
        ..controls()
    }
}

fn chamber(gas: RealGas, temperature: f64) -> RealGasEuler1d {
    let q = state(&gas, 10.0, 0.0, temperature);
    RealGasEuler1d::new(
        gas,
        Length::new(1.0 / 3.0),
        0.07,
        vec![q; 3],
        RealGasBoundaries::Ends {
            left: RealGasBoundary::Reflecting,
            right: RealGasBoundary::Reflecting,
        },
        || false,
    )
    .unwrap()
}

fn reservoir(temperature: f64, conductance: f64) -> RealGasWallContact<'static> {
    contact(WallThermal::Reservoir {
        temperature_k: temperature,
        conductance_w_per_k: conductance,
    })
}

fn reservoir_temperature(rho: f64, volume: f64, initial: f64, bath: f64, integrated: f64) -> f64 {
    root(initial.min(bath), initial.max(bath), |temperature| {
        rho * volume * (internal_energy(rho, temperature) - internal_energy(rho, initial))
            - integrated * (bath - temperature)
    })
}

#[test]
fn hot_original_finite_walls_drive_a_balanced_piston_during_the_same_interval() {
    let mut flow = chamber(gas(), 450.0);
    let mut walls = RealGasWallOperator::new(
        &flow,
        (0..3).map(|_| finite_wall(800.0, 1e5)).collect(),
        || false,
    )
    .unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 0.0,
    };
    let mut piston = Piston::new(5.0, 0.0).unwrap();
    let before = flow.totals(|| false).unwrap();
    let initial_wall: f64 = (0..3).map(|i| 0.75 * specific_wall_energy(&walls, i)).sum();
    let mut exterior_work = 0.0;
    for _ in 0..4 {
        let report = flow
            .step_with_piston_and_walls(
                &mut piston,
                Time::new(2e-5),
                muscl(),
                load,
                piston_controls(),
                &mut walls,
                wall_controls(),
                || false,
            )
            .unwrap();
        let thermal = walls.last_report().unwrap();
        assert!(piston.velocity_m_s() > 0.0);
        assert!(report.displacement_m > 0.0);
        assert!(report.gas.wall_heat > 0.0);
        assert_eq!(thermal.external_heat_to_gas_j, 0.0);
        assert_eq!(report.gas.wall_heat, thermal.heat_to_gas_j);
        assert_eq!(thermal.elapsed, report.gas.elapsed);
        near(report.total_energy_residual_j, 0.0, 5e-6);
        near(report.gas.conservation_residual.energy, 0.0, 5e-6);
        near(thermal.gas_wall_storage_residual_j, 0.0, 5e-6);
        exterior_work += report.exterior_pressure_work_j;
    }
    let after = flow.totals(|| false).unwrap();
    let final_wall: f64 = (0..3).map(|i| 0.75 * specific_wall_energy(&walls, i)).sum();
    assert!(final_wall < initial_wall);
    near(after.mass, before.mass, 1e-13);
    near(
        after.energy - before.energy
            + 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2)
            + final_wall
            - initial_wall
            + exterior_work,
        0.0,
        2e-5,
    );
}

#[test]
fn compression_heat_halves_match_independent_pr_energy_at_the_actual_moved_volumes() {
    let source = gas();
    let mut flow = chamber(source.clone(), 450.0);
    let dt = 2e-5;
    let conductance = 1e5;
    let bath = 700.0;
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: -200.0,
    };
    let first_temperature = reservoir_temperature(
        10.0,
        flow.mesh().volumes()[0],
        450.0,
        bath,
        0.5 * dt * conductance,
    );
    let first_heat: f64 = flow
        .mesh()
        .volumes()
        .iter()
        .map(|volume| {
            10.0 * volume
                * (internal_energy(10.0, first_temperature) - internal_energy(10.0, 450.0))
        })
        .sum();
    let mut manual = chamber(source, first_temperature);
    let mut manual_piston = Piston::new(50.0, -10.0).unwrap();
    let pressure = manual
        .step_with_piston(
            &mut manual_piston,
            Time::new(dt),
            muscl(),
            load,
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(pressure.gas.elapsed.value(), dt);
    let expected: Vec<_> = manual
        .primitives()
        .iter()
        .zip(manual.mesh().volumes())
        .map(|(p, &volume)| {
            let temperature = reservoir_temperature(
                p.density,
                volume,
                p.temperature,
                bath,
                0.5 * dt * conductance,
            );
            let heat = p.density
                * volume
                * (internal_energy(p.density, temperature)
                    - internal_energy(p.density, p.temperature));
            (temperature, heat)
        })
        .collect();
    let mut walls =
        RealGasWallOperator::new(&flow, vec![reservoir(bath, conductance); 3], || false).unwrap();
    let mut piston = Piston::new(50.0, -10.0).unwrap();
    let actual = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(dt),
            muscl(),
            load,
            piston_controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(actual.gas.elapsed.value(), dt);
    assert!(flow.mesh().positions()[3] < 1.0);
    near(piston.velocity_m_s(), manual_piston.velocity_m_s(), 2e-8);
    near(
        actual.gas.boundary_change.energy,
        pressure.gas.boundary_change.energy,
        2e-5,
    );
    for (p, &(temperature, _)) in flow.primitives().iter().zip(&expected) {
        near(p.temperature, temperature, 2e-7);
    }
    near(
        actual.gas.wall_heat,
        first_heat + expected.iter().map(|(_, heat)| heat).sum::<f64>(),
        2e-5,
    );
    near(actual.total_energy_residual_j, 0.0, 5e-6);
}

#[test]
fn sliding_drive_heat_spring_and_damping_keep_their_distinct_work_accounts() {
    let mut flow = chamber(gas(), 450.0);
    let mut input = reservoir(650.0, 2e4);
    let wall_speed = 25.0;
    input.contact.drag = WallDrag::MovingLinear {
        rate_per_s: 3000.0,
        wall_velocity_m_s: wall_speed,
    };
    let mut walls = RealGasWallOperator::new(&flow, vec![input; 3], || false).unwrap();
    let mut piston = Piston::new(3.0, -1.0).unwrap();
    let spring = PistonSuspension::new(4000.0, 30.0, 0.95).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 200.0,
    };
    let before = flow.totals(|| false).unwrap();
    let report = flow
        .step_with_suspended_piston_and_walls(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            load,
            spring,
            piston_controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    let thermal = walls.last_report().unwrap();
    assert!(thermal.drag_heating_j > 0.0);
    assert!(report.damping_energy_j > 0.0);
    assert!(thermal.heat_to_gas_j > 0.0);
    near(
        thermal.mechanical_work_to_gas_j,
        -wall_speed * thermal.wall_impulse_ns,
        2e-9,
    );
    near(
        report.piston.gas.wall_momentum_impulse,
        -thermal.wall_impulse_ns,
        1e-14,
    );
    let after = flow.totals(|| false).unwrap();
    near(
        after.energy - before.energy
            + report.piston.kinetic_energy_change_j
            + report.spring_energy_change_j
            + report.damping_energy_j
            + report.piston.exterior_pressure_work_j
            - report.applied_force_work_j
            - thermal.external_heat_to_gas_j
            - thermal.mechanical_work_to_gas_j,
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
                * (load.external_force_n - 0.07 * load.exterior_pressure_pa),
        0.0,
        2e-9,
    );
    near(report.total_energy_residual_j, 0.0, 5e-6);
}

#[test]
fn original_datum_shifts_and_an_inactive_wall_cannot_change_piston_dynamics() {
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    let mut a = chamber(gas(), 450.0);
    let mut b = chamber(shifted, 450.0);
    let mut pa = Piston::new(4.0, -4.0).unwrap();
    let mut pb = pa;
    let mut wa = RealGasWallOperator::new(&a, vec![reservoir(700.0, 5e4); 3], || false).unwrap();
    let mut wb = RealGasWallOperator::new(&b, vec![reservoir(700.0, 5e4); 3], || false).unwrap();
    let x = a
        .step_with_piston_and_walls(
            &mut pa,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            &mut wa,
            wall_controls(),
            || false,
        )
        .unwrap();
    let y = b
        .step_with_piston_and_walls(
            &mut pb,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            &mut wb,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(x, y);
    assert_eq!(a.states(), b.states());
    assert_eq!(a.primitives(), b.primitives());
    assert_eq!(wa.last_report(), wb.last_report());

    let source = gas();
    let mut dry = chamber(source.clone(), 450.0);
    let mut idle = chamber(source, 450.0);
    let mut pd = Piston::new(4.0, -4.0).unwrap();
    let mut pi = pd;
    let mut walls =
        RealGasWallOperator::new(&idle, vec![contact(WallThermal::Adiabatic); 3], || false)
            .unwrap();
    let dry_report = dry
        .step_with_piston(
            &mut pd,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            || false,
        )
        .unwrap();
    let idle_report = idle
        .step_with_piston_and_walls(
            &mut pi,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(dry.states(), idle.states());
    assert_eq!(dry.primitives(), idle.primitives());
    assert_eq!(dry.mesh(), idle.mesh());
    assert_eq!(pd, pi);
    assert_eq!(
        dry_report.gas.boundary_change,
        idle_report.gas.boundary_change
    );
    assert_eq!(idle_report.gas.wall_heat, 0.0);
    assert_eq!(idle_report.gas.wall_work, 0.0);
}

#[test]
fn a_refused_second_heat_half_retries_from_original_wall_storage_without_losing_muscl() {
    let limited = RealGasDomain {
        energy: GasEnergySolve {
            maximum_temperature_k: 451.0,
            ..domain().energy
        },
        ..domain()
    };
    let mut flow = chamber(
        bind(common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]), limited),
        450.0,
    );
    let mut walls = RealGasWallOperator::new(
        &flow,
        vec![
            contact(WallThermal::Lumped {
                temperature_k: 700.0,
                heat_capacity_j_per_k: 500.0,
                conductance_w_per_k: 1e5,
            });
            3
        ],
        || false,
    )
    .unwrap();
    let mut piston = Piston::new(5.0, 0.0).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 0.0,
    };
    let controls = PistonControl {
        midpoint_velocity_bracket_m_s: [-0.01, 0.01],
        ..piston_controls()
    };
    let before = (
        flow.states().to_vec(),
        flow.primitives().to_vec(),
        flow.mesh().clone(),
        piston,
        walls.contacts().to_vec(),
    );
    let gas_mass = 10.0 * flow.mesh().volumes()[0];
    let first = root(450.0, 700.0, |temperature| {
        let heat = gas_mass * (internal_energy(10.0, temperature) - internal_energy(10.0, 450.0));
        heat - (700.0 - heat / 500.0 - temperature)
    });
    assert!(first < 451.0);
    let first_heat = gas_mass * (internal_energy(10.0, first) - internal_energy(10.0, 450.0));
    let first_wall = 700.0 - first_heat / 500.0;
    assert!(first_wall < 700.0);
    let second = root(first, first_wall, |temperature| {
        let heat = gas_mass * (internal_energy(10.0, temperature) - internal_energy(10.0, first));
        heat - (first_wall - heat / 500.0 - temperature)
    });
    assert!(second > 451.0);
    let error = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(2e-5),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            load,
            controls,
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, PistonError::StepRetries { .. }));
    assert_eq!(flow.states(), before.0);
    assert_eq!(flow.primitives(), before.1);
    assert_eq!(flow.mesh(), &before.2);
    assert_eq!(piston, before.3);
    assert_eq!(walls.contacts(), before.4);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(walls.last_report(), None);
    let retry = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(2e-5),
            muscl(),
            load,
            controls,
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert!(retry.gas.halvings > 0);
    assert!(!retry.gas.first_order_retry);
    assert!(flow.primitives().iter().all(|p| p.temperature <= 451.0));
    assert_eq!(walls.last_report().unwrap().elapsed, retry.gas.elapsed);
}

#[test]
fn cancellation_after_prior_motion_preserves_gas_piston_wall_binding_and_reports() {
    let make = || {
        let flow = chamber(gas(), 450.0);
        let walls = RealGasWallOperator::new(
            &flow,
            (0..3).map(|_| finite_wall(750.0, 1e5)).collect(),
            || false,
        )
        .unwrap();
        (flow, Piston::new(5.0, -2.0).unwrap(), walls)
    };
    let (mut reference, mut reference_piston, mut reference_walls) = make();
    let (mut flow, mut piston, mut walls) = make();
    for (gas, p, w) in [
        (&mut reference, &mut reference_piston, &mut reference_walls),
        (&mut flow, &mut piston, &mut walls),
    ] {
        gas.step_with_piston_and_walls(
            p,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            w,
            wall_controls(),
            || false,
        )
        .unwrap();
    }
    let mut polls = 0;
    let expected = reference
        .step_with_piston_and_walls(
            &mut reference_piston,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            &mut reference_walls,
            wall_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 100);
    let before = (
        flow.states().to_vec(),
        flow.primitives().to_vec(),
        flow.mesh().clone(),
        flow.time(),
        piston,
        walls.contacts().to_vec(),
        walls.last_report(),
    );
    let mut seen = 0;
    let error = flow.step_with_piston_and_walls(
        &mut piston,
        Time::new(2e-6),
        muscl(),
        PistonLoad::default(),
        piston_controls(),
        &mut walls,
        wall_controls(),
        || {
            seen += 1;
            seen == polls
        },
    );
    assert_eq!(error, Err(PistonError::Cancelled));
    assert_eq!(flow.states(), before.0);
    assert_eq!(flow.primitives(), before.1);
    assert_eq!(flow.mesh(), &before.2);
    assert_eq!(flow.time(), before.3);
    assert_eq!(piston, before.4);
    assert_eq!(walls.contacts(), before.5);
    assert_eq!(walls.last_report(), before.6);
    let actual = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            piston_controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.primitives(), reference.primitives());
    assert_eq!(flow.mesh(), reference.mesh());
    assert_eq!(piston, reference_piston);
    assert_eq!(walls.contacts(), reference_walls.contacts());
    assert_eq!(walls.last_report(), reference_walls.last_report());
}

#[test]
fn shared_wall_and_link_inventories_continue_through_both_halves_with_one_work_budget() {
    use fs_gas::dynamics::nonideal::SharedWallLink;

    let make = |maximum_source_evaluations| {
        let flow = chamber(gas(), 450.0);
        let first = shared_wall(finite_wall(800.0, 0.0), &[1e5, 2e5]);
        let mut second = shared_wall(finite_wall(600.0, 0.0), &[8e4]);
        second.ports[0].gas = SharedWallGas::Branch { branch: 0, cell: 2 };
        let walls =
            RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
                .unwrap()
                .with_shared_walls(
                    vec![first, second],
                    SharedWallControl {
                        maximum_source_evaluations,
                    },
                    || false,
                )
                .unwrap()
                .with_shared_wall_links(
                    vec![SharedWallLink {
                        walls: [0, 1],
                        conductance_w_per_k: 2e5,
                    }],
                    || false,
                )
                .unwrap();
        (flow, Piston::new(50.0, -6.0).unwrap(), walls)
    };
    let inventory = |walls: &RealGasWallOperator<'_>| -> f64 {
        (0..2)
            .map(|i| {
                let WallThermal::Caloric {
                    mass_kg,
                    specific_energy_j_kg,
                    ..
                } = walls.shared_wall(i).unwrap().storage.thermal
                else {
                    panic!("original shared inventory");
                };
                mass_kg * specific_energy_j_kg
            })
            .sum()
    };
    let (mut flow, mut piston, mut walls) = make(1_000_000);
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: -100.0,
    };
    let mut first_work = 0;
    for interval in 0..2 {
        let before = flow.totals(|| false).unwrap();
        let wall_before = inventory(&walls);
        let report = flow
            .step_with_piston_and_walls(
                &mut piston,
                Time::new(2e-5),
                muscl(),
                load,
                piston_controls(),
                &mut walls,
                wall_controls(),
                || false,
            )
            .unwrap();
        let thermal = walls.last_report().unwrap();
        let shared = walls.last_shared_wall_report().unwrap();
        let wall_delta = inventory(&walls) - wall_before;
        near(wall_delta, thermal.wall_storage_change_j, 5e-6);
        near(
            flow.totals(|| false).unwrap().energy - before.energy
                + report.kinetic_energy_change_j
                + wall_delta
                + report.exterior_pressure_work_j
                - report.external_force_work_j,
            0.0,
            5e-6,
        );
        assert_eq!(thermal.external_heat_to_gas_j, 0.0);
        assert_eq!(shared.elapsed, report.gas.elapsed);
        assert!(shared.contact_solves >= 6);
        let link = walls.shared_wall_link_reports().unwrap()[0];
        assert_eq!(link.elapsed.value(), 0.5 * report.gas.elapsed.value());
        assert!(link.heat_a_to_b_j > 0.0);
        if interval == 0 {
            first_work = shared.source_evaluations;
        }
    }
    assert!(first_work > 10);
    let (mut refused, mut rp, mut rw) = make(first_work - 1);
    let before = (
        refused.states().to_vec(),
        refused.primitives().to_vec(),
        refused.mesh().clone(),
        rp,
        inventory(&rw),
        rw.shared_wall(0).unwrap().clone(),
        rw.shared_wall(1).unwrap().clone(),
    );
    let error = refused
        .step_with_piston_and_walls(
            &mut rp,
            Time::new(2e-5),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            load,
            piston_controls(),
            &mut rw,
            wall_controls(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, PistonError::StepRetries { .. }));
    assert_eq!(refused.states(), before.0);
    assert_eq!(refused.primitives(), before.1);
    assert_eq!(refused.mesh(), &before.2);
    assert_eq!(rp, before.3);
    assert_eq!(inventory(&rw), before.4);
    assert_eq!(rw.shared_wall(0), Some(&before.5));
    assert_eq!(rw.shared_wall(1), Some(&before.6));
    assert_eq!(refused.time().value(), 0.0);
    assert_eq!(rw.last_report(), None);
    assert_eq!(rw.last_shared_wall_report(), None);
    assert_eq!(rw.shared_wall_link_reports(), None);
}

#[test]
fn post_heat_acoustic_refusal_is_atomic_and_retries_the_complete_muscl_interval() {
    let mut flow = chamber(gas(), 450.0);
    let volume = flow.mesh().volumes()[0];
    let controls = PistonControl {
        midpoint_velocity_bracket_m_s: [-0.01, 0.01],
        ..piston_controls()
    };
    let dt = 0.9 * muscl().cfl * volume / (0.07 * (flow.primitives()[0].sound_speed + 0.01));
    let conductance = 10.0 * volume * (internal_energy(10.0, 540.0) - internal_energy(10.0, 450.0))
        / (0.5 * dt * (1000.0 - 540.0));
    let first = chamber(gas(), 540.0);
    let second_temperature =
        reservoir_temperature(10.0, volume, 540.0, 1000.0, 0.5 * dt * conductance);
    let second = chamber(gas(), second_temperature);
    assert!(first.primitives()[0].sound_speed * 0.07 * dt / volume < muscl().cfl);
    assert!(second.primitives()[0].sound_speed * 0.07 * dt / volume > muscl().cfl);
    let mut clean = chamber(gas(), 450.0);
    let mut walls =
        RealGasWallOperator::new(&flow, vec![reservoir(1000.0, conductance); 3], || false).unwrap();
    let mut clean_walls =
        RealGasWallOperator::new(&clean, vec![reservoir(1000.0, conductance); 3], || false)
            .unwrap();
    let mut piston = Piston::new(1e10, 0.0).unwrap();
    let mut clean_piston = piston;
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 0.0,
    };
    let before = (
        flow.states().to_vec(),
        flow.primitives().to_vec(),
        flow.mesh().clone(),
        piston,
        walls.contacts().to_vec(),
    );
    let error = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(dt),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            load,
            controls,
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, PistonError::StepRetries { last, .. }
        if matches!(*last, PistonError::Gas(fs_gas::dynamics::Error::CflExceeded { .. }))));
    assert_eq!(flow.states(), before.0);
    assert_eq!(flow.primitives(), before.1);
    assert_eq!(flow.mesh(), &before.2);
    assert_eq!(piston, before.3);
    assert_eq!(walls.contacts(), before.4);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(walls.last_report(), None);
    let retry = flow
        .step_with_piston_and_walls(
            &mut piston,
            Time::new(dt),
            muscl(),
            load,
            controls,
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert!(retry.gas.halvings > 0);
    assert!(!retry.gas.first_order_retry);
    assert!(retry.gas.maximum_cfl <= muscl().cfl);
    clean
        .step_with_piston_and_walls(
            &mut clean_piston,
            retry.gas.elapsed,
            muscl(),
            load,
            controls,
            &mut clean_walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(flow.states(), clean.states());
    assert_eq!(flow.primitives(), clean.primitives());
    assert_eq!(flow.mesh(), clean.mesh());
    assert_eq!(piston, clean_piston);
    assert_eq!(walls.contacts(), clean_walls.contacts());
    assert_eq!(walls.last_report(), clean_walls.last_report());
}

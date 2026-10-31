//! Moving chamber work through real chemistry, molecular transport and finite walls.
use super::molecular::operator;
use super::*;
use crate::multicomponent::dissipation::{ConstantTransport, DissipationControl};
use crate::multicomponent::network::dissipation::NetworkDissipationOperator;
use crate::multicomponent::network::reactive::NetworkSources;
use crate::multicomponent::network::valve::{NetworkValveOperator, ValveSolve};
use crate::multicomponent::network::wall_exchange::NetworkWallOperator;
use crate::multicomponent::wall_exchange::{WallContact, WallControl, WallDrag, WallThermal};

type Molecular = NetworkDissipationOperator<ConstantTransport>;

fn finite_walls(raw: &MixtureNetwork1d) -> NetworkWallOperator<'static> {
    let thermal = WallThermal::Lumped {
        temperature_k: 0.8,
        heat_capacity_j_per_k: 1.0,
        conductance_w_per_k: 0.02,
    };
    NetworkWallOperator::new(
        raw,
        raw.branches()
            .iter()
            .map(|branch| {
                vec![
                    WallContact {
                        drag: WallDrag::Linear { rate_per_s: 0.1 },
                        thermal,
                    };
                    branch.cell_count()
                ]
            })
            .collect(),
        vec![
            WallContact {
                drag: WallDrag::None,
                thermal
            };
            raw.plenums().len()
        ],
        || false,
    )
    .unwrap()
}

fn fixture(hot: bool) -> (ReactingNetwork1d, Molecular, NetworkWallOperator<'static>) {
    let mechanism = if hot {
        mechanism(40.0, 500.0, 0.0)
    } else {
        mechanism(2.0, 2.0, 3.0)
    };
    let raw = network(&mechanism, 4, !hot, if hot { 1.0 } else { 1.4 });
    let molecular = operator(&raw, 0.01, 0.002, 0.03);
    let wall = finite_walls(&raw);
    let flow = ReactingNetwork1d::new(raw, mechanism, chemistry(), || false).unwrap();
    (flow, molecular, wall)
}

fn sources<'a>(wall: &'a mut NetworkWallOperator<'static>) -> NetworkSources<'a, 'static> {
    NetworkSources {
        wall: Some((wall, WallControl::default())),
        valves: None,
    }
}

fn run(
    flow: &mut ReactingNetwork1d,
    molecular: &mut Molecular,
    wall: &mut NetworkWallOperator<'static>,
    duration: f64,
    rate: f64,
) -> ReactiveNetworkReport {
    flow.step_with_dissipation_volume_rates_and_sources(
        molecular,
        Time::new(duration),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        &[rate],
        sources(wall),
        || false,
    )
    .unwrap()
}

fn fields(flow: &ReactingNetwork1d) -> Vec<Conserved> {
    let mut result = Vec::new();
    for branch in flow.branches() {
        for i in 0..branch.cell_count() {
            result.push(owned(branch.cell(i).unwrap()));
        }
    }
    result.extend(flow.plenums().iter().map(|p| owned(p.state())));
    result
}

fn wall_contacts(flow: &ReactingNetwork1d, wall: &NetworkWallOperator) -> Vec<WallContact> {
    let mut result = Vec::new();
    for i in 0..flow.branches().len() {
        result.extend_from_slice(wall.branch_contacts(i).unwrap());
    }
    result.extend_from_slice(wall.plenum_contacts());
    result
}

// Independent physical inventory; no coupled report or stored energy-gauge sum.
fn inventory(flow: &ReactingNetwork1d, wall: &NetworkWallOperator) -> (f64, f64) {
    let (mut mass, mut energy) = (0.0, 0.0);
    for branch in flow.branches() {
        for (i, &volume) in branch.mesh().volumes().iter().enumerate() {
            let q = branch.cell(i).unwrap();
            mass += volume * q.density();
            energy += volume * flow.gas().absolute_energy_density(q).unwrap();
        }
    }
    for plenum in flow.plenums() {
        mass += plenum.volume_m3() * plenum.state().density();
        energy += plenum.volume_m3() * flow.gas().absolute_energy_density(plenum.state()).unwrap();
    }
    for contact in wall_contacts(flow, wall) {
        match contact.thermal {
            WallThermal::Lumped {
                temperature_k,
                heat_capacity_j_per_k,
                ..
            } => {
                energy += heat_capacity_j_per_k * temperature_k;
            }
            _ => panic!("fixture requires independently reconstructible finite wall storage"),
        }
    }
    (mass, energy)
}

#[test]
fn reacting_compression_closes_physical_gas_wall_energy_with_live_molecular_and_heat_paths() {
    let (mut flow, mut molecular, mut wall) = fixture(false);
    let initial = inventory(&flow, &wall);
    let initial_volume = flow.plenums()[0].volume_m3();
    let mut total_work = 0.0;
    let mut chemical_conversion = 0.0;
    let mut molecular_exchange = 0.0_f64;
    for _ in 0..8 {
        let report = run(&mut flow, &mut molecular, &mut wall, 0.003, -0.2);
        let d = report.dissipation.unwrap();
        let w = report.wall.unwrap();
        assert_eq!(d.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(d.second.elapsed, d.first.elapsed);
        assert_eq!(w.first.elapsed, d.first.elapsed);
        assert_eq!(w.second.elapsed, d.first.elapsed);
        assert_eq!(
            report.boundary_work_on_gas_j,
            report.transport.boundary_work_on_gas_j
        );
        assert!(report.boundary_work_on_gas_j > 0.0);
        let motion = flow.accepted_plenum_boundary_work().unwrap()[0];
        assert_eq!(motion.final_volume_m3, flow.plenums()[0].volume_m3());
        assert_eq!(motion.work_on_gas_j, report.boundary_work_on_gas_j);
        close(
            motion.final_volume_m3 - motion.initial_volume_m3,
            -0.2 * report.elapsed.value(),
            1e-17,
        );
        // Walls exchange heat internally; chemistry converts physical stored
        // formation energy. Neither is a second mechanical-work supply.
        close(
            report.boundary_change.energy,
            w.first.heat_to_gas_j + w.second.heat_to_gas_j,
            1e-12,
        );
        close(report.conservation_residual.energy, 0.0, 3e-8);
        close(report.mass_residual_kg, 0.0, 3e-12);
        close(report.absolute_energy_residual_j, 0.0, 3e-8);
        assert!(
            report.conservation_residual.species_masses[..2]
                .iter()
                .all(|x| x.abs() < 3e-12)
        );
        assert!(w.first.heat_to_gas_j.abs() + w.second.heat_to_gas_j.abs() > 1e-7);
        chemical_conversion += report.chemistry.species_change_kg[1];
        molecular_exchange = molecular_exchange.max(
            flow.accepted_molecular_port_exchanges()
                .unwrap()
                .iter()
                .map(|p| p.energy_to_plenum.abs())
                .fold(0.0, f64::max),
        );
        total_work += report.boundary_work_on_gas_j;
    }
    assert!(chemical_conversion > 1e-5);
    assert!(molecular_exchange > 1e-8);
    close(
        flow.plenums()[0].volume_m3(),
        initial_volume - 0.2 * flow.time().value(),
        1e-15,
    );
    let final_inventory = inventory(&flow, &wall);
    close(final_inventory.0, initial.0, 3e-12);
    close(final_inventory.1 - initial.1, total_work, 8e-8);

    // Whole-step publication must adopt the operators' endpoint volume binding.
    // Their ordinary rigid entry points remain usable at that accepted geometry.
    let mut raw = flow.into_network();
    molecular
        .advance_for(
            &mut raw,
            Time::new(0.0),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    wall.advance_for(&mut raw, Time::new(0.0), WallControl::default(), || false)
        .unwrap();
    let changed = 1.01 * raw.plenums()[0].volume_m3();
    raw.quasistatic_adiabatic_plenum_to(0, changed, || false)
        .unwrap();
    assert!(
        molecular
            .advance_for(
                &mut raw,
                Time::new(0.0),
                DissipationControl::default(),
                || false
            )
            .is_err()
    );
    assert!(
        wall.advance_for(&mut raw, Time::new(0.0), WallControl::default(), || false)
            .is_err()
    );
}

#[test]
fn exothermic_moving_composition_restarts_fields_and_geometry_at_actual_duration() {
    let (mut flow, mut molecular, mut wall) = fixture(true);
    let initial_volume = flow.plenums()[0].volume_m3();
    let report = run(&mut flow, &mut molecular, &mut wall, 0.2, -0.3);
    assert!(report.coupling_retries > 0);
    close(
        flow.plenums()[0].volume_m3(),
        initial_volume - 0.3 * report.elapsed.value(),
        1e-16,
    );
    close(report.absolute_energy_residual_j, 0.0, 5e-8);
    let (mut reference, mut reference_molecular, mut reference_wall) = fixture(true);
    let expected = run(
        &mut reference,
        &mut reference_molecular,
        &mut reference_wall,
        report.elapsed.value(),
        -0.3,
    );
    assert_eq!(expected.elapsed, report.elapsed);
    assert_eq!(
        flow.plenums()[0].volume_m3(),
        reference.plenums()[0].volume_m3()
    );
    for (actual, expected) in fields(&flow).iter().zip(fields(&reference)) {
        for k in 0..2 {
            close(
                actual.partial_densities[k],
                expected.partial_densities[k],
                3e-9,
            );
        }
        close(actual.momentum, expected.momentum, 3e-9);
        close(actual.energy, expected.energy, 5e-8);
    }
    close(
        inventory(&flow, &wall).1,
        inventory(&reference, &reference_wall).1,
        5e-8,
    );
}

#[test]
fn final_cancellation_preserves_motion_operator_bindings_and_subsequent_restart_equivalence() {
    let (mut reference, mut reference_molecular, mut reference_wall) = fixture(false);
    run(
        &mut reference,
        &mut reference_molecular,
        &mut reference_wall,
        0.003,
        -0.2,
    );
    let mut polls = 0;
    reference
        .step_with_dissipation_volume_rates_and_sources(
            &mut reference_molecular,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            &[-0.2],
            sources(&mut reference_wall),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut molecular, mut wall) = fixture(false);
    run(&mut flow, &mut molecular, &mut wall, 0.003, -0.2);
    let original_fields = fields(&flow);
    let original_volumes = flow.plenums().to_vec();
    let original_time = flow.time();
    let original_report = flow.last_report();
    let original_work = flow.accepted_plenum_boundary_work().unwrap().to_vec();
    let original_ports = flow.accepted_molecular_port_exchanges().unwrap().to_vec();
    let original_wall = wall_contacts(&flow, &wall);
    let wall_report = wall.last_report();
    let molecular_report = molecular.last_report();
    let mut seen = 0;
    assert_eq!(
        flow.step_with_dissipation_volume_rates_and_sources(
            &mut molecular,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            &[-0.2],
            sources(&mut wall),
            || {
                seen += 1;
                seen >= polls
            },
        )
        .unwrap_err(),
        ReactiveNetworkError::Cancelled
    );
    assert_eq!(fields(&flow), original_fields);
    assert_eq!(flow.plenums(), original_volumes);
    assert_eq!(flow.time(), original_time);
    assert_eq!(flow.last_report(), original_report);
    assert_eq!(flow.accepted_plenum_boundary_work().unwrap(), original_work);
    assert_eq!(
        flow.accepted_molecular_port_exchanges().unwrap(),
        original_ports
    );
    assert_eq!(wall_contacts(&flow, &wall), original_wall);
    assert_eq!(wall.last_report(), wall_report);
    assert_eq!(molecular.last_report(), molecular_report);
    run(&mut flow, &mut molecular, &mut wall, 0.003, -0.2);
    assert_eq!(
        flow.plenums()[0].volume_m3(),
        reference.plenums()[0].volume_m3()
    );
    for (actual, expected) in fields(&flow).iter().zip(fields(&reference)) {
        for k in 0..2 {
            close(
                actual.partial_densities[k],
                expected.partial_densities[k],
                2e-10,
            );
        }
        close(actual.momentum, expected.momentum, 2e-10);
        close(actual.energy, expected.energy, 2e-9);
    }
    close(
        inventory(&flow, &wall).1,
        inventory(&reference, &reference_wall).1,
        2e-9,
    );
}

#[test]
fn zero_motion_reuses_rigid_sources_and_nonzero_motion_with_valves_refuses() {
    let (mut rigid, mut rigid_molecular, mut rigid_wall) = fixture(false);
    let (mut moving, mut moving_molecular, mut moving_wall) = fixture(false);
    let expected = rigid
        .step_with_dissipation_and_sources(
            &mut rigid_molecular,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            sources(&mut rigid_wall),
            || false,
        )
        .unwrap();
    let actual = run(
        &mut moving,
        &mut moving_molecular,
        &mut moving_wall,
        0.003,
        0.0,
    );
    assert_eq!(actual, expected);
    assert_eq!(fields(&moving), fields(&rigid));
    assert_eq!(
        wall_contacts(&moving, &moving_wall),
        wall_contacts(&rigid, &rigid_wall)
    );
    assert_eq!(
        moving_molecular.last_report(),
        rigid_molecular.last_report()
    );

    let mechanism = mechanism(0.0, 2.0, 0.0);
    let raw = network(&mechanism, 2, false, 1.0);
    let mut valves = NetworkValveOperator::new(&raw, vec![], || false).unwrap();
    let mut flow = ReactingNetwork1d::new(raw, mechanism, chemistry(), || false).unwrap();
    let old = fields(&flow);
    let old_volumes = flow.plenums().to_vec();
    assert!(matches!(
        flow.step_with_volume_rates_and_sources(
            Time::new(0.001),
            ReactiveStepControl::default(),
            &[-0.1],
            NetworkSources {
                wall: None,
                valves: Some((&mut valves, ValveSolve::default()))
            },
            || false,
        ),
        Err(ReactiveNetworkError::Valve(_))
    ));
    assert_eq!(fields(&flow), old);
    assert_eq!(flow.plenums(), old_volumes);
    assert_eq!(flow.time().value(), 0.0);
    assert!(valves.last_report().is_none());
    let accepted = flow
        .step_with_volume_rates_and_sources(
            Time::new(0.001),
            ReactiveStepControl::default(),
            &[0.0],
            NetworkSources {
                wall: None,
                valves: Some((&mut valves, ValveSolve::default())),
            },
            || false,
        )
        .unwrap();
    assert_eq!(accepted.boundary_work_on_gas_j, 0.0);
    assert!(accepted.valves.is_some());
}

//! Complete wall/molecular/wave conservation and atomic coupled publication.
use super::*;
use crate::multicomponent::dissipation::ConstantTransport;
use crate::multicomponent::network::dissipating::{
    DissipatingNetwork1d, DissipatingNetworkError, DissipatingStepControl,
};
use crate::multicomponent::network::dissipation::tests::{law, owned};
use crate::multicomponent::network::{BranchEnd, PlenumConnection};
use crate::multicomponent::{Boundary, Conserved};

type Flow = DissipatingNetwork1d<ConstantTransport>;

fn coupled(g: &MixtureGas, contact: WallContact) -> (Flow, NetworkWallOperator<'static>) {
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&[1.0], 0.0, 2.0).unwrap(); 2],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    ).unwrap();
    let n = MixtureNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.1, g.conserved(&[2.0], 0.0, 4.0).unwrap()).unwrap()],
        vec![PlenumConnection { branch: 0, end: BranchEnd::Right, plenum: 0 }],
        || false,
    ).unwrap();
    let wall = uniform_operator(&n, contact);
    let flow = DissipatingNetwork1d::new(
        n,
        law(g, 0.004, 0.0, 0.02, 0.0),
        vec![Length::new(0.25)],
        || false,
    ).unwrap();
    (flow, wall)
}
fn states(flow: &Flow) -> (Vec<Vec<Conserved>>, Vec<Conserved>) {
    (
        flow.branches().iter()
            .map(|b| (0..b.cell_count()).map(|i| owned(b.cell(i).unwrap())).collect())
            .collect(),
        flow.plenums().iter().map(|p| owned(p.state())).collect(),
    )
}
// Integrate physical gas and finite-wall storage directly from fields.
fn inventory(flow: &Flow, wall: &NetworkWallOperator) -> (f64, f64) {
    let mut mass = 0.0;
    let mut energy = 0.0;
    for branch in flow.branches() {
        for (i, volume) in branch.mesh().volumes().iter().enumerate() {
            let q = branch.cell(i).unwrap();
            mass += volume * q.density();
            energy += volume * flow.gas().absolute_energy_density(q).unwrap();
        }
    }
    for plenum in flow.plenums() {
        mass += plenum.volume_m3() * plenum.state().density();
        energy += plenum.volume_m3() * flow.gas().absolute_energy_density(plenum.state()).unwrap();
    }
    for contact in wall.branch_contacts(0).unwrap().iter().chain(wall.plenum_contacts()) {
        if let WallThermal::Lumped { temperature_k, heat_capacity_j_per_k, .. } = contact.thermal {
            energy += heat_capacity_j_per_k * temperature_k;
        }
    }
    (mass, energy)
}

#[test]
fn coupled_finite_walls_close_storage_and_ports_and_zero_walls_retain_original_step() {
    let g = gas(1, false, 0.0);
    let (mut flow, mut wall) = coupled(&g, lumped());
    let initial = inventory(&flow, &wall);
    let mut maximum_advective_energy = 0.0_f64;
    let mut maximum_molecular_energy = 0.0_f64;
    for _ in 0..8 {
        let plenum_before = owned(flow.plenums()[0].state());
        let wall_before = wall_temperature(wall.plenum_contacts()[0]);
        let report = flow.step_with_wall(
            &mut wall, Time::new(0.01), DissipatingStepControl::default(),
            controls(), || false,
        ).unwrap();
        let halves = report.wall.unwrap();
        assert_eq!(halves.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(halves.second.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(wall.last_report().unwrap().elapsed, report.elapsed);
        assert_eq!(flow.last_report(), Some(report));
        close(
            report.boundary_change.energy,
            halves.first.heat_to_gas_j + halves.second.heat_to_gas_j,
            3e-12,
        );
        assert_eq!(wall.last_report().unwrap().external_heat_to_gas_j, 0.0);
        close(report.absolute_energy_residual_j, 0.0, 3e-11);
        let advective = flow.accepted_port_exchanges().unwrap()[0];
        let molecular = flow.accepted_molecular_port_exchanges().unwrap()[0];
        // The plenum wall has Cwall=1 J/K: use actual wall storage independently.
        let wall_heat = wall_before - wall_temperature(wall.plenum_contacts()[0]);
        let plenum = &flow.plenums()[0];
        close(
            plenum.volume_m3() * (plenum.state().energy - plenum_before.energy),
            advective.energy_to_plenum + molecular.energy_to_plenum + wall_heat,
            3e-11,
        );
        maximum_advective_energy = maximum_advective_energy.max(advective.energy_to_plenum.abs());
        maximum_molecular_energy = maximum_molecular_energy.max(molecular.energy_to_plenum.abs());
        let current = inventory(&flow, &wall);
        close(current.0, initial.0, 3e-12);
        close(current.1, initial.1, 3e-10);
    }
    assert!(maximum_advective_energy > 1e-8);
    assert!(maximum_molecular_energy > 0.0);
    let zero = WallContact { drag: WallDrag::None, thermal: WallThermal::Adiabatic };
    let (mut plain, _) = coupled(&g, zero);
    let (mut with_wall, mut zero_wall) = coupled(&g, zero);
    let report = with_wall.step_with_wall(
        &mut zero_wall, Time::new(0.01), DissipatingStepControl::default(),
        controls(), || false,
    ).unwrap();
    plain.step(report.elapsed, DissipatingStepControl::default(), || false).unwrap();
    let actual = states(&with_wall);
    let expected = states(&plain);
    for (a, b) in actual.0.iter().flatten().chain(&actual.1)
        .zip(expected.0.iter().flatten().chain(&expected.1))
    {
        close(a.partial_densities[0], b.partial_densities[0], 3e-12);
        close(a.momentum, b.momentum, 3e-12);
        close(a.energy, b.energy, 3e-11);
    }
    assert_eq!(with_wall.time(), plain.time());
    assert_eq!(zero_wall.last_report().unwrap().heat_to_gas_j, 0.0);
    assert_eq!(zero_wall.last_report().unwrap().drag_heating_j, 0.0);
}

#[test]
fn late_coupled_cancellation_preserves_gas_walls_and_both_port_reports_then_reuses_work() {
    let g = gas(1, false, 0.0);
    let make = || {
        let (mut flow, mut wall) = coupled(&g, lumped());
        flow.step_with_wall(
            &mut wall, Time::new(0.01), DissipatingStepControl::default(),
            controls(), || false,
        ).unwrap();
        (flow, wall)
    };
    let (mut reference, mut reference_wall) = make();
    let mut polls = 0;
    reference.step_with_wall(
        &mut reference_wall, Time::new(0.01), DissipatingStepControl::default(),
        controls(), || { polls += 1; false },
    ).unwrap();
    assert!(polls > 1);
    let (mut flow, mut wall) = make();
    let before = states(&flow);
    let clock = flow.time();
    let report = flow.last_report();
    let advective = flow.accepted_port_exchanges().unwrap().to_vec();
    let molecular = flow.accepted_molecular_port_exchanges().unwrap().to_vec();
    let branches = flow.accepted_branch_transport_reports().unwrap().to_vec();
    let wall_report = wall.last_report();
    let branch_contacts = wall.branch_contacts(0).unwrap().to_vec();
    let plenum_contacts = wall.plenum_contacts().to_vec();
    let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
    let mut visits = 0;
    let result = flow.step_with_wall(
        &mut wall, Time::new(0.01), DissipatingStepControl::default(),
        controls(), || { visits += 1; visits >= polls },
    );
    assert_eq!(result, Err(DissipatingNetworkError::Cancelled));
    assert_eq!(states(&flow), before);
    assert_eq!(flow.time(), clock);
    assert_eq!(flow.last_report(), report);
    assert_eq!(flow.accepted_port_exchanges().unwrap(), advective);
    assert_eq!(flow.accepted_molecular_port_exchanges().unwrap(), molecular);
    assert_eq!(flow.accepted_branch_transport_reports().unwrap(), branches);
    assert_eq!(wall.last_report(), wall_report);
    assert_eq!(wall.branch_contacts(0).unwrap(), branch_contacts);
    assert_eq!(wall.plenum_contacts(), plenum_contacts);
    assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
    flow.step_with_wall(
        &mut wall, Time::new(0.01), DissipatingStepControl::default(),
        controls(), || false,
    ).unwrap();
    let actual = states(&flow);
    let expected = states(&reference);
    for (a, b) in actual.0.iter().flatten().chain(&actual.1)
        .zip(expected.0.iter().flatten().chain(&expected.1))
    {
        close(a.partial_densities[0], b.partial_densities[0], 3e-12);
        close(a.momentum, b.momentum, 3e-12);
        close(a.energy, b.energy, 3e-11);
    }
    for (a, b) in wall.branch_contacts(0).unwrap().iter().chain(wall.plenum_contacts())
        .zip(reference_wall.branch_contacts(0).unwrap().iter().chain(reference_wall.plenum_contacts()))
    {
        close(wall_temperature(*a), wall_temperature(*b), 3e-11);
    }
    assert_eq!(flow.time(), reference.time());
}

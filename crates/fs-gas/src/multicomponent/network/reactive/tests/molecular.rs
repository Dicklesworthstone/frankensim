//! Independent combined diffusion/reaction oracle and coupled publication checks.
use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, TransportCoefficients, TransportDomain,
};
use crate::multicomponent::network::dissipation::NetworkDissipationOperator;
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn operator(
    network: &MixtureNetwork1d,
    diffusion: f64,
    viscosity: f64,
    conductivity: f64,
) -> NetworkDissipationOperator<ConstantTransport> {
    let law = ConstantTransport::new(
        network.gas().source().clone(),
        "reacting-network-constant-transport-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e5,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(viscosity),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        vec![0.0, diffusion, diffusion, 0.0],
    )
    .unwrap();
    NetworkDissipationOperator::new(
        network,
        law,
        vec![Length::new(0.5); network.connections().len()],
        || false,
    )
    .unwrap()
}

#[test]
fn molecular_reacting_network_matches_two_pool_analytic_diffusion_reaction() {
    let rate = 2.0;
    let diffusion = 0.1;
    let mechanism = mechanism(0.0, rate, 0.0);
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    // Two identical branch cells each have one plenum port. Symmetry leaves
    // them identical: V_cell=0.1, V_plenum=0.2, A=0.1, total distance=1.
    // Both species have R_specific=1, Cv=2.5 and identical formation energy.
    // Advection vanishes, so the composition contrast decays at 2*k + 2*D.
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(1.0),
        0.1,
        vec![gas.conserved(&[0.9, 0.1], 0.0, 1.0).unwrap(); 2],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let network = MixtureNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.2, gas.conserved(&[0.1, 0.9], 0.0, 1.0).unwrap()).unwrap()],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Left,
                plenum: 0,
            },
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap();
    let mut molecular = operator(&network, diffusion, 0.0, 0.0);
    let mut flow = ReactingNetwork1d::new(network, mechanism, chemistry(), || false).unwrap();
    for _ in 0..4 {
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(0.05),
                ReactiveStepControl::default(),
                DissipationControl {
                    maximum_step_s: 0.002,
                    ..DissipationControl::default()
                },
                || false,
            )
            .unwrap();
        let t = flow.time().value();
        let contrast = 0.4 * (-(2.0 * rate + 2.0 * diffusion) * t).exp();
        for cell in 0..2 {
            let q = flow.branches()[0].cell(cell).unwrap();
            close(q.partial_densities[0], 0.5 + contrast, 3e-5);
            let p = flow.gas().primitive(q, || false).unwrap();
            close(p.temperature, 1.0, 2e-8);
            close(p.velocity, 0.0, 2e-8);
        }
        close(
            flow.plenums()[0].state().partial_densities[0],
            0.5 - contrast,
            3e-5,
        );
        let d = report.dissipation.unwrap();
        assert_eq!(d.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(d.second.elapsed, d.first.elapsed);
        assert!(d.first.accepted_steps > 1 && d.second.accepted_steps > 1);
        close(report.mass_residual_kg, 0.0, 3e-13);
        close(report.absolute_energy_residual_j, 0.0, 3e-12);
        assert!(
            flow.accepted_molecular_port_exchanges()
                .unwrap()
                .iter()
                .all(|port| port.species_masses_to_plenum[0] > 0.0)
        );
    }
}

#[test]
fn hot_molecular_reacting_network_restarts_all_five_stages_at_actual_duration() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let initial = network(&mechanism, 4, false, 1.0);
    let mut molecular = operator(&initial, 0.1, 0.0, 0.0);
    let mut flow =
        ReactingNetwork1d::new(initial, mechanism.clone(), chemistry(), || false).unwrap();
    let report = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.coupling_retries > 0);
    let d = report.dissipation.unwrap();
    assert_eq!(d.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(d.second.elapsed, d.first.elapsed);
    assert_eq!(molecular.last_report(), Some(d.second));
    // No molecular gradients exist in this homogeneous fixture. The independently
    // owned chemistry/network path must reach the same actual-duration endpoint.
    let mut reference = reacting(&mechanism, 4, false);
    let expected = reference
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    close(report.elapsed.value(), expected.elapsed.value(), 2e-11);
    for (a, b) in flow.branches().iter().zip(reference.branches()) {
        for i in 0..a.cell_count() {
            let qa = a.cell(i).unwrap();
            let qb = b.cell(i).unwrap();
            for k in 0..2 {
                close(qa.partial_densities[k], qb.partial_densities[k], 2e-7);
            }
            close(qa.energy, qb.energy, 2e-6);
        }
    }
    for port in flow.accepted_molecular_port_exchanges().unwrap() {
        close(port.mass_to_plenum(), 0.0, 2e-14);
        close(
            port.absolute_energy_to_plenum(flow.gas()).unwrap(),
            0.0,
            2e-10,
        );
    }
    close(report.absolute_energy_residual_j, 0.0, 4e-11);
}

fn snapshot(flow: &ReactingNetwork1d) -> Vec<Conserved> {
    let mut fields = Vec::new();
    for branch in flow.branches() {
        for i in 0..branch.cell_count() {
            fields.push(owned(branch.cell(i).unwrap()));
        }
    }
    fields.extend(flow.plenums().iter().map(|p| owned(p.state())));
    fields
}

#[test]
fn late_five_stage_cancellation_preserves_operator_and_network_publication() {
    let make = || {
        let mechanism = mechanism(1.0, 2.0, 0.0);
        let network = network(&mechanism, 3, true, 1.4);
        let operator = operator(&network, 0.01, 0.002, 0.03);
        let flow = ReactingNetwork1d::new(network, mechanism, chemistry(), || false).unwrap();
        (flow, operator)
    };
    let run = |flow: &mut ReactingNetwork1d,
               op: &mut NetworkDissipationOperator<ConstantTransport>| {
        flow.step_with_dissipation(
            op,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap()
    };
    let (mut reference, mut reference_op) = make();
    run(&mut reference, &mut reference_op);
    let mut polls = 0usize;
    reference
        .step_with_dissipation(
            &mut reference_op,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 20);
    let (mut flow, mut molecular) = make();
    run(&mut flow, &mut molecular);
    let fields = snapshot(&flow);
    let time = flow.time();
    let report = flow.last_report();
    let advection = flow.accepted_port_exchanges().unwrap().to_vec();
    let diffusion = flow.accepted_molecular_port_exchanges().unwrap().to_vec();
    let op_report = molecular.last_report();
    let op_ports = molecular.accepted_port_exchanges().unwrap().to_vec();
    let mut current = 0usize;
    let error = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.003),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || {
                current += 1;
                current >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveNetworkError::Cancelled);
    assert_eq!(snapshot(&flow), fields);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_report(), report);
    assert_eq!(flow.accepted_port_exchanges().unwrap(), advection);
    assert_eq!(flow.accepted_molecular_port_exchanges().unwrap(), diffusion);
    assert_eq!(molecular.last_report(), op_report);
    assert_eq!(molecular.accepted_port_exchanges().unwrap(), op_ports);
    run(&mut flow, &mut molecular);
    for (actual, expected) in snapshot(&flow).iter().zip(snapshot(&reference)) {
        for k in 0..2 {
            close(
                actual.partial_densities[k],
                expected.partial_densities[k],
                2e-7,
            );
        }
        close(actual.momentum, expected.momentum, 2e-7);
        close(actual.energy, expected.energy, 2e-6);
    }
}

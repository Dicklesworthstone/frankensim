//! Physical bath exchange and whole-composition RKL2 work/rollback.
use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use crate::multicomponent::network::dissipation::tests::reservoir::{bath_network, bath_step};

fn bath_owner(
    gas: &MixtureGas,
) -> DissipatingNetwork1d<crate::multicomponent::dissipation::ConstantTransport> {
    DissipatingNetwork1d::new(
        bath_network(gas, false),
        law(gas, 0.0, 0.0, 0.0, 0.03),
        vec![Length::new(0.15)],
        || false,
    )
    .unwrap()
}

#[test]
fn g1_stabilized_baths_cross_both_halves_and_close_finite_plenum_inventory() {
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let mut flow = bath_owner(&gas);
    let before = inventories(&flow);
    let report = flow
        .step_stabilized_with_sources(
            Time::new(0.04),
            StabilizedDissipatingStepControl::default(),
            Some(&[[Some(bath.as_cell()), None]]),
            NetworkSources::default(),
            || false,
        )
        .unwrap();
    // At s=2 the Legendre stability polynomial is 1+z+z²/2. This
    // independently assembled three-storage ODE includes the exterior
    // half-cell and finite plenum; the isobaric advection contributes zero.
    let half = 0.5 * report.coupled.elapsed.value();
    let (first, incoming_a) = bath_step([0.4, 0.4, 0.1], half);
    let (expected, incoming_b) = bath_step(first, half);
    assert_eq!(report.dissipation.first.maximum_stages_used, 2);
    assert_eq!(report.dissipation.second.maximum_stages_used, 2);
    assert_eq!(report.coupled.coupling_retries, 0);
    assert_eq!(
        report.operator_evaluations,
        report.dissipation.first.operator_evaluations
            + report.dissipation.second.operator_evaluations
    );
    assert_eq!(
        report.dissipation.first.transport,
        report.coupled.molecular.first
    );
    assert_eq!(
        report.dissipation.second.transport,
        report.coupled.molecular.second
    );
    for (i, &value) in expected[..2].iter().enumerate() {
        close(
            flow.branches()[0].cell(i).unwrap().partial_densities[0],
            value,
            3e-11,
        );
    }
    close(
        flow.plenums()[0].state().partial_densities[0],
        expected[2],
        3e-11,
    );
    close(
        report.coupled.boundary_change.species_masses[0],
        incoming_a + incoming_b,
        3e-12,
    );
    let after = inventories(&flow);
    for k in 0..2 {
        close(
            after.0[k] - before.0[k],
            report.coupled.boundary_change.species_masses[k],
            3e-12,
        );
    }
    let advective = flow.accepted_port_exchanges().unwrap()[0];
    let molecular = flow.accepted_molecular_port_exchanges().unwrap()[0];
    close(
        0.07 * (expected[2] - 0.1),
        advective.species_masses_to_plenum[0] + molecular.species_masses_to_plenum[0],
        3e-12,
    );
    close(
        after.1 - before.1,
        report
            .coupled
            .boundary_change
            .absolute_energy(&gas)
            .unwrap(),
        3e-11,
    );
    assert_eq!(flow.last_report(), Some(report.coupled));
    assert_eq!(flow.time(), report.coupled.elapsed);
}

#[test]
fn g4_stabilized_network_budget_and_late_cancellation_publish_no_partial_half() {
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let reservoirs = [[Some(bath.as_cell()), None]];
    let controls = StabilizedDissipatingStepControl::default();
    let run = |flow: &mut DissipatingNetwork1d<_>| {
        flow.step_stabilized_with_sources(
            Time::new(0.04),
            controls,
            Some(&reservoirs),
            NetworkSources::default(),
            || false,
        )
        .unwrap()
    };
    let make = || {
        let mut flow = bath_owner(&gas);
        run(&mut flow);
        flow
    };
    let mut reference = make();
    let mut polls = 0;
    let reference_report = reference
        .step_stabilized_with_sources(
            Time::new(0.04),
            controls,
            Some(&reservoirs),
            NetworkSources::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = make();
    let original = snapshot(&flow);
    let original_operator = flow.molecular.last_report();
    let original_ports = flow.molecular.accepted_port_exchanges().unwrap().to_vec();
    let first_work = reference_report.dissipation.first.operator_evaluations;
    let error = flow
        .step_stabilized_with_sources(
            Time::new(0.04),
            StabilizedDissipatingStepControl {
                molecular: StabilizedDissipationControl {
                    maximum_operator_evaluations: first_work,
                    ..controls.molecular
                },
                ..controls
            },
            Some(&reservoirs),
            NetworkSources::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,
        DissipatingNetworkError::Molecular(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == first_work));
    assert_eq!(snapshot(&flow), original);
    assert_eq!(flow.molecular.last_report(), original_operator);
    assert_eq!(
        flow.molecular.accepted_port_exchanges().unwrap(),
        original_ports
    );
    let mut cancelled = make();
    let original = snapshot(&cancelled);
    let mut visits = 0;
    assert_eq!(
        cancelled.step_stabilized_with_sources(
            Time::new(0.04),
            controls,
            Some(&reservoirs),
            NetworkSources::default(),
            || {
                visits += 1;
                visits >= polls
            },
        ),
        Err(DissipatingNetworkError::Cancelled)
    );
    assert_eq!(snapshot(&cancelled), original);
    for candidate in [&mut flow, &mut cancelled] {
        let report = run(candidate);
        assert_eq!(
            report.operator_evaluations,
            reference_report.operator_evaluations
        );
        let actual = fields(&candidate.accepted);
        let expected = fields(&reference.accepted);
        for (a, b) in actual
            .0
            .iter()
            .flatten()
            .chain(actual.1.iter())
            .zip(expected.0.iter().flatten().chain(expected.1.iter()))
        {
            for k in 0..2 {
                close(a.partial_densities[k], b.partial_densities[k], 3e-11);
            }
            close(a.momentum, b.momentum, 3e-11);
            close(a.energy, b.energy, 3e-10);
        }
    }
}

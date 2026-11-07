//! Source entropy and atomic work budgets in the complete RKL reacting owner.
use super::*;
use crate::multicomponent::dissipation::{DissipationError, StabilizedDissipationControl};

#[test]
fn g1_rkl_composition_keeps_variable_calorics_and_physical_chemical_entropy() {
    let (mut flow, mut molecular) = variable_caloric();
    let calorics = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
    let stripped = [[3.5, 0.2, 0.0], [4.0, 0.1, 0.0]];
    let initial = source_entropy(&flow, &calorics);
    let initial_stripped = source_entropy(&flow, &stripped);
    let before = flow.totals();
    let report = flow
        .step_with_stabilized_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.002),
            controls(),
            StabilizedDissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let admitted = report.coupled;
    let step = admitted.coupled;
    let physical = step.dissipation.unwrap();
    let entropy = admitted.dissipation_entropy.unwrap();
    assert_eq!(step.elapsed.value(), 0.002);
    assert_eq!(flow.time(), step.elapsed);
    assert_eq!(report.dissipation.first.transport, physical.first);
    assert_eq!(report.dissipation.second.transport, physical.second);
    for half in [report.dissipation.first, report.dissipation.second] {
        assert_eq!(half.transport.elapsed.value(), 0.5 * step.elapsed.value());
        assert!(half.maximum_stages_used >= 2);
        assert!(half.maximum_scaled_parabolic_number <= 0.40000000000001);
    }
    assert_eq!(
        report.operator_evaluations,
        report.dissipation.first.operator_evaluations
            + report.dissipation.second.operator_evaluations
    );
    assert!(entropy.first.entropy_change_j_k > 0.0);
    assert!(entropy.second.entropy_change_j_k > 0.0);
    assert!(admitted.transport_entropy.entropy_change_j_k >= -1e-10);
    let actual = source_entropy(&flow, &calorics) - initial;
    let without_datum = source_entropy(&flow, &stripped) - initial_stripped;
    close(
        actual - without_datum,
        1.5 * (step.after.species_masses[1] - before.species_masses[1]),
        3e-13,
    );
    close(admitted.entropy.source_entropy_change_j_k, actual, 3e-10);
    close(
        admitted.entropy.entropy_change_j_k,
        actual - admitted.entropy.reference_entropy_residual_j_k,
        3e-10,
    );
    assert!((admitted.entropy.entropy_change_j_k - without_datum).abs() > 1e-4);
    assert!(admitted.entropy.maximum_relative_species_residual < 1e-11);
    close(step.mass_residual_kg, 0.0, 3e-12);
    close(step.absolute_energy_residual_j, 0.0, 5e-11);
    assert_eq!(molecular.last_report(), Some(physical.second));
    assert!(
        flow.accepted_molecular_port_exchanges()
            .unwrap()
            .iter()
            .any(|port| port.species_masses_to_plenum[0] != 0.0)
    );
}

#[test]
fn g4_heat_release_retry_cannot_reset_stabilized_source_entropy_work_allowance() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let make = || {
        let raw = network(&mechanism, 4, false, 1.0);
        let molecular = super::super::molecular::operator(&raw, 0.1, 0.0, 0.0);
        let flow = ReactingNetwork1d::new(raw, mechanism.clone(), chemistry(), || false).unwrap();
        (flow, molecular)
    };
    let (mut reference, mut reference_op) = make();
    let initial = source_entropy(&reference, &[[3.5, 0.0, 0.0]; 2]);
    let report = reference
        .step_with_stabilized_dissipation_and_entropy(
            &mut reference_op,
            Time::new(0.2),
            controls(),
            StabilizedDissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let accepted_work = report.dissipation.first.operator_evaluations
        + report.dissipation.second.operator_evaluations;
    assert!(report.coupled.coupled.coupling_retries > 0);
    assert!(report.operator_evaluations > accepted_work);
    assert!(report.coupled.coupled.after.species_masses[1] > 0.01);
    close(
        report.coupled.entropy.source_entropy_change_j_k,
        source_entropy(&reference, &[[3.5, 0.0, 0.0]; 2]) - initial,
        3e-9,
    );
    let (mut flow, mut molecular) = make();
    let before = snapshot(&flow);
    let error = flow
        .step_with_stabilized_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.2),
            controls(),
            StabilizedDissipationControl {
                maximum_operator_evaluations: accepted_work,
                ..Default::default()
            },
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,
        ReactiveNetworkError::Dissipation(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == accepted_work));
    assert_eq!(snapshot(&flow), before);
    assert!(molecular.last_report().is_none());
    let replay = flow
        .step_with_stabilized_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.2),
            controls(),
            StabilizedDissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(replay.operator_evaluations, report.operator_evaluations);
    close(
        replay.coupled.coupled.elapsed.value(),
        report.coupled.coupled.elapsed.value(),
        3e-12,
    );
    close(
        replay.coupled.entropy.source_entropy_change_j_k,
        report.coupled.entropy.source_entropy_change_j_k,
        3e-9,
    );
    close(
        replay.coupled.coupled.absolute_energy_residual_j,
        0.0,
        5e-10,
    );
}

#[test]
fn g4_final_rkl_cancellation_keeps_accepted_states_ports_and_both_entropy_receipts() {
    let entropy = MixtureEntropyControl::default();
    let rkl = StabilizedDissipationControl::default();
    let (mut probe, mut probe_op) = variable_caloric();
    probe
        .step_with_stabilized_dissipation_and_entropy(
            &mut probe_op,
            Time::new(0.001),
            controls(),
            rkl,
            entropy,
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = probe
        .step_with_stabilized_dissipation_and_entropy(
            &mut probe_op,
            Time::new(0.002),
            controls(),
            rkl,
            entropy,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut molecular) = variable_caloric();
    flow.step_with_stabilized_dissipation_and_entropy(
        &mut molecular,
        Time::new(0.001),
        controls(),
        rkl,
        entropy,
        || false,
    )
    .unwrap();
    let before = snapshot(&flow);
    let molecular_ports = flow.accepted_molecular_port_exchanges().unwrap().to_vec();
    let operator_report = molecular.last_report();
    let operator_ports = molecular.accepted_port_exchanges().unwrap().to_vec();
    let operator_branches = molecular.accepted_branch_reports().unwrap().to_vec();
    let mut visits = 0;
    let result = flow.step_with_stabilized_dissipation_and_entropy(
        &mut molecular,
        Time::new(0.002),
        controls(),
        rkl,
        entropy,
        || {
            visits += 1;
            visits == polls
        },
    );
    assert_eq!(result, Err(ReactiveNetworkError::Cancelled));
    assert_eq!(snapshot(&flow), before);
    assert_eq!(
        flow.accepted_molecular_port_exchanges().unwrap(),
        molecular_ports
    );
    assert_eq!(molecular.last_report(), operator_report);
    assert_eq!(molecular.accepted_port_exchanges().unwrap(), operator_ports);
    assert_eq!(
        molecular.accepted_branch_reports().unwrap(),
        operator_branches
    );
    let recovered = flow
        .step_with_stabilized_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.002),
            controls(),
            rkl,
            entropy,
            || false,
        )
        .unwrap();
    assert_eq!(flow.time(), probe.time());
    let actual = fields(&flow.accepted);
    let reference = fields(&probe.accepted);
    for (a, b) in actual
        .0
        .iter()
        .flatten()
        .chain(&actual.1)
        .zip(reference.0.iter().flatten().chain(&reference.1))
    {
        for k in 0..2 {
            close(a.partial_densities[k], b.partial_densities[k], 3e-10);
        }
        close(a.momentum, b.momentum, 3e-10);
        close(a.energy, b.energy, 3e-9);
    }
    close(
        recovered.coupled.entropy.source_entropy_change_j_k,
        expected.coupled.entropy.source_entropy_change_j_k,
        3e-9,
    );
    close(
        recovered.coupled.coupled.absolute_energy_residual_j,
        0.0,
        5e-11,
    );
}

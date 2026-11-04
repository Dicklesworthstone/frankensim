//! Independently solved reaction/diffusion decay and coupled RKL2 work.
use super::molecular::operator;
use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, StabilizedDissipationControl,
};
use crate::multicomponent::network::dissipation::NetworkDissipationOperator;

fn two_pool(
    rate: f64,
    diffusion: f64,
) -> (
    ReactingNetwork1d,
    NetworkDissipationOperator<ConstantTransport>,
) {
    let mechanism = mechanism(0.0, rate, 0.0);
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
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
    let raw = MixtureNetwork1d::new(
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
    let molecular = operator(&raw, diffusion, 0.0, 0.0);
    let flow = ReactingNetwork1d::new(raw, mechanism, chemistry(), || false).unwrap();
    (flow, molecular)
}

#[test]
fn g1_stabilized_reaction_diffusion_tracks_two_pool_continuous_decay() {
    let rate = 2.0;
    let diffusion = 10.0;
    let (mut flow, mut molecular) = two_pool(rate, diffusion);
    let controls = StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: 0.002,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    };
    // The two identical branch cells and the plenum have equal total volumes.
    // Their chemical/diffusive contrast obeys d(delta)/dt = -(2k+2D)delta.
    // No native stage recurrence or balance report enters this continuum oracle.
    for _ in 0..4 {
        let report = flow
            .step_with_stabilized_dissipation(
                &mut molecular,
                Time::new(0.05),
                ReactiveStepControl::default(),
                controls,
                || false,
            )
            .unwrap();
        let contrast = 0.4 * (-(2.0 * rate + 2.0 * diffusion) * flow.time().value()).exp();
        for cell in 0..2 {
            close(
                flow.branches()[0].cell(cell).unwrap().partial_densities[0],
                0.5 + contrast,
                5e-5,
            );
        }
        close(
            flow.plenums()[0].state().partial_densities[0],
            0.5 - contrast,
            5e-5,
        );
        assert_eq!(report.coupled.coupling_retries, 0);
        assert_eq!(
            report.operator_evaluations,
            report.dissipation.first.operator_evaluations
                + report.dissipation.second.operator_evaluations
        );
        let physical = report.coupled.dissipation.unwrap();
        assert_eq!(physical.first, report.dissipation.first.transport);
        assert_eq!(physical.second, report.dissipation.second.transport);
        assert_eq!(
            physical.first.elapsed.value(),
            0.5 * report.coupled.elapsed.value()
        );
        assert_eq!(physical.second.elapsed, physical.first.elapsed);
        close(report.coupled.mass_residual_kg, 0.0, 3e-12);
        close(report.coupled.absolute_energy_residual_j, 0.0, 3e-11);
        assert_eq!(molecular.last_report(), Some(physical.second));
    }
}

#[test]
fn g4_heat_release_retries_share_actual_rkl2_work_and_restore_failed_state() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let make = || {
        let raw = network(&mechanism, 4, false, 1.0);
        let molecular = operator(&raw, 0.1, 0.0, 0.0);
        let flow = ReactingNetwork1d::new(raw, mechanism.clone(), chemistry(), || false).unwrap();
        (flow, molecular)
    };
    let (mut reference, mut op) = make();
    let controls = StabilizedDissipationControl::default();
    let report = reference
        .step_with_stabilized_dissipation(
            &mut op,
            Time::new(0.2),
            ReactiveStepControl::default(),
            controls,
            || false,
        )
        .unwrap();
    assert!(report.coupled.coupling_retries > 0);
    let accepted_work = report.dissipation.first.operator_evaluations
        + report.dissipation.second.operator_evaluations;
    assert!(report.operator_evaluations > accepted_work);
    let (mut flow, mut molecular) = make();
    let before = super::super::super::dissipation::tests::fields(&flow.accepted);
    let error = flow
        .step_with_stabilized_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            StabilizedDissipationControl {
                maximum_operator_evaluations: accepted_work,
                ..controls
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,
        ReactiveNetworkError::Dissipation(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == accepted_work));
    assert_eq!(
        super::super::super::dissipation::tests::fields(&flow.accepted),
        before
    );
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_report().is_none());
    assert!(molecular.last_report().is_none());
    let replay = flow
        .step_with_stabilized_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            controls,
            || false,
        )
        .unwrap();
    close(
        replay.coupled.elapsed.value(),
        report.coupled.elapsed.value(),
        2e-11,
    );
    // The independent existing C/A/C owner sees uniform fields, so D=0 exactly.
    let mut plain = reacting(&mechanism, 4, false);
    let expected = plain
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    close(
        replay.coupled.elapsed.value(),
        expected.elapsed.value(),
        2e-11,
    );
    for (a, b) in flow.branches().iter().zip(plain.branches()) {
        for i in 0..a.cell_count() {
            let (a, b) = (a.cell(i).unwrap(), b.cell(i).unwrap());
            for k in 0..2 {
                close(a.partial_densities[k], b.partial_densities[k], 3e-7);
            }
            close(a.energy, b.energy, 3e-6);
        }
    }
}

#[test]
fn g1_adaptive_reacting_network_resolves_continuous_decay_without_manual_diffusion_steps() {
    use crate::multicomponent::network::dissipation::AdaptiveNetworkDissipationControl;
    let rate = 2.0;
    let diffusion = 10.0;
    let (mut flow, mut molecular) = two_pool(rate, diffusion);
    let controls = AdaptiveNetworkDissipationControl {
        relative_tolerance: 1e-8,
        absolute_partial_density_tolerance_kg_m3: 1e-11,
        ..AdaptiveNetworkDissipationControl::default()
    };
    let mut temporal_rejections = 0;
    for _ in 0..4 {
        let report = flow
            .step_with_adaptive_dissipation(
                &mut molecular,
                Time::new(0.05),
                ReactiveStepControl::default(),
                controls,
                || false,
            )
            .unwrap();
        // Equal-volume A<->B kinetics commutes with this isobaric two-pool
        // diffusion mode, giving an independent continuous-time solution.
        let contrast = 0.4 * (-(2.0 * rate + 2.0 * diffusion) * flow.time().value()).exp();
        for cell in 0..2 {
            close(
                flow.branches()[0].cell(cell).unwrap().partial_densities[0],
                0.5 + contrast,
                4e-6,
            );
        }
        close(
            flow.plenums()[0].state().partial_densities[0],
            0.5 - contrast,
            4e-6,
        );
        let physical = report.coupled.dissipation.unwrap();
        assert_eq!(
            physical.first,
            report.dissipation.first.stabilized.transport
        );
        assert_eq!(
            physical.second,
            report.dissipation.second.stabilized.transport
        );
        assert_eq!(
            report.operator_evaluations,
            report.dissipation.first.stabilized.operator_evaluations
                + report.dissipation.second.stabilized.operator_evaluations
        );
        for half in [report.dissipation.first, report.dissipation.second] {
            assert!(half.maximum_estimated_local_error <= 1.0);
            temporal_rejections += half.temporal_rejections;
        }
        close(report.coupled.mass_residual_kg, 0.0, 3e-12);
        close(report.coupled.absolute_energy_residual_j, 0.0, 3e-11);
    }
    assert!(
        temporal_rejections > 0,
        "the initial stable macrostep needs temporal refinement"
    );
}

#[test]
fn g4_adaptive_heat_release_retries_cannot_reset_the_molecular_work_allowance() {
    use crate::multicomponent::network::dissipation::AdaptiveNetworkDissipationControl;
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let make = || {
        let raw = network(&mechanism, 4, false, 1.0);
        let molecular = operator(&raw, 0.1, 0.0, 0.0);
        let flow = ReactingNetwork1d::new(raw, mechanism.clone(), chemistry(), || false).unwrap();
        (flow, molecular)
    };
    let controls = AdaptiveNetworkDissipationControl::default();
    let (mut reference, mut reference_operator) = make();
    let report = reference
        .step_with_adaptive_dissipation(
            &mut reference_operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            controls,
            || false,
        )
        .unwrap();
    assert!(report.coupled.coupling_retries > 0);
    let accepted = report.dissipation.first.stabilized.operator_evaluations
        + report.dissipation.second.stabilized.operator_evaluations;
    assert!(report.operator_evaluations > accepted);
    let (mut flow, mut molecular) = make();
    let original = snapshot(&flow);
    let mut limited = controls;
    limited.stabilized.maximum_operator_evaluations = accepted;
    let error = flow
        .step_with_adaptive_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            limited,
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,
        ReactiveNetworkError::Dissipation(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == accepted));
    assert_eq!(snapshot(&flow), original);
    assert!(molecular.last_report().is_none());
    let replay = flow
        .step_with_adaptive_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            controls,
            || false,
        )
        .unwrap();
    assert_eq!(replay.operator_evaluations, report.operator_evaluations);
    close(
        replay.coupled.elapsed.value(),
        report.coupled.elapsed.value(),
        3e-12,
    );
    for (a, b) in flow.branches().iter().zip(reference.branches()) {
        for i in 0..a.cell_count() {
            let (a, b) = (a.cell(i).unwrap(), b.cell(i).unwrap());
            for k in 0..2 {
                close(a.partial_densities[k], b.partial_densities[k], 3e-10);
            }
            close(a.momentum, b.momentum, 3e-10);
            close(a.energy, b.energy, 3e-9);
        }
    }
    for (a, b) in flow.plenums().iter().zip(reference.plenums()) {
        for k in 0..2 {
            close(
                a.state().partial_densities[k],
                b.state().partial_densities[k],
                3e-10,
            );
        }
        close(a.state().energy, b.state().energy, 3e-9);
    }
}

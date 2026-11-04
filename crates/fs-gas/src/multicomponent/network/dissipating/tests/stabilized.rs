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

// Independent continuum integration of the three physical capacities and
// exterior species transfer. No production stage, flux or report is used.
fn bath_continuum(duration: f64) -> [f64; 4] {
    let rhs = |q: [f64; 4]| {
        let incoming = 0.03 * (0.9 - q[0]) / 0.2;
        let interior = 0.03 * (q[0] - q[1]) / 0.4;
        let plenum = 0.03 * (q[1] - q[2]) / 0.35;
        [
            (incoming - interior) / 0.4,
            (interior - plenum) / 0.4,
            0.2 * plenum / 0.07,
            0.2 * incoming,
        ]
    };
    let mut q = [0.4, 0.4, 0.1, 0.0];
    let h = duration / 1024.0;
    for _ in 0..1024 {
        let a = rhs(q);
        let b = rhs(core::array::from_fn(|i| q[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| q[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| q[i] + h * c[i]));
        q = core::array::from_fn(|i| q[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    q
}

#[test]
fn g1_adaptive_composition_controls_bath_diffusion_error_and_physical_port_exchange() {
    use crate::multicomponent::network::dissipation::AdaptiveNetworkDissipationControl;
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1e-4, 1e-10] {
        let mut flow = bath_owner(&gas);
        let report = flow
            .step_adaptive_with_sources(
                Time::new(0.04),
                AdaptiveDissipatingStepControl {
                    molecular: AdaptiveNetworkDissipationControl {
                        relative_tolerance: tolerance,
                        absolute_partial_density_tolerance_kg_m3: 1e-13,
                        ..AdaptiveNetworkDissipationControl::default()
                    },
                    ..AdaptiveDissipatingStepControl::default()
                },
                Some(&[[Some(bath.as_cell()), None]]),
                NetworkSources::default(),
                || false,
            )
            .unwrap();
        let expected = bath_continuum(report.coupled.elapsed.value());
        let actual = [
            flow.branches()[0].cell(0).unwrap().partial_densities[0],
            flow.branches()[0].cell(1).unwrap().partial_densities[0],
            flow.plenums()[0].state().partial_densities[0],
            report.coupled.boundary_change.species_masses[0],
        ];
        let error = actual
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        work.push(report.operator_evaluations);
        for half in [report.dissipation.first, report.dissipation.second] {
            assert!(half.maximum_estimated_local_error <= 1.0);
            assert_eq!(
                half.stabilized.transport.elapsed.value(),
                0.5 * report.coupled.elapsed.value()
            );
        }
        assert_eq!(
            report.operator_evaluations,
            report.dissipation.first.stabilized.operator_evaluations
                + report.dissipation.second.stabilized.operator_evaluations
        );
        assert_eq!(
            report.coupled.molecular.first,
            report.dissipation.first.stabilized.transport
        );
        assert_eq!(
            report.coupled.molecular.second,
            report.dissipation.second.stabilized.transport
        );
        let port = flow.accepted_molecular_port_exchanges().unwrap()[0];
        let advective = flow.accepted_port_exchanges().unwrap()[0];
        close(
            0.07 * (actual[2] - 0.1),
            port.species_masses_to_plenum[0] + advective.species_masses_to_plenum[0],
            3e-12,
        );
        close(
            report.coupled.conservation_residual.species_masses[0],
            0.0,
            3e-12,
        );
        close(report.coupled.absolute_energy_residual_j, 0.0, 3e-11);
        assert_eq!(flow.last_report(), Some(report.coupled));
    }
    assert!(errors[1] < 0.1 * errors[0], "physical errors: {errors:?}");
    assert!(errors[1] < 2e-8, "tight physical error: {}", errors[1]);
    assert!(work[1] > work[0], "actual operator work: {work:?}");
}

#[test]
fn g4_adaptive_composition_charges_both_halves_and_preserves_previous_publication() {
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let reservoirs = [[Some(bath.as_cell()), None]];
    let controls = AdaptiveDissipatingStepControl::default();
    let advance = |flow: &mut DissipatingNetwork1d<_>| {
        flow.step_adaptive_with_sources(
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
        advance(&mut flow);
        flow
    };
    let mut reference = make();
    let mut polls = 0;
    let expected = reference
        .step_adaptive_with_sources(
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
    let mut candidate = make();
    let original = snapshot(&candidate);
    let original_operator = candidate.molecular.last_report();
    let original_ports = candidate
        .molecular
        .accepted_port_exchanges()
        .unwrap()
        .to_vec();
    let limit = expected.dissipation.first.stabilized.operator_evaluations;
    let mut limited = controls;
    limited.molecular.stabilized.maximum_operator_evaluations = limit;
    let error = candidate
        .step_adaptive_with_sources(
            Time::new(0.04),
            limited,
            Some(&reservoirs),
            NetworkSources::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,
        DissipatingNetworkError::Molecular(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == limit));
    assert_eq!(snapshot(&candidate), original);
    assert_eq!(candidate.molecular.last_report(), original_operator);
    assert_eq!(
        candidate.molecular.accepted_port_exchanges().unwrap(),
        original_ports
    );
    let mut visits = 0;
    assert_eq!(
        candidate.step_adaptive_with_sources(
            Time::new(0.04),
            controls,
            Some(&reservoirs),
            NetworkSources::default(),
            || {
                visits += 1;
                visits >= polls
            }
        ),
        Err(DissipatingNetworkError::Cancelled)
    );
    assert_eq!(snapshot(&candidate), original);
    assert_eq!(candidate.molecular.last_report(), original_operator);
    let replay = advance(&mut candidate);
    assert_eq!(replay.operator_evaluations, expected.operator_evaluations);
    close(
        replay.coupled.elapsed.value(),
        expected.coupled.elapsed.value(),
        3e-12,
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
            close(a.partial_densities[k], b.partial_densities[k], 3e-12);
        }
        close(a.momentum, b.momentum, 3e-12);
        close(a.energy, b.energy, 3e-11);
    }
}

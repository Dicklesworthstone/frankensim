//! Accurate source-conjugate diffusion, finite storage and atomic entropy admission.
use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

fn controls(tolerance: f64, step: f64) -> AdaptiveNetworkDissipationControl {
    AdaptiveNetworkDissipationControl {
        stabilized: StabilizedDissipationControl {
            transport: DissipationControl {
                maximum_step_s: step,
                ..Default::default()
            },
            ..Default::default()
        },
        relative_tolerance: 0.0,
        absolute_partial_density_tolerance_kg_m3: tolerance,
        absolute_temperature_tolerance_k: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
    }
}

#[test]
fn tolerance_refinement_retains_exact_fourier_storage_and_source_entropy() {
    let gas = gas(1, false, 0.0);
    let volume = 0.07;
    let duration = 0.6;
    let cb = 2.0 * AREA * WIDTH * 2.5;
    let cp = volume * 2.0 * 2.5;
    let conductance = 2.0 * AREA * 0.05 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let decay = (-duration * conductance * (1.0 / cb + 1.0 / cp)).exp();
    let mean = (cb * 2.0 + cp * 4.0) / (cb + cp);
    let expected = [mean + (2.0 - mean) * decay, mean + (4.0 - mean) * decay];
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1e-4, 1e-7] {
        let mut network = two_pools(
            &gas,
            gas.conserved(&[1.0], 0.0, 2.0).unwrap(),
            gas.conserved(&[2.0], 0.0, 8.0).unwrap(),
            volume,
        );
        let mut op = NetworkDissipationOperator::new(
            &network,
            law(&gas, 0.0, 0.0, 0.05, 0.0),
            vec![Length::new(PLENUM_DISTANCE); 2],
            || false,
        )
        .unwrap();
        let initial_entropy = source_entropy(&network, &[(1.0, 2.5, 0.0)]);
        let (report, entropy) = op
            .advance_adaptive_entropy_consistent_for(
                &mut network,
                Time::new(duration),
                controls(tolerance, duration),
                MolecularEntropyControl::default(),
                || false,
            )
            .unwrap();
        let values = [
            gas.primitive(network.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .temperature,
            gas.primitive(network.plenums()[0].state(), || false)
                .unwrap()
                .temperature,
        ];
        errors.push(
            values.into_iter().zip(expected)
                .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max),
        );
        work.push(report.stabilized.operator_evaluations);
        assert!(report.maximum_estimated_local_error <= 1.0);
        assert!(entropy.minimum_step_entropy_change_j_k > 0.0);
        close(
            entropy.entropy_change_j_k,
            cb * (values[0] / 2.0).ln() + cp * (values[1] / 4.0).ln(),
            4e-12,
        );
        close(
            entropy.source_entropy_change_j_k,
            source_entropy(&network, &[(1.0, 2.5, 0.0)]) - initial_entropy,
            4e-12,
        );
        close(report.stabilized.transport.absolute_energy_residual_j, 0.0, 4e-12);
        for port in op.accepted_port_exchanges().unwrap() {
            close(port.energy_to_plenum, 0.5 * cp * (values[1] - 4.0), 4e-12);
            assert_eq!(port.mass_to_plenum(), 0.0);
        }
        assert_eq!(network.time().value(), 0.0);
        assert_eq!(op.last_report(), Some(report.stabilized.transport));
    }
    assert!(errors[1] < 2e-6, "errors={errors:?}");
    assert!(errors[1] < 0.04 * errors[0], "errors={errors:?}");
    assert!(work[1] > work[0]);
}

#[test]
fn accurate_conjugate_mixing_repairs_engineering_entropy_and_retains_source_gauge() {
    let duration = 0.01;
    let controls = controls(1e-8, duration);
    let calorics = [(1.0, 2.5, 0.0), (1.0, 5.5, 0.0), (1.0, 3.0, 0.0)];
    let mut native = contrasting_pools(0.0, true, false);
    let initial_entropy = source_entropy(&native, &calorics);
    diffusion(&native)
        .advance_adaptive_for(&mut native, Time::new(duration), controls, || false)
        .unwrap();
    assert!(source_entropy(&native, &calorics) < initial_entropy - 1e-6);

    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut network = contrasting_pools(formation, true, false);
        let initial_entropy = source_entropy(&network, &calorics);
        let mut op = diffusion(&network);
        let (report, entropy) = op
            .advance_adaptive_entropy_consistent_for(
                &mut network,
                Time::new(duration),
                controls,
                MolecularEntropyControl::default(),
                || false,
            )
            .unwrap();
        assert!(entropy.entropy_change_j_k > 0.0);
        assert!(entropy.minimum_step_entropy_change_j_k > 0.0);
        close(
            entropy.source_entropy_change_j_k,
            source_entropy(&network, &calorics) - initial_entropy,
            4e-12,
        );
        assert!(report.maximum_estimated_local_error <= 1.0);
        assert_eq!(report.stabilized.transport.after.species_masses[2], 0.0);
        assert_eq!(report.stabilized.transport.boundary_change.energy, 0.0);
        assert!(entropy.maximum_relative_species_residual < 1e-13);
        if let Some((state, measured)) = &reference {
            assert_eq!(&fields(&network), state);
            assert_eq!(&entropy, measured);
        } else {
            close(report.stabilized.transport.absolute_energy_residual_j, 0.0, 4e-12);
            reference = Some((fields(&network), entropy));
        }
    }
}

fn variable_mixture() -> (MixtureNetwork1d, NetworkDissipationOperator<ConstantTransport>) {
    let gas = gas(2, true, 0.0);
    let network = one_port(&gas, false);
    let op = NetworkDissipationOperator::new(
        &network,
        law(&gas, 0.02, 0.01, 0.04, 0.015),
        vec![Length::new(PLENUM_DISTANCE)],
        || false,
    )
    .unwrap();
    (network, op)
}

#[test]
fn unequal_masses_variable_calorics_and_port_impulse_keep_complete_balances() {
    let (mut network, mut op) = variable_mixture();
    let before = network.totals();
    let calorics = [(1.0, 2.5, 0.1), (0.5, 3.3, 0.05)];
    let initial_entropy = source_entropy(&network, &calorics);
    let (report, entropy) = op
        .advance_adaptive_entropy_consistent_for(
            &mut network,
            Time::new(0.2),
            controls(1e-8, 0.2),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.temporal_rejections > 0);
    assert!(report.stabilized.transport.accepted_steps > 1);
    assert!(report.maximum_estimated_local_error <= 1.0);
    assert!(entropy.minimum_step_entropy_change_j_k > 0.0);
    close(
        entropy.source_entropy_change_j_k,
        source_entropy(&network, &calorics) - initial_entropy,
        4e-12,
    );
    for k in 0..2 {
        close(report.stabilized.transport.after.species_masses[k], before.species_masses[k], 4e-14);
    }
    close(report.stabilized.transport.after.energy, before.energy, 4e-12);
    let port = op.accepted_port_exchanges().unwrap()[0];
    assert_ne!(port.momentum_to_plenum, 0.0);
    assert_eq!(network.plenums()[0].state().momentum, 0.0);
    for k in 0..2 {
        close(
            port.species_masses_to_plenum[k],
            0.07 * (network.plenums()[0].state().partial_densities[k] - [0.2, 0.8][k]),
            4e-14,
        );
    }
    assert_eq!(report.stabilized.transport.boundary_change.energy, 0.0);
    assert_eq!(network.time().value(), 0.0);
}

#[test]
fn entropy_path_preserves_accepted_state_and_receipts_on_late_refusal_and_replay() {
    let controls = controls(1e-8, 0.2);
    let make = || {
        let (mut network, mut op) = variable_mixture();
        op.advance_adaptive_entropy_consistent_for(
            &mut network, Time::new(0.001), controls,
            MolecularEntropyControl::default(), || false,
        ).unwrap();
        (network, op)
    };
    let (mut reference, mut reference_op) = make();
    let mut polls = 0;
    let expected = reference_op.advance_adaptive_entropy_consistent_for(
        &mut reference, Time::new(0.2), controls, MolecularEntropyControl::default(), || {
            polls += 1;
            false
        },
    ).unwrap();
    for work_limit in [false, true] {
        let (mut network, mut op) = make();
        let before = fields(&network);
        let prior = op.last_report();
        let ports = op.accepted_port_exchanges().unwrap().to_vec();
        let branches = op.accepted_branch_reports().unwrap().to_vec();
        let mut limited = controls;
        if work_limit {
            limited.stabilized.maximum_operator_evaluations = expected.0.stabilized.operator_evaluations - 1;
        }
        let mut visits = 0;
        let error = op.advance_adaptive_entropy_consistent_for(
            &mut network, Time::new(0.2), limited, MolecularEntropyControl::default(), || {
                visits += 1;
                !work_limit && visits >= polls
            },
        ).unwrap_err();
        if work_limit {
            assert!(matches!(error, DissipationError::OperatorWorkLimit { .. }));
        } else {
            assert_eq!(error, DissipationError::Cancelled);
        }
        assert_eq!(fields(&network), before);
        assert_eq!(op.last_report(), prior);
        assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
        assert_eq!(op.accepted_branch_reports().unwrap(), branches);
        assert_eq!(network.time().value(), 0.0);
        let replay = op.advance_adaptive_entropy_consistent_for(
            &mut network, Time::new(0.2), controls, MolecularEntropyControl::default(), || false,
        ).unwrap();
        assert_eq!(replay.0.stabilized.operator_evaluations, expected.0.stabilized.operator_evaluations);
        let actual = fields(&network);
        let wanted = fields(&reference);
        for (a, b) in actual.0.iter().flatten().chain(&actual.1)
            .zip(wanted.0.iter().flatten().chain(&wanted.1))
        {
            for (&a, &b) in a.partial_densities.iter().zip(&b.partial_densities) {
                close(a, b, 3e-12);
            }
            close(a.momentum, b.momentum, 3e-12);
            close(a.energy, b.energy, 3e-12);
        }
        close(replay.1.entropy_change_j_k, expected.1.entropy_change_j_k, 4e-12);
    }
}

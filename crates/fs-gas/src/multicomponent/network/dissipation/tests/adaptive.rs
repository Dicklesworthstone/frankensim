//! Accuracy against independent finite-storage modes and boundary quadrature.
use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

fn control(density: f64, temperature: f64, step: f64) -> AdaptiveNetworkDissipationControl {
    AdaptiveNetworkDissipationControl {
        stabilized: StabilizedDissipationControl {
            transport: DissipationControl {
                maximum_step_s: step,
                ..DissipationControl::default()
            },
            ..StabilizedDissipationControl::default()
        },
        relative_tolerance: 0.0,
        absolute_partial_density_tolerance_kg_m3: density,
        absolute_temperature_tolerance_k: temperature,
        absolute_velocity_tolerance_m_s: 1e-6,
    }
}

// Closed RKL2 stability polynomial, independent of the conservative stages.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut previous, mut current) = (1.0, x);
    for j in 2..=stages {
        let jf = j as f64;
        let next = ((2.0 * jf - 1.0) * x * current - (jf - 1.0) * previous) / jf;
        previous = current;
        current = next;
    }
    1.0 - b + b * current
}

fn binary(
    volume: f64,
) -> (
    MixtureNetwork1d,
    NetworkDissipationOperator<ConstantTransport>,
) {
    let gas = gas(2, false, 0.0);
    let network = two_pools(
        &gas,
        gas.conserved(&[0.55, 0.45], 0.0, 2.0).unwrap(),
        gas.conserved(&[0.45, 0.55], 0.0, 2.0).unwrap(),
        volume,
    );
    let op = NetworkDissipationOperator::new(
        &network,
        law(&gas, 0.0, 0.0, 0.0, 0.03),
        vec![Length::new(PLENUM_DISTANCE); 2],
        || false,
    )
    .unwrap();
    (network, op)
}

fn binary_mode(volume: f64, decay: f64) -> [f64; 2] {
    let branch_volume = 2.0 * AREA * WIDTH;
    let mean = (branch_volume * 0.55 + volume * 0.45) / (branch_volume + volume);
    [mean + (0.55 - mean) * decay, mean + (0.45 - mean) * decay]
}

fn binary_rate(volume: f64) -> f64 {
    let conductance = 2.0 * AREA * 0.03 / (0.5 * WIDTH + PLENUM_DISTANCE);
    conductance * (1.0 / (2.0 * AREA * WIDTH) + 1.0 / volume)
}

#[test]
fn accepted_probe_uses_same_stage_count_and_conservative_half_step_polynomial() {
    let volume = 0.0002;
    let duration = 0.2;
    let (mut network, mut op) = binary(volume);
    // Loose physical accuracy admits one long, linearly stable probe. Its
    // independently known two-half polynomial distinguishes stage-count changes
    // and Richardson extrapolation from the published physical trajectory.
    let report = op
        .advance_adaptive_for(
            &mut network,
            Time::new(duration),
            control(1.0, 1.0, duration),
            || false,
        )
        .unwrap();
    let s = report.stabilized.maximum_stages_used;
    assert!(s > 8);
    assert_eq!(report.stabilized.transport.accepted_steps, 1);
    assert_eq!(report.stabilized.transport.rejected_steps, 0);
    assert_eq!(report.stabilized.operator_evaluations, 3 * s + 1);
    assert!(report.maximum_estimated_local_error <= 1.0);
    let decay = amplification(s, -0.5 * duration * binary_rate(volume)).powi(2);
    let expected = binary_mode(volume, decay);
    close(
        network.branches()[0].cell(0).unwrap().partial_densities[0],
        expected[0],
        3e-12,
    );
    close(
        network.plenums()[0].state().partial_densities[0],
        expected[1],
        3e-12,
    );
    for port in op.accepted_port_exchanges().unwrap() {
        close(
            port.species_masses_to_plenum[0],
            0.5 * volume * (expected[1] - 0.45),
            2e-14,
        );
        close(port.mass_to_plenum(), 0.0, 3e-14);
    }
    for residual in &report
        .stabilized
        .transport
        .conservation_residual
        .species_masses[..2]
    {
        close(*residual, 0.0, 3e-14);
    }
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(op.last_report(), Some(report.stabilized.transport));
}

#[test]
fn requested_accuracy_resolves_stiff_finite_inventory_exponential() {
    let volume = 0.001;
    let duration = 0.1;
    let expected = binary_mode(volume, (-duration * binary_rate(volume)).exp());
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for accuracy in [1e-4, 1e-7] {
        let (mut network, mut op) = binary(volume);
        let report = op
            .advance_adaptive_for(
                &mut network,
                Time::new(duration),
                control(accuracy, 1.0, duration),
                || false,
            )
            .unwrap();
        let values = [
            network.branches()[0].cell(0).unwrap().partial_densities[0],
            network.plenums()[0].state().partial_densities[0],
        ];
        errors.push(
            values
                .into_iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
        work.push(report.stabilized.operator_evaluations);
        assert!(report.temporal_rejections > 0);
        assert!(report.maximum_estimated_local_error <= 1.0);
        assert!(report.stabilized.maximum_scaled_parabolic_number <= 0.40000000000001);
    }
    assert!(errors[1] < 3e-6, "errors={errors:?}");
    assert!(errors[1] < 0.04 * errors[0], "errors={errors:?}");
    assert!(work[1] > work[0]);
    // The original stability-only macrostep remains available, but its error
    // does not satisfy this requested accuracy merely because it is stable.
    let (mut network, mut op) = binary(volume);
    op.advance_stabilized_for(
        &mut network,
        Time::new(duration),
        control(1.0, 1.0, duration).stabilized,
        || false,
    )
    .unwrap();
    let uncontrolled = (network.plenums()[0].state().partial_densities[0] - expected[1]).abs();
    assert!(uncontrolled > 100.0 * errors[1]);
}

#[test]
fn temperature_accuracy_uses_actual_unequal_heat_capacities() {
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
        let report = op
            .advance_adaptive_for(
                &mut network,
                Time::new(duration),
                control(1.0, tolerance, duration),
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
            values
                .into_iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
        assert!(report.maximum_estimated_local_error <= 1.0);
        close(
            op.accepted_port_exchanges()
                .unwrap()
                .iter()
                .map(|p| p.energy_to_plenum)
                .sum(),
            cp * (values[1] - 4.0),
            4e-12,
        );
        close(
            report.stabilized.transport.absolute_energy_residual_j,
            0.0,
            4e-12,
        );
    }
    assert!(errors[1] < 2e-6, "errors={errors:?}");
    assert!(errors[1] < 0.04 * errors[0], "errors={errors:?}");
}

// Taylor evaluation of exp(t A) for the independently specified three-storage
// affine diffusion system, shifted by its uniform .9 bath equilibrium.
fn bath_reference(duration: f64) -> [f64; 3] {
    let mut term = [-0.5, -0.5, -0.8];
    let mut sum = term;
    for order in 1..=64 {
        let exterior = -0.03 * term[0] / 0.2;
        let interior = 0.03 * (term[0] - term[1]) / 0.4;
        let port = 0.03 * (term[1] - term[2]) / 0.35;
        let derivative = [
            (exterior - interior) / 0.4,
            (interior - port) / 0.4,
            0.2 * port / 0.07,
        ];
        term = derivative.map(|x| duration * x / order as f64);
        for k in 0..3 {
            sum[k] += term[k];
        }
    }
    sum.map(|x| x + 0.9)
}

#[test]
fn bath_quadrature_and_finite_ports_refine_without_formation_or_orientation_dependence() {
    let duration = 2.0;
    let expected = bath_reference(duration);
    let expected_external =
        AREA * WIDTH * (expected[0] + expected[1] - 0.8) + 0.07 * (expected[2] - 0.1);
    for reversed in [false, true] {
        let mut original = None;
        for formation in [0.0, 1e18] {
            let gas = gas(2, false, formation);
            let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
            let reservoirs = if reversed {
                [[None, Some(bath.as_cell())]]
            } else {
                [[Some(bath.as_cell()), None]]
            };
            let mut network = reservoir::bath_network(&gas, reversed);
            let mut op = NetworkDissipationOperator::new(
                &network,
                law(&gas, 0.0, 0.0, 0.0, 0.03),
                vec![Length::new(PLENUM_DISTANCE)],
                || false,
            )
            .unwrap();
            let report = op
                .advance_adaptive_with_reservoirs_for(
                    &mut network,
                    Time::new(duration),
                    control(1e-8, 1e-6, duration),
                    &reservoirs,
                    || false,
                )
                .unwrap();
            for (i, value) in expected[..2].iter().enumerate() {
                close(
                    network.branches()[0]
                        .cell(if reversed { 1 - i } else { i })
                        .unwrap()
                        .partial_densities[0],
                    *value,
                    4e-6,
                );
            }
            close(
                network.plenums()[0].state().partial_densities[0],
                expected[2],
                4e-6,
            );
            close(
                report.stabilized.transport.boundary_change.species_masses[0],
                expected_external,
                8e-7,
            );
            let port = op.accepted_port_exchanges().unwrap()[0];
            close(
                port.species_masses_to_plenum[0],
                0.07 * (network.plenums()[0].state().partial_densities[0] - 0.1),
                4e-14,
            );
            let branch = op.accepted_branch_reports().unwrap()[0];
            let inner = if reversed {
                branch.left_boundary_change
            } else {
                branch.right_boundary_change
            };
            assert_eq!(inner.species_masses[0], -port.species_masses_to_plenum[0]);
            close(
                report
                    .stabilized
                    .transport
                    .conservation_residual
                    .species_masses[0],
                0.0,
                5e-14,
            );
            close(
                report.stabilized.transport.boundary_change.mass(),
                0.0,
                5e-14,
            );
            assert!(report.maximum_estimated_local_error <= 1.0);
            if let Some((state, boundary)) = &original {
                assert_eq!(&fields(&network), state);
                assert_eq!(&report.stabilized.transport.boundary_change, boundary);
            } else {
                original = Some((
                    fields(&network),
                    report.stabilized.transport.boundary_change,
                ));
            }
        }
    }
}

#[test]
fn all_probes_share_work_and_late_failures_preserve_fields_reports_and_replay() {
    let controls = control(1e-7, 1e-6, 0.1);
    let make = || {
        let (mut network, mut op) = binary(0.02);
        op.advance_adaptive_for(&mut network, Time::new(0.001), controls, || false)
            .unwrap();
        network
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        (network, op)
    };
    let (mut probe, mut probe_op) = make();
    let mut polls = 0;
    let expected = probe_op
        .advance_adaptive_for(&mut probe, Time::new(0.3), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(expected.temporal_rejections > 0);
    assert!(expected.stabilized.transport.accepted_steps > 2);
    for budget in [false, true] {
        let (mut network, mut op) = make();
        let before = fields(&network);
        let time = network.time();
        let flow_report = network.last_report();
        let flow_ports = network.accepted_port_exchanges().unwrap().to_vec();
        let report = op.last_report();
        let ports = op.accepted_port_exchanges().unwrap().to_vec();
        let branches = op.accepted_branch_reports().unwrap().to_vec();
        let mut limited = controls;
        if budget {
            limited.stabilized.maximum_operator_evaluations =
                expected.stabilized.operator_evaluations - 1;
        }
        let mut visits = 0;
        let result = op.advance_adaptive_for(&mut network, Time::new(0.3), limited, || {
            visits += 1;
            !budget && visits >= polls
        });
        if budget {
            assert!(
                matches!(result, Err(DissipationError::OperatorWorkLimit { operator_evaluations, .. })
                if operator_evaluations == limited.stabilized.maximum_operator_evaluations)
            );
        } else {
            assert_eq!(result, Err(DissipationError::Cancelled));
        }
        assert_eq!(fields(&network), before);
        assert_eq!(network.time(), time);
        assert_eq!(network.last_report(), flow_report);
        assert_eq!(network.accepted_port_exchanges().unwrap(), flow_ports);
        assert_eq!(op.last_report(), report);
        assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
        assert_eq!(op.accepted_branch_reports().unwrap(), branches);
        let replay = op
            .advance_adaptive_for(&mut network, Time::new(0.3), controls, || false)
            .unwrap();
        assert_eq!(
            replay.stabilized.operator_evaluations,
            expected.stabilized.operator_evaluations
        );
        let actual = fields(&network);
        let wanted = fields(&probe);
        for (a, b) in actual
            .0
            .iter()
            .flatten()
            .chain(&actual.1)
            .zip(wanted.0.iter().flatten().chain(&wanted.1))
        {
            for (a, b) in a.partial_densities.iter().zip(&b.partial_densities) {
                close(*a, *b, 3e-12);
            }
            close(a.momentum, b.momentum, 3e-12);
            close(a.energy, b.energy, 3e-12);
        }
    }
}

#[test]
fn accuracy_admission_and_zero_interval_have_no_hidden_work() {
    let (mut network, mut op) = binary(0.02);
    let before = fields(&network);
    let mut invalid = control(1e-7, 1e-6, 0.1);
    invalid.relative_tolerance = f64::NAN;
    assert!(matches!(
        op.advance_adaptive_for(&mut network, Time::new(0.1), invalid, || false),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(fields(&network), before);
    let report = op
        .advance_adaptive_for(
            &mut network,
            Time::new(0.0),
            control(1e-7, 1e-6, 0.1),
            || false,
        )
        .unwrap();
    assert_eq!(report.stabilized.operator_evaluations, 0);
    assert_eq!(report.maximum_estimated_local_error, 0.0);
    assert_eq!(report.temporal_rejections, 0);
    assert_eq!(fields(&network), before);
}

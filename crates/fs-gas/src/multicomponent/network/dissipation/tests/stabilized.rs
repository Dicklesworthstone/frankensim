//! Independent finite-storage eigenmodes and original network publication.
use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    }
}

// Closed stability polynomial 1-b_s+b_s P_s(1+w_1 z), evaluated with the
// ordinary Legendre polynomial. This does not reproduce the stage recurrence.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut old, mut previous) = (1.0, x);
    for j in 2..=stages {
        let jf = j as f64;
        let next = ((2.0 * jf - 1.0) * x * previous - (jf - 1.0) * old) / jf;
        old = previous;
        previous = next;
    }
    1.0 - b + b * previous
}

#[test]
fn small_finite_plenum_uses_legendre_eigenmode_and_fewer_spatial_calls() {
    let g = gas(2, false, 0.0);
    let volume = 0.0002;
    let duration = 0.2;
    let make = || {
        let network = two_pools(
            &g,
            g.conserved(&[0.55, 0.45], 0.0, 2.0).unwrap(),
            g.conserved(&[0.45, 0.55], 0.0, 2.0).unwrap(),
            volume,
        );
        let op = NetworkDissipationOperator::new(
            &network,
            law(&g, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(PLENUM_DISTANCE); 2],
            || false,
        )
        .unwrap();
        (network, op)
    };
    let (mut network, mut op) = make();
    let report = op
        .advance_stabilized_for(
            &mut network,
            Time::new(duration),
            controls(duration),
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 8);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.transport.maximum_parabolic_number > 10.0);
    let branch_volume = 2.0 * AREA * WIDTH;
    let conductance = 2.0 * AREA * 0.03 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let rate = conductance * (1.0 / branch_volume + 1.0 / volume);
    let mean = (branch_volume * 0.55 + volume * 0.45) / (branch_volume + volume);
    let decay = amplification(report.maximum_stages_used, -rate * duration);
    let expected_branch = mean + (0.55 - mean) * decay;
    let expected_plenum = mean + (0.45 - mean) * decay;
    for i in 0..2 {
        let q = network.branches()[0].cell(i).unwrap();
        close(q.partial_densities[0], expected_branch, 3e-13);
        close(q.partial_densities.iter().sum(), 1.0, 5e-13);
        close(g.primitive(q, || false).unwrap().temperature, 2.0, 3e-11);
    }
    close(
        network.plenums()[0].state().partial_densities[0],
        expected_plenum,
        3e-12,
    );
    for port in op.accepted_port_exchanges().unwrap() {
        close(
            port.species_masses_to_plenum[0],
            0.5 * volume * (expected_plenum - 0.45),
            2e-14,
        );
        close(port.mass_to_plenum(), 0.0, 2e-14);
    }
    for value in report.transport.conservation_residual.species_masses[..2].iter() {
        close(*value, 0.0, 5e-14);
    }
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-12);
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(op.last_report(), Some(report.transport));
    let (mut explicit, mut explicit_op) = make();
    let explicit_report = explicit_op
        .advance_for(
            &mut explicit,
            Time::new(duration),
            controls(duration).transport,
            || false,
        )
        .unwrap();
    // SSP evaluates two spatial operators per accepted step in this constant
    // linear fixture. This compares stability work, not equal temporal error.
    assert!(4 * report.operator_evaluations < 2 * explicit_report.accepted_steps);
}

fn heat_error(step: f64) -> f64 {
    let g = gas(1, false, 0.0);
    let volume = 0.07;
    let mut network = two_pools(
        &g,
        g.conserved(&[1.0], 0.0, 2.0).unwrap(),
        g.conserved(&[2.0], 0.0, 8.0).unwrap(),
        volume,
    );
    let mut op = NetworkDissipationOperator::new(
        &network,
        law(&g, 0.0, 0.0, 0.05, 0.0),
        vec![Length::new(PLENUM_DISTANCE); 2],
        || false,
    )
    .unwrap();
    let duration = 0.6;
    let report = op
        .advance_stabilized_for(&mut network, Time::new(duration), controls(step), || false)
        .unwrap();
    let cb = (2.0 * AREA * WIDTH) * 2.5;
    let cp = volume * 2.0 * 2.5;
    let conductance = 2.0 * AREA * 0.05 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let rate = conductance * (1.0 / cb + 1.0 / cp);
    let mean = (cb * 2.0 + cp * 4.0) / (cb + cp);
    let decay = (-rate * duration).exp();
    let expected = [mean + (2.0 - mean) * decay, mean + (4.0 - mean) * decay];
    let mut error = 0.0_f64;
    for i in 0..2 {
        let q = network.branches()[0].cell(i).unwrap();
        assert_eq!(q.partial_densities, [1.0]);
        error = error.max((g.primitive(q, || false).unwrap().temperature - expected[0]).abs());
    }
    let q = network.plenums()[0].state();
    assert_eq!(q.partial_densities, [2.0]);
    let final_temperature = g.primitive(q, || false).unwrap().temperature;
    error = error.max((final_temperature - expected[1]).abs());
    let ports = op.accepted_port_exchanges().unwrap();
    close(
        ports.iter().map(|p| p.energy_to_plenum).sum(),
        cp * (final_temperature - 4.0),
        3e-12,
    );
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-12);
    error
}

#[test]
fn finite_heat_capacities_refine_to_independent_fourier_exponential() {
    let coarse = heat_error(0.06);
    let fine = heat_error(0.03);
    assert!(fine < 5e-6, "fine={fine:e}");
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

// Independent affine three-storage ODE with a prescribed .9 bath. Shift by
// equilibrium and apply the Legendre matrix polynomial, including finite Vp.
fn bath_polynomial(stages: usize, duration: f64) -> [f64; 3] {
    let anchor = [-0.5, -0.5, -0.8];
    let s = stages as f64;
    let weight = 4.0 * duration / (s * s + s - 2.0);
    let apply = |q: [f64; 3]| {
        let exterior = -0.03 * q[0] / 0.2;
        let interior = 0.03 * (q[0] - q[1]) / 0.4;
        let port = 0.03 * (q[1] - q[2]) / 0.35;
        [
            q[0] + weight * (exterior - interior) / 0.4,
            q[1] + weight * (interior - port) / 0.4,
            q[2] + weight * 0.2 * port / 0.07,
        ]
    };
    let mut old = anchor;
    let mut previous = apply(anchor);
    for j in 2..=stages {
        let jf = j as f64;
        let applied = apply(previous);
        let next =
            core::array::from_fn(|k| ((2.0 * jf - 1.0) * applied[k] - (jf - 1.0) * old[k]) / jf);
        old = previous;
        previous = next;
    }
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    core::array::from_fn(|k| 0.9 + (1.0 - b) * anchor[k] + b * previous[k])
}

#[test]
fn exterior_bath_and_finite_inventory_follow_one_polynomial_and_signed_ledger() {
    let duration = 2.0;
    for reversed in [false, true] {
        let mut baseline = None;
        for reference in [0.0, 1e18] {
            let g = gas(2, false, reference);
            let bath = g.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
            let reservoirs = if reversed {
                [[None, Some(bath.as_cell())]]
            } else {
                [[Some(bath.as_cell()), None]]
            };
            let mut network = reservoir::bath_network(&g, reversed);
            let mut op = NetworkDissipationOperator::new(
                &network,
                law(&g, 0.0, 0.0, 0.0, 0.03),
                vec![Length::new(PLENUM_DISTANCE)],
                || false,
            )
            .unwrap();
            let report = op
                .advance_stabilized_with_reservoirs_for(
                    &mut network,
                    Time::new(duration),
                    controls(duration),
                    &reservoirs,
                    || false,
                )
                .unwrap();
            assert_eq!(report.transport.accepted_steps, 1);
            assert_eq!(report.transport.rejected_steps, 0);
            assert!(report.maximum_stages_used > 2);
            let expected = bath_polynomial(report.maximum_stages_used, duration);
            let exterior =
                AREA * WIDTH * (expected[0] + expected[1] - 0.8) + 0.07 * (expected[2] - 0.1);
            for (i, value) in expected[..2].iter().enumerate() {
                let q = network.branches()[0]
                    .cell(if reversed { 1 - i } else { i })
                    .unwrap();
                close(q.partial_densities[0], *value, 5e-13);
                close(q.partial_densities[1], 1.0 - value, 5e-13);
                close(g.primitive(q, || false).unwrap().temperature, 2.0, 3e-11);
            }
            close(
                network.plenums()[0].state().partial_densities[0],
                expected[2],
                5e-13,
            );
            close(
                report.transport.boundary_change.species_masses[0],
                exterior,
                5e-14,
            );
            close(
                report.transport.boundary_change.species_masses[1],
                -exterior,
                5e-14,
            );
            close(report.transport.boundary_change.mass(), 0.0, 5e-14);
            close(report.transport.conservation_residual.energy, 0.0, 3e-12);
            let port = op.accepted_port_exchanges().unwrap()[0];
            close(
                port.species_masses_to_plenum[0],
                0.07 * (expected[2] - 0.1),
                5e-14,
            );
            let branch = op.accepted_branch_reports().unwrap()[0];
            let (outer, inner) = if reversed {
                (branch.right_boundary_change, branch.left_boundary_change)
            } else {
                (branch.left_boundary_change, branch.right_boundary_change)
            };
            close(outer.species_masses[0], exterior, 5e-14);
            assert_eq!(inner.species_masses[0], -port.species_masses_to_plenum[0]);
            if let Some(previous) = &baseline {
                assert_eq!(&fields(&network), previous);
            } else {
                baseline = Some(fields(&network));
            }
        }
    }
}

#[test]
fn bounded_work_and_late_cancellation_preserve_accepted_network_and_replay() {
    let g = gas(2, true, 5.0);
    let make = || {
        let mut network = one_port(&g, false);
        let mut op = NetworkDissipationOperator::new(
            &network,
            law(&g, 0.04, 0.02, 0.05, 0.03),
            vec![Length::new(0.1)],
            || false,
        )
        .unwrap();
        op.advance_stabilized_for(&mut network, Time::new(0.01), controls(0.005), || false)
            .unwrap();
        network
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        (network, op)
    };
    let (mut probe, mut probe_op) = make();
    let mut polls = 0;
    let expected = probe_op
        .advance_stabilized_for(&mut probe, Time::new(0.04), controls(0.005), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(expected.operator_evaluations >= 24);
    let (mut network, mut op) = make();
    let before = fields(&network);
    let time = network.time();
    let wave_report = network.last_report();
    let wave_ports = network.accepted_port_exchanges().unwrap().to_vec();
    let wave_branches = network.accepted_branch_reports().unwrap().to_vec();
    let report = op.last_report();
    let ports = op.accepted_port_exchanges().unwrap().to_vec();
    let branches = op.accepted_branch_reports().unwrap().to_vec();
    for budget_refusal in [false, true] {
        let mut visits = 0;
        let mut control = controls(0.005);
        if budget_refusal {
            // Exhaust after one completed but unpublished macrostep, proving
            // the work budget is interval-wide rather than reset per step.
            control.maximum_operator_evaluations = 4;
        }
        let result = op.advance_stabilized_for(&mut network, Time::new(0.04), control, || {
            visits += 1;
            !budget_refusal && visits >= polls
        });
        if budget_refusal {
            assert!(matches!(
                result,
                Err(DissipationError::OperatorWorkLimit {
                    operator_evaluations: 4,
                    ..
                })
            ));
        } else {
            assert_eq!(result, Err(DissipationError::Cancelled));
        }
        assert_eq!(fields(&network), before);
        assert_eq!(network.time(), time);
        assert_eq!(network.last_report(), wave_report);
        assert_eq!(network.accepted_port_exchanges().unwrap(), wave_ports);
        assert_eq!(network.accepted_branch_reports().unwrap(), wave_branches);
        assert_eq!(op.last_report(), report);
        assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
        assert_eq!(op.accepted_branch_reports().unwrap(), branches);
    }
    let replay = op
        .advance_stabilized_for(&mut network, Time::new(0.04), controls(0.005), || false)
        .unwrap();
    assert_eq!(
        replay.transport.accepted_steps,
        expected.transport.accepted_steps
    );
    let replay_fields = fields(&network);
    let expected_fields = fields(&probe);
    for (a, b) in replay_fields
        .0
        .iter()
        .flatten()
        .chain(&replay_fields.1)
        .zip(expected_fields.0.iter().flatten().chain(&expected_fields.1))
    {
        for k in 0..2 {
            close(a.partial_densities[k], b.partial_densities[k], 3e-12);
        }
        close(a.momentum, b.momentum, 3e-12);
        close(a.energy, b.energy, 3e-11);
    }
}

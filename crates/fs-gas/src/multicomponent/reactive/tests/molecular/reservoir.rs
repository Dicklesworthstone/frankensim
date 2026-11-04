//! Physical reservoir exchange through the complete original reacting owner.
use super::*;
use crate::multicomponent::Cell;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

fn open_field(
    mechanism: &IdealGasKineticsV1,
    count: usize,
    velocity: f64,
    initial: impl Fn(f64) -> [f64; 2],
) -> MixtureEuler1d {
    let mut flow = field(mechanism, count, velocity, initial);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    };
    flow
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}",
    );
}

fn balance(report: &ReactiveStepReport) {
    for k in 0..report.before.species_count + 2 {
        close(report.conservation_residual.get(k), 0.0, 5e-12);
    }
    close(report.mass_residual_kg, 0.0, 5e-12);
    close(report.absolute_energy_residual_j, 0.0, 5e-11);
}

fn step(
    stabilized: bool,
    flow: &mut ReactingEuler1d,
    operator: &mut DissipationOperator<ConstantTransport>,
    duration: f64,
    baths: [Option<Cell<'_>>; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveStepReport, ReactiveError> {
    let transport = DissipationControl {
        maximum_step_s: 0.001,
        ..Default::default()
    };
    if stabilized {
        flow.step_with_stabilized_dissipation_and_reservoirs(
            operator,
            Time::new(duration),
            ReactiveStepControl::default(),
            StabilizedDissipationControl {
                transport,
                ..Default::default()
            },
            baths,
            cancelled,
        )
        .map(|report| report.coupled)
    } else {
        flow.step_with_dissipation_and_reservoirs(
            operator,
            Time::new(duration),
            ReactiveStepControl::default(),
            transport,
            baths,
            cancelled,
        )
    }
}

#[test]
fn g1_isobaric_reacting_baths_match_independent_dirichlet_diffusion_mode() {
    let count = 16;
    let diffusion = 0.08;
    let reaction = 1.7;
    let end = 0.07;
    let pi = core::f64::consts::PI;
    // A <-> B has equal forward/reverse rates. Its linear relaxation and the
    // cell-centered Dirichlet diffusion operator share this exact eigenmode.
    let decay = 2.0 * reaction
        + 4.0 * diffusion * (count as f64 * (pi / (2.0 * count as f64)).sin()).powi(2);
    for stabilized in [false, true] {
        let mechanism = mechanism(0.0, reaction, 0.0);
        let initial = open_field(&mechanism, count, 0.0, |x| {
            let a = 0.5 + 0.15 * (pi * x).sin();
            [a, 1.0 - a]
        });
        let bath = initial.gas().conserved(&[0.5, 0.5], 0.0, 1.0).unwrap();
        let mut operator = molecular(&initial, diffusion, 0.0, 0.0);
        let chemical_controls = IsochoricKineticsSolve {
            relative_tolerance: 2e-9,
            absolute_amount_tolerance_mol: 1e-16,
            ..chemistry()
        };
        let mut flow =
            ReactingEuler1d::new(initial, mechanism, chemical_controls, || false).unwrap();
        let mut chemical_a = 0.0;
        let mut molecular_a = 0.0;
        let mut steps = 0;
        while flow.flow().time().value() < end {
            let remaining = end - flow.flow().time().value();
            let report = step(
                stabilized,
                &mut flow,
                &mut operator,
                remaining,
                [Some(bath.as_cell()); 2],
                || false,
            )
            .unwrap();
            balance(&report);
            chemical_a += report.chemistry.species_change_kg[0];
            let halves = report.dissipation.unwrap();
            for half in [halves.first, halves.second] {
                close(half.elapsed.value(), 0.5 * report.elapsed.value(), 1e-15);
                molecular_a += half.boundary_change.species_masses[0];
                assert!(half.left_boundary_change.species_masses[0] < 0.0);
                assert!(half.right_boundary_change.species_masses[0] < 0.0);
                close(half.boundary_change.mass(), 0.0, 2e-13);
            }
            close(report.transport.boundary_change.mass(), 0.0, 2e-12);
            steps += 1;
            assert!(steps < 100);
        }
        assert!(chemical_a < 0.0 && molecular_a < 0.0);
        assert_eq!(flow.flow().time().value(), end);
        let amplitude = 0.15 * (-decay * end).exp();
        for i in 0..count {
            let x = (i as f64 + 0.5) / count as f64;
            let cell = flow.flow().cell(i).unwrap();
            close(
                cell.partial_densities[0],
                0.5 + amplitude * (pi * x).sin(),
                5e-6,
            );
            close(cell.density(), 1.0, 2e-11);
            close(
                flow.flow()
                    .gas()
                    .primitive(cell, || false)
                    .unwrap()
                    .temperature,
                1.0,
                2e-10,
            );
        }
    }
}

#[test]
fn g1_inert_reacting_inlet_preserves_characteristic_species_and_enthalpy_transport() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut reference = open_field(&mechanism, 8, 0.0, |_| [0.3, 0.7]);
    let bath = reference.gas().conserved(&[0.9, 0.1], 0.0, 1.5).unwrap();
    let initial = reference.duplicate(&mut || false).unwrap();
    let mut flow = ReactingEuler1d::new(initial, mechanism, chemistry(), || false).unwrap();
    let duration = Time::new(0.0004);
    let controls = ReactiveStepControl::default();
    let expected = reference
        .step_with_reservoirs(
            duration,
            controls.transport,
            [Some(bath.as_cell()), None],
            || false,
        )
        .unwrap();
    let report = flow
        .step_with_reservoirs(duration, controls, [Some(bath.as_cell()), None], || false)
        .unwrap();
    assert_eq!(report.elapsed, expected.elapsed);
    assert_eq!(flow.flow().time(), reference.time());
    assert!(report.transport.left_boundary_change.mass() > 0.0);
    assert!(report.transport.left_boundary_change.energy > 0.0);
    assert_eq!(report.dissipation, None);
    close(
        report.transport.left_boundary_change.species_masses[0],
        0.9 * report.transport.left_boundary_change.mass(),
        2e-13,
    );
    balance(&report);
    for (&actual, &expected) in flow.accepted.states.iter().zip(&reference.states) {
        close(actual, expected, 3e-12);
    }
}

#[test]
fn g4_complete_bath_window_rolls_back_late_failure_and_preserves_global_work_budget() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let make = || {
        let initial = open_field(&mechanism, 16, 0.0, |_| [0.5, 0.5]);
        let operator = molecular(&initial, 0.2, 0.0, 0.0);
        (
            ReactingEuler1d::new(initial, mechanism.clone(), chemistry(), || false).unwrap(),
            operator,
        )
    };
    let (mut reference, mut reference_operator) = make();
    let left = reference
        .flow()
        .gas()
        .conserved(&[0.8, 0.2], 0.0, 1.0)
        .unwrap();
    let right = reference
        .flow()
        .gas()
        .conserved(&[0.2, 0.8], 0.0, 1.0)
        .unwrap();
    let baths = [Some(left.as_cell()), Some(right.as_cell())];
    let controls = StabilizedReactiveIntervalControl::default();
    let end = Time::new(0.06);
    let mut polls = 0;
    let expected = reference
        .advance_stabilized_with_reservoirs_to(
            &mut reference_operator,
            end,
            controls,
            baths,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.coupled_steps > 1 && expected.operator_evaluations > 1);
    close(expected.mass_residual_kg, 0.0, 5e-12);
    close(expected.absolute_energy_residual_j, 0.0, 5e-11);
    assert!(expected.left_boundary_change.species_masses[0] > 0.0);
    assert!(expected.right_boundary_change.species_masses[0] < 0.0);
    let (mut flow, mut operator) = make();
    let initial = flow.accepted.states.clone();
    let mut seen = 0;
    assert_eq!(
        flow.advance_stabilized_with_reservoirs_to(&mut operator, end, controls, baths, || {
            seen += 1;
            seen == polls
        },)
            .unwrap_err(),
        ReactiveError::Cancelled
    );
    assert_eq!(flow.accepted.states, initial);
    assert_eq!(flow.flow().time().value(), 0.0);

    for limited in [
        StabilizedReactiveIntervalControl {
            maximum_coupled_attempts: 1,
            ..controls
        },
        StabilizedReactiveIntervalControl {
            dissipation: StabilizedDissipationControl {
                maximum_operator_evaluations: expected.operator_evaluations - 1,
                ..controls.dissipation
            },
            ..controls
        },
    ] {
        let error = flow
            .advance_stabilized_with_reservoirs_to(&mut operator, end, limited, baths, || false)
            .unwrap_err();
        assert!(matches!(
            error,
            ReactiveError::IntervalWorkLimit { .. }
                | ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { .. })
        ));
        assert_eq!(flow.accepted.states, initial);
        assert_eq!(flow.flow().time().value(), 0.0);
    }
    let invalid = Cell {
        momentum: 0.1,
        ..left.as_cell()
    };
    assert!(
        flow.advance_stabilized_with_reservoirs_to(
            &mut operator,
            Time::new(0.0),
            controls,
            [Some(invalid), None],
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.accepted.states, initial);
    let repeated = flow
        .advance_stabilized_with_reservoirs_for(&mut operator, end, controls, baths, || false)
        .unwrap();
    assert_eq!(flow.accepted.states, reference.accepted.states);
    assert_eq!(repeated, expected);
}

#[test]
fn g4_bath_heat_release_retries_replay_both_halves_and_do_not_reset_rkl2_work() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let make = || {
        let initial = open_field(&mechanism, 4, 0.7, |_| [1.0, 0.0]);
        let operator = molecular(&initial, 0.05, 0.01, 0.1);
        (
            ReactingEuler1d::new(initial, mechanism.clone(), chemistry(), || false).unwrap(),
            operator,
        )
    };
    let (mut flow, mut operator) = make();
    let bath = flow.flow().gas().conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let baths = [Some(bath.as_cell()), None];
    let controls = ReactiveStepControl::default();
    let transport = StabilizedDissipationControl::default();
    let report = flow
        .step_with_stabilized_dissipation_and_reservoirs(
            &mut operator,
            Time::new(0.2),
            controls,
            transport,
            baths,
            || false,
        )
        .unwrap();
    assert!(report.coupled.coupling_retries > 0);
    balance(&report.coupled);
    let accepted_calls = report.dissipation.first.operator_evaluations
        + report.dissipation.second.operator_evaluations;
    assert!(report.operator_evaluations > accepted_calls);
    for half in [report.dissipation.first, report.dissipation.second] {
        assert_eq!(
            half.transport.elapsed.value(),
            0.5 * report.coupled.elapsed.value()
        );
    }
    let (mut reference, mut reference_operator) = make();
    let retried = reference
        .step_with_stabilized_dissipation_and_reservoirs(
            &mut reference_operator,
            report.coupled.elapsed,
            controls,
            transport,
            baths,
            || false,
        )
        .unwrap();
    assert_eq!(retried.coupled.coupling_retries, 0);
    // Discarded trials may leave different caloric inversion guesses. The
    // original EOS accepts guesses within its energy tolerance, so compare
    // physical trajectories and accounts at that accuracy, not cached history.
    for (&actual, &expected) in flow.accepted.states.iter().zip(&reference.accepted.states) {
        close(actual, expected, 3e-10);
    }
    assert_eq!(report.coupled.elapsed, retried.coupled.elapsed);
    assert_eq!(
        report.coupled.transport.elapsed,
        retried.coupled.transport.elapsed
    );
    for k in 0..report.coupled.before.species_count + 2 {
        close(
            report.coupled.transport.boundary_change.get(k),
            retried.coupled.transport.boundary_change.get(k),
            3e-10,
        );
        for (actual, expected) in [
            (report.dissipation.first, retried.dissipation.first),
            (report.dissipation.second, retried.dissipation.second),
        ] {
            close(
                actual.transport.boundary_change.get(k),
                expected.transport.boundary_change.get(k),
                3e-10,
            );
        }
    }
    for k in 0..report.coupled.before.species_count {
        close(
            report.coupled.chemistry.species_change_kg[k],
            retried.coupled.chemistry.species_change_kg[k],
            3e-10,
        );
    }
    close(
        report.coupled.chemistry.referenced_energy_change_j,
        retried.coupled.chemistry.referenced_energy_change_j,
        3e-10,
    );
    let (mut limited, mut limited_operator) = make();
    let initial = limited.accepted.states.clone();
    let error = limited
        .step_with_stabilized_dissipation_and_reservoirs(
            &mut limited_operator,
            Time::new(0.2),
            controls,
            StabilizedDissipationControl {
                maximum_operator_evaluations: accepted_calls,
                ..transport
            },
            baths,
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { .. })
    ));
    assert_eq!(limited.accepted.states, initial);
    assert_eq!(limited.flow().time().value(), 0.0);
}

#[test]
fn g0_absent_baths_keep_original_reacting_states_reports_and_cancellation_sequence() {
    for stabilized in [false, true] {
        let mechanism = mechanism(0.0, 2.0, 0.0);
        let make = || {
            let initial = field(&mechanism, 8, 0.3, |_| [0.6, 0.4]);
            let operator = molecular(&initial, 0.08, 0.01, 0.02);
            (
                ReactingEuler1d::new(initial, mechanism.clone(), chemistry(), || false).unwrap(),
                operator,
            )
        };
        let (mut a, mut oa) = make();
        let (mut b, mut ob) = make();
        let (mut pa, mut pb) = (0, 0);
        let transport = DissipationControl {
            maximum_step_s: 0.001,
            ..Default::default()
        };
        let original = if stabilized {
            a.step_with_stabilized_dissipation(
                &mut oa,
                Time::new(0.005),
                ReactiveStepControl::default(),
                StabilizedDissipationControl {
                    transport,
                    ..Default::default()
                },
                || {
                    pa += 1;
                    false
                },
            )
            .unwrap()
            .coupled
        } else {
            a.step_with_dissipation(
                &mut oa,
                Time::new(0.005),
                ReactiveStepControl::default(),
                transport,
                || {
                    pa += 1;
                    false
                },
            )
            .unwrap()
        };
        let extended = step(stabilized, &mut b, &mut ob, 0.005, [None; 2], || {
            pb += 1;
            false
        })
        .unwrap();
        assert_eq!(a.accepted.states, b.accepted.states);
        assert_eq!(original, extended);
        assert_eq!(pa, pb);
        let initial = b.accepted.states.clone();
        let bath = b.flow().gas().conserved(&[0.5, 0.5], 0.0, 1.0).unwrap();
        assert!(
            step(
                stabilized,
                &mut b,
                &mut ob,
                0.005,
                [Some(bath.as_cell()), None],
                || false
            )
            .is_err()
        );
        assert_eq!(b.accepted.states, initial);
    }
}

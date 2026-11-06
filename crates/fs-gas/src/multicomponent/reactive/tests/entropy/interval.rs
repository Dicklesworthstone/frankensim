//! Complete closed reacting windows retain one original entropy allowance.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use crate::multicomponent::dissipation::entropy::reactive::ReactiveEntropyAccount;

fn interval_controls() -> StabilizedReactiveIntervalControl {
    StabilizedReactiveIntervalControl {
        step: controls(),
        dissipation: StabilizedDissipationControl {
            maximum_operator_evaluations: 4096,
            ..StabilizedDissipationControl::default()
        },
        maximum_coupled_attempts: 128,
    }
}

fn warmed_variable_caloric() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>) {
    let (mut simulation, mut operator) = variable_caloric();
    simulation
        .step_entropy_consistent_with_stabilized_dissipation(
            &mut operator,
            Time::new(0.001),
            controls(),
            interval_controls().dissipation,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    (simulation, operator)
}

fn uniform_isomers() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>) {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let flow = field(&mechanism, 8, 0.7, |_| [0.8, 0.2]);
    let operator = molecular(&flow, 0.02, 0.0, 0.0);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap(),
        operator,
    )
}

#[test]
fn g1_complete_reacting_window_retains_original_caloric_entropy_and_chemical_sources() {
    let (mut simulation, mut operator) = warmed_variable_caloric();
    let start = simulation.flow().time();
    let duration = Time::new(0.025);
    let target = Time::new(start.value() + duration.value());
    let calorics = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
    let initial_entropy = source_entropy(simulation.flow(), &calorics);
    let initial = simulation.flow().totals();
    let entropy_control = MixtureEntropyControl::default();
    let report = simulation
        .advance_stabilized_entropy_consistent_for(
            &mut operator,
            duration,
            interval_controls(),
            entropy_control,
            || false,
        )
        .unwrap();
    let coupled = report.coupled;
    assert!(coupled.coupled_steps > 1);
    assert_eq!(coupled.initial_time, start);
    assert_eq!(coupled.final_time, target);
    assert_eq!(simulation.flow().time(), target);
    assert!((coupled.elapsed.value() - duration.value()).abs() < 3e-15);
    assert_eq!(
        coupled.chemistry.cell_solves,
        2 * simulation.flow().cell_count() * coupled.coupled_steps
    );
    assert!(coupled.operator_evaluations >= 6 * coupled.coupled_steps);
    assert!(coupled.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(coupled.chemistry.species_change_kg[1] > 1e-4);
    for k in 0..2 {
        assert!(
            (coupled.after.species_masses[k]
                - initial.species_masses[k]
                - coupled.chemistry.species_change_kg[k])
                .abs()
                < 3e-11
        );
    }
    assert!(coupled.absolute_energy_residual_j.abs() < 3e-10);
    assert!(coupled.mass_residual_kg.abs() < 3e-11);
    assert!(coupled.conservation_residual.momentum.abs() < 3e-11);
    assert_eq!(coupled.advection_boundary_change.mass(), 0.0);
    assert_eq!(coupled.molecular_boundary_change.mass(), 0.0);
    assert_eq!(coupled.advection_boundary_change.energy, 0.0);
    assert_eq!(coupled.molecular_boundary_change.energy, 0.0);

    // Reconstruct S(final)-S(original) directly. B's nonzero a9 is physical
    // chemical entropy; a sum that removes all source reference terms is wrong.
    let measured = source_entropy(simulation.flow(), &calorics) - initial_entropy;
    let chemical_datum = 1.5 * (coupled.after.species_masses[1] - initial.species_masses[1]);
    assert!(chemical_datum > 1e-4 && measured > chemical_datum);
    assert!((report.entropy.source_entropy_change_j_k - measured).abs() < 5e-10);
    assert!((report.entropy.entropy_change_j_k - measured).abs() < 5e-10);
    assert!(
        (report.entropy.source_entropy_change_j_k
            - report.entropy.reference_entropy_residual_j_k
            - report.entropy.entropy_change_j_k)
            .abs()
            < 3e-12
    );
    assert!(
        report.entropy.maximum_relative_species_residual
            <= entropy_control.relative_species_tolerance
    );
    assert!(
        report.maximum_step_relative_species_residual <= entropy_control.relative_species_tolerance
    );
    assert!(report.minimum_step_entropy_change_j_k >= -entropy_control.absolute_tolerance_j_k);
}

#[test]
fn g4_interval_source_and_rhs_budgets_restore_a_completed_prefix_then_reuse() {
    let target = Time::new(0.07);
    let controls = interval_controls();
    let entropy = MixtureEntropyControl::default();
    let (mut probe, mut probe_operator) = uniform_isomers();
    let first = probe
        .step_entropy_consistent_with_stabilized_dissipation(
            &mut probe_operator,
            target,
            controls.step,
            controls.dissipation,
            entropy,
            || false,
        )
        .unwrap();
    assert_eq!(first.coupled.coupled.coupling_retries, 0);
    assert!(first.coupled.coupled.elapsed.value() < target.value());
    assert!(first.coupled.coupled.chemistry.species_change_kg[1] > 0.0);
    let first_calls = first.operator_evaluations;
    assert!(first_calls >= 6);

    let (mut simulation, mut operator) = uniform_isomers();
    let original = simulation.accepted.states.clone();
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            target,
            StabilizedReactiveIntervalControl {
                maximum_coupled_attempts: 1,
                ..controls
            },
            entropy,
            || false,
        )
        .unwrap_err();
    match error {
        ReactiveError::IntervalWorkLimit {
            elapsed_time_s,
            completed_steps,
            attempts,
        } => {
            assert_eq!((completed_steps, attempts), (1, 1));
            assert!((elapsed_time_s - first.coupled.coupled.elapsed.value()).abs() < 3e-14);
        }
        other => panic!("expected interval refusal after one accepted step, got {other:?}"),
    }
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);

    // The first complete coupled step consumes this measured budget exactly.
    // A later step cannot obtain a fresh molecular allowance.
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            target,
            StabilizedReactiveIntervalControl {
                dissipation: StabilizedDissipationControl {
                    maximum_operator_evaluations: first_calls,
                    ..controls.dissipation
                },
                ..controls
            },
            entropy,
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ReactiveError::Dissipation(DissipationError::OperatorWorkLimit {
            operator_evaluations, ..
        }) if operator_evaluations == first_calls
    ));
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let recovered = simulation
        .advance_stabilized_entropy_consistent_to(&mut operator, target, controls, entropy, || {
            false
        })
        .unwrap();
    assert_eq!(simulation.flow().time(), target);
    assert!(recovered.coupled.coupled_steps > 1);
    assert!(recovered.coupled.operator_evaluations > first_calls);
    assert!(recovered.entropy.entropy_change_j_k > 0.0);
    assert!(recovered.coupled.absolute_energy_residual_j.abs() < 3e-11);
}

#[test]
fn g4_late_interval_entropy_cancellation_restores_the_original_field_clock_and_source_anchor() {
    let (mut reference, mut reference_operator) = warmed_variable_caloric();
    let start = reference.flow().time();
    let target = Time::new(start.value() + 0.025);
    let mut polls = 0;
    let expected = reference
        .advance_stabilized_entropy_consistent_to(
            &mut reference_operator,
            target,
            interval_controls(),
            MixtureEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.coupled.coupled_steps > 1);
    let (mut simulation, mut operator) = warmed_variable_caloric();
    let original = simulation.accepted.states.clone();
    let original_mesh = simulation.accepted.mesh.clone();
    let mut visits = 0;
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            target,
            interval_controls(),
            MixtureEntropyControl::default(),
            || {
                visits += 1;
                visits >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.accepted.mesh, original_mesh);
    assert_eq!(simulation.flow().time(), start);
    let recovered = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            target,
            interval_controls(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(simulation.flow().time(), target);
    assert_eq!(
        recovered.coupled.coupled_steps,
        expected.coupled.coupled_steps
    );
    assert!(
        simulation
            .accepted
            .states
            .iter()
            .zip(&reference.accepted.states)
            .all(|(a, b)| (a - b).abs() < 3e-9)
    );
    assert!(
        (recovered.entropy.source_entropy_change_j_k - expected.entropy.source_entropy_change_j_k)
            .abs()
            < 3e-9
    );
    assert!(recovered.coupled.absolute_energy_residual_j.abs() < 3e-10);
}

#[test]
fn g0_entropy_window_noop_and_represented_clock_contracts_remain_native() {
    let (mut simulation, mut operator) = uniform_isomers();
    let controls = interval_controls();
    let entropy = MixtureEntropyControl::default();
    let original = simulation.accepted.states.clone();
    let noop = simulation
        .advance_stabilized_entropy_consistent_for(
            &mut operator,
            Time::new(0.0),
            controls,
            entropy,
            || false,
        )
        .unwrap();
    assert_eq!(noop.coupled.elapsed.value(), 0.0);
    assert_eq!(noop.coupled.coupled_steps, 0);
    assert_eq!(noop.coupled.operator_evaluations, 0);
    assert_eq!(noop.coupled.initial_time, noop.coupled.final_time);
    assert_eq!(noop.minimum_step_entropy_change_j_k, 0.0);
    assert_eq!(noop.maximum_step_relative_species_residual, 0.0);
    assert_eq!(noop.entropy.entropy_change_j_k, 0.0);
    assert_eq!(noop.entropy.source_entropy_change_j_k, 0.0);
    assert_eq!(simulation.accepted.states, original);

    // Even a zero-duration request must admit the original entropy scope.
    let mut invalid = controls;
    invalid.step.transport.reconstruction = Reconstruction::Muscl;
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            Time::new(0.0),
            invalid,
            entropy,
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ReactiveError::Flow(FlowError::EntropyDomain(_))
    ));
    assert_eq!(simulation.accepted.states, original);

    simulation.accepted.time = 1e16;
    let error = simulation
        .advance_stabilized_entropy_consistent_for(
            &mut operator,
            Time::new(3.0),
            controls,
            entropy,
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("reacting duration is not representable at this clock")
    );
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            Time::new(1e16 - 2.0),
            controls,
            entropy,
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("finite nondecreasing reacting target required")
    );
    let error = simulation
        .advance_stabilized_entropy_consistent_to(
            &mut operator,
            Time::new(1e16 + 2.0),
            controls,
            entropy,
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("reacting time increment cannot advance the clock")
    );
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 1e16);
}

#[test]
fn g4_two_locally_admitted_entropy_losses_do_not_multiply_the_whole_window_allowance() {
    // This directly exercises the production state account on source-admitted
    // manufactured endpoints. It is not a claim that a flux evolved these states.
    // M=R, rho=1, Cv=2.5, total volume=0.1 gives dS=0.25*ln(Tfinal/Tinitial).
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let original = field(&mechanism, 2, 0.0, |_| [1.0, 0.0]);
    let endpoint = |temperature: f64| {
        let gas = original.gas().clone();
        let cell = gas.conserved(&[1.0, 0.0], 0.0, temperature).unwrap();
        MixtureEuler1d::new(
            gas,
            Length::new(0.5),
            0.1,
            vec![cell.clone(), cell],
            Boundaries::Periodic,
            || false,
        )
        .unwrap()
    };
    let epsilon = 1e-5;
    let middle = endpoint((-0.6_f64 * epsilon / 0.25).exp());
    let final_state = endpoint((-1.2_f64 * epsilon / 0.25).exp());
    let controls = MixtureEntropyControl {
        absolute_tolerance_j_k: epsilon,
        ..MixtureEntropyControl::default()
    };
    let mut first = ReactiveEntropyAccount::new(&original, controls, &mut || false).unwrap();
    first.begin_attempt();
    let a = first
        .check_candidate(&original, &middle, None, &mut || false)
        .unwrap();
    let mut second = ReactiveEntropyAccount::new(&middle, controls, &mut || false).unwrap();
    second.begin_attempt();
    let b = second
        .check_candidate(&middle, &final_state, None, &mut || false)
        .unwrap();
    assert!((a.entropy_change_j_k + 0.6 * epsilon).abs() < 1e-11);
    assert!((b.entropy_change_j_k + 0.6 * epsilon).abs() < 1e-11);
    assert!(a.entropy_change_j_k >= -epsilon && b.entropy_change_j_k >= -epsilon);

    let mut window = ReactiveEntropyAccount::new(&original, controls, &mut || false).unwrap();
    window.begin_attempt();
    window
        .merge_inventory(&first.accepted_inventory(), &mut || false)
        .unwrap();
    window
        .merge_inventory(&second.accepted_inventory(), &mut || false)
        .unwrap();
    let error = window
        .check_candidate_rows(&original.states, &final_state, None, &mut || false)
        .unwrap_err();
    match error {
        FlowError::EntropyDecrease {
            change_j_k,
            tolerance_j_k,
        } => {
            assert_eq!(tolerance_j_k, epsilon);
            assert!((change_j_k + 1.2 * epsilon).abs() < 1e-11);
        }
        other => panic!("expected complete-window entropy refusal, got {other:?}"),
    }
}

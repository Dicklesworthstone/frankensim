//! Complete open reacting windows use one original source/boundary entropy basis.

use super::*;
use super::exterior_molecular::{CALORICS, WITHOUT_DATUM, fixture, supply};
use crate::multicomponent::Conserved;
use crate::multicomponent::dissipation::{DissipationError, StabilizedDissipationControl};

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

fn warmed() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>, [Conserved; 2]) {
    let (mut simulation, mut operator, baths) = fixture();
    simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
        &mut operator, Time::new(0.001), controls(), interval_controls().dissipation,
        MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
    ).unwrap();
    (simulation, operator, baths)
}

fn uniform_isomers() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>, Conserved) {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let mut flow = field(&mechanism, 8, 0.0, |_| [0.8, 0.2]);
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    let bath = flow.gas().conserved(&[0.8, 0.2], 0.0, 1.0).unwrap();
    let operator = molecular(&flow, 0.02, 0.0, 0.0);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap(),
        operator, bath,
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 6e-10 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}",
    );
}

fn inputs(report: &StabilizedReactiveIntervalReport, calorics: &[[f64; 3]; 2]) -> [f64; 2] {
    // The native interval combines all accepted advective and both molecular
    // end transfers. Fixed stationary baths make their physical entropy effort
    // linear in these SAME transfers, with no momentum or advective psi term.
    [
        supply([0.8, 0.2], 1.3, calorics, &report.left_boundary_change),
        supply([0.3, 0.7], 0.9, calorics, &report.right_boundary_change),
    ]
}

#[test]
fn g1_warm_open_interval_integrates_exact_remainder_and_original_source_boundary_entropy() {
    let (mut simulation, mut operator, baths) = warmed();
    let start = simulation.flow().time();
    let duration = Time::new(0.025);
    let target = Time::new(start.value() + duration.value());
    let before = source_entropy(simulation.flow(), &CALORICS);
    let stripped_before = source_entropy(simulation.flow(), &WITHOUT_DATUM);
    let initial = simulation.flow().totals();
    let entropy = MixtureEntropyControl::default();
    let report = simulation.advance_stabilized_entropy_consistent_with_exterior_states_for(
        &mut operator, duration, interval_controls(), entropy,
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
    ).unwrap();
    let coupled = report.coupled;
    assert!(coupled.coupled_steps > 1);
    assert_eq!(coupled.initial_time, start);
    assert_eq!(coupled.final_time, target);
    assert_eq!(simulation.flow().time(), target);
    assert!((coupled.elapsed.value() - duration.value()).abs() < 3e-15);
    assert_eq!(coupled.chemistry.cell_solves, 2 * simulation.flow().cell_count() * coupled.coupled_steps);
    assert!(coupled.operator_evaluations >= 6 * coupled.coupled_steps);
    assert!(coupled.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(coupled.chemistry.species_change_kg[1] > 1e-4);
    for k in 0..2 {
        near(
            coupled.after.species_masses[k] - initial.species_masses[k],
            coupled.chemistry.species_change_kg[k]
                + coupled.advection_boundary_change.species_masses[k]
                + coupled.molecular_boundary_change.species_masses[k],
        );
        near(
            coupled.left_boundary_change.species_masses[k] + coupled.right_boundary_change.species_masses[k],
            coupled.advection_boundary_change.species_masses[k] + coupled.molecular_boundary_change.species_masses[k],
        );
    }
    near(coupled.absolute_energy_residual_j, 0.0);
    near(coupled.mass_residual_kg, 0.0);
    near(coupled.conservation_residual.momentum, 0.0);
    let measured = source_entropy(simulation.flow(), &CALORICS) - before;
    let boundary = inputs(&coupled, &CALORICS);
    for end in 0..2 { near(report.entropy.boundary_entropy_input_j_k[end], boundary[end]); }
    near(report.entropy.source_entropy_change_j_k, measured);
    near(report.entropy.source_entropy_generation_j_k, measured - boundary.iter().sum::<f64>());
    near(report.entropy.entropy_generation_j_k,
        report.entropy.source_entropy_generation_j_k - report.entropy.reference_entropy_residual_j_k);
    assert!(report.entropy.entropy_generation_j_k > 0.0);

    // The physical B datum follows chemical conversion after ALL boundary
    // transports have been removed, including the first/last molecular halves.
    let stripped = source_entropy(simulation.flow(), &WITHOUT_DATUM) - stripped_before
        - inputs(&coupled, &WITHOUT_DATUM).iter().sum::<f64>();
    let chemical_datum = 1.5 * coupled.chemistry.species_change_kg[1];
    assert!(chemical_datum > 1e-4);
    near(report.entropy.source_entropy_generation_j_k - stripped, chemical_datum);
    assert!(report.entropy.maximum_relative_species_residual <= entropy.relative_species_tolerance);
    assert!(report.maximum_step_relative_species_residual <= entropy.relative_species_tolerance);
    assert!(report.minimum_step_entropy_generation_j_k >= -entropy.absolute_tolerance_j_k);
}

#[test]
fn g4_open_interval_global_work_caps_restore_a_changed_accepted_prefix_and_allow_reuse() {
    let target = Time::new(0.07);
    let controls = interval_controls();
    let entropy = MixtureEntropyControl::default();
    let (mut probe, mut probe_operator, bath) = uniform_isomers();
    let probe_initial = probe.accepted.states.clone();
    let first = probe.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
        &mut probe_operator, target, controls.step, controls.dissipation, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap();
    assert_eq!(first.coupled.coupled.coupling_retries, 0);
    assert!(first.coupled.coupled.elapsed.value() < target.value());
    assert!(first.coupled.coupled.chemistry.species_change_kg[1] > 0.0);
    assert_ne!(probe.accepted.states, probe_initial);
    let first_calls = first.operator_evaluations;
    assert!(first_calls >= 6);

    let (mut simulation, mut operator, bath) = uniform_isomers();
    let initial = simulation.accepted.states.clone();
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, target,
        StabilizedReactiveIntervalControl { maximum_coupled_attempts: 1, ..controls },
        entropy, [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap_err();
    match error {
        ReactiveError::IntervalWorkLimit { elapsed_time_s, completed_steps, attempts } => {
            assert_eq!((completed_steps, attempts), (1, 1));
            near(elapsed_time_s, first.coupled.coupled.elapsed.value());
        }
        other => panic!("expected refusal after an accepted open prefix, got {other:?}"),
    }
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);

    // The measured first complete step spends the complete interval allowance.
    // A subsequent step must not create a fresh RHS budget.
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, target,
        StabilizedReactiveIntervalControl {
            dissipation: StabilizedDissipationControl {
                maximum_operator_evaluations: first_calls, ..controls.dissipation
            },
            ..controls
        },
        entropy, [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap_err();
    assert!(matches!(error, ReactiveError::Dissipation(
        DissipationError::OperatorWorkLimit { operator_evaluations, .. }
    ) if operator_evaluations == first_calls));
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let recovered = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, target, controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap();
    assert_eq!(simulation.flow().time(), target);
    assert!(recovered.coupled.coupled_steps > 1);
    assert!(recovered.coupled.operator_evaluations > first_calls);
    assert!(recovered.entropy.entropy_generation_j_k > 0.0);
    near(recovered.coupled.absolute_energy_residual_j, 0.0);
}

#[test]
fn g4_final_open_window_callback_cancellation_restores_original_state_clock_and_entropy_anchor() {
    let (mut reference, mut reference_operator, baths) = warmed();
    let start = reference.flow().time();
    let target = Time::new(start.value() + 0.025);
    let mut polls = 0;
    let expected = reference.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut reference_operator, target, interval_controls(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
        || { polls += 1; false },
    ).unwrap();
    assert!(expected.coupled.coupled_steps > 1);
    let (mut simulation, mut operator, baths) = warmed();
    let original = simulation.accepted.states.clone();
    let original_mesh = simulation.accepted.mesh.clone();
    let mut visits = 0;
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, target, interval_controls(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
        || { visits += 1; visits >= polls },
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.accepted.mesh, original_mesh);
    assert_eq!(simulation.flow().time(), start);
    let recovered = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, target, interval_controls(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
    ).unwrap();
    assert_eq!(simulation.flow().time(), target);
    assert_eq!(recovered.coupled.coupled_steps, expected.coupled.coupled_steps);
    assert!(simulation.accepted.states.iter().zip(&reference.accepted.states)
        .all(|(a, b)| (a - b).abs() < 3e-9));
    assert!((recovered.entropy.source_entropy_change_j_k
        - expected.entropy.source_entropy_change_j_k).abs() < 3e-9);
    for end in 0..2 {
        assert!((recovered.entropy.boundary_entropy_input_j_k[end]
            - expected.entropy.boundary_entropy_input_j_k[end]).abs() < 3e-9);
    }
    near(recovered.coupled.absolute_energy_residual_j, 0.0);
}

#[test]
fn g0_open_noop_still_admits_scope_stationary_baths_and_original_operator_source() {
    let (mut simulation, mut operator, bath) = uniform_isomers();
    let controls = interval_controls();
    let entropy = MixtureEntropyControl::default();
    let initial = simulation.accepted.states.clone();
    let noop = simulation.advance_stabilized_entropy_consistent_with_exterior_states_for(
        &mut operator, Time::new(0.0), controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap();
    assert_eq!(noop.coupled.elapsed.value(), 0.0);
    assert_eq!(noop.coupled.coupled_steps, 0);
    assert_eq!(noop.coupled.operator_evaluations, 0);
    assert_eq!(noop.coupled.initial_time, noop.coupled.final_time);
    assert_eq!(noop.minimum_step_entropy_generation_j_k, 0.0);
    assert_eq!(noop.maximum_step_relative_species_residual, 0.0);
    assert_eq!(noop.entropy.entropy_generation_j_k, 0.0);
    assert_eq!(noop.entropy.source_entropy_change_j_k, 0.0);
    assert_eq!(noop.entropy.boundary_entropy_input_j_k, [0.0; 2]);
    assert_eq!(simulation.accepted.states, initial);

    let mut invalid = controls;
    invalid.step.transport.reconstruction = Reconstruction::Muscl;
    assert!(simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, Time::new(0.0), invalid, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).is_err());
    let moving = simulation.flow().gas().conserved(&[0.8, 0.2], 0.1, 1.0).unwrap();
    assert!(simulation.advance_stabilized_entropy_consistent_with_exterior_states_for(
        &mut operator, Time::new(0.0), controls, entropy,
        [Some(moving.as_cell()), Some(bath.as_cell())], || false,
    ).is_err());

    // Same dimensions and physical mesh, different ORIGINAL source bank.
    // A zero window cannot bypass the original molecular owner/source check.
    let other = mechanism(0.0, 2.0, 1.0);
    let mut foreign_flow = field(&other, 8, 0.0, |_| [0.8, 0.2]);
    foreign_flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    let mut foreign_operator = molecular(&foreign_flow, 0.02, 0.0, 0.0);
    assert!(simulation.advance_stabilized_entropy_consistent_with_exterior_states_for(
        &mut foreign_operator, Time::new(0.0), controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).is_err());
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);
}

#[test]
fn g4_open_interval_preserves_represented_duration_target_and_stagnant_clock_refusals() {
    let (mut simulation, mut operator, bath) = uniform_isomers();
    let initial = simulation.accepted.states.clone();
    let controls = interval_controls();
    let entropy = MixtureEntropyControl::default();
    simulation.accepted.time = 1e16;
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_for(
        &mut operator, Time::new(3.0), controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Invalid("reacting duration is not representable at this clock"));
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, Time::new(1e16 - 2.0), controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Invalid("finite nondecreasing reacting target required"));
    let error = simulation.advance_stabilized_entropy_consistent_with_exterior_states_to(
        &mut operator, Time::new(1e16 + 2.0), controls, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], || false,
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Invalid("reacting time increment cannot advance the clock"));
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 1e16);
}

#[test]
fn g4_two_open_step_allowances_do_not_accumulate_into_the_complete_window() {
    // A direct production-account test on source-admitted manufactured
    // endpoints, NOT a claim that any flux evolved these cooling states.
    // The zero boundary inventory is produced by an actual zero-duration native
    // molecular call, never fabricated from private receipt fields.
    use crate::multicomponent::dissipation::entropy::reactive::ReactiveEntropyAccount;
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut original = field(&mechanism, 2, 0.0, |_| [1.0, 0.0]);
    let boundaries = Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive };
    original.set_boundaries(boundaries, || false).unwrap();
    let bath = original.gas().conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let epsilon = 1e-5;
    let entropy = MixtureEntropyControl {
        absolute_tolerance_j_k: epsilon, ..MixtureEntropyControl::default()
    };
    let mut zero_state = original.duplicate(&mut || false).unwrap();
    let mut operator = molecular(&zero_state, 0.02, 0.0, 0.0);
    let mut evaluations = 0;
    let (_, _, zero) = operator.advance_stabilized_entropy_consistent_with_reservoir_inventory_budget(
        &mut zero_state, Time::new(0.0), interval_controls().dissipation, entropy,
        [Some(bath.as_cell()), Some(bath.as_cell())], &mut evaluations, &mut || false,
    ).unwrap();
    assert_eq!(evaluations, 0);
    assert_eq!(zero_state.states, original.states);
    let endpoint = |temperature: f64| {
        let gas = original.gas().clone();
        let cell = gas.conserved(&[1.0, 0.0], 0.0, temperature).unwrap();
        MixtureEuler1d::new(
            gas, Length::new(0.5), 0.1, vec![cell.clone(), cell], boundaries, || false,
        ).unwrap()
    };
    // M=R, rho=1, Cv=2.5, total volume=.1:
    // Delta S=.25 ln(T_final/T_initial).
    let middle = endpoint((-0.6_f64 * epsilon / 0.25).exp());
    let final_state = endpoint((-1.2_f64 * epsilon / 0.25).exp());
    let mut window = ReactiveEntropyAccount::new_with_exterior_states(
        &original, entropy, [Some(bath.as_cell()), Some(bath.as_cell())], &mut || false,
    ).unwrap();
    let mut first = ReactiveEntropyAccount::new_step_in_basis(
        &original, &window, entropy, &mut || false,
    ).unwrap();
    first.begin_attempt();
    let a = first.check_candidate_with_exterior(
        &original, &middle, None, &zero, &mut || false,
    ).unwrap();
    window.merge_exterior_step(&first, &mut || false).unwrap();
    let mut second = ReactiveEntropyAccount::new_step_in_basis(
        &middle, &window, entropy, &mut || false,
    ).unwrap();
    second.begin_attempt();
    let b = second.check_candidate_with_exterior(
        &middle, &final_state, None, &zero, &mut || false,
    ).unwrap();
    window.merge_exterior_step(&second, &mut || false).unwrap();
    assert!((a.entropy_generation_j_k + 0.6 * epsilon).abs() < 1e-11);
    assert!((b.entropy_generation_j_k + 0.6 * epsilon).abs() < 1e-11);
    assert!(a.entropy_generation_j_k >= -epsilon && b.entropy_generation_j_k >= -epsilon);
    let error = window.check_candidate_exterior_rows(
        &original.states, &final_state, None, &mut || false,
    ).unwrap_err();
    match error {
        FlowError::EntropyDecrease { change_j_k, tolerance_j_k } => {
            assert_eq!(tolerance_j_k, epsilon);
            assert!((change_j_k + 1.2 * epsilon).abs() < 1e-11);
        }
        other => panic!("expected original-window entropy refusal, got {other:?}"),
    }
}

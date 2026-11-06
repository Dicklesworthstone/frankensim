//! Complete open chemistry/advection/chemistry entropy and source balances.
//! The original physical B entropy constant belongs to chemical conversion;
//! only the contribution of a numerical species continuity defect is removed.

use super::*;
use crate::multicomponent::Conserved;

const CALORICS: [[f64; 3]; 2] = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
const WITHOUT_DATUM: [[f64; 3]; 2] = [[3.5, 0.2, 0.0], [4.0, 0.1, 0.0]];
const T0: f64 = 0.001;
const AREA: f64 = 0.1;

fn exterior(gas: &MixtureGas, rho: [f64; 2], velocity: f64, temperature: f64) -> Conserved {
    // Both original fixture species have M=R and p0=1.
    gas.conserved(&rho, velocity, (rho[0] + rho[1]) * temperature).unwrap()
}

fn open_fixture() -> (ReactingEuler1d, [Conserved; 2]) {
    let (simulation, _) = variable_caloric();
    let mechanism = simulation.mechanism().clone();
    let chemistry = simulation.chemistry_controls();
    let mut flow = simulation.into_flow();
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    let outside = [
        exterior(flow.gas(), [0.8, 0.2], 0.25, 1.3),
        exterior(flow.gas(), [0.3, 0.7], -0.15, 0.9),
    ];
    (ReactingEuler1d::new(flow, mechanism, chemistry, || false).unwrap(), outside)
}

fn hot_fixture() -> (ReactingEuler1d, [Conserved; 2]) {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    // Both species are present in the exterior as well as after the chemistry
    // half step; no ideal-mixture activity floor is introduced.
    let mut flow = field(&mechanism, 4, 0.7, |_| [0.99, 0.01]);
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    let outside = exterior(flow.gas(), [0.99, 0.01], 0.7, 1.0);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap(),
        [outside.clone(), outside],
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-10 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}",
    );
}

// Differentiate the original polynomial ideal-mixture entropy directly.
// eta=-S has v_i=e_i,ref/T+1-s_i-u²/(2T), v_m=u/T, v_E=-1/T,
// psi=u*sum rho_i. Transfers already carry their inward orientation.
#[allow(clippy::too_many_arguments)]
fn boundary_input(
    rho: [f64; 2],
    velocity: f64,
    temperature: f64,
    calorics: &[[f64; 3]; 2],
    transfer: &Totals,
    sign: f64,
    duration: f64,
) -> f64 {
    let mut result = transfer.energy / temperature - velocity * transfer.momentum / temperature;
    for (k, &[a, b, c]) in calorics.iter().enumerate() {
        let internal = (a - 1.0) * (temperature - T0)
            + 0.5 * b * (temperature * temperature - T0 * T0);
        let entropy = (a - 1.0) * temperature.ln() + b * temperature + c - rho[k].ln();
        let effort = internal / temperature + 1.0 - entropy
            - 0.5 * velocity * velocity / temperature;
        result -= effort * transfer.species_masses[k];
    }
    result + sign * AREA * duration * velocity * (rho[0] + rho[1])
}

fn variable_inputs(step: &ReactiveStepReport, calorics: &[[f64; 3]; 2]) -> [f64; 2] {
    [
        boundary_input(
            [0.8, 0.2], 0.25, 1.3, calorics,
            &step.transport.left_boundary_change, 1.0, step.elapsed.value(),
        ),
        boundary_input(
            [0.3, 0.7], -0.15, 0.9, calorics,
            &step.transport.right_boundary_change, -1.0, step.elapsed.value(),
        ),
    ]
}

#[test]
fn g1_open_variable_caloric_reaction_keeps_true_species_entropy_and_actual_boundary_supply() {
    let (mut simulation, outside) = open_fixture();
    let initial_entropy = source_entropy(simulation.flow(), &CALORICS);
    let without_datum = source_entropy(simulation.flow(), &WITHOUT_DATUM);
    let initial = simulation.flow().totals();
    let original = simulation.accepted.states.clone();
    let report = simulation.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).unwrap();
    let step = report.coupled;
    assert_eq!(step.elapsed.value(), 0.002);
    assert_eq!(step.elapsed, step.transport.elapsed);
    assert_eq!(simulation.flow().time(), step.elapsed);
    assert_eq!(step.chemistry.cell_solves, 2 * simulation.flow().cell_count());
    assert!(step.dissipation.is_none());
    assert_eq!(step.transport.hll_fallback_faces, 0);
    assert_ne!(simulation.accepted.states, original);

    // Open exchange and chemical conversion each enter species continuity once.
    assert!(step.chemistry.species_change_kg[1] > 1e-4);
    for k in 0..2 {
        near(
            step.after.species_masses[k] - initial.species_masses[k],
            step.transport.boundary_change.species_masses[k] + step.chemistry.species_change_kg[k],
        );
    }
    assert!(report.entropy.maximum_relative_species_residual < 1e-12);
    near(step.absolute_energy_residual_j, 0.0);
    near(step.mass_residual_kg, 0.0);
    near(step.conservation_residual.momentum, 0.0);

    let measured = source_entropy(simulation.flow(), &CALORICS) - initial_entropy;
    let input = variable_inputs(&step, &CALORICS);
    for end in 0..2 {
        near(report.entropy.boundary_entropy_input_j_k[end], input[end]);
        near(report.transport_entropy.boundary_entropy_input_j_k[end], input[end]);
    }
    near(report.entropy.source_entropy_change_j_k, measured);
    near(report.entropy.source_entropy_generation_j_k, measured - input.iter().sum::<f64>());
    near(
        report.entropy.entropy_generation_j_k,
        report.entropy.source_entropy_generation_j_k - report.entropy.reference_entropy_residual_j_k,
    );
    assert!(report.entropy.entropy_generation_j_k > 0.0);
    assert!(report.transport_entropy.entropy_generation_j_k >= -1e-10);
    assert!((measured - step.chemistry.entropy_change_j_k
        - report.transport_entropy.source_entropy_change_j_k).abs() < 2e-8);

    // Remove B's 1.5*R datum only in the INDEPENDENT oracle. The difference
    // between these two state-minus-boundary values is actual reaction entropy,
    // not 1.5 times all B inventory change and not a numerical drift correction.
    let stripped = source_entropy(simulation.flow(), &WITHOUT_DATUM) - without_datum
        - variable_inputs(&step, &WITHOUT_DATUM).iter().sum::<f64>();
    let chemical_datum = 1.5 * step.chemistry.species_change_kg[1];
    assert!(chemical_datum > 1e-4);
    near(report.entropy.source_entropy_generation_j_k - stripped, chemical_datum);
}

#[test]
fn g1_hot_open_flow_restarts_both_chemical_halves_and_keeps_only_accepted_boundary_exchange() {
    let (mut simulation, outside) = hot_fixture();
    let original_entropy = source_entropy(simulation.flow(), &[[3.5, 0.0, 0.0]; 2]);
    let report = simulation.step_entropy_consistent_with_exterior_states(
        Time::new(0.2), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).unwrap();
    let step = report.coupled;
    assert!(step.coupling_retries > 0, "reaction heat must shorten the common acoustic interval");
    assert!(step.elapsed.value() < 0.2);
    assert_eq!(step.elapsed, step.transport.elapsed);
    assert_eq!(step.chemistry.cell_solves, 8);
    assert_eq!(simulation.flow().time(), step.elapsed);

    let inputs = [
        boundary_input(
            [0.99, 0.01], 0.7, 1.0, &[[3.5, 0.0, 0.0]; 2],
            &step.transport.left_boundary_change, 1.0, step.elapsed.value(),
        ),
        boundary_input(
            [0.99, 0.01], 0.7, 1.0, &[[3.5, 0.0, 0.0]; 2],
            &step.transport.right_boundary_change, -1.0, step.elapsed.value(),
        ),
    ];
    for end in 0..2 { near(report.entropy.boundary_entropy_input_j_k[end], inputs[end]); }
    let measured = source_entropy(simulation.flow(), &[[3.5, 0.0, 0.0]; 2]) - original_entropy;
    near(report.entropy.source_entropy_generation_j_k, measured - inputs.iter().sum::<f64>());
    for k in 0..2 {
        near(
            step.after.species_masses[k] - step.before.species_masses[k],
            step.transport.boundary_change.species_masses[k] + step.chemistry.species_change_kg[k],
        );
    }

    // Repeating the accepted duration from the SAME original state must agree.
    // Keeping a chemistry half from a longer rejected attempt cannot pass this.
    let (mut replay, outside) = hot_fixture();
    let replayed = replay.step_entropy_consistent_with_exterior_states(
        step.elapsed, controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).unwrap();
    assert_eq!(replayed.coupled.elapsed, step.elapsed);
    for (a, b) in replay.accepted.states.iter().zip(&simulation.accepted.states) { near(*a, *b); }
    near(replayed.entropy.entropy_generation_j_k, report.entropy.entropy_generation_j_k);
    for end in 0..2 {
        near(replayed.entropy.boundary_entropy_input_j_k[end], report.entropy.boundary_entropy_input_j_k[end]);
    }
}

#[test]
fn g0_absent_exterior_states_keep_original_closed_reacting_dispatch() {
    let (mut closed, _) = variable_caloric();
    let (mut through_open_api, _) = variable_caloric();
    let expected = closed.step_entropy_consistent(
        Time::new(0.002), controls(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    let actual = through_open_api.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(), [None; 2], || false,
    ).unwrap();
    assert_eq!(actual.coupled, expected.coupled);
    assert_eq!(through_open_api.accepted.states, closed.accepted.states);
    assert_eq!(through_open_api.flow().time(), closed.flow().time());
    assert_eq!(actual.entropy.boundary_entropy_input_j_k, [0.0; 2]);
    assert_eq!(actual.entropy.entropy_generation_j_k, expected.entropy.entropy_change_j_k);
    assert_eq!(actual.entropy.source_entropy_change_j_k, expected.entropy.source_entropy_change_j_k);
    assert_eq!(actual.entropy.reference_entropy_residual_j_k, expected.entropy.reference_entropy_residual_j_k);
    assert_eq!(actual.transport_entropy.boundary_entropy_input_j_k, [0.0; 2]);
    assert_eq!(actual.transport_entropy.entropy_generation_j_k, expected.transport_entropy.entropy_change_j_k);
}

#[test]
fn g4_complete_open_reacting_cancellation_is_atomic_and_the_original_owners_are_reusable() {
    let (mut reference, outside) = open_fixture();
    let mut polls = 0;
    let expected = reference.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())],
        || { polls += 1; false },
    ).unwrap();
    assert!(polls > 100);
    let (mut simulation, outside) = open_fixture();
    let initial = simulation.accepted.states.clone();
    let mut visits = 0;
    let error = simulation.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())],
        || { visits += 1; visits >= polls },
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let actual = simulation.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).unwrap();
    assert_eq!(actual.coupled.elapsed, expected.coupled.elapsed);
    assert_eq!(simulation.flow().time(), reference.flow().time());
    for (a, b) in simulation.accepted.states.iter().zip(&reference.accepted.states) { near(*a, *b); }
    near(actual.entropy.entropy_generation_j_k, expected.entropy.entropy_generation_j_k);
    for end in 0..2 {
        near(actual.entropy.boundary_entropy_input_j_k[end], expected.entropy.boundary_entropy_input_j_k[end]);
    }
}

#[test]
fn g4_open_reacting_scope_and_exterior_refusals_leave_the_accepted_field_unchanged() {
    let (mut simulation, outside) = open_fixture();
    let initial = simulation.accepted.states.clone();
    let wrong_axis = Conserved { partial_densities: vec![1.0], momentum: 0.0, energy: 2.0 };
    let one_sided_zero = exterior(simulation.flow().gas(), [1.0, 0.0], 0.0, 1.0);
    let mut bounded = controls();
    bounded.maximum_coupling_retries = 0;
    for ends in [
        [None, Some(outside[1].as_cell())],
        [Some(wrong_axis.as_cell()), Some(outside[1].as_cell())],
        [Some(one_sided_zero.as_cell()), Some(outside[1].as_cell())],
    ] {
        assert!(simulation.step_entropy_consistent_with_exterior_states(
            Time::new(0.002), bounded, MixtureEntropyControl::default(), ends, || false,
        ).is_err());
        assert_eq!(simulation.accepted.states, initial);
        assert_eq!(simulation.flow().time().value(), 0.0);
    }
    let mut muscl = controls();
    muscl.transport.reconstruction = Reconstruction::Muscl;
    assert!(simulation.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), muscl, MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).is_err());
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);
}

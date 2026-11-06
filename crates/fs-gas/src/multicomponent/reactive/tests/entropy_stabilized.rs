//! Original source-conjugate RKL2 inside the complete reacting entropy owner.

use super::*;
use super::molecular::molecular;
use crate::multicomponent::{MixtureEntropyControl, Reconstruction};
use crate::multicomponent::dissipation::{
    DissipationError, StabilizedDissipationControl,
};

fn controls() -> ReactiveStepControl {
    ReactiveStepControl {
        transport: StepControl {
            reconstruction: Reconstruction::FirstOrder,
            ..StepControl::default()
        },
        ..ReactiveStepControl::default()
    }
}

// Direct polynomial source entropy for this equal-mass, constant-Cv fixture.
// Source constants cancel with its individually admitted closed species/source
// balances; no entropy report or divided flux mean enters the oracle.
fn state_entropy(flow: &MixtureEuler1d) -> f64 {
    (0..flow.cell_count()).map(|i| {
        let cell = flow.cell(i).unwrap();
        let t = flow.gas().primitive(cell, || false).unwrap().temperature;
        let density: f64 = cell.partial_densities.iter().filter(|&&rho| rho > 0.0)
            .map(|&rho| rho * (2.5 * t.ln() - rho.ln())).sum();
        flow.mesh().volumes()[i] * density
    }).sum()
}

#[test]
fn g1_conjugate_rkl_coupled_mode_retains_actual_entropy_and_reduces_diffusive_rhs_work() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let initial = |x: f64| {
        let perturbation = 0.02 * (2.0 * core::f64::consts::PI * x).sin();
        [0.5 + perturbation, 0.5 - perturbation]
    };
    let original = field(&mechanism, 64, 0.1, initial);
    let before = state_entropy(&original);
    let mut rkl_operator = molecular(&original, 0.5, 0.02, 0.4);
    let mut rkl = ReactingEuler1d::new(
        original, mechanism.clone(), chemistry(), || false,
    ).unwrap();
    let explicit_original = field(&mechanism, 64, 0.1, initial);
    let mut explicit_operator = molecular(&explicit_original, 0.5, 0.02, 0.4);
    let mut explicit = ReactingEuler1d::new(
        explicit_original, mechanism, chemistry(), || false,
    ).unwrap();
    let rkl_report = rkl.step_entropy_consistent_with_stabilized_dissipation(
        &mut rkl_operator, Time::new(0.002), controls(),
        StabilizedDissipationControl::default(),
        MixtureEntropyControl::default(), || false,
    ).unwrap();
    let explicit_report = explicit.step_entropy_consistent_with_dissipation(
        &mut explicit_operator, Time::new(0.002), controls(),
        DissipationControl::default(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    let coupled = &rkl_report.coupled.coupled;
    assert_eq!(coupled.elapsed.value(), 0.002);
    assert_eq!(coupled.elapsed, explicit_report.coupled.elapsed);
    assert_eq!(rkl.flow().time(), coupled.elapsed);
    assert_eq!(coupled.coupling_retries, 0);
    assert_eq!(rkl_report.dissipation.first.transport.elapsed.value(), 0.001);
    assert_eq!(rkl_report.dissipation.second.transport.elapsed.value(), 0.001);
    assert!(rkl_report.dissipation.first.maximum_stages_used > 2);
    let explicit_halves = explicit_report.coupled.dissipation.unwrap();
    let explicit_accepted_rhs = 2 * (
        explicit_halves.first.accepted_steps + explicit_halves.second.accepted_steps
    );
    assert!(
        rkl_report.operator_evaluations < explicit_accepted_rhs,
        "RKL={}, explicit accepted RHS={explicit_accepted_rhs}",
        rkl_report.operator_evaluations,
    );
    assert_eq!(
        rkl_report.operator_evaluations,
        rkl_report.dissipation.first.operator_evaluations
            + rkl_report.dissipation.second.operator_evaluations,
    );
    let halves = rkl_report.coupled.dissipation_entropy.unwrap();
    assert!(halves.first.entropy_change_j_k >= -1e-10);
    assert!(halves.second.entropy_change_j_k >= -1e-10);
    let measured = state_entropy(rkl.flow()) - before;
    assert!(measured > 0.0);
    assert!((rkl_report.coupled.entropy.source_entropy_change_j_k - measured).abs() < 2e-10);
    assert!(rkl_report.coupled.entropy.entropy_change_j_k > 0.0);
    assert!(rkl_report.coupled.entropy.maximum_relative_species_residual < 1e-12);
    assert!(coupled.absolute_energy_residual_j.abs() < 2e-11);
    assert!(coupled.mass_residual_kg.abs() < 2e-12);
    let maximum_difference = rkl.accepted.states.iter().zip(&explicit.accepted.states)
        .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
    assert!(maximum_difference < 4e-7, "RK split difference={maximum_difference:e}");
}

#[test]
fn g4_conjugate_rkl_budget_counts_discarded_heat_release_attempts_and_reuse_is_atomic() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut simulation = uniform(&mechanism);
    let original = simulation.accepted.states.clone();
    let mut operator = molecular(simulation.flow(), 0.05, 0.01, 0.1);
    let error = simulation.step_entropy_consistent_with_stabilized_dissipation(
        &mut operator, Time::new(0.2), controls(),
        StabilizedDissipationControl {
            maximum_stages: 2,
            maximum_operator_evaluations: 3,
            ..StabilizedDissipationControl::default()
        },
        MixtureEntropyControl::default(), || false,
    ).unwrap_err();
    // The first uniform molecular half consumes exactly three native RHS calls.
    // Exothermic heat release shortens A and discards that entire attempt.
    // Its work remains spent: the next half cannot acquire a fresh counter.
    assert!(matches!(error, ReactiveError::Dissipation(
        DissipationError::OperatorWorkLimit { operator_evaluations: 3, .. }
    )));
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);

    let recovered = simulation.step_entropy_consistent_with_stabilized_dissipation(
        &mut operator, Time::new(0.2), controls(),
        StabilizedDissipationControl::default(),
        MixtureEntropyControl::default(), || false,
    ).unwrap();
    let coupled = &recovered.coupled.coupled;
    assert!(coupled.coupling_retries > 0);
    assert!(recovered.operator_evaluations
        > recovered.dissipation.first.operator_evaluations
            + recovered.dissipation.second.operator_evaluations);
    assert_eq!(simulation.flow().time(), coupled.elapsed);
    assert!(recovered.coupled.entropy.entropy_change_j_k > 0.0);
    assert!(coupled.absolute_energy_residual_j.abs() < 2e-11);

    let mut reference = uniform(&mechanism);
    let mut reference_operator = molecular(reference.flow(), 0.05, 0.01, 0.1);
    let expected = reference.step_entropy_consistent_with_stabilized_dissipation(
        &mut reference_operator, Time::new(0.2), controls(),
        StabilizedDissipationControl::default(),
        MixtureEntropyControl::default(), || false,
    ).unwrap();
    assert_eq!(coupled.elapsed, expected.coupled.coupled.elapsed);
    assert_eq!(recovered.operator_evaluations, expected.operator_evaluations);
    for (actual, expected) in simulation.accepted.states.iter().zip(&reference.accepted.states) {
        assert!((actual - expected).abs() < 3e-11);
    }
    assert!((recovered.coupled.entropy.entropy_change_j_k
        - expected.coupled.entropy.entropy_change_j_k).abs() < 3e-11);
}

use super::*;
use crate::chemistry::advance_adiabatic;
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::{Boundaries, CaloricDomain, MixtureGas};
use fs_qty::Length;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn mechanism(heat: f64, rate: f64, reference: f64) -> IdealGasKineticsV1 {
    let gas = bank(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", R, 3.5, 0.0, reference - heat),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(rate),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    IdealGasKineticsV1::new(gas, "reactive-flow-isomers", 1, vec![reaction], 1e-12).unwrap()
}
fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10,
        relative_tolerance: 2e-7,
        absolute_amount_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 1.0,
        maximum_steps: 30000,
        ..IsochoricKineticsSolve::default()
    }
}
fn field(
    mechanism: &IdealGasKineticsV1,
    n: usize,
    velocity: f64,
    initial: impl Fn(f64) -> [f64; 2],
) -> MixtureEuler1d {
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let states = (0..n)
        .map(|i| {
            gas.conserved(&initial((i as f64 + 0.5) / n as f64), velocity, 1.0)
                .unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        0.1,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap()
}
fn uniform(mechanism: &IdealGasKineticsV1) -> ReactingEuler1d {
    ReactingEuler1d::new(
        field(mechanism, 4, 0.7, |_| [1.0, 0.0]),
        mechanism.clone(),
        chemistry(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_moving_exothermic_field_restarts_at_actual_duration_and_matches_homogeneous_reactor() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = uniform(&mechanism);
    let before = flow.flow().totals();
    let step = flow
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    assert!(
        step.coupling_retries > 0,
        "heat release must shorten acoustic step: {step:?}"
    );
    assert_eq!(step.elapsed, step.transport.elapsed);
    assert_eq!(flow.flow().time(), step.elapsed);
    assert_eq!(step.chemistry.cell_solves, 8);
    assert!(step.chemistry.entropy_change_j_k > 0.0);
    assert!(step.chemistry.referenced_energy_change_j > 0.0);
    assert!(step.absolute_energy_residual_j.abs() < 2e-12);
    assert!(step.mass_residual_kg.abs() < 1e-13);
    assert!((step.after.momentum - before.momentum).abs() < 1e-13);
    let volume = flow.flow().mesh().volumes()[0];
    let reference = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        step.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..4 {
        let q = flow.flow().cell(i).unwrap();
        let p = flow.flow().gas().primitive(q, || false).unwrap();
        assert!((p.velocity - 0.7).abs() < 1e-12);
        assert!(
            (p.temperature - reference.temperature().value()).abs() < 2e-4,
            "cellT={},reactorT={}",
            p.temperature,
            reference.temperature().value()
        );
        assert!(
            (q.partial_densities[1] - reference.amounts()[1].value() * R / volume).abs() < 2e-5
        );
        assert!(p.temperature > 5.0);
    }
}

fn reacting_contact(n: usize) -> f64 {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let initial = |x: f64| {
        let perturbation = 0.2 * (2.0 * core::f64::consts::PI * x).sin();
        [0.5 + perturbation, 0.5 - perturbation]
    };
    let flow = field(&mechanism, n, 0.7, initial);
    let mut flow = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let end = 0.1;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let step = flow
            .step(
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                || false,
            )
            .unwrap();
        assert!(step.absolute_energy_residual_j.abs() < 2e-12);
        assert!(step.mass_residual_kg.abs() < 2e-13);
        assert!(
            step.conservation_residual.species_masses[..2]
                .iter()
                .all(|x| x.abs() < 2e-13)
        );
        steps += 1;
        assert!(steps < 1000);
    }
    let mut error = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let exact = 0.5
            + 0.2 * (-4.0_f64 * end).exp() * (2.0 * core::f64::consts::PI * (x - 0.7 * end)).sin();
        let cell = flow.flow().cell(i).unwrap();
        error += (cell.partial_densities[0] - exact).abs() / n as f64;
        let p = flow.flow().gas().primitive(cell, || false).unwrap();
        assert!((p.pressure - 1.0).abs() < 2e-9);
        assert!((p.temperature - 1.0).abs() < 2e-9);
    }
    error
}

#[test]
fn g1_advection_reaction_matches_independent_decay_and_refines() {
    let coarse = reacting_contact(40);
    let fine = reacting_contact(80);
    eprintln!("reacting contact L1: 40={coarse:e},80={fine:e}");
    assert!(coarse < 1e-3);
    assert!(fine < coarse * 0.45, "coarse={coarse:e},fine={fine:e}");
}

#[test]
fn g4_late_cancel_keeps_all_cells_and_clock_and_next_step_recovers() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let mut reference = uniform(&mechanism);
    let mut polls = 0;
    let expected = reference
        .step(Time::new(0.01), ReactiveStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 100);
    let mut flow = uniform(&mechanism);
    let original = flow.accepted.states.clone();
    let mut attempted = 0;
    let error = flow
        .step(Time::new(0.01), ReactiveStepControl::default(), || {
            attempted += 1;
            attempted >= polls - 2
        })
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(flow.accepted.states, original);
    assert_eq!(flow.flow().time().value(), 0.0);
    let recovered = flow
        .step(Time::new(0.01), ReactiveStepControl::default(), || false)
        .unwrap();
    assert_eq!(recovered.elapsed, expected.elapsed);
    for (a, b) in flow.accepted.states.iter().zip(&reference.accepted.states) {
        assert!((a - b).abs() < 2e-11, "reused={a},fresh={b}");
    }
}

#[test]
fn g0_mismatched_sources_and_lost_chemical_energy_precision_are_refused_atomically() {
    let original = mechanism(0.0, 2.0, 0.0);
    let other = mechanism(0.0, 2.0, 1.0);
    assert!(matches!(
        ReactingEuler1d::new(
            field(&original, 4, 0.0, |_| [1.0, 0.0]),
            other,
            chemistry(),
            || false
        ),
        Err(ReactiveError::Invalid(_))
    ));
    let shifted = mechanism(0.0, 2.0, 1e18);
    let mut flow = uniform(&shifted);
    let before = flow.accepted.states.clone();
    let error = flow
        .step(
            Time::new(0.01),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            || false,
        )
        .unwrap_err();
    match error {
        ReactiveError::RetryLimit { last, .. } => {
            assert!(matches!(*last, ReactiveError::EnergyConsistency { .. }))
        }
        e => panic!("expected explicit energy precision refusal, got {e:?}"),
    }
    assert_eq!(flow.accepted.states, before);
    assert_eq!(flow.flow().time().value(), 0.0);
}

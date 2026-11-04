use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

fn mode_error(n: usize) -> f64 {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let flow = field(&mechanism, n, 0.7, |x| {
        let a = 0.6 + 0.2 * (core::f64::consts::TAU * x).sin();
        [a, 1.0 - a]
    });
    let mut molecular = molecular(&flow, 0.3, 0.0, 0.0);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let end = 0.04;
    let report = simulation
        .advance_stabilized_to(
            &mut molecular,
            Time::new(end),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(simulation.flow().time().value(), end);
    assert_eq!(report.final_time.value(), end);
    assert!((report.elapsed.value() - end).abs() < 3e-16);
    assert!(report.coupled_steps > 1);
    assert_eq!(report.chemistry.cell_solves, 2 * n * report.coupled_steps);
    assert!(report.chemistry.accepted_steps > 0);
    assert!(report.chemistry.entropy_change_j_k > 0.0);
    assert!(report.maximum_stages_used > 2);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.mass_residual_kg.abs() < 3e-13);
    assert!(report.absolute_energy_residual_j.abs() < 3e-12);
    assert!(
        report.conservation_residual.species_masses[..2]
            .iter()
            .all(|x| x.abs() < 3e-13)
    );
    // Equal calorics/masses make A<->B reaction and molecular diffusion
    // commute. The mean approaches 1/2 at k_f+k_r=4; the nonzero Fourier
    // mode also decays at the independently known discrete diffusion rate.
    let mean = 0.5 + 0.1 * (-4.0_f64 * end).exp();
    let laplacian = 4.0 * (n * n) as f64 * (core::f64::consts::PI / n as f64).sin().powi(2);
    let amplitude = 0.2 * (-(4.0 + 0.3 * laplacian) * end).exp();
    let mut error = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let q = simulation.flow().cell(i).unwrap();
        let expected = mean + amplitude * (core::f64::consts::TAU * (x - 0.7 * end)).sin();
        error += (q.partial_densities[0] - expected).abs() / n as f64;
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        assert!((p.temperature - 1.0).abs() < 3e-9);
        assert!((p.density - 1.0).abs() < 3e-11);
    }
    error
}

#[test]
fn g1_complete_rkl2_reacting_window_matches_independent_advected_reaction_diffusion_mode() {
    let coarse = mode_error(40);
    let fine = mode_error(80);
    eprintln!("coupled RKL2 Fourier L1: 40={coarse:e}, 80={fine:e}");
    assert!(coarse < 1e-3);
    assert!(fine < 0.45 * coarse);
}

#[test]
fn g1_rkl2_uses_the_original_coupled_path_and_reduces_molecular_rhs_work() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let make = || {
        let flow = field(&mechanism, 96, 0.7, |x| {
            let a = 0.6 + 0.2 * (core::f64::consts::TAU * x).sin();
            [a, 1.0 - a]
        });
        let operator = molecular(&flow, 1.0, 0.0, 0.0);
        (
            ReactingEuler1d::new(flow, mechanism.clone(), chemistry(), || false).unwrap(),
            operator,
        )
    };
    let (mut stabilized, mut op) = make();
    let (mut coarse, mut coarse_operator) = make();
    let (mut original, mut explicit) = make();
    let coarse_report = coarse
        .step_with_stabilized_dissipation(
            &mut coarse_operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    let report = stabilized
        .step_with_stabilized_dissipation(
            &mut op,
            Time::new(0.001),
            ReactiveStepControl::default(),
            StabilizedDissipationControl {
                transport: DissipationControl {
                    maximum_step_s: 0.00025,
                    ..DissipationControl::default()
                },
                ..StabilizedDissipationControl::default()
            },
            || false,
        )
        .unwrap();
    let reference = original
        .step_with_dissipation(
            &mut explicit,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.coupled.elapsed, reference.elapsed);
    assert_eq!(
        report.coupled.chemistry.cell_solves,
        reference.chemistry.cell_solves
    );
    let reference_halves = reference.dissipation.unwrap();
    let explicit_calls =
        2 * (reference_halves.first.accepted_steps + reference_halves.second.accepted_steps);
    assert_eq!(reference_halves.first.rejected_steps, 0);
    assert_eq!(reference_halves.second.rejected_steps, 0);
    assert!(report.operator_evaluations * 2 < explicit_calls);
    let difference = |flow: &ReactingEuler1d| {
        flow.accepted
            .states
            .iter()
            .zip(&original.accepted.states)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
    };
    let coarse_error = difference(&coarse);
    let fine_error = difference(&stabilized);
    let mean_error = stabilized
        .accepted
        .states
        .iter()
        .zip(&original.accepted.states)
        .map(|(a, b)| (a - b).abs())
        .sum::<f64>()
        / stabilized.accepted.states.len() as f64;
    // The maximum macrostep is an accuracy input. Reconstruction introduces
    // fine modes after the first D half, so acceleration alone need not match
    // a finely subcycled SSP step at an arbitrarily small pointwise tolerance.
    eprintln!(
        "coupled SSP discrepancy: coarse max={coarse_error:e}, refined max={fine_error:e}, refined mean={mean_error:e}; RHS coarse={}, refined={}, SSP={explicit_calls}",
        coarse_report.operator_evaluations, report.operator_evaluations
    );
    assert!(fine_error < coarse_error);
    assert!(fine_error < 1e-5);
    assert!(mean_error < 1e-6);
    assert_eq!(
        report.dissipation.first.transport.elapsed.value(),
        report.coupled.elapsed.value() * 0.5
    );
    assert_eq!(
        report.dissipation.second.transport.elapsed.value(),
        report.coupled.elapsed.value() * 0.5
    );
    assert_eq!(
        report.coupled.dissipation.unwrap().first,
        report.dissipation.first.transport
    );
    assert!(report.coupled.absolute_energy_residual_j.abs() < 3e-12);
    let (mut limited, mut limited_operator) = make();
    let initial = limited.accepted.states.clone();
    let first_half_calls = report.dissipation.first.operator_evaluations;
    let error = limited
        .step_with_stabilized_dissipation(
            &mut limited_operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            StabilizedDissipationControl {
                transport: DissipationControl {
                    maximum_step_s: 0.00025,
                    ..DissipationControl::default()
                },
                maximum_operator_evaluations: first_half_calls,
                ..StabilizedDissipationControl::default()
            },
            || false,
        )
        .unwrap_err();
    // The first D half completes; the second cannot reset its RHS budget.
    assert!(
        matches!(error, ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { operator_evaluations, .. }) if operator_evaluations == first_half_calls)
    );
    assert_eq!(limited.accepted.states, initial);
    assert_eq!(limited.flow().time().value(), 0.0);
    eprintln!(
        "coupled RKL2: {} RHS calls vs {explicit_calls} SSP RK2 RHS calls",
        report.operator_evaluations
    );
}

#[test]
fn g4_heat_release_restarts_keep_actual_rkl2_work_and_cannot_reset_the_budget() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut reference = uniform(&mechanism);
    let mut reference_operator = molecular(reference.flow(), 0.05, 0.01, 0.1);
    let report = reference
        .step_with_stabilized_dissipation(
            &mut reference_operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.coupled.coupling_retries > 0);
    let accepted_calls = report.dissipation.first.operator_evaluations
        + report.dissipation.second.operator_evaluations;
    assert!(report.operator_evaluations > accepted_calls);
    assert_eq!(
        report.dissipation.first.transport.elapsed.value(),
        0.5 * report.coupled.elapsed.value()
    );
    assert_eq!(
        report.dissipation.second.transport.elapsed.value(),
        0.5 * report.coupled.elapsed.value()
    );
    let volume = reference.flow().mesh().volumes()[0];
    let reactor = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.coupled.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..reference.flow().cell_count() {
        let p = reference
            .flow()
            .gas()
            .primitive(reference.flow().cell(i).unwrap(), || false)
            .unwrap();
        assert!((p.temperature - reactor.temperature().value()).abs() < 2e-4);
    }
    let mut simulation = uniform(&mechanism);
    let original = simulation.accepted.states.clone();
    let mut operator = molecular(simulation.flow(), 0.05, 0.01, 0.1);
    let error = simulation
        .step_with_stabilized_dissipation(
            &mut operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            StabilizedDissipationControl {
                maximum_operator_evaluations: accepted_calls,
                ..StabilizedDissipationControl::default()
            },
            || false,
        )
        .unwrap_err();
    assert!(
        matches!(error, ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { operator_evaluations, .. }) if operator_evaluations == accepted_calls)
    );
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let recovered = simulation
        .step_with_stabilized_dissipation(
            &mut operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(recovered.coupled.elapsed, report.coupled.elapsed);
    assert!(
        simulation
            .accepted
            .states
            .iter()
            .zip(&reference.accepted.states)
            .all(|(a, b)| (a - b).abs() < 3e-10)
    );
}

#[test]
fn g4_complete_window_rolls_back_late_cancel_and_both_aggregate_budgets_then_recovers() {
    let mechanism = mechanism(2.0, 4.0, 0.0);
    let make = || {
        let flow = field(&mechanism, 24, 0.7, |x| {
            let a = 0.6 + 0.1 * (core::f64::consts::TAU * x).sin();
            [a, 1.0 - a]
        });
        let mut operator = molecular(&flow, 0.2, 0.01, 0.03);
        let mut simulation =
            ReactingEuler1d::new(flow, mechanism.clone(), chemistry(), || false).unwrap();
        simulation
            .step_with_stabilized_dissipation(
                &mut operator,
                Time::new(0.001),
                ReactiveStepControl::default(),
                StabilizedDissipationControl::default(),
                || false,
            )
            .unwrap();
        (simulation, operator)
    };
    let (mut reference, mut reference_op) = make();
    let start = reference.flow().time();
    let target = Time::new(0.031);
    let mut polls = 0;
    let expected = reference
        .advance_stabilized_to(
            &mut reference_op,
            target,
            StabilizedReactiveIntervalControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.coupled_steps > 2);
    let (mut simulation, mut op) = make();
    let original = simulation.accepted.states.clone();
    let original_mesh = simulation.accepted.mesh.clone();
    let error = simulation
        .advance_stabilized_to(
            &mut op,
            target,
            StabilizedReactiveIntervalControl {
                maximum_coupled_attempts: 1,
                ..StabilizedReactiveIntervalControl::default()
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ReactiveError::IntervalWorkLimit {
            attempts: 1,
            completed_steps: 1,
            ..
        }
    ));
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time(), start);
    let limit = expected.operator_evaluations / 2;
    let error = simulation
        .advance_stabilized_to(
            &mut op,
            target,
            StabilizedReactiveIntervalControl {
                dissipation: StabilizedDissipationControl {
                    maximum_operator_evaluations: limit,
                    ..StabilizedDissipationControl::default()
                },
                ..StabilizedReactiveIntervalControl::default()
            },
            || false,
        )
        .unwrap_err();
    assert!(
        matches!(error, ReactiveError::Dissipation(DissipationError::OperatorWorkLimit { operator_evaluations, .. }) if operator_evaluations == limit)
    );
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time(), start);
    let mut calls = 0;
    let error = simulation
        .advance_stabilized_to(
            &mut op,
            target,
            StabilizedReactiveIntervalControl::default(),
            || {
                calls += 1;
                calls >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time(), start);
    assert_eq!(simulation.accepted.mesh, original_mesh);
    let recovered = simulation
        .advance_stabilized_to(
            &mut op,
            target,
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(simulation.flow().time(), target);
    assert_eq!(recovered.coupled_steps, expected.coupled_steps);
    assert!(
        simulation
            .accepted
            .states
            .iter()
            .zip(&reference.accepted.states)
            .all(|(a, b)| (a - b).abs() < 3e-10)
    );
    assert!(recovered.absolute_energy_residual_j.abs() < 3e-11);
}

#[test]
fn g0_complete_window_retains_both_boundary_ledgers_and_original_chemical_energy() {
    let mechanism = mechanism(0.0, 3.0, 0.0);
    let mut flow = field(&mechanism, 16, 0.5, |x| [0.6 + 0.1 * x, 0.4 - 0.1 * x]);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let mut op = molecular(&flow, 0.3, 0.05, 0.1);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let report = simulation
        .advance_stabilized_for(
            &mut op,
            Time::new(0.03),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.coupled_steps > 1);
    assert!(report.chemistry.species_change_kg[0] < 0.0);
    assert!(report.chemistry.species_change_kg[1] > 0.0);
    assert!(report.molecular_boundary_change.momentum < 0.0);
    assert_eq!(report.advection_boundary_change.mass(), 0.0);
    assert_eq!(report.molecular_boundary_change.mass(), 0.0);
    assert_eq!(report.advection_boundary_change.energy, 0.0);
    assert_eq!(report.molecular_boundary_change.energy, 0.0);
    for k in 0..4 {
        let combined =
            report.advection_boundary_change.get(k) + report.molecular_boundary_change.get(k);
        assert!(
            (combined - report.left_boundary_change.get(k) - report.right_boundary_change.get(k))
                .abs()
                < 3e-13
        );
    }
    assert!(report.conservation_residual.momentum.abs() < 3e-12);
    assert!(
        report.conservation_residual.species_masses[..2]
            .iter()
            .all(|x| x.abs() < 3e-13)
    );
    assert!(report.mass_residual_kg.abs() < 3e-13);
    assert!(report.absolute_energy_residual_j.abs() < 3e-12);
}

#[test]
fn g4_final_rkl2_heating_cannot_leave_the_original_chemistry_interval() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let make = || {
        let gas = MixtureGas::new(
            mechanism.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || false,
        )
        .unwrap();
        let states = [0.5, -0.5]
            .into_iter()
            .map(|u| gas.conserved(&[0.6, 0.4], u, 1.0).unwrap())
            .collect();
        let flow = MixtureEuler1d::new(
            gas,
            Length::new(0.5),
            0.1,
            states,
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let op = molecular(&flow, 0.03, 1.0, 0.0);
        (flow, op)
    };
    let temperature = |flow: &MixtureEuler1d| {
        (0..flow.cell_count())
            .map(|i| {
                flow.gas()
                    .primitive(flow.cell(i).unwrap(), || false)
                    .unwrap()
                    .temperature
            })
            .fold(0.0, f64::max)
    };
    let (mut wide, mut wide_op) = make();
    wide_op
        .advance_stabilized_for(
            &mut wide,
            Time::new(0.005),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    wide.step(
        Time::new(0.01),
        ReactiveStepControl::default().transport,
        || false,
    )
    .unwrap();
    let before_last = temperature(&wide);
    wide_op
        .advance_stabilized_for(
            &mut wide,
            Time::new(0.005),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    let after_last = temperature(&wide);
    assert!(after_last > before_last + 1e-5);
    let ceiling = 0.5 * (before_last + after_last);
    let (flow, mut op) = make();
    let mut simulation = ReactingEuler1d::new(
        flow,
        mechanism,
        IsochoricKineticsSolve {
            minimum_temperature_k: 0.5,
            maximum_temperature_k: ceiling,
            ..chemistry()
        },
        || false,
    )
    .unwrap();
    let original = simulation.accepted.states.clone();
    let error = simulation
        .step_with_stabilized_dissipation(
            &mut op,
            Time::new(0.01),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap_err();
    match error {
        ReactiveError::RetryLimit { last, .. } => assert_eq!(
            *last,
            ReactiveError::Invalid("final molecular state is outside chemistry interval")
        ),
        other => panic!("expected final molecular chemistry-domain refusal, got {other:?}"),
    }
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
}

#[test]
fn g0_clock_targets_solve_remainders_and_refuse_material_duration_rounding_or_stagnation() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut simulation = uniform(&mechanism);
    let mut op = molecular(simulation.flow(), 0.01, 0.0, 0.0);
    simulation.accepted.time = 0.1;
    let report = simulation
        .advance_stabilized_for(
            &mut op,
            Time::new(0.2),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(simulation.flow().time().value(), 0.1 + 0.2);
    assert!((report.elapsed.value() - 0.2).abs() < 2e-15);
    let original = simulation.accepted.states.clone();
    let noop = simulation
        .advance_stabilized_to(
            &mut op,
            report.final_time,
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(noop.operator_evaluations, 0);
    assert_eq!(noop.coupled_steps, 0);
    assert_eq!(simulation.accepted.states, original);
    simulation.accepted.time = 1e16;
    let error = simulation
        .advance_stabilized_for(
            &mut op,
            Time::new(3.0),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("reacting duration is not representable at this clock")
    );
    assert_eq!(simulation.flow().time().value(), 1e16);
    assert_eq!(simulation.accepted.states, original);
    let error = simulation
        .advance_stabilized_to(
            &mut op,
            Time::new(1e16 + 2.0),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("reacting time increment cannot advance the clock")
    );
    assert_eq!(simulation.flow().time().value(), 1e16);
    assert_eq!(simulation.accepted.states, original);
    // A coarse clock can progress without representing the actual accepted dt:
    // two approximately 3 s physical steps each round the clock ahead by 4 s.
    // Reaching an 8 s represented target must not publish only 6 s of physics.
    let mut coarse = field(&mechanism, 4, 0.7, |_| [0.6, 0.4]);
    let p = coarse
        .gas()
        .primitive(coarse.cell(0).unwrap(), || false)
        .unwrap();
    let width = 3.0 * (p.velocity.abs() + p.sound_speed)
        / (0.8 * ReactiveStepControl::default().transport.cfl);
    coarse.mesh = crate::multicomponent::Mesh1d::uniform(4, Length::new(width), 0.1).unwrap();
    let mut coarse_op = molecular(&coarse, 0.01, 0.0, 0.0);
    let mut coarse = ReactingEuler1d::new(coarse, mechanism, chemistry(), || false).unwrap();
    coarse.accepted.time = 1e16;
    let original = coarse.accepted.states.clone();
    let error = coarse
        .advance_stabilized_to(
            &mut coarse_op,
            Time::new(1e16 + 8.0),
            StabilizedReactiveIntervalControl::default(),
            || false,
        )
        .unwrap_err();
    assert_eq!(
        error,
        ReactiveError::Invalid("reacting clock span differs from integrated durations")
    );
    assert_eq!(coarse.flow().time().value(), 1e16);
    assert_eq!(coarse.accepted.states, original);
}

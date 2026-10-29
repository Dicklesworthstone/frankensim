use super::*;
use crate::multicomponent::Boundary;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationOperator, TransportCoefficients,
    TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn molecular(
    flow: &MixtureEuler1d,
    diffusion: f64,
    viscosity: f64,
    conductivity: f64,
) -> DissipationOperator<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "reacting-molecular-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(viscosity),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        vec![0.0, diffusion, diffusion, 0.0],
    )
    .unwrap();
    DissipationOperator::new(flow, law, || false).unwrap()
}

fn coupled_mode(n: usize) -> f64 {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let initial = |x: f64| {
        let perturbation = 0.2 * (2.0 * core::f64::consts::PI * x).sin();
        [0.5 + perturbation, 0.5 - perturbation]
    };
    let flow = field(&mechanism, n, 0.7, initial);
    let mut operator = molecular(&flow, 0.05, 0.0, 0.0);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let end = 0.1;
    let mut steps = 0;
    while simulation.flow().time().value() < end {
        let report = simulation
            .step_with_dissipation(
                &mut operator,
                Time::new(end - simulation.flow().time().value()),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let molecular = report.dissipation.expect("both molecular halves retained");
        assert_eq!(
            molecular.first.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert_eq!(
            molecular.second.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert!(report.absolute_energy_residual_j.abs() < 2e-12);
        assert!(report.mass_residual_kg.abs() < 2e-13);
        assert!(
            report.conservation_residual.species_masses[..2]
                .iter()
                .all(|x| x.abs() < 2e-13)
        );
        steps += 1;
        assert!(steps < 1000);
    }
    let mut error = 0.0;
    let wave_number = 2.0 * core::f64::consts::PI;
    let amplitude = 0.2 * (-(4.0 + 0.05 * wave_number * wave_number) * end).exp();
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let exact = 0.5 + amplitude * (wave_number * (x - 0.7 * end)).sin();
        let q = simulation.flow().cell(i).unwrap();
        error += (q.partial_densities[0] - exact).abs() / n as f64;
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        assert!((p.temperature - 1.0).abs() < 3e-9);
    }
    error
}

#[test]
fn g1_full_advection_diffusion_reaction_matches_analytic_mode_and_refines() {
    let coarse = coupled_mode(40);
    let fine = coupled_mode(80);
    eprintln!("advection/diffusion/reaction L1:40={coarse:e},80={fine:e}");
    assert!(coarse < 1e-3);
    assert!(fine < 0.45 * coarse, "coarse={coarse:e},fine={fine:e}");
}

#[test]
fn g1_heat_release_restart_replays_both_molecular_half_durations() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut simulation = uniform(&mechanism);
    let mut operator = molecular(simulation.flow(), 0.05, 0.01, 0.1);
    let report = simulation
        .step_with_dissipation(
            &mut operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.coupling_retries > 0);
    assert_eq!(simulation.flow().time(), report.elapsed);
    let halves = report.dissipation.unwrap();
    assert_eq!(halves.first.elapsed.value(), report.elapsed.value() * 0.5);
    assert_eq!(halves.second.elapsed.value(), report.elapsed.value() * 0.5);
    let volume = simulation.flow().mesh().volumes()[0];
    let reference = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..simulation.flow().cell_count() {
        let q = simulation.flow().cell(i).unwrap();
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        assert!((p.temperature - reference.temperature().value()).abs() < 2e-4);
    }
    assert!(report.absolute_energy_residual_j.abs() < 2e-12);
}

#[test]
fn g0_closed_wall_molecular_traction_enters_coupled_momentum_balance_without_work() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut flow = field(&mechanism, 8, 0.5, |_| [0.6, 0.4]);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let mut operator = molecular(&flow, 0.03, 0.02, 0.04);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let report = simulation
        .step_with_dissipation(
            &mut operator,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let halves = report.dissipation.unwrap();
    assert!(halves.first.boundary_change.momentum < 0.0);
    assert_eq!(halves.first.boundary_change.energy, 0.0);
    assert_eq!(halves.second.boundary_change.energy, 0.0);
    assert!(report.conservation_residual.momentum.abs() < 2e-12);
    assert!(report.absolute_energy_residual_j.abs() < 2e-12);
    assert!(report.mass_residual_kg.abs() < 2e-13);
}

#[test]
fn g4_final_viscous_heating_respects_a_narrow_chemistry_interval() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let make = || {
        let gas = MixtureGas::new(
            mechanism.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || false,
        )
        .unwrap();
        let states = [0.5, -0.5]
            .iter()
            .map(|u| gas.conserved(&[0.6, 0.4], *u, 1.0).unwrap())
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
        let operator = molecular(&flow, 0.03, 1.0, 0.0);
        (flow, operator)
    };
    let maximum_temperature = |flow: &MixtureEuler1d| {
        (0..flow.cell_count())
            .map(|i| {
                flow.gas()
                    .primitive(flow.cell(i).unwrap(), || false)
                    .unwrap()
                    .temperature
            })
            .fold(0.0, f64::max)
    };
    let duration = Time::new(0.01);
    let (mut wide, mut wide_operator) = make();
    // The inert mechanism makes both chemistry halves exact no-ops. Isolate
    // the extra heat produced by the final D half using these physical stages.
    wide_operator
        .advance_for(
            &mut wide,
            Time::new(0.5 * duration.value()),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let transport = wide
        .step(duration, ReactiveStepControl::default().transport, || false)
        .unwrap();
    assert_eq!(transport.elapsed, duration);
    let before_temperature = maximum_temperature(&wide);
    wide_operator
        .advance_for(
            &mut wide,
            Time::new(0.5 * duration.value()),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let final_temperature = maximum_temperature(&wide);
    assert!(final_temperature > before_temperature + 1e-5);
    // Admit both chemistry endpoints while excluding the additional heat
    // produced by the final molecular half step of the original duration.
    let ceiling = 0.5 * (before_temperature + final_temperature);
    let controls = IsochoricKineticsSolve {
        minimum_temperature_k: 0.5,
        maximum_temperature_k: ceiling,
        ..chemistry()
    };
    let (flow, mut operator) = make();
    let mut simulation = ReactingEuler1d::new(flow, mechanism.clone(), controls, || false).unwrap();
    let original = simulation.accepted.states.clone();
    let error = simulation
        .step_with_dissipation(
            &mut operator,
            duration,
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl::default(),
            || false,
        )
        .unwrap_err();
    match error {
        ReactiveError::RetryLimit { attempts, last } => {
            assert_eq!(attempts, 1);
            assert_eq!(
                *last,
                ReactiveError::Invalid("final molecular state is outside chemistry interval")
            );
        }
        other => panic!("expected final molecular domain refusal, got {other:?}"),
    }
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let recovered = simulation
        .step_with_dissipation(
            &mut operator,
            duration,
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert!(recovered.coupling_retries > 0);
    assert!(recovered.elapsed.value() < duration.value());
    assert_eq!(simulation.flow().time(), recovered.elapsed);
    let temperature = maximum_temperature(simulation.flow());
    assert!(temperature <= ceiling);
    assert!(recovered.absolute_energy_residual_j.abs() < 2e-12);
    assert!(recovered.conservation_residual.momentum.abs() < 2e-12);
}

#[test]
fn g4_late_molecular_cancellation_keeps_the_entire_coupled_field_and_clock() {
    let mechanism = mechanism(4.0, 30.0, 0.0);
    let make = || {
        let flow = field(&mechanism, 8, 0.7, |x| {
            let y = 0.6 + 0.1 * (2.0 * core::f64::consts::PI * x).sin();
            [y, 1.0 - y]
        });
        let operator = molecular(&flow, 0.04, 0.005, 0.01);
        (
            ReactingEuler1d::new(flow, mechanism.clone(), chemistry(), || false).unwrap(),
            operator,
        )
    };
    let (mut reference, mut reference_operator) = make();
    let mut polls = 0;
    let expected = reference
        .step_with_dissipation(
            &mut reference_operator,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 100);
    let (mut simulation, mut operator) = make();
    let original = simulation.accepted.states.clone();
    let mut calls = 0;
    let error = simulation
        .step_with_dissipation(
            &mut operator,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || {
                calls += 1;
                calls >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let recovered = simulation
        .step_with_dissipation(
            &mut operator,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(recovered.elapsed, expected.elapsed);
    for (a, b) in simulation
        .accepted
        .states
        .iter()
        .zip(&reference.accepted.states)
    {
        assert!((a - b).abs() < 2e-10, "reused={a},fresh={b}");
    }
}

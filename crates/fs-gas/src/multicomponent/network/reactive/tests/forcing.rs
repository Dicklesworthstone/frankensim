//! Physical reaction/acceleration and molecular-source composition checks.

use super::*;
use super::molecular::operator;
use crate::body_force::{BodyForce, ForceControl};
use crate::multicomponent::dissipation::{DissipationControl, DissipationError, StabilizedDissipationControl};
use crate::multicomponent::network::{NetworkBodyForces, reactive::NetworkSources};

const FORCE: ForceControl = ForceControl { maximum_mach: 100.0 };

fn energy(flow: &ReactingNetwork1d) -> f64 {
    let mut result = 0.0;
    for branch in flow.branches() {
        for (i, volume) in branch.mesh().volumes().iter().enumerate() {
            result += volume * flow.gas().absolute_energy_density(branch.cell(i).unwrap()).unwrap();
        }
    }
    for plenum in flow.plenums() {
        result += plenum.volume_m3() * flow.gas().absolute_energy_density(plenum.state()).unwrap();
    }
    result
}

#[test]
fn physical_force_history_and_chemistry_replay_the_shortened_common_duration() {
    let rate = 2.0;
    let mechanism = mechanism(0.0, rate, 0.0);
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)), || false,
    ).unwrap();
    let initial = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let branch = MixtureEuler1d::new(
        gas, Length::new(0.125), 0.1, vec![initial; 8],
        Boundaries::Periodic, || false,
    ).unwrap();
    let network = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let mut flow = ReactingNetwork1d::new(network, mechanism, chemistry(), || false).unwrap();
    let forces = [BodyForce::UniformAccelerationRamp {
        reference_time: 0.0, acceleration: 1000.0, jerk: -10_000.0,
    }];
    let initial_energy = energy(&flow);
    let mut work = 0.0;
    for call in 0..2 {
        let report = flow.step_with_body_forces(
            Time::new(0.02), ReactiveStepControl::default(),
            NetworkBodyForces { branches: &forces, controls: FORCE }, || false,
        ).unwrap();
        if call == 0 { assert!(report.coupling_retries > 0); }
        assert_eq!(report.elapsed, report.transport.elapsed);
        assert_eq!(report.chemistry.cell_solves, 16);
        let t = flow.time().value();
        let expected_velocity = 1000.0 * t - 5000.0 * t * t;
        let expected_a = 0.5 + 0.5 * (-2.0 * rate * t).exp();
        for branch in flow.branches() {
            assert_eq!(branch.time(), flow.time());
            for i in 0..branch.cell_count() {
                let q = branch.cell(i).unwrap();
                let p = flow.gas().primitive(q, || false).unwrap();
                close(p.velocity, expected_velocity, 3e-10);
                close(p.temperature, 1.0, 3e-8);
                close(q.partial_densities[0], expected_a, 3e-5);
                close(q.density(), 1.0, 3e-12);
            }
        }
        work += report.transport.body_force_work_j;
        close(energy(&flow) - initial_energy, work, 3e-9);
        close(report.absolute_energy_residual_j, 0.0, 3e-9);
    }
}

#[test]
fn accelerated_reacting_molecular_ports_baths_and_compression_share_one_physical_balance() {
    for stabilized in [false, true] {
        let mechanism = mechanism(5.0, 20.0, 0.0);
        let raw = network(&mechanism, 4, true, 1.4);
        let bath = raw.gas().conserved(&[1.1, 0.2], 0.0, 1.5).unwrap();
        let reservoirs = [[Some(bath.as_cell()), None], [None; 2]];
        let mut diffusion = operator(&raw, 0.02, 0.015, 0.01);
        let mut flow = ReactingNetwork1d::new(raw, mechanism, chemistry(), || false).unwrap();
        let before = energy(&flow);
        let forces = [BodyForce::UniformAcceleration(1.2), BodyForce::UniformForceDensity(-0.7)];
        let forcing = NetworkBodyForces { branches: &forces, controls: FORCE };
        let report = if stabilized {
            let accepted = flow.step_with_stabilized_dissipation_and_body_forces(
                &mut diffusion, Time::new(0.001), ReactiveStepControl::default(),
                StabilizedDissipationControl::default(), forcing, NetworkSources::default(),
                Some(&reservoirs), Some(&[-0.001]), || false,
            ).unwrap();
            assert!(accepted.operator_evaluations >= accepted.dissipation.first.operator_evaluations
                + accepted.dissipation.second.operator_evaluations);
            accepted.coupled
        } else {
            flow.step_with_dissipation_and_body_forces(
                &mut diffusion, Time::new(0.001), ReactiveStepControl::default(),
                DissipationControl::default(), forcing, NetworkSources::default(),
                Some(&reservoirs), Some(&[-0.001]), || false,
            ).unwrap()
        };
        let actual = energy(&flow) - before;
        let expected = report.boundary_change.absolute_energy(flow.gas()).unwrap()
            + report.transport.boundary_work_on_gas_j + report.transport.body_force_work_j;
        close(actual, expected, 3e-9);
        close(report.absolute_energy_residual_j, 0.0, 3e-9);
        close(report.mass_residual_kg, 0.0, 3e-12);
        assert!(report.transport.body_force_work_j.abs() > 1e-10);
        assert!(report.transport.boundary_work_on_gas_j > 0.0);
        assert!(report.chemistry.referenced_energy_change_j.abs() > 1e-10);
        let d = report.dissipation.unwrap();
        assert_eq!(d.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(d.second.elapsed, d.first.elapsed);
        assert_eq!(diffusion.last_report(), Some(d.second));
        assert_eq!(flow.time(), report.elapsed);
    }
}

#[test]
fn forced_reacting_rkl_budget_and_late_cancellation_preserve_every_accepted_owner() {
    let mechanism = mechanism(5.0, 20.0, 0.0);
    let make = || {
        let raw = network(&mechanism, 4, true, 1.4);
        let op = operator(&raw, 0.02, 0.015, 0.01);
        (ReactingNetwork1d::new(raw, mechanism.clone(), chemistry(), || false).unwrap(), op)
    };
    let forces = [BodyForce::UniformAcceleration(1.2); 2];
    let forcing = NetworkBodyForces { branches: &forces, controls: FORCE };
    let (mut probe, mut op) = make();
    let mut count = 0;
    let expected = probe.step_with_stabilized_dissipation_and_body_forces(
        &mut op, Time::new(0.001), ReactiveStepControl::default(),
        StabilizedDissipationControl::default(), forcing, NetworkSources::default(),
        None, None, || { count += 1; false },
    ).unwrap();
    let (mut refused, mut op) = make();
    let before = snapshot(&refused);
    let operator_before = op.last_report();
    let error = refused.step_with_stabilized_dissipation_and_body_forces(
        &mut op, Time::new(0.001), ReactiveStepControl::default(),
        StabilizedDissipationControl {
            maximum_operator_evaluations: expected.dissipation.first.operator_evaluations,
            ..StabilizedDissipationControl::default()
        }, forcing, NetworkSources::default(), None, None, || false,
    ).unwrap_err();
    assert!(matches!(error, ReactiveNetworkError::Dissipation(DissipationError::OperatorWorkLimit { .. })));
    assert_eq!(snapshot(&refused), before);
    assert_eq!(op.last_report(), operator_before);
    let (mut cancelled, mut op) = make();
    let before = snapshot(&cancelled);
    let mut visits = 0;
    let result = cancelled.step_with_stabilized_dissipation_and_body_forces(
        &mut op, Time::new(0.001), ReactiveStepControl::default(),
        StabilizedDissipationControl::default(), forcing, NetworkSources::default(),
        None, None, || { visits += 1; visits == count },
    );
    assert_eq!(result, Err(ReactiveNetworkError::Cancelled));
    assert_eq!(snapshot(&cancelled), before);
    assert_eq!(op.last_report(), None);
}

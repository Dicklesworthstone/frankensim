//! Physical heat-mode and atomic-refusal checks for the existing planar owner.
use super::*;
use crate::compressible::molecular::MolecularStabilization;
use fs_gas::multicomponent::dissipation::DissipationError;

#[test]
fn rkl_fourier_mode_converges_in_time_and_reduces_actual_diffusion_work() {
    let n = 16;
    let w = 2.0 * std::f64::consts::PI;
    let mode = |x: [f64; 2]| (w * x[0]).sin() * (w * x[1]).cos();
    let initial = |x| [1.0, 0.0, 0.0, 2.0 + 0.2 * mode(x)];
    let end = 0.2;
    // Independent eigenvalue of the original centered Cartesian Laplacian;
    // exact time evolution isolates temporal error from spatial truncation.
    let lambda = 8.0 * (0.05 / 2.5) * (n as f64).powi(2)
        * (std::f64::consts::PI / n as f64).sin().powi(2);
    let exact_amplitude = 0.2 * (-lambda * end).exp();
    let mut errors = Vec::new();
    let mut first_work = 0;
    for step in [0.2, 0.1, 0.05] {
        let mut flow = planar(n, initial);
        let mut op = MolecularTransport2d::new(
            &flow, law(&flow, 0.0, 0.0, 0.05), || false,
        ).unwrap();
        let report = op.advance_stabilized(
            &mut flow, Time::new(end),
            MolecularControl { maximum_step_s: step, ..MolecularControl::default() },
            MolecularStabilization::default(), || false,
        ).unwrap();
        let mut projection = 0.0;
        let mut norm = 0.0;
        for (i, q) in flow.states().iter().enumerate() {
            let phi = mode(flow.grid().cell_center(i).unwrap());
            projection += (q.primitive(flow.gas(), || false).unwrap().temperature - 2.0) * phi;
            norm += phi * phi;
            assert_eq!(q.density, 1.0);
            assert_eq!(q.momentum, [0.0; 2]);
        }
        errors.push((projection / norm - exact_amplitude).abs());
        if first_work == 0 { first_work = report.operator_evaluations; }
        assert!(report.maximum_stages_used > 2);
        assert!(report.transport.maximum_parabolic_number > 0.4);
        assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
        assert_eq!(report.transport.rejected_steps, 0);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(flow.last_step(), None);
        for residual in report.transport.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
        assert!(op.stabilized_workspace_bytes() >= op.required_stabilized_workspace_bytes().unwrap());
    }
    assert!(errors[0] < 4e-4, "{errors:?}");
    assert!(errors[1] < 0.30 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.30 * errors[1], "{errors:?}");
    let mut reference = planar(n, initial);
    let mut op = MolecularTransport2d::new(
        &reference, law(&reference, 0.0, 0.0, 0.05), || false,
    ).unwrap();
    let ssp = op.advance(&mut reference, Time::new(end), MolecularControl::default(), || false).unwrap();
    // Even two calls per SSP step is a lower bound on its actual source work.
    assert!(first_work * 2 < 2 * ssp.accepted_steps);
}

#[test]
fn rkl_work_storage_and_late_cancellation_preserve_the_complete_planar_state() {
    let initial = |x: [f64; 2]| [1.0, 0.0, 0.0, 2.0 + 0.2 * (6.0 * x[0]).sin()];
    let mut flow = planar(6, initial);
    let mut op = MolecularTransport2d::new(&flow, law(&flow, 0.0, 0.0, 0.05), || false).unwrap();
    let original = flow.states().to_vec();
    let controls = MolecularControl::default();
    let selected = MolecularStabilization::default();
    let duration = Time::new(0.02);
    let required = op.required_stabilized_workspace_bytes().unwrap();
    assert!(matches!(
        op.advance_stabilized(
            &mut flow, duration, controls,
            MolecularStabilization { maximum_extra_workspace_bytes: required - 1, ..selected },
            || false,
        ),
        Err(MolecularError::StabilizedWorkspace { .. })
    ));
    assert_eq!(op.stabilized_workspace_bytes(), 0);
    assert_eq!(flow.states(), original);
    assert!(matches!(
        op.advance_stabilized(
            &mut flow, duration, controls,
            MolecularStabilization { maximum_operator_evaluations: 2, ..selected }, || false,
        ),
        Err(MolecularError::Law(DissipationError::OperatorWorkLimit {
            operator_evaluations: 2, ..
        }))
    ));
    assert_eq!(flow.states(), original);
    assert_eq!(flow.time().value(), 0.0);
    let impulses = SIDES.map(|s| op.boundary_impulses(s).to_vec());

    let mut probe = planar(6, initial);
    let mut native = MolecularTransport2d::new(&probe, law(&probe, 0.0, 0.0, 0.05), || false).unwrap();
    let mut polls = 0;
    let expected = native.advance_stabilized(&mut probe, duration, controls, selected, || {
        polls += 1;
        false
    }).unwrap();
    let mut calls = 0;
    assert_eq!(
        op.advance_stabilized(&mut flow, duration, controls, selected, || {
            calls += 1;
            calls == polls
        }),
        Err(MolecularError::Cancelled)
    );
    assert_eq!(flow.states(), original);
    for (i, side) in SIDES.into_iter().enumerate() {
        assert_eq!(op.boundary_impulses(side), impulses[i]);
    }
    let report = op.advance_stabilized(&mut flow, duration, controls, selected, || false).unwrap();
    assert_eq!(report, expected);
    assert_eq!(flow.states(), probe.states());
}

#[test]
fn rkl_coupled_heat_uses_the_common_physical_interval_and_complete_balance() {
    let initial = planar(12, |x| [
        1.0, 0.0, 0.0,
        2.0 + 0.1 * (2.0 * std::f64::consts::PI * x[0]).sin(),
    ]);
    let source = law(&initial, 0.01, 0.0, 5.0);
    let mut flow = MolecularEuler2d::new(initial, source, || false).unwrap();
    let report = flow.step_stabilized(
        Time::new(0.001), MolecularStepControl::default(),
        MolecularStabilization::default(), || false,
    ).unwrap();
    let physical = report.transport;
    assert!(report.maximum_stages_used > 2);
    assert_eq!(physical.coupling_retries, 0);
    assert_eq!(physical.first.before, physical.before);
    assert_eq!(physical.first.after, physical.advection.before);
    assert_eq!(physical.second.before, physical.advection.after);
    assert_eq!(physical.second.after, physical.after);
    assert_eq!(physical.first.elapsed.value() * 2.0, physical.elapsed.value());
    assert_eq!(physical.second.elapsed, physical.first.elapsed);
    assert_eq!(physical.advection.elapsed, physical.elapsed);
    assert_eq!(flow.flow().time(), physical.elapsed);
    assert_eq!(flow.last_step().copied(), Some(physical));
    for residual in physical.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
}

#[test]
fn rkl_coupled_work_is_shared_between_halves_and_refusal_preserves_history() {
    let mut flow = coupled();
    let controls = MolecularStepControl::default();
    let selected = MolecularStabilization::default();
    flow.step_stabilized(Time::new(0.001), controls, selected, || false).unwrap();
    let states = flow.flow().states().to_vec();
    let time = flow.flow().time();
    let previous = flow.last_step().copied();
    let loads = SIDES.map(|s| flow.boundary_impulses(s).to_vec());
    // Each s=2 half requires three evaluations, so five cannot finish both.
    assert!(matches!(
        flow.step_stabilized(
            Time::new(0.001), controls,
            MolecularStabilization { maximum_operator_evaluations: 5, ..selected }, || false,
        ),
        Err(MolecularError::Law(DissipationError::OperatorWorkLimit {
            operator_evaluations: 5, ..
        }))
    ));
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step().copied(), previous);
    for (i, side) in SIDES.into_iter().enumerate() {
        assert_eq!(flow.boundary_impulses(side), loads[i]);
    }
    let report = flow.step_stabilized(
        Time::new(0.001), controls,
        MolecularStabilization { maximum_operator_evaluations: 6, ..selected }, || false,
    ).unwrap();
    assert_eq!(report.operator_evaluations, 6);
    assert_eq!(report.transport.coupling_retries, 0);
}

#[test]
fn rkl_forced_uniform_flow_recovers_exact_body_impulse_and_kinetic_work() {
    let velocity = [0.2, -0.1];
    let acceleration = [0.5, -0.25];
    let initial = planar(4, |_| [1.0, velocity[0], velocity[1], 2.0]);
    let mass = initial.totals().mass;
    let source = law(&initial, 0.02, 0.003, 0.05);
    let mut flow = MolecularEuler2d::new(initial, source, || false).unwrap();
    let mut workspace = forcing_workspace(flow.flow());
    let report = flow.step_stabilized_with_forcing(
        Time::new(0.01), MolecularStepControl::default(),
        MolecularStabilization { maximum_operator_evaluations: 6, ..MolecularStabilization::default() },
        BodyForce2d::UniformAcceleration(acceleration), ForceControl2d { maximum_mach: 2.0 },
        &mut workspace, || false,
    ).unwrap();
    let dt = report.transport.elapsed.value();
    assert_eq!(dt, 0.01);
    assert_eq!(report.operator_evaluations, 6);
    let mut kinetic_change = 0.0;
    for axis in 0..2 {
        let next = velocity[axis] + acceleration[axis] * dt;
        close(report.transport.source_change.impulse[axis], mass * acceleration[axis] * dt, 3e-13);
        kinetic_change += 0.5 * mass * (next * next - velocity[axis] * velocity[axis]);
    }
    close(report.transport.source_change.work, kinetic_change, 3e-13);
    for q in flow.flow().states() {
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.temperature, 2.0, 3e-13);
        for axis in 0..2 {
            close(p.velocity[axis], velocity[axis] + acceleration[axis] * dt, 3e-13);
        }
    }
    for residual in report.transport.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
}

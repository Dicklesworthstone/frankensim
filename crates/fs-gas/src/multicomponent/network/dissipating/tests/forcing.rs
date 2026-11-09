//! Whole molecular/acceleration energy and exact inactive-path checks.

use super::*;
use crate::body_force::{BodyForce, ForceControl};
use crate::multicomponent::network::NetworkBodyForces;

const FORCE: ForceControl = ForceControl { maximum_mach: 100.0 };

#[test]
fn explicit_and_rkl_accelerated_networks_preserve_species_and_close_physical_work() {
    let g = gas(2, true, 5.0);
    let forces = [BodyForce::UniformAccelerationRamp {
        reference_time: 0.0, acceleration: 2.0, jerk: -0.4,
    }];
    let forcing = NetworkBodyForces { branches: &forces, controls: FORCE };
    for stabilized in [false, true] {
        let mut flow = owner(&g);
        let initial = inventories(&flow);
        let mut work = 0.0;
        for _ in 0..6 {
            let report = if stabilized {
                let report = flow.step_stabilized_with_body_forces(
                    Time::new(0.002), StabilizedDissipatingStepControl::default(),
                    forcing, None, NetworkSources::default(), || false,
                ).unwrap();
                assert!(report.operator_evaluations >= report.dissipation.first.operator_evaluations
                    + report.dissipation.second.operator_evaluations);
                report.coupled
            } else {
                flow.step_with_body_forces(
                    Time::new(0.002), DissipatingStepControl::default(), forcing, || false,
                ).unwrap()
            };
            work += report.transport.body_force_work_j;
            assert_eq!(report.elapsed, report.transport.elapsed);
            assert_eq!(report.molecular.first.elapsed.value(), 0.5 * report.elapsed.value());
            assert_eq!(report.molecular.second.elapsed, report.molecular.first.elapsed);
            assert_eq!(report.boundary_change.energy, 0.0);
            close(report.absolute_energy_residual_j, 0.0, 3e-10);
            let after = inventories(&flow);
            for k in 0..2 { close(after.0[k], initial.0[k], 3e-12); }
            close(after.1 - initial.1, work, 3e-10);
        }
        assert!(work.abs() > 1e-8);
    }
}

#[test]
fn inactive_force_returns_to_original_stabilized_owner_after_nonzero_use() {
    let g = gas(2, true, 5.0);
    let mut ordinary = owner(&g);
    let mut forced = owner(&g);
    let active = [BodyForce::UniformAcceleration(1.1)];
    for flow in [&mut ordinary, &mut forced] {
        flow.step_stabilized_with_body_forces(
            Time::new(0.002), StabilizedDissipatingStepControl::default(),
            NetworkBodyForces { branches: &active, controls: FORCE },
            None, NetworkSources::default(), || false,
        ).unwrap();
    }
    let zero = [BodyForce::UniformAcceleration(0.0)];
    let a = ordinary.step_stabilized(
        Time::new(0.002), StabilizedDissipatingStepControl::default(), || false,
    ).unwrap();
    let b = forced.step_stabilized_with_body_forces(
        Time::new(0.002), StabilizedDissipatingStepControl::default(),
        NetworkBodyForces { branches: &zero, controls: FORCE },
        None, NetworkSources::default(), || false,
    ).unwrap();
    assert_eq!(a, b);
    assert_eq!(snapshot(&ordinary), snapshot(&forced));
    assert_eq!(b.coupled.transport.body_force_work_j, 0.0);
    assert!(forced.accepted_branch_transport_reports().unwrap().iter()
        .all(|report| report.source_change == crate::body_force::ForceChange::default()));
}

//! Independent mass/force and finite-stretch oracles for the declared split.
use super::*;
use fs_couple::finite_solid::impulse::{FiniteSolidImpulseControls, FiniteSolidPoint};

fn vector(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        near(actual[axis], expected[axis], tolerance);
    }
}

#[test]
fn original_force_kick_and_kinetic_drift_match_finite_stretch_equations() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let mut u = [[0.0; 3]; 4];
        let mut v = [[0.0; 3]; 4];
        u[1][0] = 0.3;
        v[1][0] = 0.1;
        let force = [(3, 0.2)];
        let body = FiniteSolidBody::new(
            problem(&material, &fixed, &force),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let input = body.native_state().clone();
        let h = 1e-3;
        let mass = 6.0 * VOLUME / 10.0;
        let elastic = VOLUME * nominal_stress(1.3, 1.0);
        let expected_v = 0.1 + h * (0.2 - elastic) / mass;
        let expected_q = 0.3 + h * expected_v;
        let (kicked, kick) = body
            .force_kick_candidate(
                h,
                &force,
                [0.0; 3],
                [0.2, -0.4, 0.1],
                Default::default(),
                cx,
            )
            .unwrap();
        near(kicked.native_state().v[0], expected_v, 2e-12);
        assert_eq!(kicked.native_state().q, input.q);
        assert_eq!(kicked.native_state().a, input.a);
        assert_eq!(kicked.native_state().t, input.t);
        assert_eq!(kicked.native_state().history, input.history);
        near(
            kick.source_impulse_work_j,
            h * 0.2 * 0.5 * (0.1 + expected_v),
            2e-14,
        );
        near(
            kick.internal_impulse_work_j,
            -h * elastic * 0.5 * (0.1 + expected_v),
            2e-13,
        );
        let kicked_state = kicked.native_state().clone();
        let (next, drift) = kicked
            .contact_drift_candidate(&kick, controls(1e-3), cx)
            .unwrap();
        near(next.native_state().q[0], expected_q, 2e-14);
        near(next.native_state().v[0], expected_v, 2e-12);
        near(
            next.native_state().a[0],
            (0.2 - VOLUME * nominal_stress(1.0 + expected_q, 1.0)) / mass,
            2e-10,
        );
        near(
            drift.strain_after_j,
            strain_energy(1.0 + expected_q, 1.0),
            2e-11,
        );
        near(
            drift.source_displacement_work_j,
            0.2 * h * expected_v,
            2e-14,
        );
        near(
            drift.source_splitting_work_j,
            h * 0.2 * 0.5 * (0.1 - expected_v),
            2e-14,
        );
        near(
            drift.elastic_splitting_work_j,
            strain_energy(1.0 + expected_q, 1.0)
                - strain_energy(1.3, 1.0)
                - h * elastic * 0.5 * (0.1 + expected_v),
            2e-11,
        );
        near(
            drift.integration_defect_j,
            0.5 * mass * (expected_v * expected_v - 0.01) + strain_energy(1.0 + expected_q, 1.0)
                - strain_energy(1.3, 1.0)
                - 0.2 * h * expected_v,
            2e-11,
        );
        assert!(drift.absolute_integration_defect_j >= drift.integration_defect_j.abs());
        assert!(drift.time.newton.converged);
        assert!(drift.endpoint_acceleration.converged);
        assert_eq!(next.native_state().history.last(), Some(&drift.time));
        assert_eq!(next.native_state().steps, 1);
        near(next.native_state().t, h, 0.0);
        assert_eq!(body.native_state(), &input);
        assert_eq!(kicked.native_state(), &kicked_state);
    });
}

#[test]
fn uniform_force_drift_retains_consistent_mass_orbital_momentum_and_source_work() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let origin = [0.5, -0.7, 0.3];
        let initial_v = [0.2, -0.3, 0.4];
        let h = 0.01;
        let density_force = [0.6, -1.2, 1.8];
        let acceleration = density_force.map(|force| force / 6.0);
        let velocity: [f64; 3] =
            std::array::from_fn(|axis| initial_v[axis] + h * acceleration[axis]);
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[initial_v; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (kicked, kick) = body
            .force_kick_candidate(h, &[], density_force, origin, Default::default(), cx)
            .unwrap();
        let (next, drift) = kicked
            .contact_drift_candidate(&kick, controls(1e-3), cx)
            .unwrap();
        for node in 0..4 {
            vector(next.velocity_m_s()[node], velocity, 3e-13);
            vector(
                next.displacement_m()[node],
                velocity.map(|value| h * value),
                3e-15,
            );
        }
        vector(
            kick.source_impulse.linear_kg_m_s,
            density_force.map(|value| h * VOLUME * value),
            2e-15,
        );
        let arm = [0.25 - origin[0], 0.25 - origin[1], 0.25 - origin[2]];
        let impulse = density_force.map(|value| h * VOLUME * value);
        let torque = [
            arm[1] * impulse[2] - arm[2] * impulse[1],
            arm[2] * impulse[0] - arm[0] * impulse[2],
            arm[0] * impulse[1] - arm[1] * impulse[0],
        ];
        vector(kick.source_impulse.angular_kg_m2_s, torque, 2e-15);
        vector(drift.after.linear_kg_m_s, velocity, 3e-13);
        vector(drift.linear_residual_kg_m_s, [0.0; 3], 1e-14);
        vector(drift.angular_residual_kg_m2_s, [0.0; 3], 1e-14);
        near(
            drift.source_displacement_work_j,
            h * VOLUME
                * (0..3)
                    .map(|axis| density_force[axis] * velocity[axis])
                    .sum::<f64>(),
            2e-14,
        );
        near(
            drift.source_splitting_work_j,
            -0.5 * h * h * acceleration.iter().map(|a| a * a).sum::<f64>(),
            2e-14,
        );
    });
}

#[test]
fn supported_contact_impulse_can_cancel_load_then_advance_real_time_without_bounce() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let body = FiniteSolidBody::new(
            problem(&material, &fixed, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let force = [(3, -2.0)];
        let h = 0.01;
        let (kicked, kick) = body
            .force_kick_candidate(h, &force, [0.0; 3], [0.0; 3], Default::default(), cx)
            .unwrap();
        // Independent one-DOF zero-restitution wall law: J = -h F. The full
        // coupled contact solver has separate fixtures; this isolates ordering.
        let (contacted, contact) = kicked
            .point_impulse_candidate(
                FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap(),
                [2.0 * h, 0.0, 0.0],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        let (next, drift) = contacted
            .contact_drift_candidate(&kick, controls(0.01), cx)
            .unwrap();
        near(next.native_state().v[0], 0.0, 1e-14);
        near(next.native_state().q[0], 0.0, 1e-16);
        near(next.native_state().t, h, 0.0);
        assert_eq!(next.native_state().steps, 1);
        // Smooth acceleration continues to represent the held load, not an
        // invented zero total acceleration absorbing the impulsive contact law.
        near(next.native_state().a[0], -20.0, 2e-12);
        near(drift.source_displacement_work_j, 0.0, 2e-16);
        near(
            kick.impulse.kinetic_change_j + contact.kinetic_change_j,
            0.0,
            2e-15,
        );
        near(drift.integration_defect_j, 0.002, 2e-15);
        near(drift.elastic_splitting_work_j, 0.0, 1e-14);
    });
}

#[test]
fn late_split_refusal_and_foreign_material_binding_preserve_accepted_state() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let make = || {
            FiniteSolidBody::new(
                problem(&material, &[], &[]),
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap()
        };
        let original = make();
        let (kicked, mut kick) = original
            .force_kick_candidate(
                0.1,
                &[],
                [0.0, -6.0, 0.0],
                [0.0; 3],
                FiniteSolidImpulseControls::default(),
                cx,
            )
            .unwrap();
        let state = kicked.native_state().clone();
        assert!(matches!(
            kicked.contact_drift_candidate(&kick, controls(1e-8), cx),
            Err(FiniteSolidError::Energy { .. })
        ));
        assert_eq!(kicked.native_state(), &state);
        assert!(matches!(
            make().contact_drift_candidate(&kick, controls(1.0), cx),
            Err(FiniteSolidError::Invalid(_))
        ));
        // Public diagnostics are observations, never mutable solver inputs.
        kick.duration_s = f64::NAN;
        kick.source_impulse_work_j = 0.0;
        kick.strain_before_j = 1e99;
        let (next, report) = kicked
            .contact_drift_candidate(&kick, controls(1.0), cx)
            .unwrap();
        near(next.native_state().t, 0.1, 0.0);
        near(report.source_splitting_work_j, -0.005, 2e-14);
        assert_eq!(kicked.native_state(), &state);
        assert!(
            next.contact_drift_candidate(&kick, controls(1.0), cx)
                .is_err()
        );
        gate.request();
        assert!(
            kicked
                .contact_drift_candidate(&kick, controls(1.0), cx)
                .is_err()
        );
        assert_eq!(kicked.native_state(), &state);
    });
}

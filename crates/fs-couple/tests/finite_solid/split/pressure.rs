//! Native uniform pressure in the original supported force/contact/drift map.
use super::*;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use std::sync::Arc;

fn surface(cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap(),
    )
}

#[test]
fn pressure_kick_and_drift_keep_original_force_work_and_endpoint_acceleration() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let mut u = [[0.0; 3]; 4];
        u[1][0] = 0.3;
        let mut v = [[0.0; 3]; 4];
        v[1][0] = 0.1;
        let force = [(3, 0.2)];
        let body = FiniteSolidBody::new_with_pressure(
            problem(&material, &fixed, &force),
            surface(cx),
            3.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let initial = body.native_state().clone();
        let h = 1e-3;
        let mass = 6.0 * VOLUME / 10.0;
        // This uniaxial tet has constant dV/dq=1/6 at its only free DOF.
        let pressure_force = -3.0 / 6.0;
        let elastic = VOLUME * nominal_stress(1.3, 1.0);
        let expected_v = 0.1 + h * (0.2 + pressure_force - elastic) / mass;
        let expected_q = 0.3 + h * expected_v;
        let (kicked, mut kick) = body
            .force_kick_candidate_with_pressure(
                h,
                &force,
                [0.0; 3],
                3.0,
                [0.2, -0.4, 0.1],
                Default::default(),
                cx,
            )
            .unwrap();
        near(kicked.native_state().v[0], expected_v, 2e-12);
        assert_eq!(kicked.native_state().q, initial.q);
        assert_eq!(kicked.native_state().a, initial.a);
        assert_eq!(kicked.native_state().t, initial.t);
        assert_eq!(kicked.native_state().history, initial.history);
        let pressure_kick_work = h * pressure_force * 0.5 * (0.1 + expected_v);
        near(kick.pressure_impulse_work_j, pressure_kick_work, 2e-14);
        near(
            kick.source_impulse_work_j,
            h * (0.2 + pressure_force) * 0.5 * (0.1 + expected_v),
            2e-14,
        );
        vector(
            kick.source_impulse.linear_kg_m_s,
            [0.2 * h, 0.0, 0.0],
            1e-14,
        );
        let (next, drift) = kicked
            .contact_drift_candidate(&kick, controls(1e-3), cx)
            .unwrap();
        near(next.native_state().q[0], expected_q, 2e-14);
        near(next.native_state().v[0], expected_v, 2e-12);
        near(
            next.native_state().a[0],
            (0.2 + pressure_force - VOLUME * nominal_stress(1.0 + expected_q, 1.0)) / mass,
            2e-10,
        );
        let pressure_work = pressure_force * h * expected_v;
        near(drift.pressure_work_j, pressure_work, 2e-14);
        near(
            drift.source_displacement_work_j,
            (0.2 + pressure_force) * h * expected_v,
            2e-14,
        );
        near(
            drift.pressure_splitting_work_j,
            pressure_kick_work - pressure_work,
            2e-14,
        );
        near(
            drift.integration_defect_j,
            0.5 * mass * (expected_v * expected_v - 0.01) + strain_energy(1.0 + expected_q, 1.0)
                - strain_energy(1.3, 1.0)
                - (0.2 + pressure_force) * h * expected_v,
            2e-11,
        );
        assert!(drift.absolute_integration_defect_j >= drift.pressure_splitting_work_j.abs());
        vector(drift.linear_residual_kg_m_s, [0.0; 3], 1e-13);
        vector(drift.angular_residual_kg_m2_s, [0.0; 3], 1e-13);
        assert!(drift.endpoint_acceleration.converged);
        assert_eq!(next.native_state().history.last(), Some(&drift.time));
        // Public diagnostics never choose the pressure retained by continuation.
        kick.pressure_impulse_work_j = f64::NAN;
        let replay = kicked
            .contact_drift_candidate(&kick, controls(1e-3), cx)
            .unwrap();
        assert_eq!(replay.0.native_state(), next.native_state());
        assert_eq!(
            replay.1.pressure_work_j.to_bits(),
            drift.pressure_work_j.to_bits()
        );
        assert_eq!(body.native_state(), &initial);
    });
}

#[test]
fn zero_pressure_preserves_legacy_numerics_and_invalid_pressure_is_atomic() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let u = NODES.map(|x| x.map(|v| 0.01 * v));
        let v = [[0.2, -0.1, 0.3]; 4];
        let plain = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let bound = FiniteSolidBody::new_with_pressure(
            problem(&material, &[], &[]),
            surface(cx),
            0.0,
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = bound.native_state().clone();
        let old = plain
            .force_kick_candidate(
                1e-4,
                &[],
                [0.1, -0.2, 0.3],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        let zero = bound
            .force_kick_candidate_with_pressure(
                1e-4,
                &[],
                [0.1, -0.2, 0.3],
                0.0,
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(old.0.native_state(), zero.0.native_state());
        let old = old
            .0
            .contact_drift_candidate(&old.1, controls(1e-3), cx)
            .unwrap();
        let zero = zero
            .0
            .contact_drift_candidate(&zero.1, controls(1e-3), cx)
            .unwrap();
        assert_eq!(old.0.native_state(), zero.0.native_state());
        assert_eq!(
            old.1.source_displacement_work_j.to_bits(),
            zero.1.source_displacement_work_j.to_bits()
        );
        assert_eq!(
            old.1.absolute_integration_defect_j.to_bits(),
            zero.1.absolute_integration_defect_j.to_bits()
        );
        assert!(
            bound
                .force_kick_candidate(1e-4, &[], [0.0; 3], [0.0; 3], Default::default(), cx)
                .is_err()
        );
        assert!(
            bound
                .force_kick_candidate_with_pressure(
                    1e-4,
                    &[],
                    [0.0; 3],
                    f64::NAN,
                    [0.0; 3],
                    Default::default(),
                    cx
                )
                .is_err()
        );
        assert!(
            plain
                .force_kick_candidate_with_pressure(
                    1e-4,
                    &[],
                    [0.0; 3],
                    3.0,
                    [0.0; 3],
                    Default::default(),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(
            bound
                .force_kick_candidate_with_pressure(
                    1e-4,
                    &[],
                    [0.0; 3],
                    3.0,
                    [0.0; 3],
                    Default::default(),
                    cx
                )
                .is_err()
        );
        assert_eq!(bound.native_state(), &before);
    });
}

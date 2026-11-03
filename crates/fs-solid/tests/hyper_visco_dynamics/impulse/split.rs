//! Native Piola kick, intervening impulse, and retained Maxwell drift.
use super::*;
use fs_solid::hyper3::visco::dynamics::FiniteViscoContactDriftControls;

fn drift_controls() -> FiniteViscoContactDriftControls {
    FiniteViscoContactDriftControls {
        material_balance_tolerance_j: 2e-9,
        maximum_integration_defect_j: 0.1,
    }
}

// Unit-tetrahedron virtual work, using the ACTUAL relaxed Piola field.
fn internal(motion: &FiniteViscoTetMotion<'_>) -> [[f64; 3]; 4] {
    let gradient = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let piola = motion.solid().element_first_piola_pa()[0];
    std::array::from_fn(|node| {
        std::array::from_fn(|axis| {
            VOLUME
                * (0..3)
                    .map(|j| piola[3 * axis + j] * gradient[node][j])
                    .sum::<f64>()
        })
    })
}

#[test]
fn supported_graded_force_contact_drift_retains_native_memory_and_private_phase_binding() {
    with_cx(|cx| {
        let law = material(true);
        let motion = FiniteViscoTetMotion::try_new_with_nodal_density(
            relaxed(&law, cx),
            &RHO,
            &[[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let old_states = motion.solid().element_states().to_vec();
        let old_u = motion.solid().displacement_m().to_vec();
        let held: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
        let dt = 0.01;
        let origin = [0.4, -0.3, 0.2];
        let piola_force = internal(&motion);
        let body = [0.4, -0.2, 0.1];
        let mut external = [body.map(|f| VOLUME * f / 4.0); 4];
        let load = piola_force[1][0] + 0.03 / dt - external[1][0];
        external[1][0] += load;
        let (kicked, kick) = motion
            .force_kick_candidate(
                dt,
                &held,
                &[(3, 0.5 * load), (3, 0.5 * load)],
                body,
                origin,
                impulse_controls(),
                cx,
            )
            .unwrap();
        // M_11=1.5, from cubic graded-density integration, not mean density.
        close(kicked.velocity_mps()[1][0], 0.03 / mass(1, 1), 1e-13);
        for &(dof, support) in &kick.impulse.support_impulses_kg_m_s {
            let node = dof / 3;
            let axis = dof % 3;
            let momentum = if axis == 0 { mass(node, 1) * 0.02 } else { 0.0 };
            close(
                support,
                momentum - dt * (external[node][axis] - piola_force[node][axis]),
                1e-12,
            );
        }
        close(
            kick.source_impulse_work_j,
            dt * external[1][0] * 0.01,
            1e-13,
        );
        close(
            kick.internal_impulse_work_j,
            -dt * piola_force[1][0] * 0.01,
            1e-13,
        );
        assert_eq!(kicked.solid().element_states(), old_states);
        assert_eq!(kicked.solid().time_s(), HOLD);
        let large = 18014398509481984.0; // 2^54: ordinary addition loses 0.1.
        let (_, resultant) = motion
            .force_kick_candidate(
                dt,
                &held,
                &[(3, large), (3, 0.1), (3, -large)],
                [0.0; 3],
                origin,
                impulse_controls(),
                cx,
            )
            .unwrap();
        close(resultant.source_impulse.linear_kg_m_s[0], 0.1 * dt, 1e-15);
        assert!(
            motion
                .contact_drift_candidate(&kick, drift_controls(), cx)
                .is_err()
        );
        let (second_kick, _) = kicked
            .force_kick_candidate(dt, &held, &[], [0.0; 3], origin, impulse_controls(), cx)
            .unwrap();
        assert!(
            second_kick
                .contact_drift_candidate(&kick, drift_controls(), cx)
                .is_err()
        );

        let target_v = 0.015;
        let contact = kicked
            .impulse_candidate(
                &held,
                &[(3, mass(1, 1) * (target_v - 0.02))],
                origin,
                impulse_controls(),
                cx,
            )
            .unwrap();
        let (next, drift) = contact
            .motion
            .contact_drift_candidate(&kick, drift_controls(), cx)
            .unwrap();
        close(next.solid().time_s(), HOLD + dt, 1e-15);
        close(
            next.solid().displacement_m()[1][0],
            dt.mul_add(target_v, old_u[1][0]),
            1e-14,
        );
        for &dof in &held {
            assert_eq!(
                next.solid().displacement_m()[dof / 3][dof % 3],
                old_u[dof / 3][dof % 3]
            );
        }
        assert_eq!(next.velocity_mps(), contact.motion.velocity_mps());
        assert_eq!(next.reference_nodal_density_kg_m3(), Some(RHO.as_slice()));
        let mut f = old_states[0].deformation_gradient;
        f[0] += dt * target_v;
        let point = law.step_response(&old_states[0], &f, dt).unwrap();
        close(
            next.solid().stored_energy_j(),
            VOLUME * point.state.stored,
            2e-11,
        );
        close(drift.dissipated_j, VOLUME * point.dissipated_pa, 2e-11);
        close(
            drift.native_loading_work_j,
            VOLUME * point.loading_work_pa,
            2e-11,
        );
        assert!(drift.dissipated_j > 0.0);
        for (actual, expected) in next.solid().element_states()[0]
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(point.state.inverse_viscous_cauchy_green.iter().flatten())
        {
            close(*actual, *expected, 1e-12);
        }
        close(
            drift.source_displacement_work_j,
            external[1][0] * dt * target_v,
            1e-13,
        );
        close(
            drift.source_splitting_work_j,
            kick.source_impulse_work_j - drift.source_displacement_work_j,
            1e-13,
        );
        close(
            drift.integration_defect_j,
            kick.impulse.kinetic_change_j + drift.stored_after_j - drift.stored_before_j
                + drift.dissipated_j
                - drift.source_displacement_work_j,
            2e-11,
        );
        assert!(drift.absolute_integration_defect_j >= drift.integration_defect_j.abs() - 1e-13);

        // Editable receipts are diagnostics, not authority over continuation.
        let mut edited = kick.clone();
        edited.duration_s = 99.0;
        edited.source_impulse_work_j = f64::NAN;
        edited.internal_impulse_work_j = 1e9;
        edited.stored_before_j = f64::NAN;
        edited.impulse.work_residual_j = 1e9;
        let (replayed, replay) = contact
            .motion
            .contact_drift_candidate(&edited, drift_controls(), cx)
            .unwrap();
        assert_eq!(
            replayed.solid().element_states(),
            next.solid().element_states()
        );
        assert_eq!(
            replayed.solid().displacement_m(),
            next.solid().displacement_m()
        );
        assert_eq!(replay.integration_defect_j, drift.integration_defect_j);
        assert!(
            next.contact_drift_candidate(&kick, drift_controls(), cx)
                .is_err()
        );
        let moving_support = contact
            .motion
            .impulse_candidate(&[], &[(0, 1e-3)], origin, impulse_controls(), cx)
            .unwrap();
        assert!(
            moving_support
                .motion
                .contact_drift_candidate(&kick, drift_controls(), cx)
                .is_err()
        );
        assert!(matches!(
            contact.motion.contact_drift_candidate(
                &kick,
                FiniteViscoContactDriftControls {
                    maximum_integration_defect_j: 0.0,
                    ..drift_controls()
                },
                cx
            ),
            Err(fs_solid::hyper3::visco::FiniteViscoTetError::EnergyGate { .. })
        ));
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().displacement_m(), old_u);
        assert_eq!(motion.velocity_mps(), [[0.0; 3]; 4]);
    });
}

#[test]
fn free_graded_material_drift_matches_physical_inventory_without_numerical_kinetic_loss() {
    with_cx(|cx| {
        let law = material(true);
        let velocity = [
            [0.001, -0.002, 0.003],
            [0.02, 0.01, 0.0],
            [-0.01, 0.03, 0.0],
            [0.0, 0.005, -0.01],
        ];
        let motion = FiniteViscoTetMotion::try_new_with_nodal_density(
            relaxed(&law, cx),
            &RHO,
            &velocity,
            cx,
        )
        .unwrap();
        let force = internal(&motion);
        let loads: Vec<_> = force.as_flattened().iter().copied().enumerate().collect();
        let origin = [0.4, -0.3, 0.2];
        let (kicked, kick) = motion
            .force_kick_candidate(0.003, &[], &loads, [0.0; 3], origin, impulse_controls(), cx)
            .unwrap();
        let (next, drift) = kicked
            .contact_drift_candidate(&kick, drift_controls(), cx)
            .unwrap();
        let x0 = std::array::from_fn(|i| {
            std::array::from_fn(|c| NODES[i][c] + kicked.solid().displacement_m()[i][c])
        });
        let x1 = std::array::from_fn(|i| {
            std::array::from_fn(|c| NODES[i][c] + next.solid().displacement_m()[i][c])
        });
        let v = std::array::from_fn(|i| kicked.velocity_mps()[i]);
        let before = inventory(&x0, &v, origin);
        let after = inventory(&x1, &v, origin);
        for axis in 0..3 {
            close(drift.before.linear_kg_m_s[axis], before.0[axis], 1e-12);
            close(drift.after.linear_kg_m_s[axis], after.0[axis], 1e-12);
            close(drift.before.angular_kg_m2_s[axis], before.1[axis], 1e-12);
            close(drift.after.angular_kg_m2_s[axis], after.1[axis], 1e-12);
            close(
                drift.angular_residual_kg_m2_s[axis],
                after.1[axis] - before.1[axis],
                1e-12,
            );
            close(drift.angular_residual_kg_m2_s[axis], 0.0, 1e-12);
            assert_eq!(drift.linear_residual_kg_m_s[axis], 0.0);
        }
        close(drift.kinetic_before_j, before.2, 1e-13);
        close(drift.kinetic_after_j, after.2, 1e-13);
        assert_eq!(drift.kinetic_change_j, 0.0);
        assert_eq!(next.velocity_mps(), kicked.velocity_mps());
        assert!(drift.dissipated_j > 0.0);
        assert_eq!(motion.velocity_mps(), velocity);
    });
}

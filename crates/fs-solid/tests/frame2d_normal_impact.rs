//! G0/G1/G3: coupled frame impact, original-mass work and atomic continuation.
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy, MultiImpactError};
use fs_solid::frame2d::{
    FrameMassFormulation, FrameNormalContact, FrameNormalImpact, FrameNormalImpactError,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const UNDAMPED: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:.3e}"
    );
}
fn vector_close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        close(a, b, tolerance);
    }
}
fn contact(row: &[f64], obstacle_speed: f64, restitution: f64) -> FrameNormalContact<'_> {
    FrameNormalContact {
        normal_covector: row,
        gap: 0.0,
        activation_tolerance: 0.0,
        obstacle_normal_speed: obstacle_speed,
        restitution,
    }
}
fn event(
    frame: &mut Frame2d,
    contacts: &[FrameNormalContact<'_>],
) -> Result<FrameNormalImpact, FrameNormalImpactError> {
    let load = vec![0.0; frame.displacements().len()];
    frame.apply_normal_impacts(
        contacts,
        &load,
        UNDAMPED,
        0.01,
        MultiImpactConfig::default(),
        || false,
    )
}
fn point_mass(mass: f64, velocity: f64) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_mass(0, mass);
    frame.fix(0, [false, true, true]);
    if velocity != 0.0 {
        frame
            .apply_impulse(&[mass * velocity, 0.0, 0.0], &[0.0; 3], UNDAMPED, 0.01)
            .unwrap();
    }
    frame
}
fn elastic() -> SectionModel {
    SectionModel::Elastic {
        ea: 4000.0,
        ei: 20.0,
    }
}
fn axial_pair() -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [false, true, true]);
    frame.fix(1, [false, true, true]);
    frame.add_element(0, 1, 3, 6.0, &elastic).unwrap();
    // Independent consistent axial mass [[2,1],[1,2]].
    // This gives v_before = [-1, 1/4].
    frame
        .apply_impulse(
            &[-1.75, 0.0, 0.0, -0.5, 0.0, 0.0],
            &[0.0; 6],
            UNDAMPED,
            0.01,
        )
        .unwrap();
    frame
}
fn transverse(formulation: FrameMassFormulation) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(2.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [true, false, false]);
    frame.add_element(0, 1, 5, 3.0, &elastic).unwrap();
    frame
}

#[test]
fn g1_moving_obstacle_rebound_separates_input_work_loss_and_right_limit_acceleration() {
    let mut frame = point_mass(2.0, -2.0);
    let row = [1.0, 0.0, 0.0];
    let impact = frame
        .apply_normal_impacts(
            &[contact(&row, 1.0, 0.5)],
            &[3.0, 0.0, 0.0],
            Rayleigh { a0: 0.4, a1: 0.0 },
            0.01,
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
    // Relative approach=-3, effective mass=2, p=2*(1+1/2)*3=9.
    close(impact.contact.normal_impulses[0], 9.0, 1e-13);
    close(frame.velocities()[0], 2.5, 1e-13);
    close(impact.relative_normal_speeds_after[0], 1.5, 1e-13);
    close(frame.accelerations()[0], 0.5, 1e-13);
    close(impact.impulse.kinetic_energy_change_j, 2.25, 1e-13);
    close(impact.obstacle_work_j, 9.0, 1e-13);
    close(impact.relative_work_j, -6.75, 1e-13);
    close(impact.dissipated_energy_j, 6.75, 1e-13);
    close(
        impact.impulse.kinetic_energy_change_j,
        impact.obstacle_work_j + impact.relative_work_j,
        1e-13,
    );
    assert_eq!(frame.displacements(), &[0.0; 3]);
}

#[test]
fn g1_simultaneous_consistent_contacts_prevent_reclosing_an_initially_separating_contact() {
    let mut frame = axial_pair();
    let row0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let row1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let before = format!("{:?}", frame.element_forces(0));
    let result = event(
        &mut frame,
        &[contact(&row0, 0.0, 0.0), contact(&row1, 0.0, 0.0)],
    )
    .unwrap();
    // M^-1=[[2,-1],[-1,2]]/3. Solving both rows gives p=[7/4,1/2].
    // A first-contact-only kick would drive the separating second point closed.
    vector_close(&result.contact.normal_impulses, &[1.75, 0.5], 2e-13);
    vector_close(frame.velocities(), &[0.0; 6], 2e-13);
    close(result.dissipated_energy_j, 13.0 / 16.0, 2e-13);
    close(result.obstacle_work_j, 0.0, 0.0);
    assert_eq!(format!("{:?}", frame.element_forces(0)), before);
    assert!(matches!(
        result.contact.determinacy,
        MultiImpactDeterminacy::UniqueImpulses
    ));
}

#[test]
fn g1_inactive_contact_still_reports_the_coupled_endpoint_without_an_impulse() {
    let mut frame = axial_pair();
    let row0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let row1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let mut distant = contact(&row1, 0.0, 0.0);
    distant.gap = 1.0;
    let result = event(&mut frame, &[contact(&row0, 0.0, 0.0), distant]).unwrap();
    vector_close(&result.contact.normal_impulses, &[1.5, 0.0], 2e-13);
    vector_close(&result.relative_normal_speeds_after, &[0.0, -0.25], 2e-13);
    assert_eq!(result.contact.active_indices, vec![0]);
}

#[test]
fn g1_consistent_tip_impact_retains_support_cross_mass_and_rotational_response() {
    let mut frame = transverse(FrameMassFormulation::ConsistentReference);
    let row = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let result = event(&mut frame, &[contact(&row, 2.0, 0.0)]).unwrap();
    // Tip inertia inverse=[[2,11/2],[11/2,39/2]], so p=1.
    close(result.contact.normal_impulses[0], 1.0, 2e-13);
    vector_close(frame.velocities(), &[0.0, 0.0, 0.0, 0.0, 2.0, 5.5], 5e-13);
    vector_close(
        &result.impulse.support_impulse,
        &[0.0, -0.5, -0.2, 0.0, 0.0, 0.0],
        5e-13,
    );
    close(result.obstacle_work_j, 2.0, 3e-13);
    close(result.impulse.kinetic_energy_change_j, 1.0, 3e-13);
    close(result.dissipated_energy_j, 1.0, 3e-13);
}

#[test]
fn g1_member_contact_uses_original_hermite_trace_and_its_conjugate_moment() {
    let mut frame = transverse(FrameMassFormulation::ConsistentReference);
    let (dofs, weights, length) = frame.reference_transverse_trace(0, 0.5).unwrap();
    close(length, 2.0, 0.0);
    let mut row = [0.0; 6];
    for (dof, weight) in dofs.into_iter().zip(weights) {
        row[dof] += weight;
    }
    let before_u = frame.displacements().to_vec();
    let before_sections = frame.section_deformations(0).to_vec();
    let result = event(&mut frame, &[contact(&row, 1.0, 0.0)]).unwrap();
    // Free trace=[1/2,-1/4]. Its mobility is 11/32, including rotation.
    close(result.contact.normal_impulses[0], 32.0 / 11.0, 2e-12);
    vector_close(
        frame.velocities(),
        &[0.0, 0.0, 0.0, 0.0, -12.0 / 11.0, -68.0 / 11.0],
        2e-12,
    );
    close(result.relative_normal_speeds_after[0], 0.0, 3e-13);
    close(result.impulse.support_impulse[1], 0.0, 2e-12);
    close(result.impulse.support_impulse[2], -4.0 / 55.0, 2e-12);
    close(result.obstacle_work_j, 32.0 / 11.0, 2e-12);
    close(result.dissipated_energy_j, 16.0 / 11.0, 2e-12);
    assert_eq!(frame.displacements(), before_u);
    assert_eq!(frame.section_deformations(0), before_sections);
}

#[test]
fn g3_internal_contact_conserves_total_momentum_and_is_invariant_under_uniform_velocity_boost() {
    let run = |boost: f64| {
        let mut frame = Frame2d::new(Geometry::Linear);
        for (x, mass) in [(0.0, 2.0), (1.0, 3.0)] {
            let i = frame.add_node(x, 0.0);
            frame.add_mass(i, mass);
            frame.fix(i, [false, true, true]);
        }
        let v = [2.0 + boost, -1.0 + boost];
        frame
            .apply_impulse(
                &[2.0 * v[0], 0.0, 0.0, 3.0 * v[1], 0.0, 0.0],
                &[0.0; 6],
                UNDAMPED,
                0.01,
            )
            .unwrap();
        let row = [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let result = event(&mut frame, &[contact(&row, 0.0, 0.5)]).unwrap();
        close(result.contact.normal_impulses[0], 5.4, 3e-13);
        close(frame.velocities()[0], boost - 0.7, 3e-13);
        close(frame.velocities()[3], boost + 0.8, 3e-13);
        close(
            2.0 * frame.velocities()[0] + 3.0 * frame.velocities()[3],
            1.0 + 5.0 * boost,
            3e-13,
        );
        close(result.dissipated_energy_j, 4.05, 4e-13);
        close(result.impulse.kinetic_energy_change_j, -4.05, 4e-13);
        close(result.obstacle_work_j, 0.0, 0.0);
        result
    };
    let rest = run(0.0);
    let boosted = run(7.0);
    vector_close(
        &rest.contact.normal_impulses,
        &boosted.contact.normal_impulses,
        4e-13,
    );
    vector_close(
        &rest.impulse.velocity_change,
        &boosted.impulse.velocity_change,
        4e-13,
    );
}

#[test]
fn g1_impact_continues_with_the_actual_elastic_force_damping_and_acceleration() {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame.add_element(0, 1, 3, 6.0, &elastic).unwrap();
    frame.add_mass(1, 2.0);
    frame
        .static_load(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0], 1)
        .unwrap();
    let before_u = frame.displacements().to_vec();
    let before_sections = frame.section_deformations(0).to_vec();
    let before_force = frame.element_forces(0);
    let row = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let damping = Rayleigh { a0: 0.2, a1: 0.001 };
    let load = [0.0, 0.0, 0.0, 8.0, 0.0, 0.0];
    let h = 0.001;
    let result = frame
        .apply_normal_impacts(
            &[contact(&row, 1.0, 0.0)],
            &load,
            damping,
            h,
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
    let (m, k, c) = (4.0, 4000.0, 4.8);
    let q0 = before_u[3];
    let v0 = 1.0;
    let a0 = (8.0 - k * q0 - c * v0) / m;
    close(frame.velocities()[3], v0, 2e-14);
    close(frame.accelerations()[3], a0, 2e-12);
    assert_eq!(frame.displacements(), before_u);
    assert_eq!(frame.section_deformations(0), before_sections);
    assert_eq!(frame.element_forces(0), before_force);
    let a1 = (8.0 - k * (q0 + h * v0 + 0.25 * h * h * a0) - c * (v0 + 0.5 * h * a0))
        / (m + 0.5 * h * c + 0.25 * h * h * k);
    let q1 = q0 + h * v0 + 0.25 * h * h * (a0 + a1);
    let v1 = v0 + 0.5 * h * (a0 + a1);
    let history = frame
        .newmark_with_load(&[0.0, 0.0], h, (0.0, 0.0), damping, &load)
        .unwrap();
    close(
        history.kinetic_energy[0],
        result.impulse.kinetic_energy_after_j,
        3e-13,
    );
    close(frame.displacements()[3], q1, 3e-15);
    close(frame.velocities()[3], v1, 3e-13);
    close(frame.accelerations()[3], a1, 2e-10);
}

#[test]
fn g3_redundant_contacts_preserve_velocity_only_determinacy_and_refuse_unresolved_source_work() {
    let mut frame = point_mass(2.0, -2.0);
    let row = [1.0, 0.0, 0.0];
    let result = event(&mut frame, &[contact(&row, 0.0, 0.0); 2]).unwrap();
    assert!(matches!(
        result.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    close(result.contact.normal_impulses.iter().sum(), 4.0, 1e-12);
    close(frame.velocities()[0], 0.0, 1e-12);
    close(result.dissipated_energy_j, 4.0, 1e-12);
    let mut unresolved = point_mass(2.0, -2.0);
    let before = format!("{unresolved:?}");
    assert!(matches!(
        event(
            &mut unresolved,
            &[contact(&row, 0.0, 0.0), contact(&row, 1.0, 0.0)]
        ),
        Err(FrameNormalImpactError::Contact(
            MultiImpactError::UnresolvedWork { .. }
        ))
    ));
    assert_eq!(format!("{unresolved:?}"), before);
}

#[test]
fn g0_globally_active_restitution_is_refused_without_a_hidden_inelastic_repair() {
    let mut frame = axial_pair();
    let before = format!("{frame:?}");
    let row0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let row1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let refusal = event(
        &mut frame,
        &[contact(&row0, 0.0, 1.0), contact(&row1, 0.0, 1.0)],
    )
    .unwrap_err();
    match refusal {
        FrameNormalImpactError::Contact(MultiImpactError::EnergyIncrease { work, .. }) => {
            close(work, 3.0 / 16.0, 3e-13)
        }
        other => panic!("wrong refusal: {other:?}"),
    }
    assert_eq!(format!("{frame:?}"), before);
    let retry = event(
        &mut frame,
        &[contact(&row0, 0.0, 0.0), contact(&row1, 0.0, 0.0)],
    )
    .unwrap();
    close(retry.dissipated_energy_j, 13.0 / 16.0, 3e-13);
}

#[test]
fn g0_late_cancellation_preserves_the_complete_frame_and_retry_matches_uninterrupted_event() {
    let original = point_mass(2.0, 0.0); // No damping reference has been captured.
    let row = [1.0, 0.0, 0.0];
    let contacts = [contact(&row, 1.0, 0.0)];
    let mut expected = original.clone();
    let mut polls = 0;
    let expected_event = expected
        .apply_normal_impacts(
            &contacts,
            &[0.0; 3],
            UNDAMPED,
            0.01,
            MultiImpactConfig::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 5);
    for cancel_at in [1, 3, polls] {
        let mut frame = original.clone();
        let before = format!("{frame:?}");
        let mut count = 0;
        let refusal = frame.apply_normal_impacts(
            &contacts,
            &[0.0; 3],
            UNDAMPED,
            0.01,
            MultiImpactConfig::default(),
            || {
                count += 1;
                count == cancel_at
            },
        );
        assert!(matches!(
            refusal,
            Err(FrameNormalImpactError::Contact(MultiImpactError::Cancelled))
        ));
        assert_eq!(format!("{frame:?}"), before);
        let retry = event(&mut frame, &contacts).unwrap();
        assert_eq!(retry, expected_event);
        assert_eq!(format!("{frame:?}"), format!("{expected:?}"));
    }
}

#[test]
fn g0_geometry_limits_invalid_continuation_and_massless_or_locked_motion_refuse_atomically() {
    let row = [1.0, 0.0, 0.0];
    let original = point_mass(2.0, 0.0);
    for case in 0..6 {
        let mut frame = original.clone();
        let before = format!("{frame:?}");
        let mut normal = contact(&row, 1.0, 0.0);
        let mut config = MultiImpactConfig::default();
        let mut dt = 0.01;
        match case {
            0 => normal.gap = -1e-9,
            1 => normal.normal_covector = &[1.0],
            2 => normal.obstacle_normal_speed = f64::NAN,
            3 => normal.restitution = 1.1,
            4 => config.solver.max_variables = 0,
            _ => dt = 0.0,
        }
        assert!(
            frame
                .apply_normal_impacts(&[normal], &[0.0; 3], UNDAMPED, dt, config, || false)
                .is_err()
        );
        assert_eq!(format!("{frame:?}"), before);
    }
    let mut frame = transverse(FrameMassFormulation::Lumped);
    let before = format!("{frame:?}");
    let row = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    assert!(event(&mut frame, &[contact(&row, 1.0, 0.0)]).is_err());
    assert_eq!(format!("{frame:?}"), before);
    let mut locked = original;
    locked.fix(0, [true; 3]);
    let before = format!("{locked:?}");
    assert!(event(&mut locked, &[contact(&[1.0, 0.0, 0.0], 1.0, 0.0)]).is_err());
    assert_eq!(format!("{locked:?}"), before);
}

#[test]
fn g0_a_contact_solution_cannot_publish_an_unrepresentable_stored_velocity_jump() {
    let mut frame = point_mass(1.0, 1e16);
    let before = format!("{frame:?}");
    let row = [1.0, 0.0, 0.0];
    // The law requests delta-v=3. The stored spacing at 1e16 is 2:
    // adding that kick changes the represented velocity by 4, not 3.
    let refusal = event(&mut frame, &[contact(&row, 1e16 + 2.0, 0.5)]);
    assert!(matches!(refusal, Err(FrameNormalImpactError::Frame(_))));
    assert_eq!(format!("{frame:?}"), before);
}

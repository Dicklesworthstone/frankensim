//! G1/G3/G4: original changing-area torsion acts on actual finite bodies.
use super::*;
use fs_mbd::tangential_friction::adhesive::torsion::normal::*;

fn loaded() -> (SourceJkrTorsionContact, Fixture, RigidAdhesiveTorsionState) {
    let source = source(false);
    let dt = 1e-4;
    let mut f = fixture(
        0.3 * source.law().parameters().critical_twist_angle_rad / dt,
        0.0,
        true,
    );
    let mut state = f
        .contact(&source)
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let contact = f.contact(&source);
    apply_rigid_adhesive_torsion(
        &mut f.bodies,
        &f.properties,
        &mut state,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(state.elastic_twist_rad() > 0.0);
    (source, f, state)
}

#[test]
fn g1_changing_area_preserves_preload_and_uses_actual_twisting_inertia() {
    let (source, f, state) = loaded();
    let normal = source.normal_source();
    let pull = normal.force_controlled_pull_off_n().abs();
    let dt = 1e-4;
    let c = RigidAdhesiveTorsionControls {
        energy_tolerance_j: 1e-20,
        ..controls()
    };
    for load in [-0.5 * pull, 0.75 * pull] {
        let target = normal.evaluate_force_controlled(load).unwrap();
        let next = prepare_rigid_adhesive_torsion_with_normal_update(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&source),
            target.indentation_m,
            [Wrench::ZERO; 2],
            dt,
            c,
            &CancelGate::new(),
        )
        .unwrap();
        let ratio = target.contact_radius_m / source.normal_response().contact_radius_m;
        let k = source.law().stiffness_nm_per_rad() * ratio.powi(3);
        let phi0 = state.elastic_twist_rad();
        let rate = next.report.mechanical.twisting_rate_before_rad_s;
        // Solve the independent linear spring/midpoint equation with the
        // original loaded twist and the changed cubic stiffness.
        let expected_angle = (dt * rate - 0.5 * dt * dt * f.mobility() * k * phi0)
            / (1.0 + 0.25 * dt * dt * f.mobility() * k);
        close(
            next.report.mechanical.angle_increment_rad,
            expected_angle,
            1e-18,
        );
        close(
            next.state_after.elastic_twist_rad(),
            phi0 + expected_angle,
            1e-18,
        );
        close(
            next.state_after.accumulated_angle_rad(),
            state.accumulated_angle_rad() + expected_angle,
            1e-18,
        );
        close(
            next.report.torsion_parameter_work_j,
            0.5 * (k - source.law().stiffness_nm_per_rad()) * phi0 * phi0,
            1e-28,
        );
        close(
            next.report.normal_contact_work_j,
            target.reversible_energy_j - source.normal_response().reversible_energy_j,
            1e-20,
        );
        let mechanical = &next.report.mechanical;
        close(
            mechanical.kinetic_energy_change_j,
            -0.5 * k * ((phi0 + expected_angle).powi(2) - phi0.powi(2)),
            1e-25,
        );
        let actual_kinetic_change: f64 = (0..2)
            .map(|i| {
                next.bodies_after[i]
                    .kinetic_energy(f.properties[i])
                    .unwrap()
                    - f.bodies[i].kinetic_energy(f.properties[i]).unwrap()
            })
            .sum();
        close(
            actual_kinetic_change,
            mechanical.kinetic_energy_change_j,
            1e-25,
        );
        vec_close(
            mechanical.angular_impulse_on_a_world_nms,
            vector(f.frame.normal()).scale(-dt * k * (phi0 + 0.5 * expected_angle)),
            1e-25,
        );
        close(
            next.report.represented_storage_change_j
                + mechanical.contact.torsion().irreversible_loss_j,
            next.report.applied_work_j,
            c.energy_tolerance_j,
        );
        assert_eq!(mechanical.contact.torsion().irreversible_loss_j, 0.0);
        assert_eq!(next.source_after.normal_source(), normal);
        assert_eq!(next.source_after.controls(), source.controls());
        assert_eq!(
            next.report.torsion_parameter_work_j.is_sign_positive(),
            load > 0.0
        );
        for i in 0..2 {
            assert_eq!(next.bodies_after[i].pose(), f.bodies[i].pose());
            assert_eq!(
                next.bodies_after[i].linear_momentum_world(),
                f.bodies[i].linear_momentum_world()
            );
        }
    }
}

#[test]
fn g3_unchanged_area_replays_the_original_driven_torsion_phase() {
    let (source, f, state) = loaded();
    let body_axis = f.bodies[0]
        .pose()
        .orientation()
        .rotate_world_to_body(vector(f.frame.normal()));
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-9, 0.0, -2e-9),
            torque_body: body_axis.scale(source.law().critical_moment_n_m()),
        },
        Wrench::ZERO,
    ];
    let ordinary = prepare_rigid_adhesive_torsion_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&source),
        wrenches,
        1e-4,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let continued = prepare_rigid_adhesive_torsion_with_normal_update(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&source),
        source.normal_response().indentation_m,
        wrenches,
        1e-4,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(continued.report.mechanical, ordinary.report);
    assert_eq!(continued.bodies_after, ordinary.bodies_after);
    assert_eq!(continued.state_after, ordinary.state_after);
    assert_eq!(continued.source_after, source);
    assert_eq!(continued.report.normal_contact_work_j, 0.0);
    assert_eq!(continued.report.torsion_parameter_work_j, 0.0);
    assert_ne!(continued.report.mechanical.external_work_j, 0.0);
}

#[test]
fn g4_normal_work_twist_memory_and_native_momenta_publish_together() {
    let (original_source, f, original_state) = loaded();
    let target = original_source
        .normal_source()
        .evaluate_force_controlled(0.001)
        .unwrap()
        .indentation_m;
    let c = RigidAdhesiveTorsionControls {
        energy_tolerance_j: 1e-20,
        ..controls()
    };
    let mut source = original_source.clone();
    let mut state = original_state;
    let mut bodies = f.bodies;
    for failure in 0..4 {
        let gate = CancelGate::new();
        let mut indentation = target;
        let mut wrenches = [Wrench::ZERO; 2];
        let mut frame = f.frame;
        match failure {
            0 => gate.request(),
            1 => indentation = 1.1 * source.normal_source().minimum_attached_indentation_m(),
            2 => wrenches[0].torque_body = Vec3::new(0.0, 0.0, f64::NAN),
            _ => frame = TangentFrame::new(f.frame.normal(), f.frame.lateral()).unwrap(),
        }
        assert!(
            apply_rigid_adhesive_torsion_with_normal_update(
                &mut bodies,
                &f.properties,
                &mut state,
                &mut source,
                f.binding,
                frame,
                indentation,
                wrenches,
                1e-4,
                c,
                &gate,
            )
            .is_err()
        );
        assert_eq!(source, original_source);
        assert_eq!(state, original_state);
        assert_eq!(bodies, f.bodies);
    }
    let candidate = prepare_rigid_adhesive_torsion_with_normal_update(
        &bodies,
        &f.properties,
        &state,
        f.contact(&source),
        target,
        [Wrench::ZERO; 2],
        1e-4,
        c,
        &CancelGate::new(),
    )
    .unwrap();
    let report = apply_rigid_adhesive_torsion_with_normal_update(
        &mut bodies,
        &f.properties,
        &mut state,
        &mut source,
        f.binding,
        f.frame,
        target,
        [Wrench::ZERO; 2],
        1e-4,
        c,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, candidate.report);
    assert_eq!(bodies, candidate.bodies_after);
    assert_eq!(state, candidate.state_after);
    assert_eq!(source, candidate.source_after);
}

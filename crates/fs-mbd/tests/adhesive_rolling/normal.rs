//! G1/G3/G4: real recoil after a load-dependent adhesive area change.
use super::*;
use fs_mbd::tangential_friction::adhesive::rolling::normal::*;

fn loaded() -> (SourceJkrRollingContact, Fixture, RigidAdhesiveRollingState) {
    let source = source(false);
    let dt = 1e-4;
    let mut f = fixture(
        0.3 * source.law().critical_rolling_angle_rad() / dt,
        0.0,
        true,
    );
    let mut state = f
        .contact(&source)
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let contact = f.contact(&source);
    apply_rigid_adhesive_rolling(
        &mut f.bodies,
        &f.properties,
        &mut state,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(state.elastic_shift_m() > 0.0);
    (source, f, state)
}

#[test]
fn g1_changed_area_retains_preload_and_drives_the_actual_inertial_response() {
    let (source, f, state) = loaded();
    let normal = source.normal_source();
    let pull_off = normal.force_controlled_pull_off_n().abs();
    let dt = 1e-4;
    for load in [-0.5 * pull_off, 0.75 * pull_off] {
        let endpoint = normal.evaluate_force_controlled(load).unwrap();
        let c = RigidAdhesiveRollingControls {
            energy_tolerance_j: 1e-20,
            ..controls()
        };
        let next = prepare_rigid_adhesive_rolling_with_normal_update(
            &f.bodies,
            &f.properties,
            &state,
            f.contact(&source),
            endpoint.indentation_m,
            [Wrench::ZERO; 2],
            dt,
            c,
            &CancelGate::new(),
        )
        .unwrap();
        let r = normal.effective_radius_m();
        let k = source.law().stiffness_n_per_m()
            * (endpoint.contact_radius_m / source.normal_response().contact_radius_m).powf(1.5);
        let rotational_k = k * r * r;
        let rate = next.report.mechanical.rolling_rate_before_rad_s;
        let initial_elastic_angle = state.elastic_shift_m() / r;
        let expected_angle = (dt * rate
            - 0.5 * dt * dt * f.mobility() * rotational_k * initial_elastic_angle)
            / (1.0 + 0.25 * dt * dt * f.mobility() * rotational_k);
        close(
            next.report.mechanical.angle_increment_rad,
            expected_angle,
            1e-18,
        );
        close(
            next.state_after.elastic_shift_m(),
            state.elastic_shift_m() + r * expected_angle,
            1e-22,
        );
        close(
            next.state_after.accumulated_angle_rad(),
            state.accumulated_angle_rad() + expected_angle,
            1e-18,
        );
        close(
            next.report.rolling_parameter_work_j,
            0.5 * (k - source.law().stiffness_n_per_m()) * state.elastic_shift_m().powi(2),
            1e-28,
        );
        close(
            next.report.normal_contact_work_j,
            endpoint.reversible_energy_j - source.normal_response().reversible_energy_j,
            1e-20,
        );
        close(
            next.report.represented_storage_change_j
                + next.report.mechanical.contact.rolling().irreversible_loss_j,
            next.report.applied_work_j,
            c.energy_tolerance_j,
        );
        assert!(next.report.absolute_stage_residual_j < c.energy_tolerance_j);
        assert_eq!(next.source_after.normal_source(), normal);
        assert_eq!(next.source_after.controls(), source.controls());
        assert_eq!(
            next.report.mechanical.contact.rolling().irreversible_loss_j,
            0.0
        );
        for i in 0..2 {
            assert_eq!(next.bodies_after[i].pose(), f.bodies[i].pose());
            assert_eq!(
                next.bodies_after[i].linear_momentum_world(),
                f.bodies[i].linear_momentum_world(),
            );
        }
        assert_eq!(
            next.report.rolling_parameter_work_j.is_sign_positive(),
            load > 0.0
        );
    }
}

#[test]
fn g3_unchanged_normal_area_replays_the_original_native_phase_exactly() {
    let (source, f, state) = loaded();
    let ordinary = prepare_rigid_adhesive_rolling(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&source),
        1e-4,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let continued = prepare_rigid_adhesive_rolling_with_normal_update(
        &f.bodies,
        &f.properties,
        &state,
        f.contact(&source),
        source.normal_response().indentation_m,
        [Wrench::ZERO; 2],
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
    assert_eq!(continued.report.rolling_parameter_work_j, 0.0);
}

#[test]
fn g4_normal_update_and_recoil_publish_all_owners_only_after_acceptance() {
    let (original_source, f, original_state) = loaded();
    let target = original_source
        .normal_source()
        .evaluate_force_controlled(0.001)
        .unwrap()
        .indentation_m;
    let c = RigidAdhesiveRollingControls {
        energy_tolerance_j: 1e-20,
        ..controls()
    };
    let mut source = original_source.clone();
    let mut state = original_state;
    let mut bodies = f.bodies;
    for failure in 0..4 {
        let gate = CancelGate::new();
        let mut wrenches = [Wrench::ZERO; 2];
        let mut indentation = target;
        let mut frame = f.frame;
        match failure {
            0 => gate.request(),
            1 => {
                indentation = original_source
                    .normal_source()
                    .minimum_attached_indentation_m()
                    - 1e-8;
            }
            2 => wrenches[0].torque_body = Vec3::new(f64::NAN, 0.0, 0.0),
            _ => frame = TangentFrame::new(f.frame.normal(), f.frame.lateral()).unwrap(),
        }
        assert!(
            apply_rigid_adhesive_rolling_with_normal_update(
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
    let expected = prepare_rigid_adhesive_rolling_with_normal_update(
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
    let report = apply_rigid_adhesive_rolling_with_normal_update(
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
    assert_eq!(report, expected.report);
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
    assert_eq!(source, expected.source_after);
}

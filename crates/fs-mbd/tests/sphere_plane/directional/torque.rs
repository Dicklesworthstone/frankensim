//! G1/G3/G4: real motor work through anisotropic support, impact and source replay.

use super::*;
use fs_mbd::sphere_plane::directional::impact::advance_sphere_plane_directional_impacts_with_torque;
use fs_tribo::tangential::elliptic::impulse::EllipticImpulseControls;
use std::cell::Cell;

fn torque_balances(
    scene: &SpherePlane,
    initial: SpherePlaneState,
    duration: f64,
    torque: Vec3,
    result: &SpherePlaneDirectionalAdvance,
) {
    let motion = &result.motion;
    close(
        kinetic(motion.state_after.body()) - kinetic(initial.body()),
        motion.external_force_work + motion.external_torque_work + motion.plane_work
            - motion.dissipated_energy,
    );
    close_vector(motion.external_torque_impulse_world, torque.scale(duration));
    close_vector(
        motion
            .state_after
            .body()
            .linear_momentum_world()
            .sub(initial.body().linear_momentum_world()),
        scene
            .force_world
            .scale(duration)
            .add(motion.contact_impulse_world),
    );
    close(motion.energy_balance_residual, 0.0);
    close_vector(motion.linear_momentum_residual, Vec3::ZERO);
    close_vector(motion.angular_momentum_residual, Vec3::ZERO);
    close(
        motion.dissipated_energy,
        result.sustained_friction_heat_j + result.impact_loss_j,
    );
}

fn impulse_controls() -> EllipticImpulseControls {
    EllipticImpulseControls {
        maximum_iterations: 128,
        maximum_impulse_ns: 1e4,
        impulse_absolute_tolerance_ns: 1e-12,
        impulse_relative_tolerance: 1e-13,
        velocity_tolerance_m_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}

#[test]
fn the_same_motor_breaks_weak_texture_but_drives_lossless_rolling_on_strong_texture() {
    let scene = scene();
    let support = EllipticCoulomb::new([0.2, 0.4], [0.1, 0.3], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.01, 0.02], [0.01, 0.02], [1.0, 0.0]).unwrap();
    let initial = supported(Vec3::ZERO, Vec3::ZERO);
    let mut history = support.initial_state();
    let weak_torque = Vec3::new(0.0, 0.0, -4.0);
    let weak = advance_sphere_plane_directional_impacts_with_torque(
        &scene,
        initial,
        weak_torque,
        &support,
        &mut history,
        &impact,
        &impact.initial_state(),
        frame(),
        0.2,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    // Static demand is 40/7 N, above the weak 4 N limit. Kinetic traction
    // is 2 N: a_x=1, alpha_z=-15 and contact-slip acceleration is -6.5.
    assert_eq!(weak.motion.state_after.mode(), SpherePlaneMode::Sliding);
    close_vector(
        velocity(weak.motion.state_after.body()),
        Vec3::new(0.2, 0.0, 0.0),
    );
    close_vector(
        spin(weak.motion.state_after.body()),
        Vec3::new(0.0, 0.0, -3.0),
    );
    close_vector(
        weak.motion.state_after.body().pose().position_world(),
        Vec3::new(0.02, RADIUS, 0.0),
    );
    close(weak.motion.external_torque_work, 1.2);
    close(weak.sustained_friction_heat_j, 0.26);
    close(weak.support_distance_m, 0.13);
    close(history.slip_m()[0], -0.13);
    close(history.slip_m()[1], 0.0);
    torque_balances(&scene, initial, 0.2, weak_torque, &weak);

    let mut history = support.initial_state();
    let strong_torque = Vec3::new(4.0, 0.0, 0.0);
    let strong = advance_sphere_plane_directional_impacts_with_torque(
        &scene,
        initial,
        strong_torque,
        &support,
        &mut history,
        &impact,
        &impact.initial_state(),
        frame(),
        0.2,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    // The lateral static limit is 8 N. Traction 40/7 N produces
    // a_z=20/7 and alpha_x=40/7, so the actual contact point remains at rest.
    assert_eq!(strong.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(
        velocity(strong.motion.state_after.body()),
        Vec3::new(0.0, 0.0, 4.0 / 7.0),
    );
    close_vector(
        spin(strong.motion.state_after.body()),
        Vec3::new(8.0 / 7.0, 0.0, 0.0),
    );
    close_vector(
        strong.motion.state_after.body().pose().position_world(),
        Vec3::new(0.0, RADIUS, 2.0 / 35.0),
    );
    close(strong.motion.external_torque_work, 16.0 / 35.0);
    close(strong.sustained_friction_heat_j, 0.0);
    close(strong.support_distance_m, 0.0);
    assert_eq!(history, support.initial_state());
    torque_balances(&scene, initial, 0.2, strong_torque, &strong);
}

#[test]
fn cancelling_contact_drive_keeps_distinct_force_and_couple_work() {
    let mut scene = scene();
    scene.force_world.x = 4.0;
    let support = EllipticCoulomb::new([0.2, 0.4], [0.1, 0.3], [1.0, 0.0]).unwrap();
    let initial = supported(Vec3::ZERO, Vec3::ZERO);
    let torque = Vec3::new(0.0, 0.0, -0.8);
    let mut history = support.initial_state();
    let result = advance_sphere_plane_directional_with_torque(
        &scene,
        initial,
        torque,
        &mut history,
        frame(),
        0.2,
        controls(),
        |history, phase, visit, cancel| {
            // F/m and (T/I) cross r cancel in slip space. They still
            // accelerate both physical momenta and perform independent work.
            close(phase.input.driving_acceleration_m_s2[0], 0.0);
            close(phase.input.driving_acceleration_m_s2[1], 0.0);
            let report = support.advance_driven_segments(
                history,
                &phase.input,
                controls().driven,
                &mut || cancel(),
                &mut |segment| visit(segment),
            )?;
            let next = report.state.clone();
            Ok::<_, EllipticDrivenSegmentError<SpherePlaneError>>((report, next))
        },
        |_, _, _| Ok(None),
        || false,
    )
    .unwrap();
    assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close_vector(
        velocity(result.motion.state_after.body()),
        Vec3::new(0.4, 0.0, 0.0),
    );
    close_vector(
        spin(result.motion.state_after.body()),
        Vec3::new(0.0, 0.0, -0.8),
    );
    close_vector(
        result.motion.contact_impulse_world,
        Vec3::new(0.0, 4.0, 0.0),
    );
    close(result.motion.external_force_work, 0.16);
    close(result.motion.external_torque_work, 0.064);
    close(result.sustained_friction_heat_j, 0.0);
    close(result.support_distance_m, 0.0);
    assert_eq!(history, support.initial_state());
    torque_balances(&scene, initial, 0.2, torque, &result);
}

#[test]
fn flight_motor_changes_directional_impact_and_preserves_world_and_plane_work() {
    let scene = scene();
    let support = EllipticCoulomb::new([0.2, 0.4], [0.1, 0.3], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.01, 0.02], [0.01, 0.02], [1.0, 0.0]).unwrap();
    let initial = state(
        Vec3::new(0.0, RADIUS + 0.2, 0.0),
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let torque = Vec3::new(0.0, 0.0, 0.4);
    let mut history = support.initial_state();
    let impact_history = impact.initial_state();
    let result = advance_sphere_plane_directional_impacts_with_torque(
        &scene,
        initial,
        torque,
        &support,
        &mut history,
        &impact,
        &impact_history,
        frame(),
        0.3,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    let event = result
        .motion
        .events
        .iter()
        .find(|event| event.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(event.time_s, 0.2);
    close(event.normal_impulse, 6.0);
    // Flight motor spin makes pre-impact slip .2 m/s. The weak impulse
    // ellipse caps J_x at -.06 N s, leaving slip .095 m/s and omega_z=.25.
    close_vector(event.tangent_impulse_world, Vec3::new(-0.06, 0.0, 0.0));
    close_vector(
        velocity(result.motion.state_after.body()),
        Vec3::new(-0.03, 0.0, 0.0),
    );
    close_vector(
        spin(result.motion.state_after.body()),
        Vec3::new(0.0, 0.0, 0.45),
    );
    close_vector(
        result.motion.state_after.body().pose().position_world(),
        Vec3::new(-0.003, 0.55, 0.0),
    );
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.075).unwrap();
    close_vector(
        result
            .motion
            .state_after
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
    );
    close(result.motion.external_torque_work, 0.03);
    close(result.impact_loss_j, 3.0 + 0.06 * (0.2 + 0.095) / 2.0);
    close(result.sustained_friction_heat_j, 0.0);
    assert_eq!(result.support_phases, 0);
    assert_eq!(history, support.initial_state());
    assert_eq!(impact_history, impact.initial_state());
    torque_balances(&scene, initial, 0.3, torque, &result);

    // Transform the complete physical problem, including the world torque
    // and plane texture. Common translation changes only the prescribed work.
    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let boost = Vec3::new(0.3, -0.2, 0.4);
    let transformed_scene = SpherePlane {
        normal_world: rotate(scene.normal_world),
        force_world: rotate(scene.force_world),
        velocity_world: boost,
        ..scene
    };
    let transformed_initial = state(
        rotate(initial.body().pose().position_world()),
        boost,
        Vec3::ZERO,
        rotation,
    );
    let n = transformed_scene.normal_world;
    let t = rotate(Vec3::new(1.0, 0.0, 0.0));
    let transformed_frame = TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap();
    let mut transformed_history = support.initial_state();
    let transformed = advance_sphere_plane_directional_impacts_with_torque(
        &transformed_scene,
        transformed_initial,
        rotate(torque),
        &support,
        &mut transformed_history,
        &impact,
        &impact_history,
        transformed_frame,
        0.3,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    close_vector(
        transformed
            .motion
            .state_after
            .body()
            .pose()
            .position_world(),
        rotate(result.motion.state_after.body().pose().position_world()).add(boost.scale(0.3)),
    );
    close_vector(
        velocity(transformed.motion.state_after.body()),
        rotate(velocity(result.motion.state_after.body())).add(boost),
    );
    close_vector(
        spin(transformed.motion.state_after.body()),
        rotate(spin(result.motion.state_after.body())),
    );
    close_vector(
        transformed.motion.contact_impulse_world,
        rotate(result.motion.contact_impulse_world),
    );
    close(
        transformed.motion.external_torque_work,
        result.motion.external_torque_work,
    );
    close(
        transformed.motion.plane_work,
        transformed.motion.contact_impulse_world.dot(boost),
    );
    close(transformed.impact_loss_j, result.impact_loss_j);
    assert_eq!(transformed_history, history);
    torque_balances(
        &transformed_scene,
        transformed_initial,
        0.3,
        rotate(torque),
        &transformed,
    );
}

#[test]
fn driven_source_failure_keeps_history_private_until_body_and_work_accept() {
    let scene = scene();
    let support = EllipticCoulomb::new([0.2, 0.4], [0.1, 0.3], [1.0, 0.0]).unwrap();
    let initial = supported(Vec3::ZERO, Vec3::ZERO);
    let original = (support.initial_state(), 17usize);
    let mut control = controls();
    control.maximum_support_phase_s = 0.02;
    for failure in [1, 2, 3, 4, 0] {
        let torque = Vec3::new(0.0, 0.0, if failure == 4 { f64::NAN } else { -4.0 });
        let mut history = original.clone();
        let visits = Cell::new(0usize);
        let cancel = Cell::new(false);
        let result = advance_sphere_plane_directional_with_torque(
            &scene,
            initial,
            torque,
            &mut history,
            frame(),
            0.1,
            control,
            |history, phase, visit, poll| {
                close(phase.input.driving_acceleration_m_s2[0], -10.0);
                close(phase.input.driving_acceleration_m_s2[1], 0.0);
                if phase.time_s >= 0.02 {
                    if failure == 1 {
                        return Err("late source refusal");
                    }
                    if failure == 2 {
                        cancel.set(true);
                    }
                }
                let mut report = support
                    .advance_driven_segments(
                        &history.0,
                        &phase.input,
                        control.driven,
                        &mut || poll(),
                        &mut |segment| {
                            visits.set(visits.get() + 1);
                            visit(segment)
                        },
                    )
                    .map_err(|_| "native replay refused")?;
                if failure == 3 {
                    report.dissipated_j += 1.0;
                }
                let next = (report.state.clone(), history.1 + 1);
                Ok((report, next))
            },
            |_, _, _| Ok(None),
            || cancel.get(),
        );
        if failure != 0 {
            assert!(result.is_err(), "failure {failure} unexpectedly published");
            assert_eq!(history, original);
            if failure == 4 {
                assert_eq!(visits.get(), 0);
            } else {
                assert!(visits.get() > 0, "must stage actual torque-driven motion");
            }
        } else {
            let result = result.unwrap();
            assert_eq!(history.1, 17 + result.support_phases);
            close(history.0.slip_m()[0], -0.0325);
            close(result.sustained_friction_heat_j, 0.065);
            close(result.motion.external_torque_work, 0.3);
            torque_balances(&scene, initial, 0.1, torque, &result);
        }
    }
    close_vector(
        initial.body().pose().position_world(),
        Vec3::new(0.0, RADIUS, 0.0),
    );
    close_vector(velocity(initial.body()), Vec3::ZERO);
    close_vector(spin(initial.body()), Vec3::ZERO);
}

//! G1/G4: motor work and source provenance on the same accepted contact trajectory.
use super::*;

fn motor_balance(
    scene: &SpherePlane,
    before: SpherePlaneState,
    report: &SpherePlaneDirectionalAdvance,
) {
    let motion = &report.motion;
    close(
        motion
            .state_after
            .body()
            .kinetic_energy(scene.properties)
            .unwrap()
            - before.body().kinetic_energy(scene.properties).unwrap()
            + report.sustained_friction_heat_j
            + report.impact_loss_j,
        motion.external_force_work + motion.external_torque_work + motion.plane_work,
    );
    close(motion.energy_balance_residual, 0.0);
    for residual in [
        motion.linear_momentum_residual,
        motion.angular_momentum_residual,
    ] {
        for x in [residual.x, residual.y, residual.z] {
            close(x, 0.0);
        }
    }
}

#[test]
fn motor_rolling_uses_original_static_capacity_without_advancing_material_path() {
    let (card, law, query) = source_fixture(scalar(0.2));
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = resting(0.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let dt = 0.125;
    let report = advance_source_sphere_plane_with_torque(
        &scene,
        before,
        Vec3::new(0.0, 7.0, 0.0),
        &law,
        &mut history,
        frame(),
        dt,
        controls(),
        |_, phase| {
            close(phase.input.driving_acceleration_m_s2[0], -17.5);
            Ok::<_, &'static str>(support_conditions(&query, phase))
        },
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        || false,
    )
    .unwrap();
    let motion = &report.trajectory;
    // Rolling needs 10 N, above the 4 N kinetic limit but below 18 N static.
    assert_eq!(motion.motion.state_after.mode(), SpherePlaneMode::Rolling);
    assert!(report.impacts.is_empty());
    close(motion.motion.contact_impulse_world.x, 10.0 * dt);
    close(
        motion.motion.state_after.body().pose().position_world().x,
        2.5 * dt * dt,
    );
    close(
        physical_slip(motion.motion.state_after.body(), scene.properties),
        0.0,
    );
    close(motion.motion.external_torque_work, 35.0 * dt * dt);
    close(motion.motion.external_torque_impulse_world.y, 7.0 * dt);
    close(motion.support_distance_m, 0.0);
    close(motion.sustained_friction_heat_j, 0.0);
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25);
        close(source.state_coordinate("speed").unwrap(), 0.0);
    }
    original_claims(&card, &original, &history);
    motor_balance(&scene, before, motion);
}

#[test]
fn motor_breakaway_closes_original_sources_against_actual_slip_and_heat() {
    let (card, law, query) = source_fixture(scalar(0.2));
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = resting(0.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let dt = 0.125;
    let report = advance_source_sphere_plane_with_torque(
        &scene,
        before,
        Vec3::new(0.0, -14.0, 0.0),
        &law,
        &mut history,
        frame(),
        dt,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        || false,
    )
    .unwrap();
    let motion = &report.trajectory;
    // The 20 N sticking demand exceeds 18 N. With 4 N kinetic traction,
    // actual slip accelerates at 35 - 1.75*4 = 28 m/s².
    let distance = 14.0 * dt * dt;
    assert_eq!(motion.motion.state_after.mode(), SpherePlaneMode::Sliding);
    assert!(report.support_source_evaluations > 0);
    assert!(report.impacts.is_empty());
    close(
        physical_slip(motion.motion.state_after.body(), scene.properties),
        28.0 * dt,
    );
    close(motion.support_distance_m, distance);
    close(motion.sustained_friction_heat_j, 4.0 * distance);
    close(motion.motion.external_torque_work, 420.0 * dt * dt);
    close(motion.motion.external_torque_impulse_world.y, -14.0 * dt);
    close(history.constitutive_state().slip_m()[0], distance);
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25 + distance);
    }
    original_claims(&card, &original, &history);
    motor_balance(&scene, before, motion);
}

#[test]
fn flight_motor_spin_changes_the_sourced_impact_coefficient_and_bounce_work() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.7)]));
    let mut scene = sphere(Vec3::ZERO);
    scene.restitution = 0.5;
    let before = incoming(1.0, 0.0, -1.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let report = advance_source_sphere_plane_with_torque(
        &scene,
        before,
        Vec3::new(0.0, -0.8, 0.0),
        &law,
        &mut history,
        frame(),
        0.625,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |accepted, phase| {
            assert_eq!(accepted, &original);
            close(phase.time_s, 0.5);
            close(phase.input.velocity_before_m_s[0], 1.0);
            Ok::<_, &'static str>(conditions(&query))
        },
        || false,
    )
    .unwrap();
    let motion = &report.trajectory;
    assert_eq!(report.impacts.len(), 1);
    assert_eq!(motion.support_phases, 0);
    assert_eq!(report.support_source_evaluations, 0);
    let impact = &report.impacts[0].response;
    // Pure motor spin creates 1 m/s preimpact slip; source mu=.1+.05*1=.15.
    close(impact.queried_speed_m_s, 1.0);
    close(impact.coefficients_mu[0], 0.15);
    close(impact.impulse.impulse_ns[0], -0.45);
    close(impact.impulse.dissipated_j, 0.2728125);
    close(motion.impact_loss_j, 1.0228125);
    close(motion.motion.external_torque_work, 0.5125);
    close(motion.motion.external_torque_impulse_world.y, -0.5);
    let after = motion.motion.state_after.body();
    close(after.pose().position_world().x, -0.028125);
    close(after.pose().position_world().z, 0.5625);
    close(
        after
            .point_kinematics(scene.properties, Vec3::ZERO)
            .unwrap()
            .angular_velocity_world
            .y,
        -1.375,
    );
    close(motion.sustained_friction_heat_j, 0.0);
    close(motion.support_distance_m, 0.0);
    assert_eq!(history.constitutive_state().slip_m(), [0.0; 2]);
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25);
        close(source.state_coordinate("speed").unwrap(), 1.0);
    }
    original_claims(&card, &original, &impact.state);
    original_claims(&card, &original, &history);
    motor_balance(&scene, before, motion);
}

#[test]
fn failed_source_after_private_motor_motion_and_impact_preserves_all_history() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = incoming(0.5, 1.0, -1.0);
    let original = law.initial_state();
    for cancel_late in [false, true] {
        let mut history = original.clone();
        let impacts = Cell::new(0);
        let supports = Cell::new(0);
        let stop = Cell::new(false);
        let result = advance_source_sphere_plane_with_torque(
            &scene,
            before,
            Vec3::new(0.0, -14.0, 0.0),
            &law,
            &mut history,
            frame(),
            0.125,
            controls(),
            |accepted, phase| {
                supports.set(supports.get() + 1);
                assert_eq!(impacts.get(), 1);
                if supports.get() == 3 {
                    assert!(phase.time_s > 0.0);
                    assert!(accepted.constitutive_state().slip_m()[0] > 0.0);
                    if cancel_late {
                        stop.set(true);
                    } else {
                        return Err("motor-path source unavailable");
                    }
                }
                Ok::<_, &'static str>(support_conditions(&query, phase))
            },
            |_, _| {
                impacts.set(impacts.get() + 1);
                Ok::<_, &'static str>(conditions(&query))
            },
            || stop.get(),
        );
        assert!(result.is_err());
        if !cancel_late {
            assert!(matches!(
                result,
                Err(SpherePlaneDirectionalError::Source(
                    SourceSpherePlaneError::SupportConditions("motor-path source unavailable")
                ))
            ));
        }
        assert_eq!(supports.get(), 3);
        assert_eq!(impacts.get(), 1);
        assert_eq!(history, original);
        assert_eq!(before.body(), incoming(0.5, 1.0, -1.0).body());
    }
}

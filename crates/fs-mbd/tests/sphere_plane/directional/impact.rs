//! G1/G3/G4: instantaneous directional friction through actual contact events.
use super::*;
use fs_mbd::sphere_plane::directional::impact::*;
use fs_tribo::tangential::elliptic::impulse::EllipticImpulseControls;
use std::cell::Cell;

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
fn principal_impact_has_exact_spin_loss_pose_and_isotropic_recovery() {
    let support = EllipticCoulomb::new([0.4; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.05, 0.2], [0.05, 0.2], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.force_world = Vec3::ZERO;
    scene.friction_coefficient = 9.0;
    let input = supported(Vec3::new(1.0, -1.0, 0.0), Vec3::ZERO);
    let mut history = support.initial_state();
    let original = history.clone();
    let impact_history = impact.initial_state();
    let impact_original = impact_history.clone();
    let result = advance_sphere_plane_directional_impacts(
        &scene,
        input,
        &support,
        &mut history,
        &impact,
        &impact_history,
        frame(),
        0.2,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    close_vector(
        result.motion.contact_impulse_world,
        Vec3::new(-0.15, 3.0, 0.0),
    );
    close_vector(
        velocity(result.motion.state_after.body()),
        Vec3::new(0.925, 0.5, 0.0),
    );
    close_vector(
        spin(result.motion.state_after.body()),
        Vec3::new(0.0, 0.0, -0.375),
    );
    close_vector(
        result.motion.state_after.body().pose().position_world(),
        Vec3::new(0.185, 0.6, 0.0),
    );
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), -0.075).unwrap();
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
    close(result.impact_loss_j, 0.75 + 0.15 * (1.0 + 0.7375) / 2.0);
    assert_eq!(result.sustained_friction_heat_j, 0.0);
    assert_eq!(result.support_distance_m, 0.0);
    assert_eq!(history, original);
    assert_eq!(impact_history, impact_original);
    assert_eq!(result.motion.events.len(), 1);
    assert_eq!(result.motion.events[0].kind, SpherePlaneEventKind::Impact);
    close(result.motion.events[0].time_s, 0.0);
    balances(&scene, input, 0.2, &result.motion);

    scene.friction_coefficient = 0.05;
    let isotropic = EllipticCoulomb::new([0.05; 2], [0.05; 2], [0.6, 0.8]).unwrap();
    let isotropic_history = isotropic.initial_state();
    let mut native_history = support.initial_state();
    let native = advance_sphere_plane_directional_impacts(
        &scene,
        input,
        &support,
        &mut native_history,
        &isotropic,
        &isotropic_history,
        frame(),
        0.2,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    let mut legacy_history = support.initial_state();
    let legacy = advance_sphere_plane_directional(
        &scene,
        input,
        &support,
        &mut legacy_history,
        frame(),
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    close_vector(
        native.motion.state_after.body().pose().position_world(),
        legacy.motion.state_after.body().pose().position_world(),
    );
    close_vector(
        velocity(native.motion.state_after.body()),
        velocity(legacy.motion.state_after.body()),
    );
    close_vector(
        spin(native.motion.state_after.body()),
        spin(legacy.motion.state_after.body()),
    );
    close(native.impact_loss_j, legacy.impact_loss_j);
    assert_eq!(native_history, legacy_history);
}

#[test]
fn detected_directional_impact_flows_into_sustained_arrest_without_invented_distance() {
    let support = EllipticCoulomb::new([0.3, 0.4], [0.2, 0.3], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.01, 0.02], [0.01, 0.02], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.restitution = 0.0;
    let input = state(
        Vec3::new(0.0, 1.75, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
    );
    let mut history = support.initial_state();
    let impact_history = impact.initial_state();
    let result = advance_sphere_plane_directional_impacts(
        &scene,
        input,
        &support,
        &mut history,
        &impact,
        &impact_history,
        frame(),
        0.75,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    let impact_event = result
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(impact_event.time_s, 0.5);
    close(impact_event.normal_impulse, 10.0);
    close_vector(
        impact_event.tangent_impulse_world,
        Vec3::new(-0.1, 0.0, 0.0),
    );
    let post_impact_slip: f64 = 0.825;
    let arrest = result
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Rolling)
        .unwrap();
    close(arrest.time_s, 0.5 + post_impact_slip / 7.0);
    close(
        result.impact_loss_j,
        25.0 + 0.1 * (1.0 + post_impact_slip) / 2.0,
    );
    close(
        result.sustained_friction_heat_j,
        post_impact_slip.powi(2) / 3.5,
    );
    close(result.support_distance_m, post_impact_slip.powi(2) / 14.0);
    close(history.slip_m()[0], result.support_distance_m);
    close(history.slip_m()[1], 0.0);
    assert_eq!(impact_history.slip_m(), [0.0; 2]);
    assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close(velocity(result.motion.state_after.body()).x, 5.0 / 7.0);
    balances(&scene, input, 0.75, &result.motion);
}

#[test]
fn off_axis_impact_preserves_rotated_texture_and_prescribed_plane_work() {
    let support = EllipticCoulomb::new([0.4; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.1, 0.4], [0.1, 0.4], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.force_world = Vec3::ZERO;
    let input = supported(Vec3::new(2.0, -1.0, -1.0), Vec3::ZERO);
    let mut history = support.initial_state();
    let baseline = advance_sphere_plane_directional_impacts(
        &scene,
        input,
        &support,
        &mut history,
        &impact,
        &impact.initial_state(),
        frame(),
        0.125,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    let impulse = baseline.motion.events[0].tangent_impulse_world;
    assert!((impulse.x * -1.0 - impulse.z * 2.0).abs() > 0.1);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), 0.71).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let boost = Vec3::new(0.3, -0.2, 0.4);
    scene.normal_world = rotate(scene.normal_world);
    scene.velocity_world = boost;
    // Rotating exact contact can round inside the strictly nonpenetrating plane.
    // Give the original event detector a positive, roundoff-scale flight gap.
    let clearance = 64.0 * f64::EPSILON * scene.radius;
    let transformed = state(
        rotate(input.body().pose().position_world()).add(scene.normal_world.scale(clearance)),
        rotate(velocity(input.body())).add(boost),
        Vec3::ZERO,
        rotation,
    );
    let n = scene.normal_world;
    let t = rotate(Vec3::new(1.0, 0.0, 0.0));
    let transformed_frame = TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap();
    let mut transformed_history = support.initial_state();
    let result = advance_sphere_plane_directional_impacts(
        &scene,
        transformed,
        &support,
        &mut transformed_history,
        &impact,
        &impact.initial_state(),
        transformed_frame,
        0.125,
        controls(),
        impulse_controls(),
        || false,
    )
    .unwrap();
    close_vector(
        result.motion.state_after.body().pose().position_world(),
        rotate(baseline.motion.state_after.body().pose().position_world()).add(boost.scale(0.125)),
    );
    close_vector(
        velocity(result.motion.state_after.body()),
        rotate(velocity(baseline.motion.state_after.body())).add(boost),
    );
    close_vector(
        spin(result.motion.state_after.body()),
        rotate(spin(baseline.motion.state_after.body())),
    );
    close_vector(
        result.motion.contact_impulse_world,
        rotate(baseline.motion.contact_impulse_world),
    );
    close(result.impact_loss_j, baseline.impact_loss_j);
    close(
        result.motion.plane_work,
        result.motion.contact_impulse_world.dot(boost),
    );
    assert_eq!(transformed_history, history);
    balances(&scene, transformed, 0.125, &result.motion);
}

#[derive(Clone, Debug, PartialEq)]
struct CoupledHistory {
    support: EllipticCoulombState,
    impacts: usize,
}

#[test]
fn late_callback_cancel_and_forged_impact_work_roll_back_then_retry_exactly() {
    let support = EllipticCoulomb::new([0.4; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let impact = EllipticCoulomb::new([0.05, 0.2], [0.05, 0.2], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.restitution = 0.0;
    let input = supported(Vec3::new(1.0, -1.0, 0.0), Vec3::ZERO);
    let original = CoupledHistory {
        support: support.initial_state(),
        impacts: 0,
    };
    let mut successful = None;
    for failure in [1, 2, 3, 0, 0] {
        let mut history = original.clone();
        let phases = Cell::new(0usize);
        let visited = Cell::new(0usize);
        let cancel_now = Cell::new(false);
        let result = advance_sphere_plane_directional_with_impacts(
            &scene,
            input,
            &mut history,
            frame(),
            0.2,
            controls(),
            |history, phase, visit, cancel| {
                assert_eq!(history.impacts, 1);
                phases.set(phases.get() + 1);
                if failure == 1 && phases.get() == 2 {
                    return Err("late source");
                }
                let report = support
                    .advance_driven_segments(
                        &history.support,
                        &phase.input,
                        controls().driven,
                        &mut || cancel(),
                        &mut |segment| {
                            visit(segment)?;
                            visited.set(visited.get() + 1);
                            if failure == 2 && visited.get() == 2 {
                                cancel_now.set(true);
                            }
                            Ok::<(), SpherePlaneError>(())
                        },
                    )
                    .map_err(|_| "support")?;
                let mut next = history.clone();
                next.support = report.state.clone();
                Ok((report, next))
            },
            |history, phase, cancel| {
                assert_eq!(history.impacts, 0);
                close(phase.time_s, 0.0);
                close(phase.normal_velocity_before_m_s, -1.0);
                close(phase.normal_velocity_after_m_s, 0.0);
                let mut report = impact
                    .solve_impulse(
                        &impact.initial_state(),
                        &phase.input,
                        impulse_controls(),
                        &mut || cancel(),
                    )
                    .map_err(|_| "impact")?;
                if failure == 3 {
                    report.dissipated_j += 1.0;
                }
                let mut next = history.clone();
                next.impacts += 1;
                Ok((report, next))
            },
            || cancel_now.get(),
        );
        if failure != 0 {
            assert!(result.is_err());
            assert_eq!(history, original);
            if failure != 3 {
                assert!(visited.get() > 0);
            }
        } else {
            let report = result.unwrap();
            assert_eq!(history.impacts, 1);
            assert!(history.support.slip_m()[0] > 0.0);
            balances(&scene, input, 0.2, &report.motion);
            let accepted = (
                report.motion.state_after,
                history,
                report.impact_loss_j,
                report.sustained_friction_heat_j,
            );
            if let Some(previous) = &successful {
                assert_eq!(&accepted, previous);
            } else {
                successful = Some(accepted);
            }
        }
    }
    close_vector(
        input.body().pose().position_world(),
        Vec3::new(0.0, RADIUS, 0.0),
    );
}

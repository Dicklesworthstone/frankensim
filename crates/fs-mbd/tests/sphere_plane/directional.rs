//! G1/G3/G4: original moving geometry, anisotropic support and accepted history.
use super::*;
use fs_mbd::sphere_plane::advance_sphere_plane_with_support_friction;
use fs_mbd::sphere_plane::directional::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::driven::*;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;
use fs_tribo::tangential::elliptic::{EllipticCoulomb, EllipticCoulombState};
use fs_tribo::{FrictionLaw, TangentialLaw, TangentialSolveConfig};

fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap()
}
fn controls() -> SpherePlaneDirectionalControls {
    SpherePlaneDirectionalControls {
        motion: SpherePlaneConfig {
            max_attitude_step_s: 0.0005,
            ..SpherePlaneConfig::default()
        },
        maximum_support_phase_s: 0.03,
        maximum_support_phases: 1000,
        driven: EllipticDrivenControls {
            maximum_step_s: 0.002,
            minimum_step_s: 1e-10,
            maximum_substeps: 100_000,
            rest_velocity_tolerance_m_s: 1e-10,
            solve: TangentialSolveConfig {
                max_iterations: 128,
                maximum_force_n: 1e5,
                force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-13,
            },
            relaxation: EllipticRelaxationControls {
                time_tolerance_s: 1e-11,
                travel_tolerance_m: 1e-11,
                maximum_evaluations: 40_000,
                maximum_root_iterations: 96,
            },
        },
        velocity_tolerance_m_s: 1e-9,
        work_tolerance_j: 1e-9,
    }
}
fn supported(velocity: Vec3, omega: Vec3) -> SpherePlaneState {
    state(
        Vec3::new(0.0, RADIUS, 0.0),
        velocity,
        omega,
        UnitQuaternion::IDENTITY,
    )
}

#[test]
fn exact_arrest_advances_original_pose_attitude_and_boosted_plane_work() {
    let law = EllipticCoulomb::new([0.4; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    for boost in [Vec3::ZERO, Vec3::new(0.3, 0.2, -0.1)] {
        let mut scene = scene();
        scene.velocity_world = boost;
        let before = supported(boost.add(Vec3::new(1.0, 0.0, 0.0)), Vec3::ZERO);
        let mut history = law.initial_state();
        let result = advance_sphere_plane_directional(
            &scene,
            before,
            &law,
            &mut history,
            frame(),
            0.4,
            controls(),
            || false,
        )
        .unwrap();
        let legacy = advance_sphere_plane_with_support_friction(
            &scene,
            before,
            &FrictionLaw::Coulomb {
                static_mu: 0.4,
                kinetic_mu: 0.2,
            },
            0.4,
            controls().motion,
            || false,
        )
        .unwrap();
        let after = result.motion.state_after.body();
        close_vector(
            after.pose().position_world(),
            Vec3::new(15.0 / 49.0, RADIUS, 0.0).add(boost.scale(0.4)),
        );
        close_vector(velocity(after), boost.add(Vec3::new(5.0 / 7.0, 0.0, 0.0)));
        close_vector(spin(after), Vec3::new(0.0, 0.0, -10.0 / 7.0));
        let expected_orientation =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), -23.0 / 49.0).unwrap();
        close_vector(
            after
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            expected_orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        );
        close_vector(
            after.pose().position_world(),
            legacy.state_after.body().pose().position_world(),
        );
        close_vector(velocity(after), velocity(legacy.state_after.body()));
        close(result.sustained_friction_heat_j, 2.0 / 7.0);
        close(result.impact_loss_j, 0.0);
        close(result.support_distance_m, 1.0 / 14.0);
        close(result.support_slip_m[0], 1.0 / 14.0);
        close(history.slip_m()[0], 1.0 / 14.0);
        assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Rolling);
        let arrest = result
            .motion
            .events
            .iter()
            .find(|e| e.kind == SpherePlaneEventKind::Rolling)
            .unwrap();
        close(arrest.time_s, 1.0 / 7.0);
        assert!(result.support_phases > 1);
        assert!(result.segments >= result.support_phases);
        balances(&scene, before, 0.4, &result.motion);
    }
}

#[test]
fn original_normal_support_holds_above_kinetic_limit_then_driven_breakaway_moves_body() {
    let law = EllipticCoulomb::new([0.4, 0.5], [0.2, 0.3], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.force_world.x = 20.0;
    let before = supported(Vec3::ZERO, Vec3::ZERO);
    let mut history = law.initial_state();
    let held = advance_sphere_plane_directional(
        &scene,
        before,
        &law,
        &mut history,
        frame(),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(held.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close(held.sustained_friction_heat_j, 0.0);
    close(held.support_distance_m, 0.0);
    close(held.motion.contact_impulse_world.y, 2.0);
    close(held.motion.contact_impulse_world.x, -4.0 / 7.0);
    close(held.motion.external_force_work, 5.0 / 7.0);
    close(velocity(held.motion.state_after.body()).x, 5.0 / 7.0);
    close(
        held.motion.state_after.body().pose().position_world().x,
        1.0 / 28.0,
    );
    balances(&scene, before, 0.1, &held.motion);

    scene.force_world.x = 40.0;
    let before = held.motion.state_after;
    let result = advance_sphere_plane_directional(
        &scene,
        before,
        &law,
        &mut history,
        frame(),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Sliding);
    close(result.support_distance_m, 0.065);
    close(result.sustained_friction_heat_j, 0.26);
    close(
        result.motion.external_force_work,
        40.0 * (0.1 * 5.0 / 7.0 + 0.09),
    );
    close(
        velocity(result.motion.state_after.body()).x,
        5.0 / 7.0 + 1.8,
    );
    close(spin(result.motion.state_after.body()).z, -10.0 / 7.0 - 1.0);
    close(
        result.motion.state_after.body().pose().position_world().x,
        1.0 / 28.0 + 0.1 * 5.0 / 7.0 + 0.09,
    );
    close(history.slip_m()[0], 0.065);
    let onset = result
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Sliding)
        .unwrap();
    close(onset.time_s, 0.1);
    balances(&scene, before, 0.1, &result.motion);
}

#[test]
fn directional_curve_keeps_plane_texture_gauge_and_common_translation_covariance() {
    let mut reference = None;
    for (rotated, boost) in [
        (false, Vec3::ZERO),
        (true, Vec3::ZERO),
        (true, Vec3::new(0.7, 0.2, -0.3)),
    ] {
        let frame = if rotated {
            TangentFrame::new([0.0, 1.0, 0.0], [0.6, 0.0, -0.8]).unwrap()
        } else {
            frame()
        };
        let axis = if rotated { [1.0, 0.0] } else { [0.6, 0.8] };
        let law = EllipticCoulomb::new([0.5, 0.6], [0.15, 0.45], axis).unwrap();
        let mut history = law.initial_state();
        let mut scene = scene();
        scene.force_world = Vec3::new(0.6, -20.0, -0.8);
        scene.velocity_world = boost;
        let before = supported(
            boost.add(Vec3::new(1.0, 0.0, 0.4)),
            Vec3::new(0.1, 0.3, -0.2),
        );
        let result = advance_sphere_plane_directional(
            &scene,
            before,
            &law,
            &mut history,
            frame,
            0.08,
            controls(),
            || false,
        )
        .unwrap();
        let after = result.motion.state_after.body();
        assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Sliding);
        assert!(result.motion.sliding_substeps > result.support_phases);
        assert!(
            result.support_distance_m
                > result.support_slip_m[0].hypot(result.support_slip_m[1]) + 1e-6
        );
        let values = (
            after.pose().position_world().sub(boost.scale(0.08)),
            velocity(after).sub(boost),
            spin(after),
            after
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            result.support_distance_m,
            result.sustained_friction_heat_j,
        );
        if let Some((x, v, omega, direction, path, heat)) = reference {
            close_vector(values.0, x);
            close_vector(values.1, v);
            close_vector(values.2, omega);
            close_vector(values.3, direction);
            close(values.4, path);
            close(values.5, heat);
        } else {
            reference = Some(values);
        }
        balances(&scene, before, 0.08, &result.motion);
    }
}

#[test]
fn moving_segment_discretization_converges_to_independent_curved_trajectory() {
    let mut scene = scene();
    scene.force_world = Vec3::new(0.0, -20.0, 2.0);
    let law = EllipticCoulomb::new([0.2; 2], [2.0 / 35.0; 2], [1.0, 0.0]).unwrap();
    let before = supported(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO);
    // v'=a-2v/|v| has the exact w=tan(theta/2) solution used by the
    // original curve test: T=19/48, with independently integrated COM/spin.
    let duration = 19.0 / 48.0;
    let expected_v = Vec3::new(11.0 / 14.0, 0.0, 113.0 / 336.0);
    let expected_x = Vec3::new(1183.0 / 3360.0, RADIUS, 2291.0 / 32256.0);
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let mut history = law.initial_state();
        let mut control = controls();
        control.maximum_support_phase_s = duration;
        control.driven.maximum_step_s = duration / f64::from(count);
        let result = advance_sphere_plane_directional(
            &scene,
            before,
            &law,
            &mut history,
            frame(),
            duration,
            control,
            || false,
        )
        .unwrap();
        errors.push([
            velocity(result.motion.state_after.body())
                .sub(expected_v)
                .norm_squared()
                .sqrt(),
            result
                .motion
                .state_after
                .body()
                .pose()
                .position_world()
                .sub(expected_x)
                .norm_squared()
                .sqrt(),
            (result.sustained_friction_heat_j - 129.0 / 448.0).abs(),
        ]);
        assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Sliding);
        assert_eq!(result.support_phases, 1);
        balances(&scene, before, duration, &result.motion);
    }
    for pair in errors.windows(2) {
        for metric in 0..3 {
            let ratio = pair[0][metric] / pair[1][metric];
            assert!(
                (1.8..2.2).contains(&ratio),
                "moving curve errors: {errors:?}"
            );
        }
    }
    assert!(errors[2][0] < 0.002 && errors[2][1] < 0.0005 && errors[2][2] < 0.001);
}

#[test]
fn original_impacts_and_flight_keep_loss_separate_from_sustained_friction() {
    let law = EllipticCoulomb::new([0.9; 2], [0.7; 2], [1.0, 0.0]).unwrap();
    for restitution in [0.5, 0.0] {
        let mut scene = scene();
        scene.restitution = restitution;
        let before = state(
            Vec3::new(0.0, 0.9, 0.0),
            Vec3::new(1.0, -1.0, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        );
        let mut history = law.initial_state();
        let initial = history.clone();
        let result = advance_sphere_plane_directional(
            &scene,
            before,
            &law,
            &mut history,
            frame(),
            0.3,
            controls(),
            || false,
        )
        .unwrap();
        close(
            result.impact_loss_j,
            9.0 * (1.0 - restitution * restitution) + 2.0 / 7.0,
        );
        close(result.sustained_friction_heat_j, 0.0);
        close(result.support_distance_m, 0.0);
        assert_eq!(history, initial);
        let impact = result
            .motion
            .events
            .iter()
            .find(|e| e.kind == SpherePlaneEventKind::Impact)
            .unwrap();
        close(impact.time_s, 0.2);
        close_vector(
            impact.tangent_impulse_world,
            Vec3::new(-4.0 / 7.0, 0.0, 0.0),
        );
        if restitution == 0.5 {
            assert_eq!(result.support_phases, 0);
            assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Flight);
        } else {
            assert!(result.support_phases > 0);
            assert_eq!(result.motion.state_after.mode(), SpherePlaneMode::Rolling);
        }
        balances(&scene, before, 0.3, &result.motion);
    }
}

#[test]
fn late_source_cancel_budget_and_report_refusals_preserve_original_history() {
    let law = EllipticCoulomb::new([0.4; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let mut scene = scene();
    scene.force_world.x = 40.0;
    let before = supported(Vec3::ZERO, Vec3::ZERO);
    let initial = (law.initial_state(), 17usize);
    for failure in 0..4 {
        let mut history = initial.clone();
        let mut control = controls();
        control.maximum_support_phase_s = 0.02;
        if failure == 2 {
            control.maximum_support_phases = 1;
        }
        let visits = std::cell::Cell::new(0usize);
        let cancel = std::cell::Cell::new(false);
        let result = advance_sphere_plane_directional_with(
            &scene,
            before,
            &mut history,
            frame(),
            0.1,
            control,
            |state: &(EllipticCoulombState, usize), phase, visit, poll| {
                if phase.time_s >= 0.02 {
                    if failure == 0 {
                        return Err("later source refusal");
                    }
                    if failure == 1 {
                        cancel.set(true);
                    }
                }
                let mut report = law
                    .advance_driven_segments(
                        &state.0,
                        &phase.input,
                        control.driven,
                        &mut || poll(),
                        &mut |segment| {
                            visits.set(visits.get() + 1);
                            visit(segment)
                        },
                    )
                    .map_err(|_| "native/visitor refusal")?;
                if failure == 3 {
                    report.dissipated_j += 1.0;
                }
                let next = (report.state.clone(), state.1 + 1);
                Ok((report, next))
            },
            || cancel.get(),
        );
        assert!(result.is_err(), "failure {failure} unexpectedly published");
        assert!(visits.get() > 0, "must stage original sphere motion first");
        assert_eq!(history, initial);
    }
    let mut history = initial;
    let result = advance_sphere_plane_directional_with(
        &scene,
        before,
        &mut history,
        frame(),
        0.1,
        controls(),
        |state, phase, visit, poll| {
            close(phase.input.normal_n, 20.0);
            close(phase.input.mobility_per_kg[0][0], 7.0 / 4.0);
            let report = law.advance_driven_segments(
                &state.0,
                &phase.input,
                controls().driven,
                &mut || poll(),
                &mut |segment| visit(segment),
            )?;
            let next = (report.state.clone(), state.1 + 1);
            Ok::<_, EllipticDrivenSegmentError<SpherePlaneError>>((report, next))
        },
        || false,
    )
    .unwrap();
    assert_eq!(history.1, 17 + result.support_phases);
    close(result.support_distance_m, 0.065);
    close(result.sustained_friction_heat_j, 0.26);
    balances(&scene, before, 0.1, &result.motion);
}

//! Prescribed rigid motion through the original panel/Kelvin wake owner.

use fs_bem::BemError;
use fs_bem::panel2d::{PreparedPanel2d, naca4_symmetric};
use fs_bem::wake2d::{PanelFrame2d, PanelWakeAdaptiveOptions, PanelWakeTransportOptions, WakeSim};

fn transport(substeps: usize) -> PanelWakeTransportOptions {
    PanelWakeTransportOptions {
        substeps,
        max_coupled_solves: 2 * substeps + 1,
        max_displacement_fraction: None,
    }
}

fn adaptive(tolerance: f64) -> PanelWakeAdaptiveOptions {
    PanelWakeAdaptiveOptions {
        initial_substeps: 1,
        max_substeps: 256,
        max_coupled_solves: 2048,
        absolute_position_tolerance: tolerance,
        relative_position_tolerance: 0.0,
        max_displacement_fraction: None,
    }
}

fn assert_same_body(a: &WakeSim, b: &WakeSim) {
    assert_eq!(
        a.steady_circulation().to_bits(),
        b.steady_circulation().to_bits()
    );
    assert_eq!(a.history().len(), b.history().len());
    for (a, b) in a.history().iter().zip(b.history()) {
        assert_eq!(a.t.to_bits(), b.t.to_bits());
        assert_eq!(a.bound.to_bits(), b.bound.to_bits());
        assert_eq!(a.peak_speed.to_bits(), b.peak_speed.to_bits());
        assert_eq!(a.vortices, b.vortices);
    }
    assert_eq!(a.wake().len(), b.wake().len());
    for (a, b) in a.wake().iter().zip(b.wake()) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.gamma.to_bits(), b.gamma.to_bits());
    }
    let a = a.panel_state().unwrap();
    let b = b.panel_state().unwrap();
    assert_eq!(a.sources, b.sources);
    assert_eq!(a.tangential_velocity, b.tangential_velocity);
    assert_eq!(a.bound_circulation.to_bits(), b.bound_circulation.to_bits());
    assert_eq!(
        a.max_normal_residual.to_bits(),
        b.max_normal_residual.to_bits()
    );
    assert_eq!(a.kutta_residual.to_bits(), b.kutta_residual.to_bits());
}

fn assert_unchanged(a: &WakeSim, b: &WakeSim) {
    assert_same_body(a, b);
    assert_eq!(a.panel_frame(), b.panel_frame());
    assert_eq!(a.trace_json(1).unwrap(), b.trace_json(1).unwrap());
}

#[test]
fn identity_frame_preserves_fixed_and_adaptive_native_results() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut fixed = WakeSim::new_panel_coupled(&foil, 0.0, 0.02, 0.04).unwrap();
    let mut moving =
        WakeSim::new_panel_coupled_with_motion(&foil, 0.02, 0.04, PanelFrame2d::default()).unwrap();
    for substeps in [1, 2, 4] {
        let a = fixed
            .step_panel_coupled_with_transport(0.0, transport(substeps), || false)
            .unwrap();
        let b = moving
            .step_panel_coupled_with_motion(
                transport(substeps),
                |_| Ok(PanelFrame2d::default()),
                |_, _| Ok([1.0, 0.0]),
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        assert_same_body(&fixed, &moving);
    }
    let a = fixed
        .step_panel_coupled_adaptive(0.0, adaptive(1e-7), || false)
        .unwrap();
    let b = moving
        .step_panel_coupled_adaptive_with_motion(
            adaptive(1e-7),
            |_| Ok(PanelFrame2d::default()),
            |_, _| Ok([1.0, 0.0]),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_same_body(&fixed, &moving);
    let trace = moving.trace_json(1).unwrap();
    assert!(trace.contains("\"model\":\"inviscid-panel-kelvin-moving-frame\""));
    assert!(trace.contains("\"coordinate_frame\":\"moving-body-chord\""));
    assert!(trace.contains("\"alpha_role\":\"nominal-reference-diagnostic\""));
    assert!(trace.contains("\"peak_body_relative_speed_reference_ratio\":"));
    assert!(!trace.contains("\"peak_speed_freestream_ratio\":"));
    assert!(
        !fixed
            .trace_json(1)
            .unwrap()
            .contains("\"coordinate_frame\"")
    );
    for (body, world) in moving.wake().iter().zip(moving.wake_in_world().unwrap()) {
        assert_eq!(body.pos.map(f64::to_bits), world.pos.map(f64::to_bits));
        assert_eq!(body.gamma.to_bits(), world.gamma.to_bits());
    }
}

fn translation(time: f64, boost: [f64; 2]) -> PanelFrame2d {
    let velocity = [0.125 + boost[0], 0.0625 + boost[1]];
    PanelFrame2d {
        origin_world: velocity.map(|v| v * time),
        origin_velocity_world: velocity,
        ..PanelFrame2d::default()
    }
}

#[test]
fn translating_heaving_body_matches_relative_air_and_is_galilean_invariant() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let (dt, core, boost) = (0.025, 0.04, [64.0, -32.0]);
    let mut fixed = WakeSim::new_panel_coupled(&foil, 0.0, dt, core).unwrap();
    let mut moving =
        WakeSim::new_panel_coupled_with_motion(&foil, dt, core, translation(0.0, [0.0; 2]))
            .unwrap();
    let mut boosted =
        WakeSim::new_panel_coupled_with_motion(&foil, dt, core, translation(0.0, boost)).unwrap();
    for step in 1..=5 {
        fixed
            .step_panel_coupled_with_incident(
                0.0,
                transport(4),
                |_, _| Ok([-0.125, -0.0625]),
                || false,
            )
            .unwrap();
        for (wake, boost) in [(&mut moving, [0.0; 2]), (&mut boosted, boost)] {
            wake.step_panel_coupled_with_motion(
                transport(4),
                |t| Ok(translation(t, boost)),
                |_, _| Ok([1.0 + boost[0], boost[1]]),
                || false,
            )
            .unwrap();
        }
        assert_same_body(&moving, &fixed);
        assert_same_body(&boosted, &moving);
        let time = step as f64 * dt;
        assert_eq!(moving.panel_frame(), Some(translation(time, [0.0; 2])));
        let exported = moving.wake_in_world().unwrap();
        let boosted_export = boosted.wake_in_world().unwrap();
        for ((body, world), changed) in moving.wake().iter().zip(exported).zip(boosted_export) {
            for axis in 0..2 {
                let expected = body.pos[axis] + translation(time, [0.0; 2]).origin_world[axis];
                assert!((world.pos[axis] - expected).abs() < 2e-15);
                assert!((changed.pos[axis] - world.pos[axis] - boost[axis] * time).abs() < 3e-15);
            }
            assert_eq!(body.gamma.to_bits(), world.gamma.to_bits());
            assert_eq!(world.gamma.to_bits(), changed.gamma.to_bits());
        }
    }
}

// Quarter-chord pitch plus heave, with analytic velocity of the LE origin.
fn pitch_heave(time: f64) -> PanelFrame2d {
    let angle = 0.06 * (4.0 * time).sin();
    let omega = 0.24 * (4.0 * time).cos();
    let (s, c) = angle.sin_cos();
    PanelFrame2d {
        origin_world: [0.25 - 0.25 * c, 0.025 * (3.0 * time).sin() - 0.25 * s],
        angle_rad: angle,
        origin_velocity_world: [
            0.25 * s * omega,
            0.075 * (3.0 * time).cos() - 0.25 * c * omega,
        ],
        angular_velocity: omega,
    }
}

// Smooth potential flow: spatial strain has zero divergence and zero curl.
fn world_air(point: [f64; 2], time: f64) -> [f64; 2] {
    let strain = 0.02 * (2.0 * time).sin();
    [
        1.0 + 0.025 * (5.0 * time).sin() + strain * point[0],
        0.04 * (3.0 * time).sin() - strain * point[1],
    ]
}

fn assert_moving_boundary(wake: &WakeSim, panels: &PreparedPanel2d, core: f64) {
    let time = wake.history().last().unwrap().t;
    let frame = pitch_heave(time);
    assert_eq!(wake.panel_frame(), Some(frame));
    let (s, c) = frame.angle_rad.sin_cos();
    let particles: Vec<_> = wake.wake().iter().map(|v| v.as_particle()).collect();
    let incident: Vec<_> = panels
        .midpoints()
        .iter()
        .map(|&r| {
            let world = [
                frame.origin_world[0] + c * r[0] - s * r[1],
                frame.origin_world[1] + s * r[0] + c * r[1],
            ];
            let air = world_air(world, time);
            let slip = [
                air[0] - frame.origin_velocity_world[0],
                air[1] - frame.origin_velocity_world[1],
            ];
            let induced = fs_vpm::induced_velocity(&particles, r, core);
            [
                c * slip[0] + s * slip[1] + frame.angular_velocity * r[1] + induced[0],
                -s * slip[0] + c * slip[1] - frame.angular_velocity * r[0] + induced[1],
            ]
        })
        .collect();
    // Independent Cartesian transform and ordinary solve of the FULL endpoint
    // wake, without the native coupled-newborn response or motion helper.
    let expected = panels.solve_incident(&incident).unwrap();
    let actual = wake.panel_state().unwrap();
    let tolerance = 6e-11 * (1.0 + expected.bound_circulation.abs());
    assert!((actual.bound_circulation - expected.bound_circulation).abs() < tolerance);
    for (a, b) in actual.sources.iter().zip(expected.sources) {
        assert!((a - b).abs() < tolerance);
    }
    assert!(actual.max_normal_residual < tolerance);
    assert!(actual.kutta_residual.abs() < tolerance);
    // The frame's -omega*perp(r) is kinematic, NOT another vortex inventory.
    let kelvin = actual.bound_circulation + wake.wake().iter().map(|v| v.gamma).sum::<f64>();
    assert!(kelvin.abs() < tolerance);
    for (body, world) in wake.wake().iter().zip(wake.wake_in_world().unwrap()) {
        let expected = [
            frame.origin_world[0] + c * body.pos[0] - s * body.pos[1],
            frame.origin_world[1] + s * body.pos[0] + c * body.pos[1],
        ];
        assert!((world.pos[0] - expected[0]).hypot(world.pos[1] - expected[1]) < 2e-14);
        assert_eq!(body.gamma.to_bits(), world.gamma.to_bits());
    }
}

fn advance_motion(wake: &mut WakeSim, substeps: usize) {
    wake.step_panel_coupled_with_motion(
        transport(substeps),
        |time| Ok(pitch_heave(time)),
        |point, time| Ok(world_air(point, time)),
        || false,
    )
    .unwrap();
}

#[test]
fn pitch_and_heave_drive_local_boundary_velocity_shedding_and_world_trajectories() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let (dt, core) = (0.02, 0.04);
    let mut wake =
        WakeSim::new_panel_coupled_with_motion(&foil, dt, core, pitch_heave(0.0)).unwrap();
    let mut samples = Vec::new();
    wake.step_panel_coupled_with_motion(
        transport(2),
        |time| {
            samples.push(time);
            Ok(pitch_heave(time))
        },
        |point, time| Ok(world_air(point, time)),
        || false,
    )
    .unwrap();
    // A single prepared frame serves every panel and wake point at a stage.
    assert_eq!(samples, vec![0.0, 0.25 * dt, 0.5 * dt, 0.75 * dt, dt]);
    assert_moving_boundary(&wake, &panels, core);
    let first = wake.panel_state().unwrap().bound_circulation;
    for _ in 0..4 {
        let old: Vec<_> = wake.wake().iter().map(|v| v.gamma.to_bits()).collect();
        advance_motion(&mut wake, 4);
        assert_moving_boundary(&wake, &panels, core);
        assert_eq!(wake.wake().len(), old.len() + 1);
        for (v, bits) in wake.wake().iter().zip(old) {
            assert_eq!(v.gamma.to_bits(), bits);
        }
    }
    assert!((wake.panel_state().unwrap().bound_circulation - first).abs() > 1e-4);
    // A stationary foil with the same changing freestream angle cannot model
    // pitch/heave: its boundary omits the rotating point's material velocity.
    let mut angle_only = WakeSim::new_panel_coupled(&foil, 0.0, dt, core).unwrap();
    for step in 1..=5 {
        angle_only
            .step_panel_coupled_at_angle(-pitch_heave(step as f64 * dt).angle_rad)
            .unwrap();
    }
    assert!(
        (wake.panel_state().unwrap().bound_circulation
            - angle_only.panel_state().unwrap().bound_circulation)
            .abs()
            > 1e-3
    );
}

fn world_position_error(a: &WakeSim, b: &WakeSim) -> f64 {
    let a = a.wake_in_world().unwrap();
    let b = b.wake_in_world().unwrap();
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(a, b)| (a.pos[0] - b.pos[0]).hypot(a.pos[1] - b.pos[1]))
        .fold(0.0, f64::max)
}

#[test]
fn moving_transport_converges_and_adaptive_retries_publish_one_endpoint_frame() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let panels = PreparedPanel2d::new(&foil).unwrap();
    let mut start =
        WakeSim::new_panel_coupled_with_motion(&foil, 0.02, 0.04, pitch_heave(0.0)).unwrap();
    for _ in 0..3 {
        advance_motion(&mut start, 4);
    }
    let mut reference = start.clone();
    advance_motion(&mut reference, 128);
    let mut errors = Vec::new();
    let mut trials = Vec::new();
    for substeps in [1, 2, 4, 8] {
        let mut result = start.clone();
        advance_motion(&mut result, substeps);
        assert_moving_boundary(&result, &panels, 0.04);
        assert_eq!(result.wake().len(), start.wake().len() + 1);
        assert_eq!(result.history().len(), start.history().len() + 1);
        errors.push(world_position_error(&result, &reference));
        trials.push(result);
    }
    assert!(
        errors[0] > 1e-9,
        "fixture must resolve temporal error: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.4 * pair[0], "midpoint refinement: {errors:?}");
    }
    let estimate = world_position_error(&trials[0], &trials[1]) / 3.0;
    let mut result = start.clone();
    let options = adaptive(estimate / 64.0);
    let report = result
        .step_panel_coupled_adaptive_with_motion(
            options,
            |time| Ok(pitch_heave(time)),
            |point, time| Ok(world_air(point, time)),
            || false,
        )
        .unwrap();
    assert!(report.rejected_comparisons > 0);
    assert!(report.estimated_max_position_error <= options.absolute_position_tolerance);
    assert_eq!(result.history().len(), start.history().len() + 1);
    assert_eq!(result.wake().len(), start.wake().len() + 1);
    assert_moving_boundary(&result, &panels, 0.04);
    let mut fixed = start.clone();
    advance_motion(&mut fixed, report.substeps);
    assert_unchanged(&result, &fixed);
    assert!(world_position_error(&result, &reference) < 0.1 * errors[0]);
    for (old, accepted) in start.wake().iter().zip(result.wake()) {
        assert_eq!(old.gamma.to_bits(), accepted.gamma.to_bits());
    }
}

#[test]
fn frame_mode_mismatch_and_all_fixed_step_paths_refuse_atomically() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut start =
        WakeSim::new_panel_coupled_with_motion(&foil, 0.02, 0.04, pitch_heave(0.0)).unwrap();
    advance_motion(&mut start, 2);
    let mut wake = start.clone();
    assert!(wake.step().is_err());
    assert_unchanged(&wake, &start);
    assert!(wake.step_panel_coupled_at_angle(0.1).is_err());
    assert_unchanged(&wake, &start);
    assert!(
        wake.step_panel_coupled_with_transport(0.1, transport(2), || false)
            .is_err()
    );
    assert_unchanged(&wake, &start);
    assert!(
        wake.step_panel_coupled_with_incident(0.1, transport(2), |_, _| Ok([0.0; 2]), || false)
            .is_err()
    );
    assert_unchanged(&wake, &start);
    assert!(
        wake.step_panel_coupled_adaptive(0.1, adaptive(1e-5), || false)
            .is_err()
    );
    assert_unchanged(&wake, &start);
    assert!(
        wake.step_panel_coupled_adaptive_with_incident(
            0.1,
            adaptive(1e-5),
            |_, _| Ok([0.0; 2]),
            || false
        )
        .is_err()
    );
    assert_unchanged(&wake, &start);

    let mut fixed = WakeSim::new_panel_coupled(&foil, 0.0, 0.02, 0.04).unwrap();
    fixed.step().unwrap();
    let before = fixed.clone();
    for adaptive_mode in [false, true] {
        let calls = std::cell::Cell::new(0);
        let mut motion = |_| {
            calls.set(calls.get() + 1);
            Ok(PanelFrame2d::default())
        };
        let result = if adaptive_mode {
            fixed
                .step_panel_coupled_adaptive_with_motion(
                    adaptive(1e-5),
                    &mut motion,
                    |_, _| Ok([1.0, 0.0]),
                    || false,
                )
                .map(|_| ())
        } else {
            fixed
                .step_panel_coupled_with_motion(
                    transport(2),
                    &mut motion,
                    |_, _| Ok([1.0, 0.0]),
                    || false,
                )
                .map(|_| ())
        };
        assert!(result.is_err());
        assert_eq!(calls.get(), 0);
        assert_unchanged(&fixed, &before);
    }
}

#[test]
fn discontinuous_pose_endpoint_failure_and_motion_cancellation_keep_accepted_state() {
    let foil = naca4_symmetric(0.12, 24).unwrap();
    let mut start =
        WakeSim::new_panel_coupled_with_motion(&foil, 0.02, 0.04, pitch_heave(0.0)).unwrap();
    advance_motion(&mut start, 2);
    for refusal in 0..5 {
        let mut wake = start.clone();
        let cancelled = std::cell::Cell::new(false);
        let mut visits = 0;
        let result = wake.step_panel_coupled_with_motion(
            transport(2),
            |time| {
                visits += 1;
                let mut frame = pitch_heave(time);
                if refusal == 0 {
                    frame.origin_world[0] += 0.1;
                } else if refusal == 1 {
                    frame.angle_rad += 2.0 * std::f64::consts::PI;
                } else if time == 0.04 {
                    if refusal == 2 {
                        return Err(BemError::NonFiniteResult {
                            operation: "test endpoint motion refusal",
                        });
                    } else if refusal == 3 {
                        frame.angular_velocity = f64::NAN;
                    } else {
                        cancelled.set(true);
                    }
                }
                Ok(frame)
            },
            |point, time| Ok(world_air(point, time)),
            || cancelled.get(),
        );
        assert!(result.is_err(), "refusal {refusal}");
        assert_eq!(visits, if refusal < 2 { 1 } else { 5 });
        assert_unchanged(&wake, &start);
    }
    let mut wake = start.clone();
    assert!(
        wake.step_panel_coupled_adaptive_with_motion(
            adaptive(1e-6),
            |time| Ok(pitch_heave(time)),
            |point, time| {
                if time == 0.04 {
                    Err(BemError::NonFiniteResult {
                        operation: "test endpoint air refusal",
                    })
                } else {
                    Ok(world_air(point, time))
                }
            },
            || false,
        )
        .is_err()
    );
    assert_unchanged(&wake, &start);
}

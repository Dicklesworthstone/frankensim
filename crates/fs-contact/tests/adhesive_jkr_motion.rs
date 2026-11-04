//! Original JKR dynamics: capture, rebound, smooth accuracy and bounded refusal.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{motion::*, *};
use fs_exec::CancelGate;

fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.05).normal_state().clone();
    SourceJkrSphereContact::new(
        normal,
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}

const MASS: f64 = 1.2e-5;
fn scales(normal: &SourceJkrSphereContact) -> (f64, f64, f64, f64) {
    let length = normal.indentation_scale_m();
    let force = PI * 0.05 * 0.01 / 6.0;
    let speed = (force * length / MASS).sqrt();
    (length, force, speed, length / speed)
}
fn budget(normal: &SourceJkrSphereContact) -> JkrNormalMotionBudget {
    let (length, force, speed, time) = scales(normal);
    JkrNormalMotionBudget {
        maximum_step_s: 10.0 * time,
        maximum_substep_s: 0.01 * time,
        maximum_substeps: 4096,
        maximum_root_iterations: 48,
        maximum_events: MAX_JKR_NORMAL_EVENTS,
        minimum_indentation_m: -100.0 * length,
        maximum_accumulated_travel_m: 200.0 * length,
        maximum_normal_speed_m_s: 100.0 * speed,
        maximum_force_n: 1e4 * force,
        coordinate_absolute_tolerance_m: 1e-14 * length,
        coordinate_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * length,
        energy_relative_tolerance: 1e-10,
        impulse_absolute_tolerance_n_s: 1e-12 * MASS * speed,
        impulse_relative_tolerance: 1e-11,
    }
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, tolerance={:e}",
        tolerance * scale
    );
}

#[test]
fn detached_flight_is_exact_with_held_load_and_no_contact_moments() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, budget(&normal)).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-8.0 * length).unwrap(), -speed)
        .unwrap();
    let dt = 0.5 * time;
    let load = 2.0 * force;
    let step = motion
        .advance(initial, load, dt, &CancelGate::new())
        .unwrap();
    close(
        step.contact.indentation_m,
        -8.0 * length - dt * speed + 0.5 * load / MASS * dt * dt,
        length,
        2e-13,
    );
    close(
        step.state().velocity_m_s(),
        -speed + load / MASS * dt,
        speed,
        2e-13,
    );
    assert_eq!(step.normal_impulse_n_s, 0.0);
    assert_eq!(step.normal_time_moment_n_s2, 0.0);
    assert_eq!(step.released_energy_j, 0.0);
    assert_eq!(step.event_count(), 0);
    assert_eq!(step.initial_state(), initial);
    assert_eq!(step.state().time_s(), dt);
    assert_eq!(step.state().steps(), 1);
}

#[test]
fn touching_at_step_end_has_one_original_snap_and_continuous_velocity() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let dt = 0.4 * time;
    let mut controls = budget(&normal);
    controls.maximum_substep_s = dt;
    let motion = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let velocity = 0.5 * speed;
    let initial = motion
        .initialize(normal.detached_state(-velocity * dt).unwrap(), velocity)
        .unwrap();
    let step = motion
        .advance(initial, 0.0, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(step.contact.branch, JkrBranch::Attached);
    assert_eq!(step.event_count(), 1);
    let event = step.events().next().unwrap();
    assert_eq!(event.before().branch(), JkrBranch::Detached);
    assert_eq!(event.after().branch(), JkrBranch::Attached);
    close(event.offset_s, dt, time, 1e-14);
    close(event.velocity_m_s, velocity, speed, 1e-14);
    assert_eq!(event.before().indentation_m(), 0.0);
    assert_eq!(event.after().indentation_m(), 0.0);
    assert_eq!(step.normal_impulse_n_s, 0.0);
    assert_eq!(step.normal_time_moment_n_s2, 0.0);
    close(
        step.released_energy_j,
        -normal.evaluate_attached(0.0).unwrap().reversible_energy_j,
        force * length,
        1e-12,
    );
    let next = motion
        .advance(step.state(), 0.0, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(next.event_count(), 0);
    assert_eq!(next.released_energy_j, 0.0);
}

#[test]
fn collision_rebounds_with_original_snap_energy_and_work_conjugate_impulses() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, budget(&normal)).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-0.25 * length).unwrap(), 10.0 * speed)
        .unwrap();
    let dt = 8.0 * time;
    let step = motion
        .advance(initial, 0.0, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(step.contact.branch, JkrBranch::Detached);
    assert_eq!(step.event_count(), 2);
    let mut events = step.events();
    let touch = events.next().unwrap();
    let escape = events.next().unwrap();
    assert!(touch.offset_s < escape.offset_s);
    assert_eq!(
        escape.before().indentation_m(),
        normal.minimum_attached_indentation_m()
    );
    assert_eq!(escape.before().branch(), JkrBranch::Attached);
    assert_eq!(escape.after().branch(), JkrBranch::Detached);
    let original_loss = -normal.evaluate_attached(0.0).unwrap().reversible_energy_j
        + normal
            .evaluate_attached(normal.minimum_attached_indentation_m())
            .unwrap()
            .reversible_energy_j;
    close(step.released_energy_j, original_loss, force * length, 1e-12);
    let expected_velocity = -(initial.velocity_m_s().powi(2) - 2.0 * original_loss / MASS).sqrt();
    close(step.state().velocity_m_s(), expected_velocity, speed, 2e-9);
    close(escape.velocity_m_s, expected_velocity, speed, 2e-9);
    close(
        step.normal_impulse_n_s,
        MASS * (initial.velocity_m_s() - expected_velocity),
        MASS * speed,
        2e-9,
    );
    let restored = dt * initial.velocity_m_s()
        - (dt * step.normal_impulse_n_s - step.normal_time_moment_n_s2) / MASS;
    close(step.normal_increment_m, restored, length, 2e-11);
    close(
        step.total_energy_j + step.released_energy_j,
        0.5 * MASS * initial.velocity_m_s().powi(2),
        force * length,
        2e-9,
    );
    assert!(step.accumulated_travel_m > step.normal_increment_m.abs());
}

#[test]
fn lower_energy_collision_is_captured_without_velocity_reset_or_extra_loss() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, budget(&normal)).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-0.25 * length).unwrap(), 5.0 * speed)
        .unwrap();
    let step = motion
        .advance(initial, 0.0, 8.0 * time, &CancelGate::new())
        .unwrap();
    assert_eq!(step.contact.branch, JkrBranch::Attached);
    assert_eq!(step.event_count(), 1);
    close(
        step.events().next().unwrap().velocity_m_s,
        5.0 * speed,
        speed,
        1e-14,
    );
    close(
        step.released_energy_j,
        -normal.evaluate_attached(0.0).unwrap().reversible_energy_j,
        force * length,
        1e-12,
    );
    assert!(step.accumulated_travel_m > 3.0 * step.normal_increment_m.abs());
    assert!(step.state().velocity_m_s().abs() > 0.01 * speed);
    close(
        step.total_energy_j + step.released_energy_j,
        0.5 * MASS * initial.velocity_m_s().powi(2),
        force * length,
        2e-9,
    );
}

#[test]
fn near_fold_outward_motion_localizes_release_without_a_tangent_at_the_fold() {
    let normal = source();
    let (length, _, speed, time) = scales(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, budget(&normal)).unwrap();
    let fold = normal.minimum_attached_indentation_m();
    let initial = motion
        .initialize(
            normal
                .evaluate_attached(fold + 1e-8 * length)
                .unwrap()
                .state(),
            -0.2 * speed,
        )
        .unwrap();
    let step = motion
        .advance(initial, 0.0, 0.02 * time, &CancelGate::new())
        .unwrap();
    assert_eq!(step.event_count(), 1);
    let event = step.events().next().unwrap();
    assert!(event.offset_s > 0.0 && event.offset_s < 1e-6 * time);
    assert_eq!(event.before().indentation_m(), fold);
    assert_eq!(event.after().indentation_m(), fold);
    assert_eq!(step.contact.branch, JkrBranch::Detached);
    assert_eq!(event.velocity_m_s, step.state().velocity_m_s());
    assert!(
        motion
            .initialize(normal.evaluate_attached(fold).unwrap().state(), -speed)
            .is_err()
    );
}

#[test]
fn force_equilibrium_has_nonzero_reaction_moments_and_zero_work() {
    let normal = source();
    let (length, force, _, time) = scales(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, budget(&normal)).unwrap();
    let response = normal.evaluate_attached(length).unwrap();
    let initial = motion.initialize(response.state(), 0.0).unwrap();
    let dt = 0.5 * time;
    let step = motion
        .advance(initial, response.normal_force_n, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(step.normal_increment_m, 0.0);
    assert_eq!(step.state().velocity_m_s(), 0.0);
    assert_eq!(step.normal_work_j, 0.0);
    assert_eq!(step.applied_work_j, 0.0);
    close(
        step.normal_impulse_n_s,
        response.normal_force_n * dt,
        force * time,
        1e-13,
    );
    close(
        step.normal_time_moment_n_s2,
        0.5 * response.normal_force_n * dt * dt,
        force * time * time,
        1e-13,
    );
}

#[test]
fn a_coincident_root_on_the_wrong_slope_is_refined_before_acceptance() {
    let normal = source();
    let (length, force, _, _) = scales(&normal);
    let from = normal.minimum_attached_indentation_m() + 0.01 * length;
    let response = normal.evaluate_attached(from).unwrap();
    let tangent = response.tangent_n_per_m.unwrap();
    assert!(tangent < 0.0);
    let dt = 2.5 * (MASS / -tangent).sqrt();
    let load = response.normal_force_n + 1e-5 * force;
    let velocity = -0.5 * dt * ((load - response.normal_force_n) / MASS);
    // This coarse coincidence solves the residual while its Jacobian is
    // negative. Its tiny turning arc remains inside the source support, so
    // source-domain checks alone cannot distinguish the wrong root.
    assert!(1.0 + dt * dt * tangent / (4.0 * MASS) < 0.0);
    let mut controls = budget(&normal);
    controls.maximum_substep_s = dt;
    let motion = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let state = motion.initialize(response.state(), velocity).unwrap();
    let step = motion.advance(state, load, dt, &CancelGate::new()).unwrap();
    assert!(step.substeps > 1 && step.trials > step.substeps);
    assert_eq!(step.event_count(), 0);
    assert!(step.contact.indentation_m < from - 1e-8 * length);
}

// An independent force-sampled Runge-Kutta trajectory checks time accuracy;
// it does not reuse the discrete-gradient residual or its path mean.
fn rk4(normal: &SourceJkrSphereContact, mut d: f64, mut v: f64, dt: f64) -> (f64, f64) {
    let h = dt / 4096.0;
    let a = |x| -normal.evaluate_attached(x).unwrap().normal_force_n / MASS;
    for _ in 0..4096 {
        let (k1d, k1v) = (v, a(d));
        let (k2d, k2v) = (v + 0.5 * h * k1v, a(d + 0.5 * h * k1d));
        let (k3d, k3v) = (v + 0.5 * h * k2v, a(d + 0.5 * h * k2d));
        let (k4d, k4v) = (v + h * k3v, a(d + h * k3d));
        d += h / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d);
        v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    }
    (d, v)
}

#[test]
fn smooth_trajectory_converges_at_second_order_to_independent_force_integration() {
    let normal = source();
    let (length, _, speed, time) = scales(&normal);
    let dt = 0.4 * time;
    let exact = rk4(&normal, 0.0, 0.3 * speed, dt);
    let mut errors = Vec::new();
    for resolution in [0.08, 0.04, 0.02] {
        let mut controls = budget(&normal);
        controls.maximum_substep_s = resolution * time;
        let motion = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
        let initial = motion
            .initialize(normal.evaluate_attached(0.0).unwrap().state(), 0.3 * speed)
            .unwrap();
        let step = motion
            .advance(initial, 0.0, dt, &CancelGate::new())
            .unwrap();
        assert_eq!(step.event_count(), 0);
        errors.push(
            ((step.contact.indentation_m - exact.0) / length)
                .hypot((step.state().velocity_m_s() - exact.1) / speed),
        );
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    assert!(
        errors[0] / errors[1] < 4.2 && errors[1] / errors[2] < 4.2,
        "{errors:?}"
    );
}

#[test]
fn cancellation_and_finite_budgets_leave_history_reusable() {
    let normal = source();
    let (length, _, speed, time) = scales(&normal);
    let controls = budget(&normal);
    let motion = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-0.25 * length).unwrap(), 10.0 * speed)
        .unwrap();
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        motion.advance(initial, 0.0, time, &cancelled),
        Err(JkrNormalMotionError::Cancelled)
    );
    let retry = motion
        .advance(initial, 0.0, time, &CancelGate::new())
        .unwrap();
    assert_eq!(
        retry,
        motion
            .advance(initial, 0.0, time, &CancelGate::new())
            .unwrap()
    );
    let foreign = JkrNormalMotion::new(&normal, 2.0 * MASS, controls).unwrap();
    assert_eq!(
        foreign.advance(initial, 0.0, time, &CancelGate::new()),
        Err(JkrNormalMotionError::StateMismatch)
    );
    let mut controls = controls;
    controls.maximum_substeps = 1;
    let limited = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let state = limited
        .initialize(initial.contact_state(), initial.velocity_m_s())
        .unwrap();
    assert_eq!(
        limited.advance(state, 0.0, time, &CancelGate::new()),
        Err(JkrNormalMotionError::TrialBudget)
    );
    controls.maximum_substeps = 4096;
    controls.maximum_events = 1;
    let limited = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let state = limited
        .initialize(initial.contact_state(), initial.velocity_m_s())
        .unwrap();
    assert_eq!(
        limited.advance(state, 0.0, 8.0 * time, &CancelGate::new()),
        Err(JkrNormalMotionError::EventBudget)
    );
}

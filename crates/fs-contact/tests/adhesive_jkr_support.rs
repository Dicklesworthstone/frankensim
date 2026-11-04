//! Finite machine stiffness: conservative support and original JKR events.
#[path = "source_point/common.rs"]
mod common;
#[path = "adhesive_jkr_support/damping.rs"]
mod damping;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{motion::*, *};
use fs_exec::CancelGate;

const MASS: f64 = 1.2e-5;
fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
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
fn accounts(motion: &JkrNormalMotion<'_>, step: &JkrNormalMotionStep, load: f64, dt: f64) {
    let (length, force, speed, _) = scales(motion.normal());
    let initial = step.initial_state();
    let response = motion.response(initial).unwrap();
    let support_before = motion.support_energy_j(initial).unwrap();
    let initial_energy =
        0.5 * MASS * initial.velocity_m_s().powi(2) + response.reversible_energy_j + support_before;
    close(
        step.total_energy_j + step.released_energy_j,
        initial_energy + load * step.normal_increment_m,
        force * length,
        2e-9,
    );
    close(
        step.support_energy_j - support_before,
        step.support_work_j,
        force * length,
        2e-9,
    );
    close(
        step.contact.reversible_energy_j - response.reversible_energy_j + step.released_energy_j,
        step.normal_work_j,
        force * length,
        2e-9,
    );
    close(
        MASS * (step.state().velocity_m_s() - initial.velocity_m_s()),
        load * dt - step.reaction_impulse_n_s(),
        MASS * speed,
        2e-9,
    );
    let restored = dt * initial.velocity_m_s()
        + (0.5 * load * dt * dt - dt * step.reaction_impulse_n_s()
            + step.reaction_time_moment_n_s2())
            / MASS;
    close(step.normal_increment_m, restored, length, 2e-10);
}

#[test]
fn detached_spring_oscillator_preserves_energy_and_converges_to_exact_solution() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let k = 4.0 * force / length;
    let anchor = -10.0 * length;
    let dt = time;
    let expected_d = anchor + 2.0 * length * 2.0_f64.cos();
    let expected_v = -4.0 * speed * 2.0_f64.sin();
    let mut errors = Vec::new();
    for resolution in [0.1, 0.05, 0.025] {
        let mut controls = budget(&normal);
        controls.maximum_substep_s = resolution * time;
        let motion = JkrNormalMotion::new_with_support(&normal, MASS, controls, k, anchor).unwrap();
        let initial = motion
            .initialize(normal.detached_state(anchor + 2.0 * length).unwrap(), 0.0)
            .unwrap();
        let step = motion
            .advance(initial, 0.0, dt, &CancelGate::new())
            .unwrap();
        assert_eq!(step.event_count(), 0);
        assert_eq!(step.normal_impulse_n_s, 0.0);
        assert_eq!(step.normal_work_j, 0.0);
        assert_eq!(step.released_energy_j, 0.0);
        accounts(&motion, &step, 0.0, dt);
        close(step.support_force_bound_n, 8.0 * force, force, 1e-12);
        errors.push(
            ((step.contact.indentation_m - expected_d) / length)
                .hypot((step.state().velocity_m_s() - expected_v) / speed),
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
fn attached_equilibrium_retains_distinct_contact_and_machine_reactions() {
    let normal = source();
    let (length, force, _, time) = scales(&normal);
    let k = 2.0 * force / length;
    let anchor = -length;
    let response = normal.evaluate_attached(length).unwrap();
    let support = k * (length - anchor);
    let load = response.normal_force_n + support;
    let motion =
        JkrNormalMotion::new_with_support(&normal, MASS, budget(&normal), k, anchor).unwrap();
    let initial = motion.initialize(response.state(), 0.0).unwrap();
    let dt = 0.5 * time;
    let step = motion
        .advance(initial, load, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(step.normal_increment_m, 0.0);
    assert_eq!(step.state().velocity_m_s(), 0.0);
    assert_eq!(step.normal_work_j, 0.0);
    assert_eq!(step.support_work_j, 0.0);
    close(
        step.normal_impulse_n_s,
        response.normal_force_n * dt,
        force * time,
        1e-13,
    );
    close(step.support_impulse_n_s, support * dt, force * time, 1e-13);
    close(
        step.support_time_moment_n_s2,
        0.5 * support * dt * dt,
        force * time * time,
        1e-13,
    );
    close(
        step.support_energy_j,
        0.5 * k * (length - anchor).powi(2),
        force * length,
        1e-13,
    );
    accounts(&motion, &step, load, dt);
}

#[test]
fn machine_stiffness_and_preload_change_pull_off_without_changing_snap_law() {
    let normal = source();
    let (length, force, _, time) = scales(&normal);
    let dt = 4.0 * time;
    let load = -10.0 * force;
    let initial_contact = normal.evaluate_attached(0.0).unwrap().state();
    let mut branches = Vec::new();
    for (k, anchor) in [
        (0.0, 0.0),
        (4.0 * force / length, 0.0),
        (4.0 * force / length, -4.0 * length),
    ] {
        let motion =
            JkrNormalMotion::new_with_support(&normal, MASS, budget(&normal), k, anchor).unwrap();
        let initial = motion.initialize(initial_contact, 0.0).unwrap();
        let step = motion
            .advance(initial, load, dt, &CancelGate::new())
            .unwrap();
        branches.push(step.contact.branch);
        accounts(&motion, &step, load, dt);
        if step.contact.branch == JkrBranch::Detached {
            assert_eq!(step.event_count(), 1);
            let event = step.events().next().unwrap();
            assert_eq!(
                event.before().indentation_m(),
                normal.minimum_attached_indentation_m()
            );
            close(
                step.released_energy_j,
                normal
                    .evaluate_attached(normal.minimum_attached_indentation_m())
                    .unwrap()
                    .reversible_energy_j,
                force * length,
                1e-12,
            );
        } else {
            assert_eq!(step.event_count(), 0);
            assert_eq!(step.released_energy_j, 0.0);
        }
    }
    assert_eq!(
        branches,
        [
            JkrBranch::Detached,
            JkrBranch::Attached,
            JkrBranch::Detached
        ]
    );
}

#[test]
fn coarse_motion_cannot_tunnel_through_the_supported_adhesion_barrier() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let k = 2.0 * force / length;
    let load = -12.0 * force; // Outside the unsupported law's [-9,-5] force interval.
    let mut lo = normal.minimum_attached_indentation_m();
    let mut hi = -2.8 * length;
    for _ in 0..60 {
        let mid = lo.midpoint(hi);
        if normal.evaluate_attached(mid).unwrap().normal_force_n + k * mid > load {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let barrier = lo.midpoint(hi);
    let from = -2.4 * length;
    let potential = |d: f64| {
        normal.evaluate_attached(d).unwrap().reversible_energy_j + 0.5 * k * d * d - load * d
    };
    let energy_to_barrier = potential(barrier) - potential(from);
    let velocity = -0.5 * speed;
    assert!(0.5 * MASS * velocity * velocity < 0.5 * energy_to_barrier);
    let mut controls = budget(&normal);
    controls.maximum_substep_s = 8.0 * time;
    let motion = JkrNormalMotion::new_with_support(&normal, MASS, controls, k, 0.0).unwrap();
    let initial = motion
        .initialize(normal.evaluate_attached(from).unwrap().state(), velocity)
        .unwrap();
    let step = motion
        .advance(initial, load, 8.0 * time, &CancelGate::new())
        .unwrap();
    assert_eq!(step.event_count(), 0);
    assert_eq!(step.contact.branch, JkrBranch::Attached);
    assert!(step.contact.indentation_m > barrier);
    assert!(step.trials > step.substeps);
    accounts(&motion, &step, load, 8.0 * time);
}

#[test]
fn support_force_budget_includes_interior_reversals() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let mut controls = budget(&normal);
    controls.maximum_substep_s = 4.0 * time;
    controls.maximum_force_n = 2.1 * force;
    let motion =
        JkrNormalMotion::new_with_support(&normal, MASS, controls, force / length, -10.0 * length)
            .unwrap();
    let initial = motion
        .initialize(normal.detached_state(-10.0 * length).unwrap(), 2.0 * speed)
        .unwrap();
    let step = motion
        .advance(initial, 0.0, 4.0 * time, &CancelGate::new())
        .unwrap();
    // The rejected single midpoint arc reaches 2.5 Fs internally, although
    // its endpoint reaction is only 1.6 Fs. Bounded refinement resolves it.
    assert!(step.trials > step.substeps);
    assert!(step.support_force_bound_n <= controls.maximum_force_n);
    assert!(step.support_force_bound_n > 1.9 * force);
    assert_eq!(step.event_count(), 0);
    accounts(&motion, &step, 0.0, 4.0 * time);
}

#[test]
fn support_parameters_are_validated_and_bound_to_the_sealed_history() {
    let normal = source();
    let (length, force, _, time) = scales(&normal);
    let controls = budget(&normal);
    assert!(JkrNormalMotion::new_with_support(&normal, MASS, controls, -1.0, 0.0).is_err());
    assert!(JkrNormalMotion::new_with_support(&normal, MASS, controls, 1.0, f64::NAN).is_err());
    let motion =
        JkrNormalMotion::new_with_support(&normal, MASS, controls, force / length, 0.0).unwrap();
    let initial = motion
        .initialize(normal.evaluate_attached(0.0).unwrap().state(), 0.0)
        .unwrap();
    for (k, anchor) in [(2.0 * force / length, 0.0), (force / length, length)] {
        let other = JkrNormalMotion::new_with_support(&normal, MASS, controls, k, anchor).unwrap();
        assert_eq!(
            other.advance(initial, 0.0, time, &CancelGate::new()),
            Err(JkrNormalMotionError::StateMismatch)
        );
    }
    let free = JkrNormalMotion::new(&normal, MASS, controls).unwrap();
    let explicit = JkrNormalMotion::new_with_support(&normal, MASS, controls, 0.0, 0.0).unwrap();
    let contact = normal.detached_state(-length).unwrap();
    let state = free.initialize(contact, 0.0).unwrap();
    assert_eq!(state, explicit.initialize(contact, 0.0).unwrap());
    assert_eq!(
        free.advance(state, 0.0, time, &CancelGate::new()),
        explicit.advance(state, 0.0, time, &CancelGate::new())
    );
}

//! Passive midpoint viscosity, physical events and optional sealed loss timing.
use super::*;

fn audit(motion: &JkrNormalMotion<'_>, trace: &JkrNormalMotionTrace, load: f64, dt: f64) {
    let step = trace.step();
    let (length, force, speed, _) = scales(motion.normal());
    let initial = step.initial_state();
    let initial_energy = 0.5 * MASS * initial.velocity_m_s().powi(2)
        + motion.response(initial).unwrap().reversible_energy_j
        + motion.support_energy_j(initial).unwrap();
    close(
        step.total_energy_j + step.released_energy_j + step.damping_dissipation_j,
        initial_energy + load * step.normal_increment_m,
        force * length,
        3e-9,
    );
    close(
        MASS * (step.state().velocity_m_s() - initial.velocity_m_s()),
        load * dt - step.reaction_impulse_n_s(),
        MASS * speed,
        3e-10,
    );
    close(
        step.normal_increment_m,
        dt * initial.velocity_m_s()
            + (0.5 * load * dt * dt - dt * step.reaction_impulse_n_s()
                + step.reaction_time_moment_n_s2())
                / MASS,
        length,
        3e-10,
    );
    let (mut work, mut impulse, mut moment, mut previous_end) = (0.0, 0.0, 0.0, 0.0);
    for segment in trace.segments() {
        assert!(segment.start_offset_s() >= previous_end);
        assert!(segment.end_offset_s() <= dt);
        assert!(segment.dissipated_work_j() > 0.0);
        let c = motion.support_damping_n_s_per_m();
        let v = segment.midpoint_velocity_m_s();
        let h = segment.duration_s();
        close(
            segment.dissipated_work_j(),
            c * v * v * h,
            force * length,
            1e-12,
        );
        close(
            segment.power_w() * h,
            segment.dissipated_work_j(),
            force * length,
            1e-12,
        );
        work += segment.dissipated_work_j();
        impulse += c * v * h;
        moment += c * v * h * (segment.start_offset_s() + 0.5 * h);
        previous_end = segment.end_offset_s();
    }
    close(step.damping_dissipation_j, work, force * length, 2e-11);
    close(step.damping_impulse_n_s, impulse, MASS * speed, 2e-11);
    close(step.damping_time_moment_n_s2, moment, MASS * length, 2e-11);
    assert_eq!(
        *step,
        motion
            .advance(initial, load, dt, &CancelGate::new())
            .unwrap()
    );
}

#[test]
fn damped_oscillator_loss_and_motion_converge_through_velocity_reversals() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let (k, c, anchor) = (4.0 * force / length, 0.8 * MASS / time, -10.0 * length);
    let dt = 3.0 * time;
    let decay = c / (2.0 * MASS);
    let frequency = (k / MASS - decay * decay).sqrt();
    let amplitude = 2.0 * length;
    let expected_x = amplitude
        * (-decay * dt).exp()
        * ((frequency * dt).cos() + decay / frequency * (frequency * dt).sin());
    let expected_v =
        -amplitude * (k / MASS) / frequency * (-decay * dt).exp() * (frequency * dt).sin();
    let expected_loss = 0.5 * k * amplitude * amplitude
        - 0.5 * k * expected_x * expected_x
        - 0.5 * MASS * expected_v * expected_v;
    let mut errors = Vec::new();
    for resolution in [0.1, 0.05, 0.025] {
        let mut controls = budget(&normal);
        controls.maximum_substep_s = resolution * time;
        let motion =
            JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, k, anchor, c)
                .unwrap();
        let initial = motion
            .initialize(normal.detached_state(anchor + amplitude).unwrap(), 0.0)
            .unwrap();
        let trace = motion
            .advance_with_dissipation_trace(initial, 0.0, dt, &CancelGate::new())
            .unwrap();
        let step = trace.step();
        assert_eq!(step.event_count(), 0);
        assert!(
            trace
                .segments()
                .iter()
                .any(|s| s.midpoint_velocity_m_s() > 0.0)
        );
        assert!(
            trace
                .segments()
                .iter()
                .any(|s| s.midpoint_velocity_m_s() < 0.0)
        );
        audit(&motion, &trace, 0.0, dt);
        errors.push(
            ((step.contact.indentation_m - anchor - expected_x) / length)
                .hypot((step.state().velocity_m_s() - expected_v) / speed)
                .hypot((step.damping_dissipation_j - expected_loss) / (force * length)),
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
fn pure_viscous_decay_refines_instead_of_reversing_or_crossing_an_unreachable_touch() {
    let normal = source();
    let (length, _, speed, time) = scales(&normal);
    let mut controls = budget(&normal);
    controls.maximum_substep_s = 4.0 * time;
    let c = MASS / time;
    let motion =
        JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, 0.0, 0.0, c).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-1.5 * length).unwrap(), speed)
        .unwrap();
    let trace = motion
        .advance_with_dissipation_trace(initial, 0.0, 4.0 * time, &CancelGate::new())
        .unwrap();
    let step = trace.step();
    assert!(step.trials > step.substeps);
    assert!(step.state().velocity_m_s() > 0.0 && step.state().velocity_m_s() < speed);
    assert_eq!(step.event_count(), 0);
    assert!(step.contact.indentation_m < -0.5 * length);
    close(
        step.normal_increment_m,
        MASS / c * (speed - step.state().velocity_m_s()),
        length,
        2e-13,
    );
    audit(&motion, &trace, 0.0, 4.0 * time);
}

#[test]
fn viscous_touch_time_uses_the_actual_event_velocity_and_leaves_snap_impulse_zero() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let c = MASS / time;
    let distance = 0.25 * length;
    let dt = distance / (speed - 0.5 * c * distance / MASS);
    let mut controls = budget(&normal);
    controls.maximum_substep_s = dt;
    let motion =
        JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, 0.0, 0.0, c).unwrap();
    let initial = motion
        .initialize(normal.detached_state(-distance).unwrap(), speed)
        .unwrap();
    let trace = motion
        .advance_with_dissipation_trace(initial, 0.0, dt, &CancelGate::new())
        .unwrap();
    let step = trace.step();
    assert_eq!(step.event_count(), 1);
    let event = step.events().next().unwrap();
    close(event.offset_s, dt, time, 1e-14);
    close(
        event.velocity_m_s,
        speed - c * distance / MASS,
        speed,
        1e-13,
    );
    // Event localization may leave one ulp of attached smooth time. Its force
    // contributes that tiny impulse; the original snap still preserves v.
    let tail_bound = 9.0 * force * (dt - event.offset_s);
    assert!(step.normal_impulse_n_s.abs() <= tail_bound);
    assert!(step.normal_time_moment_n_s2.abs() <= dt * tail_bound);
    close(
        step.released_energy_j,
        -normal.evaluate_attached(0.0).unwrap().reversible_energy_j,
        force * length,
        1e-12,
    );
    assert!(
        trace
            .segments()
            .iter()
            .any(|s| s.end_offset_s() == event.offset_s)
    );
    audit(&motion, &trace, 0.0, dt);
}

#[test]
fn viscous_loss_changes_capture_and_pull_off_without_becoming_snap_loss() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let dt = 8.0 * time;
    for (c, branch, events) in [
        (0.0, JkrBranch::Detached, 2),
        (MASS / time, JkrBranch::Attached, 1),
    ] {
        let motion =
            JkrNormalMotion::new_with_damped_support(&normal, MASS, budget(&normal), 0.0, 0.0, c)
                .unwrap();
        let initial = motion
            .initialize(normal.detached_state(-0.25 * length).unwrap(), 10.0 * speed)
            .unwrap();
        let trace = motion
            .advance_with_dissipation_trace(initial, 0.0, dt, &CancelGate::new())
            .unwrap();
        let step = trace.step();
        assert_eq!(step.contact.branch, branch);
        assert_eq!(step.event_count(), events);
        let snap_sum: f64 = step.events().map(|e| e.jump.released_energy_j).sum();
        close(step.released_energy_j, snap_sum, force * length, 1e-12);
        assert_eq!(step.damping_dissipation_j > 0.0, c > 0.0);
        for event in step.events() {
            assert!(trace.segments().iter().all(
                |s| s.end_offset_s() <= event.offset_s || s.start_offset_s() >= event.offset_s
            ));
        }
        audit(&motion, &trace, 0.0, dt);
    }
    let motion = JkrNormalMotion::new_with_damped_support(
        &normal,
        MASS,
        budget(&normal),
        4.0 * force / length,
        -4.0 * length,
        MASS / time,
    )
    .unwrap();
    let initial = motion
        .initialize(normal.evaluate_attached(0.0).unwrap().state(), 0.0)
        .unwrap();
    let trace = motion
        .advance_with_dissipation_trace(initial, -10.0 * force, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(trace.step().contact.branch, JkrBranch::Detached);
    assert_eq!(trace.step().event_count(), 1);
    assert!(trace.step().damping_dissipation_j > 0.0);
    audit(&motion, &trace, -10.0 * force, dt);
}

#[test]
fn combined_spring_dashpot_force_includes_its_interior_extremum() {
    let normal = source();
    let (length, force, speed, time) = scales(&normal);
    let (k, c, anchor, dt) = (
        force / length,
        0.2 * MASS / time,
        -10.0 * length,
        4.0 * time,
    );
    let mut controls = budget(&normal);
    controls.maximum_substep_s = dt;
    let motion =
        JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, k, anchor, c).unwrap();
    let initial = motion
        .initialize(normal.detached_state(anchor).unwrap(), 2.0 * speed)
        .unwrap();
    let trace = motion
        .advance_with_dissipation_trace(initial, 0.0, dt, &CancelGate::new())
        .unwrap();
    let step = trace.step();
    assert_eq!(step.substeps, 1);
    let v0 = initial.velocity_m_s();
    let acceleration = (step.state().velocity_m_s() - v0) / dt;
    let at = -v0 / acceleration - c / k;
    assert!(at > 0.0 && at < dt);
    let expected = k * (v0 * at + 0.5 * acceleration * at * at) + c * (v0 + acceleration * at);
    close(step.support_reaction_force_bound_n, expected, force, 1e-12);
    assert!(step.support_reaction_force_bound_n > step.support_force_bound_n);
    close(step.damping_force_bound_n, c * v0, force, 1e-12);
    audit(&motion, &trace, 0.0, dt);
}

#[test]
fn zero_damping_is_exact_and_damped_history_refuses_foreign_policy_and_cancellation() {
    let normal = source();
    let (length, force, _, time) = scales(&normal);
    let controls = budget(&normal);
    let elastic =
        JkrNormalMotion::new_with_support(&normal, MASS, controls, force / length, 0.0).unwrap();
    let zero =
        JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, force / length, 0.0, 0.0)
            .unwrap();
    let state = elastic
        .initialize(normal.evaluate_attached(0.0).unwrap().state(), 0.0)
        .unwrap();
    assert_eq!(state, zero.initialize(state.contact_state(), 0.0).unwrap());
    let trace = zero
        .advance_with_dissipation_trace(state, 0.0, time, &CancelGate::new())
        .unwrap();
    assert!(trace.segments().is_empty());
    assert_eq!(
        *trace.step(),
        elastic
            .advance(state, 0.0, time, &CancelGate::new())
            .unwrap()
    );
    let damped = JkrNormalMotion::new_with_damped_support(
        &normal,
        MASS,
        controls,
        force / length,
        0.0,
        MASS / time,
    )
    .unwrap();
    assert!(matches!(
        damped.advance_with_dissipation_trace(state, 0.0, time, &CancelGate::new()),
        Err(JkrNormalMotionError::StateMismatch)
    ));
    let state = damped.initialize(state.contact_state(), 0.0).unwrap();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        damped.advance_with_dissipation_trace(state, 0.0, time, &gate),
        Err(JkrNormalMotionError::Cancelled)
    ));
    assert!(
        JkrNormalMotion::new_with_damped_support(&normal, MASS, controls, 0.0, 0.0, -1.0).is_err()
    );
}

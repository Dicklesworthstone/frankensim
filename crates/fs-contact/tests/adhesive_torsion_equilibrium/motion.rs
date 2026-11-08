//! G1/G3/G4: actual two-coordinate inertia, source work and whole-path refusal.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::motion::*;

const MASS: f64 = 1e-5;

// Translation and rotation have different units. The fixture chooses the
// inertia so both coordinates evolve on comparable physical time scales.
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let length = source.normal_source().indentation_scale_m();
    let angle = oracle(source).2;
    let force = PI
        * source
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2()
        * source.normal_source().effective_radius_m();
    let time = (MASS * length / force).sqrt();
    (length, angle, force, time)
}

fn mobility(source: &SourceJkrMixedMode) -> [[f64; 2]; 2] {
    let (length, angle, _, _) = scales(source);
    let inverse_inertia = (angle / length).powi(2) / MASS;
    let cross = 0.3 * (inverse_inertia / MASS).sqrt();
    [[1.0 / MASS, cross], [cross, inverse_inertia]]
}

fn budget(source: &SourceJkrMixedMode) -> JkrTorsionalMotionBudget {
    let (length, angle, force, time) = scales(source);
    JkrTorsionalMotionBudget {
        maximum_step_s: 0.5 * time,
        maximum_coordinate_increment: [2.0 * length, 2.0 * angle],
        maximum_rate: [100.0 * length / time, 100.0 * angle / time],
        maximum_effort: [1e4 * force, 1e4 * force * length / angle],
        maximum_iterations: 48,
        maximum_evaluations: 512,
        maximum_backtracks: 24,
        coordinate_absolute_tolerance: [1e-14 * length, 1e-14 * angle],
        coordinate_relative_tolerance: 1e-13,
        rate_absolute_tolerance: [1e-14 * length / time, 1e-14 * angle / time],
        rate_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * length,
        energy_relative_tolerance: 1e-10,
    }
}

fn mul(w: [[f64; 2]; 2], value: [f64; 2]) -> [f64; 2] {
    w.map(|row| row[0] * value[0] + row[1] * value[1])
}

fn force(source: &SourceJkrMixedMode, q: [f64; 2]) -> [f64; 2] {
    let r = source
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap();
    [r.normal_force_n, r.torsional_moment_n_m]
}

// Independent inverse of the full physical mixed-unit mobility, including its
// off-diagonal term. A diagonal effective-mass surrogate would fail this audit.
fn kinetic(w: [[f64; 2]; 2], v: [f64; 2]) -> f64 {
    let determinant = w[0][0] * w[1][1] - w[0][1] * w[1][0];
    0.5 * (w[1][1] * v[0] * v[0] - 2.0 * w[0][1] * v[0] * v[1] + w[0][0] * v[1] * v[1])
        / determinant
}

#[test]
fn g1_full_mobility_retains_source_work_native_impulses_and_quadratic_force_bounds() {
    let source = source();
    let gate = CancelGate::new();
    let (length, angle, force_scale, time) = scales(&source);
    let (_, delta, _) = oracle(&source);
    let w = mobility(&source);
    assert!(w[1][1] / w[0][0] > 1e6);
    let motion = JkrTorsionalMotion::new(&source, w, budget(&source)).unwrap();
    let initial = motion
        .initial_state(
            [delta, 0.6 * angle],
            [0.12 * length / time, -0.2 * angle / time],
            &gate,
        )
        .unwrap();
    let applied = [0.3 * force_scale, -0.15 * force_scale * length / angle];
    let external = mul(w, applied);
    let dt = 0.08 * time;
    let step = motion.advance(initial, dt, external, &gate).unwrap();
    assert_eq!(motion.advance(initial, dt, external, &gate).unwrap(), step);
    assert_eq!(step.initial_state(), initial);
    assert_eq!(step.initial_contact, initial.contact());
    assert_eq!(step.contact, step.state().contact());
    assert_eq!(step.state().steps(), 1);
    assert_eq!(step.state().time_s(), dt);
    assert_eq!(step.contact.model_identity(), source.identity());
    assert!(step.iterations > 0);
    assert!(step.evaluations <= budget(&source).maximum_evaluations);
    let reaction = mul(w, step.mean_effort);
    let dq: [f64; 2] =
        core::array::from_fn(|i| step.state().coordinates()[i] - initial.coordinates()[i]);
    for i in 0..2 {
        let qscale = [length, angle][i];
        work_close(
            dq[i],
            0.5 * dt * (initial.rates()[i] + step.state().rates()[i]),
            qscale,
            2e-12,
        );
        work_close(
            step.state().rates()[i] - initial.rates()[i],
            dt * (external[i] - reaction[i]),
            qscale / time,
            2e-12,
        );
        close(step.effort_impulse[i], step.mean_effort[i] * dt, 2e-15);
        close(
            step.effort_time_moment[i],
            0.5 * dt * step.effort_impulse[i],
            2e-15,
        );
    }
    let energy = force_scale * length;
    work_close(
        step.contact_work_j,
        step.contact.reversible_energy_j - initial.contact().reversible_energy_j,
        energy,
        2e-12,
    );
    work_close(
        step.kinetic_energy_change_j,
        kinetic(w, step.state().rates()) - kinetic(w, initial.rates()),
        energy,
        2e-12,
    );
    work_close(
        step.applied_work_j,
        applied[0] * dq[0] + applied[1] * dq[1],
        energy,
        2e-12,
    );
    work_close(
        step.kinetic_energy_change_j + step.contact_work_j,
        step.applied_work_j,
        energy,
        2e-11,
    );
    assert!(step.energy_residual_j.abs() < 2e-11 * energy);
    assert!(step.source_work_residual_j.abs() < 2e-12 * energy);
    for k in 0..=32 {
        let t = dt * k as f64 / 32.0;
        let q = core::array::from_fn(|i| {
            initial.coordinates()[i]
                + t * initial.rates()[i]
                + 0.5 * t * t * (external[i] - reaction[i])
        });
        for (i, effort) in force(&source, q).into_iter().enumerate() {
            assert!(effort.abs() <= step.maximum_source_effort[i]);
        }
    }
}

// Independent explicit RK4 consults only the pointwise original source. It
// shares no discrete gradient, Jacobian, quadrature or production integrator.
fn rk4(
    source: &SourceJkrMixedMode,
    w: [[f64; 2]; 2],
    initial: [f64; 4],
    external: [f64; 2],
    duration: f64,
) -> [f64; 4] {
    let rhs = |y: [f64; 4]| {
        let reaction = mul(w, force(source, [y[0], y[1]]));
        [
            y[2],
            y[3],
            external[0] - reaction[0],
            external[1] - reaction[1],
        ]
    };
    let plus = |y: [f64; 4], k: [f64; 4], dt: f64| core::array::from_fn(|i| y[i] + dt * k[i]);
    let dt = duration / 512.0;
    let mut y = initial;
    for _ in 0..512 {
        let a = rhs(y);
        let b = rhs(plus(y, a, 0.5 * dt));
        let c = rhs(plus(y, b, 0.5 * dt));
        let d = rhs(plus(y, c, dt));
        y = core::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn g1_coupled_recoil_and_radius_relaxation_converge_to_original_force_ode() {
    let source = source();
    let gate = CancelGate::new();
    let (length, angle, force_scale, time) = scales(&source);
    let (_, delta, _) = oracle(&source);
    let initial = [
        delta,
        0.6 * angle,
        0.12 * length / time,
        0.25 * angle / time,
    ];
    let q = [initial[0], initial[1]];
    let w = mobility(&source);
    let external = mul(w, [force(&source, q)[0], 0.0]);
    let duration = 0.4 * time;
    let reference = rk4(&source, w, initial, external, duration);
    let motion = JkrTorsionalMotion::new(&source, w, budget(&source)).unwrap();
    let start = motion
        .initial_state(q, [initial[2], initial[3]], &gate)
        .unwrap();
    let run = |count: usize| {
        let mut state = start;
        for _ in 0..count {
            state = motion
                .advance(state, duration / count as f64, external, &gate)
                .unwrap()
                .state();
        }
        state
    };
    let coarse = run(4);
    let fine = run(8);
    let error = |state: JkrTorsionalMotionState| -> f64 {
        (0..2)
            .map(|i| {
                (state.coordinates()[i] - reference[i]).abs() / [length, angle][i]
                    + time * (state.rates()[i] - reference[i + 2]).abs() / [length, angle][i]
            })
            .sum()
    };
    assert!(
        error(coarse) > 3.0 * error(fine),
        "coarse={}, fine={}",
        error(coarse),
        error(fine)
    );
    assert!(error(fine) < 2e-3);
    assert_eq!(fine.steps(), 8);
    close(fine.time_s(), duration, 2e-15);
    assert!(
        (fine.contact().contact_radius_m - start.contact().contact_radius_m).abs()
            > 1e-5 * start.contact().contact_radius_m
    );
    assert!(
        (fine.contact().normal_force_n - start.contact().normal_force_n).abs() > 1e-5 * force_scale
    );
    // The held-force reversible phase must also retrace under time reversal.
    let reversed = motion
        .initial_state(fine.coordinates(), fine.rates().map(|v| -v), &gate)
        .unwrap();
    let mut returned = reversed;
    for _ in 0..8 {
        returned = motion
            .advance(returned, duration / 8.0, external, &gate)
            .unwrap()
            .state();
    }
    for i in 0..2 {
        work_close(
            returned.coordinates()[i],
            start.coordinates()[i],
            [length, angle][i],
            2e-10,
        );
        work_close(
            returned.rates()[i],
            -start.rates()[i],
            [length, angle][i] / time,
            2e-10,
        );
    }
}

#[test]
fn g3_hidden_twist_turning_point_refuses_even_when_both_endpoints_are_identical() {
    let source = source();
    let gate = CancelGate::new();
    let (_, delta, _) = oracle(&source);
    let (_, _, _, time) = scales(&source);
    let limit = source
        .torsional_displacement_limit(delta)
        .unwrap()
        .unwrap()
        .twist_angle_rad;
    let q = [delta, 0.0];
    source.torsional_discrete_path(q, q, &gate).unwrap();
    assert!(
        source
            .evaluate_torsion_displacement_controlled(delta, 1.05 * limit)
            .is_err()
    );
    let mut bounds = budget(&source);
    bounds.maximum_coordinate_increment[1] = 2.0 * limit;
    bounds.maximum_rate[1] = 100.0 * limit / time;
    bounds.maximum_effort = bounds.maximum_effort.map(|f| f * 100.0);
    let motion = JkrTorsionalMotion::new(&source, mobility(&source), bounds).unwrap();
    let dt = 0.2 * time;
    let rates = [0.0, 4.0 * 1.05 * limit / dt];
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let reaction = mul(mobility(&source), force(&source, q));
    // Exactly returning q1=q0 solves midpoint inertia: v1=-v0. The actual
    // parabola reaches 1.05*limit at h/2, invisible to a straight endpoint path.
    let external = core::array::from_fn(|i| reaction[i] - 2.0 * rates[i] / dt);
    assert!(matches!(
        motion.advance(initial, dt, external, &gate),
        Err(JkrTorsionalMotionError::Path(
            JkrTorsionalPathError::UnresolvedEnvelope
        ))
    ));
    assert_eq!(initial.coordinates(), q);
    assert_eq!(initial.rates(), rates);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.steps(), 0);
}

#[test]
fn g3_nonreciprocal_unresolved_or_ambiguous_inertia_is_not_a_surrogate_motion() {
    let source = source();
    let gate = CancelGate::new();
    let b = budget(&source);
    for w in [
        [[1.0, 0.1], [0.2, 1.0]],
        [[1.0, 1.0], [1.0, 1.0]],
        [[1.0, 1.1], [1.1, 1.0]],
        [[0.0, 0.0], [0.0, 1.0]],
    ] {
        assert!(matches!(
            JkrTorsionalMotion::new(&source, w, b),
            Err(JkrTorsionalMotionError::Invalid(_))
        ));
    }
    // Validity is checked in a unit-normalized metric, not with a threshold
    // that compares a translational diagonal directly to a rotational one.
    let extreme = [[1e-120, 0.2], [0.2, 1e120]];
    assert_eq!(
        JkrTorsionalMotion::new(&source, extreme, b)
            .unwrap()
            .mobility(),
        extreme
    );
    let q = [
        0.99 * source.normal_source().minimum_attached_indentation_m(),
        0.0,
    ];
    let contact = source
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap();
    let stiffness = contact.coupled_stiffness.unwrap()[0][0];
    assert!(stiffness < 0.0);
    let mut w = mobility(&source);
    w[0][1] = 0.0;
    w[1][0] = 0.0;
    let dt = 4.0 / (-w[0][0] * stiffness).sqrt();
    let mut b = b;
    b.maximum_step_s = dt;
    let motion = JkrTorsionalMotion::new(&source, w, b).unwrap();
    let initial = motion.initial_state(q, [0.0; 2], &gate).unwrap();
    assert_eq!(
        motion.advance(initial, dt, mul(w, force(&source, q)), &gate),
        Err(JkrTorsionalMotionError::InertialTangent)
    );
}

#[test]
fn g4_budget_cancel_and_owner_mismatch_leave_the_original_candidate_retryable() {
    let source = source();
    let gate = CancelGate::new();
    let (length, angle, _, time) = scales(&source);
    let (_, delta, _) = oracle(&source);
    let q = [delta, 0.5 * angle];
    let rates = [0.1 * length / time, 0.2 * angle / time];
    let w = mobility(&source);
    let external = mul(w, force(&source, q));
    let motion = JkrTorsionalMotion::new(&source, w, budget(&source)).unwrap();
    let state = motion.initial_state(q, rates, &gate).unwrap();
    let successful = motion.advance(state, 0.05 * time, external, &gate).unwrap();
    let mut one_query = budget(&source);
    one_query.maximum_evaluations = 1;
    let limited = JkrTorsionalMotion::new(&source, w, one_query).unwrap();
    let limited_state = limited.initial_state(q, rates, &gate).unwrap();
    assert_eq!(
        limited.advance(limited_state, 0.05 * time, external, &gate),
        Err(JkrTorsionalMotionError::Budget)
    );
    assert_eq!(
        limited.advance(state, 0.05 * time, external, &gate),
        Err(JkrTorsionalMotionError::StateMismatch)
    );
    assert_eq!(
        motion.advance(state, 0.05 * time, external, &gate).unwrap(),
        successful
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        motion.advance(state, 0.05 * time, external, &cancelled),
        Err(JkrTorsionalMotionError::Cancelled)
    );
    assert_eq!(state.coordinates(), q);
    assert_eq!(state.rates(), rates);
    assert_eq!(state.steps(), 0);
}

#[test]
fn g1_projection_handoff_preserves_original_endpoint_and_clock_without_work() {
    let source = source();
    let gate = CancelGate::new();
    let (length, angle, force_scale, time) = scales(&source);
    let (_, delta, _) = oracle(&source);
    let w = mobility(&source);
    let motion = JkrTorsionalMotion::new(&source, w, budget(&source)).unwrap();
    let initial = motion
        .initial_state(
            [delta, 0.6 * angle],
            [0.1 * length / time, -0.2 * angle / time],
            &gate,
        )
        .unwrap();
    let step = motion
        .advance(initial, 0.05 * time, [0.0; 2], &gate)
        .unwrap();
    let previous = step.state();
    let changed = w.map(|row| row.map(|v| 1.4 * v));
    let rates = [0.2 * length / time, 0.3 * angle / time];
    let transfer = motion
        .rebind_projection(previous, changed, rates, &gate)
        .unwrap();
    assert_eq!(transfer.previous_state(), previous);
    assert_eq!(transfer.state().coordinates(), previous.coordinates());
    assert_eq!(transfer.state().contact(), previous.contact());
    assert_eq!(transfer.state().time_s(), previous.time_s());
    assert_eq!(transfer.state().steps(), previous.steps());
    assert_eq!(transfer.state().rates(), rates);
    assert_eq!(transfer.motion().source().identity(), source.identity());
    let energy = force_scale * length;
    work_close(
        transfer.report().projected_kinetic_energy_before_j,
        kinetic(w, previous.rates()),
        energy,
        2e-13,
    );
    work_close(
        transfer.report().projected_kinetic_energy_after_j,
        kinetic(changed, rates),
        energy,
        2e-13,
    );
    work_close(
        transfer.report().projected_kinetic_repartition_j,
        kinetic(changed, rates) - kinetic(w, previous.rates()),
        energy,
        2e-13,
    );
    let (next, state, _) = transfer.into_parts();
    assert_eq!(
        next.advance(previous, 0.01 * time, [0.0; 2], &gate),
        Err(JkrTorsionalMotionError::StateMismatch)
    );
    let resumed = next.advance(state, 0.01 * time, [0.0; 2], &gate).unwrap();
    assert_eq!(resumed.initial_state(), state);
    assert_eq!(resumed.state().steps(), previous.steps() + 1);
    close(
        resumed.state().time_s(),
        previous.time_s() + 0.01 * time,
        2e-15,
    );
}

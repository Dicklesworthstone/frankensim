//! Objective I/II/III history and independent mixed-unit kinetic observations.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::combined::{
    JkrCombinedResponse, motion::transport::JkrCombinedProjectionReport,
};
use fs_tribo::{partial_slip::TangentFrame, tangential::transport::TangentialHistoryRotation};

fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}

fn identity() -> TangentialHistoryRotation {
    TangentialHistoryRotation::between(frame(), frame(), 0.0).unwrap()
}

fn projection(rotation: TangentialHistoryRotation) -> [[f64; 4]; 4] {
    let [c, s] = rotation.cosine_sine();
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, c, -s, 0.0],
        [0.0, s, c, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn congruence(p: [[f64; 4]; 4], w: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut result = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            let value = (0..4)
                .flat_map(|k| (0..4).map(move |l| p[i][k] * w[k][l] * p[j][l]))
                .sum();
            result[i][j] = value;
            result[j][i] = value;
        }
    }
    result
}

fn dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

fn storage(response: JkrCombinedResponse) -> [f64; 4] {
    [
        response.normal_elastic_energy_j,
        response.tangential_elastic_energy_j,
        response.torsional_elastic_energy_j,
        response.surface_energy_j,
    ]
}

fn check_storage(report: &JkrCombinedProjectionReport, energy: f64) {
    let before = storage(report.before);
    let after = storage(report.after);
    for i in 0..4 {
        close(
            report.storage_component_residual_j[i],
            after[i] - before[i],
            energy,
            4e-15,
        );
        close(after[i], before[i], energy, 4e-12);
    }
    close(
        report.absolute_storage_residual_j,
        (0..4).map(|i| (after[i] - before[i]).abs()).sum(),
        energy,
        4e-15,
    );
    close(
        report.source_energy_residual_j,
        report.after.reversible_energy_j - report.before.reversible_energy_j,
        energy,
        4e-15,
    );
    assert!(report.source_energy_residual_j.abs() <= report.source_energy_tolerance_j);
    assert!(report.absolute_storage_residual_j <= report.source_energy_tolerance_j);
    close(
        report.after.contact_area_m2,
        report.before.contact_area_m2,
        report.before.contact_area_m2,
        4e-12,
    );
}

#[test]
fn warm_gauge_rebind_preserves_full_metric_work_and_following_coupled_motion() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (scale, force, time) = scales(&source);
    let energy = force * scale[0];
    let w = mobility(&source);
    let mut policy = budget(&source);
    policy.maximum_coordinate_increment = scale.map(|x| 0.01 * x);
    let motion = JkrCombinedMotion::new(&source, w, policy).unwrap();
    let (q, rates) = start(&source);
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let warm = motion
        .advance(initial, 0.01 * time, [0.0; 4], &gate)
        .unwrap()
        .state();
    assert_eq!(warm.steps(), 1);
    assert!(warm.time_s() > 0.0);
    assert!(warm.coordinates().into_iter().all(|x| x != 0.0));

    let unchanged = motion
        .rebind_projection(warm, identity(), w, warm.rates(), &gate)
        .unwrap();
    assert_eq!(unchanged.previous_state(), warm);
    assert_eq!(unchanged.state(), warm);
    assert_eq!(unchanged.report().before, unchanged.report().after);
    assert_eq!(unchanged.report().storage_component_residual_j, [0.0; 4]);
    assert_eq!(unchanged.report().source_energy_residual_j, 0.0);
    assert_eq!(unchanged.report().absolute_storage_residual_j, 0.0);
    assert_eq!(unchanged.report().projected_kinetic_repartition_j, 0.0);
    assert_eq!(unchanged.report().evaluations, 0);

    let rotation =
        TangentialHistoryRotation::between(frame(), frame().rotated(0.713).unwrap(), 0.0).unwrap();
    let p = projection(rotation);
    let rotated_w = congruence(p, w);
    let transfer = motion
        .rebind_projection(warm, rotation, rotated_w, mul(p, warm.rates()), &gate)
        .unwrap();
    let state = transfer.state();
    let report = transfer.report();
    assert_eq!(transfer.previous_state(), warm);
    assert_eq!(report.before, warm.contact());
    assert_eq!(report.after, state.contact());
    assert_eq!(report.rotation, rotation);
    assert_eq!(report.evaluations, 1);
    assert_eq!(state.coordinates()[0], warm.coordinates()[0]);
    assert_eq!(state.coordinates()[3], warm.coordinates()[3]);
    assert_eq!(state.time_s(), warm.time_s());
    assert_eq!(state.steps(), warm.steps());
    assert_eq!(transfer.motion().budget(), policy);
    assert_eq!(transfer.motion().mobility(), rotated_w);
    assert!(core::ptr::eq(transfer.motion().source(), &source));
    assert_eq!(state.contact().model_identity(), source.identity());
    assert_eq!(
        state.contact(),
        source
            .evaluate_combined_displacement_controlled(state.coordinates())
            .unwrap()
    );
    // Re-expressing a retained vector is not a physical coordinate increment.
    assert!(
        (state.coordinates()[1] - warm.coordinates()[1]).abs()
            > policy.maximum_coordinate_increment[1]
    );
    let expected_q = mul(p, warm.coordinates());
    let expected_effort = mul(p, warm.contact().force);
    let expected_stiffness = congruence(p, warm.contact().coupled_stiffness.unwrap());
    let stiffness = state.contact().coupled_stiffness.unwrap();
    for i in 0..4 {
        close(state.coordinates()[i], expected_q[i], scale[i], 3e-15);
        close(
            state.contact().force[i],
            expected_effort[i],
            energy / scale[i],
            4e-12,
        );
        for j in 0..4 {
            close(
                stiffness[i][j],
                expected_stiffness[i][j],
                energy / (scale[i] * scale[j]),
                4e-12,
            );
        }
    }
    let virtual_increment = core::array::from_fn(|i| [0.2, -0.3, 0.4, 0.5][i] * scale[i]);
    close(
        dot(state.contact().force, mul(p, virtual_increment)),
        dot(warm.contact().force, virtual_increment),
        energy,
        4e-12,
    );
    check_storage(report, energy);
    let expected_kinetic = kinetic(&source, warm.rates());
    close(
        report.projected_kinetic_energy_before_j,
        expected_kinetic,
        energy,
        2e-14,
    );
    close(
        report.projected_kinetic_energy_after_j,
        expected_kinetic,
        energy,
        2e-14,
    );
    close(report.projected_kinetic_repartition_j, 0.0, energy, 2e-14);

    let (rebound, state, _) = transfer.into_parts();
    let applied = core::array::from_fn(|i| [0.1, -0.08, 0.03, 0.04][i] * energy / scale[i]);
    let external = mul(w, applied);
    let reference = motion.advance(warm, 0.01 * time, external, &gate).unwrap();
    let next = rebound
        .advance(state, 0.01 * time, mul(p, external), &gate)
        .unwrap();
    assert_eq!(next.initial_state(), state);
    assert_eq!(next.state().steps(), warm.steps() + 1);
    assert_eq!(next.state().time_s(), reference.state().time_s());
    let expected_q = mul(p, reference.state().coordinates());
    let expected_rates = mul(p, reference.state().rates());
    let expected_mean = mul(p, reference.mean_effort);
    for i in 0..4 {
        close(
            next.state().coordinates()[i],
            expected_q[i],
            scale[i],
            3e-10,
        );
        close(
            next.state().rates()[i],
            expected_rates[i],
            scale[i] / time,
            3e-10,
        );
        close(
            next.mean_effort[i],
            expected_mean[i],
            energy / scale[i],
            3e-10,
        );
    }
    close(next.contact_work_j, reference.contact_work_j, energy, 3e-10);
    close(next.applied_work_j, reference.applied_work_j, energy, 3e-10);
    close(
        next.kinetic_energy_change_j,
        reference.kinetic_energy_change_j,
        energy,
        3e-10,
    );
    assert!(matches!(
        motion.advance(state, 0.01 * time, external, &gate),
        Err(JkrCombinedMotionError::StateMismatch)
    ));
}

#[test]
fn changed_native_observation_repartitions_fixed_body_energy_without_source_work() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (scale, force, time) = scales(&source);
    let energy = force * scale[0];
    let d = metric_scale(&source);
    let diagonal =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] * d[i] } else { 0.0 }));
    let motion = JkrCombinedMotion::new(&source, diagonal, budget(&source)).unwrap();
    let (q, rates) = start(&source);
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let warm = motion
        .advance(initial, 0.005 * time, [0.0; 4], &gate)
        .unwrap()
        .state();

    // Five native mass-normalized velocities are unchanged. The old contact
    // projection is D[I|0], and the new observation is D[I|b]. Twist has its
    // own mixed-unit scale; no translational mass is reused as an inertia.
    let b = [0.4, -0.8, 0.3, 0.6];
    let unseen = 2.0 * energy.sqrt();
    let old_modes = core::array::from_fn(|i| warm.rates()[i] / d[i]);
    let new_modes = core::array::from_fn(|i| old_modes[i] + b[i] * unseen);
    let new_rates = core::array::from_fn(|i| d[i] * new_modes[i]);
    let new_w = core::array::from_fn(|i| {
        core::array::from_fn(|j| d[i] * d[j] * ((if i == j { 1.0 } else { 0.0 }) + b[i] * b[j]))
    });
    let native_energy = 0.5 * (dot(old_modes, old_modes) + unseen * unseen);
    let expected_old = 0.5 * dot(old_modes, old_modes);
    // Sherman-Morrison inverse of I+b b^T, independent of production factors.
    let expected_new =
        0.5 * (dot(new_modes, new_modes) - dot(b, new_modes).powi(2) / (1.0 + dot(b, b)));
    assert!(expected_old < native_energy && expected_new < native_energy);
    assert!((expected_new - expected_old).abs() > 0.1 * energy);
    let transfer = motion
        .rebind_projection(warm, identity(), new_w, new_rates, &gate)
        .unwrap();
    let report = transfer.report();
    close(
        report.projected_kinetic_energy_before_j,
        expected_old,
        native_energy,
        4e-15,
    );
    close(
        report.projected_kinetic_energy_after_j,
        expected_new,
        native_energy,
        4e-15,
    );
    close(
        report.projected_kinetic_repartition_j,
        expected_new - expected_old,
        native_energy,
        4e-15,
    );
    assert_eq!(report.before, report.after);
    assert_eq!(report.source_energy_residual_j, 0.0);
    assert_eq!(report.storage_component_residual_j, [0.0; 4]);
    assert_eq!(report.absolute_storage_residual_j, 0.0);
    assert_eq!(report.evaluations, 0);
    assert_eq!(transfer.state().coordinates(), warm.coordinates());
    assert_eq!(transfer.state().rates(), new_rates);
    assert_eq!(transfer.state().time_s(), warm.time_s());
    assert_eq!(transfer.state().steps(), warm.steps());
    assert_eq!(transfer.motion().mobility(), new_w);
    let (rebound, state, _) = transfer.into_parts();
    let next = rebound
        .advance(state, 0.002 * time, [0.0; 4], &gate)
        .unwrap();
    assert_eq!(next.initial_contact, warm.contact());
    assert_eq!(next.state().steps(), warm.steps() + 1);
    assert!(next.state().time_s() > warm.time_s());
}

#[test]
fn normal_rotation_and_common_spin_preserve_relative_twist_without_a_shear_chord() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (scale, force, time) = scales(&source);
    let energy = force * scale[0];
    let motion = JkrCombinedMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let (q, rates) = start(&source);
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let warm = motion
        .advance(initial, 0.01 * time, [0.0; 4], &gate)
        .unwrap()
        .state();
    let old = warm.coordinates();
    let tilted = TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    // Minimum z->x rotation sends old x->-z and old y->y. A subsequent
    // common quarter-turn about x sends old x->y and old y->z.
    for (current, spin, expected_world) in [
        (frame(), 0.5 * PI, [-old[2], old[1], 0.0]),
        (tilted, 0.0, [0.0, old[2], -old[1]]),
        (tilted, 0.5 * PI, [0.0, old[1], old[2]]),
    ] {
        let rotation = TangentialHistoryRotation::between(frame(), current, spin).unwrap();
        let p = projection(rotation);
        let transfer = motion
            .rebind_projection(
                warm,
                rotation,
                congruence(p, motion.mobility()),
                mul(p, warm.rates()),
                &gate,
            )
            .unwrap();
        let state = transfer.state();
        let u = state.coordinates();
        let world: [f64; 3] = core::array::from_fn(|i| {
            u[1] * current.longitudinal()[i] + u[2] * current.lateral()[i]
        });
        for i in 0..3 {
            close(world[i], expected_world[i], scale[1], 4e-15);
        }
        assert_eq!(u[0], old[0]);
        assert_eq!(u[3], old[3]);
        assert_eq!(state.time_s(), warm.time_s());
        assert_eq!(state.steps(), warm.steps());
        assert_eq!(transfer.report().evaluations, usize::from(u != old));
        check_storage(transfer.report(), energy);
        close(
            transfer.report().projected_kinetic_repartition_j,
            0.0,
            energy,
            2e-14,
        );
        if current == frame() {
            // A linear chord would partially unload shear, enlarge the shared
            // disk, and change the other stores. Transport admits only the
            // isometric endpoint and leaves physical time and twist untouched.
            let middle = core::array::from_fn(|i| old[i].midpoint(u[i]));
            let chord = source
                .evaluate_combined_displacement_controlled(middle)
                .unwrap();
            assert!(chord.contact_area_m2 > warm.contact().contact_area_m2);
            assert!(
                chord.tangential_elastic_energy_j
                    < 0.75 * warm.contact().tangential_elastic_energy_j
            );
            assert_eq!(transfer.report().evaluations, 1);
        }
        let reverse = TangentialHistoryRotation::between(current, frame(), -spin).unwrap();
        let restored = transfer
            .motion()
            .rebind_projection(
                state,
                reverse,
                motion.mobility(),
                mul(projection(reverse), state.rates()),
                &gate,
            )
            .unwrap();
        for i in 0..4 {
            close(restored.state().coordinates()[i], old[i], scale[i], 5e-15);
            close(
                restored.state().rates()[i],
                warm.rates()[i],
                scale[i] / time,
                5e-15,
            );
        }
        assert_eq!(restored.state().time_s(), warm.time_s());
        assert_eq!(restored.state().steps(), warm.steps());
    }
    let antipodal = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
    assert!(TangentialHistoryRotation::between(frame(), antipodal, 0.0).is_err());
}

#[test]
fn rates_metric_ownership_effort_storage_and_cancellation_refuse_without_losing_retry() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let w = mobility(&source);
    let (q, rates) = start(&source);
    let rotation = TangentialHistoryRotation::between(frame(), frame(), 0.5 * PI).unwrap();
    let mut policy = budget(&source);
    policy.maximum_evaluations = 1;
    let motion = JkrCombinedMotion::new(&source, w, policy).unwrap();
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let expected = motion
        .rebind_projection(initial, rotation, w, rates, &gate)
        .unwrap();
    assert_eq!(expected.report().evaluations, 1);
    assert_eq!(expected.motion().budget(), policy);

    for i in 0..4 {
        for value in [f64::NAN, 1.01 * policy.maximum_rate[i]] {
            let mut invalid = rates;
            invalid[i] = value;
            assert!(matches!(
                motion.rebind_projection(initial, rotation, w, invalid, &gate),
                Err(JkrCombinedMotionError::Range(_))
            ));
        }
    }
    let mut zero_twist_mobility = w;
    zero_twist_mobility[3][3] = 0.0;
    let mut nonreciprocal = w;
    nonreciprocal[0][3] *= 1.1;
    let mut indefinite = w;
    indefinite[2][3] = 2.0 * w[2][2].sqrt() * w[3][3].sqrt();
    indefinite[3][2] = indefinite[2][3];
    let mut nonfinite = w;
    nonfinite[3][0] = f64::INFINITY;
    for invalid in [zero_twist_mobility, nonreciprocal, indefinite, nonfinite] {
        assert!(matches!(
            motion.rebind_projection(initial, rotation, invalid, rates, &gate),
            Err(JkrCombinedMotionError::Invalid(_))
        ));
    }
    let mut other_policy = policy;
    other_policy.maximum_iterations += 1;
    let foreign = JkrCombinedMotion::new(&source, w, other_policy).unwrap();
    assert!(matches!(
        foreign.rebind_projection(initial, rotation, w, rates, &gate),
        Err(JkrCombinedMotionError::StateMismatch)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.rebind_projection(initial, rotation, w, rates, &cancelled),
        Err(JkrCombinedMotionError::Cancelled)
    ));
    let retried = motion
        .rebind_projection(initial, rotation, w, rates, &gate)
        .unwrap();
    assert_eq!(retried.previous_state(), initial);
    assert_eq!(retried.state(), expected.state());
    assert_eq!(retried.report(), expected.report());
    assert_eq!(
        motion
            .rebind_projection(initial, identity(), w, rates, &gate)
            .unwrap()
            .state(),
        initial
    );

    // Rotating an accepted force beyond a component limit must be refused;
    // source attachment and rate admission cannot substitute for that limit.
    let mut effort_policy = policy;
    effort_policy.maximum_effort[2] = 1.1 * initial.contact().force[2].abs();
    let limited = JkrCombinedMotion::new(&source, w, effort_policy).unwrap();
    let limited_state = limited.initial_state(q, rates, &gate).unwrap();
    assert!(matches!(
        limited.rebind_projection(limited_state, rotation, w, rates, &gate),
        Err(JkrCombinedMotionError::Range(_))
    ));
    assert_eq!(
        limited
            .rebind_projection(limited_state, identity(), w, rates, &gate)
            .unwrap()
            .state(),
        limited_state
    );

    // Exact-zero allowance must respect defects of individual original
    // stores, independently of whether their rounded total changes.
    let mut strict_policy = policy;
    strict_policy.energy_absolute_tolerance_j = 0.0;
    strict_policy.energy_relative_tolerance = 0.0;
    let strict = JkrCombinedMotion::new(&source, w, strict_policy).unwrap();
    let strict_state = strict.initial_state(q, rates, &gate).unwrap();
    let mut observed_component_defect = false;
    for angle in [0.713, 0.371, 1.137, -2.281, 0.17, 1.89, 2.73, -0.93] {
        let r = TangentialHistoryRotation::between(frame(), frame(), angle).unwrap();
        let shear = r.apply([q[1], q[2]]).unwrap();
        let after = source
            .evaluate_combined_displacement_controlled([q[0], shear[0], shear[1], q[3]])
            .unwrap();
        if storage(after) != storage(strict_state.contact()) {
            observed_component_defect = true;
            assert!(matches!(
                strict.rebind_projection(strict_state, r, w, rates, &gate),
                Err(JkrCombinedMotionError::EnergyBalance {
                    tolerance_j: 0.0,
                    ..
                })
            ));
            break;
        }
    }
    assert!(
        observed_component_defect,
        "fixture must resolve original storage roundoff"
    );
    assert_eq!(
        strict
            .rebind_projection(strict_state, identity(), w, rates, &gate)
            .unwrap()
            .state(),
        strict_state
    );
}

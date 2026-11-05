//! Warm objective history and independently checked native kinetic partition.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::{motion::*, path::JkrMixedModePathError};
use fs_exec::CancelGate;
use fs_tribo::{partial_slip::TangentFrame, tangential::transport::TangentialHistoryRotation};

const MASS: f64 = 1e-5;
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let h = source.normal_source().indentation_scale_m();
    let force = PI * 0.05 * 0.01;
    let speed = (force * h / MASS).sqrt();
    (h, force, speed, h / speed)
}
fn budget(source: &SourceJkrMixedMode) -> JkrMixedModeMotionBudget {
    let (h, force, speed, time) = scales(source);
    JkrMixedModeMotionBudget {
        maximum_step_s: 0.5 * time,
        maximum_coordinate_increment_m: 2.0 * h,
        maximum_speed_m_s: 100.0 * speed,
        maximum_force_n: 1e4 * force,
        maximum_iterations: 48,
        maximum_evaluations: 512,
        maximum_backtracks: 24,
        coordinate_absolute_tolerance_m: 1e-14 * h,
        coordinate_relative_tolerance: 1e-13,
        velocity_absolute_tolerance_m_s: 1e-14 * speed,
        velocity_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * h,
        energy_relative_tolerance: 1e-10,
    }
}
fn mobility() -> [[f64; 3]; 3] {
    [[1.0, 0.2, -0.15], [0.2, 0.8, 0.12], [-0.15, 0.12, 1.1]].map(|row| row.map(|v| v / MASS))
}
fn mul(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    a.map(|row| (0..3).map(|i| row[i] * x[i]).sum())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn projection(rotation: TangentialHistoryRotation) -> [[f64; 3]; 3] {
    let [c, s] = rotation.cosine_sine();
    [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
}
fn congruence(p: [[f64; 3]; 3], w: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .flat_map(|k| (0..3).map(move |l| p[i][k] * w[k][l] * p[j][l]))
                .sum()
        })
    })
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn identity() -> TangentialHistoryRotation {
    TangentialHistoryRotation::between(frame(), frame(), 0.0).unwrap()
}

#[test]
fn warm_gauge_transfer_preserves_chronology_and_following_full_mobility_dynamics() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, force, speed, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let mut policy = budget(&source);
    policy.maximum_coordinate_increment_m = 0.05 * h;
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), policy).unwrap();
    let initial = motion
        .initialize(
            [2.0 * h, 0.3 * h, -0.2 * h],
            [0.1, 0.4, -0.3].map(|v| v * speed),
            &gate,
        )
        .unwrap();
    let warm = motion
        .advance(initial, [0.0; 3], 0.01 * time, &gate)
        .unwrap()
        .state();
    assert_eq!(warm.steps(), 1);
    assert!(warm.time_s() > 0.0);

    let unchanged = motion
        .rebind_native_projection(warm, identity(), mobility(), warm.velocity_m_s(), &gate)
        .unwrap();
    assert_eq!(unchanged.state(), warm);
    assert_eq!(unchanged.previous_state(), warm);
    assert_eq!(unchanged.report().before, unchanged.report().after);
    assert_eq!(unchanged.report().source_energy_residual_j, 0.0);
    assert_eq!(unchanged.report().projected_kinetic_repartition_j, 0.0);
    assert_eq!(unchanged.report().evaluations, 1);

    let rotation =
        TangentialHistoryRotation::between(frame(), frame().rotated(0.713).unwrap(), 0.0).unwrap();
    let p = projection(rotation);
    let rotated_w = congruence(p, mobility());
    let transfer = motion
        .rebind_native_projection(
            warm,
            rotation,
            rotated_w,
            mul(p, warm.velocity_m_s()),
            &gate,
        )
        .unwrap();
    let transported = transfer.state();
    assert_eq!(transported.position_m()[0], warm.position_m()[0]);
    assert_eq!(transported.time_s(), warm.time_s());
    assert_eq!(transported.steps(), warm.steps());
    assert_eq!(transfer.motion().budget(), policy);
    assert!(core::ptr::eq(transfer.motion().source(), &source));
    assert!(core::ptr::eq(transfer.motion().shear(), &shear));
    assert_eq!(transfer.motion().mobility_per_kg(), rotated_w);
    assert_eq!(transfer.report().evaluations, 2);
    // A large change of components is a gauge change, not a physical increment.
    assert!(
        (transported.position_m()[1] - warm.position_m()[1]).abs()
            > policy.maximum_coordinate_increment_m
    );
    for i in 0..3 {
        close(
            transported.position_m()[i],
            mul(p, warm.position_m())[i],
            h,
            2e-15,
        );
    }
    close(
        transfer.report().source_energy_residual_j,
        0.0,
        force * h,
        2e-13,
    );
    close(
        transfer.report().projected_kinetic_repartition_j,
        0.0,
        force * h,
        2e-14,
    );
    let old_t = transfer.report().before.tangential_force_n;
    let new_t = rotation.apply(old_t).unwrap();
    for i in 0..2 {
        close(
            transfer.report().after.tangential_force_n[i],
            new_t[i],
            force,
            2e-13,
        );
    }

    let (mut rebound, state, _) = transfer.into_parts();
    let external = [0.2, -0.15, 0.1].map(|v| v * force / MASS);
    let reference = motion.advance(warm, external, 0.01 * time, &gate).unwrap();
    let next = rebound
        .advance(state, mul(p, external), 0.01 * time, &gate)
        .unwrap();
    assert_eq!(next.state().steps(), 2);
    assert_eq!(next.state().time_s(), reference.state().time_s());
    for i in 0..3 {
        close(
            next.state().position_m()[i],
            mul(p, reference.state().position_m())[i],
            h,
            3e-11,
        );
        close(
            next.state().velocity_m_s()[i],
            mul(p, reference.state().velocity_m_s())[i],
            speed,
            3e-11,
        );
    }
    close(
        next.contact_work_j,
        reference.contact_work_j,
        force * h,
        3e-11,
    );
    close(
        next.applied_work_j,
        reference.applied_work_j,
        force * h,
        3e-11,
    );
    assert!(matches!(
        motion.advance(state, external, 0.01 * time, &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));
}

#[test]
fn changed_observation_repartitions_fixed_native_kinetic_energy_without_work() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, force, speed, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let diagonal =
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]].map(|row| row.map(|v| v / MASS));
    let mut motion = JkrMixedModeMotion::new(&source, &shear, diagonal, budget(&source)).unwrap();
    let initial = motion
        .initialize(
            [2.0 * h, 0.3 * h, -0.2 * h],
            [0.1, 0.4, -0.3].map(|v| v * speed),
            &gate,
        )
        .unwrap();
    let warm = motion
        .advance(initial, [0.0; 3], 0.01 * time, &gate)
        .unwrap()
        .state();

    // Independent four-DOF native system M=m I. Old J=[I|0]; new J=[I|b].
    // The SAME four native velocities carry the same total energy in both views.
    let b = [1.0, 0.5, -0.25];
    let unseen_velocity = 4.0 * speed;
    let old_velocity = warm.velocity_m_s();
    let new_velocity = core::array::from_fn(|i| old_velocity[i] + b[i] * unseen_velocity);
    let new_w = core::array::from_fn(|i| {
        core::array::from_fn(|j| ((if i == j { 1.0 } else { 0.0 }) + b[i] * b[j]) / MASS)
    });
    let native_energy = 0.5 * MASS * (dot(old_velocity, old_velocity) + unseen_velocity.powi(2));
    let expected_old = 0.5 * MASS * dot(old_velocity, old_velocity);
    // Sherman-Morrison inverse of I+b b^T, independent of the production LU.
    let expected_new = 0.5
        * MASS
        * (dot(new_velocity, new_velocity) - dot(b, new_velocity).powi(2) / (1.0 + dot(b, b)));
    assert!(expected_old < native_energy && expected_new < native_energy);
    assert!((expected_new - expected_old).abs() > force * h);

    let transfer = motion
        .rebind_native_projection(warm, identity(), new_w, new_velocity, &gate)
        .unwrap();
    let report = transfer.report();
    close(
        report.projected_kinetic_energy_before_j,
        expected_old,
        native_energy,
        2e-15,
    );
    close(
        report.projected_kinetic_energy_after_j,
        expected_new,
        native_energy,
        2e-15,
    );
    close(
        report.projected_kinetic_repartition_j,
        expected_new - expected_old,
        native_energy,
        2e-15,
    );
    assert_eq!(report.source_energy_residual_j, 0.0);
    assert_eq!(report.before, report.after);
    assert_eq!(report.evaluations, 1);
    assert_eq!(transfer.state().position_m(), warm.position_m());
    assert_eq!(transfer.state().velocity_m_s(), new_velocity);
    assert_eq!(transfer.state().time_s(), warm.time_s());
    assert_eq!(transfer.state().steps(), warm.steps());
    assert_eq!(transfer.motion().mobility_per_kg(), new_w);
    let (mut rebound, state, _) = transfer.into_parts();
    let next = rebound
        .advance(state, [0.0; 3], 0.005 * time, &gate)
        .unwrap();
    assert_eq!(next.state().steps(), warm.steps() + 1);
    assert!(next.state().time_s() > warm.time_s());
    assert!(next.state().position_m()[0] > warm.position_m()[0]);
}

#[test]
fn normal_rotation_and_common_material_spin_transport_the_original_elastic_vector() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, force, speed, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let initial = motion
        .initialize([2.0 * h, 0.3 * h, -0.2 * h], [0.0; 3], &gate)
        .unwrap();
    let warm = motion
        .advance(initial, [0.0; 3], 0.01 * time, &gate)
        .unwrap()
        .state();
    let tilted = TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    for spin in [0.0, 0.5 * PI] {
        let rotation = TangentialHistoryRotation::between(frame(), tilted, spin).unwrap();
        let p = projection(rotation);
        let transfer = motion
            .rebind_native_projection(
                warm,
                rotation,
                congruence(p, mobility()),
                mul(p, warm.velocity_m_s()),
                &gate,
            )
            .unwrap();
        let u = transfer.state().position_m();
        let world: [f64; 3] =
            core::array::from_fn(|i| u[1] * tilted.longitudinal()[i] + u[2] * tilted.lateral()[i]);
        let old = warm.position_m();
        // Minimum rotation z->x sends old x->-z and old y->y. A subsequent
        // positive quarter-turn about x sends old x->y and old y->z.
        let expected = if spin == 0.0 {
            [0.0, old[2], -old[1]]
        } else {
            [0.0, old[1], old[2]]
        };
        for i in 0..3 {
            close(world[i], expected[i], h, 3e-15);
        }
        close(
            transfer.report().source_energy_residual_j,
            0.0,
            force * h,
            2e-13,
        );
        let reverse = TangentialHistoryRotation::between(tilted, frame(), -spin).unwrap();
        let restored_velocity = mul(projection(reverse), transfer.state().velocity_m_s());
        let restored = transfer
            .motion()
            .rebind_native_projection(
                transfer.state(),
                reverse,
                mobility(),
                restored_velocity,
                &gate,
            )
            .unwrap();
        for i in 0..3 {
            close(
                restored.state().position_m()[i],
                warm.position_m()[i],
                h,
                4e-15,
            );
            close(
                restored.state().velocity_m_s()[i],
                warm.velocity_m_s()[i],
                speed,
                4e-15,
            );
        }
        assert_eq!(restored.state().time_s(), warm.time_s());
        assert_eq!(restored.state().steps(), warm.steps());
    }
    let antipodal = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
    assert!(TangentialHistoryRotation::between(frame(), antipodal, 0.0).is_err());
}

#[test]
fn original_owner_capacity_energy_work_budget_and_cancellation_remain_authoritative() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, _, speed, _) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let q = [2.0 * h, 0.731 * h, -0.419 * h];
    let mut policy = budget(&source);
    policy.maximum_evaluations = 5;
    let motion = JkrMixedModeMotion::new(&source, &shear, mobility(), policy).unwrap();
    let initial = motion.initialize(q, [0.0; 3], &gate).unwrap();
    let rotation =
        TangentialHistoryRotation::between(frame(), frame().rotated(0.713).unwrap(), 0.0).unwrap();
    let transferred = motion
        .rebind_native_projection(initial, rotation, mobility(), [0.0; 3], &gate)
        .unwrap();
    assert_eq!(transferred.report().evaluations, 2);
    assert_eq!(transferred.motion().budget(), policy);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.rebind_native_projection(initial, rotation, mobility(), [0.0; 3], &cancelled),
        Err(JkrMixedModeMotionError::Cancelled)
    ));
    assert_eq!(
        motion
            .rebind_native_projection(initial, identity(), mobility(), [0.0; 3], &gate)
            .unwrap()
            .state(),
        initial
    );
    let mut other_policy = policy;
    other_policy.maximum_evaluations += 1;
    let foreign = JkrMixedModeMotion::new(&source, &shear, mobility(), other_policy).unwrap();
    assert!(matches!(
        foreign.rebind_native_projection(initial, rotation, mobility(), [0.0; 3], &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));
    let mut indefinite = mobility();
    indefinite[0][0] = -indefinite[0][0];
    assert!(
        motion
            .rebind_native_projection(initial, rotation, indefinite, [0.0; 3], &gate)
            .is_err()
    );
    for bad_velocity in [[f64::NAN, 0.0, 0.0], [101.0 * speed, 0.0, 0.0]] {
        assert!(matches!(
            motion.rebind_native_projection(initial, rotation, mobility(), bad_velocity, &gate),
            Err(JkrMixedModeMotionError::Range(_))
        ));
    }

    let response = source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let stress = response.tangential_force_n[0].hypot(response.tangential_force_n[1])
        / response.contact_area_m2;
    let (near_source, near_shear) = super::path::sourced_law(1.01 * stress);
    let near = JkrMixedModeMotion::new(&near_source, &near_shear, mobility(), budget(&near_source))
        .unwrap();
    let near_state = near.initialize(q, [0.0; 3], &gate).unwrap();
    let near_transfer = near
        .rebind_native_projection(near_state, rotation, mobility(), [0.0; 3], &gate)
        .unwrap();
    let after = near_transfer.report().after;
    let capacity = near_shear
        .law()
        .static_capacity_n(after.friction_contact_state().unwrap())
        .unwrap();
    close(
        after.tangential_force_n[0].hypot(after.tangential_force_n[1]) / capacity,
        1.0 / 1.01,
        1.0,
        2e-13,
    );
    let (weak_source, weak_shear) = super::path::sourced_law(0.99 * stress);
    let weak = JkrMixedModeMotion::new(&weak_source, &weak_shear, mobility(), budget(&weak_source))
        .unwrap();
    assert!(matches!(
        weak.initialize(q, [0.0; 3], &gate),
        Err(JkrMixedModeMotionError::Path(
            JkrMixedModePathError::StaticCapacity { .. }
        ))
    ));
    assert!(matches!(
        weak.rebind_native_projection(near_state, rotation, mobility(), [0.0; 3], &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));

    // Exact-zero storage allowance must refuse real floating-point drift;
    // the production operation must not silently repair the source energy.
    policy.energy_absolute_tolerance_j = 0.0;
    policy.energy_relative_tolerance = 0.0;
    let strict = JkrMixedModeMotion::new(&source, &shear, mobility(), policy).unwrap();
    // At this higher shear fraction the original total energy resolves the
    // last-bit changes of vector norm instead of rounding them all away.
    let limit = source.critical_tangential_displacement_m(q[0]).unwrap();
    let high_shear = [q[0], 0.72 * limit, 0.54 * limit];
    let reference = source
        .evaluate_displacement_controlled(high_shear[0], [high_shear[1], high_shear[2]])
        .unwrap();
    let exact = strict.initialize(high_shear, [0.0; 3], &gate).unwrap();
    let mut observed_roundoff = false;
    for angle in [0.713, 0.371, 1.137, -2.281, 0.17, 1.89, 2.73, -0.93] {
        let r = TangentialHistoryRotation::between(frame(), frame().rotated(angle).unwrap(), 0.0)
            .unwrap();
        let u = r.apply([high_shear[1], high_shear[2]]).unwrap();
        let changed = source
            .evaluate_displacement_controlled(high_shear[0], u)
            .unwrap();
        if changed.reversible_energy_j != reference.reversible_energy_j {
            observed_roundoff = true;
            assert!(matches!(
                strict.rebind_native_projection(exact, r, mobility(), [0.0; 3], &gate),
                Err(JkrMixedModeMotionError::EnergyBalance {
                    tolerance_j: 0.0,
                    ..
                })
            ));
            break;
        }
    }
    assert!(
        observed_roundoff,
        "fixture must exercise a nonzero original source roundoff defect"
    );
    assert_eq!(
        strict
            .rebind_native_projection(exact, identity(), mobility(), [0.0; 3], &gate)
            .unwrap()
            .state(),
        exact
    );
}

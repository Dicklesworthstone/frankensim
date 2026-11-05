//! Thermal work on live coupled trajectories, through original source receipts.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::mixed_mode::{motion::*, path::JkrMixedModePathError};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;
use fs_material::state_point::MaterialPropertySelection;

#[path = "../adhesive_temperature/common.rs"]
mod original;

const MASS: f64 = 1e-5;
fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}
fn sources(data: &original::Fixture) -> (SourceJkrMixedMode, SourceAdhesiveShear) {
    let normal = data.jkr();
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        SourceJkrMixedMode::new(normal, shear_envelope()).unwrap(),
        shear,
    )
}
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let normal = source.normal_source();
    let h = normal.indentation_scale_m();
    let force = PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m();
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
    [[1.0, 0.2, -0.15], [0.2, 0.8, 0.12], [-0.15, 0.12, 1.1]].map(|r| r.map(|v| v / MASS))
}
fn refresh(
    source: &SourceJkrMixedMode,
    shear: &SourceAdhesiveShear,
    q: [f64; 3],
    data: &original::Fixture,
    t: AdhesiveTemperatures,
) -> (
    SourceJkrThermalRefresh,
    SourceJkrMixedMode,
    SourceAdhesiveShear,
) {
    let normal = source.normal_source();
    let update = normal
        .refresh_at_held_indentation(
            normal.evaluate_attached(q[0]).unwrap().state(),
            data.cards(),
            t,
        )
        .unwrap();
    let next = SourceJkrMixedMode::new(update.source().clone(), source.envelope()).unwrap();
    let next_shear = shear.rebind_jkr_temperature(normal, &update).unwrap();
    (update, next, next_shear)
}

#[test]
fn thermal_transfer_preserves_live_history_and_accounts_for_both_elastic_stores() {
    let data = original::Fixture::with_options(1.0, 1.0, true, 320.0);
    let (source, shear) = sources(&data);
    let (h, force, speed, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let start = motion
        .initialize(
            [2.0 * h, 0.3 * h, -0.2 * h],
            [0.1 * speed, 0.05 * speed, -0.04 * speed],
            &gate,
        )
        .unwrap();
    let previous = motion
        .advance(start, [0.0; 3], 0.02 * time, &gate)
        .unwrap()
        .state();
    assert!(previous.steps() > 0 && previous.time_s() > 0.0);
    let (normal_refresh, successor, next_shear) = refresh(
        &source,
        &shear,
        previous.position_m(),
        &data,
        temperatures(310.0, 290.0, 315.0),
    );
    let candidate = motion
        .continue_at_temperature(previous, &normal_refresh, &successor, &next_shear, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(
        candidate.state().position_m().map(f64::to_bits),
        previous.position_m().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().velocity_m_s().map(f64::to_bits),
        previous.velocity_m_s().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().time_s().to_bits(),
        previous.time_s().to_bits()
    );
    assert_eq!(candidate.state().steps(), previous.steps());
    assert_eq!(
        candidate.motion().mobility_per_kg(),
        motion.mobility_per_kg()
    );
    assert_eq!(
        candidate.motion().effective_mass_matrix_kg(),
        motion.effective_mass_matrix_kg()
    );
    assert_eq!(candidate.motion().budget(), motion.budget());
    let report = *candidate.report();
    assert_ne!(report.before.contact_area_m2, report.after.contact_area_m2);
    assert_ne!(report.normal_elastic_change_j, 0.0);
    assert_ne!(report.tangential_elastic_change_j, 0.0);
    assert_ne!(report.surface_energy_change_j, 0.0);
    close(
        report.parameter_work_j,
        report.normal_elastic_change_j
            + report.tangential_elastic_change_j
            + report.surface_energy_change_j,
        force * h,
        2e-14,
    );
    close(
        report.storage_change_j,
        report.after.reversible_energy_j - report.before.reversible_energy_j,
        force * h,
        0.0,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert!(
        (report.parameter_work_j - normal_refresh.parameter_work_j()).abs() > 1e-6 * force * h,
        "pure-normal parameter work must not replace the actual sheared-contact supply"
    );
    assert_eq!(report.evaluations, 7);
    for (old, new) in shear
        .resolved()
        .properties()
        .iter()
        .zip(next_shear.resolved().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    close(
        successor
            .normal_source()
            .normal_state()
            .elastic()
            .surface_a_state()
            .young_modulus_pa(),
        9e5,
        1e6,
        1e-14,
    );
    close(
        successor
            .normal_source()
            .normal_state()
            .elastic()
            .surface_b_state()
            .young_modulus_pa(),
        2.2e6,
        2e6,
        1e-14,
    );
    assert_eq!(next_shear.resolved().state_coordinate("T"), Some(315.0));
    assert_eq!(next_shear.resolved().state_coordinate("p"), Some(2e5));

    let (mut next, transferred, _) = candidate.into_parts();
    let resumed = next
        .advance(transferred, [0.0; 3], 0.02 * time, &gate)
        .unwrap();
    let frozen = motion
        .advance(previous, [0.0; 3], 0.02 * time, &gate)
        .unwrap();
    assert_ne!(
        resumed.state().velocity_m_s(),
        frozen.state().velocity_m_s()
    );
    assert_eq!(resumed.state().steps(), previous.steps() + 1);
    assert_eq!(resumed.state().time_s(), previous.time_s() + 0.02 * time);
    let mut independent =
        JkrMixedModeMotion::new(&successor, &next_shear, mobility(), motion.budget()).unwrap();
    let independent_start = independent
        .initialize(previous.position_m(), previous.velocity_m_s(), &gate)
        .unwrap();
    let expected = independent
        .advance(independent_start, [0.0; 3], 0.02 * time, &gate)
        .unwrap();
    assert_eq!(resumed.state().position_m(), expected.state().position_m());
    assert_eq!(
        resumed.state().velocity_m_s(),
        expected.state().velocity_m_s()
    );
    assert_eq!(resumed.force_impulse_n_s, expected.force_impulse_n_s);
    assert_eq!(resumed.contact_work_j, expected.contact_work_j);
}

#[test]
fn no_op_keeps_exact_token_and_a_thermal_cycle_returns_signed_storage_work() {
    let data = original::Fixture::new();
    let (source, shear) = sources(&data);
    let (h, force, speed, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let initial = motion
        .initialize([h, 0.3 * h, -0.15 * h], [0.05 * speed; 3], &gate)
        .unwrap();
    let previous = motion
        .advance(initial, [0.0; 3], 0.02 * time, &gate)
        .unwrap()
        .state();
    let (same, same_source, same_shear) = refresh(
        &source,
        &shear,
        previous.position_m(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let no_op = motion
        .continue_at_temperature(previous, &same, &same_source, &same_shear, &gate)
        .unwrap();
    assert_eq!(no_op.state(), previous);
    assert_eq!(no_op.report().before, no_op.report().after);
    assert_eq!(no_op.report().parameter_work_j, 0.0);
    assert_eq!(no_op.report().storage_change_j, 0.0);
    let (mut unchanged, same_state, _) = no_op.into_parts();
    assert_eq!(
        motion
            .advance(previous, [0.0; 3], 0.02 * time, &gate)
            .unwrap(),
        unchanged
            .advance(same_state, [0.0; 3], 0.02 * time, &gate)
            .unwrap()
    );

    let (hot, hot_source, hot_shear) = refresh(
        &source,
        &shear,
        previous.position_m(),
        &data,
        temperatures(305.0, 295.0, 310.0),
    );
    let warmed = motion
        .continue_at_temperature(previous, &hot, &hot_source, &hot_shear, &gate)
        .unwrap();
    let warm_work = warmed.report().parameter_work_j;
    let (warm_motion, warm_state, _) = warmed.into_parts();
    let (cold, cold_source, cold_shear) = refresh(
        &hot_source,
        &hot_shear,
        warm_state.position_m(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let restored = warm_motion
        .continue_at_temperature(warm_state, &cold, &cold_source, &cold_shear, &gate)
        .unwrap();
    assert_eq!(restored.state().position_m(), previous.position_m());
    assert_eq!(restored.state().velocity_m_s(), previous.velocity_m_s());
    assert_eq!(restored.state().time_s(), previous.time_s());
    assert_eq!(restored.state().steps(), previous.steps());
    assert!(warm_work * restored.report().parameter_work_j < 0.0);
    close(
        warm_work + restored.report().parameter_work_j,
        0.0,
        force * h,
        3e-14,
    );
}

#[test]
fn refreshed_static_capacity_or_mixed_mode_stability_loss_refuses_without_reset() {
    let gate = CancelGate::new_clock_free();
    let data = original::Fixture::new();
    let (source, shear) = sources(&data);
    let (h, _, _, _) = scales(&source);
    let delta = 2.0 * h;
    let mut low: f64 = 0.0;
    let mut high: f64 = 0.99 * source.critical_tangential_displacement_m(delta).unwrap();
    for _ in 0..48 {
        let u = low.midpoint(high);
        let response = source
            .evaluate_displacement_controlled(delta, [u, 0.0])
            .unwrap();
        if response.tangential_force_n[0] / response.contact_area_m2
            < 0.8 * shear.law().static_shear_strength_pa()
        {
            low = u;
        } else {
            high = u;
        }
    }
    let motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let state = motion
        .initialize([delta, low.midpoint(high), 0.0], [0.0; 3], &gate)
        .unwrap();
    let (hot, successor, next_shear) = refresh(
        &source,
        &shear,
        state.position_m(),
        &data,
        temperatures(280.0, 280.0, 320.0),
    );
    let q = state.position_m();
    let after = successor
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    assert!(
        after.tangential_force_n[0]
            > next_shear.law().static_shear_strength_pa() * after.contact_area_m2
    );
    assert!(matches!(
        motion.continue_at_temperature(state, &hot, &successor, &next_shear, &gate),
        Err(JkrMixedModeMotionError::Path(
            JkrMixedModePathError::StaticCapacity { .. }
        ))
    ));
    assert_eq!(motion.initialize(q, [0.0; 3], &gate).unwrap(), state);

    let weak_data = original::Fixture::with_options(1.0, 0.01, false, 320.0);
    let (weak, weak_shear) = sources(&weak_data);
    let weak_motion =
        JkrMixedModeMotion::new(&weak, &weak_shear, mobility(), budget(&weak)).unwrap();
    let u = 0.75 * weak.critical_tangential_displacement_m(0.0).unwrap();
    let attached = weak_motion
        .initialize([0.0, u, 0.0], [0.0; 3], &gate)
        .unwrap();
    let (hot, successor, next_shear) = refresh(
        &weak,
        &weak_shear,
        attached.position_m(),
        &weak_data,
        temperatures(300.0, 300.0, 320.0),
    );
    assert!(matches!(
        weak_motion.continue_at_temperature(attached, &hot, &successor, &next_shear, &gate),
        Err(JkrMixedModeMotionError::Path(
            JkrMixedModePathError::Source(JkrMixedModeError::NoFullStickEquilibrium { .. })
        ))
    ));
    assert_eq!(
        weak_motion
            .initialize(attached.position_m(), [0.0; 3], &gate)
            .unwrap(),
        attached
    );
}

#[test]
fn foreign_history_refresh_shear_cancellation_and_shared_budget_are_atomic() {
    let data = original::Fixture::new();
    let (source, shear) = sources(&data);
    let (h, _, _, _) = scales(&source);
    let q = [h, 0.2 * h, -0.1 * h];
    let gate = CancelGate::new_clock_free();
    let motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let state = motion.initialize(q, [0.0; 3], &gate).unwrap();
    let (hot, successor, next_shear) =
        refresh(&source, &shear, q, &data, temperatures(305.0, 295.0, 310.0));
    assert!(matches!(
        motion.continue_at_temperature(state, &hot, &successor, &shear, &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));
    assert!(matches!(
        motion.continue_at_temperature(state, &hot, &source, &next_shear, &gate),
        Err(JkrMixedModeMotionError::Invalid(_))
    ));
    let unrelated = JkrMixedModeMotion::new(
        &source,
        &shear,
        mobility().map(|r| r.map(|v| 2.0 * v)),
        budget(&source),
    )
    .unwrap();
    let wrong_state = unrelated.initialize(q, [0.0; 3], &gate).unwrap();
    assert!(matches!(
        motion.continue_at_temperature(wrong_state, &hot, &successor, &next_shear, &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));
    let mut other_q = q;
    other_q[0] += 0.1 * h;
    let (wrong, wrong_source, wrong_shear) = refresh(
        &source,
        &shear,
        other_q,
        &data,
        temperatures(305.0, 295.0, 310.0),
    );
    assert!(matches!(
        motion.continue_at_temperature(state, &wrong, &wrong_source, &wrong_shear, &gate),
        Err(JkrMixedModeMotionError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.continue_at_temperature(state, &hot, &successor, &next_shear, &cancelled),
        Err(JkrMixedModeMotionError::Cancelled)
    ));
    let mut limited_budget = budget(&source);
    limited_budget.maximum_evaluations = 6;
    let limited = JkrMixedModeMotion::new(&source, &shear, mobility(), limited_budget).unwrap();
    let limited_state = limited.initialize(q, [0.0; 3], &gate).unwrap();
    assert!(matches!(
        limited.continue_at_temperature(limited_state, &hot, &successor, &next_shear, &gate),
        Err(JkrMixedModeMotionError::Path(JkrMixedModePathError::Budget))
    ));
    assert_eq!(motion.initialize(q, [0.0; 3], &gate).unwrap(), state);
}

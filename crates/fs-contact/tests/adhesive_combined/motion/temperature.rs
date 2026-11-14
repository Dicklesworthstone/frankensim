//! Original four-mode storage, independent temperatures and retained chronology.
#[path = "../../adhesive_temperature/common.rs"]
mod original;

use super::{budget, close, mobility, scales, start};
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, SourceJkrSphereContact,
    mixed_mode::{
        SourceJkrMixedMode,
        combined::{JkrCombinedResponse, motion::*},
        torsion::motion::{JkrTorsionalMotion, JkrTorsionalMotionBudget},
    },
    temperature::{AdhesiveTemperatures, SourceJkrThermalRefresh},
};
use fs_exec::CancelGate;

fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}
fn thermal_source(data: &original::Fixture) -> SourceJkrMixedMode {
    let template = super::source();
    SourceJkrMixedMode::new(
        SourceJkrSphereContact::new(
            data.normal.clone(),
            template.normal_source().effective_radius_m(),
            template.normal_source().interaction_range_m(),
            template.normal_source().envelope(),
        )
        .unwrap(),
        template.envelope(),
    )
    .unwrap()
}
fn refresh(
    source: &SourceJkrMixedMode,
    q: [f64; 4],
    data: &original::Fixture,
    t: AdhesiveTemperatures,
) -> (SourceJkrThermalRefresh, SourceJkrMixedMode) {
    let normal = source.normal_source();
    let update = normal
        .refresh_at_held_indentation(
            normal.evaluate_attached(q[0]).unwrap().state(),
            data.cards(),
            t,
        )
        .unwrap();
    let successor = SourceJkrMixedMode::new(update.source().clone(), source.envelope()).unwrap();
    (update, successor)
}

// Independent evaluation of the original four-term potential. Only the
// admitted radius comes from the source; no reported energy is used here.
fn components(source: &SourceJkrMixedMode, response: JkrCombinedResponse) -> [f64; 4] {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let g = source.reduced_shear_modulus_pa();
    let gt = source.torsional_reduced_shear_modulus_pa();
    let q = response.coordinates;
    let a = response.contact_radius_m;
    let hertz = a * a / normal.effective_radius_m();
    [
        e * a * ((q[0] - hertz / 3.0).powi(2) + 4.0 * hertz * hertz / 45.0),
        4.0 * g * a * (q[1] * q[1] + q[2] * q[2]),
        (8.0 / 3.0) * gt * a.powi(3) * q[3] * q[3],
        -PI * w * a * a,
    ]
}

#[test]
fn g1_live_all_mode_temperature_uses_one_relaxed_area_and_complete_potential() {
    let data = original::Fixture::with_options(10.0, 0.05, true, 320.0);
    let source = thermal_source(&data);
    let (scale, force, time) = scales(&source);
    let energy = force * scale[0];
    let gate = CancelGate::new_clock_free();
    let motion = JkrCombinedMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let (q, rates) = start(&source);
    let initial = motion.initial_state(q, rates, &gate).unwrap();
    let previous = motion
        .advance(initial, 0.02 * time, [0.0; 4], &gate)
        .unwrap()
        .state();
    let (normal, successor) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(310.0, 305.0, 315.0),
    );
    let candidate = motion
        .continue_at_temperature(previous, &normal, &successor, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(candidate.state().coordinates(), previous.coordinates());
    assert_eq!(candidate.state().rates(), previous.rates());
    assert_eq!(candidate.state().time_s(), previous.time_s());
    assert_eq!(candidate.state().steps(), previous.steps());
    assert!(previous.steps() > 0 && previous.time_s() > 0.0);
    assert_eq!(candidate.motion().mobility(), motion.mobility());
    assert_eq!(candidate.motion().budget(), motion.budget());
    let report = *candidate.report();
    assert_eq!(report.before, previous.contact());
    assert_eq!(report.after, candidate.state().contact());
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, successor.identity());
    assert_ne!(
        report.before.contact_radius_m,
        report.after.contact_radius_m
    );
    let changes = [
        report.normal_elastic_change_j,
        report.tangential_elastic_change_j,
        report.torsional_elastic_change_j,
        report.surface_energy_change_j,
    ];
    let before = components(&source, report.before);
    let after = components(&successor, report.after);
    for i in 0..4 {
        assert_ne!(changes[i], 0.0);
        close(changes[i], after[i] - before[i], energy, 3e-13);
    }
    close(
        report.parameter_work_j,
        after.into_iter().sum::<f64>() - before.into_iter().sum::<f64>(),
        energy,
        4e-13,
    );
    close(
        report.parameter_work_j,
        changes.into_iter().sum(),
        energy,
        3e-14,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert_eq!(report.evaluations, 2);
    assert!((report.parameter_work_j - normal.parameter_work_j()).abs() > 1e-8 * energy);
    // The two original bulk queries and the interface query remain independent,
    // with pinned original claims despite competing hot-only alternatives.
    let material = successor.normal_source().normal_state();
    close(
        material.elastic().surface_a_state().young_modulus_pa(),
        9e6,
        1e7,
        2e-14,
    );
    close(
        material.elastic().surface_b_state().young_modulus_pa(),
        19e6,
        2e7,
        2e-14,
    );
    close(material.adhesion_energy_j_per_m2(), 0.0275, 0.05, 2e-14);
    assert_eq!(
        material.interface_state().state_coordinate("T"),
        Some(315.0)
    );
    for (old, new) in source
        .normal_source()
        .normal_state()
        .interface_state()
        .properties()
        .iter()
        .zip(material.interface_state().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    let (next, transferred, _) = candidate.into_parts();
    let resumed = next
        .advance(transferred, 0.02 * time, [0.0; 4], &gate)
        .unwrap();
    let frozen = motion
        .advance(previous, 0.02 * time, [0.0; 4], &gate)
        .unwrap();
    assert_ne!(resumed.state().rates(), frozen.state().rates());
    assert_eq!(resumed.state().steps(), previous.steps() + 1);
    assert_eq!(resumed.state().time_s(), previous.time_s() + 0.02 * time);
    let independent =
        JkrCombinedMotion::new(&successor, motion.mobility(), motion.budget()).unwrap();
    let start = independent
        .initial_state(previous.coordinates(), previous.rates(), &gate)
        .unwrap();
    let expected = independent
        .advance(start, 0.02 * time, [0.0; 4], &gate)
        .unwrap();
    assert_eq!(
        resumed.state().coordinates(),
        expected.state().coordinates()
    );
    assert_eq!(resumed.state().rates(), expected.state().rates());
    assert_eq!(resumed.effort_impulse, expected.effort_impulse);
    assert_eq!(resumed.contact_work_j, expected.contact_work_j);
}

#[test]
fn g1_exact_no_op_round_trip_and_original_shear_twist_reductions() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = thermal_source(&data);
    let (scale, force, time) = scales(&source);
    let energy = force * scale[0];
    let gate = CancelGate::new_clock_free();
    let motion = JkrCombinedMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let (q, rates) = start(&source);
    let previous = motion
        .advance(
            motion.initial_state(q, rates, &gate).unwrap(),
            0.02 * time,
            [0.0; 4],
            &gate,
        )
        .unwrap()
        .state();
    let (same, same_source) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let unchanged = motion
        .continue_at_temperature(previous, &same, &same_source, &gate)
        .unwrap();
    assert_eq!(unchanged.state(), previous);
    assert_eq!(unchanged.report().parameter_work_j, 0.0);
    assert_eq!(unchanged.report().storage_change_j, 0.0);
    let (warm, warm_source) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(305.0, 295.0, 310.0),
    );
    let warmed = motion
        .continue_at_temperature(previous, &warm, &warm_source, &gate)
        .unwrap();
    let (back, original) = refresh(
        &warm_source,
        previous.coordinates(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let restored = warmed
        .motion()
        .continue_at_temperature(warmed.state(), &back, &original, &gate)
        .unwrap();
    // The way back re-resolves the cards with pinned-claim receipts, so the
    // identities differ from the original's single-claim-only ones; the
    // chronology and every physical value must return exactly.
    original::same_after_refresh(&restored.state(), &previous);
    close(
        warmed.report().parameter_work_j + restored.report().parameter_work_j,
        0.0,
        energy,
        5e-14,
    );
    original::same_after_refresh(
        &restored
            .motion()
            .advance(restored.state(), 0.02 * time, [0.0; 4], &gate)
            .unwrap(),
        &motion
            .advance(previous, 0.02 * time, [0.0; 4], &gate)
            .unwrap(),
    );
    // Exact zero-twist reduction retains the original sheared contact energies.
    let shear_q = [q[0], q[1], q[2], 0.0];
    let shear_state = motion.initial_state(shear_q, [0.0; 4], &gate).unwrap();
    let (refresh_shear, successor) =
        refresh(&source, shear_q, &data, temperatures(310.0, 305.0, 315.0));
    let transfer = motion
        .continue_at_temperature(shear_state, &refresh_shear, &successor, &gate)
        .unwrap();
    let before = source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let after = successor
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    assert_eq!(
        transfer.report().before.contact_radius_m,
        before.contact_radius_m
    );
    assert_eq!(
        transfer.report().after.contact_radius_m,
        after.contact_radius_m
    );
    assert_eq!(transfer.report().torsional_elastic_change_j, 0.0);
    close(
        transfer.report().parameter_work_j,
        after.reversible_energy_j - before.reversible_energy_j,
        energy,
        5e-14,
    );
    // Exact zero-shear reduction agrees with the original coupled I/III driver.
    let twist_q = [q[0], 0.0, 0.0, q[3]];
    let twist_state = motion.initial_state(twist_q, [0.0; 4], &gate).unwrap();
    let (refresh_twist, successor) =
        refresh(&source, twist_q, &data, temperatures(310.0, 305.0, 315.0));
    let transfer = motion
        .continue_at_temperature(twist_state, &refresh_twist, &successor, &gate)
        .unwrap();
    let b = motion.budget();
    let pair = |values: [f64; 4]| [values[0], values[3]];
    let w = motion.mobility();
    let torsion = JkrTorsionalMotion::new(
        &source,
        [[w[0][0], w[0][3]], [w[3][0], w[3][3]]],
        JkrTorsionalMotionBudget {
            maximum_step_s: b.maximum_step_s,
            maximum_coordinate_increment: pair(b.maximum_coordinate_increment),
            maximum_rate: pair(b.maximum_rate),
            maximum_effort: pair(b.maximum_effort),
            maximum_iterations: b.maximum_iterations,
            maximum_evaluations: b.maximum_evaluations,
            maximum_backtracks: b.maximum_backtracks,
            coordinate_absolute_tolerance: pair(b.coordinate_absolute_tolerance),
            coordinate_relative_tolerance: b.coordinate_relative_tolerance,
            rate_absolute_tolerance: pair(b.rate_absolute_tolerance),
            rate_relative_tolerance: b.rate_relative_tolerance,
            energy_absolute_tolerance_j: b.energy_absolute_tolerance_j,
            energy_relative_tolerance: b.energy_relative_tolerance,
        },
    )
    .unwrap();
    let original_state = torsion
        .initial_state([q[0], q[3]], [0.0; 2], &gate)
        .unwrap();
    let original_transfer = torsion
        .continue_at_temperature(original_state, &refresh_twist, &successor, &gate)
        .unwrap();
    assert_eq!(transfer.report().tangential_elastic_change_j, 0.0);
    assert_eq!(
        transfer.report().after.contact_radius_m,
        original_transfer.report().after.contact_radius_m
    );
    assert_eq!(
        transfer.report().parameter_work_j,
        original_transfer.report().parameter_work_j
    );
}

#[test]
fn g3_normal_thermal_admission_does_not_hide_combined_instability_or_effort_limits() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = thermal_source(&data);
    let (scale, _, _) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let motion = JkrCombinedMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let twist_limit = source
        .torsional_displacement_limit(scale[0])
        .unwrap()
        .unwrap()
        .twist_angle_rad;
    let q = [
        scale[0],
        0.001 * scale[1],
        -0.001 * scale[2],
        0.85 * twist_limit,
    ];
    let previous = motion.initial_state(q, [0.0; 4], &gate).unwrap();
    let (hot, successor) = refresh(&source, q, &data, temperatures(300.0, 300.0, 320.0));
    assert_eq!(hot.response().branch, JkrBranch::Attached);
    assert!(
        successor
            .torsional_displacement_limit(q[0])
            .unwrap()
            .unwrap()
            .twist_angle_rad
            < q[3]
    );
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &successor, &gate),
        Err(JkrCombinedMotionError::Path(_))
    ));
    assert_eq!(motion.initial_state(q, [0.0; 4], &gate).unwrap(), previous);
    let (q, _) = start(&source);
    let (cold, stronger) = refresh(&source, q, &data, temperatures(300.0, 300.0, 280.0));
    let old = source.evaluate_combined_displacement_controlled(q).unwrap();
    let new = stronger
        .evaluate_combined_displacement_controlled(q)
        .unwrap();
    assert!(new.force[3] > old.force[3]);
    let mut limited = budget(&source);
    limited.maximum_effort[3] = old.force[3].midpoint(new.force[3]);
    let constrained = JkrCombinedMotion::new(&source, mobility(&source), limited).unwrap();
    let state = constrained.initial_state(q, [0.0; 4], &gate).unwrap();
    assert!(matches!(
        constrained.continue_at_temperature(state, &cold, &stronger, &gate),
        Err(JkrCombinedMotionError::Range(_))
    ));
    assert_eq!(
        constrained.initial_state(q, [0.0; 4], &gate).unwrap(),
        state
    );
}

#[test]
fn g4_foreign_history_source_control_budget_and_cancellation_leave_history_unchanged() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = thermal_source(&data);
    let gate = CancelGate::new_clock_free();
    let motion = JkrCombinedMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let (q, _) = start(&source);
    let previous = motion.initial_state(q, [0.0; 4], &gate).unwrap();
    let (hot, successor) = refresh(&source, q, &data, temperatures(305.0, 295.0, 310.0));
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &source, &gate),
        Err(JkrCombinedMotionError::Invalid(_))
    ));
    let wrong_q = [1.1 * q[0], q[1], q[2], q[3]];
    let (wrong, wrong_source) = refresh(&source, wrong_q, &data, temperatures(305.0, 295.0, 310.0));
    assert!(matches!(
        motion.continue_at_temperature(previous, &wrong, &wrong_source, &gate),
        Err(JkrCombinedMotionError::Invalid(_))
    ));
    let foreign = JkrCombinedMotion::new(
        &source,
        mobility(&source).map(|row| row.map(|x| 2.0 * x)),
        budget(&source),
    )
    .unwrap();
    let foreign_state = foreign.initial_state(q, [0.0; 4], &gate).unwrap();
    assert!(matches!(
        motion.continue_at_temperature(foreign_state, &hot, &successor, &gate),
        Err(JkrCombinedMotionError::StateMismatch)
    ));
    let normal = source.normal_source();
    let load = normal.evaluate_attached(q[0]).unwrap().normal_force_n;
    let held_load = normal
        .refresh_at_held_load(load, data.cards(), temperatures(300.0, 300.0, 300.0))
        .unwrap();
    let held_source =
        SourceJkrMixedMode::new(held_load.source().clone(), source.envelope()).unwrap();
    assert!(matches!(
        motion.continue_at_temperature(previous, &held_load, &held_source, &gate),
        Err(JkrCombinedMotionError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &successor, &cancelled),
        Err(JkrCombinedMotionError::Cancelled)
    ));
    let mut b = budget(&source);
    b.maximum_evaluations = 1;
    let limited = JkrCombinedMotion::new(&source, mobility(&source), b).unwrap();
    let state = limited.initial_state(q, [0.0; 4], &gate).unwrap();
    assert!(matches!(
        limited.continue_at_temperature(state, &hot, &successor, &gate),
        Err(JkrCombinedMotionError::Budget)
    ));
    b.maximum_evaluations = 2;
    let exact = JkrCombinedMotion::new(&source, mobility(&source), b).unwrap();
    let state = exact.initial_state(q, [0.0; 4], &gate).unwrap();
    assert_eq!(
        exact
            .continue_at_temperature(state, &hot, &successor, &gate)
            .unwrap()
            .report()
            .evaluations,
        2
    );
    assert_eq!(motion.initial_state(q, [0.0; 4], &gate).unwrap(), previous);
    motion
        .continue_at_temperature(previous, &hot, &successor, &gate)
        .unwrap();
}

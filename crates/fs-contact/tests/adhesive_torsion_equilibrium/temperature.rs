//! Original temperature sources, live coupled history and signed parameter work.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::motion::*;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};

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
fn thermal_source(data: &original::Fixture) -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(
        data.jkr(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let normal = source.normal_source();
    let length = normal.indentation_scale_m();
    let angle = 0.2
        * source
            .torsional_displacement_limit(length)
            .unwrap()
            .unwrap()
            .twist_angle_rad;
    let force = PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m();
    (length, angle, force, (MASS * length / force).sqrt())
}
fn mobility(source: &SourceJkrMixedMode) -> [[f64; 2]; 2] {
    let (length, angle, _, _) = scales(source);
    let inverse_inertia = (angle / length).powi(2) / MASS;
    let cross = 0.25 * (inverse_inertia / MASS).sqrt();
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
fn refresh(
    source: &SourceJkrMixedMode,
    q: [f64; 2],
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
    let next = SourceJkrMixedMode::new(update.source().clone(), source.envelope()).unwrap();
    (update, next)
}

// Independent original I/III potential; the thermal continuation does not get
// to replace this changed-area storage by the pure-normal refresh's work.
fn potential(source: &SourceJkrMixedMode, response: JkrTorsionalResponse) -> f64 {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let g = source.torsional_reduced_shear_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let a = response.contact_radius_m;
    let hertz = a * a / normal.effective_radius_m();
    e * a * ((response.indentation_m - hertz / 3.0).powi(2) + (4.0 / 45.0) * hertz * hertz)
        + (8.0 / 3.0) * g * a.powi(3) * response.twist_angle_rad.powi(2)
        - PI * w * a * a
}

#[test]
fn g1_live_thermal_transfer_retains_twist_and_uses_all_original_energy_components() {
    let data = original::Fixture::with_options(1.0, 1.0, true, 320.0);
    let source = thermal_source(&data);
    let (length, angle, force, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let motion = JkrTorsionalMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let initial = motion
        .initial_state(
            [length, angle],
            [0.05 * length / time, 0.1 * angle / time],
            &gate,
        )
        .unwrap();
    let previous = motion
        .advance(initial, 0.02 * time, [0.0; 2], &gate)
        .unwrap()
        .state();
    assert!(previous.steps() > 0 && previous.time_s() > 0.0);
    let (normal_refresh, successor) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(310.0, 290.0, 315.0),
    );
    let candidate = motion
        .continue_at_temperature(previous, &normal_refresh, &successor, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(
        candidate.state().coordinates().map(f64::to_bits),
        previous.coordinates().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().rates().map(f64::to_bits),
        previous.rates().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().time_s().to_bits(),
        previous.time_s().to_bits()
    );
    assert_eq!(candidate.state().steps(), previous.steps());
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
    assert_ne!(report.normal_elastic_change_j, 0.0);
    assert_ne!(report.torsional_elastic_change_j, 0.0);
    assert_ne!(report.surface_energy_change_j, 0.0);
    let energy = force * length;
    work_close(
        report.parameter_work_j,
        potential(&successor, report.after) - potential(&source, report.before),
        energy,
        3e-13,
    );
    work_close(
        report.parameter_work_j,
        report.normal_elastic_change_j
            + report.torsional_elastic_change_j
            + report.surface_energy_change_j,
        energy,
        3e-14,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert!((report.parameter_work_j - normal_refresh.parameter_work_j()).abs() > 1e-6 * energy);
    assert_eq!(report.evaluations, 2);
    // Independent temperatures and pinned selected claims survive even when
    // the fixture offers conflicting hot-only material/interface alternatives.
    let normal = successor.normal_source().normal_state();
    close(
        normal.elastic().surface_a_state().young_modulus_pa(),
        9e5,
        2e-14,
    );
    close(
        normal.elastic().surface_b_state().young_modulus_pa(),
        2.2e6,
        2e-14,
    );
    close(normal.adhesion_energy_j_per_m2(), 0.55, 2e-14);
    assert_eq!(normal.interface_state().state_coordinate("T"), Some(315.0));
    assert_eq!(normal.interface_state().state_coordinate("p"), Some(2e5));
    for (old, new) in source
        .normal_source()
        .normal_state()
        .interface_state()
        .properties()
        .iter()
        .zip(normal.interface_state().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    let (next, transferred, _) = candidate.into_parts();
    let resumed = next
        .advance(transferred, 0.02 * time, [0.0; 2], &gate)
        .unwrap();
    let frozen = motion
        .advance(previous, 0.02 * time, [0.0; 2], &gate)
        .unwrap();
    assert_ne!(resumed.state().rates(), frozen.state().rates());
    assert_eq!(resumed.state().steps(), previous.steps() + 1);
    assert_eq!(resumed.state().time_s(), previous.time_s() + 0.02 * time);
    let independent =
        JkrTorsionalMotion::new(&successor, motion.mobility(), motion.budget()).unwrap();
    let expected_start = independent
        .initial_state(previous.coordinates(), previous.rates(), &gate)
        .unwrap();
    let expected = independent
        .advance(expected_start, 0.02 * time, [0.0; 2], &gate)
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
fn g1_no_op_is_exact_and_warming_cooling_returns_signed_parameter_work_without_motion() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (length, angle, force, time) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let motion = JkrTorsionalMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let start = motion
        .initial_state(
            [length, angle],
            [0.03 * length / time, -0.05 * angle / time],
            &gate,
        )
        .unwrap();
    let previous = motion
        .advance(start, 0.02 * time, [0.0; 2], &gate)
        .unwrap()
        .state();
    let (same, same_source) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let no_op = motion
        .continue_at_temperature(previous, &same, &same_source, &gate)
        .unwrap();
    assert_eq!(no_op.state(), previous);
    assert_eq!(no_op.report().before, no_op.report().after);
    assert_eq!(no_op.report().parameter_work_j, 0.0);
    assert_eq!(no_op.report().storage_change_j, 0.0);
    assert_eq!(
        motion
            .advance(previous, 0.02 * time, [0.0; 2], &gate)
            .unwrap(),
        no_op
            .motion()
            .advance(no_op.state(), 0.02 * time, [0.0; 2], &gate)
            .unwrap()
    );
    let (hot, hot_source) = refresh(
        &source,
        previous.coordinates(),
        &data,
        temperatures(305.0, 295.0, 310.0),
    );
    let warmed = motion
        .continue_at_temperature(previous, &hot, &hot_source, &gate)
        .unwrap();
    let warm_work = warmed.report().parameter_work_j;
    let (cold, cold_source) = refresh(
        &hot_source,
        warmed.state().coordinates(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let restored = warmed
        .motion()
        .continue_at_temperature(warmed.state(), &cold, &cold_source, &gate)
        .unwrap();
    assert_eq!(restored.state().coordinates(), previous.coordinates());
    assert_eq!(restored.state().rates(), previous.rates());
    assert_eq!(restored.state().time_s(), previous.time_s());
    assert_eq!(restored.state().steps(), previous.steps());
    assert!(warm_work * restored.report().parameter_work_j < 0.0);
    work_close(
        warm_work + restored.report().parameter_work_j,
        0.0,
        force * length,
        5e-14,
    );
    work_close(
        restored.state().contact().reversible_energy_j,
        previous.contact().reversible_energy_j,
        force * length,
        5e-14,
    );
}

#[test]
fn g3_torsional_attachment_and_effort_limits_are_rechecked_after_normal_refresh() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (length, angle, _, _) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let motion = JkrTorsionalMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let limit = source
        .torsional_displacement_limit(length)
        .unwrap()
        .unwrap()
        .twist_angle_rad;
    let q = [length, 0.85 * limit];
    let previous = motion.initial_state(q, [0.0; 2], &gate).unwrap();
    let (hot, successor) = refresh(&source, q, &data, temperatures(300.0, 300.0, 320.0));
    assert_eq!(
        hot.response().branch,
        fs_contact::normal_patch::adhesive::JkrBranch::Attached
    );
    assert!(
        successor
            .torsional_displacement_limit(length)
            .unwrap()
            .unwrap()
            .twist_angle_rad
            < q[1]
    );
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &successor, &gate),
        Err(JkrTorsionalMotionError::Path(
            JkrTorsionalPathError::Source(JkrTorsionalError::BeyondAttachedLimit { .. })
        ))
    ));
    assert_eq!(motion.initial_state(q, [0.0; 2], &gate).unwrap(), previous);

    let q = [length, angle];
    let (cold, stronger) = refresh(&source, q, &data, temperatures(300.0, 300.0, 280.0));
    let old = source
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap();
    let new = stronger
        .evaluate_torsion_displacement_controlled(q[0], q[1])
        .unwrap();
    assert!(new.torsional_moment_n_m > old.torsional_moment_n_m);
    let mut limited = budget(&source);
    limited.maximum_effort[1] = old.torsional_moment_n_m.midpoint(new.torsional_moment_n_m);
    let constrained = JkrTorsionalMotion::new(&source, mobility(&source), limited).unwrap();
    let state = constrained.initial_state(q, [0.0; 2], &gate).unwrap();
    assert!(matches!(
        constrained.continue_at_temperature(state, &cold, &stronger, &gate),
        Err(JkrTorsionalMotionError::Range(_))
    ));
    assert_eq!(
        constrained.initial_state(q, [0.0; 2], &gate).unwrap(),
        state
    );
}

#[test]
fn g4_foreign_controls_sources_history_budget_and_missing_claim_support_refuse_atomically() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (length, angle, _, _) = scales(&source);
    let gate = CancelGate::new_clock_free();
    let q = [length, angle];
    let motion = JkrTorsionalMotion::new(&source, mobility(&source), budget(&source)).unwrap();
    let previous = motion.initial_state(q, [0.0; 2], &gate).unwrap();
    let (hot, successor) = refresh(&source, q, &data, temperatures(305.0, 295.0, 310.0));
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &source, &gate),
        Err(JkrTorsionalMotionError::Invalid(_))
    ));
    let (other, other_source) = refresh(
        &source,
        [1.1 * length, angle],
        &data,
        temperatures(305.0, 295.0, 310.0),
    );
    assert!(matches!(
        motion.continue_at_temperature(previous, &other, &other_source, &gate),
        Err(JkrTorsionalMotionError::Invalid(_))
    ));
    let foreign = JkrTorsionalMotion::new(
        &source,
        mobility(&source).map(|row| row.map(|v| 2.0 * v)),
        budget(&source),
    )
    .unwrap();
    let foreign_state = foreign.initial_state(q, [0.0; 2], &gate).unwrap();
    assert!(matches!(
        motion.continue_at_temperature(foreign_state, &hot, &successor, &gate),
        Err(JkrTorsionalMotionError::StateMismatch)
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
        Err(JkrTorsionalMotionError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.continue_at_temperature(previous, &hot, &successor, &cancelled),
        Err(JkrTorsionalMotionError::Cancelled)
    ));
    let mut limited = budget(&source);
    limited.maximum_evaluations = 1;
    let limited = JkrTorsionalMotion::new(&source, mobility(&source), limited).unwrap();
    let state = limited.initial_state(q, [0.0; 2], &gate).unwrap();
    assert!(matches!(
        limited.continue_at_temperature(state, &hot, &successor, &gate),
        Err(JkrTorsionalMotionError::Budget)
    ));
    let mut exact = budget(&source);
    exact.maximum_evaluations = 2;
    let exact = JkrTorsionalMotion::new(&source, mobility(&source), exact).unwrap();
    let state = exact.initial_state(q, [0.0; 2], &gate).unwrap();
    assert_eq!(
        exact
            .continue_at_temperature(state, &hot, &successor, &gate)
            .unwrap()
            .report()
            .evaluations,
        2
    );
    // Original selected claims expire at 305 K. Hot-only alternatives exist
    // but may not silently replace the selected source to continue adhesion.
    let unsupported = original::Fixture::with_options(1.0, 1.0, true, 305.0);
    let original = thermal_source(&unsupported);
    let h = original.normal_source().indentation_scale_m();
    let contact = original.normal_source().evaluate_attached(h).unwrap();
    assert!(
        original
            .normal_source()
            .refresh_at_held_indentation(
                contact.state(),
                unsupported.cards(),
                temperatures(315.0, 300.0, 315.0),
            )
            .is_err()
    );
    assert_eq!(motion.initial_state(q, [0.0; 2], &gate).unwrap(), previous);
    motion
        .continue_at_temperature(previous, &hot, &successor, &gate)
        .unwrap();
}

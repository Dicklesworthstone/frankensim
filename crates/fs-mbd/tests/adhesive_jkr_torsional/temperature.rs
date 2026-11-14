//! Original thermal coefficient work with native torsion and chronology retained.
use super::*;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod original;

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

fn thermal_scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64, f64) {
    let normal = source.normal_source();
    let length = normal.indentation_scale_m();
    let force =
        PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m() / 6.0;
    let speed = (force * length / (1.2 * MASS)).sqrt();
    let stiffness = source
        .evaluate_torsion_displacement_controlled(0.0, 0.0)
        .unwrap()
        .coupled_stiffness
        .unwrap()[1][1];
    (
        length,
        speed,
        length / speed,
        (force * length / stiffness).sqrt(),
        force,
    )
}

fn thermal_controls(source: &SourceJkrMixedMode) -> RigidJkrTorsionalControls {
    let (length, speed, time, twist, force) = thermal_scales(source);
    let torque = force * length / twist;
    RigidJkrTorsionalControls {
        motion: JkrTorsionalMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_coordinate_increment: [20.0 * length, 20.0 * twist],
            maximum_rate: [100.0 * speed, 100.0 * twist / time],
            maximum_effort: [1e4 * force, 1e4 * torque],
            maximum_iterations: 64,
            maximum_evaluations: 2048,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [1e-13 * length, 1e-13 * twist],
            coordinate_relative_tolerance: 2e-12,
            rate_absolute_tolerance: [1e-12 * speed, 1e-12 * twist / time],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-10 * force * length,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 4,
        // These COMs lie near +/-ARM; pose rounding is distinct from the
        // original constitutive solve tolerance and cannot reset approach.
        geometry_tolerance_m: 256.0 * f64::EPSILON * ARM,
        displacement_tolerance_m: (1e-9 * length).max(256.0 * f64::EPSILON * ARM),
        twist_tolerance_rad: 1e-9 * twist,
        velocity_tolerance_m_s: 1e-9 * speed,
        angular_velocity_tolerance_rad_s: 1e-9 * twist / time,
        linear_momentum_tolerance_ns: 1e-11 * MASS * speed,
        angular_momentum_tolerance_nms: 1e-11 * MASS * speed * ARM,
        energy_tolerance_j: 1e-8 * force * length,
    }
}

fn thermal_contact<'a>(
    f: &Fixture,
    source: &'a SourceJkrMixedMode,
) -> RigidJkrTorsionalContact<'a> {
    let mut contact = f.contact(source);
    contact.reference_approach_m = 2.0 * source.normal_source().indentation_scale_m();
    contact
}

fn prepare<'a>(f: &Fixture, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrTorsionalPair<'a> {
    PreparedRigidJkrTorsionalPair::new(
        &f.bodies,
        &f.properties,
        thermal_contact(f, source),
        thermal_controls(source),
        &CancelGate::new(),
    )
    .unwrap()
}

fn refresh(
    data: &original::Fixture,
    source: &SourceJkrMixedMode,
    previous: RigidJkrTorsionalState,
    temperatures: [f64; 3],
) -> (SourceJkrThermalRefresh, SourceJkrMixedMode) {
    let normal = source.normal_source();
    let update = normal
        .refresh_at_held_indentation(
            normal
                .evaluate_attached(previous.coordinates()[0])
                .unwrap()
                .state(),
            data.cards(),
            AdhesiveTemperatures {
                surface_a_k: temperatures[0],
                surface_b_k: temperatures[1],
                interface_k: temperatures[2],
            },
        )
        .unwrap();
    let successor = SourceJkrMixedMode::new(update.source().clone(), source.envelope()).unwrap();
    (update, successor)
}

#[test]
fn unequal_temperatures_preserve_native_motion_and_supply_full_twisted_storage_work() {
    let data = original::Fixture::with_options(1.0, 1.0, true, 320.0);
    let source = thermal_source(&data);
    let (_, speed, time, twist, _) = thermal_scales(&source);
    let limits = thermal_controls(&source);
    let gate = CancelGate::new();
    let mut f = axial();
    f.set_rates(
        0,
        Vec3::new(0.0, 0.0, 0.03 * speed),
        Vec3::new(0.0, 0.0, 0.02 * twist / time),
    );
    let pair = prepare(&f, &source);
    let previous = pair.initialize(0.3 * twist, &gate).unwrap();
    let (update, successor) = refresh(&data, &source, previous, [310.0, 290.0, 315.0]);
    let candidate = pair
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(candidate.bodies(), &f.bodies);
    assert_ne!(candidate.bodies()[0].linear_momentum_world(), Vec3::ZERO);
    assert_ne!(candidate.bodies()[0].angular_momentum_body(), Vec3::ZERO);
    assert_eq!(candidate.state().coordinates(), previous.coordinates());
    assert_eq!(candidate.state().rates(), previous.rates());
    assert_eq!(candidate.state().time_s(), previous.time_s());
    assert_eq!(candidate.state().steps(), previous.steps());
    assert_eq!(candidate.pair().mobility(), pair.mobility());
    assert_eq!(
        candidate.pair().reaction_arms_body(),
        pair.reaction_arms_body()
    );
    assert_eq!(candidate.pair().approach_m(), pair.approach_m());
    assert_eq!(candidate.pair().contact().frame, pair.contact().frame);
    assert_eq!(
        candidate.pair().contact().reference_approach_m,
        pair.contact().reference_approach_m
    );
    let old = pair.contact().binding;
    let new = candidate.pair().contact().binding;
    assert_eq!(
        (new.body_a, new.body_b, new.arm_a_body, new.arm_b_body),
        (old.body_a, old.body_b, old.arm_a_body, old.arm_b_body)
    );
    assert_eq!(candidate.pair().contact().source, &successor);
    let report = *candidate.report();
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, successor.identity());
    assert_ne!(
        report.original_model_identity,
        report.refreshed_model_identity
    );
    assert_eq!(candidate.state().motion_state().contact(), report.after);
    assert_ne!(report.before.contact_area_m2, report.after.contact_area_m2);
    assert_ne!(report.normal_elastic_change_j, 0.0);
    assert_ne!(report.torsional_elastic_change_j, 0.0);
    assert_ne!(report.surface_energy_change_j, 0.0);
    close(
        report.parameter_work_j,
        report.normal_elastic_change_j
            + report.torsional_elastic_change_j
            + report.surface_energy_change_j,
        limits.energy_tolerance_j,
    );
    close(
        report.parameter_work_j,
        report.after.reversible_energy_j - report.before.reversible_energy_j,
        limits.energy_tolerance_j,
    );
    assert!(report.energy_residual_j.abs() <= limits.energy_tolerance_j);
    assert!(
        (report.parameter_work_j - update.parameter_work_j()).abs()
            > 100.0 * limits.energy_tolerance_j
    );

    let (hot, state, _) = candidate.into_parts();
    assert!(hot.advance(previous, 0.002 * time, &gate).is_err());
    let advanced = hot.advance(state, 0.002 * time, &gate).unwrap();
    let frozen = pair.advance(previous, 0.002 * time, &gate).unwrap();
    assert_ne!(advanced.bodies_after, frozen.bodies_after);
    assert_eq!(advanced.state_after.steps(), previous.steps() + 1);
    assert_eq!(
        advanced.state_after.time_s(),
        previous.time_s() + 0.002 * time
    );
    let kinetic_change: f64 = (0..2)
        .map(|i| {
            f.bodies[i]
                .kinetic_energy_change_to(f.properties[i], advanced.bodies_after[i])
                .unwrap()
        })
        .sum();
    close(
        kinetic_change + advanced.report.motion.contact.reversible_energy_j
            - report.before.reversible_energy_j
            - report.parameter_work_j,
        0.0,
        2.0 * limits.energy_tolerance_j,
    );
}

#[test]
fn cooling_reverses_warming_parameter_work_without_an_impulse_or_clock_reset() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (_, _, _, twist, _) = thermal_scales(&source);
    let gate = CancelGate::new();
    let f = axial();
    let pair = prepare(&f, &source);
    let previous = pair.initialize(0.25 * twist, &gate).unwrap();
    let (warming, warm_source) = refresh(&data, &source, previous, [315.0; 3]);
    let warm = pair
        .continue_at_temperature(previous, &warming, &warm_source, &gate)
        .unwrap();
    let work = warm.report().parameter_work_j;
    let (warm_pair, warm_state, _) = warm.into_parts();
    let (cooling, cold_source) = refresh(&data, &warm_source, warm_state, [300.0; 3]);
    let cold = warm_pair
        .continue_at_temperature(warm_state, &cooling, &cold_source, &gate)
        .unwrap();
    assert!(work * cold.report().parameter_work_j < 0.0);
    close(
        work + cold.report().parameter_work_j,
        0.0,
        2.0 * thermal_controls(&source).energy_tolerance_j,
    );
    same_physical_source(&cold_source, &source);
    same_values(cold.state(), previous);
    assert_eq!(cold.bodies(), &f.bodies);
    assert_eq!(cold.pair().mobility(), pair.mobility());
    same_values(cold.report().after, previous.motion_state().contact());
}

#[test]
fn thermal_transfer_requires_actual_normal_drift_and_retains_nonzero_clock() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (_, _, time, twist, _) = thermal_scales(&source);
    let gate = CancelGate::new();
    let f = axial();
    let pair = prepare(&f, &source);
    let original = pair.initialize(0.2 * twist, &gate).unwrap();
    let accepted = pair.advance(original, 0.002 * time, &gate).unwrap();
    let (update, successor) = refresh(&data, &source, accepted.state_after, [300.0; 3]);
    let held = PreparedRigidJkrTorsionalPair::new(
        &accepted.bodies_after,
        &f.properties,
        thermal_contact(&f, &source),
        thermal_controls(&source),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        held.continue_at_temperature(accepted.state_after, &update, &successor, &gate),
        Err(RigidJkrTorsionalError::GeometryMismatch { .. })
    ));
    let geometry = pair
        .continue_axial_geometry(accepted.state_after, &gate)
        .unwrap();
    let native_bodies = *geometry.bodies();
    let (current, previous, _) = geometry.into_parts();
    let continuation = current
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
    assert_eq!(continuation.state(), previous);
    assert_eq!(continuation.state().steps(), 1);
    assert_eq!(continuation.state().time_s(), 0.002 * time);
    assert_eq!(continuation.bodies(), &native_bodies);
    assert_eq!(continuation.report().parameter_work_j, 0.0);
    assert_eq!(continuation.report().storage_change_j, 0.0);
    let (same, state, _) = continuation.into_parts();
    assert_eq!(
        same.advance(state, 0.002 * time, &gate).unwrap(),
        current.advance(previous, 0.002 * time, &gate).unwrap()
    );
}

#[test]
fn held_twist_travel_survives_noop_refresh_before_compatible_attitude_drift() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (_, _, time, twist, _) = thermal_scales(&source);
    let gate = CancelGate::new();
    let dt = 0.002 * time;
    let mut f = axial();
    f.set_rates(0, Vec3::ZERO, Vec3::new(0.0, 0.0, 0.5 * twist / time));
    let pair = prepare(&f, &source);
    let previous = pair.initialize(0.15 * twist, &gate).unwrap();
    let q0 = previous.coordinates();
    let q1 = [q0[0], q0[1] + dt * previous.rates()[1]];
    let effort = source
        .torsional_discrete_path(q0, q1, &gate)
        .unwrap()
        .mean_effort;
    // These source-owned discrete reactions balance the body impulses while
    // the genuine initial spin still advances the elastic twist and attitude.
    let wrenches = core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        Wrench {
            force_world: Vec3::new(0.0, 0.0, sign * effort[0]),
            torque_body: Vec3::new(0.0, 0.0, sign * effort[1]),
        }
    });
    let accepted = pair
        .advance_with_wrenches(previous, wrenches, dt, &gate)
        .unwrap();
    assert!(accepted.report.angular_displacement_body_rad[0].z > 0.0);
    let held = PreparedRigidJkrTorsionalPair::new(
        &accepted.bodies_after,
        &f.properties,
        thermal_contact(&f, &source),
        thermal_controls(&source),
        &gate,
    )
    .unwrap();
    let (update, successor) = refresh(&data, &source, accepted.state_after, [300.0; 3]);
    let thermal = held
        .continue_at_temperature(accepted.state_after, &update, &successor, &gate)
        .unwrap();
    assert_eq!(thermal.state(), accepted.state_after);
    let without_refresh = held
        .continue_axial_geometry(accepted.state_after, &gate)
        .unwrap();
    let with_refresh = thermal
        .pair()
        .continue_axial_geometry(thermal.state(), &gate)
        .unwrap();
    assert_eq!(with_refresh.bodies(), without_refresh.bodies());
    assert_eq!(with_refresh.state(), without_refresh.state());
    assert_ne!(
        with_refresh.bodies()[0].pose().orientation(),
        f.bodies[0].pose().orientation()
    );
    assert_eq!(with_refresh.state().steps(), 1);
    assert_eq!(with_refresh.state().time_s(), dt);
}

#[test]
fn foreign_sources_native_staleness_and_cancellation_leave_temperature_retry_unchanged() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let (_, speed, time, twist, _) = thermal_scales(&source);
    let gate = CancelGate::new();
    let f = axial();
    let pair = prepare(&f, &source);
    let previous = pair.initialize(0.2 * twist, &gate).unwrap();
    let (update, successor) = refresh(&data, &source, previous, [310.0, 290.0, 315.0]);
    let baseline = pair
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
    let trial = pair.advance(previous, 0.002 * time, &gate).unwrap();
    assert!(
        pair.continue_at_temperature(previous, &update, &source, &gate)
            .is_err()
    );
    let foreign_data = original::Fixture::with_options(1.2, 1.0, false, 320.0);
    let foreign_source = thermal_source(&foreign_data);
    let (foreign_update, foreign_successor) =
        refresh(&foreign_data, &foreign_source, previous, [310.0; 3]);
    assert!(
        pair.continue_at_temperature(previous, &foreign_update, &foreign_successor, &gate)
            .is_err()
    );
    let mut changed = axial();
    changed.set_rates(0, Vec3::new(0.0, 0.0, 0.01 * speed), Vec3::ZERO);
    assert!(matches!(
        prepare(&changed, &source).continue_at_temperature(previous, &update, &successor, &gate),
        Err(RigidJkrTorsionalError::StateMismatch)
    ));
    let mut changed_contact = thermal_contact(&f, &source);
    changed_contact.frame = TangentFrame::new([0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]).unwrap();
    let changed = PreparedRigidJkrTorsionalPair::new(
        &f.bodies,
        &f.properties,
        changed_contact,
        thermal_controls(&source),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        changed.continue_at_temperature(previous, &update, &successor, &gate),
        Err(RigidJkrTorsionalError::StateMismatch)
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        pair.continue_at_temperature(previous, &update, &successor, &cancelled)
            .is_err()
    );
    let retry = pair
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
    assert_eq!(retry.bodies(), &f.bodies);
    assert_eq!(retry.previous_state(), previous);
    assert_eq!(retry.state(), baseline.state());
    assert_eq!(retry.report(), baseline.report());
    assert_eq!(pair.advance(previous, 0.002 * time, &gate).unwrap(), trial);
}

/// Back at the original temperatures, the cooled source was refreshed, so its
/// receipts record pinned-claim selection where the original resolution
/// recorded SingleClaimOnly, and its identity differs by that receipt alone
/// (fs-material's rule: a refresh pins, an independent resolution does not).
/// Every card, query, selected claim and value is the original's.
fn same_physical_source(refreshed: &SourceJkrMixedMode, original: &SourceJkrMixedMode) {
    let (a, b) = (
        refreshed.normal_source().normal_state(),
        original.normal_source().normal_state(),
    );
    for (x, y) in [
        (
            a.elastic().surface_a_state().resolved(),
            b.elastic().surface_a_state().resolved(),
        ),
        (
            a.elastic().surface_b_state().resolved(),
            b.elastic().surface_b_state().resolved(),
        ),
    ] {
        assert_eq!(x.card_identity(), y.card_identity());
        assert_eq!(x.query(), y.query());
        same_claims(x.properties(), y.properties());
    }
    let (x, y) = (a.interface_state(), b.interface_state());
    assert_eq!(x.card_identity(), y.card_identity());
    assert_eq!(x.query(), y.query());
    same_claims(x.properties(), y.properties());
    assert_eq!(
        a.adhesion_energy_j_per_m2().to_bits(),
        b.adhesion_energy_j_per_m2().to_bits()
    );
    assert_eq!(refreshed.envelope(), original.envelope());
}

fn same_claims(
    refreshed: &[fs_material::state_point::ResolvedScalarProperty],
    original: &[fs_material::state_point::ResolvedScalarProperty],
) {
    assert_eq!(refreshed.len(), original.len());
    for (x, y) in refreshed.iter().zip(original) {
        assert_eq!(x.answer().receipt.selected, y.answer().receipt.selected);
        assert_eq!(x.value_si().to_bits(), y.value_si().to_bits());
    }
}

/// The cooled state holds the original's values; its content hashes are
/// derived from the refreshed source's identity, so they differ by the same
/// receipt (see `same_physical_source`). Every field is compared except those
/// hashes. Debug prints each f64 in its shortest round-trip form, so equal
/// text means equal bits.
fn same_values<T: core::fmt::Debug>(cooled: T, original: T) {
    assert_eq!(
        without_hashes(&format!("{cooled:?}")),
        without_hashes(&format!("{original:?}"))
    );
}

fn without_hashes(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(start) = rest.find("ContentHash(") {
        let open = start + "ContentHash(".len();
        let close = open + rest[open..].find(')').unwrap();
        out.push_str(&rest[..open]);
        out.push_str("..)");
        rest = &rest[close + 1..];
    }
    out.push_str(rest);
    out
}

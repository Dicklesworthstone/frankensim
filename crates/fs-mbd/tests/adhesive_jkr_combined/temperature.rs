//! G1/G3/G4: the original coupled area, actual native history and thermal work.
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

fn thermal_scales(source: &SourceJkrMixedMode) -> [f64; 6] {
    let normal = source.normal_source();
    let length = normal.indentation_scale_m();
    let force =
        PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m() / 6.0;
    let speed = (force * length / (1.2 * MASS)).sqrt();
    let stiffness = source
        .evaluate_combined_displacement_controlled([0.0; 4])
        .unwrap()
        .coupled_stiffness
        .unwrap();
    [
        length,
        speed,
        length / speed,
        (force * length / stiffness[1][1]).sqrt(),
        (force * length / stiffness[3][3]).sqrt(),
        force,
    ]
}

fn thermal_controls(source: &SourceJkrMixedMode) -> RigidJkrCombinedControls {
    let [length, speed, time, shear, twist, force] = thermal_scales(source);
    let torque = force * length / twist;
    RigidJkrCombinedControls {
        motion: JkrCombinedMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_coordinate_increment: [20.0 * length, 20.0 * shear, 20.0 * shear, 20.0 * twist],
            maximum_rate: [
                100.0 * speed,
                100.0 * shear / time,
                100.0 * shear / time,
                100.0 * twist / time,
            ],
            maximum_effort: [1e4 * force, 1e4 * force, 1e4 * force, 1e4 * torque],
            maximum_iterations: 64,
            maximum_evaluations: 8192,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [
                1e-13 * length,
                1e-13 * shear,
                1e-13 * shear,
                1e-13 * twist,
            ],
            coordinate_relative_tolerance: 2e-12,
            rate_absolute_tolerance: [
                1e-12 * speed,
                1e-12 * shear / time,
                1e-12 * shear / time,
                1e-12 * twist / time,
            ],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-10 * force * length,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 4,
        geometry_tolerance_m: 256.0 * f64::EPSILON * ARM,
        displacement_tolerance_m: (1e-9 * length.max(shear)).max(256.0 * f64::EPSILON * ARM),
        twist_tolerance_rad: 1e-9 * twist,
        velocity_tolerance_m_s: 1e-9 * speed.max(shear / time),
        angular_velocity_tolerance_rad_s: 1e-9 * twist / time,
        linear_momentum_tolerance_ns: 1e-11 * MASS * speed,
        angular_momentum_tolerance_nms: 1e-11 * MASS * speed * ARM,
        energy_tolerance_j: 1e-8 * force * length,
    }
}

fn thermal_contact<'a>(f: &Fixture, source: &'a SourceJkrMixedMode) -> RigidJkrCombinedContact<'a> {
    let mut contact = f.contact(source);
    contact.reference_approach_m = 2.0 * source.normal_source().indentation_scale_m();
    contact
}

fn prepare<'a>(f: &Fixture, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrCombinedPair<'a> {
    PreparedRigidJkrCombinedPair::new(
        &f.bodies,
        &f.properties,
        thermal_contact(f, source),
        thermal_controls(source),
        &CancelGate::new(),
    )
    .unwrap()
}

fn preload(source: &SourceJkrMixedMode) -> [f64; 3] {
    let [_, _, _, shear, twist, _] = thermal_scales(source);
    [0.18 * shear, -0.12 * shear, 0.2 * twist]
}

fn refresh(
    data: &original::Fixture,
    source: &SourceJkrMixedMode,
    previous: RigidJkrCombinedState,
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
fn independent_temperatures_reprice_the_complete_loaded_area_without_native_impulse() {
    let data = original::Fixture::with_options(1.0, 1.0, true, 320.0);
    let source = thermal_source(&data);
    let [_, speed, time, shear, twist, _] = thermal_scales(&source);
    let limits = thermal_controls(&source);
    let gate = CancelGate::new();
    let mut f = facing();
    f.set_rates(
        0,
        Vec3::new(0.02 * shear / time, -0.01 * shear / time, 0.03 * speed),
        Vec3::new(0.0, 0.0, 0.02 * twist / time),
    );
    let pair = prepare(&f, &source);
    let previous = pair.initialize(preload(&source), &gate).unwrap();
    assert!(previous.coordinates().into_iter().all(|q| q != 0.0));
    assert!(previous.rates().into_iter().all(|rate| rate != 0.0));
    let (update, successor) = refresh(&data, &source, previous, [310.0, 290.0, 315.0]);
    let candidate = pair
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(candidate.bodies(), &f.bodies);
    assert_eq!(candidate.properties(), &f.properties);
    assert_eq!(candidate.pair().bodies(), pair.bodies());
    assert_eq!(candidate.pair().properties(), pair.properties());
    assert_eq!(candidate.state().coordinates(), previous.coordinates());
    assert_eq!(candidate.state().rates(), previous.rates());
    assert_eq!(candidate.state().time_s(), previous.time_s());
    assert_eq!(candidate.state().steps(), previous.steps());
    assert_eq!(candidate.pair().mobility(), pair.mobility());
    assert_eq!(
        candidate.pair().reaction_arms_body(),
        pair.reaction_arms_body()
    );
    assert_eq!(candidate.pair().contact().frame, pair.contact().frame);
    assert_eq!(candidate.pair().approach_m(), pair.approach_m());
    assert_eq!(candidate.pair().contact().source, &successor);
    let report = *candidate.report();
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, successor.identity());
    assert_ne!(
        report.original_model_identity,
        report.refreshed_model_identity
    );
    assert_ne!(report.before.contact_area_m2, report.after.contact_area_m2);
    let direct = successor
        .evaluate_combined_displacement_controlled(previous.coordinates())
        .unwrap();
    assert_eq!(report.after, direct);
    assert_eq!(candidate.state().motion_state().contact(), direct);
    let changes = [
        report.normal_elastic_change_j,
        report.tangential_elastic_change_j,
        report.torsional_elastic_change_j,
        report.surface_energy_change_j,
    ];
    assert!(changes.into_iter().all(|energy| energy != 0.0));
    close(
        report.parameter_work_j,
        changes.into_iter().sum(),
        limits.energy_tolerance_j,
    );
    close(
        report.parameter_work_j,
        direct.reversible_energy_j - report.before.reversible_energy_j,
        limits.energy_tolerance_j,
    );
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
    let kinetic: f64 = (0..2)
        .map(|i| {
            f.bodies[i]
                .kinetic_energy_change_to(f.properties[i], advanced.bodies_after[i])
                .unwrap()
        })
        .sum();
    close(
        kinetic + advanced.report.motion.contact.reversible_energy_j
            - report.before.reversible_energy_j
            - report.parameter_work_j,
        0.0,
        2.0 * limits.energy_tolerance_j,
    );
}

#[test]
fn full_mobility_nonzero_chronology_and_pending_native_travel_survive_round_trip() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let [_, speed, time, _, twist, _] = thermal_scales(&source);
    let gate = CancelGate::new();
    let mut f = fixture([
        Vec3::new(ARM, 0.3 * ARM, 0.2 * ARM),
        Vec3::new(0.4 * ARM, -0.2 * ARM, -0.5 * ARM),
    ]);
    for i in 0..2 {
        f.properties[i] = MassProperties::from_principal_frame(
            f.properties[i].mass(),
            Vec3::ZERO,
            f.properties[i].principal_inertia_body(),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.45 + 0.1 * i as f64)
                .unwrap(),
        )
        .unwrap();
    }
    let arms = prepare(&f, &source).reaction_arms_body();
    let common = Vec3::new(0.01 * speed, -0.02 * speed, 0.015 * speed);
    let omega = Vec3::new(0.0, 0.0, 0.02 * twist / time);
    for (i, arm) in arms.into_iter().enumerate() {
        f.set_rates(i, common.sub(omega.cross(arm)), omega);
    }
    let pair = prepare(&f, &source);
    let previous = pair.initialize(preload(&source), &gate).unwrap();
    let mobility = pair.mobility();
    for row in 0..4 {
        for column in 0..row {
            assert!(
                mobility[row][column].abs()
                    > 1e-5 * (mobility[row][row] * mobility[column][column]).sqrt()
            );
        }
    }
    let response = previous.motion_state().contact();
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    let force = axes[0]
        .scale(response.force[0])
        .add(axes[1].scale(response.force[1]))
        .add(axes[2].scale(response.force[2]));
    // Actual COM forces and body torques balance the complete original source.
    // Common native motion still contributes nonzero sealed pending travel.
    let loads = core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        Wrench {
            force_world: force.scale(sign),
            torque_body: arms[i]
                .cross(force.scale(sign))
                .add(axes[0].scale(sign * response.force[3])),
        }
    });
    let dt = 0.002 * time;
    let advanced = pair
        .advance_with_wrenches(previous, loads, dt, &gate)
        .unwrap();
    let mut current = f;
    current.bodies = advanced.bodies_after;
    let held = prepare(&current, &source);
    let previous = advanced.state_after;
    assert_eq!(previous.steps(), 1);
    assert_eq!(previous.time_s(), dt);
    assert!(
        previous
            .pending_com_travel_m()
            .into_iter()
            .all(|travel| travel != Vec3::ZERO)
    );
    assert!(
        previous
            .pending_angular_displacement_body_rad()
            .into_iter()
            .all(|travel| travel != Vec3::ZERO)
    );
    let (warming, warm_source) = refresh(&data, &source, previous, [310.0, 290.0, 315.0]);
    let warm = held
        .continue_at_temperature(previous, &warming, &warm_source, &gate)
        .unwrap();
    assert_eq!(warm.bodies(), &current.bodies);
    assert_eq!(warm.pair().mobility(), mobility);
    assert_eq!(
        warm.state().pending_com_travel_m(),
        previous.pending_com_travel_m()
    );
    assert_eq!(
        warm.state().pending_angular_displacement_body_rad(),
        previous.pending_angular_displacement_body_rad()
    );
    assert_eq!(warm.state().steps(), 1);
    assert_eq!(warm.state().time_s(), dt);
    let work = warm.report().parameter_work_j;
    let (warm_pair, warm_state, _) = warm.into_parts();
    let (cooling, cold_source) = refresh(&data, &warm_source, warm_state, [300.0; 3]);
    let cold = warm_pair
        .continue_at_temperature(warm_state, &cooling, &cold_source, &gate)
        .unwrap();
    same_physical_source(&cold_source, &source);
    same_values(cold.state(), previous);
    assert_eq!(cold.bodies(), &current.bodies);
    close(
        work + cold.report().parameter_work_j,
        0.0,
        2.0 * thermal_controls(&source).energy_tolerance_j,
    );
}

#[test]
fn source_body_staleness_geometry_and_cancellation_refuse_without_changing_retry() {
    let data = original::Fixture::new();
    let source = thermal_source(&data);
    let [_, speed, time, _, _, _] = thermal_scales(&source);
    let gate = CancelGate::new();
    let f = facing();
    let pair = prepare(&f, &source);
    let previous = pair.initialize(preload(&source), &gate).unwrap();
    let (update, successor) = refresh(&data, &source, previous, [310.0, 290.0, 315.0]);
    let baseline = pair
        .continue_at_temperature(previous, &update, &successor, &gate)
        .unwrap();
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
    let mut changed = facing();
    changed.set_rates(0, Vec3::new(0.0, 0.0, 0.01 * speed), Vec3::ZERO);
    assert!(matches!(
        prepare(&changed, &source).continue_at_temperature(previous, &update, &successor, &gate),
        Err(RigidJkrCombinedError::StateMismatch)
    ));
    let trial = pair.advance(previous, 0.002 * time, &gate).unwrap();
    let held = PreparedRigidJkrCombinedPair::new(
        &trial.bodies_after,
        &f.properties,
        thermal_contact(&f, &source),
        thermal_controls(&source),
        &gate,
    )
    .unwrap();
    let (stale_update, stale_successor) = refresh(&data, &source, trial.state_after, [300.0; 3]);
    assert!(matches!(
        held.continue_at_temperature(trial.state_after, &stale_update, &stale_successor, &gate),
        Err(RigidJkrCombinedError::GeometryMismatch { .. })
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
    assert_eq!(retry.state(), baseline.state());
    assert_eq!(retry.report(), baseline.report());
    assert_eq!(retry.bodies(), &f.bodies);
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

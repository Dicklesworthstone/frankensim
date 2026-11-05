//! Source parameter work with the captured physical pair and chronology intact.
use super::*;
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_material::state_point::MaterialPropertySelection;

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod original;

fn thermal_sources(data: &original::Fixture) -> (SourceJkrMixedMode, SourceAdhesiveShear) {
    let normal = data.jkr();
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = SourceJkrMixedMode::new(
        normal,
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap();
    (source, shear)
}

fn thermal_scales(source: &SourceJkrMixedMode) -> (f64, f64, f64) {
    let normal = source.normal_source();
    let length = normal.indentation_scale_m();
    let force =
        PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m() / 6.0;
    let speed = (force * length / (1.2 * MASS)).sqrt();
    (length, speed, length / speed)
}

fn thermal_controls(source: &SourceJkrMixedMode) -> RigidJkrMixedModeControls {
    // The original thermal card has a different R, W and stiffness from this
    // target's stiff-sphere fixture. Retain its dimensionless solve budgets.
    let normal = source.normal_source();
    let ratio = PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m()
        / (6.0 * FORCE);
    let rate = ratio.sqrt();
    let mut result = controls(source);
    result.motion.maximum_step_s /= rate;
    result.motion.maximum_speed_m_s *= rate;
    result.motion.maximum_force_n *= ratio;
    result.motion.velocity_absolute_tolerance_m_s *= rate;
    result.motion.energy_absolute_tolerance_j *= ratio;
    result.velocity_tolerance_m_s *= rate;
    result.angular_velocity_tolerance_rad_s *= rate;
    result.linear_momentum_tolerance_ns *= rate;
    result.angular_momentum_tolerance_nms *= rate;
    result.energy_tolerance_j *= ratio;
    result
}

fn thermal_contact<'a>(
    f: &Fixture,
    source: &'a SourceJkrMixedMode,
    shear: &'a SourceAdhesiveShear,
) -> RigidJkrMixedModeContact<'a> {
    let mut contact = f.contact(source, shear);
    contact.reference_approach_m = 2.0 * source.normal_source().indentation_scale_m();
    contact
}

fn prepare<'a>(
    f: &Fixture,
    source: &'a SourceJkrMixedMode,
    shear: &'a SourceAdhesiveShear,
) -> PreparedRigidJkrMixedModePair<'a> {
    PreparedRigidJkrMixedModePair::new(
        &f.bodies,
        &f.properties,
        thermal_contact(f, source, shear),
        thermal_controls(source),
        &CancelGate::new(),
    )
    .unwrap()
}

fn refresh(
    data: &original::Fixture,
    source: &SourceJkrMixedMode,
    shear: &SourceAdhesiveShear,
    previous: RigidJkrMixedModeState,
    temperatures: [f64; 3],
) -> (
    SourceJkrThermalRefresh,
    SourceJkrMixedMode,
    SourceAdhesiveShear,
) {
    let normal = source.normal_source();
    let update = normal
        .refresh_at_held_indentation(
            normal
                .evaluate_attached(previous.position_m()[0])
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
    let next_shear = shear.rebind_jkr_temperature(normal, &update).unwrap();
    (update, successor, next_shear)
}

#[test]
fn unequal_temperatures_preserve_native_momenta_and_supply_complete_mixed_work() {
    let data = original::Fixture::with_options(1.0, 1.0, true, 320.0);
    let (source, shear) = thermal_sources(&data);
    let (length, speed, time) = thermal_scales(&source);
    let limits = thermal_controls(&source);
    let gate = CancelGate::new();
    let mut f = axial();
    f.set_velocity(0, Vec3::new(0.03, -0.02, 0.01).scale(speed));
    let mut pair = prepare(&f, &source, &shear);
    let previous = pair
        .initialize([0.3 * length, -0.2 * length], &gate)
        .unwrap();
    let (update, successor, next_shear) =
        refresh(&data, &source, &shear, previous, [310.0, 290.0, 315.0]);
    let candidate = pair
        .continue_at_temperature(previous, &update, &successor, &next_shear, &gate)
        .unwrap();
    assert_eq!(candidate.previous_state(), previous);
    assert_eq!(candidate.bodies(), &f.bodies);
    assert_ne!(candidate.bodies()[0].linear_momentum_world(), Vec3::ZERO);
    assert_eq!(candidate.state().position_m(), previous.position_m());
    assert_eq!(candidate.state().velocity_m_s(), previous.velocity_m_s());
    assert_eq!(candidate.state().time_s(), previous.time_s());
    assert_eq!(candidate.state().steps(), previous.steps());
    assert_eq!(candidate.pair().mobility_per_kg(), pair.mobility_per_kg());
    assert_eq!(
        candidate.pair().reaction_arms_body(),
        pair.reaction_arms_body()
    );
    assert_eq!(candidate.pair().approach_m(), pair.approach_m());
    let binding = candidate.pair().contact().binding;
    let original_binding = pair.contact().binding;
    assert_eq!(binding.body_a, original_binding.body_a);
    assert_eq!(binding.body_b, original_binding.body_b);
    assert_eq!(binding.arm_a_body, original_binding.arm_a_body);
    assert_eq!(binding.arm_b_body, original_binding.arm_b_body);
    assert_eq!(candidate.pair().contact().frame, pair.contact().frame);
    assert_eq!(
        candidate.pair().contact().reference_approach_m,
        pair.contact().reference_approach_m
    );
    assert_eq!(candidate.pair().contact().source, &successor);
    assert_eq!(candidate.pair().contact().shear, &next_shear);
    let report = *candidate.report();
    assert_ne!(
        report.original_model_identity,
        report.refreshed_model_identity
    );
    assert_ne!(report.before.contact_area_m2, report.after.contact_area_m2);
    assert_ne!(report.normal_elastic_change_j, 0.0);
    assert_ne!(report.tangential_elastic_change_j, 0.0);
    assert_ne!(report.surface_energy_change_j, 0.0);
    close(
        report.parameter_work_j,
        report.normal_elastic_change_j
            + report.tangential_elastic_change_j
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

    let (mut hot_pair, hot_state, _) = candidate.into_parts();
    assert!(hot_pair.advance(previous, 0.02 * time, &gate).is_err());
    let advanced = hot_pair.advance(hot_state, 0.02 * time, &gate).unwrap();
    let frozen = pair.advance(previous, 0.02 * time, &gate).unwrap();
    assert_ne!(advanced.bodies_after, frozen.bodies_after);
    assert_eq!(advanced.state_after.steps(), previous.steps() + 1);
    assert_eq!(
        advanced.state_after.time_s(),
        previous.time_s() + 0.02 * time
    );
    let kinetic_change: f64 = (0..2)
        .map(|i| {
            advanced.bodies_after[i]
                .kinetic_energy(f.properties[i])
                .unwrap()
                - f.bodies[i].kinetic_energy(f.properties[i]).unwrap()
        })
        .sum();
    close(
        kinetic_change + advanced.report.motion.contact.reversible_energy_j
            - report.before.reversible_energy_j
            - report.parameter_work_j,
        0.0,
        2.0 * limits.energy_tolerance_j,
    );
    assert!(advanced.report.audit.total_energy_residual_j.abs() <= limits.energy_tolerance_j);
}

#[test]
fn no_op_temperature_transfer_retains_accepted_clock_and_exact_future_native_trial() {
    let data = original::Fixture::new();
    let (source, shear) = thermal_sources(&data);
    let (length, _, time) = thermal_scales(&source);
    let gate = CancelGate::new();
    let f = axial();
    let mut pair = prepare(&f, &source, &shear);
    let initial = pair
        .initialize([0.3 * length, -0.2 * length], &gate)
        .unwrap();
    let q = initial.position_m();
    let response = source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let force = Vec3::new(
        response.tangential_force_n[0],
        response.tangential_force_n[1],
        response.normal_force_n,
    );
    let hold = core::array::from_fn(|i| {
        let force = force.scale(if i == 0 { 1.0 } else { -1.0 });
        Wrench {
            force_world: force,
            torque_body: f.arms[i].cross(force),
        }
    });
    let first = pair
        .advance_with_wrenches(initial, hold, 0.02 * time, &gate)
        .unwrap();
    assert_eq!(first.state_after.steps(), 1);
    assert!(first.state_after.time_s() > 0.0);
    let next = Fixture {
        bodies: first.bodies_after,
        ..f
    };
    let mut current = prepare(&next, &source, &shear);
    let (update, successor, next_shear) =
        refresh(&data, &source, &shear, first.state_after, [300.0; 3]);
    let candidate = current
        .continue_at_temperature(first.state_after, &update, &successor, &next_shear, &gate)
        .unwrap();
    assert_eq!(candidate.bodies(), &next.bodies);
    assert_eq!(candidate.state(), first.state_after);
    assert_eq!(candidate.report().parameter_work_j, 0.0);
    assert_eq!(candidate.report().storage_change_j, 0.0);
    let (mut unchanged, state, _) = candidate.into_parts();
    assert_eq!(
        unchanged
            .advance_with_wrenches(state, hold, 0.02 * time, &gate)
            .unwrap(),
        current
            .advance_with_wrenches(first.state_after, hold, 0.02 * time, &gate)
            .unwrap()
    );
}

#[test]
fn native_staleness_and_cancellation_refuse_temperature_transfer_before_any_publication() {
    let data = original::Fixture::new();
    let (source, shear) = thermal_sources(&data);
    let (length, speed, _) = thermal_scales(&source);
    let limits = thermal_controls(&source);
    let gate = CancelGate::new();
    let f = axial();
    let pair = prepare(&f, &source, &shear);
    let previous = pair
        .initialize([0.2 * length, -0.1 * length], &gate)
        .unwrap();
    let (update, successor, next_shear) =
        refresh(&data, &source, &shear, previous, [310.0, 290.0, 315.0]);
    let transfer = |other: &PreparedRigidJkrMixedModePair<'_>, gate: &CancelGate| {
        other.continue_at_temperature(previous, &update, &successor, &next_shear, gate)
    };
    let mut changed = thermal_contact(&f, &source, &shear);
    changed.frame = TangentFrame::new([0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]).unwrap();
    let rotated =
        PreparedRigidJkrMixedModePair::new(&f.bodies, &f.properties, changed, limits, &gate)
            .unwrap();
    assert!(matches!(
        transfer(&rotated, &gate),
        Err(RigidJkrMixedModeError::StateMismatch)
    ));

    let mut changed = axial();
    changed.set_velocity(0, Vec3::new(0.01 * speed, 0.0, 0.0));
    assert!(matches!(
        transfer(&prepare(&changed, &source, &shear), &gate),
        Err(RigidJkrMixedModeError::StateMismatch)
    ));

    let mut changed = axial();
    let old = changed.bodies[0];
    changed.bodies[0] = RigidBodyState::new(
        Pose::new(
            old.pose()
                .position_world()
                .add(Vec3::new(0.0, 0.0, 0.01 * length)),
            old.pose().orientation(),
        )
        .unwrap(),
        old.linear_momentum_world(),
        old.angular_momentum_body(),
    )
    .unwrap();
    assert!(matches!(
        transfer(&prepare(&changed, &source, &shear), &gate),
        Err(RigidJkrMixedModeError::GeometryMismatch { .. })
    ));

    let mut changed = axial();
    changed.properties[0] = MassProperties::new(
        2.1 * MASS,
        Vec3::ZERO,
        changed.properties[0].principal_inertia_body(),
    )
    .unwrap();
    assert!(transfer(&prepare(&changed, &source, &shear), &gate).is_err());

    let mut changed = axial();
    changed.bodies[0] = RigidBodyState::new(
        changed.bodies[0].pose(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, INERTIA * speed / ARM),
    )
    .unwrap();
    assert!(matches!(
        transfer(&prepare(&changed, &source, &shear), &gate),
        Err(RigidJkrMixedModeError::UnsupportedTorsion { .. })
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(transfer(&pair, &cancelled).is_err());

    let retry = transfer(&pair, &gate).unwrap();
    let again = transfer(&pair, &gate).unwrap();
    assert_eq!(retry.previous_state(), previous);
    assert_eq!(retry.bodies(), &f.bodies);
    assert_eq!(retry.bodies(), again.bodies());
    assert_eq!(retry.state(), again.state());
    assert_eq!(retry.report(), again.report());
}

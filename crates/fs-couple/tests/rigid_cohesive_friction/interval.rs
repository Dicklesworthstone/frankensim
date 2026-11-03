//! Native interval stopping, complete path accounts and outer rollback.
use super::*;

pub(super) fn controls() -> RigidCohesiveIntervalConfig {
    RigidCohesiveIntervalConfig {
        minimum_step_s: 1e-14,
        reversal_time_tolerance_s: 1e-12,
        maximum_trials: 128,
        maximum_substeps: 32,
        maximum_total_evaluations: 2_000_000,
        maximum_total_material_evaluations: 12_000_000,
        maximum_total_iterations: 12_000,
    }
}

// Finite, very large rotational inertias isolate the analytic translation limit
// without suppressing native rotational equations or their actual work.
pub(super) fn translation_fixture<'a>(
    bond: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
    incoming: f64,
) -> RigidCohesiveSystem<'a> {
    let seed = system(
        bond,
        None,
        failed(bond),
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, Vec3::new(0.0, incoming, 0.0)],
        zero(),
        0.02,
        config(),
    );
    let bodies = core::array::from_fn(|i| {
        let old = &seed.bodies()[i];
        let properties = MassProperties::from_principal_frame(
            old.properties().mass(),
            Vec3::ZERO,
            Vec3::new(1e6, 1e6, 1e6),
            UnitQuaternion::IDENTITY,
        )
        .unwrap();
        RigidCohesiveBody::new(
            old.material().clone(),
            properties,
            old.state(),
            old.limits(),
        )
        .unwrap()
    });
    let gate = CancelGate::new();
    let mut owner = RigidCohesiveSystem::new(
        bodies,
        seed.sites().to_vec(),
        seed.gravity(),
        seed.step_duration_s(),
        config(),
        &gate,
    )
    .unwrap();
    owner
        .enable_damaged_area_friction(
            vec![RigidDamageFrictionSite {
                source: friction,
                activation: ReclosureActivation::AtRest,
            }],
            friction_config(),
            &gate,
        )
        .unwrap();
    owner
}

#[test]
fn unaligned_native_stop_preserves_current_temperature_then_resticks_and_breaks_away() {
    let card = source_card([0.6, 0.3], false, Some([0.4, 0.1]));
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    let bond = resolve_bilinear_cohesive(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
    )
    .unwrap();
    let residual =
        resolve_residual_coulomb(&card, &query, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let incoming = 3.7e-4;
    let reduced_mass = 1.0 / (1.0 / 0.4 + 1.0 / 0.7);
    for hot in [false, true] {
        let mut owner = translation_fixture(&bond, &residual, incoming);
        let pressure = if hot { 0.6 } else { 0.4 };
        let friction_force = pressure * if hot { 0.1 } else { 0.3 };
        if hot {
            let jump = v(owner.sites()[0].point().numerical().material().jump_m());
            let loads = balanced_load(
                Vec3::new(-pressure, friction_force, 0.0),
                jump.add(Vec3::new(0.0, incoming * 0.01, 0.0)),
            );
            owner
                .step_with_temperature_and_friction(loads, &[&hot_bond], &[&hot_residual], "T")
                .unwrap();
        }
        let jump_before = owner.sites()[0].point().numerical().material().jump_m();
        let time_before = owner.time_s();
        let mut loads = balanced_load(Vec3::new(-pressure, 0.0, 0.0), v(jump_before));
        // No imposed shear or projection: relative kinetic energy is dissipated
        // by the original finite-pressure Coulomb forces on both finite masses.
        let interval = owner.advance_interval(loads, 0.02, controls()).unwrap();
        assert!(interval.attempted_substeps > interval.accepted_substeps);
        assert!(interval.accepted_substeps >= 2);
        // At least one rejected coarse trial is charged at its full native
        // allowance; discarded candidate work is not hidden by successful stops.
        assert!(interval.charged_evaluations >= config().maximum_evaluations);
        assert!(interval.charged_material_evaluations >= config().maximum_material_evaluations);
        assert!(interval.charged_iterations >= friction_config().maximum_iterations);
        assert!(interval.charged_evaluations <= controls().maximum_total_evaluations);
        assert!(
            interval.charged_material_evaluations <= controls().maximum_total_material_evaluations
        );
        assert!(interval.charged_iterations <= controls().maximum_total_iterations);
        close(interval.time_s, time_before + 0.02, 0.0);
        let slip = owner.sites()[0].point().numerical().material().jump_m()[1] - jump_before[1];
        close(
            slip,
            0.5 * reduced_mass * incoming * incoming / friction_force,
            2e-10,
        );
        close(
            interval.friction_dissipation_j,
            0.5 * reduced_mass * incoming * incoming,
            2e-12,
        );
        close(
            interval.sites[0].friction_dissipation_j,
            interval.friction_dissipation_j,
            1e-20,
        );
        assert_eq!(interval.fracture_dissipation_j, 0.0);
        assert_eq!(interval.viscous_dissipation_j, 0.0);
        assert!(interval.energy_residual_j.abs() <= interval.energy_tolerance_j);
        assert!(
            norm(interval.linear_momentum_residual_n_s) <= interval.linear_momentum_tolerance_n_s
        );
        assert_eq!(
            owner.damaged_area_friction_states().unwrap()[0].occupancy(),
            DamageFrictionRegime::Sticking
        );
        let expected_bond = if hot { &hot_bond } else { &bond };
        let expected_friction = if hot { &hot_residual } else { &residual };
        assert!(core::ptr::eq(
            owner.sites()[0].point().source(),
            expected_bond
        ));
        assert!(core::ptr::eq(
            owner.damaged_area_friction_sources().unwrap()[0].source,
            expected_friction
        ));
        // .20 N exceeds both kinetic disks but stays inside both static disks.
        let jump = v(owner.sites()[0].point().numerical().material().jump_m());
        loads = balanced_load(Vec3::new(-pressure, 0.2, 0.0), jump);
        let held = owner.advance_interval(loads, 0.007, controls()).unwrap();
        assert_eq!(held.friction_dissipation_j, 0.0);
        close(
            owner.next_step_time_s().unwrap(),
            owner.time_s() + 0.02,
            0.0,
        );
        let held_clock = owner.time_s();
        let ordinary = *owner.step(loads).unwrap();
        close(ordinary.time_s, held_clock + 0.02, 0.0);
        assert_eq!(ordinary.friction_dissipation_j, 0.0);
        let jump = v(owner.sites()[0].point().numerical().material().jump_m());
        let broken = owner
            .advance_interval(
                balanced_load(Vec3::new(-pressure, 0.3, 0.0), jump),
                0.002,
                controls(),
            )
            .unwrap();
        assert!(broken.friction_dissipation_j > 0.0);
        assert_eq!(
            owner.damaged_area_friction_states().unwrap()[0].occupancy(),
            DamageFrictionRegime::Sliding
        );
    }
}

pub(super) fn rotating_fixture<'a>(
    bond: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
) -> RigidCohesiveSystem<'a> {
    system(
        bond,
        Some(friction),
        partial(bond),
        &[
            [0.0, -0.002, 0.001],
            [0.0, 0.001, -0.002],
            [0.0, 0.002, 0.002],
        ],
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), 0.8).unwrap(),
        [Vec3::new(0.0, -0.025, 0.015), Vec3::new(0.0, 0.020, -0.010)],
        [Vec3::new(0.3, -0.2, 0.4), Vec3::new(-0.1, 0.15, 0.2)],
        1e-4,
        config(),
    )
}

#[test]
fn interval_retains_spread_rotating_sites_and_sums_original_native_paths() {
    let (bond, friction) = sources([0.6, 0.3], true);
    let mut interval_owner = rotating_fixture(&bond, &friction);
    let mut ordinary = rotating_fixture(&bond, &friction);
    let mut sums = [RigidCohesiveIntervalSite::default(); 3];
    for (sum, site) in sums.iter_mut().zip(ordinary.sites()) {
        sum.maximum_normal_jump_m = site.point().numerical().material().jump_m()[0];
    }
    for _ in 0..2 {
        ordinary.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        for (sum, site) in sums.iter_mut().zip(ordinary.site_frames()) {
            let site = site.unwrap();
            sum.fracture_dissipation_j += site.fracture_dissipation_j;
            sum.viscous_dissipation_j += site.viscous_dissipation_j;
            sum.friction_dissipation_j += site.friction.unwrap().friction_dissipation_j;
            sum.cohesive_spatial_work_j += site.spatial_work_j;
            sum.maximum_normal_jump_m = sum.maximum_normal_jump_m.max(site.jump_m[0]);
        }
    }
    let interval = interval_owner
        .advance_interval([RigidCohesiveLoad::ZERO; 2], 2e-4, controls())
        .unwrap();
    assert_eq!(interval.accepted_substeps, 2);
    assert_eq!(interval.attempted_substeps, 2);
    assert_eq!(interval_owner.bodies(), ordinary.bodies());
    assert_eq!(interval_owner.sites(), ordinary.sites());
    assert_eq!(
        interval_owner.damaged_area_friction_states(),
        ordinary.damaged_area_friction_states()
    );
    assert_eq!(interval.sites, sums);
    assert!(interval.friction_dissipation_j > 0.0 && interval.viscous_dissipation_j > 0.0);
    assert!(
        norm(interval.angular_momentum_residual_n_m_s) <= interval.angular_momentum_tolerance_n_m_s
    );
    close(
        ordinary.time_s(),
        ordinary.steps() as f64 * ordinary.step_duration_s(),
        0.0,
    );
}

#[test]
fn late_interval_budget_and_cancel_preserve_complete_native_state_and_retry() {
    let (bond, friction) = sources([0.6, 0.3], true);
    let mut owner = rotating_fixture(&bond, &friction);
    let mut control = rotating_fixture(&bond, &friction);
    let loads = [RigidCohesiveLoad::ZERO; 2];
    owner.step(loads).unwrap();
    control.step(loads).unwrap();
    let bodies = owner.bodies().clone();
    let sites = owner.sites().to_vec();
    let states = owner.damaged_area_friction_states().unwrap().to_vec();
    let frames = owner.site_frames().to_vec();
    let last = owner.last_frame().copied();
    let next_time = owner.next_step_time_s().unwrap();
    let limits = RigidCohesiveIntervalConfig {
        maximum_substeps: 1,
        ..controls()
    };
    assert_eq!(
        owner.advance_interval(loads, 2e-4, limits),
        Err(RigidCohesiveError::Budget("rigid interval substeps"))
    );
    let work_limit = RigidCohesiveIntervalConfig {
        maximum_total_evaluations: 1,
        ..controls()
    };
    assert!(matches!(
        owner.advance_interval(loads, 2e-4, work_limit),
        Err(RigidCohesiveError::Budget(_))
    ));
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        owner.advance_interval_under_gate(loads, 2e-4, controls(), &gate),
        Err(RigidCohesiveError::Cancelled)
    );
    assert_eq!(owner.bodies(), &bodies);
    assert_eq!(owner.sites(), sites);
    assert_eq!(owner.damaged_area_friction_states().unwrap(), states);
    assert_eq!(owner.site_frames(), frames);
    assert_eq!(owner.last_frame().copied(), last);
    assert_eq!(owner.next_step_time_s().unwrap(), next_time);
    assert_eq!(owner.config(), config());
    assert_eq!(owner.step_duration_s(), 1e-4);
    let retry = owner.advance_interval(loads, 2e-4, controls()).unwrap();
    let first = control.advance_interval(loads, 2e-4, controls()).unwrap();
    assert_eq!(retry, first);
    assert_eq!(owner.bodies(), control.bodies());
    assert_eq!(owner.sites(), control.sites());
}

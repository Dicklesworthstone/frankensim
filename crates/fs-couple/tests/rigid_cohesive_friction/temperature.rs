//! G1/G3/G4: paired temperature sources on the original rigid Coulomb mixture.
use super::*;

fn cold(card: &InterfaceSystemCard) -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    (
        resolve_bilinear_cohesive(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            0.05,
        )
        .unwrap(),
        resolve_residual_coulomb(card, &query, MaterialPropertySelection::SingleClaimOnly).unwrap(),
    )
}

fn check_sources(
    owner: &RigidCohesiveSystem<'_>,
    cohesive: &ResolvedBilinearCohesive,
    friction: &ResidualCoulombSource,
) {
    for ((site, state), source) in owner
        .sites()
        .iter()
        .zip(owner.damaged_area_friction_states().unwrap())
        .zip(owner.damaged_area_friction_sources().unwrap())
    {
        assert!(core::ptr::eq(site.point().source(), cohesive));
        assert!(core::ptr::eq(source.source, friction));
        assert!(core::ptr::eq(state.point().law(), cohesive.law()));
        assert_eq!(state.point(), site.point().numerical().material());
    }
}

#[test]
fn residual_temperature_refresh_pins_claims_and_rejects_invalid_coefficient_order() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot = residual.at_temperature(&card, "T", 320.0).unwrap();
    close(hot.law().static_coefficient().unwrap(), 0.3, 1e-15);
    close(
        hot.law().zero_speed_kinetic_coefficient().unwrap(),
        0.1,
        1e-15,
    );
    residual.require_temperature_successor(&hot, "T").unwrap();
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    hot.require_cohesive_source(&hot_bond).unwrap();
    assert!(hot.require_cohesive_source(&bond).is_err());
    for (old, new) in residual
        .resolved()
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(old.requirement(), new.requirement());
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    assert_eq!(
        residual.resolved().card_identity(),
        hot.resolved().card_identity()
    );
    assert_eq!(residual.surfaces(), hot.surfaces());
    for (axis, value) in [("T", 321.0), ("T", 0.0), ("pressure", 300.0)] {
        assert!(residual.at_temperature(&card, axis, value).is_err());
    }
    let foreign_card = source_card([0.7, 0.3], false, Some([0.3, 0.1]));
    assert!(residual.at_temperature(&foreign_card, "T", 320.0).is_err());
    let (_, foreign) = cold(&foreign_card);
    assert!(
        residual
            .require_temperature_successor(&foreign, "T")
            .is_err()
    );
    // A valid cold source may become invalid at another point in its table.
    // The new static and kinetic values must be admitted together.
    let crossed = source_card([0.6, 0.3], false, Some([0.1, 0.2]));
    let (_, valid_cold) = cold(&crossed);
    assert!(valid_cold.at_temperature(&crossed, "T", 320.0).is_err());
}

#[test]
fn failed_interface_holds_new_static_disks_through_heating_and_cooling_without_heat() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let point = failed(&bond);
    let mut owner = system(
        &bond,
        Some(&residual),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let before = owner.bodies().clone();
    let energy = owner.total_energy_j().unwrap();
    let original_fracture = point.dissipated_j_m2();
    // Hot: |T|=.15 < mu_s(T)*Kc(T)*|g|*A=.18, above its kinetic .06.
    let force = Vec3::new(-0.6, 0.12, -0.09);
    let frame = *owner
        .step_with_temperature_and_friction(
            balanced_load(force, v(point.jump_m())),
            &[&hot_bond],
            &[&hot_residual],
            "T",
        )
        .unwrap();
    close(frame.parameter_work_j, 2e-5, 2e-14);
    close(frame.thermal_stored_change_j, 2e-5, 2e-14);
    assert_eq!(frame.thermal_damage_dissipation_j, 0.0);
    assert_eq!(frame.fracture_dissipation_j, 0.0);
    assert_eq!(frame.friction_dissipation_j, 0.0);
    let actual = owner.site_frames()[0].unwrap().friction.unwrap();
    assert_eq!(actual.regime, DamageFrictionRegime::Sticking);
    close(actual.damaged_compression_force_n, 0.6, 2e-9);
    close(actual.local_force_n[0], 0.12, 2e-9);
    close(actual.local_force_n[1], -0.09, 2e-9);
    check_sources(&owner, &hot_bond, &hot_residual);
    // Cold: |T|=.2 fits the restored static .24, while stale hot mu would
    // admit only .12 at the restored compression stiffness.
    let frame = *owner
        .step_with_temperature_and_friction(
            balanced_load(Vec3::new(-0.4, 0.16, -0.12), v(point.jump_m())),
            &[&bond],
            &[&residual],
            "T",
        )
        .unwrap();
    close(frame.parameter_work_j, -2e-5, 2e-13);
    assert_eq!(frame.friction_dissipation_j, 0.0);
    check_sources(&owner, &bond, &residual);
    close(owner.total_energy_j().unwrap(), energy, 2e-12);
    for i in 0..2 {
        close_vec(
            owner.bodies()[i].state().pose().position_world(),
            before[i].state().pose().position_world(),
            2e-11,
        );
        close_vec(
            owner.bodies()[i].state().linear_momentum_world(),
            before[i].state().linear_momentum_world(),
            2e-10,
        );
    }
    let retained = owner.damaged_area_friction_states().unwrap()[0];
    assert_eq!(retained.point().damage(), 1.0);
    assert_eq!(retained.point().dissipated_j_m2(), original_fracture);
    assert_eq!(retained.occupancy(), DamageFrictionRegime::Sticking);
}

#[test]
fn fixed_jump_thermal_damage_activates_only_the_newly_damaged_friction_area() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let jump = [-2e-4, 1.8e-4, 0.0];
    let point = bond.law().bonded_point().trial(jump, 1e-8).unwrap().point;
    assert_eq!(point.damage(), 0.0);
    let mut owner = system(
        &bond,
        Some(&residual),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let original_energy = owner.total_energy_j().unwrap();
    let damage = (0.00125 / (0.00125 - 0.00016)) * (1.0 - 0.00016 / 0.00018);
    let thermal_release = AREA * 8.1 * damage;
    let force = Vec3::new(-0.6, (1.0 - damage) * 0.09 + 0.02, 0.0);
    let frame = *owner
        .step_with_temperature_and_friction(
            balanced_load(force, v(jump)),
            &[&hot_bond],
            &[&hot_residual],
            "T",
        )
        .unwrap();
    // At fixed jump, compression adds 20 J/m² while the intact shear
    // stiffness change removes 8.1 J/m². New damage releases 8.1*d J/m².
    close(frame.parameter_work_j, 11.9e-6, 2e-13);
    close(frame.thermal_damage_dissipation_j, thermal_release, 2e-13);
    close(frame.fracture_dissipation_j, thermal_release, 2e-12);
    assert_eq!(frame.friction_dissipation_j, 0.0);
    let local = owner.site_frames()[0].unwrap().friction.unwrap();
    assert_eq!(local.regime, DamageFrictionRegime::Sticking);
    close(local.damaged_compression_force_n, 0.6 * damage, 2e-9);
    close(local.local_force_n[0], 0.02, 2e-9);
    let state = owner.damaged_area_friction_states().unwrap()[0];
    close(state.point().damage(), damage, 2e-10);
    close(
        AREA * state.point().dissipated_j_m2(),
        thermal_release,
        2e-12,
    );
    close_vec(v(state.point().jump_m()), v(jump), 2e-11);
    close(
        owner.total_energy_j().unwrap(),
        original_energy + frame.parameter_work_j + frame.external_work_j
            - frame.fracture_dissipation_j,
        2e-10,
    );
    check_sources(&owner, &hot_bond, &hot_residual);
}

#[test]
fn thermal_vector_sliding_keeps_damage_and_separates_parameter_work_from_actual_loss() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let point = partial(&bond);
    let delta = Vec3::new(0.0, 6e-5, 8e-5);
    let dt = 0.01;
    let mut owner = system(
        &bond,
        Some(&residual),
        point,
        &[[0.0; 3]; 2],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let mut accounted = owner.total_energy_j().unwrap();
    // Independent constant-damage calculation: d=.625, Kt_hot=5e8,
    // mean bond shear=(1-d)*Kt*delta/2 and tau_f=mu_k*d*Kc*|g_n|.
    let hot_force = Vec3::new(-0.6, 0.028125, 0.0375);
    let first = *owner
        .step_with_temperature_and_friction(
            balanced_load(hot_force, v(point.jump_m()).add(delta.scale(0.5))),
            &[&hot_bond; 2],
            &[&hot_residual; 2],
            "T",
        )
        .unwrap();
    close(first.parameter_work_j, 2e-5, 2e-13);
    close(first.friction_dissipation_j, 3.75e-6, 2e-11);
    close_vec(first.contact_loads[0].force_world_n, hot_force, 3e-8);
    assert_eq!(first.fracture_dissipation_j, 0.0);
    assert_eq!(first.viscous_dissipation_j, 0.0);
    assert_eq!(first.material_evaluations, 2 * (first.evaluations + 2));
    accounted += first.parameter_work_j + first.external_work_j - first.friction_dissipation_j;
    close(owner.total_energy_j().unwrap(), accounted, 2e-10);
    check_sources(&owner, &hot_bond, &hot_residual);
    for (site, history) in owner
        .site_frames()
        .iter()
        .zip(owner.damaged_area_friction_states().unwrap())
    {
        let site = site.unwrap();
        close_vec(v(site.jump_m), v(point.jump_m()).add(delta), 2e-10);
        let f = site.friction.unwrap();
        assert_eq!(f.regime, DamageFrictionRegime::Sliding);
        close(f.damaged_compression_force_n, 0.1875, 2e-9);
        close(f.local_force_n[0], 0.01125, 2e-9);
        close(f.local_force_n[1], 0.015, 2e-9);
        close(history.point().damage(), point.damage(), 2e-14);
        close(
            history.point().dissipated_j_m2(),
            point.dissipated_j_m2(),
            2e-12,
        );
        assert_eq!(history.occupancy(), DamageFrictionRegime::Sliding);
    }
    // Cooling at nonzero shear adds retained-bond storage while reducing
    // compression storage: A[-20 + .5*(1-d)*(1e9-5e8)*|delta|^2].
    let cold_force = Vec3::new(-0.4, 0.07875, 0.105);
    let second = *owner
        .step_with_temperature_and_friction(
            balanced_load(cold_force, v(point.jump_m()).add(delta.scale(1.5))),
            &[&bond; 2],
            &[&residual; 2],
            "T",
        )
        .unwrap();
    close(second.parameter_work_j, -19.0625e-6, 2e-11);
    close(second.friction_dissipation_j, 7.5e-6, 3e-11);
    close_vec(second.contact_loads[0].force_world_n, cold_force, 4e-8);
    close(second.thermal_damage_dissipation_j, 0.0, 1e-15);
    close(second.fracture_dissipation_j, 0.0, 1e-15);
    accounted += second.parameter_work_j + second.external_work_j
        - second.friction_dissipation_j
        - second.fracture_dissipation_j;
    close(owner.total_energy_j().unwrap(), accounted, 4e-10);
    check_sources(&owner, &bond, &residual);
    for history in owner.damaged_area_friction_states().unwrap() {
        close(history.point().damage(), point.damage(), 2e-14);
        assert_eq!(history.occupancy(), DamageFrictionRegime::Sliding);
        close_vec(
            v(history.point().jump_m()),
            v(point.jump_m()).add(delta.scale(2.0)),
            3e-10,
        );
    }
}

#[test]
fn partial_source_and_cancel_refusals_preserve_accepted_histories_and_identical_retries() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let foreign_card = source_card([0.7, 0.3], false, Some([0.3, 0.1]));
    let (_, foreign) = cold(&foreign_card);
    let point = partial(&bond);
    let make = || {
        system(
            &bond,
            Some(&residual),
            point,
            &[[0.0; 3]; 2],
            UnitQuaternion::IDENTITY,
            zero(),
            zero(),
            0.002,
            config(),
        )
    };
    let mut owner = make();
    let mut control = make();
    let cold_load = balanced_load(Vec3::new(-0.4, 0.0, 0.0), v(point.jump_m()));
    owner.step(cold_load).unwrap();
    control.step(cold_load).unwrap();
    let bodies = owner.bodies().clone();
    let states = owner.damaged_area_friction_states().unwrap().to_vec();
    let sites = owner.sites().to_vec();
    let frames = owner.site_frames().to_vec();
    let frame = *owner.last_frame().unwrap();
    let hot_load = balanced_load(Vec3::new(-0.6, 0.0, 0.0), v(point.jump_m()));
    for friction in [
        &[&hot_residual, &foreign][..],
        &[&hot_residual, &residual][..],
        &[&hot_residual][..],
    ] {
        assert!(
            owner
                .step_with_temperature_and_friction(hot_load, &[&hot_bond; 2], friction, "T",)
                .is_err()
        );
        assert_eq!(owner.bodies(), &bodies);
        assert_eq!(owner.sites(), sites);
        assert_eq!(owner.damaged_area_friction_states().unwrap(), states);
        assert_eq!(owner.site_frames(), frames);
        assert_eq!(*owner.last_frame().unwrap(), frame);
        assert_eq!(owner.steps(), 1);
        check_sources(&owner, &bond, &residual);
    }
    let gate = CancelGate::new();
    gate.request();
    assert!(
        owner
            .step_with_temperature_and_friction_under_gate(
                hot_load,
                &[&hot_bond; 2],
                &[&hot_residual; 2],
                "T",
                &gate,
            )
            .is_err()
    );
    assert_eq!(owner.bodies(), &bodies);
    assert_eq!(owner.damaged_area_friction_states().unwrap(), states);
    // Ordinary evaluation must ignore every partially staged thermal source.
    assert_eq!(
        *owner.step(cold_load).unwrap(),
        *control.step(cold_load).unwrap()
    );
    let a = *owner
        .step_with_temperature_and_friction(hot_load, &[&hot_bond; 2], &[&hot_residual; 2], "T")
        .unwrap();
    let b = *control
        .step_with_temperature_and_friction(hot_load, &[&hot_bond; 2], &[&hot_residual; 2], "T")
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(owner.bodies(), control.bodies());
    assert_eq!(
        owner.damaged_area_friction_states(),
        control.damaged_area_friction_states()
    );
    check_sources(&owner, &hot_bond, &hot_residual);
}

#[test]
fn one_transfer_per_site_shares_native_budget_and_late_refusal_restores_cold_solver() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold(&card);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let cfg = RigidCohesiveConfig {
        maximum_material_evaluations: 4,
        ..config()
    };
    let make = || {
        system(
            &bond,
            Some(&residual),
            bond.law().bonded_point(),
            &[[0.0; 3]; 2],
            UnitQuaternion::IDENTITY,
            zero(),
            zero(),
            0.002,
            cfg,
        )
    };
    let mut owner = make();
    let mut control = make();
    let before = owner.bodies().clone();
    let states = owner.damaged_area_friction_states().unwrap().to_vec();
    assert!(matches!(
        owner.step_with_temperature_and_friction(
            [RigidCohesiveLoad::ZERO; 2],
            &[&hot_bond; 2],
            &[&hot_residual; 2],
            "T"
        ),
        Err(RigidCohesiveError::Budget(
            "original friction audit evaluations"
        ))
    ));
    assert_eq!(owner.bodies(), &before);
    assert_eq!(owner.damaged_area_friction_states().unwrap(), states);
    assert_eq!(owner.steps(), 0);
    check_sources(&owner, &bond, &residual);
    let a = *owner.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    let b = *control.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.material_evaluations, 4);
    assert_eq!(a.parameter_work_j, 0.0);
    assert_eq!(
        owner.damaged_area_friction_states(),
        control.damaged_area_friction_states()
    );
}

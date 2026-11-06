//! Original DMT temperatures on actual native normal/shear body motion.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_contact::normal_patch::dmt::DmtRegime;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::adhesive::normal_motion::temperature::RigidDmtSourceRefreshError;

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod thermal_fixture;

fn temperatures(t: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0 + 0.5 * (t - 300.0),
        surface_b_k: 300.0 - 0.25 * (t - 300.0),
        interface_k: t,
    }
}
fn sources(f: &thermal_fixture::Fixture) -> (SourceDmtSphereContact, SourceAdhesiveShear) {
    let normal = f.dmt(common::Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (normal, shear)
}
fn pair<'a>(
    f: &Fixture,
    normal: &'a SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    dt: f64,
) -> PreparedRigidDmtPair<'a> {
    PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        f.contact(normal, shear),
        dt,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn all_original_phases_refresh_at_fixed_native_state_and_drive_the_successor_step() {
    let thermal = thermal_fixture::Fixture::with_options(1e5, 0.01, true, 320.0);
    let (normal, shear) = sources(&thermal);
    let dt = 5e-10;
    for delta in [
        0.2 * common::Z,
        0.0,
        -0.5 * common::Z,
        -common::Z,
        -2.0 * common::Z,
    ] {
        // Central reaction points preserve exact piecewise corners and permit
        // genuine attached translation without an unrepresented torsion mode.
        let mut f = fixture(delta, [Vec3::ZERO; 2]);
        f.set_a_velocity(Vec3::new(0.25, 0.0, 0.0), 0.0);
        for t in [300.0, 285.0, 315.0] {
            let mut original = pair(&f, &normal, &shear, dt);
            assert_eq!(original.approach_m(), delta);
            let state = original
                .initialize([0.3 * common::Z, -0.1 * common::Z])
                .unwrap();
            let refresh = normal
                .refresh_at_held_indentation(
                    original.approach_m(),
                    thermal.cards(),
                    temperatures(t),
                )
                .unwrap();
            let next_shear = shear.rebind_dmt_temperature(&normal, &refresh).unwrap();
            let (mut next, transferred, report) = original
                .with_temperature_refresh(
                    state,
                    &refresh,
                    &next_shear,
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            assert_eq!(original.contact().normal.identity(), normal.identity());
            assert_eq!(original.contact().shear, &shear);
            assert!(core::ptr::eq(next.contact().normal, refresh.source()));
            assert!(core::ptr::eq(next.contact().shear, &next_shear));
            let binding = next.contact().binding;
            let original_binding = original.contact().binding;
            assert_eq!(
                (binding.body_a, binding.body_b),
                (original_binding.body_a, original_binding.body_b)
            );
            assert_eq!(
                [binding.arm_a_body, binding.arm_b_body],
                [original_binding.arm_a_body, original_binding.arm_b_body]
            );
            assert_eq!(next.contact().frame, original.contact().frame);
            assert_eq!(
                next.contact().reference_approach_m,
                original.contact().reference_approach_m
            );
            assert_eq!(next.mobility_per_kg(), original.mobility_per_kg());
            assert_eq!(next.compliance_m_per_n(), original.compliance_m_per_n());
            assert_eq!(next.reaction_arms_body(), original.reaction_arms_body());
            assert_eq!(next.approach_m(), delta);
            assert_eq!(transferred.position_m(), state.position_m());
            assert_eq!(transferred.regime(), state.regime());
            assert_eq!(report.position_m, state.position_m());
            assert_eq!(report.regime, state.regime());
            assert_eq!(report.normal_before, *refresh.before());
            assert_eq!(report.normal_after, *refresh.response());
            assert_eq!(report.geometry_residual_m, 0.0);
            assert_eq!(report.parameter_work_j, refresh.parameter_work_j());
            close(
                report.normal_storage_change_j,
                refresh.potential_change_j(),
                1e-30,
            );
            assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
            assert_eq!(report.energy_tolerance_j, controls().energy_tolerance_j);
            assert_eq!(refresh.mechanical_work_j(), 0.0);
            assert_eq!(refresh.before().regime, refresh.response().regime);
            if delta <= 0.0 {
                assert_eq!(transferred.regime(), DmtCoupledRegime::Released);
                assert_eq!(report.normal_after.contact_area_m2, 0.0);
                assert_eq!(report.normal_after.repulsive_normal_force_n, 0.0);
            } else {
                assert_eq!(transferred.regime(), DmtCoupledRegime::Sliding);
            }
            if delta <= -common::Z {
                assert_eq!(report.normal_after.regime, DmtRegime::Separated);
                assert_eq!(report.parameter_work_j, 0.0);
            }
            // Compare the complete native solve with independently binding
            // the same actual bodies to the original successor source owners.
            let mut direct = pair(&f, refresh.source(), &next_shear, dt);
            let initial = direct
                .initialize([state.position_m()[1], state.position_m()[2]])
                .unwrap();
            assert_eq!(transferred, initial);
            let wrenches = [
                Wrench {
                    force_world: Vec3::new(1e-9, 0.0, 2e-8),
                    torque_body: Vec3::new(0.0, 2e-18, 0.0),
                },
                Wrench::ZERO,
            ];
            let advanced = next
                .advance_with_wrenches(transferred, wrenches, &CancelGate::new_clock_free())
                .unwrap();
            let expected = direct
                .advance_with_wrenches(initial, wrenches, &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(advanced, expected);
            assert!(
                advanced.report.audit.total_energy_residual_j.abs()
                    <= controls().energy_tolerance_j
            );
            if delta > 0.0 {
                assert_eq!(advanced.report.contact.regime, DmtCoupledRegime::Sliding);
                assert!(advanced.report.dissipated_energy_j > 0.0);
            } else if advanced.state_after.position_m()[0] <= 0.0 {
                assert_eq!(advanced.report.dissipated_energy_j, 0.0);
                assert_eq!(advanced.report.applied_mean_force_n[1..], [0.0; 2]);
            }
            let unchanged = original
                .advance_with_wrenches(state, wrenches, &CancelGate::new_clock_free())
                .unwrap();
            if t == 300.0 {
                assert_eq!(report.parameter_work_j, 0.0);
                assert_eq!(advanced, unchanged);
            } else if delta > -common::Z {
                assert_ne!(
                    advanced.report.applied_mean_force_n[0].to_bits(),
                    unchanged.report.applied_mean_force_n[0].to_bits()
                );
            }
        }
    }
}

#[test]
fn noncentral_mobility_and_small_thermal_work_survive_a_large_common_body_motion() {
    let thermal = thermal_fixture::Fixture::with_options(1e5, 0.01, false, 320.0);
    let (normal, shear) = sources(&thermal);
    let dt = 5e-10;
    let mut f = fixture(
        -0.5 * common::Z,
        [
            Vec3::new(0.2, 0.1, 0.4).scale(LENGTH),
            Vec3::new(-0.12, 0.15, -0.3).scale(LENGTH),
        ],
    );
    f.set_a_velocity(Vec3::new(0.12, -0.04, 0.02), 0.0);
    let original = pair(&f, &normal, &shear, dt);
    let state = original.initialize([0.0; 2]).unwrap();
    let refresh = normal
        .refresh_at_held_indentation(original.approach_m(), thermal.cards(), temperatures(315.0))
        .unwrap();
    let next_shear = shear.rebind_dmt_temperature(&normal, &refresh).unwrap();
    let (mut successor, history, report) = original
        .with_temperature_refresh(state, &refresh, &next_shear, &CancelGate::new_clock_free())
        .unwrap();
    let w = original.mobility_per_kg();
    assert!(w[0][1].abs() > 1e12 && w[0][2].abs() > 1e12 && w[1][2].abs() > 1e12);
    assert_eq!(successor.mobility_per_kg(), w);
    let mut direct = pair(&f, refresh.source(), &next_shear, dt);
    let direct_state = direct.initialize([0.0; 2]).unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.2, -0.1, 0.3).scale(1e-8),
            torque_body: Vec3::new(0.4, 0.2, -0.3).scale(1e-15),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.4, -0.1).scale(1e-8),
            torque_body: Vec3::new(-0.2, 0.1, 0.4).scale(1e-15),
        },
    ];
    let a = successor
        .advance_with_wrenches(history, wrenches, &CancelGate::new_clock_free())
        .unwrap();
    let b = direct
        .advance_with_wrenches(direct_state, wrenches, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(a, b);
    assert!(a.state_after.position_m()[0] < 0.0);
    assert_eq!(a.report.dissipated_energy_j, 0.0);
    assert_eq!(a.report.contact.path.mean_real_area_m2, 0.0);

    let boost = Vec3::new(0.0, 2.0_f64.powi(30), 0.0);
    for i in 0..2 {
        f.bodies[i] = RigidBodyState::new(
            f.bodies[i].pose(),
            f.bodies[i]
                .linear_momentum_world()
                .add(boost.scale(f.properties[i].mass())),
            f.bodies[i].angular_momentum_body(),
        )
        .unwrap();
    }
    let captured = f.bodies;
    let mut boosted = pair(&f, &normal, &shear, dt);
    let mut boosted_history = boosted.initialize([0.0; 2]).unwrap();
    let updated = boosted
        .apply_temperature_refresh(
            &mut boosted_history,
            &refresh,
            &next_shear,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(updated, report);
    assert_eq!(f.bodies, captured);
    assert_eq!(boosted_history, history);
    assert_ne!(report.parameter_work_j, 0.0);
    let kinetic = kinetic(f.bodies, f.properties);
    assert_eq!(kinetic + report.parameter_work_j, kinetic);
    assert_eq!(updated.energy_tolerance_j, controls().energy_tolerance_j);
}

#[test]
fn foreign_stale_held_load_and_late_receiver_refusals_preserve_both_owners() {
    let thermal = thermal_fixture::Fixture::with_options(1e5, 0.01, false, 320.0);
    let (normal, shear) = sources(&thermal);
    let dt = 5e-10;
    let mut f = fixture(0.2 * common::Z, [Vec3::ZERO; 2]);
    f.set_a_velocity(Vec3::new(0.25, 0.0, 0.0), 0.0);
    let mut owner = pair(&f, &normal, &shear, dt);
    let mut history = owner.initialize([0.0; 2]).unwrap();
    let initial = history;
    let refresh = normal
        .refresh_at_held_indentation(owner.approach_m(), thermal.cards(), temperatures(315.0))
        .unwrap();
    let next_shear = shear.rebind_dmt_temperature(&normal, &refresh).unwrap();
    let stale = normal
        .refresh_at_held_indentation(0.3 * common::Z, thermal.cards(), temperatures(315.0))
        .unwrap();
    let foreign = SourceDmtSphereContact::new(
        normal.normal_state().clone(),
        1.1 * normal.effective_radius_m(),
        normal.interaction_range_m(),
        normal.envelope(),
    )
    .unwrap();
    let foreign = foreign
        .refresh_at_held_indentation(owner.approach_m(), thermal.cards(), temperatures(315.0))
        .unwrap();
    let force = normal.evaluate(owner.approach_m()).unwrap().normal_force_n;
    let held = normal
        .refresh_at_held_load(force, thermal.cards(), temperatures(315.0))
        .unwrap();
    let held_shear = shear.rebind_dmt_temperature(&normal, &held).unwrap();
    let mobility = owner.mobility_per_kg();
    for (invalid, junction) in [
        (&stale, &next_shear),
        (&foreign, &next_shear),
        (&held, &held_shear),
        (&refresh, &shear),
    ] {
        assert!(
            owner
                .apply_temperature_refresh(
                    &mut history,
                    invalid,
                    junction,
                    &CancelGate::new_clock_free(),
                )
                .is_err()
        );
        assert_eq!(history, initial);
        assert_eq!(owner.contact().normal.identity(), normal.identity());
        assert_eq!(owner.contact().shear, &shear);
        assert_eq!(owner.mobility_per_kg(), mobility);
    }
    let refused: Result<(), RigidDmtSourceRefreshError> = owner
        .apply_temperature_refresh_with_admission(
            &mut history,
            &refresh,
            &next_shear,
            &CancelGate::new_clock_free(),
            |report, candidate, state| {
                assert_eq!(report.regime, DmtCoupledRegime::Sliding);
                assert_eq!(state.position_m(), initial.position_m());
                assert_eq!(
                    candidate.contact().normal.identity(),
                    refresh.source().identity()
                );
                Err(RigidDmtError::Invalid("receiving thermal owner refused").into())
            },
        );
    assert!(refused.is_err());
    let cancelled = CancelGate::new_clock_free();
    let result: Result<(), RigidDmtSourceRefreshError> = owner
        .apply_temperature_refresh_with_admission(
            &mut history,
            &refresh,
            &next_shear,
            &cancelled,
            |_, _, _| {
                cancelled.request();
                Ok(())
            },
        );
    assert!(result.is_err());
    assert!(
        owner
            .apply_temperature_refresh(&mut history, &refresh, &next_shear, &cancelled)
            .is_err()
    );
    assert_eq!(history, initial);
    assert_eq!(owner.contact().normal.identity(), normal.identity());
    assert_eq!(owner.contact().shear, &shear);
    let mut reference = pair(&f, &normal, &shear, dt);
    assert_eq!(
        owner
            .advance(history, &CancelGate::new_clock_free())
            .unwrap(),
        reference
            .advance(initial, &CancelGate::new_clock_free())
            .unwrap(),
    );
    owner
        .apply_temperature_refresh(
            &mut history,
            &refresh,
            &next_shear,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let accepted = history;
    assert!(
        owner
            .apply_temperature_refresh(
                &mut history,
                &refresh,
                &next_shear,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert_eq!(history, accepted);
    assert_eq!(
        owner.contact().normal.identity(),
        refresh.source().identity()
    );
}

#[test]
fn an_advanced_normal_history_requires_the_pose_owner_to_consume_actual_travel() {
    let thermal = thermal_fixture::Fixture::with_options(1e5, 0.01, false, 320.0);
    let (normal, shear) = sources(&thermal);
    let dt = 5e-10;
    let mut f = fixture(-0.5 * common::Z, [Vec3::ZERO; 2]);
    let mut captured = pair(&f, &normal, &shear, dt);
    let state = captured.initialize([0.0; 2]).unwrap();
    let next = captured
        .advance(state, &CancelGate::new_clock_free())
        .unwrap();
    let refresh = normal
        .refresh_at_held_indentation(
            next.state_after.position_m()[0],
            thermal.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let next_shear = shear.rebind_dmt_temperature(&normal, &refresh).unwrap();
    assert!(matches!(
        captured.with_temperature_refresh(
            next.state_after,
            &refresh,
            &next_shear,
            &CancelGate::new_clock_free()
        ),
        Err(RigidDmtSourceRefreshError::Mechanical(
            RigidDmtError::GeometryMismatch { .. }
        ))
    ));
    for i in 0..2 {
        let after = next.bodies_after[i];
        f.bodies[i] = RigidBodyState::new(
            Pose::new(
                after
                    .pose()
                    .position_world()
                    .add(next.report.body_displacement_world_m[i]),
                after.pose().orientation(),
            )
            .unwrap(),
            after.linear_momentum_world(),
            after.angular_momentum_body(),
        )
        .unwrap();
    }
    let updated_geometry = pair(&f, &normal, &shear, dt);
    let (mut updated, history, report) = updated_geometry
        .with_temperature_refresh(
            next.state_after,
            &refresh,
            &next_shear,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(history.position_m(), next.state_after.position_m());
    assert!(report.geometry_residual_m.abs() <= controls().geometry_tolerance_m);
    let next = updated
        .advance(history, &CancelGate::new_clock_free())
        .unwrap();
    assert!(next.state_after.position_m()[0] > history.position_m()[0]);
    assert_eq!(next.report.dissipated_energy_j, 0.0);
}

#[test]
fn stronger_original_tail_refuses_at_both_corners_before_receiver_publication() {
    let thermal = thermal_fixture::Fixture::with_options(1e5, 0.01, false, 320.0);
    let (normal, shear) = sources(&thermal);
    let base = fixture(0.0, [Vec3::ZERO; 2]);
    let inverse_mass = 1.0 / base.properties[0].mass() + 1.0 / base.properties[1].mass();
    let tail =
        2.0 * PI * normal.effective_radius_m() * normal.normal_state().adhesion_energy_j_per_m2()
            / normal.interaction_range_m();
    // The source's original 1 - S*k/2 margin is 0.2. Cooling increases W by
    // 45%, so the same complete native mechanical response must now refuse.
    let dt = (3.2 / (inverse_mass * tail)).sqrt();
    for delta in [0.0, -common::Z] {
        let f = fixture(delta, [Vec3::ZERO; 2]);
        let mut owner = pair(&f, &normal, &shear, dt);
        let mut history = owner.initialize([0.0; 2]).unwrap();
        let initial = history;
        let no_change = normal
            .refresh_at_held_indentation(delta, thermal.cards(), temperatures(300.0))
            .unwrap();
        let same_shear = shear.rebind_dmt_temperature(&normal, &no_change).unwrap();
        owner
            .apply_temperature_refresh(
                &mut history,
                &no_change,
                &same_shear,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(history, initial);
        let cold = normal
            .refresh_at_held_indentation(delta, thermal.cards(), temperatures(285.0))
            .unwrap();
        let cold_shear = shear.rebind_dmt_temperature(&normal, &cold).unwrap();
        let mut called = false;
        let refused: Result<(), RigidDmtSourceRefreshError> = owner
            .apply_temperature_refresh_with_admission(
                &mut history,
                &cold,
                &cold_shear,
                &CancelGate::new_clock_free(),
                |_, _, _| {
                    called = true;
                    Ok(())
                },
            );
        assert!(matches!(
            refused,
            Err(RigidDmtSourceRefreshError::Mechanical(
                RigidDmtError::Contact(DmtCoupledError::MechanicalStability)
            ))
        ));
        assert!(!called);
        assert_eq!(history, initial);
        assert_eq!(owner.contact().normal.identity(), normal.identity());
        assert_eq!(owner.contact().shear, &shear);
    }
}

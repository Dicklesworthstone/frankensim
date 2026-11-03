//! Actual stopping-path heat, closed-path admission and coupled rollback.
use super::*;
use fs_couple::friction_heat::rigid_cohesive::interval::*;

#[test]
fn stopping_then_static_hold_and_breakaway_keep_all_native_path_heat() {
    let (source, friction) = sources([0.6, 0.3], false);
    let incoming = 3.7e-4;
    let mut owner = crate::interval::translation_fixture(&source, &friction, incoming);
    let mut reference = crate::interval::translation_fixture(&source, &friction, incoming);
    let jump = v(owner.sites()[0].point().numerical().material().jump_m());
    let loads = balanced_load(Vec3::new(-0.4, 0.0, 0.0), jump);
    let duration = 0.02;
    let plain = reference
        .advance_interval(loads, duration, crate::interval::controls())
        .unwrap();
    let initial = [10.0, 0.0];
    let mut heat = thermal(initial);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let report = apply_rigid_cohesive_interval_heat(
        &mut owner,
        loads,
        duration,
        crate::interval::controls(),
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(report.mechanical, plain);
    assert_eq!(owner.bodies(), reference.bodies());
    assert_eq!(owner.sites(), reference.sites());
    assert_eq!(
        owner.damaged_area_friction_states(),
        reference.damaged_area_friction_states()
    );
    assert!(report.mechanical.attempted_substeps > report.mechanical.accepted_substeps);
    assert!(report.mechanical.accepted_substeps >= 2);
    let reduced_mass = 1.0 / (1.0 / 0.4 + 1.0 / 0.7);
    let dissipated = 0.5 * reduced_mass * incoming * incoming;
    close(report.heat.dissipated_j, dissipated, 2e-12);
    close(
        report.heat.dissipated_j,
        report.mechanical.sites[0].friction_dissipation_j,
        1e-20,
    );
    assert_eq!(owner.last_frame().unwrap().friction_dissipation_j, 0.0);
    assert_eq!(report.fracture_work_j, 0.0);
    assert_eq!(report.mechanical.viscous_dissipation_j, 0.0);
    let (enthalpy, face) = thermal_reference(initial, duration, dissipated);
    for (i, expected) in enthalpy.into_iter().enumerate() {
        close(
            heat.cells()[&CellId(i as u64 + 1)].specific_h_j_kg,
            expected,
            2e-9,
        );
    }
    close(
        report.heat.interfaces[0].interface_temperature_k,
        face,
        2e-8,
    );
    assert!(report.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(report.heat.interfaces[0].body_heat_j[1] > report.heat.dissipated_j);
    close(
        report.heat.entropy.unwrap().friction_generation_j_k,
        report.heat.dissipated_j / face,
        1e-15,
    );
    close(report.combined_energy_residual_j, 0.0, 3e-10);
    close(report.summed_energy_residual_j, 0.0, 3e-10);
    assert_eq!(owner.time_s(), heat.time_s());
    assert_eq!(
        owner.damaged_area_friction_states().unwrap()[0].occupancy(),
        DamageFrictionRegime::Sticking
    );

    // A later static interval conducts stored heat but invents no friction heat.
    let jump = v(owner.sites()[0].point().numerical().material().jump_m());
    let held = apply_rigid_cohesive_interval_heat(
        &mut owner,
        balanced_load(Vec3::new(-0.4, 0.2, 0.0), jump),
        0.007,
        crate::interval::controls(),
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(held.heat.dissipated_j, 0.0);
    close(held.heat.thermal.stored_energy_change_j, 0.0, 1e-10);
    assert_eq!(owner.time_s(), heat.time_s());
    close(
        owner.next_step_time_s().unwrap(),
        owner.time_s() + owner.step_duration_s(),
        0.0,
    );

    // Crossing the original static disk resumes actual sliding and heat.
    let jump = v(owner.sites()[0].point().numerical().material().jump_m());
    let broken = apply_rigid_cohesive_interval_heat(
        &mut owner,
        balanced_load(Vec3::new(-0.4, 0.3, 0.0), jump),
        0.002,
        crate::interval::controls(),
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(broken.heat.dissipated_j > 0.0);
    close(
        broken.heat.dissipated_j,
        broken.mechanical.friction_dissipation_j,
        1e-20,
    );
    assert_eq!(
        owner.damaged_area_friction_states().unwrap()[0].occupancy(),
        DamageFrictionRegime::Sliding
    );
    assert_eq!(owner.time_s(), heat.time_s());
}

#[test]
fn thermal_source_failure_after_native_stopping_restores_both_original_owners() {
    let (source, friction) = sources([0.6, 0.3], false);
    let mut owner = crate::interval::translation_fixture(&source, &friction, 3.7e-4);
    let jump = v(owner.sites()[0].point().numerical().material().jump_m());
    let loads = balanced_load(Vec3::new(-0.4, 0.0, 0.0), jump);
    let mut heat = thermal([0.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let old_bodies = owner.bodies().clone();
    let old_sites = owner.sites().to_vec();
    let old_friction = owner.damaged_area_friction_states().unwrap().to_vec();
    let old_cells = heat.cells().clone();
    let old_next = owner.next_step_time_s().unwrap();
    for cancel in [false, true] {
        let gate = CancelGate::new_clock_free();
        let calls = Cell::new(0);
        let refused = apply_rigid_cohesive_interval_heat(
            &mut owner,
            loads,
            0.02,
            crate::interval::controls(),
            &mut heat,
            bodies(&source),
            &contact(&route),
            controls(),
            &[0.0; 2],
            |_, ids, _| {
                calls.set(calls.get() + 1);
                if cancel {
                    gate.request();
                    Ok(inputs(ids))
                } else {
                    Err("late interval conductivity source")
                }
            },
            &gate,
        );
        if cancel {
            assert!(refused.is_err());
        } else {
            assert!(matches!(
                refused,
                Err(RigidCohesiveHeatError::Heat(CohesiveHeatError::Source(
                    "late interval conductivity source"
                )))
            ));
        }
        assert!(calls.get() > 0); // Native stopping and resticking already succeeded.
        assert_eq!(owner.bodies(), &old_bodies);
        assert_eq!(owner.sites(), old_sites);
        assert_eq!(owner.damaged_area_friction_states().unwrap(), old_friction);
        assert!(owner.last_frame().is_none());
        assert!(owner.site_frames().iter().all(Option::is_none));
        assert_eq!(owner.time_s(), 0.0);
        assert_eq!(owner.next_step_time_s().unwrap(), old_next);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(heat.cells(), &old_cells);
    }
    let accepted = apply_rigid_cohesive_interval_heat(
        &mut owner,
        loads,
        0.02,
        crate::interval::controls(),
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(accepted.heat.dissipated_j > 0.0);
    assert_eq!(owner.time_s(), heat.time_s());
    assert_eq!(
        owner.damaged_area_friction_states().unwrap()[0].occupancy(),
        DamageFrictionRegime::Sticking
    );
}

fn opening_fixture<'a>(
    source: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
) -> RigidCohesiveSystem<'a> {
    system(
        source,
        Some(friction),
        failed(source),
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, Vec3::new(0.02, 0.0, 0.0)],
        zero(),
        0.002,
        config(),
    )
}

#[test]
fn an_intermediate_open_contact_cannot_hide_behind_a_closed_final_endpoint() {
    let (source, friction) = sources([0.6, 0.3], false);
    let mut owner = opening_fixture(&source, &friction);
    let mut reference = opening_fixture(&source, &friction);
    let jump = v(owner.sites()[0].point().numerical().material().jump_m());
    let loads = balanced_load(Vec3::new(-0.4, 0.0, 0.0), jump);
    let native = reference
        .advance_interval(loads, 0.04, crate::interval::controls())
        .unwrap();
    assert!(native.sites[0].maximum_normal_jump_m > 0.0);
    assert!(reference.sites()[0].point().numerical().material().jump_m()[0] <= 0.0);
    let mut heat = thermal([0.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let old_bodies = owner.bodies().clone();
    let old_sites = owner.sites().to_vec();
    let old_cells = heat.cells().clone();
    let calls = Cell::new(0);
    let refused = apply_rigid_cohesive_interval_heat(
        &mut owner,
        loads,
        0.04,
        crate::interval::controls(),
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(inputs(ids))
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(refused, Err(RigidCohesiveHeatError::Invalid(_))));
    assert_eq!(calls.get(), 0);
    assert_eq!(owner.bodies(), &old_bodies);
    assert_eq!(owner.sites(), old_sites);
    assert!(owner.last_frame().is_none());
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(heat.cells(), &old_cells);
}

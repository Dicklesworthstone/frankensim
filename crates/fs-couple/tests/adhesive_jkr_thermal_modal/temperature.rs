//! Prescribed source temperature changes retain original native/thermal history.
use super::*;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod temperature_fixture;

fn thermal_source(fixture: &temperature_fixture::Fixture) -> SourceJkrSphereContact {
    SourceJkrSphereContact::new(fixture.normal.clone(), 0.01, 1e-9, source().envelope()).unwrap()
}
fn temperatures(t: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0 + 0.5 * (t - 300.0),
        surface_b_k: 300.0 - 0.25 * (t - 300.0),
        interface_k: t,
    }
}

#[test]
fn temperature_changes_subsequent_native_force_without_impulse_heat_or_history_reset() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let (d, e, v) = scales(&law);
    let k = 0.02 * e / (d * d);
    let c = 0.4 * (k * MU).sqrt();
    let state = law.evaluate_attached(0.0).unwrap().state();
    let make = || mechanical_with_damping(&law, state, 0.2 * v, Some((k, 0.0, c)));
    let caloric = curves(&law, [1.0; 2]);
    let origin = 7.0 * DT;
    let mut frozen = make();
    let mut owner = joint_with_support_heat(
        make(),
        timeline(&law, &caloric, [1.0; 2], origin),
        &caloric,
        origin,
    );
    for _ in 0..2 {
        owner.step(&[0.0; 2]).unwrap();
        frozen.step(&[0.0; 2]).unwrap();
    }
    let before = owner.mechanical().state();
    let native = states(owner.mechanical());
    let cells = owner.thermal().cells().clone();
    let time = owner.thermal().time_s();
    let heat = *owner.heat_history().ledger();
    let stamp = owner.heat_history().stamp();
    let energy = owner.stored_energy_j().unwrap();
    let last = owner.last_frame().cloned();
    let refresh = law
        .refresh_at_held_indentation(before.contact_state(), fixture.cards(), temperatures(315.0))
        .unwrap();
    let report = owner
        .continue_at_temperature(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_ne!(
        report.source.before.normal_force_n,
        report.source.after.normal_force_n
    );
    assert_ne!(
        report.source.before.contact_area_m2,
        report.source.after.contact_area_m2
    );
    assert_ne!(
        report.source.original_model_identity,
        report.source.refreshed_model_identity
    );
    assert_ne!(report.source.parameter_work_j, 0.0);
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert_eq!(states(owner.mechanical()), native);
    assert_eq!(owner.mechanical().state().time_s(), before.time_s());
    assert_eq!(owner.mechanical().state().steps(), before.steps());
    assert_eq!(
        owner.mechanical().state().velocity_m_s(),
        before.velocity_m_s()
    );
    assert_eq!(owner.mechanical().samples_rendered(), 2);
    assert_eq!(owner.thermal().time_s(), time);
    assert_eq!(owner.thermal().cells(), &cells);
    assert_eq!(*owner.heat_history().ledger(), heat);
    assert_ne!(owner.heat_history().stamp(), stamp);
    assert_eq!(owner.last_frame(), last.as_ref());
    close(
        owner.stored_energy_j().unwrap() - energy,
        report.source.parameter_work_j,
        e,
        1e-8,
    );
    let next = owner.step(&[0.0; 2]).unwrap().clone();
    let old = frozen.step(&[0.0; 2]).unwrap().clone();
    assert_eq!(next.mechanical.path.initial_state().steps(), before.steps());
    assert_eq!(
        next.mechanical.path.initial_state().time_s(),
        before.time_s()
    );
    assert_ne!(
        next.mechanical.path.state().velocity_m_s(),
        old.path.state().velocity_m_s()
    );
    assert_eq!(
        next.mechanical.path.state().contact_state(),
        owner.heat_history().state()
    );
    assert!(next.combined_energy_residual_j.abs() <= next.combined_energy_tolerance_j);
}

#[test]
fn unchanged_temperature_preserves_snap_ledger_and_later_path_bit_for_bit() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let (d, _, v) = scales(&law);
    let state = law.detached_state(-0.25 * d).unwrap();
    let caloric = curves(&law, [1.0; 2]);
    let make = || {
        joint(
            mechanical(&law, state, 10.0 * v, None),
            timeline(&law, &caloric, [1.0; 2], 0.0),
            &caloric,
            0.0,
        )
    };
    let mut owner = make();
    let mut reference = make();
    for _ in 0..2 {
        owner.step(&[0.0; 2]).unwrap();
        reference.step(&[0.0; 2]).unwrap();
    }
    assert!(owner.heat_history().ledger().dissipated_work_j() > 0.0);
    let motion = owner.mechanical().state();
    let ledger = *owner.heat_history().ledger();
    let stamp = owner.heat_history().stamp();
    let refresh = law
        .refresh_at_held_indentation(motion.contact_state(), fixture.cards(), temperatures(300.0))
        .unwrap();
    let report = owner
        .continue_at_temperature(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(report.source.parameter_work_j, 0.0);
    assert_eq!(owner.mechanical().state(), motion);
    assert_eq!(*owner.heat_history().ledger(), ledger);
    assert_ne!(owner.heat_history().stamp(), stamp);
    assert_eq!(
        owner.step(&[0.0; 2]).unwrap(),
        reference.step(&[0.0; 2]).unwrap()
    );
    assert_eq!(states(owner.mechanical()), states(reference.mechanical()));
    assert_eq!(owner.thermal().cells(), reference.thermal().cells());
}

#[test]
fn refused_temperature_candidate_and_late_cancellation_publish_no_native_change() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let state = law.evaluate_attached(0.0).unwrap().state();
    let mut owner = mechanical(&law, state, 0.0, None);
    owner.step(&[0.0; 2]).unwrap();
    let native = states(&owner);
    let motion = owner.state();
    let frame = owner.last_frame().cloned();
    let refresh = law
        .refresh_at_held_indentation(motion.contact_state(), fixture.cards(), temperatures(315.0))
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let refused = owner.continue_at_temperature_with_admission(&refresh, &gate, |_, _| {
        Err::<(), _>(JkrNormalEventModalError::Compatibility {
            what: "receiving history refused",
            residual: 1.0,
            tolerance: 0.0,
        })
    });
    assert!(refused.is_err());
    let cancelled = owner.continue_at_temperature_with_admission(&refresh, &gate, |_, _| {
        gate.request();
        Ok::<_, JkrNormalEventModalError>(())
    });
    assert!(cancelled.is_err());
    assert_eq!(states(&owner), native);
    assert_eq!(owner.state(), motion);
    assert_eq!(owner.last_frame(), frame.as_ref());
    assert_eq!(owner.normal().identity(), law.identity());
    owner
        .continue_at_temperature(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(owner.state().steps(), motion.steps());
    assert_eq!(owner.normal().identity(), refresh.source().identity());
}

#[test]
fn foreign_refresh_and_cancelled_joint_refresh_preserve_heat_and_source() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let state = law.evaluate_attached(0.0).unwrap().state();
    let caloric = curves(&law, [1.0; 2]);
    let mut owner = joint(
        mechanical(&law, state, 0.0, None),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    let foreign =
        SourceJkrSphereContact::new(fixture.normal.clone(), 0.011, 1e-9, law.envelope()).unwrap();
    let wrong = foreign
        .refresh_at_held_indentation(
            foreign.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let cells = owner.thermal().cells().clone();
    let stamp = owner.heat_history().stamp();
    let motion = owner.mechanical().state();
    assert!(
        owner
            .continue_at_temperature(&wrong, &CancelGate::new_clock_free())
            .is_err()
    );
    let refresh = law
        .refresh_at_held_indentation(state, fixture.cards(), temperatures(315.0))
        .unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(owner.continue_at_temperature(&refresh, &gate).is_err());
    assert_eq!(owner.mechanical().state(), motion);
    assert_eq!(owner.heat_history().stamp(), stamp);
    assert_eq!(owner.thermal().cells(), &cells);
    assert_eq!(owner.mechanical().normal().identity(), law.identity());
    assert_eq!(owner.thermal().time_s(), 0.0);
}

//! Original-source changes preserve the complete native mechanical history.
use super::*;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;

#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod temperature_fixture;

fn temperatures(t: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0 + 0.5 * (t - 300.0),
        surface_b_k: 300.0 - 0.25 * (t - 300.0),
        interface_k: t,
    }
}
fn thermal_source(fixture: &temperature_fixture::Fixture) -> SourceJkrSphereContact {
    SourceJkrSphereContact::new(fixture.normal.clone(), 0.01, 1e-9, source().envelope()).unwrap()
}
fn thermal_owner(law: &SourceJkrSphereContact) -> SourceNormalModalSystem {
    let config = controls(law);
    SourceNormalModalSystem::new(
        network(law, true, 0.0),
        BoundNormalContact::new(
            NormalLawSource::JkrAttached(Box::new(law.clone())),
            column(),
            0.0,
        )
        .unwrap(),
        config,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn native_bits(owner: &SourceNormalModalSystem) -> Vec<u64> {
    states(owner.network())
        .iter()
        .flat_map(|state| {
            [
                state.displacement_m_sqrt_kg.to_bits(),
                state.velocity_m_sqrt_kg_per_s.to_bits(),
            ]
        })
        .chain(
            owner
                .network()
                .material_memory()
                .iter()
                .map(|x| x.to_bits()),
        )
        .collect()
}

#[test]
fn source_temperature_changes_force_and_signed_work_without_native_history_or_heat_changes() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let mut work = Vec::new();
    for t in [300.0, 285.0, 315.0] {
        let mut owner = thermal_owner(&law);
        let mut reference = thermal_owner(&law);
        for _ in 0..2 {
            owner.step(&[0.0; 4]).unwrap();
            reference.step(&[0.0; 4]).unwrap();
        }
        let approach = owner.last_frame().unwrap().approach_after_m;
        let refresh = law
            .refresh_at_held_indentation(
                law.evaluate_attached(approach).unwrap().state(),
                fixture.cards(),
                temperatures(t),
            )
            .unwrap();
        let body_and_memory = native_bits(&owner);
        let network_energy = owner.network().total_energy_j().unwrap().to_bits();
        let total = owner.total_energy_j().unwrap();
        let last = owner.last_frame().copied();
        let native = owner.network().last_frame().cloned();
        let period = owner.sample_period_s().to_bits();
        let config = owner.config();
        let row = owner.contact().work_column().to_vec();
        let gap = owner.contact().gap_m().to_bits();
        if t == 315.0 {
            // Heating moves the original fold inside the old authored bracket.
            // The actual retained contact is still stable, so continuation is valid.
            assert!(config.minimum_approach_m < refresh.source().minimum_attached_indentation_m());
        }
        let report = owner
            .apply_jkr_source_refresh(&refresh, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(native_bits(&owner), body_and_memory);
        assert_eq!(
            owner.network().total_energy_j().unwrap().to_bits(),
            network_energy
        );
        assert_eq!(owner.samples_rendered(), 2);
        assert_eq!(report.sample, 2);
        assert_eq!(report.approach_m.to_bits(), approach.to_bits());
        assert_eq!(owner.sample_period_s().to_bits(), period);
        assert_eq!(owner.last_frame(), last.as_ref());
        assert_eq!(owner.network().last_frame(), native.as_ref());
        assert_eq!(owner.config(), config);
        assert_eq!(owner.contact().work_column(), row);
        assert_eq!(owner.contact().gap_m().to_bits(), gap);
        assert_eq!(
            owner.contact().source().identity(),
            refresh.source().identity()
        );
        assert_eq!(report.parameter_work_j, refresh.parameter_work_j());
        assert_eq!(refresh.mechanical_work_j(), 0.0);
        close(
            owner.total_energy_j().unwrap() - total,
            report.parameter_work_j,
            1e-9,
            1e-9,
        );
        close(
            report.contact_energy_after_j - report.contact_energy_before_j,
            report.parameter_work_j,
            1e-9,
            1e-10,
        );
        assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
        // The unchanged native frame retains its original physical losses;
        // parameter work was not added to either connection or material heat.
        assert!(native.as_ref().unwrap().material_dissipation_j > 0.0);
        // A subsequent outward native load needs the lower trial bracket;
        // it must use the successor's source domain without changing config.
        let outward: Vec<_> = column().into_iter().map(|b| -0.02 * b).collect();
        let next = *owner.step(&outward).unwrap();
        let old = *reference.step(&outward).unwrap();
        assert!(next.approach_after_m < next.approach_before_m);
        if t == 300.0 {
            assert_eq!(report.parameter_work_j, 0.0);
            assert_eq!(next, old);
            assert_eq!(native_bits(&owner), native_bits(&reference));
        } else {
            assert_ne!(next.normal_force_n.to_bits(), old.normal_force_n.to_bits());
            assert_ne!(native_bits(&owner), native_bits(&reference));
        }
        work.push(report.parameter_work_j);
    }
    assert!(work[1] < 0.0);
    assert!(work[2] > 0.0);
}

#[test]
fn foreign_and_late_refused_source_refresh_never_publish_partial_native_change() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let mut owner = thermal_owner(&law);
    owner.step(&[0.0; 4]).unwrap();
    let approach = owner.last_frame().unwrap().approach_after_m;
    let refresh = law
        .refresh_at_held_indentation(
            law.evaluate_attached(approach).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let body_and_memory = native_bits(&owner);
    let energy = owner.total_energy_j().unwrap().to_bits();
    let frame = owner.last_frame().copied();
    let native = owner.network().last_frame().cloned();
    let source_identity = owner.contact().source().identity();
    let foreign =
        SourceJkrSphereContact::new(fixture.normal.clone(), 0.011, 1e-9, law.envelope()).unwrap();
    let foreign_refresh = foreign
        .refresh_at_held_indentation(
            foreign.evaluate_attached(approach).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    assert!(
        owner
            .apply_jkr_source_refresh(&foreign_refresh, &CancelGate::new_clock_free())
            .is_err()
    );
    let stale = law
        .refresh_at_held_indentation(
            law.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    assert!(
        owner
            .apply_jkr_source_refresh(&stale, &CancelGate::new_clock_free())
            .is_err()
    );
    let refused: Result<(), SourceNormalModalError> = owner
        .apply_jkr_source_refresh_with_admission(
            &refresh,
            &CancelGate::new_clock_free(),
            |_, _| {
                Err(SourceNormalModalError::Unresolved(
                    "parameter-work receiver refused",
                ))
            },
        );
    assert!(refused.is_err());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), SourceNormalModalError> = owner
        .apply_jkr_source_refresh_with_admission(&refresh, &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(cancelled.is_err());
    assert_eq!(native_bits(&owner), body_and_memory);
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.samples_rendered(), 1);
    assert_eq!(owner.last_frame(), frame.as_ref());
    assert_eq!(owner.network().last_frame(), native.as_ref());
    assert_eq!(owner.contact().source().identity(), source_identity);
    owner
        .apply_jkr_source_refresh(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(native_bits(&owner), body_and_memory);
    assert_eq!(owner.samples_rendered(), 1);
    let accepted_energy = owner.total_energy_j().unwrap().to_bits();
    // The old-source token cannot debit the same signed parameter work twice.
    assert!(
        owner
            .apply_jkr_source_refresh(&refresh, &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), accepted_energy);
}

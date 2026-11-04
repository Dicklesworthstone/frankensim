//! Thermal finite-zone sources continue the original native modal/Maxwell owner.
use super::*;
use fs_contact::normal_patch::maugis::continuation::{
    MaugisContinuationBudget, MaugisPhase,
};
use fs_contact::normal_patch::maugis::mechanical::MaugisAttachedPotential;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, MaugisWorkBudget, SourceMaugisSphereContact,
};

fn maugis_source(
    fixture: &temperature_fixture::Fixture,
    lambda: f64,
) -> SourceMaugisSphereContact {
    let r: f64 = 0.01;
    let w = fixture.normal.adhesion_energy_j_per_m2();
    let e = fixture.normal.elastic().reduced_modulus_pa();
    let a0 = (core::f64::consts::PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    SourceMaugisSphereContact::new(
        fixture.normal.clone(),
        r,
        2.0 * d0 / (core::f64::consts::PI * lambda),
        MaugisEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}
fn maugis_potential(law: &SourceMaugisSphereContact) -> MaugisAttachedPotential {
    let r = law.effective_radius_m();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let a0 = (core::f64::consts::PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    let f0 = core::f64::consts::PI * w * r;
    MaugisAttachedPotential::new(
        law.clone(),
        MaugisContinuationBudget {
            work: MaugisWorkBudget {
                absolute_tolerance_j: 1e-12 * f0 * d0,
                relative_tolerance: 1e-11,
                maximum_evaluations: 10000,
            },
            energy_absolute_tolerance_j: 1e-10 * f0 * d0,
            energy_relative_tolerance: 1e-10,
        },
    )
    .unwrap()
}
fn initial_approach(law: &SourceMaugisSphereContact, connected: bool) -> f64 {
    if connected {
        let limits = law.noncontact_limits().unwrap().unwrap();
        0.25 * limits.minimum_indentation_m + 0.75 * limits.first_contact_indentation_m
    } else {
        0.0
    }
}
fn maugis_owner(
    fixture: &temperature_fixture::Fixture,
    law: &SourceMaugisSphereContact,
    connected: bool,
) -> SourceNormalModalSystem {
    let reference_normal = thermal_source(fixture);
    let mut potential = maugis_potential(law);
    if connected {
        potential = potential.with_noncontact_continuation().unwrap();
    }
    let d0 = potential.indentation_scale_m();
    let minimum = potential.noncontact_continuation().map_or_else(
        || law.minimum_attached_indentation_m(),
        |limits| limits.minimum_indentation_m,
    );
    let mut config = controls(&reference_normal);
    config.minimum_approach_m = minimum + 1e-3 * d0;
    config.maximum_approach_m = 4.0 * d0;
    SourceNormalModalSystem::new(
        network(&reference_normal, true, 0.0),
        BoundNormalContact::new(
            potential.into(),
            column(),
            -initial_approach(law, connected),
        )
        .unwrap(),
        config,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn original_potential(owner: &SourceNormalModalSystem) -> &MaugisAttachedPotential {
    let NormalLawSource::MaugisAttached(source) = owner.contact().source() else {
        panic!("original finite-zone source missing");
    };
    source
}

#[test]
fn thermal_contact_and_noncontact_sources_preserve_native_history_and_continue_after_shifted_fold() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    for (lambda, connected) in [(5.0, false), (0.2, true)] {
        let law = maugis_source(&fixture, lambda);
        let mut work = Vec::new();
        for t in [300.0, 285.0, 315.0] {
            let mut owner = maugis_owner(&fixture, &law, connected);
            let mut reference = maugis_owner(&fixture, &law, connected);
            for _ in 0..2 {
                owner.step(&[0.0; 4]).unwrap();
                reference.step(&[0.0; 4]).unwrap();
            }
            let approach = owner.last_frame().unwrap().approach_after_m;
            let before = original_potential(&owner).evaluate(approach).unwrap();
            assert_eq!(
                before.phase(),
                if connected { MaugisPhase::Noncontact } else { MaugisPhase::Contact }
            );
            let refresh = law
                .refresh_at_held_indentation(before, fixture.cards(), temperatures(t))
                .unwrap();
            let original_budget = original_potential(&owner).budget();
            let original_policy = original_potential(&owner).noncontact_continuation();
            let body_and_memory = native_bits(&owner);
            let network_energy = owner.network().total_energy_j().unwrap().to_bits();
            let total = owner.total_energy_j().unwrap();
            let last = owner.last_frame().copied();
            let native = owner.network().last_frame().cloned();
            let period = owner.sample_period_s().to_bits();
            let config = owner.config();
            let row = owner.contact().work_column().to_vec();
            let gap = owner.contact().gap_m().to_bits();
            let report = owner
                .apply_maugis_source_refresh(&refresh, &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(native_bits(&owner), body_and_memory);
            assert_eq!(owner.network().total_energy_j().unwrap().to_bits(), network_energy);
            assert_eq!(owner.samples_rendered(), 2);
            assert_eq!(report.sample, 2);
            assert_eq!(report.approach_m.to_bits(), approach.to_bits());
            assert_eq!(owner.sample_period_s().to_bits(), period);
            assert_eq!(owner.last_frame(), last.as_ref());
            assert_eq!(owner.network().last_frame(), native.as_ref());
            assert_eq!(owner.config(), config);
            assert_eq!(owner.contact().work_column(), row);
            assert_eq!(owner.contact().gap_m().to_bits(), gap);
            assert_eq!(owner.contact().source().identity(), refresh.source().identity());
            let accepted = original_potential(&owner);
            assert_eq!(accepted.budget(), original_budget);
            assert_eq!(accepted.noncontact_continuation().is_some(), connected);
            assert_eq!(accepted.evaluate(approach).unwrap().phase(), before.phase());
            if connected {
                assert_eq!(
                    accepted.evaluate(approach).unwrap().friction_contact_state().unwrap().area_m2(),
                    0.0
                );
                if t != 300.0 {
                    assert_ne!(accepted.noncontact_continuation(), original_policy);
                }
            }
            let lower = accepted.noncontact_continuation().map_or_else(
                || accepted.source().minimum_attached_indentation_m(),
                |limits| limits.minimum_indentation_m,
            );
            if t == 315.0 {
                assert!(lower > config.minimum_approach_m);
            }
            assert_eq!(report.parameter_work_j, refresh.parameter_work_j());
            assert_eq!(refresh.mechanical_work_j(), 0.0);
            close(owner.total_energy_j().unwrap() - total, report.parameter_work_j, 1e-9, 1e-9);
            close(
                report.contact_energy_after_j - report.contact_energy_before_j,
                report.parameter_work_j,
                1e-9,
                1e-10,
            );
            assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
            assert!(native.as_ref().unwrap().material_dissipation_j > 0.0);

            // An outward solve needs the refreshed original lower trial domain.
            // It retains the complete original material memory and work ports.
            let outward: Vec<_> = column().into_iter().map(|b| -0.02 * b).collect();
            let next = *owner.step(&outward).unwrap();
            let old = *reference.step(&outward).unwrap();
            assert!(next.approach_after_m < next.approach_before_m);
            assert!(next.contact_work_residual_j.abs() <= next.contact_work_tolerance_j);
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
        assert!(work[1] < 0.0 && work[2] > 0.0);
    }
}

#[test]
fn maugis_refresh_refusal_and_late_cancellation_preserve_the_original_complete_network() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = maugis_source(&fixture, 5.0);
    let mut owner = maugis_owner(&fixture, &law, false);
    owner.step(&[0.0; 4]).unwrap();
    let approach = owner.last_frame().unwrap().approach_after_m;
    let refresh = law
        .refresh_at_held_indentation(
            original_potential(&owner).evaluate(approach).unwrap(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let body_and_memory = native_bits(&owner);
    let energy = owner.total_energy_j().unwrap().to_bits();
    let frame = owner.last_frame().copied();
    let native = owner.network().last_frame().cloned();
    let source_identity = owner.contact().source().identity();
    let foreign = maugis_source(&fixture, 4.5);
    let foreign_refresh = foreign
        .refresh_at_held_indentation(
            foreign.evaluate_contact_energy(approach).unwrap(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    assert!(owner.apply_maugis_source_refresh(&foreign_refresh, &CancelGate::new_clock_free()).is_err());
    let stale = law
        .refresh_at_held_indentation(
            law.evaluate_contact_energy(0.0).unwrap(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    assert!(owner.apply_maugis_source_refresh(&stale, &CancelGate::new_clock_free()).is_err());
    let load = law.refresh_at_held_load(0.0, fixture.cards(), temperatures(315.0)).unwrap();
    assert!(owner.apply_maugis_source_refresh(&load, &CancelGate::new_clock_free()).is_err());
    let refused: Result<(), SourceNormalModalError> = owner
        .apply_maugis_source_refresh_with_admission(
            &refresh,
            &CancelGate::new_clock_free(),
            |_, _| Err(SourceNormalModalError::Unresolved("parameter-work receiver refused")),
        );
    assert!(refused.is_err());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), SourceNormalModalError> = owner
        .apply_maugis_source_refresh_with_admission(&refresh, &gate, |_, _| {
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
    owner.apply_maugis_source_refresh(&refresh, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(native_bits(&owner), body_and_memory);
    assert_eq!(owner.samples_rendered(), 1);
    let accepted_energy = owner.total_energy_j().unwrap().to_bits();
    assert!(owner.apply_maugis_source_refresh(&refresh, &CancelGate::new_clock_free()).is_err());
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), accepted_energy);
}

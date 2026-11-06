//! Finite-range thermal adhesion on the original scalar and shared modal owners.
use super::*;
use fs_contact::normal_patch::dmt::{DmtEnvelope, DmtRegime, SourceDmtSphereContact};
use fs_couple::render::schedule::force::coupled::contact::{
    multiple::MultiContactConfig,
    source_normal::multiple::{MultipleSourceNormalModalSystem, temperature::NormalSourceRefresh},
};

fn dmt(fixture: &temperature_fixture::Fixture) -> SourceDmtSphereContact {
    SourceDmtSphereContact::new(
        fixture.normal.clone(),
        0.01,
        2e-5,
        DmtEnvelope {
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
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}
fn dmt_controls(
    jkr: &SourceJkrSphereContact,
    law: &SourceDmtSphereContact,
) -> SourceNormalModalConfig {
    let mut c = controls(jkr);
    c.minimum_approach_m = -2.0 * law.interaction_range_m();
    c.maximum_approach_m = law.envelope().maximum_indentation_m;
    c
}
fn dmt_owner(
    fixture: &temperature_fixture::Fixture,
    law: &SourceDmtSphereContact,
    approach: f64,
    scale: f64,
) -> SourceNormalModalSystem {
    let jkr = thermal_source(fixture);
    SourceNormalModalSystem::new(
        network(&jkr, true, 0.0),
        BoundNormalContact::new(
            law.clone().into(),
            column().into_iter().map(|x| x * scale).collect(),
            -approach,
        )
        .unwrap(),
        dmt_controls(&jkr, law),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn batch_controls() -> MultiContactConfig {
    MultiContactConfig {
        max_contacts: 4,
        max_sweeps: 96,
        max_setup_terms: 32768,
    }
}
fn batch_bits(owner: &MultipleSourceNormalModalSystem) -> Vec<u64> {
    states(owner.network())
        .iter()
        .flat_map(|s| {
            [
                s.displacement_m_sqrt_kg.to_bits(),
                s.velocity_m_sqrt_kg_per_s.to_bits(),
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
fn contact_gap_cutoff_and_separated_refresh_retain_native_state_and_original_work() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = dmt(&fixture);
    let range = law.interaction_range_m();
    for initial in [1e-6, 0.0, -0.5 * range, -range, -1.5 * range] {
        for t in [300.0, 285.0, 315.0] {
            let mut owner = dmt_owner(&fixture, &law, initial, 1.0);
            let mut reference = dmt_owner(&fixture, &law, initial, 1.0);
            // Preserve exact corner fixtures; interior fixtures retain an
            // already accepted native frame and evolved material memory.
            if initial != 0.0 && initial != -range {
                for _ in 0..2 {
                    owner.step(&[0.0; 4]).unwrap();
                    reference.step(&[0.0; 4]).unwrap();
                }
            }
            let approach = owner.last_frame().map_or(initial, |f| f.approach_after_m);
            let refresh = law
                .refresh_at_held_indentation(approach, fixture.cards(), temperatures(t))
                .unwrap();
            let bits = native_bits(&owner);
            let native = owner.network().last_frame().cloned();
            let frame = owner.last_frame().copied();
            let sample = owner.samples_rendered();
            let energy = owner.total_energy_j().unwrap();
            let report = owner
                .apply_dmt_source_refresh(&refresh, &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(native_bits(&owner), bits);
            assert_eq!(owner.network().last_frame(), native.as_ref());
            assert_eq!(owner.last_frame(), frame.as_ref());
            assert_eq!(owner.samples_rendered(), sample);
            assert_eq!(report.sample, sample);
            assert_eq!(report.approach_m, approach);
            assert_eq!(
                owner.contact().source().identity(),
                refresh.source().identity()
            );
            assert_eq!(report.parameter_work_j, refresh.parameter_work_j());
            assert_eq!(refresh.before().regime, refresh.response().regime);
            assert_eq!(refresh.mechanical_work_j(), 0.0);
            assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
            close(
                owner.total_energy_j().unwrap() - energy,
                report.parameter_work_j,
                1e-18,
                1e-9,
            );
            if refresh.response().regime != DmtRegime::Contact {
                assert_eq!(refresh.response().contact_area_m2, 0.0);
                assert_eq!(refresh.response().repulsive_normal_force_n, 0.0);
            }
            if refresh.response().regime == DmtRegime::Separated {
                assert_eq!(report.parameter_work_j, 0.0);
            }
            // The new source, not a copied endpoint force, drives the next
            // actual native mechanical step. Identity refresh is exact.
            let next = *owner.step(&[0.0; 4]).unwrap();
            let previous = *reference.step(&[0.0; 4]).unwrap();
            if t == 300.0 {
                assert_eq!(next, previous);
                assert_eq!(native_bits(&owner), native_bits(&reference));
            } else if refresh.response().regime != DmtRegime::Separated {
                assert_ne!(
                    next.normal_force_n.to_bits(),
                    previous.normal_force_n.to_bits()
                );
            }
        }
    }
}

#[test]
fn one_contact_batch_matches_scalar_dmt_refresh_and_preserves_same_receiver() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = dmt(&fixture);
    let jkr = thermal_source(&fixture);
    for approach in [1e-6, -0.5 * law.interaction_range_m()] {
        let mut single = dmt_owner(&fixture, &law, approach, 1.0);
        let mut batch = MultipleSourceNormalModalSystem::new(
            network(&jkr, true, 0.0),
            vec![(
                BoundNormalContact::new(law.clone().into(), column(), -approach).unwrap(),
                dmt_controls(&jkr, &law),
            )],
            batch_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let refresh = law
            .refresh_at_held_indentation(approach, fixture.cards(), temperatures(315.0))
            .unwrap();
        let a = single
            .apply_dmt_source_refresh(&refresh, &CancelGate::new_clock_free())
            .unwrap();
        let b = batch
            .apply_source_refreshes(
                &[NormalSourceRefresh::Dmt {
                    contact: 0,
                    refreshed: &refresh,
                }],
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(a.parameter_work_j, b.parameter_work_j);
        assert_eq!(
            a.contact_energy_before_j,
            b.contacts[0].change.contact_energy_before_j
        );
        assert_eq!(
            a.contact_energy_after_j,
            b.contacts[0].change.contact_energy_after_j
        );
        assert_eq!(a.energy_residual_j, b.energy_residual_j);
        assert_eq!(native_bits(&single), batch_bits(&batch));
        let a = *single.step(&[0.0; 4]).unwrap();
        let b = batch.step(&[0.0; 4]).unwrap().clone();
        close(a.normal_force_n, b.contacts[0].normal_force_n, 1e-12, 1e-9);
        close(
            a.approach_after_m,
            b.contacts[0].approach_after_m,
            1e-14,
            1e-9,
        );
    }
}

#[test]
fn mixed_dmt_jkr_subset_is_order_independent_and_keeps_untouched_source() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = dmt(&fixture);
    let jkr = thermal_source(&fixture);
    let gap = 0.4 * law.interaction_range_m();
    let make = || {
        MultipleSourceNormalModalSystem::new(
            network(&jkr, true, 0.0),
            vec![
                (
                    BoundNormalContact::new(jkr.clone().into(), column(), 0.0).unwrap(),
                    controls(&jkr),
                ),
                (
                    BoundNormalContact::new(
                        law.clone().into(),
                        column().into_iter().map(|x| 0.7 * x).collect(),
                        gap,
                    )
                    .unwrap(),
                    dmt_controls(&jkr, &law),
                ),
                (
                    BoundNormalContact::new(
                        law.clone().into(),
                        column().into_iter().map(|x| 0.2 * x).collect(),
                        -1e-6,
                    )
                    .unwrap(),
                    dmt_controls(&jkr, &law),
                ),
            ],
            batch_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let mut owner = make();
    let mut permuted = make();
    let d = law
        .refresh_at_held_indentation(-gap, fixture.cards(), temperatures(285.0))
        .unwrap();
    let j = jkr
        .refresh_at_held_indentation(
            jkr.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let dr = NormalSourceRefresh::Dmt {
        contact: 1,
        refreshed: &d,
    };
    let jr = NormalSourceRefresh::Jkr {
        contact: 0,
        refreshed: &j,
    };
    let bits = batch_bits(&owner);
    let untouched = owner.contact(2).unwrap().clone();
    let a = owner
        .apply_source_refreshes(&[dr, jr], &CancelGate::new_clock_free())
        .unwrap();
    let b = permuted
        .apply_source_refreshes(&[jr, dr], &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(
        a.contacts.iter().map(|c| c.contact).collect::<Vec<_>>(),
        [0, 1]
    );
    assert_eq!(batch_bits(&owner), bits);
    assert_eq!(owner.contact(2).unwrap(), &untouched);
    assert_eq!(
        owner.contact(1).unwrap().source().identity(),
        d.source().identity()
    );
    assert_eq!(
        owner.contact(0).unwrap().source().identity(),
        j.source().identity()
    );
    assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
    assert_eq!(
        owner.step(&[0.0; 4]).unwrap(),
        permuted.step(&[0.0; 4]).unwrap()
    );
    let current = batch_bits(&owner);
    let wrong_law = owner.apply_source_refreshes(
        &[NormalSourceRefresh::Dmt {
            contact: 0,
            refreshed: &d,
        }],
        &CancelGate::new_clock_free(),
    );
    assert!(wrong_law.is_err());
    assert_eq!(batch_bits(&owner), current);
}

#[test]
fn invalid_original_receipts_receiver_refusal_and_late_cancellation_are_atomic() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = dmt(&fixture);
    let approach = -0.4 * law.interaction_range_m();
    let mut owner = dmt_owner(&fixture, &law, approach, 1.0);
    let refresh = law
        .refresh_at_held_indentation(approach, fixture.cards(), temperatures(315.0))
        .unwrap();
    let stale = law
        .refresh_at_held_indentation(0.0, fixture.cards(), temperatures(315.0))
        .unwrap();
    let foreign = SourceDmtSphereContact::new(
        fixture.normal.clone(),
        0.011,
        law.interaction_range_m(),
        law.envelope(),
    )
    .unwrap();
    let foreign = foreign
        .refresh_at_held_indentation(approach, fixture.cards(), temperatures(315.0))
        .unwrap();
    let bits = native_bits(&owner);
    let identity = owner.contact().source().identity();
    let energy = owner.total_energy_j().unwrap().to_bits();
    for bad in [&stale, &foreign] {
        assert!(
            owner
                .apply_dmt_source_refresh(bad, &CancelGate::new_clock_free())
                .is_err()
        );
    }
    let refused: Result<(), SourceNormalModalError> = owner
        .apply_dmt_source_refresh_with_admission(
            &refresh,
            &CancelGate::new_clock_free(),
            |_, _| Err(SourceNormalModalError::Unresolved("thermal driver refused")),
        );
    assert!(refused.is_err());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), SourceNormalModalError> = owner
        .apply_dmt_source_refresh_with_admission(&refresh, &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(cancelled.is_err());
    assert_eq!(native_bits(&owner), bits);
    assert_eq!(owner.contact().source().identity(), identity);
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.samples_rendered(), 0);
    let held_load = law
        .refresh_at_held_load(0.01, fixture.cards(), temperatures(315.0))
        .unwrap();
    let mut loaded = dmt_owner(&fixture, &law, held_load.before().indentation_m, 1.0);
    assert!(
        loaded
            .apply_dmt_source_refresh(&held_load, &CancelGate::new_clock_free())
            .is_err()
    );
    owner
        .apply_dmt_source_refresh(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    let accepted = owner.total_energy_j().unwrap().to_bits();
    assert!(
        owner
            .apply_dmt_source_refresh(&refresh, &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), accepted);
}

#[test]
fn stronger_tail_refuses_at_touch_and_cutoff_before_any_receiver_publication() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = dmt(&fixture);
    let jkr = thermal_source(&fixture);
    let base = dmt_owner(&fixture, &law, 0.0, 1.0).compliance_m_per_n();
    let negative = -0.5 * law.evaluate(0.0).unwrap().one_sided_tangents_n_per_m[0];
    // Original tail margin is 0.2; cooling raises W by 45%, making the
    // original global finite-range load line inadmissible for this step.
    let scale = (0.8 / (base * negative)).sqrt();
    for approach in [0.0, -law.interaction_range_m()] {
        let mut owner = dmt_owner(&fixture, &law, approach, scale);
        let no_change = law
            .refresh_at_held_indentation(approach, fixture.cards(), temperatures(300.0))
            .unwrap();
        owner
            .apply_dmt_source_refresh(&no_change, &CancelGate::new_clock_free())
            .unwrap();
        let refresh = law
            .refresh_at_held_indentation(approach, fixture.cards(), temperatures(285.0))
            .unwrap();
        let bits = native_bits(&owner);
        let mut called = false;
        let refused: Result<(), SourceNormalModalError> = owner
            .apply_dmt_source_refresh_with_admission(
                &refresh,
                &CancelGate::new_clock_free(),
                |_, _| {
                    called = true;
                    Ok(())
                },
            );
        assert!(matches!(
            refused,
            Err(SourceNormalModalError::Unresolved(_))
        ));
        assert!(!called);
        assert_eq!(native_bits(&owner), bits);
        assert_eq!(owner.contact().source().identity(), law.identity());
        let mut batch = MultipleSourceNormalModalSystem::new(
            network(&jkr, true, 0.0),
            vec![(
                BoundNormalContact::new(
                    law.clone().into(),
                    column().into_iter().map(|x| x * scale).collect(),
                    -approach,
                )
                .unwrap(),
                dmt_controls(&jkr, &law),
            )],
            batch_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let bits = batch_bits(&batch);
        assert!(
            batch
                .apply_source_refreshes(
                    &[NormalSourceRefresh::Dmt {
                        contact: 0,
                        refreshed: &refresh
                    }],
                    &CancelGate::new_clock_free()
                )
                .is_err()
        );
        assert_eq!(batch_bits(&batch), bits);
        assert_eq!(
            batch.contact(0).unwrap().source().identity(),
            law.identity()
        );
    }
}

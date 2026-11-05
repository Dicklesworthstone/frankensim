//! G0/G3: atomic source work on the actual shared modal and material histories.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::{
    multiple::MultiContactConfig,
    source_normal::multiple::{MultipleSourceNormalModalSystem, temperature::NormalSourceRefresh},
};

fn batch_controls() -> MultiContactConfig {
    MultiContactConfig {
        max_contacts: 8,
        max_sweeps: 64,
        max_setup_terms: 32768,
    }
}
fn batch_native_bits(owner: &MultipleSourceNormalModalSystem) -> Vec<u64> {
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
fn pair_owner(
    fixture: &temperature_fixture::Fixture,
    law: &SourceMaugisSphereContact,
    connected: bool,
) -> MultipleSourceNormalModalSystem {
    let jkr = thermal_source(fixture);
    let mut maugis = maugis_potential(law);
    if connected {
        maugis = maugis.with_noncontact_continuation().unwrap();
    }
    let mut limits = controls(&jkr);
    let d0 = maugis.indentation_scale_m();
    limits.minimum_approach_m = maugis.noncontact_continuation().map_or_else(
        || law.minimum_attached_indentation_m(),
        |bounds| bounds.minimum_indentation_m,
    ) + 1e-3 * d0;
    limits.maximum_approach_m = 4.0 * d0;
    MultipleSourceNormalModalSystem::new(
        network(&jkr, true, 0.0),
        vec![
            (
                BoundNormalContact::new(jkr.clone().into(), column(), 0.0).unwrap(),
                controls(&jkr),
            ),
            (
                BoundNormalContact::new(
                    maugis.into(),
                    vec![
                        0.45 / det::sqrt(MA),
                        -0.6 / det::sqrt(MA),
                        -0.45 / det::sqrt(MB),
                        0.5 / det::sqrt(MB),
                    ],
                    -initial_approach(law, connected),
                )
                .unwrap(),
                limits,
            ),
        ],
        batch_controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn pair_maugis(owner: &MultipleSourceNormalModalSystem) -> &MaugisAttachedPotential {
    let NormalLawSource::MaugisAttached(source) = owner.contact(1).unwrap().source() else {
        panic!("missing original finite-zone source");
    };
    source
}

#[test]
fn mixed_refreshes_preserve_native_history_and_publish_canonical_original_work() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let jkr = thermal_source(&fixture);
    for (lambda, connected) in [(5.0, false), (0.2, true)] {
        let law = maugis_source(&fixture, lambda);
        for (jkr_t, maugis_t) in [(300.0, 300.0), (315.0, 285.0)] {
            let mut owner = pair_owner(&fixture, &law, connected);
            let mut ordered = pair_owner(&fixture, &law, connected);
            let mut reference = pair_owner(&fixture, &law, connected);
            for _ in 0..2 {
                owner.step(&[0.0; 4]).unwrap();
                ordered.step(&[0.0; 4]).unwrap();
                reference.step(&[0.0; 4]).unwrap();
            }
            let last = owner.last_frame().cloned().unwrap();
            let before = pair_maugis(&owner)
                .evaluate(last.contacts[1].approach_after_m)
                .unwrap();
            let jkr_refresh = jkr
                .refresh_at_held_indentation(
                    jkr.evaluate_attached(last.contacts[0].approach_after_m)
                        .unwrap()
                        .state(),
                    fixture.cards(),
                    temperatures(jkr_t),
                )
                .unwrap();
            let maugis_refresh = law
                .refresh_at_held_indentation(before, fixture.cards(), temperatures(maugis_t))
                .unwrap();
            let j = NormalSourceRefresh::Jkr {
                contact: 0,
                refreshed: &jkr_refresh,
            };
            let m = NormalSourceRefresh::Maugis {
                contact: 1,
                refreshed: &maugis_refresh,
            };
            let bits = batch_native_bits(&owner);
            let native = owner.network().last_frame().cloned();
            let network_energy = owner.network().total_energy_j().unwrap().to_bits();
            let total = owner.total_energy_j().unwrap().to_bits();
            let budget = pair_maugis(&owner).budget();
            let rows: Vec<_> = (0..2)
                .map(|i| {
                    (
                        owner.contact(i).unwrap().work_column().to_vec(),
                        owner.contact(i).unwrap().gap_m(),
                    )
                })
                .collect();
            let report = owner
                .apply_source_refreshes(&[m, j], &CancelGate::new_clock_free())
                .unwrap();
            let canonical = ordered
                .apply_source_refreshes(&[j, m], &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(report, canonical);
            assert_eq!(
                report
                    .contacts
                    .iter()
                    .map(|c| c.contact)
                    .collect::<Vec<_>>(),
                [0, 1]
            );
            assert_eq!(batch_native_bits(&owner), bits);
            assert_eq!(owner.network().last_frame(), native.as_ref());
            assert_eq!(owner.last_frame(), Some(&last));
            assert_eq!(
                owner.network().total_energy_j().unwrap().to_bits(),
                network_energy
            );
            assert_eq!(owner.samples_rendered(), 2);
            assert_eq!(report.sample, 2);
            assert_eq!(report.energy_before_j.to_bits(), total);
            assert_eq!(
                report.energy_after_j.to_bits(),
                owner.total_energy_j().unwrap().to_bits()
            );
            assert_eq!(
                owner.contact(0).unwrap().source().identity(),
                jkr_refresh.source().identity()
            );
            assert_eq!(
                owner.contact(1).unwrap().source().identity(),
                maugis_refresh.source().identity()
            );
            assert_eq!(
                report.parameter_work_j,
                jkr_refresh.parameter_work_j() + maugis_refresh.parameter_work_j()
            );
            assert_eq!(jkr_refresh.mechanical_work_j(), 0.0);
            assert_eq!(maugis_refresh.mechanical_work_j(), 0.0);
            assert_eq!(pair_maugis(&owner).budget(), budget);
            assert_eq!(
                pair_maugis(&owner).noncontact_continuation().is_some(),
                connected
            );
            assert_eq!(
                pair_maugis(&owner)
                    .evaluate(before.indentation_m())
                    .unwrap()
                    .phase(),
                before.phase()
            );
            for (i, (row, gap)) in rows.iter().enumerate() {
                assert_eq!(owner.contact(i).unwrap().work_column(), row);
                assert_eq!(owner.contact(i).unwrap().gap_m().to_bits(), gap.to_bits());
                assert!(
                    report.contacts[i].change.energy_residual_j.abs()
                        <= report.contacts[i].change.energy_tolerance_j
                );
                assert_eq!(
                    report.contacts[i].change.approach_m,
                    last.contacts[i].approach_after_m
                );
            }
            assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
            if connected {
                assert_eq!(
                    pair_maugis(&owner)
                        .evaluate(before.indentation_m())
                        .unwrap()
                        .friction_contact_state()
                        .unwrap()
                        .area_m2(),
                    0.0
                );
            }
            // The next actual native solve consumes the new laws while keeping
            // every original mechanical memory and work-conjugate row.
            let next = owner.step(&[0.0; 4]).unwrap().clone();
            let old = reference.step(&[0.0; 4]).unwrap().clone();
            assert!(next.contact_work_residual_j.abs() <= next.contact_work_tolerance_j);
            if jkr_t == 300.0 {
                assert_eq!(report.parameter_work_j, 0.0);
                assert_eq!(next, old);
                assert_eq!(batch_native_bits(&owner), batch_native_bits(&reference));
            } else {
                assert!(jkr_refresh.parameter_work_j() > 0.0);
                assert!(maugis_refresh.parameter_work_j() < 0.0);
                assert_ne!(
                    next.contacts[0].normal_force_n.to_bits(),
                    old.contacts[0].normal_force_n.to_bits()
                );
                assert_ne!(batch_native_bits(&owner), batch_native_bits(&reference));
            }
        }
    }
}

#[test]
fn single_and_subset_refreshes_retain_original_scalar_work_and_untouched_neighbours() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let jkr = thermal_source(&fixture);
    let mut scalar = thermal_owner(&jkr);
    let mut one = MultipleSourceNormalModalSystem::new(
        network(&jkr, true, 0.0),
        vec![(scalar.contact().clone(), scalar.config())],
        batch_controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let refresh = jkr
        .refresh_at_held_indentation(
            jkr.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let first = scalar
        .apply_jkr_source_refresh(&refresh, &CancelGate::new_clock_free())
        .unwrap();
    let batch = one
        .apply_source_refreshes(
            &[NormalSourceRefresh::Jkr {
                contact: 0,
                refreshed: &refresh,
            }],
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(batch.parameter_work_j, first.parameter_work_j);
    assert_eq!(batch.energy_before_j, first.energy_before_j);
    assert_eq!(batch.energy_after_j, first.energy_after_j);
    assert_eq!(
        batch.contacts[0].change.energy_residual_j,
        first.energy_residual_j
    );
    assert_eq!(batch.energy_tolerance_j, first.energy_tolerance_j);
    assert_eq!(batch_native_bits(&one), native_bits(&scalar));

    for (lambda, connected) in [(5.0, false), (0.2, true)] {
        let law = maugis_source(&fixture, lambda);
        let mut scalar = maugis_owner(&fixture, &law, connected);
        let mut one = MultipleSourceNormalModalSystem::new(
            network(&jkr, true, 0.0),
            vec![(scalar.contact().clone(), scalar.config())],
            batch_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let before = original_potential(&scalar)
            .evaluate(initial_approach(&law, connected))
            .unwrap();
        let refresh = law
            .refresh_at_held_indentation(before, fixture.cards(), temperatures(315.0))
            .unwrap();
        let first = scalar
            .apply_maugis_source_refresh(&refresh, &CancelGate::new_clock_free())
            .unwrap();
        let batch = one
            .apply_source_refreshes(
                &[NormalSourceRefresh::Maugis {
                    contact: 0,
                    refreshed: &refresh,
                }],
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(batch.parameter_work_j, first.parameter_work_j);
        assert_eq!(batch.energy_before_j, first.energy_before_j);
        assert_eq!(batch.energy_after_j, first.energy_after_j);
        assert_eq!(
            batch.contacts[0].change.energy_residual_j,
            first.energy_residual_j
        );
        assert_eq!(batch.energy_tolerance_j, first.energy_tolerance_j);
        assert_eq!(batch_native_bits(&one), native_bits(&scalar));
    }
    let law = maugis_source(&fixture, 0.2);
    let mut owner = pair_owner(&fixture, &law, true);
    let untouched = owner.contact(1).unwrap().clone();
    let before = owner.total_energy_j().unwrap();
    let report = owner
        .apply_source_refreshes(
            &[NormalSourceRefresh::Jkr {
                contact: 0,
                refreshed: &refresh,
            }],
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(owner.contact(1).unwrap(), &untouched);
    assert_eq!(report.contacts.len(), 1);
    assert_eq!(report.contacts[0].contact, 0);
    assert_eq!(report.parameter_work_j, refresh.parameter_work_j());
    assert_eq!(report.energy_before_j, before);
    assert_eq!(report.energy_after_j, owner.total_energy_j().unwrap());
    assert_eq!(
        report.energy_tolerance_j,
        report.contacts[0].change.energy_tolerance_j
    );
}

#[test]
fn stale_duplicate_wrong_owner_and_receiver_failures_publish_no_partial_batch() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let jkr = thermal_source(&fixture);
    let law = maugis_source(&fixture, 5.0);
    let mut owner = pair_owner(&fixture, &law, false);
    owner.step(&[0.0; 4]).unwrap();
    let last = owner.last_frame().cloned().unwrap();
    let native = owner.network().last_frame().cloned();
    let bits = batch_native_bits(&owner);
    let energy = owner.total_energy_j().unwrap().to_bits();
    let originals: Vec<_> = (0..2).map(|i| owner.contact(i).unwrap().clone()).collect();
    let jkr_refresh = jkr
        .refresh_at_held_indentation(
            jkr.evaluate_attached(last.contacts[0].approach_after_m)
                .unwrap()
                .state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let maugis_refresh = law
        .refresh_at_held_indentation(
            pair_maugis(&owner)
                .evaluate(last.contacts[1].approach_after_m)
                .unwrap(),
            fixture.cards(),
            temperatures(285.0),
        )
        .unwrap();
    let stale = jkr
        .refresh_at_held_indentation(
            jkr.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let foreign = maugis_source(&fixture, 4.5);
    let foreign_refresh = foreign
        .refresh_at_held_indentation(
            foreign
                .evaluate_contact_energy(last.contacts[1].approach_after_m)
                .unwrap(),
            fixture.cards(),
            temperatures(285.0),
        )
        .unwrap();
    let j = NormalSourceRefresh::Jkr {
        contact: 0,
        refreshed: &jkr_refresh,
    };
    let m = NormalSourceRefresh::Maugis {
        contact: 1,
        refreshed: &maugis_refresh,
    };
    for bad in [
        vec![],
        vec![j, j],
        vec![
            j,
            NormalSourceRefresh::Maugis {
                contact: 2,
                refreshed: &maugis_refresh,
            },
        ],
        vec![
            j,
            NormalSourceRefresh::Jkr {
                contact: 1,
                refreshed: &jkr_refresh,
            },
        ],
        vec![
            NormalSourceRefresh::Jkr {
                contact: 0,
                refreshed: &stale,
            },
            m,
        ],
        vec![
            j,
            NormalSourceRefresh::Maugis {
                contact: 1,
                refreshed: &foreign_refresh,
            },
        ],
    ] {
        assert!(
            owner
                .apply_source_refreshes(&bad, &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(batch_native_bits(&owner), bits);
        assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
        assert_eq!(owner.last_frame(), Some(&last));
        assert_eq!(owner.network().last_frame(), native.as_ref());
        for (i, original) in originals.iter().enumerate() {
            assert_eq!(owner.contact(i).unwrap(), original);
        }
    }
    let refused: Result<(), SourceNormalModalError> = owner.apply_source_refreshes_with_admission(
        &[m, j],
        &CancelGate::new_clock_free(),
        |report, candidates| {
            assert_eq!(report.contacts.len(), 2);
            assert_eq!(candidates[0].0, 0);
            assert_eq!(candidates[1].0, 1);
            assert_eq!(
                candidates[0].1.source().identity(),
                jkr_refresh.source().identity()
            );
            assert_eq!(
                candidates[1].1.source().identity(),
                maugis_refresh.source().identity()
            );
            Err(SourceNormalModalError::Unresolved(
                "temperature work receiver refused",
            ))
        },
    );
    assert!(refused.is_err());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), SourceNormalModalError> = owner
        .apply_source_refreshes_with_admission(&[j, m], &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(cancelled.is_err());
    assert!(owner.apply_source_refreshes(&[j, m], &gate).is_err());
    assert_eq!(batch_native_bits(&owner), bits);
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.last_frame(), Some(&last));
    assert_eq!(owner.network().last_frame(), native.as_ref());
    for (i, original) in originals.iter().enumerate() {
        assert_eq!(owner.contact(i).unwrap(), original);
    }
    owner
        .apply_source_refreshes(&[j, m], &CancelGate::new_clock_free())
        .unwrap();
    let accepted = owner.total_energy_j().unwrap().to_bits();
    assert!(
        owner
            .apply_source_refreshes(&[j, m], &CancelGate::new_clock_free())
            .is_err()
    );
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), accepted);
    assert_eq!(batch_native_bits(&owner), bits);
    assert_eq!(owner.samples_rendered(), 1);
}

#[test]
fn complete_successor_set_must_pass_joint_curvature_before_any_source_publishes() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let warm_at_zero = law
        .refresh_at_held_indentation(
            law.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let approach =
        warm_at_zero.source().minimum_attached_indentation_m() + 1e-4 * law.indentation_scale_m();
    let refresh = law
        .refresh_at_held_indentation(
            law.evaluate_attached(approach).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let old: NormalLawSource = law.clone().into();
    let next: NormalLawSource = refresh.source().clone().into();
    let old_negative = (-old
        .mean_force_end_derivative(approach, approach)
        .unwrap()
        .unwrap())
    .max(0.0);
    let new_negative = -next
        .mean_force_end_derivative(approach, approach)
        .unwrap()
        .unwrap();
    assert!(new_negative > 0.0 && old_negative < 0.1 * new_negative);
    let probe = SourceNormalModalSystem::new(
        network(&law, false, 0.0),
        BoundNormalContact::new(law.clone().into(), column(), -approach).unwrap(),
        controls(&law),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let scale = (0.6 / (new_negative * probe.compliance_m_per_n())).sqrt();
    let row: Vec<_> = column().iter().map(|b| b * scale).collect();
    let make = || {
        MultipleSourceNormalModalSystem::new(
            network(&law, false, 0.0),
            (0..2)
                .map(|_| {
                    (
                        BoundNormalContact::new(law.clone().into(), row.clone(), -approach)
                            .unwrap(),
                        controls(&law),
                    )
                })
                .collect(),
            batch_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let mut partial = make();
    let first = NormalSourceRefresh::Jkr {
        contact: 0,
        refreshed: &refresh,
    };
    let second = NormalSourceRefresh::Jkr {
        contact: 1,
        refreshed: &refresh,
    };
    // Each isolated replacement is admissible with the unchanged neighbour.
    partial
        .apply_source_refreshes(&[first], &CancelGate::new_clock_free())
        .unwrap();
    let mut other = make();
    other
        .apply_source_refreshes(&[second], &CancelGate::new_clock_free())
        .unwrap();
    let mut owner = make();
    let s = owner.compliance_m_per_n();
    assert!(1.0 - new_negative * s[0] > 0.3);
    assert!(new_negative * (s[0] + s[1]) > 1.1);
    let bits = batch_native_bits(&owner);
    let energy = owner.total_energy_j().unwrap().to_bits();
    assert!(matches!(
        owner.apply_source_refreshes(&[first, second], &CancelGate::new_clock_free()),
        Err(SourceNormalModalError::Unresolved(
            "coupled adhesive negative curvature needs a smaller native step"
        ))
    ));
    for i in 0..2 {
        assert_eq!(
            owner.contact(i).unwrap().source().identity(),
            law.identity()
        );
    }
    assert_eq!(batch_native_bits(&owner), bits);
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
}

#[test]
fn shifted_source_fold_refuses_the_complete_batch_without_inventing_a_snap() {
    let fixture = temperature_fixture::Fixture::with_options(10.0, 0.05, false, 320.0);
    let law = thermal_source(&fixture);
    let refresh = law
        .refresh_at_held_indentation(
            law.evaluate_attached(0.0).unwrap().state(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    let held = controls(&law).minimum_approach_m + 1e-3 * law.indentation_scale_m();
    assert!(held < refresh.source().minimum_attached_indentation_m());
    let mut owner = MultipleSourceNormalModalSystem::new(
        network(&law, true, 0.0),
        [0.0, held]
            .into_iter()
            .map(|approach| {
                (
                    BoundNormalContact::new(law.clone().into(), column(), -approach).unwrap(),
                    controls(&law),
                )
            })
            .collect(),
        batch_controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let bits = batch_native_bits(&owner);
    let energy = owner.total_energy_j().unwrap().to_bits();
    // The first replacement is valid. The second original contact still
    // exists, but the successor has no attached solution at its held approach.
    assert!(matches!(
        owner.apply_source_refreshes(
            &[
                NormalSourceRefresh::Jkr {
                    contact: 0,
                    refreshed: &refresh
                },
                NormalSourceRefresh::Jkr {
                    contact: 1,
                    refreshed: &refresh
                },
            ],
            &CancelGate::new_clock_free(),
        ),
        Err(SourceNormalModalError::Source(_))
    ));
    for i in 0..2 {
        assert_eq!(
            owner.contact(i).unwrap().source().identity(),
            law.identity()
        );
    }
    assert_eq!(batch_native_bits(&owner), bits);
    assert_eq!(owner.total_energy_j().unwrap().to_bits(), energy);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
}

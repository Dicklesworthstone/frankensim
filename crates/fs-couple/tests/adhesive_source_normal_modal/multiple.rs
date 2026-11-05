//! Joint native contact regressions: source work, reciprocal cross response,
//! redundant rows, material memory, unstable attraction and rollback (G0/G3).
use super::*;
use fs_couple::render::schedule::force::coupled::contact::{
    multiple::MultiContactConfig, source_normal::multiple::MultipleSourceNormalModalSystem,
};

fn batch_config() -> MultiContactConfig {
    MultiContactConfig {
        max_contacts: 8,
        max_sweeps: 64,
        max_setup_terms: 32768,
    }
}
fn second_column() -> Vec<f64> {
    vec![
        0.45 / det::sqrt(MA),
        -0.6 / det::sqrt(MA),
        -0.45 / det::sqrt(MB),
        0.5 / det::sqrt(MB),
    ]
}
fn mixed_sources() -> Vec<(BoundNormalContact, SourceNormalModalConfig)> {
    let jkr = source();
    let dmt = normal_fixtures::dmt();
    let mut dmt_config = controls(&jkr);
    dmt_config.minimum_approach_m = -2.0 * dmt.interaction_range_m();
    dmt_config.maximum_approach_m = dmt.envelope().maximum_indentation_m;
    vec![
        (
            BoundNormalContact::new(
                jkr.clone().into(),
                column(),
                0.05 * jkr.indentation_scale_m(),
            )
            .unwrap(),
            controls(&jkr),
        ),
        (
            BoundNormalContact::new(
                dmt.clone().into(),
                second_column(),
                0.4 * dmt.interaction_range_m(),
            )
            .unwrap(),
            dmt_config,
        ),
    ]
}
fn mixed_owner() -> MultipleSourceNormalModalSystem {
    MultipleSourceNormalModalSystem::new(
        network(&source(), true, 0.0),
        mixed_sources(),
        batch_config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn mixed_adhesion_uses_the_entire_cross_response_and_original_native_memory() {
    let mut owner = mixed_owner();
    let mut reference = network(&source(), true, 0.0);
    let original_memory = owner.network().material_memory().to_vec();
    let s = owner.compliance_m_per_n();
    assert!(s[1].abs() > 0.1 * (s[0] * s[3]).sqrt());
    close(s[1], s[2], s[0].max(s[3]), 1e-12);
    let initial = owner.total_energy_j().unwrap();
    let external = [1e-6, -2e-7, 3e-7, 1e-7];
    let mut losses = 0.0;
    let mut work = 0.0;
    for _ in 0..12 {
        let frame = owner.step(&external).unwrap().clone();
        assert!(frame.contacts.iter().all(|c| c.normal_force_n < 0.0));
        let mut held = external;
        for (i, contact) in frame.contacts.iter().enumerate() {
            let original = owner.contact(i).unwrap();
            let expected = original
                .source()
                .mean_force_between(contact.approach_before_m, contact.approach_after_m)
                .unwrap()
                .unwrap();
            close(contact.normal_force_n, expected, 0.01, 1e-10);
            for (f, b) in held.iter_mut().zip(original.work_column()) {
                *f -= b * contact.normal_force_n;
            }
            assert!(contact.constitutive_residual_n.abs() <= contact.force_tolerance_n);
            assert!(contact.contact_work_residual_j.abs() <= contact.contact_work_tolerance_j);
        }
        let native = reference.step(&held).unwrap().clone();
        for (actual, expected) in states(owner.network()).iter().zip(states(&reference)) {
            close(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                source().indentation_scale_m() * det::sqrt(MA),
                2e-9,
            );
            close(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
                source().indentation_scale_m() * det::sqrt(MA) * RATE as f64,
                2e-9,
            );
        }
        for (a, b) in owner
            .network()
            .material_memory()
            .iter()
            .zip(reference.material_memory())
        {
            close(*a, *b, source().indentation_scale_m() * 50.0, 2e-9);
        }
        close(
            frame.network_dissipation_j,
            native.total_dissipation_j(),
            1e-9,
            1e-10,
        );
        assert!(native.material_dissipation_j > 0.0);
        losses += frame.network_dissipation_j;
        work += frame.external_work_j;
        close(
            owner.total_energy_j().unwrap() + losses - work,
            initial,
            1e-8,
            2e-9,
        );
        assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert_ne!(owner.network().material_memory(), original_memory);
}

#[test]
fn redundant_adhesive_rows_are_legal_and_contact_permutation_preserves_motion() {
    let source = source();
    let h = source.indentation_scale_m();
    let contacts = || {
        vec![
            (
                BoundNormalContact::new(source.clone().into(), column(), 0.05 * h).unwrap(),
                controls(&source),
            ),
            (
                BoundNormalContact::new(source.clone().into(), column(), -0.03 * h).unwrap(),
                controls(&source),
            ),
        ]
    };
    let gate = CancelGate::new_clock_free();
    let mut a = MultipleSourceNormalModalSystem::new(
        network(&source, false, 0.0),
        contacts(),
        batch_config(),
        &gate,
    )
    .unwrap();
    let mut reversed = contacts();
    reversed.reverse();
    let mut b = MultipleSourceNormalModalSystem::new(
        network(&source, false, 0.0),
        reversed,
        batch_config(),
        &gate,
    )
    .unwrap();
    let s = a.compliance_m_per_n();
    assert_eq!(s[0], s[1]);
    assert_eq!(s[0], s[3]);
    for _ in 0..10 {
        let fa = a.step(&[0.0; 4]).unwrap().clone();
        let fb = b.step(&[0.0; 4]).unwrap().clone();
        for (pa, pb) in fa.contacts.iter().zip(fb.contacts.iter().rev()) {
            close(pa.normal_force_n, pb.normal_force_n, 0.01, 1e-10);
        }
        for (pa, pb) in states(a.network()).iter().zip(states(b.network())) {
            close(
                pa.displacement_m_sqrt_kg,
                pb.displacement_m_sqrt_kg,
                h * det::sqrt(MA),
                2e-9,
            );
            close(
                pa.velocity_m_sqrt_kg_per_s,
                pb.velocity_m_sqrt_kg_per_s,
                h * det::sqrt(MA) * RATE as f64,
                2e-9,
            );
        }
        close(
            fa.network_energy_j + fa.contact_energy_j,
            fb.network_energy_j + fb.contact_energy_j,
            1e-9,
            2e-9,
        );
    }
}

#[test]
fn one_contact_reduces_to_the_original_scalar_native_owner() {
    let source = source();
    let initial = -0.05 * source.indentation_scale_m();
    let mut single = owner(&source, true, 0.0, initial);
    let mut multiple = MultipleSourceNormalModalSystem::new(
        network(&source, true, 0.0),
        vec![(
            BoundNormalContact::new(source.clone().into(), column(), -initial).unwrap(),
            controls(&source),
        )],
        batch_config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..8 {
        let a = single.step(&[0.0; 4]).unwrap().clone();
        let b = multiple.step(&[0.0; 4]).unwrap().clone();
        close(a.normal_force_n, b.contacts[0].normal_force_n, 0.01, 1e-11);
        close(
            a.contact_work_on_network_j,
            b.contact_work_on_network_j,
            1e-9,
            1e-11,
        );
        for (pa, pb) in states(single.network())
            .iter()
            .zip(states(multiple.network()))
        {
            close(
                pa.displacement_m_sqrt_kg,
                pb.displacement_m_sqrt_kg,
                source.indentation_scale_m() * det::sqrt(MA),
                2e-10,
            );
        }
    }
}

#[test]
fn deformed_maugis_gap_attraction_coexists_with_a_compressive_hertz_contact() {
    use fs_contact::normal_patch::maugis::{
        MaugisEnvelope, MaugisWorkBudget, SourceMaugisSphereContact,
        continuation::MaugisContinuationBudget, mechanical::MaugisAttachedPotential,
    };
    let normal = source();
    let radius = normal.effective_radius_m();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let a0 = (core::f64::consts::PI * w * radius * radius / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    let f0 = core::f64::consts::PI * w * radius;
    let maugis = MaugisAttachedPotential::new(
        SourceMaugisSphereContact::new(
            normal.normal_state().clone(),
            radius,
            2.0 * d0 / (core::f64::consts::PI * 0.2),
            MaugisEnvelope {
                maximum_indentation_m: 4.0 * d0,
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
        .unwrap(),
        MaugisContinuationBudget {
            work: MaugisWorkBudget {
                absolute_tolerance_j: 1e-12 * f0 * d0,
                relative_tolerance: 1e-11,
                maximum_evaluations: 2049,
            },
            energy_absolute_tolerance_j: 1e-10 * f0 * d0,
            energy_relative_tolerance: 1e-10,
        },
    )
    .unwrap()
    .with_noncontact_continuation()
    .unwrap();
    let limits = maugis.noncontact_continuation().unwrap();
    let approach = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let mut config = controls(&normal);
    config.minimum_approach_m = limits.minimum_indentation_m
        + 0.1 * (limits.first_contact_indentation_m - limits.minimum_indentation_m);
    config.maximum_approach_m = 3.0 * d0;
    let mut hertz_config = config;
    hertz_config.minimum_approach_m = -d0;
    let mut owner = MultipleSourceNormalModalSystem::new(
        network(&normal, false, 0.0),
        vec![
            (
                BoundNormalContact::new(maugis.clone().into(), column(), -approach).unwrap(),
                config,
            ),
            (
                BoundNormalContact::new(
                    normal_fixtures::point(None).into(),
                    second_column(),
                    -0.2 * d0,
                )
                .unwrap(),
                hertz_config,
            ),
        ],
        batch_config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let initial = owner.total_energy_j().unwrap();
    let frame = owner.step(&[0.0; 4]).unwrap().clone();
    let attractive = &frame.contacts[0];
    assert!(attractive.normal_force_n < 0.0);
    assert!(attractive.approach_after_m < limits.first_contact_indentation_m);
    assert!(frame.contacts[1].normal_force_n > 0.0);
    let path = maugis
        .path(attractive.approach_before_m, attractive.approach_after_m)
        .unwrap();
    close(attractive.normal_force_n, path.mean_force_n, f0, 1e-9);
    close(attractive.source_work_j, path.work.work_j, f0 * d0, 1e-9);
    close(
        owner.total_energy_j().unwrap() + frame.network_dissipation_j,
        initial,
        f0 * d0,
        1e-8,
    );
}

#[test]
fn coupled_attraction_refuses_even_when_each_individual_diagonal_is_stable() {
    let dmt = normal_fixtures::dmt();
    let jkr = source();
    let mut config = controls(&jkr);
    config.minimum_approach_m = -2.0 * dmt.interaction_range_m();
    config.maximum_approach_m = dmt.envelope().maximum_indentation_m;
    let gap = 0.4 * dmt.interaction_range_m();
    let gate = CancelGate::new_clock_free();
    let probe = SourceNormalModalSystem::new(
        network(&jkr, false, 0.0),
        BoundNormalContact::new(dmt.clone().into(), column(), gap).unwrap(),
        config,
        &gate,
    )
    .unwrap();
    let negative = -0.5 * dmt.evaluate(0.0).unwrap().one_sided_tangents_n_per_m[0];
    let scale = (0.6 / (negative * probe.compliance_m_per_n())).sqrt();
    let row: Vec<_> = column().iter().map(|b| b * scale).collect();
    let mut owner = MultipleSourceNormalModalSystem::new(
        network(&jkr, false, 0.0),
        (0..2)
            .map(|_| {
                (
                    BoundNormalContact::new(dmt.clone().into(), row.clone(), gap).unwrap(),
                    config,
                )
            })
            .collect(),
        batch_config(),
        &gate,
    )
    .unwrap();
    let s = owner.compliance_m_per_n();
    assert!(1.0 - negative * s[0] > 0.3);
    assert!(negative * (s[0] + s[1]) > 1.1);
    let before = states(owner.network());
    assert!(matches!(
        owner.step(&[0.0; 4]),
        Err(SourceNormalModalError::Unresolved(
            "coupled adhesive negative curvature needs a smaller native step"
        ))
    ));
    assert_eq!(states(owner.network()), before);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
}

#[test]
fn late_receiver_failure_and_cancellation_keep_all_native_contact_histories() {
    let mut owner = mixed_owner();
    owner.step(&[0.0; 4]).unwrap();
    let before = states(owner.network());
    let memory = owner.network().material_memory().to_vec();
    let last = owner.last_frame().unwrap().clone();
    let gate = CancelGate::new_clock_free();
    let rejected: Result<(), SourceNormalModalError> =
        owner.step_with_admission(&[0.0; 4], &gate, |_, _| {
            Err(SourceNormalModalError::Unresolved(
                "receiving thermal state refused",
            ))
        });
    assert!(rejected.is_err());
    let cancelled: Result<(), SourceNormalModalError> =
        owner.step_with_admission(&[0.0; 4], &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(cancelled.is_err());
    assert_eq!(states(owner.network()), before);
    assert_eq!(owner.network().material_memory(), memory);
    assert_eq!(owner.last_frame(), Some(&last));
    assert_eq!(owner.samples_rendered(), 1);
    // Retrying with a fresh gate advances the retained original owners once.
    owner.step(&[0.0; 4]).unwrap();
    assert_eq!(owner.samples_rendered(), 2);
}

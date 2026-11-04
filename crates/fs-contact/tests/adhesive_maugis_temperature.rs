//! Finite-zone thermal changes preserve actual mechanical geometry and branch policy.
#[path = "adhesive_temperature/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_contact::normal_patch::maugis::continuation::{MaugisContinuationBudget, MaugisPhase};
use fs_contact::normal_patch::maugis::mechanical::MaugisAttachedPotential;
use fs_contact::normal_patch::maugis::{
    MaugisEnvelope, MaugisWorkBudget, SourceMaugisSphereContact,
};

fn temperatures(t: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 300.0 + 0.5 * (t - 300.0),
        surface_b_k: 300.0 - 0.25 * (t - 300.0),
        interface_k: t,
    }
}
fn scales(law: &SourceMaugisSphereContact) -> (f64, f64) {
    let r = law.effective_radius_m();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let w = law.normal_state().adhesion_energy_j_per_m2();
    let a0 = (PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    (a0 * a0 / r, PI * w * r)
}
fn model(fixture: &common::Fixture, lambda: f64) -> SourceMaugisSphereContact {
    let r: f64 = 0.01;
    let e = fixture.normal.elastic().reduced_modulus_pa();
    let w = fixture.normal.adhesion_energy_j_per_m2();
    let a0 = (PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    SourceMaugisSphereContact::new(
        fixture.normal.clone(),
        r,
        2.0 * d0 / (PI * lambda),
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
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
fn potential(source: &SourceMaugisSphereContact) -> MaugisAttachedPotential {
    let (d0, f0) = scales(source);
    MaugisAttachedPotential::new(
        source.clone(),
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
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * scale, "{a:e} != {b:e}");
}

#[test]
fn held_source_refresh_preserves_work_budgets_geometry_policy_and_exact_noop() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    for (lambda, connected) in [(0.2, false), (0.2, true), (5.0, false)] {
        let source = model(&fixture, lambda);
        let (d0, f0) = scales(&source);
        let mut owner = potential(&source);
        if connected {
            owner = owner.with_noncontact_continuation().unwrap();
        }
        let original = owner.clone();
        let approach = 0.3 * d0;
        let before = owner.evaluate(approach).unwrap();
        for t in [300.0, 285.0, 315.0] {
            let refresh = source
                .refresh_at_held_indentation(before, fixture.cards(), temperatures(t))
                .unwrap();
            let next = owner.with_held_indentation_refresh(approach, &refresh).unwrap();
            let after = next.evaluate(approach).unwrap();
            assert_eq!(next.source(), refresh.source());
            assert_eq!(next.budget(), owner.budget());
            assert_eq!(next.noncontact_continuation().is_some(), connected);
            assert_eq!(after.phase(), MaugisPhase::Contact);
            assert_eq!(next.source().effective_radius_m(), source.effective_radius_m());
            assert_eq!(next.source().interaction_range_m(), source.interaction_range_m());
            assert_eq!(next.source().envelope(), source.envelope());
            close(
                after.potential_j() - before.potential_j(),
                refresh.parameter_work_j(),
                f0 * d0,
                1e-10,
            );
            assert_eq!(refresh.mechanical_work_j(), 0.0);
            let end = approach + 0.1 * d0;
            let path = next.path(approach, end).unwrap();
            close(
                next.evaluate(end).unwrap().potential_j() - after.potential_j(),
                path.work.work_j,
                f0 * d0,
                1e-10,
            );
            assert!(next.path_end_derivative(approach, end).unwrap().is_finite());
            if t == 300.0 {
                assert_eq!(next, owner);
                assert_eq!(after, before);
                assert_eq!(refresh.parameter_work_j(), 0.0);
            } else {
                assert_ne!(next.identity(), owner.identity());
                assert_ne!(after.normal_force_n(), before.normal_force_n());
            }
            assert_eq!(owner, original);
        }
    }
}

#[test]
fn heated_noncontact_disk_keeps_its_original_force_gap_and_zero_solid_area() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    let source = model(&fixture, 0.2);
    let (d0, f0) = scales(&source);
    let owner = potential(&source).with_noncontact_continuation().unwrap();
    let limits = owner.noncontact_continuation().unwrap();
    let approach = 0.25 * limits.minimum_indentation_m
        + 0.75 * limits.first_contact_indentation_m;
    let before = owner.evaluate(approach).unwrap();
    for t in [300.0, 315.0] {
        let refresh = source
            .refresh_at_held_indentation(
                before,
                fixture.cards(),
                AdhesiveTemperatures {
                    surface_a_k: 300.0,
                    surface_b_k: 300.0,
                    interface_k: t,
                },
            )
            .unwrap();
        let next = owner.with_held_indentation_refresh(approach, &refresh).unwrap();
        let after = next.evaluate(approach).unwrap();
        assert_eq!(after.phase(), MaugisPhase::Noncontact);
        assert!(after.central_gap_m() > 0.0);
        assert!(after.normal_force_n() < 0.0);
        assert_eq!(after.friction_contact_state().unwrap().area_m2(), 0.0);
        assert_eq!(next.budget(), owner.budget());
        close(
            after.potential_j() - before.potential_j(),
            refresh.parameter_work_j(),
            f0 * d0,
            1e-10,
        );
        assert!(next.path_end_derivative(approach, approach + 0.01 * d0).unwrap() < 0.0);
        if t == 300.0 {
            assert_eq!(next, owner);
        } else {
            assert_ne!(next.noncontact_continuation(), owner.noncontact_continuation());
        }
        assert!(potential(&source).with_held_indentation_refresh(approach, &refresh).is_err());
    }
}

#[test]
fn loss_of_selected_noncontact_continuation_cannot_silently_reduce_to_contact_only() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    let source = model(&fixture, 0.8);
    let contacting = potential(&source);
    let connected = contacting.clone().with_noncontact_continuation().unwrap();
    let before = connected.evaluate(0.0).unwrap();
    let refresh = source
        .refresh_at_held_indentation(
            before,
            fixture.cards(),
            AdhesiveTemperatures {
                surface_a_k: 300.0,
                surface_b_k: 300.0,
                interface_k: 280.0,
            },
        )
        .unwrap();
    assert!(refresh.source().noncontact_limits().unwrap().is_none());
    assert!(contacting.with_held_indentation_refresh(0.0, &refresh).is_ok());
    let saved = connected.clone();
    assert!(connected.with_held_indentation_refresh(0.0, &refresh).is_err());
    assert_eq!(connected, saved);
}

#[test]
fn stale_foreign_replayed_or_load_controlled_refresh_cannot_change_a_mechanical_potential() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    let source = model(&fixture, 5.0);
    let (d0, _) = scales(&source);
    let owner = potential(&source);
    let saved = owner.clone();
    let approach = 0.3 * d0;
    let refresh = source
        .refresh_at_held_indentation(
            owner.evaluate(approach).unwrap(),
            fixture.cards(),
            temperatures(315.0),
        )
        .unwrap();
    assert!(owner.with_held_indentation_refresh(approach + 0.1 * d0, &refresh).is_err());
    assert!(potential(&model(&fixture, 4.5))
        .with_held_indentation_refresh(approach, &refresh)
        .is_err());
    let load = source.refresh_at_held_load(0.0, fixture.cards(), temperatures(315.0)).unwrap();
    assert!(owner.with_held_indentation_refresh(load.before().indentation_m(), &load).is_err());
    let next = owner.with_held_indentation_refresh(approach, &refresh).unwrap();
    assert!(next.with_held_indentation_refresh(approach, &refresh).is_err());
    let detached = source.evaluate_detached(-2.0 * source.interaction_range_m()).unwrap();
    let detached_refresh = source
        .refresh_at_held_indentation(detached, fixture.cards(), temperatures(315.0))
        .unwrap();
    assert!(owner
        .with_held_indentation_refresh(detached.indentation_m(), &detached_refresh)
        .is_err());
    assert_eq!(owner, saved);
}

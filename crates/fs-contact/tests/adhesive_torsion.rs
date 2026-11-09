//! G0/G1/G3: original unequal bulk compliances, tensile JKR area and twist work.
#[path = "source_point/common.rs"]
mod common;

#[path = "adhesive_torsion/events.rs"]
mod events;

#[path = "adhesive_torsion/temperature.rs"]
mod temperature;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::torsion::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_tribo::InputAuthority;

fn normal(label: &str) -> SourceJkrSphereContact {
    let a = common::solid(label, 1000.0, 1e7);
    let b = common::solid("torsion-B", 2000.0, 2e9);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn controls() -> JkrTorsionControls {
    JkrTorsionControls {
        critical_twist_angle_rad: 1e-5,
        maximum_critical_twist_angle_rad: 0.01,
        maximum_mode_iii_energy_release_to_adhesion: 0.01,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_tensile_area_and_both_unequal_bulk_compliances_set_torsion() {
    let normal = normal("torsion-A");
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let fc = 1.5 * PI * w * normal.effective_radius_m();
    let elastic = normal.normal_state().elastic();
    let g = [elastic.surface_a_state(), elastic.surface_b_state()]
        .map(|s| s.elastic_state().shear_modulus_pa());
    assert!(g[1] > 100.0 * g[0]);
    let mut stiffness = Vec::new();
    for load in [0.0, -0.5 * fc, -fc] {
        let response = normal.evaluate_force_controlled(load).unwrap();
        let source =
            SourceJkrTorsionContact::new(normal.clone(), response.indentation_m, controls())
                .unwrap();
        // Independent pair of circular-punch compliances in series.
        let a = response.contact_radius_m;
        let compliance = 3.0 / (16.0 * a.powi(3)) * (1.0 / g[0] + 1.0 / g[1]);
        let k = compliance.recip();
        close(source.law().stiffness_nm_per_rad(), k);
        stiffness.push(k);
        assert_eq!(
            source.normal_response().normal_state_identity(),
            normal.normal_state().identity()
        );
        assert_eq!(
            source.normal_response().authority,
            InputAuthority::SyntheticFixture
        );
        assert!((source.normal_response().normal_force_n - load).abs() < 1e-11 * fc);
        let phi = 0.25 * controls().critical_twist_angle_rad;
        let step = source.advance(&source.initial_state(), phi).unwrap();
        close(step.torsion().endpoint_moment_n_m, -k * phi);
        close(step.torsion().storage_after_j, 0.5 * k * phi * phi);
        assert_eq!(step.torsion().irreversible_loss_j, 0.0);
        assert!(step.torsion().storage_after_j > 0.0);
    }
    // JKR force pull-off has one quarter of the zero-load contact volume a^3.
    close(stiffness[2] / stiffness[0], 0.25);
}

#[test]
fn original_history_yields_unloads_and_cannot_be_rebound_to_changed_area() {
    let original = normal("torsion-A");
    let indentation = original
        .evaluate_force_controlled(0.0)
        .unwrap()
        .indentation_m;
    let source = SourceJkrTorsionContact::new(original.clone(), indentation, controls()).unwrap();
    let c = controls().critical_twist_angle_rad;
    let moment = source.law().critical_moment_n_m();
    let step = source.advance(&source.initial_state(), 3.0 * c).unwrap();
    close(step.state().elastic_twist_rad(), c);
    close(step.torsion().irreversible_loss_j, 2.0 * c * moment);
    close(
        step.torsion().applied_work_j,
        -3.0 * c * step.torsion().mean_moment_n_m,
    );
    let unload = source.advance(&step.state(), -c).unwrap();
    assert_eq!(unload.torsion().irreversible_loss_j, 0.0);
    close(
        unload.torsion().applied_work_j,
        -step.torsion().storage_after_j,
    );
    let reverse = source.advance(&step.state(), -3.0 * c).unwrap();
    close(reverse.torsion().irreversible_loss_j, c * moment);
    let saved = step.state();
    for other in [
        SourceJkrTorsionContact::new(original.clone(), indentation * 1.1, controls()).unwrap(),
        SourceJkrTorsionContact::new(normal("foreign-A"), indentation, controls()).unwrap(),
        SourceJkrTorsionContact::new(
            original,
            indentation,
            JkrTorsionControls {
                critical_twist_angle_rad: 2.0 * c,
                ..controls()
            },
        )
        .unwrap(),
    ] {
        assert!(matches!(
            other.advance(&saved, 0.0),
            Err(SourceJkrTorsionError::StateMismatch)
        ));
        assert!(other.stored_energy_j(&saved).is_err());
    }
    assert!(source.advance(&saved, f64::NAN).is_err());
    assert_eq!(saved, step.state());
    assert_eq!(source.advance(&saved, 0.0).unwrap().state(), saved);
}

#[test]
fn frozen_area_release_screen_is_explicit_and_not_an_added_energy_channel() {
    let normal = normal("torsion-A");
    let indentation = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrTorsionContact::new(normal.clone(), indentation, controls()).unwrap();
    let at_yield = source
        .advance(&source.initial_state(), controls().critical_twist_angle_rad)
        .unwrap();
    // U is proportional to a^3, so dU/dA=3U/(2A) at fixed elastic twist.
    let release =
        1.5 * at_yield.torsion().storage_after_j / source.normal_response().contact_area_m2;
    close(source.maximum_mode_iii_release_j_m2(), release);
    close(
        at_yield.torsion().applied_work_j,
        at_yield.torsion().storage_after_j,
    );
    assert_eq!(at_yield.torsion().irreversible_loss_j, 0.0);
    let ratio = release / normal.normal_state().adhesion_energy_j_per_m2();
    assert!(ratio > 0.0 && ratio < controls().maximum_mode_iii_energy_release_to_adhesion);
    assert!(matches!(
        SourceJkrTorsionContact::new(
            normal.clone(),
            indentation,
            JkrTorsionControls {
                maximum_mode_iii_energy_release_to_adhesion: 0.5 * ratio,
                ..controls()
            }
        ),
        Err(SourceJkrTorsionError::OutsideFrozenArea { .. })
    ));
    for limit in [0.0, 1.0, f64::NAN] {
        assert!(
            SourceJkrTorsionContact::new(
                normal.clone(),
                indentation,
                JkrTorsionControls {
                    maximum_mode_iii_energy_release_to_adhesion: limit,
                    ..controls()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn prescribed_normal_unloading_preserves_loaded_twist_and_returns_parameter_work() {
    let normal = normal("torsion-normal-A");
    let zero = normal.evaluate_force_controlled(0.0).unwrap();
    let pull = normal
        .evaluate_force_controlled(normal.force_controlled_pull_off_n())
        .unwrap();
    let source =
        SourceJkrTorsionContact::new(normal.clone(), zero.indentation_m, controls()).unwrap();
    let c = controls().critical_twist_angle_rad;
    for angle in [0.4 * c, -3.0 * c] {
        let loaded = source.advance(&source.initial_state(), angle).unwrap();
        let step = source
            .advance_normal(&loaded.state(), pull.indentation_m)
            .unwrap();
        assert_eq!(step.source().normal_source().identity(), normal.identity());
        assert_eq!(step.source().controls(), controls());
        assert_eq!(
            step.state().elastic_twist_rad(),
            loaded.state().elastic_twist_rad()
        );
        assert_eq!(
            step.state().accumulated_angle_rad(),
            loaded.state().accumulated_angle_rad()
        );
        assert!(step.source().normal_response().normal_force_n < 0.0);
        close(
            step.source().law().stiffness_nm_per_rad(),
            0.25 * source.law().stiffness_nm_per_rad(),
        );
        close(
            step.torsion_storage_before_j(),
            loaded.torsion().storage_after_j,
        );
        close(
            step.torsion_storage_after_j(),
            0.25 * loaded.torsion().storage_after_j,
        );
        close(
            step.torsion_parameter_work_j(),
            -0.75 * loaded.torsion().storage_after_j,
        );
        close(
            step.normal_contact_work_j(),
            normal
                .reversible_work_between(zero.indentation_m, pull.indentation_m)
                .unwrap(),
        );
        close(
            step.applied_work_j(),
            step.normal_contact_work_j() + step.torsion_parameter_work_j(),
        );
        // An unchanged source cannot silently adopt the changed-radius history.
        assert!(matches!(
            source.advance(&step.state(), 0.0),
            Err(SourceJkrTorsionError::StateMismatch)
        ));
        let restore = step
            .source()
            .advance_normal(&step.state(), zero.indentation_m)
            .unwrap();
        assert_eq!(restore.source(), &source);
        assert_eq!(restore.state(), loaded.state());
        assert_eq!(
            restore.torsion_parameter_work_j(),
            -step.torsion_parameter_work_j()
        );
        close(
            restore.normal_contact_work_j(),
            -step.normal_contact_work_j(),
        );
        // Only actual twist at the new area can produce additional plastic loss.
        let next = step
            .source()
            .advance(&step.state(), angle.signum() * 2.0 * c)
            .unwrap();
        assert!(next.torsion().irreversible_loss_j > 0.0);
        close(
            next.torsion().endpoint_moment_n_m.abs(),
            step.source().law().critical_moment_n_m(),
        );
    }
}

#[test]
fn normal_continuation_is_subdivision_independent_and_keeps_weak_twist_history() {
    let normal = normal("torsion-history-A");
    let zero = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let end = normal
        .evaluate_force_controlled(normal.force_controlled_pull_off_n())
        .unwrap()
        .indentation_m;
    let middle = 0.5 * (zero + end);
    let source = SourceJkrTorsionContact::new(normal, zero, controls()).unwrap();
    let c = controls().critical_twist_angle_rad;
    let half = source.advance(&source.initial_state(), 0.5 * c).unwrap();
    let tiny = source.advance(&half.state(), 1e-20 * c).unwrap();
    assert_ne!(tiny.state(), half.state());
    assert_eq!(
        tiny.torsion().storage_after_j,
        half.torsion().storage_after_j
    );
    let direct = source.advance_normal(&tiny.state(), end).unwrap();
    let first = source.advance_normal(&tiny.state(), middle).unwrap();
    let second = first.source().advance_normal(&first.state(), end).unwrap();
    assert_eq!(second.source(), direct.source());
    assert_eq!(second.state(), direct.state());
    close(
        first.normal_contact_work_j() + second.normal_contact_work_j(),
        direct.normal_contact_work_j(),
    );
    close(
        first.torsion_parameter_work_j() + second.torsion_parameter_work_j(),
        direct.torsion_parameter_work_j(),
    );
    let restored = direct
        .source()
        .advance_normal(&direct.state(), zero)
        .unwrap();
    assert_eq!(restored.state(), tiny.state());
    let unchanged = source.advance_normal(&tiny.state(), zero).unwrap();
    assert_eq!(unchanged.source(), &source);
    assert_eq!(unchanged.state(), tiny.state());
    assert_eq!(unchanged.applied_work_j(), 0.0);
}

#[test]
fn changed_area_rechecks_omitted_release_and_refuses_without_losing_memory() {
    let normal = normal("torsion-envelope-A");
    let zero = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let target = 2.0 * zero;
    let initial = SourceJkrTorsionContact::new(normal.clone(), zero, controls()).unwrap();
    let enlarged = SourceJkrTorsionContact::new(normal.clone(), target, controls()).unwrap();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let ratio0 = initial.maximum_mode_iii_release_j_m2() / w;
    let ratio1 = enlarged.maximum_mode_iii_release_j_m2() / w;
    assert!(ratio1 > ratio0);
    let bounded = SourceJkrTorsionContact::new(
        normal.clone(),
        zero,
        JkrTorsionControls {
            maximum_mode_iii_energy_release_to_adhesion: 0.5 * (ratio0 + ratio1),
            ..controls()
        },
    )
    .unwrap();
    let loaded = bounded
        .advance(
            &bounded.initial_state(),
            2.0 * controls().critical_twist_angle_rad,
        )
        .unwrap();
    let saved = loaded.state();
    assert!(matches!(
        bounded.advance_normal(&saved, target),
        Err(SourceJkrTorsionError::OutsideFrozenArea { .. })
    ));
    for target in [
        f64::NAN,
        f64::INFINITY,
        1.1 * normal.minimum_attached_indentation_m(),
    ] {
        assert!(bounded.advance_normal(&saved, target).is_err());
    }
    let foreign = initial
        .advance(
            &initial.initial_state(),
            controls().critical_twist_angle_rad,
        )
        .unwrap();
    assert!(matches!(
        bounded.advance_normal(&foreign.state(), zero),
        Err(SourceJkrTorsionError::StateMismatch)
    ));
    assert_eq!(loaded.state(), saved);
    assert_eq!(bounded.advance(&saved, 0.0).unwrap().state(), saved);
}

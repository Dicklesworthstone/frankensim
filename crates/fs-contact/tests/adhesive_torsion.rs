//! G0/G1/G3: original unequal bulk compliances, tensile JKR area and twist work.
#[path = "source_point/common.rs"]
mod common;
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

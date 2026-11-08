//! G0/G1/G3: changing-area I/III equilibrium, coupled work and physical limits.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::*;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_tribo::InputAuthority;

fn source() -> SourceJkrMixedMode {
    let a = common::solid("torsion-equilibrium-A", 1000.0, 1e7);
    let b = common::solid("torsion-equilibrium-B", 2000.0, 2e9);
    let normal = SourceJkrSphereContact::new(
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
            maximum_nominal_strain: 0.3,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap();
    SourceJkrMixedMode::new(
        normal,
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

fn oracle(source: &SourceJkrMixedMode) -> (f64, f64, f64) {
    let normal = source.normal_source();
    let a = normal
        .evaluate_force_controlled(0.0)
        .unwrap()
        .contact_radius_m;
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let solids = normal.normal_state().elastic();
    let ga = solids.surface_a_state().elastic_state().shear_modulus_pa();
    let gb = solids.surface_b_state().elastic_state().shear_modulus_pa();
    assert!(gb > 100.0 * ga);
    let g = (ga.recip() + gb.recip()).recip();
    let phi = (0.2 * PI * w / (4.0 * g * a)).sqrt();
    let delta = a * a / normal.effective_radius_m() - (0.8 * 2.0 * PI * w * a / e).sqrt();
    (a, delta, phi)
}

#[test]
fn original_source_and_independent_i_iii_intensity_oracle_set_area_and_efforts() {
    let source = source();
    let normal = source.normal_source();
    let (a, delta, phi) = oracle(&source);
    let response = source
        .evaluate_torsion_displacement_controlled(delta, phi)
        .unwrap();
    let pure = normal.evaluate_attached(delta).unwrap();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let g = source.torsional_reduced_shear_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    close(response.contact_radius_m, a, 2e-13);
    assert!(response.contact_area_m2 < pure.contact_area_m2);
    assert!(response.normal_force_n > pure.normal_force_n);
    close(
        response.normal_force_n,
        2.0 * e * a * (delta - a * a / (3.0 * normal.effective_radius_m())),
        2e-12,
    );
    close(
        response.torsional_moment_n_m,
        (16.0 / 3.0) * g * a.powi(3) * phi,
        2e-12,
    );
    // Independent edge intensity from torque, with two-body compliance.
    let k_iii = 3.0 * response.torsional_moment_n_m / (4.0 * PI.sqrt() * a.powf(2.5));
    close(k_iii * k_iii / (4.0 * g), 0.2 * w, 2e-12);
    close(response.normal_energy_release_rate_j_m2, 0.8 * w, 2e-12);
    assert!(response.griffith_residual_j_m2.abs() < 1e-13 * w);
    assert_eq!(
        response.normal_state_identity(),
        normal.normal_state().identity()
    );
    assert_eq!(response.model_identity(), source.identity());
    assert_eq!(response.authority, InputAuthority::SyntheticFixture);
    for delta in [0.0, delta, 0.5 * normal.minimum_attached_indentation_m()] {
        let plain = normal.evaluate_attached(delta).unwrap();
        let zero = source
            .evaluate_torsion_displacement_controlled(delta, 0.0)
            .unwrap();
        close(zero.contact_radius_m, plain.contact_radius_m, 2e-13);
        close(zero.normal_force_n, plain.normal_force_n, 2e-12);
        close(zero.reversible_energy_j, plain.reversible_energy_j, 2e-12);
        assert_eq!(zero.torsional_moment_n_m, 0.0);
    }
}

#[test]
fn relaxed_energy_gradient_and_reciprocal_tangent_match_independent_differences() {
    let source = source();
    let (_, delta, phi) = oracle(&source);
    let center = source
        .evaluate_torsion_displacement_controlled(delta, phi)
        .unwrap();
    let tangent = center.coupled_stiffness.unwrap();
    assert_eq!(tangent[0][1], tangent[1][0]);
    for (column, increment) in [delta.abs() * 1e-5, phi * 1e-5].into_iter().enumerate() {
        let mut p = [delta, phi];
        let mut m = p;
        p[column] += increment;
        m[column] -= increment;
        let plus = source
            .evaluate_torsion_displacement_controlled(p[0], p[1])
            .unwrap();
        let minus = source
            .evaluate_torsion_displacement_controlled(m[0], m[1])
            .unwrap();
        let effort = [center.normal_force_n, center.torsional_moment_n_m][column];
        close(
            (plus.reversible_energy_j - minus.reversible_energy_j) / (2.0 * increment),
            effort,
            2e-7,
        );
        for (row, derivative) in [
            (plus.normal_force_n - minus.normal_force_n) / (2.0 * increment),
            (plus.torsional_moment_n_m - minus.torsional_moment_n_m) / (2.0 * increment),
        ]
        .into_iter()
        .enumerate()
        {
            close(derivative, tangent[row][column], 2e-7);
        }
    }
    let reversed = source
        .evaluate_torsion_displacement_controlled(delta, -phi)
        .unwrap();
    assert_eq!(reversed.contact_radius_m, center.contact_radius_m);
    assert_eq!(reversed.reversible_energy_j, center.reversible_energy_j);
    assert_eq!(reversed.normal_force_n, center.normal_force_n);
    assert_eq!(reversed.torsional_moment_n_m, -center.torsional_moment_n_m);
    assert_eq!(reversed.coupled_stiffness.unwrap()[0][1], -tangent[0][1]);
}

#[test]
fn actual_twist_work_includes_normal_and_surface_relaxation() {
    let source = source();
    let (_, delta, _) = oracle(&source);
    let end = 0.6
        * source
            .torsional_displacement_limit(delta)
            .unwrap()
            .unwrap()
            .twist_angle_rad;
    let initial = source
        .evaluate_torsion_displacement_controlled(delta, 0.0)
        .unwrap();
    let final_state = source
        .evaluate_torsion_displacement_controlled(delta, end)
        .unwrap();
    let mut integral = 0.0;
    for i in 0..=128 {
        let phi = end * i as f64 / 128.0;
        let response = source
            .evaluate_torsion_displacement_controlled(delta, phi)
            .unwrap();
        let weight = if i == 0 || i == 128 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * response.torsional_moment_n_m;
    }
    integral *= end / (3.0 * 128.0);
    close(
        integral,
        final_state.reversible_energy_j - initial.reversible_energy_j,
        2e-8,
    );
    assert!((integral - final_state.torsional_elastic_energy_j).abs() > 0.01 * integral);
}

#[test]
fn radial_detachment_fold_and_compressive_opening_edge_are_distinct() {
    let source = source();
    let normal = source.normal_source();
    let delta = (2.0 / 3.0) * normal.minimum_attached_indentation_m();
    let limit = source.torsional_displacement_limit(delta).unwrap().unwrap();
    assert_eq!(limit.kind, JkrTorsionalLimitKind::RadialFold);
    let mut previous = f64::INFINITY;
    for fraction in [0.0, 0.2, 0.6, 0.9, 0.99, 1.0] {
        let state = source
            .evaluate_torsion_displacement_controlled(delta, fraction * limit.twist_angle_rad)
            .unwrap();
        assert!(state.contact_radius_m < previous);
        assert!(state.normal_force_n < 0.0);
        previous = state.contact_radius_m;
    }
    let fold = source
        .evaluate_torsion_displacement_controlled(delta, limit.twist_angle_rad)
        .unwrap();
    assert_eq!(fold.boundary, Some(JkrTorsionalLimitKind::RadialFold));
    assert_eq!(fold.displacement_radius_curvature_n_per_m, 0.0);
    assert_eq!(fold.coupled_stiffness, None);
    assert!(matches!(
        source.evaluate_torsion_displacement_controlled(delta, limit.twist_angle_rad.next_up()),
        Err(JkrTorsionalError::BeyondAttachedLimit {
            kind: JkrTorsionalLimitKind::RadialFold,
            ..
        })
    ));

    let (_, delta, _) = oracle(&source);
    let limit = source.torsional_displacement_limit(delta).unwrap().unwrap();
    assert_eq!(limit.kind, JkrTorsionalLimitKind::OpeningEdge);
    let edge = source
        .evaluate_torsion_displacement_controlled(delta, limit.twist_angle_rad)
        .unwrap();
    close(
        edge.contact_radius_m,
        (normal.effective_radius_m() * delta).sqrt(),
        2e-13,
    );
    assert_eq!(edge.normal_energy_release_rate_j_m2, 0.0);
    assert!(edge.displacement_radius_curvature_n_per_m > 0.0);
    assert!(edge.coupled_stiffness.is_some());
    assert!(matches!(
        source.evaluate_torsion_displacement_controlled(delta, limit.twist_angle_rad.next_up()),
        Err(JkrTorsionalError::BeyondAttachedLimit {
            kind: JkrTorsionalLimitKind::OpeningEdge,
            ..
        })
    ));
    assert_eq!(source.torsional_displacement_limit(0.0).unwrap(), None);
}

#[test]
fn original_envelopes_and_normal_spinodal_remain_effective() {
    let source = source();
    let minimum = source.normal_source().minimum_attached_indentation_m();
    let limit = source
        .torsional_displacement_limit(minimum)
        .unwrap()
        .unwrap();
    assert_eq!(limit.twist_angle_rad, 0.0);
    let fold = source
        .evaluate_torsion_displacement_controlled(minimum, 0.0)
        .unwrap();
    assert_eq!(fold.boundary, Some(JkrTorsionalLimitKind::RadialFold));
    assert!(
        source
            .evaluate_torsion_displacement_controlled(minimum, 1e-12)
            .is_err()
    );
    assert!(
        source
            .evaluate_torsion_displacement_controlled(minimum.next_down(), 0.0)
            .is_err()
    );
    assert!(
        source
            .evaluate_torsion_displacement_controlled(0.0, 0.2)
            .is_err()
    );
    assert!(
        source
            .evaluate_torsion_displacement_controlled(0.0, f64::NAN)
            .is_err()
    );
    assert!(
        source
            .evaluate_torsion_displacement_controlled(f64::NAN, 0.0)
            .is_err()
    );
    let near = source
        .torsional_displacement_limit(minimum.next_up())
        .unwrap()
        .unwrap();
    assert!(near.twist_angle_rad > 0.0);
    let attached = source
        .evaluate_torsion_displacement_controlled(minimum.next_up(), 0.5 * near.twist_angle_rad)
        .unwrap();
    assert!(attached.displacement_radius_curvature_n_per_m > 0.0);

    // The ideal twist boundary diverges on either side of zero approach;
    // ordinary finite twists must still converge to the zero-approach state.
    let (_, _, phi) = oracle(&source);
    let zero = source
        .evaluate_torsion_displacement_controlled(0.0, phi)
        .unwrap();
    for delta in [minimum * 1e-13, -minimum * 1e-13] {
        let near_zero = source
            .evaluate_torsion_displacement_controlled(delta, phi)
            .unwrap();
        close(near_zero.contact_radius_m, zero.contact_radius_m, 2e-12);
        close(
            near_zero.torsional_moment_n_m,
            zero.torsional_moment_n_m,
            2e-12,
        );
    }
}

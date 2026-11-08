//! Joint original I/II/III equilibrium, complete tangent and applicability.
#[path = "source_point/common.rs"]
mod common;
#[path = "adhesive_combined/motion.rs"]
mod motion;
#[path = "adhesive_combined/path.rs"]
mod path;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrEnvelope, SourceJkrSphereContact,
    mixed_mode::{
        JkrMixedModeEnvelope, JkrMixedModeError, SourceJkrMixedMode, combined::JkrCombinedError,
    },
};
use fs_tribo::InputAuthority;

fn source() -> SourceJkrMixedMode {
    let a = common::solid("combined-equilibrium-A", 1000.0, 1e7);
    let b = common::solid("combined-equilibrium-B", 2000.0, 2e9);
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

// Manufacture the shared radius from the original zero-load normal solution.
// Assign independent normal/shear/twist Griffith fractions summing to one;
// this uses no combined radius solve or relaxed tangent.
fn manufactured(source: &SourceJkrMixedMode) -> (f64, [f64; 4]) {
    let normal = source.normal_source();
    let a = normal
        .evaluate_force_controlled(0.0)
        .unwrap()
        .contact_radius_m;
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let shear = (0.12 * PI * w * a / (2.0 * source.reduced_shear_modulus_pa())).sqrt();
    let twist = (0.18 * PI * w / (4.0 * source.torsional_reduced_shear_modulus_pa() * a)).sqrt();
    let delta = a * a / normal.effective_radius_m() - (0.70 * 2.0 * PI * w * a / e).sqrt();
    (a, [delta, 0.6 * shear, 0.8 * shear, twist])
}

#[test]
fn manufactured_all_mode_equilibrium_retains_one_original_area_and_potential() {
    let source = source();
    let normal = source.normal_source();
    let (a, q) = manufactured(&source);
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let g = source.reduced_shear_modulus_pa();
    let gt = source.torsional_reduced_shear_modulus_pa();
    let point = source.evaluate_combined_displacement_controlled(q).unwrap();
    let pure = normal.evaluate_attached(q[0]).unwrap();
    close(point.contact_radius_m, a, 3e-13);
    close(point.contact_area_m2, PI * a * a, 6e-13);
    assert!(point.contact_area_m2 < pure.contact_area_m2);
    assert_eq!(point.coordinates, q);
    assert_eq!(point.boundary, None);
    assert_eq!(point.model_identity(), source.identity());
    assert_eq!(
        point.normal_state_identity(),
        normal.normal_state().identity()
    );
    assert_eq!(point.authority, InputAuthority::SyntheticFixture);
    let hertz = a * a / r;
    let mean = q[0] - hertz / 3.0;
    let efforts = [
        2.0 * e * a * mean,
        8.0 * g * a * q[1],
        8.0 * g * a * q[2],
        (16.0 / 3.0) * gt * a.powi(3) * q[3],
    ];
    for (actual, expected) in point.force.into_iter().zip(efforts) {
        close(actual, expected, 2e-12);
    }
    let normal_energy = e * a * (mean * mean + 4.0 * hertz * hertz / 45.0);
    let shear_energy = 4.0 * g * a * (q[1] * q[1] + q[2] * q[2]);
    let torsion_energy = (8.0 / 3.0) * gt * a.powi(3) * q[3] * q[3];
    let surface_energy = -PI * w * a * a;
    close(point.normal_elastic_energy_j, normal_energy, 2e-12);
    close(point.tangential_elastic_energy_j, shear_energy, 2e-12);
    close(point.torsional_elastic_energy_j, torsion_energy, 2e-12);
    close(point.surface_energy_j, surface_energy, 2e-12);
    close(
        point.reversible_energy_j,
        normal_energy + shear_energy + torsion_energy + surface_energy,
        3e-12,
    );
    close(
        point.repulsive_normal_force_n,
        (4.0 / 3.0) * e * a.powi(3) / r,
        2e-12,
    );
    close(point.normal_energy_release_rate_j_m2, 0.70 * w, 2e-12);
    close(point.shear_energy_release_rate_j_m2, 0.12 * w, 2e-12);
    close(point.torsional_energy_release_rate_j_m2, 0.18 * w, 2e-12);
    assert!(point.griffith_residual_j_m2.abs() <= 2e-13 * w);
    assert!(point.radial_residual_n.abs() <= 2e-13 * 2.0 * PI * w * a);
    assert!(point.displacement_radius_curvature_n_per_m > 0.0);
}

#[test]
fn exact_mode_reductions_retain_original_responses_and_inactive_stiffness() {
    let source = source();
    let (_, q) = manufactured(&source);
    let shear = source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let all = source
        .evaluate_combined_displacement_controlled([q[0], q[1], q[2], 0.0])
        .unwrap();
    assert_eq!(all.contact_radius_m, shear.contact_radius_m);
    assert_eq!(all.reversible_energy_j, shear.reversible_energy_j);
    assert_eq!(
        all.force,
        [
            shear.normal_force_n,
            shear.tangential_force_n[0],
            shear.tangential_force_n[1],
            0.0
        ]
    );
    assert_eq!(all.torsional_elastic_energy_j, 0.0);
    let tangent = all.coupled_stiffness.unwrap();
    let original = shear.coupled_stiffness_n_per_m.unwrap();
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(tangent[i][j], original[i][j]);
        }
        assert_eq!(tangent[3][i], 0.0);
        assert_eq!(tangent[i][3], 0.0);
    }
    close(
        tangent[3][3],
        (16.0 / 3.0) * source.torsional_reduced_shear_modulus_pa() * all.contact_radius_m.powi(3),
        4e-15,
    );

    let torsion = source
        .evaluate_torsion_displacement_controlled(q[0], q[3])
        .unwrap();
    let all = source
        .evaluate_combined_displacement_controlled([q[0], 0.0, 0.0, q[3]])
        .unwrap();
    assert_eq!(all.contact_radius_m, torsion.contact_radius_m);
    assert_eq!(all.reversible_energy_j, torsion.reversible_energy_j);
    assert_eq!(
        all.force,
        [
            torsion.normal_force_n,
            0.0,
            0.0,
            torsion.torsional_moment_n_m
        ]
    );
    assert_eq!(all.tangential_elastic_energy_j, 0.0);
    let tangent = all.coupled_stiffness.unwrap();
    let original = torsion.coupled_stiffness.unwrap();
    for (i, row) in [0, 3].into_iter().enumerate() {
        for (j, column) in [0, 3].into_iter().enumerate() {
            assert_eq!(tangent[row][column], original[i][j]);
        }
    }
    for i in [1, 2] {
        close(
            tangent[i][i],
            8.0 * source.reduced_shear_modulus_pa() * all.contact_radius_m,
            4e-15,
        );
        for j in 0..4 {
            if i != j {
                assert_eq!(tangent[i][j], 0.0);
                assert_eq!(tangent[j][i], 0.0);
            }
        }
    }
    for delta in [
        0.0,
        q[0],
        0.5 * source.normal_source().minimum_attached_indentation_m(),
    ] {
        let normal = source.normal_source().evaluate_attached(delta).unwrap();
        let all = source
            .evaluate_combined_displacement_controlled([delta, 0.0, 0.0, 0.0])
            .unwrap();
        assert_eq!(all.contact_radius_m, normal.contact_radius_m);
        close(all.force[0], normal.normal_force_n, 3e-13);
        close(all.reversible_energy_j, normal.reversible_energy_j, 3e-13);
        assert_eq!(&all.force[1..], &[0.0; 3]);
    }
}

#[test]
fn every_relaxed_derivative_matches_force_and_energy_and_is_frame_covariant() {
    let source = source();
    let (_, q) = manufactured(&source);
    let center = source.evaluate_combined_displacement_controlled(q).unwrap();
    let tangent = center.coupled_stiffness.unwrap();
    for j in 0..4 {
        let increment = 1e-5 * q[j].abs();
        let mut plus = q;
        let mut minus = q;
        plus[j] += increment;
        minus[j] -= increment;
        let p = source
            .evaluate_combined_displacement_controlled(plus)
            .unwrap();
        let m = source
            .evaluate_combined_displacement_controlled(minus)
            .unwrap();
        close(
            (p.reversible_energy_j - m.reversible_energy_j) / (2.0 * increment),
            center.force[j],
            4e-7,
        );
        for i in 0..4 {
            assert_eq!(tangent[i][j], tangent[j][i]);
            close(
                (p.force[i] - m.force[i]) / (2.0 * increment),
                tangent[i][j],
                4e-7,
            );
        }
    }
    // A quarter turn changes only shear coordinates and their conjugate rows.
    // Twist reversal changes every odd twist coupling without changing area.
    for (permutation, signs) in [
        ([0, 2, 1, 3], [1.0, -1.0, 1.0, 1.0]),
        ([0, 1, 2, 3], [1.0, 1.0, 1.0, -1.0]),
    ] {
        let transformed = core::array::from_fn(|i| signs[i] * q[permutation[i]]);
        let response = source
            .evaluate_combined_displacement_controlled(transformed)
            .unwrap();
        close(response.contact_radius_m, center.contact_radius_m, 4e-14);
        close(
            response.reversible_energy_j,
            center.reversible_energy_j,
            2e-13,
        );
        let transformed_tangent = response.coupled_stiffness.unwrap();
        for i in 0..4 {
            close(
                response.force[i],
                signs[i] * center.force[permutation[i]],
                2e-13,
            );
            for j in 0..4 {
                close(
                    transformed_tangent[i][j],
                    signs[i] * signs[j] * tangent[permutation[i]][permutation[j]],
                    3e-13,
                );
            }
        }
    }
}

#[test]
fn joint_rim_and_yield_bounds_refuse_combination_when_each_original_mode_fits() {
    let broad = source();
    let (_, q) = manufactured(&broad);
    let point = broad.evaluate_combined_displacement_controlled(q).unwrap();
    let shear = broad
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let magnitude = q[1].hypot(q[2]);
    let yield_stress = broad
        .normal_source()
        .normal_state()
        .elastic()
        .limiting_yield_stress_pa();
    let g = broad.reduced_shear_modulus_pa();
    let gt = broad.torsional_reduced_shear_modulus_pa();
    let shear_rim = magnitude / shear.contact_radius_m;
    let twist_rim = q[3].abs();
    let combined_rim = magnitude / point.contact_radius_m + twist_rim;
    let shear_mean = 8.0 * g * magnitude / (PI * shear.contact_radius_m * yield_stress);
    let twist_mean = 2.0 * gt * twist_rim / yield_stress;
    let combined_mean =
        (8.0 * g * magnitude / (PI * point.contact_radius_m) + 2.0 * gt * twist_rim) / yield_stress;
    assert!(combined_rim > shear_rim.max(twist_rim));
    assert!(combined_mean > shear_mean.max(twist_mean));
    for (envelope, expected_ratio) in [
        (
            JkrMixedModeEnvelope {
                maximum_tangential_displacement_to_radius: shear_rim
                    .max(twist_rim)
                    .midpoint(combined_rim),
                ..broad.envelope()
            },
            "combined rim displacement/radius",
        ),
        (
            JkrMixedModeEnvelope {
                maximum_mean_shear_stress_to_yield: shear_mean
                    .max(twist_mean)
                    .midpoint(combined_mean),
                ..broad.envelope()
            },
            "combined mean shear/yield",
        ),
    ] {
        let limited = SourceJkrMixedMode::new(broad.normal_source().clone(), envelope).unwrap();
        limited
            .evaluate_displacement_controlled(q[0], [q[1], q[2]])
            .unwrap();
        limited
            .evaluate_torsion_displacement_controlled(q[0], q[3])
            .unwrap();
        assert!(
            matches!(limited.evaluate_combined_displacement_controlled(q),
            Err(JkrCombinedError::Source(JkrMixedModeError::OutsideApplicability { ratio, .. }))
                if ratio == expected_ratio)
        );
    }
}

#[test]
fn simultaneous_modes_cannot_admit_the_squared_equations_closed_opening_branch() {
    let source = source();
    let normal = source.normal_source();
    let ell = normal.indentation_scale_m();
    let a_scale = (ell * normal.effective_radius_m()).sqrt();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let shear_scale = ell * (e / (4.0 * source.reduced_shear_modulus_pa())).sqrt();
    let twist_scale = (a_scale / normal.effective_radius_m())
        * (e / (8.0 * source.torsional_reduced_shear_modulus_pa())).sqrt();
    let shear = 4.1_f64.sqrt() * shear_scale;
    let twist = 12.0_f64.sqrt() * twist_scale;
    let q = [ell, 0.6 * shear, 0.8 * shear, twist];
    source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    source
        .evaluate_torsion_displacement_controlled(q[0], q[3])
        .unwrap();
    assert_eq!(
        source.evaluate_combined_displacement_controlled(q),
        Err(JkrCombinedError::OpeningEdge)
    );
}

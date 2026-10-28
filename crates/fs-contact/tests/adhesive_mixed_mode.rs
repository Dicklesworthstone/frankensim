//! Conservative full-stick Griffith peeling through the original source owners.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrEnvelope, SourceJkrSphereContact,
    mixed_mode::{
        JkrMixedModeControl, JkrMixedModeEnvelope, JkrMixedModeError, SourceJkrMixedMode,
    },
};
use fs_tribo::{InputAuthority, mindlin::full_stick_stiffness_n_per_m};

fn source(scale: f64) -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.05 * scale)
        .normal_state()
        .clone();
    SourceJkrSphereContact::new(
        normal,
        0.01 * scale,
        1e-9 * scale,
        JkrEnvelope {
            maximum_indentation_m: 1e-5 * scale,
            half_space_depth_m: 0.1 * scale,
            layer_thickness_m: 0.1 * scale,
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
fn shear_envelope() -> JkrMixedModeEnvelope {
    JkrMixedModeEnvelope {
        maximum_tangential_displacement_to_radius: 0.1,
        maximum_mean_shear_stress_to_yield: 0.1,
    }
}
fn model(scale: f64) -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(source(scale), shear_envelope()).unwrap()
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}, relative tolerance={tolerance:e}"
    );
}
fn constants(law: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let normal = law.normal_source();
    (
        normal.normal_state().elastic().reduced_modulus_pa(),
        normal.normal_state().adhesion_energy_j_per_m2(),
        normal.effective_radius_m(),
        law.reduced_shear_modulus_pa(),
    )
}

#[test]
fn zero_shear_preserves_original_jkr_source_energy_tangent_and_mindlin_elasticity() {
    let law = model(1.0);
    let normal = law.normal_source();
    let (e, w, r, _) = constants(&law);
    for delta in [
        normal.minimum_attached_indentation_m(),
        normal.minimum_attached_indentation_m().next_up(),
        -0.5e-6,
        0.0,
        2e-6,
        8e-6,
    ] {
        let original = normal.evaluate_attached(delta).unwrap();
        let result = law
            .evaluate_displacement_controlled(delta, [0.0; 2])
            .unwrap();
        assert_eq!(result.contact_radius_m, original.contact_radius_m);
        assert_eq!(
            result.normal_state_identity(),
            original.normal_state_identity()
        );
        assert_eq!(result.model_identity(), law.identity());
        assert_eq!(result.authority, InputAuthority::SyntheticFixture);
        close(
            result.normal_force_n,
            original.normal_force_n,
            PI * w * r,
            3e-13,
        );
        close(
            result.reversible_energy_j,
            original.reversible_energy_j,
            w * result.contact_area_m2 + result.normal_elastic_energy_j,
            2e-13,
        );
        match (result.coupled_stiffness_n_per_m, original.tangent_n_per_m) {
            (Some(matrix), Some(tangent)) => {
                // The immediately adjacent fold state has an ill-conditioned
                // radius subtraction. Both values remain physical and finite.
                let tolerance = if delta == normal.minimum_attached_indentation_m().next_up() {
                    1e-7
                } else {
                    3e-12
                };
                close(
                    matrix[0][0],
                    tangent,
                    tangent.abs().max(e * result.contact_radius_m),
                    tolerance,
                );
                let elastic = normal.normal_state().elastic();
                let bodies = [elastic.surface_a_state(), elastic.surface_b_state()];
                let shear =
                    bodies.map(|s| (0.5 * s.young_modulus_pa()) / (1.0 + s.poisson_ratio()));
                let poisson = bodies.map(|s| s.poisson_ratio());
                let stiffness =
                    full_stick_stiffness_n_per_m(result.contact_radius_m, shear, poisson).unwrap();
                close(matrix[1][1], stiffness, stiffness, 3e-15);
                assert_eq!(matrix[1][2], 0.0);
            }
            (None, None) => assert!(result.marginal),
            other => panic!("inconsistent zero-shear tangent: {other:?}"),
        }
    }
    for n in [normal.force_controlled_pull_off_n(), -1e-3, 0.0, 5e-3, 1e-2] {
        let result = law.evaluate_force_controlled(n, [0.0; 2]).unwrap();
        let original = normal.evaluate_attached(result.indentation_m).unwrap();
        close(
            result.contact_radius_m,
            original.contact_radius_m,
            result.contact_radius_m,
            4e-13,
        );
        close(result.normal_force_n, n, PI * w * r, 3e-13);
        close(
            result.reversible_energy_j,
            original.reversible_energy_j,
            w * result.contact_area_m2 + result.normal_elastic_energy_j,
            4e-13,
        );
    }
}

#[test]
fn force_control_matches_independent_griffith_peeling_and_real_area_reduction() {
    let law = model(1.0);
    let (e, w, r, g) = constants(&law);
    let alpha = e / (4.0 * g);
    for n in [-1e-3, 0.0, 5e-3] {
        let b = 3.0 * PI * w * r;
        let critical = ((b * b + 2.0 * b * n) / alpha).sqrt();
        close(
            law.critical_tangential_force_n(n).unwrap(),
            critical,
            critical,
            3e-15,
        );
        let mut previous_area = f64::INFINITY;
        for fraction in [0.0, 0.3, 0.7, 0.99, 1.0] {
            let t = fraction * law.critical_tangential_force_n(n).unwrap();
            let response = law.evaluate_force_controlled(n, [t, 0.0]).unwrap();
            let a = response.contact_radius_m;
            let hertz = 4.0 * e * a.powi(3) / (3.0 * r);
            close(
                (hertz - n).powi(2) + alpha * t * t,
                8.0 * PI * w * e * a.powi(3),
                b * b + 2.0 * b * n,
                8e-15,
            );
            close(
                response.normal_energy_release_rate_j_m2 + response.shear_energy_release_rate_j_m2,
                w,
                w,
                3e-14,
            );
            close(response.normal_force_n, n, b, 3e-15);
            close(response.tangential_force_n[0], t, critical, 3e-15);
            assert!(response.contact_area_m2 < previous_area);
            previous_area = response.contact_area_m2;
            let friction = response.friction_contact_state().unwrap();
            close(
                friction.area_m2(),
                response.contact_area_m2,
                response.contact_area_m2,
                0.0,
            );
            close(friction.repulsive_normal_force_n(), hertz, hertz, 3e-15);
            assert_eq!(response.marginal, fraction == 1.0);
            assert_eq!(
                response.coupled_compliance_m_per_n.is_none(),
                fraction == 1.0
            );
        }
        assert!(matches!(
            law.evaluate_force_controlled(
                n,
                [law.critical_tangential_force_n(n).unwrap().next_up(), 0.0]
            ),
            Err(JkrMixedModeError::NoFullStickEquilibrium {
                control: JkrMixedModeControl::Force,
                ..
            })
        ));
    }
    let minimum = law.normal_source().force_controlled_pull_off_n();
    assert_eq!(law.critical_tangential_force_n(minimum).unwrap(), 0.0);
    assert!(matches!(
        law.evaluate_force_controlled(minimum.next_down(), [0.0; 2]),
        Err(JkrMixedModeError::BelowNormalPullOff { .. })
    ));
}

#[test]
fn displacement_fold_resolves_the_distinct_peeling_limit_and_adjacent_states() {
    let law = model(1.0);
    let (e, w, r, g) = constants(&law);
    for delta in [
        law.normal_source()
            .minimum_attached_indentation_m()
            .next_up(),
        0.0,
        2e-6,
        8e-6,
    ] {
        let critical = law.critical_tangential_displacement_m(delta).unwrap();
        let terminal = law
            .evaluate_displacement_controlled(delta, [critical, 0.0])
            .unwrap();
        let a = terminal.contact_radius_m;
        close(
            a.powi(3) - r * delta * a,
            PI * w * r * r / (2.0 * e),
            a.powi(3),
            4e-15,
        );
        close(
            e * (delta - a * a / r).powi(2) + 4.0 * g * critical * critical,
            2.0 * PI * w * a,
            2.0 * PI * w * a,
            6e-15,
        );
        assert!(terminal.marginal);
        assert_eq!(terminal.coupled_stiffness_n_per_m, None);
        assert_eq!(terminal.coupled_compliance_m_per_n, None);
        assert!(terminal.force_radius_curvature_n_per_m < 0.0);
        let adjacent = law
            .evaluate_displacement_controlled(delta, [critical.next_down(), 0.0])
            .unwrap();
        assert!(!adjacent.marginal);
        assert!(adjacent.displacement_radius_curvature_n_per_m > 0.0);
        assert!(adjacent.coupled_stiffness_n_per_m.is_some());
        assert!(matches!(
            law.evaluate_displacement_controlled(delta, [critical.next_up(), 0.0]),
            Err(JkrMixedModeError::NoFullStickEquilibrium {
                control: JkrMixedModeControl::Displacement,
                ..
            })
        ));
        let pure = law
            .evaluate_displacement_controlled(delta, [0.0; 2])
            .unwrap();
        assert!(terminal.contact_area_m2 <= pure.contact_area_m2);
    }
    assert_eq!(
        law.critical_tangential_displacement_m(
            law.normal_source().minimum_attached_indentation_m()
        )
        .unwrap(),
        0.0
    );
}

#[test]
fn one_potential_gives_both_forces_full_reciprocal_tangent_and_inverse_compliance() {
    let law = model(1.0);
    let (e, _, _, _) = constants(&law);
    let q = [2e-6, 0.7e-6, -0.4e-6];
    let center = law
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    let efforts = [
        center.normal_force_n,
        center.tangential_force_n[0],
        center.tangential_force_n[1],
    ];
    let matrix = center.coupled_stiffness_n_per_m.unwrap();
    let inverse = center.coupled_compliance_m_per_n.unwrap();
    let h = 1e-11;
    let stiffness_scale = e * center.contact_radius_m;
    for j in 0..3 {
        let mut lower = q;
        lower[j] -= h;
        let mut upper = q;
        upper[j] += h;
        let left = law
            .evaluate_displacement_controlled(lower[0], [lower[1], lower[2]])
            .unwrap();
        let right = law
            .evaluate_displacement_controlled(upper[0], [upper[1], upper[2]])
            .unwrap();
        close(
            (right.reversible_energy_j - left.reversible_energy_j) / (2.0 * h),
            efforts[j],
            efforts[j].abs(),
            2e-9,
        );
        let low = [
            left.normal_force_n,
            left.tangential_force_n[0],
            left.tangential_force_n[1],
        ];
        let high = [
            right.normal_force_n,
            right.tangential_force_n[0],
            right.tangential_force_n[1],
        ];
        for i in 0..3 {
            close(
                (high[i] - low[i]) / (2.0 * h),
                matrix[i][j],
                stiffness_scale,
                2e-9,
            );
            close(matrix[i][j], matrix[j][i], stiffness_scale, 3e-15);
            let product: f64 = (0..3).map(|k| matrix[i][k] * inverse[k][j]).sum();
            close(product, if i == j { 1.0 } else { 0.0 }, 1.0, 3e-15);
        }
    }
    let force_center = law
        .evaluate_force_controlled(efforts[0], [efforts[1], efforts[2]])
        .unwrap();
    close(
        force_center.contact_radius_m,
        center.contact_radius_m,
        center.contact_radius_m,
        3e-15,
    );
    let force_step = 1e-8;
    for j in 0..3 {
        let mut lower = efforts;
        lower[j] -= force_step;
        let mut upper = efforts;
        upper[j] += force_step;
        let left = law
            .evaluate_force_controlled(lower[0], [lower[1], lower[2]])
            .unwrap();
        let right = law
            .evaluate_force_controlled(upper[0], [upper[1], upper[2]])
            .unwrap();
        let low = [
            left.indentation_m,
            left.tangential_displacement_m[0],
            left.tangential_displacement_m[1],
        ];
        let high = [
            right.indentation_m,
            right.tangential_displacement_m[0],
            right.tangential_displacement_m[1],
        ];
        for i in 0..3 {
            close(
                (high[i] - low[i]) / (2.0 * force_step),
                inverse[i][j],
                1.0 / stiffness_scale,
                2e-9,
            );
        }
    }
    assert!(matrix[0][1] != 0.0 && matrix[0][2] != 0.0 && matrix[1][2] != 0.0);
}

#[test]
fn independently_integrated_normal_and_vector_shear_work_close_the_conservative_potential() {
    let law = model(1.0);
    // A changing normal approach and two rotating shear components prevent a
    // scalar fixed-area law from accidentally satisfying this path integral.
    let path = |t: f64| {
        (
            [
                2e-6 + 0.3e-6 * t + 0.1e-6 * (2.0 * PI * t).sin(),
                0.5e-6 * (2.0 * PI * t).cos(),
                0.4e-6 * (2.0 * PI * t).sin(),
            ],
            [
                0.3e-6 + 0.2e-6 * PI * (2.0 * PI * t).cos(),
                -1e-6 * PI * (2.0 * PI * t).sin(),
                0.8e-6 * PI * (2.0 * PI * t).cos(),
            ],
        )
    };
    let integrand = |t: f64| {
        let (q, dq) = path(t);
        let response = law
            .evaluate_displacement_controlled(q[0], [q[1], q[2]])
            .unwrap();
        response.normal_force_n * dq[0]
            + response.tangential_force_n[0] * dq[1]
            + response.tangential_force_n[1] * dq[2]
    };
    let n = 1024;
    let mut work = integrand(0.0) + integrand(1.0);
    for i in 1..n {
        work += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(f64::from(i) / f64::from(n));
    }
    work /= 3.0 * f64::from(n);
    let (q0, _) = path(0.0);
    let (q1, _) = path(1.0);
    let start = law
        .evaluate_displacement_controlled(q0[0], [q0[1], q0[2]])
        .unwrap();
    let end = law
        .evaluate_displacement_controlled(q1[0], [q1[1], q1[2]])
        .unwrap();
    close(
        work,
        end.reversible_energy_j - start.reversible_energy_j,
        work.abs(),
        2e-9,
    );
    // A reversed query retraces precisely; no state mutation or loss is hidden.
    assert_eq!(
        law.evaluate_displacement_controlled(q0[0], [q0[1], q0[2]])
            .unwrap(),
        start
    );
}

#[test]
fn isotropic_rotation_and_source_length_similarity_preserve_the_coupled_physics() {
    let law = model(1.0);
    let center = law
        .evaluate_displacement_controlled(2e-6, [0.6e-6, -0.8e-6])
        .unwrap();
    let rotated = law
        .evaluate_displacement_controlled(2e-6, [0.8e-6, 0.6e-6])
        .unwrap();
    assert_eq!(center.contact_area_m2, rotated.contact_area_m2);
    assert_eq!(center.reversible_energy_j, rotated.reversible_energy_j);
    assert_eq!(
        rotated.tangential_force_n,
        [-center.tangential_force_n[1], center.tangential_force_n[0]]
    );
    for scale in [1e-4, 1e-2, 1e2, 1e4] {
        let scaled = model(scale);
        let response = scaled
            .evaluate_displacement_controlled(2e-6 * scale, [0.6e-6 * scale, -0.8e-6 * scale])
            .unwrap();
        close(
            response.contact_radius_m,
            center.contact_radius_m * scale,
            center.contact_radius_m * scale,
            8e-15,
        );
        close(
            response.contact_area_m2,
            center.contact_area_m2 * scale * scale,
            center.contact_area_m2 * scale * scale,
            1e-14,
        );
        close(
            response.normal_force_n,
            center.normal_force_n * scale * scale,
            center.normal_force_n.abs() * scale * scale,
            2e-14,
        );
        close(
            response.reversible_energy_j,
            center.reversible_energy_j * scale.powi(3),
            center.reversible_energy_j.abs() * scale.powi(3),
            3e-14,
        );
        let matrix = response.coupled_stiffness_n_per_m.unwrap();
        let original = center.coupled_stiffness_n_per_m.unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(
                    matrix[i][j],
                    original[i][j] * scale,
                    original[i][j].abs() * scale,
                    2e-14,
                );
            }
        }
        assert_ne!(
            scaled.normal_source().normal_state().identity(),
            law.normal_source().normal_state().identity()
        );
    }
}

#[test]
fn source_and_explicit_shear_limits_refuse_without_changing_the_original_owner() {
    let normal = source(1.0);
    let original = normal.clone();
    let mut envelope = shear_envelope();
    envelope.maximum_tangential_displacement_to_radius = 1e-4;
    let law = SourceJkrMixedMode::new(normal, envelope).unwrap();
    assert_eq!(law.normal_source(), &original);
    assert!(matches!(
        law.evaluate_displacement_controlled(2e-6, [1e-6, 0.0]),
        Err(JkrMixedModeError::OutsideApplicability { .. })
    ));
    assert!(matches!(
        law.evaluate_displacement_controlled(1e-4, [0.0; 2]),
        Err(JkrMixedModeError::Normal(_))
    ));
    for vector in [[f64::NAN, 0.0], [0.0, f64::INFINITY]] {
        assert!(matches!(
            law.evaluate_displacement_controlled(2e-6, vector),
            Err(JkrMixedModeError::Invalid { .. })
        ));
    }
    assert_eq!(law.normal_source(), &original);
    assert_ne!(law.identity(), model(1.0).identity());
}

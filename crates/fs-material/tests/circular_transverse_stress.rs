//! G1/G3: Saint-Venant equations, free contours, force/work and source admission.
//! Synthetic coefficients establish the implementation, not component life.

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::transverse_stress::CircularTransverseShearRecovery;
use fs_material::section::transverse_stress::layered::LayeredCircularTransverseShearRecovery;
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
    LayerIsotropicParameters,
};
use fs_material::state_point::{
    IsotropicElasticParameters, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Pressure};

const G: f64 = 2e8;

#[test]
fn bonded_layered_sources_keep_traction_strain_and_independent_two_layer_center_stress() {
    let original = section(0.0, 0.25, 0.8, true, true, true);
    assert!(CircularTransverseShearRecovery::try_new(Arc::clone(&original)).is_err());
    let field = LayeredCircularTransverseShearRecovery::try_new(Arc::clone(&original)).unwrap();
    assert_eq!(field.section(), original.as_ref());
    // Independent three-equation core/sleeve transmission solution, G2=2G1.
    let t = 0.6_f64;
    let c = 3.5 / 4.0;
    let b_sleeve = -c * (1.0 - t * t) * t * t / (3.0 - t * t);
    let a_core = b_sleeve - c + b_sleeve / (t * t);
    let center_factor = -2.0 * a_core / (1.25 * (2.0 - t.powi(4)));
    let center = field.evaluate(0, [0.0; 2], [2.0, -1.0], 0.01).unwrap();
    for (stress, q) in center.shear_stress_pa.into_iter().zip([2.0, -1.0]) {
        near(stress, center_factor * q / core::f64::consts::PI, 1.0);
    }
    for point in [[t, 0.0], [0.0, t]] {
        let core = field.evaluate(0, point, [2.0, -1.0], 0.01).unwrap();
        let sleeve = field.evaluate(1, point, [2.0, -1.0], 0.01).unwrap();
        let dot = |values: [f64; 2]| values[0] * point[0] + values[1] * point[1];
        let tangent = |values: [f64; 2]| -values[0] * point[1] + values[1] * point[0];
        near(dot(core.shear_stress_pa), dot(sleeve.shear_stress_pa), 1.0);
        near(
            tangent(core.engineering_shear_strain),
            tangent(sleeve.engineering_shear_strain),
            1e-8,
        );
        assert!(
            core.maximum_layer_engineering_shear
                >= core.engineering_shear_strain[0].hypot(core.engineering_shear_strain[1])
        );
    }
    let corrected =
        LayeredCircularTransverseShearRecovery::try_new(section(0.0, 0.25, 0.3, true, true, true))
            .unwrap();
    assert_eq!(
        center,
        corrected.evaluate(0, [0.0; 2], [2.0, -1.0], 0.01).unwrap()
    );
    // Querying a contour with zero traction cannot conceal an overstressed layer.
    assert!(field.evaluate(1, [1.0, 0.0], [2.0, 0.0], 1e-12).is_err());
    assert!(field.evaluate(2, [0.0; 2], [0.0; 2], 0.01).is_err());
    assert!(field.evaluate(0, [0.8, 0.0], [0.0; 2], 0.01).is_err());
    assert!(
        LayeredCircularTransverseShearRecovery::try_new(section(0.0, 0.25, 0.8, true, true, false))
            .is_err()
    );
}

#[test]
fn layered_homogeneous_limit_and_common_poisson_source_boundary_are_explicit() {
    for bore in [0.0, 0.3] {
        for nu in [-0.75, 0.0, 0.25, 0.49] {
            let source = section(bore, nu, 0.8, true, false, true);
            let homogeneous =
                CircularTransverseShearRecovery::try_new(Arc::clone(&source)).unwrap();
            let layered = LayeredCircularTransverseShearRecovery::try_new(source).unwrap();
            for (i, p) in [(0, [(bore + 0.6) / 2.0, 0.0]), (1, [0.0, 0.8])] {
                let a = homogeneous.evaluate(i, p, [1.0, -0.3], 0.01).unwrap();
                let b = layered.evaluate(i, p, [1.0, -0.3], 0.01).unwrap();
                for (x, y) in a.shear_stress_pa.into_iter().zip(b.shear_stress_pa) {
                    near(x, y, 1.0);
                }
                near(
                    a.maximum_layer_engineering_shear,
                    b.maximum_layer_engineering_shear,
                    1e-8,
                );
            }
        }
    }
    let a = section(0.0, 0.25, 0.8, true, false, true);
    let b = section(0.0, 0.3, 0.8, true, true, true);
    let layers = vec![
        a.elastic().layers()[0].clone(),
        b.elastic().layers()[1].clone(),
    ];
    let bindings = layers
        .iter()
        .map(|layer| {
            let role = |name| layer.material.property(name).unwrap().requirement().clone();
            LayerIsotropicParameters {
                density: role("density"),
                parameters: IsotropicElasticParameters::YoungShear {
                    young_modulus: role("young_modulus"),
                    shear_modulus: role("shear_modulus"),
                },
            }
        })
        .collect::<Vec<_>>();
    let mixed = ConcentricElasticSection::try_new_with_isotropic_parameters(0.0, layers, &bindings)
        .unwrap()
        .with_isotropic_shear(0.8)
        .unwrap();
    assert!(LayeredCircularTransverseShearRecovery::try_new(Arc::new(mixed)).is_err());
}

fn section(
    bore: f64,
    nu: f64,
    correction: f64,
    split: bool,
    heterogeneous: bool,
    isotropic: bool,
) -> Arc<ConcentricShearSection> {
    let radii = if split { vec![0.6, 1.0] } else { vec![1.0] };
    let mut layers = Vec::new();
    let mut bindings = Vec::new();
    for (index, radius) in radii.into_iter().enumerate() {
        let g = if heterogeneous && index == 1 {
            2.0 * G
        } else {
            G
        };
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [
            ("density", 1000.0 + index as f64, Density::DIMS),
            ("young_modulus", 2.0 * g * (1.0 + nu), Pressure::DIMS),
            ("shear_modulus", g, Pressure::DIMS),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic circular flexure".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
            requirements.push(
                ScalarPropertyRequirement::try_new(
                    name,
                    dims,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
            );
        }
        // Each layer is its own material state at revision 0; a revision is
        // a superseding version of one state, not a layer index.
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic circular material".into(),
                phase: "solid".into(),
                process: format!("synthetic shear test layer {index}"),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        layers.push(ConcentricMaterialLayer {
            outer_radius_m: radius,
            material: resolve_material_state_point(
                &card,
                &QueryPoint::new(),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        });
        bindings.push(LayerIsotropicParameters {
            density: requirements[0].clone(),
            parameters: IsotropicElasticParameters::YoungShear {
                young_modulus: requirements[1].clone(),
                shear_modulus: requirements[2].clone(),
            },
        });
    }
    let elastic = if isotropic {
        ConcentricElasticSection::try_new_with_isotropic_parameters(bore, layers, &bindings)
    } else {
        ConcentricElasticSection::try_new(bore, layers)
    }
    .unwrap();
    Arc::new(
        elastic
            .with_shear_modulus("shear_modulus", correction)
            .unwrap(),
    )
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// Independent Cartesian/Laurent form; production evaluates radial/tangential
// eigenfields with factored free-contour roots instead.
fn reference(point: [f64; 2], q: [f64; 2], bore: f64, nu: f64) -> [f64; 2] {
    let [y, z] = point;
    let r2 = y * y + z * z;
    let inertia = core::f64::consts::PI * (1.0 - bore.powi(4)) / 4.0;
    let denominator = 8.0 * inertia * (1.0 + nu);
    let a = 3.0 + 2.0 * nu;
    if r2 == 0.0 {
        return q.map(|v| a * v / denominator);
    }
    let dipole = a * bore * bore / r2;
    let isotropic = (a * (1.0 + bore * bore) - (1.0 - 2.0 * nu) * r2 + dipole) / denominator;
    let directional = -2.0 * (1.0 + 2.0 * nu + dipole / r2) / denominator;
    let dot = y * q[0] + z * q[1];
    [
        isotropic * q[0] + directional * dot * y,
        isotropic * q[1] + directional * dot * z,
    ]
}

#[test]
fn exact_circle_annulus_fields_reversal_contours_and_poisson_effect() {
    for bore in [0.0, 0.35, 0.999] {
        for nu in [-0.75, 0.0, 0.25, 0.49] {
            let recovery = CircularTransverseShearRecovery::try_new(section(
                bore, nu, 0.8, false, false, true,
            ))
            .unwrap();
            let q = [2.0, -1.0];
            for point in [
                [1.0, 0.0],
                [0.0, 1.0],
                [(1.0 + bore) / 2.0, 0.0],
                [0.0, -(1.0 + bore) / 2.0],
            ] {
                let sample = recovery.evaluate(0, point, q, 0.1).unwrap();
                let expected = reference(point, q, bore, nu);
                for i in 0..2 {
                    near(sample.shear_stress_pa[i], expected[i], 1.0 / (1.0 - bore));
                }
                let reversed = recovery.evaluate(0, point, q.map(|v| -v), 0.1).unwrap();
                assert_eq!(reversed.shear_stress_pa, sample.shear_stress_pa.map(|v| -v));
                assert_eq!(
                    reversed.stored_energy_density_j_m3,
                    sample.stored_energy_density_j_m3
                );
                let energy = sample
                    .shear_stress_pa
                    .iter()
                    .map(|t| t * t / (2.0 * G))
                    .sum::<f64>();
                near(sample.stored_energy_density_j_m3, energy, energy.max(1e-30));
            }
            for radius in [bore, 1.0] {
                if radius == 0.0 {
                    continue;
                }
                for point in [[radius, 0.0], [0.0, radius], [-radius, 0.0], [0.0, -radius]] {
                    let stress = recovery.evaluate(0, point, q, 0.1).unwrap().shear_stress_pa;
                    assert_eq!(stress[0] * point[0] + stress[1] * point[1], 0.0);
                }
            }
        }
    }
}

#[test]
fn independent_divergence_compatibility_and_rotating_force_field() {
    let nu = 0.25;
    let bore = 0.35;
    let recovery =
        CircularTransverseShearRecovery::try_new(section(bore, nu, 0.2, false, false, true))
            .unwrap();
    let q = [1.7, -0.6];
    let point = [0.51, 0.32];
    let h = 1e-5;
    let stress = |p| recovery.evaluate(0, p, q, 0.1).unwrap().shear_stress_pa;
    let center = stress(point);
    let plus_y = stress([point[0] + h, point[1]]);
    let minus_y = stress([point[0] - h, point[1]]);
    let plus_z = stress([point[0], point[1] + h]);
    let minus_z = stress([point[0], point[1] - h]);
    let div = (plus_y[0] - minus_y[0] + plus_z[1] - minus_z[1]) / (2.0 * h);
    let inertia = core::f64::consts::PI * (1.0 - bore.powi(4)) / 4.0;
    let expected_div = -(q[0] * point[0] + q[1] * point[1]) / inertia;
    assert!((div - expected_div).abs() < 2e-8);
    for i in 0..2 {
        let laplacian =
            (plus_y[i] + minus_y[i] + plus_z[i] + minus_z[i] - 4.0 * center[i]) / (h * h);
        assert!((laplacian + q[i] / ((1.0 + nu) * inertia)).abs() < 2e-5);
    }
    let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    let rotated = recovery.evaluate(0, rotate(point), rotate(q), 0.1).unwrap();
    for i in 0..2 {
        near(rotated.shear_stress_pa[i], rotate(center)[i], 1.0);
    }
}

#[test]
fn force_torque_and_solid_circle_energy_integrate_independently() {
    const GAUSS: [(f64, f64); 4] = [
        (-0.8611363115940526, 0.3478548451374538),
        (-0.3399810435848563, 0.6521451548625461),
        (0.3399810435848563, 0.6521451548625461),
        (0.8611363115940526, 0.3478548451374538),
    ];
    for bore in [0.0, 0.35] {
        let nu = 0.25;
        let recovery =
            CircularTransverseShearRecovery::try_new(section(bore, nu, 0.8, false, false, true))
                .unwrap();
        let q = [1.7, -0.6];
        let mut force = [0.0; 2];
        let mut torque = 0.0;
        let mut energy = 0.0;
        for radial in 0..32 {
            let half = (1.0 - bore) / 64.0;
            let mid = bore + (2 * radial + 1) as f64 * half;
            for (x, weight) in GAUSS {
                let r = mid + x * half;
                for angle in 0..32 {
                    let theta = core::f64::consts::TAU * (angle as f64 + 0.5) / 32.0;
                    let point = [r * theta.cos(), r * theta.sin()];
                    let area = weight * half * r * core::f64::consts::TAU / 32.0;
                    let sample = recovery.evaluate(0, point, q, 0.1).unwrap();
                    for i in 0..2 {
                        force[i] += sample.shear_stress_pa[i] * area;
                    }
                    torque += (point[0] * sample.shear_stress_pa[1]
                        - point[1] * sample.shear_stress_pa[0])
                        * area;
                    energy += sample.stored_energy_density_j_m3 * area;
                }
            }
        }
        for i in 0..2 {
            near(force[i], q[i], 1.0);
        }
        near(torque, 0.0, 1.0);
        if bore == 0.0 {
            let factor = (7.0 + 14.0 * nu + 8.0 * nu * nu) / (6.0 * (1.0 + nu).powi(2));
            let expected = factor * (q[0] * q[0] + q[1] * q[1]) / (2.0 * G * core::f64::consts::PI);
            near(energy, expected, expected);
        }
    }
}

#[test]
fn original_layer_subdivision_and_authored_correction_do_not_fabricate_stress() {
    let whole =
        CircularTransverseShearRecovery::try_new(section(0.1, 0.25, 0.2, false, false, true))
            .unwrap();
    let split =
        CircularTransverseShearRecovery::try_new(section(0.1, 0.25, 0.9, true, false, true))
            .unwrap();
    assert_ne!(
        whole.section().shear_stiffness_n(),
        split.section().shear_stiffness_n()
    );
    for (layer, p) in [(0, [0.4, 0.2]), (1, [0.0, 0.8])] {
        let a = whole.evaluate(0, p, [1.0, 2.0], 0.1).unwrap();
        let b = split.evaluate(layer, p, [1.0, 2.0], 0.1).unwrap();
        for i in 0..2 {
            near(a.shear_stress_pa[i], b.shear_stress_pa[i], 1.0);
        }
        near(
            a.stored_energy_density_j_m3,
            b.stored_energy_density_j_m3,
            a.stored_energy_density_j_m3,
        );
    }
    assert!(
        CircularTransverseShearRecovery::try_new(section(0.0, 0.25, 0.8, true, true, true))
            .is_err()
    );
    assert!(
        CircularTransverseShearRecovery::try_new(section(0.0, 0.25, 0.8, false, false, false))
            .is_err()
    );
}

#[test]
fn zero_traction_does_not_hide_whole_layer_overstrain_or_invalid_points() {
    let recovery =
        CircularTransverseShearRecovery::try_new(section(0.35, -0.75, 0.8, false, false, true))
            .unwrap();
    let q = [1.0, 0.0];
    let sample = recovery.evaluate(0, [1.0, 0.0], q, 0.1).unwrap();
    assert_eq!(sample.shear_stress_pa, [0.0; 2]);
    assert!(sample.maximum_layer_engineering_shear > 0.0);
    assert!(
        recovery
            .evaluate(
                0,
                [1.0, 0.0],
                q,
                sample.maximum_layer_engineering_shear * 0.99
            )
            .is_err()
    );
    for point in [[0.0; 2], [1.01, 0.0], [f64::NAN, 0.0]] {
        assert!(recovery.evaluate(0, point, q, 0.1).is_err());
    }
    assert!(recovery.evaluate(1, [0.5, 0.0], q, 0.1).is_err());
    assert!(
        recovery
            .evaluate(0, [0.5, 0.0], [f64::INFINITY, 0.0], 0.1)
            .is_err()
    );
}

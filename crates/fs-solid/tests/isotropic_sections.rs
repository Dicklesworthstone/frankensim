//! G1/G3: the same independent material parameters drive frames and spatial rods.
//! Synthetic data check the constitutive change of variables, not a real alloy.

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, LayerIsotropicParameters,
};
use fs_material::state_point::{
    IsotropicElasticParameters, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Pressure};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};
use fs_solid::{Frame2d, Geometry};

fn section(
    bulk: f64,
    shear: f64,
    derived: bool,
    restricted_bulk: bool,
) -> ConcentricElasticSection {
    let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
    let entries = if derived {
        vec![
            ("lot_density", 2400.0, Density::DIMS),
            ("compression_modulus", bulk, Pressure::DIMS),
            ("torsion_modulus", shear, Pressure::DIMS),
        ]
    } else {
        vec![
            ("density", 2400.0, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            ("shear_modulus", shear, Pressure::DIMS),
        ]
    };
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in entries {
        let validity = if restricted_bulk && name == "compression_modulus" {
            ValidityDomain::unconstrained().with("omega", 0.0, 100.0)
        } else {
            ValidityDomain::unconstrained()
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic independent elastic parameters".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "reference specimen".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = if restricted_bulk {
        QueryPoint::new().with("omega", 0.0).unwrap()
    } else {
        QueryPoint::new()
    };
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let layers = vec![ConcentricMaterialLayer {
        outer_radius_m: 0.02,
        material,
    }];
    if derived {
        ConcentricElasticSection::try_new_with_isotropic_parameters(
            0.0,
            layers,
            &[LayerIsotropicParameters {
                density: requirements[0].clone(),
                parameters: IsotropicElasticParameters::BulkShear {
                    bulk_modulus: requirements[1].clone(),
                    shear_modulus: requirements[2].clone(),
                },
            }],
        )
        .unwrap()
    } else {
        ConcentricElasticSection::try_new(0.0, layers).unwrap()
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

fn frame(section: ConcentricElasticSection) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.add_material_element(0, 1, 5, section).unwrap();
    frame
}

#[test]
fn independent_moduli_and_aliased_density_reach_weight_compliance_and_frame_modes() {
    let mut reference = None;
    for derived in [false, true] {
        let source = section(2.0e6, 0.8e6, derived, false);
        if derived {
            assert!(source.layers()[0].material.property("density").is_none());
            assert!(
                source.layers()[0]
                    .material
                    .property("young_modulus")
                    .is_none()
            );
        }
        let young = 9.0 * 2.0e6 * 0.8e6 / (3.0 * 2.0e6 + 0.8e6);
        let area = std::f64::consts::PI * 0.02_f64.powi(2);
        let inertia = std::f64::consts::PI * 0.02_f64.powi(4) / 4.0;
        let mut model = frame(source.clone());
        let modes = model.natural_frequencies().unwrap();
        let weight = model
            .inertial_forces(&[0.0, -9.81, 0.0, 0.0, -9.81, 0.0])
            .unwrap();
        near(weight[1] + weight[4], -2400.0 * area * 9.81, 1e-12);
        model
            .static_load(&[0.0, 0.0, 0.0, 0.02, 0.01, 0.0], 1)
            .unwrap();
        near(model.displacements()[3], 0.02 / (young * area), 1e-10);
        near(
            model.displacements()[4],
            0.01 / (3.0 * young * inertia),
            1e-10,
        );
        assert_eq!(model.material_section(0), Some(&source));
        if let Some(original) = &reference {
            for (a, b) in modes.iter().zip(original) {
                near(*a, *b, 1e-10);
            }
        } else {
            reference = Some(modes);
        }
    }
}

#[test]
fn independent_moduli_drive_rod_forces_torque_links_and_reference_inertia() {
    let gate = CancelGate::new_clock_free();
    let mut baseline = None;
    for derived in [false, true] {
        let elastic = section(2.0e6, 0.8e6, derived, false);
        let shear = if derived {
            elastic.clone().with_isotropic_shear(0.8).unwrap()
        } else {
            elastic
                .clone()
                .with_shear_modulus("shear_modulus", 0.8)
                .unwrap()
        };
        let link = MaterialTensionLink::try_new(
            elastic,
            [
                TensionLinkAttachment {
                    node: 0,
                    offset_body_m: [0.0; 3],
                },
                TensionLinkAttachment {
                    node: 1,
                    offset_body_m: [0.0; 3],
                },
            ],
            1.0,
            2,
            0.1,
        )
        .unwrap();
        let rod = SpatialMaterialRod::try_straight(1.0, vec![shear])
            .unwrap()
            .with_tension_links(vec![link], &gate)
            .unwrap();
        let model = rod.into_motion_model(RodRootBoundary::Free, &gate).unwrap();
        let mut pose = model.initial_pose();
        pose.positions_m[1][0] += 0.003;
        pose.positions_m[1][1] += 0.002;
        pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(0.01, 0.02, -0.015))).unwrap();
        let response = model.evaluate(&pose, &gate).unwrap();
        let mass = model.inertia().iter().map(|i| i.mass_kg()).sum::<f64>();
        near(
            mass,
            3.0 * 2400.0 * std::f64::consts::PI * 0.02_f64.powi(2),
            1e-12,
        );
        let values = std::iter::once(response.potential_energy_j)
            .chain(response.force_n.iter().flatten().copied())
            .chain(response.torque_body_nm.iter().flatten().copied())
            .chain(model.inertia().iter().flat_map(|i| i.principal_kg_m2()))
            .collect::<Vec<_>>();
        if let Some(original) = &baseline {
            for (a, b) in values.iter().zip(original) {
                near(*a, *b, 1e-10);
            }
        } else {
            baseline = Some(values);
        }
    }
}

#[test]
fn original_bulk_modulus_support_is_not_lost_by_derived_young_and_shear() {
    let elastic = section(2.0e6, 0.8e6, true, true);
    let mut model = Frame2d::new(Geometry::Linear);
    model.add_node(0.0, 0.0);
    model.add_node(1.0, 0.0);
    assert!(
        model
            .add_material_element(0, 1, 3, elastic.clone())
            .is_err()
    );
    let rod =
        SpatialMaterialRod::try_straight(1.0, vec![elastic.with_isotropic_shear(0.8).unwrap()])
            .unwrap();
    assert!(
        rod.into_motion_model(RodRootBoundary::Free, &CancelGate::new_clock_free())
            .is_err()
    );
}

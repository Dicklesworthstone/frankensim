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
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Pressure, QuantitySpec};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::rod::{MaterialRod, TipLoad};
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
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
    let mut reference: Option<Vec<f64>> = None;
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
    let mut baseline: Option<Vec<f64>> = None;
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


// These are separate stress-free specimens at a held environment, not a live
// update of geometry, mass or constitutive history.
fn environment_section(
    axis: &str,
    quantity: Option<QuantitySpec>,
    bounds: [f64; 2],
    value: f64,
    extra: Option<(&str, [f64; 2], f64)>,
) -> ConcentricElasticSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, values) in [
        ("density", Density::DIMS, [2400.0, 3600.0]),
        ("young_modulus", Pressure::DIMS, [2e6, 1e6]),
        ("shear_modulus", Pressure::DIMS, [0.8e6, 0.4e6]),
    ] {
        let mut validity = match quantity {
            Some(q) => ValidityDomain::unconstrained().with_quantity(axis, q, bounds[0], bounds[1]),
            None => ValidityDomain::unconstrained().with(axis, bounds[0], bounds[1]),
        };
        if name == "young_modulus" {
            if let Some((name, limits, _)) = extra {
                validity = validity.with(name, limits[0], limits[1]);
            }
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: quantity.map_or(fs_qty::Dims::NONE, |q| q.dims()),
                knots: vec![(bounds[0], values[0]), (bounds[1], values[1])],
                dims,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured held-environment section".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured solid".into(),
        phase: "solid".into(),
        process: "held environment".into(),
        revision: 0,
    }, claims, vec![]).unwrap();
    let mut point = match quantity {
        Some(q) => QueryPoint::new().with_quantity(axis, q, value).unwrap(),
        None => QueryPoint::new().with(axis, value).unwrap(),
    };
    if let Some((name, _, value)) = extra {
        point = point.with(name, value).unwrap();
    }
    let material = resolve_material_state_point(
        &card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 0.02,
        material,
    }]).unwrap()
}

fn environment_link(
    source: ConcentricElasticSection,
) -> Result<MaterialTensionLink, SpatialRodError> {
    MaterialTensionLink::try_new(source, [
        TensionLinkAttachment { node: 0, offset_body_m: [0.0; 3] },
        TensionLinkAttachment { node: 1, offset_body_m: [0.0; 3] },
    ], 1.0, 2, 0.1)
}

#[test]
fn g1_held_environment_reaches_native_frame_rod_shear_and_linked_inertia() {
    let gate = CancelGate::new_clock_free();
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static,
    ));
    let pressure = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Pressure, ValueForm::Static,
    ));
    for (axis, quantity, bounds) in [
        ("water_content", moisture, [0.0, 0.2]),
        ("ambient_pressure", pressure, [1e5, 2e5]),
    ] {
        let sources = bounds.map(|value| {
            environment_section(axis, Some(quantity), bounds, value, None)
        });
        let before = &sources[0].layers()[0].material;
        let after = &sources[1].layers()[0].material;
        assert_eq!(before.card_identity(), after.card_identity());
        assert_ne!(before.identity(), after.identity());
        let mut frequencies = Vec::new();
        for (source, [density, young]) in sources.into_iter().zip([
            [2400.0, 2e6], [3600.0, 1e6],
        ]) {
            let area = std::f64::consts::PI * 0.02_f64.powi(2);
            let moment = std::f64::consts::PI * 0.02_f64.powi(4) / 4.0;
            let mut elastic_frame = frame(source.clone());
            frequencies.push(elastic_frame.natural_frequencies().unwrap()[0]);
            let weight = elastic_frame.inertial_forces(
                &[0.0, -9.81, 0.0, 0.0, -9.81, 0.0],
            ).unwrap();
            near(weight[1] + weight[4], -density * area * 9.81, 1e-12);
            elastic_frame.static_load(&[0.0, 0.0, 0.0, 0.02, 0.01, 0.0], 1).unwrap();
            near(elastic_frame.displacements()[3], 0.02 / (young * area), 1e-10);
            near(elastic_frame.displacements()[4], 0.01 / (3.0 * young * moment), 1e-10);
            assert_eq!(elastic_frame.material_section(0), Some(&source));

            let shear = source.clone().with_shear_modulus("shear_modulus", 0.8).unwrap();
            let mut shear_frame = Frame2d::new(Geometry::Linear);
            shear_frame.add_node(0.0, 0.0);
            shear_frame.add_node(1.0, 0.0);
            shear_frame.fix(0, [true; 3]);
            shear_frame.add_material_shear_element(0, 1, 5, shear.clone()).unwrap();
            shear_frame.static_load(&[0.0, 0.0, 0.0, 0.0, 0.01, 0.0], 1).unwrap();
            near(shear_frame.displacements()[4],
                0.01 / (3.0 * young * moment) + 0.01 / shear.shear_stiffness_n(), 2e-10);
            assert_eq!(shear_frame.material_shear_section(0), Some(&shear));

            let mut uniform = MaterialRod::try_straight(1.0, 1, shear.clone()).unwrap();
            uniform.solve_static(&TipLoad {
                force: [0.01, 0.0, 0.0], moment: [0.001, 0.0, 0.0],
            }, 1, 1e-8).unwrap();
            let (strain, curvature) = uniform.rod().strains(0).unwrap();
            near(strain[0], 0.01 / source.axial_stiffness_n(), 3e-5);
            near(curvature[0], 0.001 / shear.torsional_stiffness_n_m2(), 3e-5);
            near(uniform.reference_mass_kg(), density * area, 1e-12);
            assert_eq!(uniform.material_section(), &shear);

            let motion = SpatialMaterialRod::try_straight(1.0, vec![shear]).unwrap()
                .with_tension_links(vec![environment_link(source.clone()).unwrap()], &gate).unwrap()
                .into_motion_model(RodRootBoundary::Free, &gate).unwrap();
            let mut pose = motion.initial_pose();
            pose.positions_m[1][0] += 0.001;
            let response = motion.evaluate(&pose, &gate).unwrap();
            near(response.force_n[1][0], -3.0 * young * area * 0.001, 1e-10);
            near(response.potential_energy_j, 1.5 * young * area * 0.001_f64.powi(2), 1e-10);
            near(motion.inertia().iter().map(|i| i.mass_kg()).sum(),
                3.0 * density * area, 1e-12);
        }
        near(frequencies[1] / frequencies[0], (1.0_f64 / 3.0).sqrt(), 2e-9);
    }
}

#[test]
fn g0_held_environment_rejects_untyped_and_unowned_mechanical_conditions() {
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static,
    ));
    for (quantity, extra) in [
        (None, None),
        (Some(moisture), Some(("strain", [0.0, 0.1], 0.01))),
        (Some(moisture), Some(("strain_rate", [0.0, 1.0], 0.1))),
    ] {
        let source = environment_section("water_content", quantity, [0.0, 0.2], 0.1, extra);
        let shear = source.clone().with_shear_modulus("shear_modulus", 0.8).unwrap();
        assert!(MaterialRod::try_straight(1.0, 1, shear.clone()).is_err());
        assert!(SpatialMaterialRod::try_straight(1.0, vec![shear.clone()]).is_err());
        assert!(environment_link(source.clone()).is_err());
        let mut model = Frame2d::new(Geometry::Linear);
        model.add_node(0.0, 0.0);
        model.add_node(1.0, 0.0);
        assert!(model.add_material_element(0, 1, 3, source).is_err());
        assert!(model.add_material_shear_element(0, 1, 3, shear).is_err());
        assert!(model.material_section(0).is_none());
    }
}

#[test]
fn g0_held_environment_preserves_static_and_transient_frequency_admission() {
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static,
    ));
    let gate = CancelGate::new_clock_free();
    for low in [0.0, 1.0] {
        let source = environment_section("water_content", Some(moisture), [0.0, 0.2],
            0.1, Some(("omega", [low, 100.0], low)));
        let shear = source.clone().with_shear_modulus("shear_modulus", 0.8).unwrap();
        assert_eq!(MaterialRod::try_straight(1.0, 1, shear.clone()).is_ok(), low == 0.0);
        assert_eq!(environment_link(source.clone()).is_ok(), low == 0.0);
        let mut model = Frame2d::new(Geometry::Linear);
        model.add_node(0.0, 0.0);
        model.add_node(1.0, 0.0);
        assert!(model.add_material_element(0, 1, 3, source).is_err());
        assert!(model.add_material_shear_element(0, 1, 3, shear.clone()).is_err());
        let spatial = SpatialMaterialRod::try_straight(1.0, vec![shear]);
        if low == 0.0 {
            assert!(spatial.unwrap().into_motion_model(RodRootBoundary::Free, &gate).is_err());
        } else {
            assert!(spatial.is_err());
        }
    }
}

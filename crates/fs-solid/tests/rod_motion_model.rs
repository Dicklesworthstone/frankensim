//! G0/G1/G3: source inertia, native force reciprocity and explicit release.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};

fn section(band: Option<(f64, f64)>) -> ConcentricShearSection {
    let mut layers = Vec::new();
    for (radius, density) in [(0.01, 1000.0), (0.02, 3000.0)] {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [("density", density, Density::DIMS),
            ("young_modulus", 1e6, Pressure::DIMS), ("shear_modulus", 4e5, Pressure::DIMS)] {
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: band.map_or_else(ValidityDomain::unconstrained,
                    |(lo, hi)| ValidityDomain::unconstrained().with("omega", lo, hi)),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "synthetic rod motion fixture".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims,
                ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
        let point = band.map_or_else(QueryPoint::new,
            |_| QueryPoint::new().with("omega", 0.0).unwrap());
        let material = resolve_material_state_point(&card, &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        layers.push(ConcentricMaterialLayer { outer_radius_m: radius, material });
    }
    ConcentricElasticSection::try_new(0.0, layers).unwrap()
        .with_shear_modulus("shear_modulus", 0.8).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0), "{a:e} != {b:e}");
}
fn v(a: [f64;3]) -> Vec3 { Vec3::new(a[0],a[1],a[2]) }
fn a(v: Vec3) -> [f64;3] { [v.x,v.y,v.z] }

#[test]
fn g0_original_layer_integrals_supply_mass_and_all_three_rotary_axes() {
    let gate = CancelGate::new();
    let rod = SpatialMaterialRod::try_straight(0.8, vec![section(None); 2]).unwrap();
    let reference_mass = rod.reference_mass().total_mass_kg();
    let model = rod.into_motion_model(RodRootBoundary::Free, &gate).unwrap();
    let line_mass = core::f64::consts::PI * (1000.0*0.01_f64.powi(2)
        + 3000.0*(0.02_f64.powi(2)-0.01_f64.powi(2)));
    let diameter = core::f64::consts::PI / 4.0 * (1000.0*0.01_f64.powi(4)
        + 3000.0*(0.02_f64.powi(4)-0.01_f64.powi(4)));
    for (i, weight) in [0.2, 0.4, 0.2].into_iter().enumerate() {
        near(model.inertia()[i].mass_kg(), line_mass*weight, 1e-14);
        for (actual, expected) in model.inertia()[i].principal_kg_m2()
            .into_iter().zip([2.0*diameter*weight,diameter*weight,diameter*weight]) {
            near(actual, expected, 1e-16);
        }
    }
    near(model.inertia().iter().map(|i| i.mass_kg()).sum(), reference_mass, 1e-14);
    let mut pose = model.initial_pose();
    pose.positions_m[2][0] += 0.01;
    model.evaluate(&pose, &gate).unwrap();
    near(model.inertia().iter().map(|i| i.mass_kg()).sum(), reference_mass, 1e-14);
}

#[test]
fn g1_full_native_force_includes_root_and_recovers_axial_law() {
    let gate = CancelGate::new();
    let source = section(None);
    let ea = source.elastic().axial_stiffness_n();
    let model = SpatialMaterialRod::try_straight(0.8, vec![source]).unwrap()
        .into_motion_model(RodRootBoundary::Clamped, &gate).unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.008;
    let response = model.evaluate(&pose, &gate).unwrap();
    near(response.force_n[0][0], ea*0.01, 1e-12);
    near(response.force_n[1][0], -ea*0.01, 1e-12);
    near(response.potential_energy_j, 0.5*ea*0.8*0.01_f64.powi(2), 1e-12);
    pose.positions_m[0][1] += 0.001;
    assert!(model.evaluate(&pose, &gate).is_err());
}

#[test]
fn g3_free_root_and_rigid_frame_change_preserve_native_virtual_work() {
    let gate = CancelGate::new();
    let model = SpatialMaterialRod::try_straight(0.8, vec![section(None);2]).unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate).unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1][1] = 0.002;
    pose.frames[2] = So3::exp(So3Tangent::new(Vec3::new(0.01,0.02,0.03))).unwrap();
    let before = model.evaluate(&pose, &gate).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.4,-0.2,0.1))).unwrap();
    for (point, frame) in pose.positions_m.iter_mut().zip(&mut pose.frames) {
        *point = a(rotation.rotate(v(*point)).unwrap());
        for value in point { *value += 0.7; }
        *frame = rotation.compose(*frame).unwrap();
    }
    let after = model.evaluate(&pose, &gate).unwrap();
    near(before.potential_energy_j, after.potential_energy_j, 1e-11);
    for i in 0..3 {
        let expected = a(rotation.rotate(v(before.force_n[i])).unwrap());
        for c in 0..3 {
            near(after.force_n[i][c], expected[c], 1e-10);
            near(after.torque_body_nm[i][c], before.torque_body_nm[i][c], 1e-10);
        }
    }
}

#[test]
fn g0_centerline_links_add_reference_inventory_but_offsets_require_spatial_inertia() {
    let gate = CancelGate::new();
    let source = section(None);
    for offset in [0.0, 0.01] {
        let rod = SpatialMaterialRod::try_straight(0.8, vec![source.clone()]).unwrap();
        let link = MaterialTensionLink::try_new(source.elastic().clone(), [
            TensionLinkAttachment { node:0,offset_body_m:[0.0,offset,0.0] },
            TensionLinkAttachment { node:1,offset_body_m:[0.0,offset,0.0] },
        ],0.8,2,0.1).unwrap();
        let expected_mass = rod.reference_mass().total_mass_kg()+link.reference_mass_kg();
        let rod = rod.with_tension_links(vec![link], &gate).unwrap();
        let result = rod.into_motion_model(RodRootBoundary::Free, &gate);
        if offset == 0.0 {
            near(result.unwrap().inertia().iter().map(|i| i.mass_kg()).sum(),expected_mass,1e-14);
        } else { assert!(result.is_err()); }
    }
}

#[test]
fn g4_bounded_frequency_sources_malformed_pose_and_cancellation_refuse() {
    let gate = CancelGate::new();
    assert!(SpatialMaterialRod::try_straight(1.0,vec![section(Some((0.0,100.0)))]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).is_err());
    let model = SpatialMaterialRod::try_straight(1.0,vec![section(None)]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m.pop();
    assert!(model.evaluate(&pose,&gate).is_err());
    gate.request();
    assert_eq!(model.evaluate(&model.initial_pose(),&gate),Err(SpatialRodError::Cancelled));
}

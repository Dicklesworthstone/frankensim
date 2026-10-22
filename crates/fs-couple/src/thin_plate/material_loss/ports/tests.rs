//! Synthetic source/geometry checks, not measured instrument validation.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};

fn material(eta: Option<f64>) -> ResolvedMaterialStatePoint {
    let mut values = vec![("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, 200e9), ("poisson_ratio", Dims::NONE, 0.3)];
    if let Some(eta) = eta { values.push((KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS, eta)); }
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic spatial plate ports".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "fixture".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn fixture(eta: Option<f64>) -> ResolvedPlateChart {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 4, 4);
    let inputs = ["left", "right"].into_iter().enumerate().map(|(r, name)| {
        let triangle_indices = mesh.tris.iter().enumerate().filter_map(|(e, tri)| {
            let x = tri.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0;
            ((x >= 0.2) == (r == 1)).then_some(e)
        }).collect();
        PlateRegionMaterial { region: PlateRegion { name: name.into(), triangle_indices }, material: material(eta),
            model: PlateMaterialModel::Isotropic,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.001), material_angle_rad: 0.0 }
    }).collect();
    let s = compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(4, 4), inputs).unwrap();
    if eta.is_some() { s.with_kelvin_voigt_bending_loss().unwrap() } else { s }
}
fn point(node: usize) -> Vec<f64> { let mut v = vec![0.0; 25]; v[node] = 1.0; v }
fn options() -> PlateChartRadiation {
    PlateChartRadiation { assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
        eigenvalue_window: (1.0, 1e8), n_modes: 3, damping_ratio: 0.0, unit_force_weights: point(12) }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1e-12), "{a:e} != {b:e}"); }
fn numeric(r: &CompactBody) -> [f64; 5] { [r.omega, r.zeta, r.mass_kg, r.area_m2, r.drive_participation] }

#[test]
fn reference_modes_match_original_harvesting_with_and_without_material_loss() {
    for eta in [None, Some(1e6)] {
        let s = fixture(eta); let mut o = options();
        if eta.is_none() { o.damping_ratio = 0.002; }
        let model = s.chart().assemble(&[], &o.assembly).unwrap();
        let projection = PlateProjection::with_force_weights(&s.chart().mesh, &o.unit_force_weights).unwrap();
        let loss = s.bending_viscosity().map(|laws|
            crate::plate_viscosity::DiagonalPlateViscosity::new(&s, laws, &model).unwrap());
        let expected = harvest_radiators(&model, &projection, o.eigenvalue_window, o.n_modes,
            o.damping_ratio, loss.as_ref().map(|l| l as &dyn ModalBendingLoss)).unwrap();
        let (actual, ports) = s.certified_radiators_with_force_ports(&o,
            &[o.unit_force_weights.clone(), point(11), point(13)], 10000).unwrap();
        assert_eq!(actual.iter().map(numeric).collect::<Vec<_>>(), expected.iter().map(numeric).collect::<Vec<_>>());
        for (k, r) in actual.iter().enumerate() {
            assert_eq!(ports[0][k], r.drive_participation / fs_math::det::sqrt(r.mass_kg));
        }
        let (same, permuted) = s.certified_radiators_with_force_ports(&o,
            &[point(13), point(11)], 10000).unwrap();
        assert_eq!(same.iter().map(numeric).collect::<Vec<_>>(), actual.iter().map(numeric).collect::<Vec<_>>());
        assert_eq!(permuted[0], ports[2]); assert_eq!(permuted[1], ports[1]);
    }
}

#[test]
fn spatial_ports_preserve_superposition_reciprocity_and_zero_net_force_couples() {
    let s = fixture(Some(1e6)); let o = options();
    let left = point(11); let right = point(13);
    let mixture: Vec<_> = left.iter().zip(&right).map(|(a,b)| 0.25*a + 0.75*b).collect();
    let (modes, b) = s.certified_radiators_with_force_ports(&o, &[left, right, mixture], 10000).unwrap();
    let mut opposed_compliance = 0.0;
    let (mut c00, mut c11, mut c01) = (0.0, 0.0, 0.0);
    for k in 0..modes.len() {
        near(b[2][k], 0.25*b[0][k] + 0.75*b[1][k]);
        let k_eff = modes[k].omega.powi(2); // all attachments are mass-normalized
        opposed_compliance += (b[0][k] - b[1][k]).powi(2) / k_eff;
        c00 += b[0][k].powi(2) / k_eff; c11 += b[1][k].powi(2) / k_eff;
        c01 += b[0][k] * b[1][k] / k_eff;
    }
    assert!(opposed_compliance > 1e-8, "separate forces were collapsed to net force");
    assert!(c00*c11-c01*c01 > 1e-8*c00*c11, "two spatial response ports collapsed to rank one");
    let force = [2.0, -2.0];
    let q: Vec<_> = (0..modes.len()).map(|k|
        (b[0][k]*force[0] + b[1][k]*force[1]) / modes[k].omega.powi(2)
    ).collect();
    let displacement = |port: usize| b[port].iter().zip(&q).map(|(b,q)| b*q).sum::<f64>();
    near(displacement(0), c00*force[0] + c01*force[1]);
    near(displacement(1), c01*force[0] + c11*force[1]);
    let modal_work: f64 = q.iter().enumerate().map(|(k,q)| q*(b[0][k]*force[0]+b[1][k]*force[1])).sum();
    near(modal_work, force[0]*displacement(0) + force[1]*displacement(1));
    assert!(modal_work > 0.0); // a single common footprint would give zero
}

#[test]
fn viscosity_changes_decay_without_replacing_any_spatial_port_or_elastic_mode() {
    let s = fixture(Some(1e6)); let source_identity = s.specimen_identity();
    let changed = s.with_region_material("right", &material(Some(4e6))).unwrap();
    let (a, ba) = s.certified_radiators_with_force_ports(&options(), &[point(11), point(13)], 10000).unwrap();
    let (b, bb) = changed.certified_radiators_with_force_ports(&options(), &[point(11), point(13)], 10000).unwrap();
    assert_eq!(ba, bb); assert_eq!(s.specimen_identity(), source_identity);
    assert!(a.iter().zip(&b).any(|(a,b)| b.zeta > a.zeta));
    for (a,b) in a.iter().zip(&b) {
        assert_eq!([a.omega,a.mass_kg,a.area_m2,a.drive_participation], [b.omega,b.mass_kg,b.area_m2,b.drive_participation]);
    }
}

#[test]
fn invalid_port_sets_and_projection_limits_refuse_before_publication() {
    let s = fixture(None); let o = options();
    let terms = 3*32 + 25*(1+o.n_modes);
    assert!(s.certified_radiators_with_force_ports(&o, &[point(11)], terms-1).is_err());
    s.certified_radiators_with_force_ports(&o, &[point(11)], terms).unwrap();
    s.certified_radiators_with_force_ports(&o, &[], 0).unwrap();
    for bad in [vec![0.0;25], vec![1.0;24], vec![f64::NAN;25]] {
        assert!(s.certified_radiators_with_force_ports(&o, &[point(11),bad], 10000).is_err());
    }
    let mut bad = options(); bad.n_modes = usize::MAX;
    assert!(s.certified_radiators_with_force_ports(&bad, &[point(11)], usize::MAX).is_err());
    let selected = fixture(Some(1e6)); bad = options(); bad.damping_ratio = 0.01;
    assert!(selected.certified_radiators_with_force_ports(&bad, &[point(11)], 10000).is_err());
}

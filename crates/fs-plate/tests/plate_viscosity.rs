//! Physical work on the original spatial DKT mesh, independent of Rayleigh fits.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateModel, PlateSection, assemble,
    viscosity::{PlateBendingViscosity, PlateViscosityBudget, ViscousPlateSection},
};

fn model(mesh: &PlateMesh, tension: f64, density: f64) -> PlateModel {
    assemble(
        mesh,
        &PlateSection::isotropic(3e9, 0.28, 0.008, density).unwrap(),
        &[],
        &[],
        &AssemblyOptions { pretension: tension, support: EdgeSupport::Clamped },
    ).unwrap()
}

fn budget(mesh: &PlateMesh, plate: &PlateModel) -> PlateViscosityBudget {
    PlateViscosityBudget {
        max_elements: mesh.tris.len(), max_dofs: plate.free,
        max_full_dofs: plate.dof_map.len(), max_assembly_entries: 81 * mesh.tris.len(),
    }
}

fn velocity(mesh: &PlateMesh, curvature: [f64; 3], rigid: [f64; 3]) -> Vec<f64> {
    mesh.nodes.iter().flat_map(|&(x, y)| {
        let [a, b, c] = curvature;
        [rigid[0] + rigid[1]*x + rigid[2]*y + 0.5*a*x*x + 0.5*b*y*y + 0.5*c*x*y,
         rigid[1] + a*x + 0.5*c*y, rigid[2] + b*y + 0.5*c*x]
    }).collect()
}

fn power(law: &PlateBendingViscosity<'_>, velocity: &[f64]) -> f64 {
    let mut force = vec![0.0; velocity.len()];
    law.matrix().spmv(velocity, &mut force);
    velocity.iter().zip(&force).map(|(v, f)| v*f).sum()
}

#[test]
fn anisotropic_curvature_power_is_exact_and_rigid_motion_is_lossless() {
    let mesh = PlateMesh::rectangle(1.3, 0.7, 3, 2);
    let plate = model(&mesh, 0.0, 600.0);
    let tensor = [4.0, 1.0, 0.5, 1.0, 9.0, -0.75, 0.5, -0.75, 16.0];
    let loss = PlateBendingViscosity::assemble(
        &mesh, &plate, vec![ViscousPlateSection::from_resultant(tensor).unwrap(); mesh.tris.len()],
        budget(&mesh, &plate), &mut || false,
    ).unwrap();
    let curvature = [0.7, -0.3, 0.2];
    let expected = 1.3*0.7*(0..3).flat_map(|i| (0..3).map(move |j|
        curvature[i]*tensor[3*i+j]*curvature[j])).sum::<f64>();
    for rigid in [[0.0; 3], [1.0, -0.2, 0.7]] {
        let actual = power(&loss, &velocity(&mesh, curvature, rigid));
        assert!((actual-expected).abs() < 2e-11*expected, "{actual} != {expected}");
    }
    assert!(power(&loss, &velocity(&mesh, [0.0; 3], [1.0, -0.2, 0.7])).abs() < 1e-11);
    let directional = PlateBendingViscosity::assemble(
        &mesh, &plate,
        vec![ViscousPlateSection::from_resultant([2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]).unwrap(); mesh.tris.len()],
        budget(&mesh, &plate), &mut || false,
    ).unwrap();
    assert!(power(&directional, &velocity(&mesh, [0.0, 1.0, 0.0], [0.0; 3])).abs() < 1e-11);
    assert!((power(&directional, &velocity(&mesh, [1.0, 0.0, 0.0], [0.0; 3]))-1.82).abs() < 1e-11);
}

#[test]
fn material_loss_does_not_follow_prestress_or_inertia_and_zero_stays_exact() {
    let mesh = PlateMesh::rectangle(1.0, 0.6, 2, 2);
    let plate = model(&mesh, 0.0, 600.0);
    let prestressed = model(&mesh, 2000.0, 1200.0);
    assert_ne!(plate.k.to_dense(), prestressed.k.to_dense());
    assert_ne!(plate.m.to_dense(), prestressed.m.to_dense());
    let section = ViscousPlateSection::from_resultant([4.0,1.0,0.0,1.0,9.0,0.0,0.0,0.0,2.0]).unwrap();
    let original = PlateBendingViscosity::assemble(&mesh, &plate,
        vec![section.clone();mesh.tris.len()], budget(&mesh,&plate), &mut || false).unwrap();
    let changed = PlateBendingViscosity::assemble(&mesh, &prestressed,
        vec![section;mesh.tris.len()], budget(&mesh,&prestressed), &mut || false).unwrap();
    assert_eq!(original.matrix().to_dense(), changed.matrix().to_dense());
    let zero = PlateBendingViscosity::assemble(&mesh, &plate,
        vec![ViscousPlateSection::from_resultant([0.0;9]).unwrap();mesh.tris.len()],
        budget(&mesh,&plate), &mut || false).unwrap();
    assert_eq!(zero.matrix().nnz(), 0);
}

#[test]
fn active_incomplete_overbudget_and_cancelled_viscosity_refuse() {
    assert!(ViscousPlateSection::from_resultant([-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]).is_err());
    let mesh = PlateMesh::rectangle(1.0,0.6,2,2);
    let plate = model(&mesh,0.0,600.0);
    let sections = vec![ViscousPlateSection::from_resultant([0.0;9]).unwrap();mesh.tris.len()];
    assert!(matches!(PlateBendingViscosity::assemble(&mesh,&plate,vec![],
        budget(&mesh,&plate),&mut || false),Err(PlateError::SectionCount {..})));
    let mut limit = budget(&mesh,&plate);
    limit.max_assembly_entries -= 1;
    assert!(matches!(PlateBendingViscosity::assemble(&mesh,&plate,sections.clone(),
        limit,&mut || false),Err(PlateError::Budget {..})));
    let mut polls = 0;
    assert!(matches!(PlateBendingViscosity::assemble(&mesh,&plate,sections,
        budget(&mesh,&plate),&mut || { polls += 1; polls == 3 }),Err(PlateError::Cancelled)));
}

#[test]
fn material_source_and_declared_frequency_band_survive_section_reduction() {
    use fs_evidence::ValidityDomain;
    use fs_material::{
        plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS},
        state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point},
    };
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_qty::DynViscosity;
    let mut claims=ClaimSet::new();
    let mut requirements=Vec::new();
    for (key,value) in KEYS.into_iter().zip([4e6,9e6,16e6,1e6,0.5e6,-0.75e6]) {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(key,DynViscosity::DIMS),
            value:PropertyValue::Scalar {value,dims:DynViscosity::DIMS},
            validity:ValidityDomain::unconstrained().with("omega",10.0,10000.0),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity, observations:vec![],
            provenance:Provenance {source:"synthetic tensor".into(),license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,DynViscosity::DIMS,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {
        chemistry:"synthetic".into(),phase:"solid".into(),process:"viscosity test".into(),revision:0,
    },claims,vec![]).unwrap();
    let state=resolve_material_state_point(&card,&QueryPoint::new().with("omega",100.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let law=PlaneStressViscosity::from_material(state.clone(),MaterialPlane::Axes12,KEYS).unwrap();
    let section=ViscousPlateSection::from_material(law.clone(),0.008,0.37,(20.0,7000.0)).unwrap();
    assert_eq!(section.source().unwrap().source(),&state);
    assert_eq!(section.resolved_band_rad_s(),Some((20.0,7000.0)));
    assert_eq!(section.bending_n_m_s(),&law.bending_matrix_n_m_s(0.008,0.37).unwrap());
    assert!(std::ptr::eq(section.source().unwrap(),section.clone().source().unwrap()),
        "repeated element selections share their immutable material source");
    assert!(ViscousPlateSection::from_material(law.clone(),0.008,0.0,(0.0,7000.0)).is_err());
    assert!(ViscousPlateSection::from_material(law,0.008,0.0,(20.0,12000.0)).is_err());
}

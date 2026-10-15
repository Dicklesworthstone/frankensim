//! Actual DKT pencils and modes from independently characterized sheet data.
use fs_couple::sheet_plate::*;
use fs_material::sheet::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::*;
use fs_evidence::ValidityDomain;
use fs_qty::Dims;
use fs_plate::{PlateMesh,AssemblyOptions,EdgeSupport,PlateSection,PlateSectionField};
fn sheet(mu:f64,bending:f64) -> ResolvedSheet {
    let mut claims=ClaimSet::new();
    for (name,value,dims) in [(AREAL_DENSITY,mu,Dims([-2,1,0,0,0,0])),
        (THICKNESS,0.001,Dims([1,0,0,0,0,0])),
        ("bending_warp",3e-5*bending,BENDING_DIMS),("bending_weft",8e-5*bending,BENDING_DIMS),
        ("bending_coupling",1e-5*bending,BENDING_DIMS),("bending_twist",2e-5*bending,BENDING_DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            interpolation:InterpolationPolicy::ConstantWithinValidity,
            validity:ValidityDomain::unconstrained().with("specimen",3.0,3.0),uncertainty:UncertaintyModel::Unstated,
            observations:vec![],provenance:Provenance {source:"synthetic independent sheet".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"sheet".into(),phase:"solid".into(),process:"synthetic".into(),revision:1},claims,vec![]).unwrap();
    resolve_sheet(&card,&QueryPoint::new().with("specimen",3.0).unwrap(),SheetCapabilities {membrane:false,bending:true},
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn build(s:ResolvedSheet,angle:f64) -> SheetPlate {
    let mesh=PlateMesh::rectangle(0.1,0.07,2,2);let n=mesh.tris.len();
    SheetPlate::assemble(mesh,vec![s;n],vec![angle;n],&PlateMesh::rectangle_boundary(2,2),
        AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        SheetPlateBudget {max_nodes:9,max_triangles:8},||false).unwrap()
}
#[test]
fn exact_existing_assembly_and_material_identity_are_preserved() {
    let s=sheet(0.2,1.0);let p=build(s.clone(),0.31);
    let section=PlateSection {d:*s.bending().unwrap().rotated(0.31).unwrap().coefficients(),
        thickness:s.thickness_m(),density:s.areal_density_kg_m2()/s.thickness_m()};
    let original=fs_plate::assemble_with_sections(p.mesh(),PlateSectionField::Uniform(&section),
        &PlateMesh::rectangle_boundary(2,2),&[],&AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped}).unwrap();
    assert_eq!(p.model().k.val(),original.k.val());assert_eq!(p.model().m.val(),original.m.val());
    assert_eq!(p.model().dof_map,original.dof_map);
    for source in p.sheets() {assert_eq!(source.material(),s.material());}
    assert!((p.reference_mass_kg()-0.0014).abs()<1e-16);
    assert!(p.sheets()[0].membrane().is_err()); // no bulk/membrane data needed to use measured D
}
#[test]
fn cloning_retains_native_operators_sources_and_independent_owned_storage() {
    let original = build(sheet(0.2, 1.0), 0.31);
    let cloned = original.clone();
    assert_eq!(cloned.mesh().nodes, original.mesh().nodes);
    assert_eq!(cloned.mesh().tris, original.mesh().tris);
    assert_eq!(cloned.sheets(), original.sheets());
    assert_eq!(cloned.grain_rad(), original.grain_rad());
    assert_eq!(cloned.reference_mass_kg(), original.reference_mass_kg());
    assert_eq!(cloned.model().k.val(), original.model().k.val());
    assert_eq!(cloned.model().m.val(), original.model().m.val());
    assert_eq!(cloned.model().dof_map, original.model().dof_map);
    assert_ne!(cloned.mesh().nodes.as_ptr(), original.mesh().nodes.as_ptr());
    assert_ne!(cloned.sheets().as_ptr(), original.sheets().as_ptr());
    assert_ne!(cloned.grain_rad().as_ptr(), original.grain_rad().as_ptr());
    drop(original);
    assert!((cloned.reference_mass_kg() - 0.0014).abs() < 1e-16);
}
#[test]
fn data_only_density_and_bending_changes_propagate_into_native_modes() {
    let base=build(sheet(0.2,1.0),0.0);let heavy=build(sheet(0.4,1.0),0.0);let stiff=build(sheet(0.2,4.0),0.0);
    assert_eq!(base.model().k.val(),heavy.model().k.val());assert_eq!(base.model().m.val(),stiff.model().m.val());
    let eigenvalues=|p:&SheetPlate| {
        let options=fs_modal::SliceOptions {mass_diagonal_equilibration:true,..Default::default()};
        let modes=fs_plate::modes(p.model(),(0.0,1e10),&options).unwrap();
        assert_eq!(modes.modes.len(),3);
        modes.modes.iter().map(|m|m.lambda).collect::<Vec<_>>()
    };
    let a=eigenvalues(&base);let b=eigenvalues(&heavy);let c=eigenvalues(&stiff);
    for i in 0..a.len() {
        assert!((b[i]/a[i]-0.5).abs()<1e-7);assert!((c[i]/a[i]-4.0).abs()<1e-7);
    }
}
#[test]
fn grain_changes_the_operator_without_rewriting_topology_or_sources() {
    let source=sheet(0.2,1.0);let a=build(source.clone(),0.0);let b=build(source,0.5);
    assert_ne!(a.model().k.val(),b.model().k.val());assert_eq!(a.model().m.val(),b.model().m.val());
    assert_eq!(a.mesh().tris,b.mesh().tris);assert_eq!(a.sheets(),b.sheets());
    let mesh=a.mesh().clone();let n=mesh.tris.len();
    let attempt=|sources,grains,cancel|SheetPlate::assemble(mesh.clone(),sources,grains,
        &PlateMesh::rectangle_boundary(2,2),AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        SheetPlateBudget {max_nodes:9,max_triangles:8},move ||cancel);
    assert!(attempt(a.sheets()[..n-1].to_vec(),vec![0.0;n],false).is_err());
    assert!(attempt(a.sheets().to_vec(),vec![f64::NAN;n],false).is_err());
    assert!(matches!(attempt(a.sheets().to_vec(),vec![0.0;n],true),Err(SheetPlateError::Cancelled)));
}

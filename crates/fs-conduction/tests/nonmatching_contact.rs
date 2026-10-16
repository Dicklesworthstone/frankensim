//! Different trace meshes, one physical interface: no nodal averaging oracle.
mod support;
use support::with_cx;
use fs_conduction::*;
use fs_conduction::fixtures::{box_grid,on_box_face};
use fs_conduction::interface::{NonmatchingOptions,NonmatchingSurface};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterfaceSystemCard,InterpolationPolicy,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,SelectionPolicy,SurfaceSpec,SystemContext,UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn resistance()->InterfaceResistance {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim{key:PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),value:PropertyValue::Scalar{value:0.1,dims:AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance{source:"analytic nonmatching contact fixture".into(),license:"internal-test-use".into(),artifact:None}}).unwrap();
    let side=|name:&str|SurfaceSpec{material:MaterialStateId{chemistry:name.into(),phase:"solid".into(),
        process:"test".into(),revision:0},texture_frame:"test".into()};
    let card=InterfaceSystemCard::assemble(side("a"),side("b"),SystemContext{medium:"dry".into(),
        third_body:None,environment:"test".into(),history:"test".into()},claims,vec![]).unwrap();
    InterfaceResistance::from_card("bond",&card,&QueryPoint::new(),SelectionPolicy::SingleClaimOnly).unwrap()
}
fn options()->NonmatchingOptions {
    NonmatchingOptions{plane_tolerance_m:1e-12,coverage_relative_tolerance:1e-10,
        max_pair_tests:100_000,max_overlap_triangles:10_000}
}
fn mesh(right:[usize;3],gap:f64)->(ConductionMesh,usize) {
    let (left,mut positions)=box_grid([1,1,1],[1.0;3]);
    let split=positions.len();let (right,points)=box_grid(right,[1.0;3]);
    let mut tets=left.tets;
    tets.extend(right.tets.into_iter().map(|t|t.map(|v|v+split as u32)));
    positions.extend(points.into_iter().map(|p|[p[0]+1.0+gap,p[1],p[2]]));
    (ConductionMesh::new(TetComplex::from_tets(positions.len(),tets),positions).unwrap(),split)
}
fn boundary(mesh:&ConductionMesh)->ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("cold",|f|on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("hot",|f|on_box_face(f.centroid[0],2.0),ThermalBc::dirichlet(330.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}
fn sides(mesh:&ConductionMesh,split:usize,gap:f64)->(Vec<usize>,Vec<usize>) {
    let mut a=vec![];let mut b=vec![];
    for (i,f) in mesh.boundary().iter().enumerate() {
        if f.vertices.iter().all(|&v|(v as usize)<split)
            && on_box_face(f.centroid[0],1.0) {a.push(i);}
        if f.vertices.iter().all(|&v|(v as usize)>=split)
            && on_box_face(f.centroid[0],1.0+gap) {b.push(i);}
    }
    (a,b)
}
fn declaration(a:Vec<usize>,b:Vec<usize>)->NonmatchingSurface {
    NonmatchingSurface::new("bond",a,b,resistance(),options()).unwrap()
}
fn solve_case(cx:&fs_exec::Cx<'_>,mesh:&ConductionMesh,bc:&ThermalBoundary,
    contact:&ThermalInterfaces)->ConductionSolution {
    let k=ConductivityModel::isotropic_declared(10.0).unwrap();
    let source=ScalarField::Uniform(0.0);
    let mut config=SolveConfig::default();config.linear.tolerance=1e-12;
    config.linear.max_iterations=60_000;config.stop.residual_rtol=1e-10;config.stop.step_atol=0.0;
    solve_with_interfaces(cx,ConductionProblem{mesh,boundary:bc,material:&k,
        element_materials:None,source:&source},contact,config).unwrap()
}
#[test]
fn unequal_trace_meshes_recover_the_exact_piecewise_linear_slab_and_contact_jump() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],0.0);let bc=boundary(&mesh);
        let (a,b)=sides(&mesh,split,0.0);assert_eq!(a.len(),2);assert_eq!(b.len(),12);
        let interfaces=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a,b)]).unwrap();
        let solution=solve_case(cx,&mesh,&bc,&interfaces);
        for (i,(p,&t)) in mesh.positions().iter().zip(&solution.temperature).enumerate() {
            let expected=if i<split {300.0+10.0*p[0]}else{310.0+10.0*p[0]};
            assert!((t-expected).abs()<1e-7,"{t} != {expected}");
        }
        let flux=interfaces.fluxes(&solution.temperature).unwrap();
        assert_eq!(flux.len(),1);assert_eq!(interfaces.surface_count(),1);
        assert!((flux[0].area_m2-1.0).abs()<1e-12);
        assert!((flux[0].heat_rate_a_to_b_w+100.0).abs()<1e-6);
        assert!((flux[0].mean_jump_k+10.0).abs()<1e-7);
        assert_eq!(solution.temperature.len(),mesh.vertex_count());
    });
}
#[test]
fn signed_flux_and_nonuniform_pullback_use_overlap_integrals_not_surface_means() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],0.0);let bc=boundary(&mesh);let (a,b)=sides(&mesh,split,0.0);
        let interfaces=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a,b)]).unwrap();
        let t:Vec<_>=mesh.positions().iter().enumerate().map(|(i,p)|if i<split {300.0+p[1]-0.5}else{300.0}).collect();
        let lambda:Vec<_>=mesh.positions().iter().enumerate().map(|(i,p)|if i<split {p[1]-0.5}else{0.0}).collect();
        assert!(interfaces.fluxes(&t).unwrap()[0].heat_rate_a_to_b_w.abs()<1e-10);
        let derivative=interfaces.nonmatching_log_resistance_pullback(cx,"bond",&t,&lambda).unwrap().unwrap();
        assert!((derivative-1.0/1.2).abs()<1e-11);
        let uniform=vec![300.0;mesh.vertex_count()];
        assert_eq!(interfaces.fluxes(&uniform).unwrap()[0].heat_rate_a_to_b_w,0.0);
    });
}
#[test]
fn reversing_sides_preserves_physics_and_reverses_the_reported_heat() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],0.0);let bc=boundary(&mesh);let (a,b)=sides(&mesh,split,0.0);
        let forward=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a.clone(),b.clone())]).unwrap();
        let reverse=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(b,a)]).unwrap();
        let x=solve_case(cx,&mesh,&bc,&forward);let y=solve_case(cx,&mesh,&bc,&reverse);
        for (a,b) in x.temperature.iter().zip(&y.temperature){assert!((a-b).abs()<1e-8);}
        assert!((forward.fluxes(&x.temperature).unwrap()[0].heat_rate_a_to_b_w
            +reverse.fluxes(&y.temperature).unwrap()[0].heat_rate_a_to_b_w).abs()<1e-7);
    });
}
#[test]
fn exactly_matching_subpatches_are_delegated_once_and_retain_the_legacy_solution_bits() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,1,1],0.0);let bc=boundary(&mesh);let (a,b)=sides(&mesh,split,0.0);
        let pairs=ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
        let original=ThermalInterfaces::new(&mesh,&bc,vec![InterfaceSurface::new("bond",pairs,resistance()).unwrap()]).unwrap();
        let general=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a,b)]).unwrap();
        let old=solve_case(cx,&mesh,&bc,&original);let new=solve_case(cx,&mesh,&bc,&general);
        assert_eq!(old.temperature,new.temperature);
        assert_eq!(original.fluxes(&old.temperature).unwrap()[0].area_m2,general.fluxes(&new.temperature).unwrap()[0].area_m2);
        assert_eq!(general.surface_count(),1);
    });
}
#[test]
fn missing_coverage_and_exhausted_geometry_budgets_refuse_before_assembly() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],0.0);let bc=boundary(&mesh);let (a,b)=sides(&mesh,split,0.0);
        let mut missing=b.clone();missing.pop();
        assert!(ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a.clone(),missing)]).is_err());
        for limits in [NonmatchingOptions{max_pair_tests:1,..options()},NonmatchingOptions{max_overlap_triangles:1,..options()}] {
            let surface=NonmatchingSurface::new("bond",a.clone(),b.clone(),resistance(),limits).unwrap();
            assert!(ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![surface]).is_err());
        }
        let mut duplicate=a.clone();duplicate.push(a[0]);
        assert!(NonmatchingSurface::new("bond",duplicate,b,resistance(),options()).is_err());
    });
}
#[test]
fn a_real_gap_is_not_repaired_by_nearest_face_projection() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],1e-4);let (a,b)=sides(&mesh,split,1e-4);
        let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        assert!(ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a,b)]).is_err());
    });
}

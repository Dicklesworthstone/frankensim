//! Different trace meshes, one physical interface: no nodal averaging oracle.
mod support;
#[path = "nonmatching_contact/temperature.rs"]
mod temperature;
#[path = "nonmatching_contact/steady.rs"]
mod steady;
#[path = "nonmatching_contact/steady_linearization.rs"]
mod steady_linearization;
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

// Two independent one-square-metre bondlines in the same plane. Their bulk
// resistances are each 0.1 + 0.1 K/W; distinct contact R'' gives exact spatial
// series solutions without suppressing either body's independent trace mesh.
fn parallel_contact_mesh(right: [[usize;3];2]) -> (ConductionMesh,Vec<bool>) {
    let mut positions=Vec::new();let mut tets=Vec::new();let mut on_a=Vec::new();
    for (patch,subdivision) in right.into_iter().enumerate() {
        for side in 0..2 {
            let (cells,points)=box_grid(if side==0 {[1,1,1]}else{subdivision},[1.0;3]);
            let offset=positions.len() as u32;
            tets.extend(cells.tets.into_iter().map(|tet|tet.map(|v|v+offset)));
            on_a.extend(std::iter::repeat_n(side==0,points.len()));
            positions.extend(points.into_iter().map(|[x,y,z]|[x+side as f64,y+2.0*patch as f64,z]));
        }
    }
    (ConductionMesh::new(TetComplex::from_tets(positions.len(),tets),positions).unwrap(),on_a)
}
fn parallel_sides(mesh:&ConductionMesh)->(Vec<usize>,Vec<usize>) {
    let faces:Vec<_>=mesh.boundary().iter().enumerate()
        .filter(|(_,face)|on_box_face(face.centroid[0],1.0)).collect();
    (faces.iter().filter(|(_,f)|f.outward_normal[0]>0.0).map(|(i,_)|*i).collect(),
        faces.iter().filter(|(_,f)|f.outward_normal[0]<0.0).map(|(i,_)|*i).collect())
}
fn mapped_source_values(scale:f64)->[InterfaceResistance;2] {
    let dims=AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Curve {abscissa:"normal_pressure".into(),abscissa_dims:fs_qty::Pressure::DIMS,
            knots:vec![(1e6,0.3*scale),(2e6,0.1*scale)],dims},
        validity:ValidityDomain::unconstrained().with("normal_pressure",1e6,2e6),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"manufactured pressure-dependent bondline, not measurements".into(),
            license:"internal-test-use".into(),artifact:None}}).unwrap();
    let side=|name:&str|SurfaceSpec {material:MaterialStateId {chemistry:name.into(),phase:"solid".into(),
        process:"fixture".into(),revision:0},texture_frame:"declared contact normal".into()};
    let card=InterfaceSystemCard::assemble(side("a"),side("b"),SystemContext {medium:"dry".into(),
        third_body:None,environment:"vacuum".into(),history:"unaged".into()},claims,vec![]).unwrap();
    [2e6,1e6].map(|pressure|InterfaceResistance::from_card("bond",&card,
        &QueryPoint::new().with_quantity("normal_pressure",fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS),pressure).unwrap(),
        SelectionPolicy::SingleClaimOnly).unwrap())
}
fn mapped_declaration(mesh:&ConductionMesh,scale:f64)->NonmatchingSurface {
    let (mut a,b)=parallel_sides(mesh);let values=mapped_source_values(scale);
    a.reverse(); // Maps follow explicit face slots, never positional sorting.
    let entries=a.iter().map(|&face|(face,values[usize::from(mesh.boundary()[face].centroid[1]>1.5)].clone())).collect();
    NonmatchingSurface::new_mapped("bond",a,b,entries,options()).unwrap()
}
fn near(actual:f64,expected:f64,tolerance:f64) {
    assert!((actual-expected).abs()<=tolerance,"{actual} versus {expected}; tolerance {tolerance}");
}

#[test]
fn g1_mapped_nonmatching_bondlines_close_independent_spatial_series_balances() {
    with_cx(|cx| {
        for right in [[[1,2,3],[1,3,2]],[[1,1,1],[1,1,1]]] {
            let (mesh,on_a)=parallel_contact_mesh(right);let bc=boundary(&mesh);
            let interfaces=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![mapped_declaration(&mesh,1.0)]).unwrap();
            let solution=solve_case(cx,&mesh,&bc,&interfaces);
            for (i,(p,&temperature)) in mesh.positions().iter().zip(&solution.temperature).enumerate() {
                let gradient=if p[1]<1.5 {10.0}else{6.0};
                let expected=if on_a[i] {300.0+gradient*p[0]}else{330.0-gradient*(2.0-p[0])};
                near(temperature,expected,2e-7);
            }
            let flux=interfaces.fluxes(&solution.temperature).unwrap().remove(0);
            near(flux.heat_rate_a_to_b_w,-160.0,2e-6);
            near(flux.conductance_w_per_k,1.0/0.1+1.0/0.3,2e-12);
            near(flux.area_m2,2.0,2e-12);near(flux.mean_jump_k,-12.0,2e-7);
            assert_eq!(interfaces.surface_is_mapped("bond"),Some(true));
            let sources=interfaces.nonmatching_face_resistance_sources("bond").unwrap();
            assert_eq!(sources.len(),parallel_sides(&mesh).0.len());
            assert!(sources.windows(2).all(|pair|pair[0].0<pair[1].0));
            for (face,value) in sources {
                let pressure=if mesh.boundary()[*face].centroid[1]<1.5 {2e6}else{1e6};
                assert!(value.receipt().query_point.contains(&("normal_pressure".into(),pressure)));
            }
            assert_eq!(flux.receipt,*sources[0].1.receipt());
            let material=ConductivityModel::isotropic_declared(10.0).unwrap();
            let source=ScalarField::Uniform(0.0);
            let system=assemble_operator_with_interfaces(cx,&mesh,&bc,&material,&source,&solution.temperature,&interfaces).unwrap();
            let mut reaction=vec![0.0;mesh.vertex_count()];system.operator.spmv(&solution.temperature,&mut reaction);
            near(reaction.iter().sum(),0.0,2e-8);
            for (x,expected) in [(0.0,-160.0),(2.0,160.0)] {
                near(reaction.iter().zip(mesh.positions()).filter(|(_,p)|on_box_face(p[0],x)).map(|(r,_)|r).sum(),expected,2e-6);
            }
        }
    });
}

#[test]
fn g3_constant_nonmatching_map_preserves_uniform_and_exact_delegation_results() {
    with_cx(|cx| {
        for right in [[1,2,3],[1,1,1]] {
            let (mesh,split)=mesh(right,0.0);let bc=boundary(&mesh);let (a,b)=sides(&mesh,split,0.0);
            let uniform=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![declaration(a.clone(),b.clone())]).unwrap();
            let entries=a.iter().map(|&face|(face,resistance())).collect();
            let mapped=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![
                NonmatchingSurface::new_mapped("bond",a,b,entries,options()).unwrap()]).unwrap();
            let original=solve_case(cx,&mesh,&bc,&uniform);let actual=solve_case(cx,&mesh,&bc,&mapped);
            assert_eq!(actual.temperature,original.temperature);
            assert_eq!(mapped.fluxes(&actual.temperature).unwrap(),uniform.fluxes(&actual.temperature).unwrap());
            assert!(uniform.nonmatching_face_resistance_sources("bond").is_none());
        }
    });
}

#[test]
fn g3_mixed_delegated_and_nonmatching_map_scale_matches_fresh_physical_solves() {
    with_cx(|cx| {
        let (mesh,_)=parallel_contact_mesh([[1,1,1],[1,2,3]]);let bc=boundary(&mesh);
        let interfaces=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![mapped_declaration(&mesh,1.0)]).unwrap();
        let material=ConductivityModel::isotropic_declared(10.0).unwrap();let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,element_materials:None,source:&source};
        let response=fs_conduction::adjoint::robin::RobinLinearization::new_with_interfaces(
            cx,problem,&interfaces,SolveConfig::default(),&[]).unwrap();
        let vertices:std::collections::BTreeSet<_>=parallel_sides(&mesh).0.iter()
            .flat_map(|&face|mesh.boundary()[face].vertices).collect();
        let mut weights=vec![0.0;mesh.vertex_count()];
        for vertex in &vertices {weights[*vertex as usize]=1.0/vertices.len() as f64;}
        let dual=response.pullback(cx,&weights,&[],&[]).unwrap();
        let actual=interfaces.nonmatching_log_resistance_pullback(cx,"bond",&response.primal().temperature,&dual.nodal_load).unwrap().unwrap();
        // Each patch's A temperature is 300 + 3/(0.2+sR); both get equal weight.
        let expected=-0.5*(3.0*0.1/(0.3*0.3)+3.0*0.3/(0.5*0.5));
        near(actual,expected,2e-7);
        let objective=|scale| {
            let contacts=ThermalInterfaces::with_nonmatching(cx,&mesh,&bc,vec![],vec![mapped_declaration(&mesh,scale)]).unwrap();
            solve_case(cx,&mesh,&bc,&contacts).temperature.iter().zip(&weights).map(|(t,w)|t*w).sum::<f64>()
        };
        let epsilon=2e-4;
        near(actual,(objective(1.0+epsilon)-objective(1.0-epsilon))/(2.0*epsilon),3e-5);
    });
}

#[test]
fn g0_nonmatching_maps_refuse_incomplete_duplicate_foreign_and_mixed_card_entries() {
    let (mesh,split)=mesh([1,2,3],0.0);let (a,b)=sides(&mesh,split,0.0);
    let values=mapped_source_values(1.0);
    let valid:Vec<_>=a.iter().map(|&face|(face,values[0].clone())).collect();
    let mut incomplete=valid.clone();incomplete.pop();
    let mut duplicate=valid.clone();duplicate[1].0=duplicate[0].0;
    let mut foreign=valid.clone();foreign[0].0=b[0];
    let mut mixed=valid.clone();mixed[1].1=resistance();
    for entries in [vec![],incomplete,duplicate,foreign,mixed] {
        assert!(NonmatchingSurface::new_mapped("bond",a.clone(),b.clone(),entries,options()).is_err());
    }
    assert!(NonmatchingSurface::new_mapped("bond",a,b,valid,options()).is_ok());
}

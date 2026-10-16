#![cfg(feature = "thermal-verification")]
//! Native component means on a full assembly: no artificial region boundaries.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::verification::{ConductionBoundError, FluxBudget, TetError};
use fs_conduction::verification::region::{GoalResidualLimits, RegionMeanConfig,
    bound_temperature_region_mean, solve_with_region_mean_bound};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:73,kernel_id:11,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        f(&cx,&gate)
    })
}
fn config(n: usize) -> (SolveConfig,RegionMeanConfig) {
    let mut primal=SolveConfig::default();
    primal.nonlinearity=fs_conduction::Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:8};
    primal.stop.residual_rtol=1e-11;primal.stop.residual_atol=1e-12;
    primal.linear.tolerance=1e-13;
    let region=RegionMeanConfig {dual:primal.linear,
        residual_limits:GoalResidualLimits {max_rows:n,max_nonzeros:64*n},flux:FluxBudget::default()};
    (primal,region)
}
fn cells(mesh:&ConductionMesh,end:f64)->Vec<usize> {
    mesh.complex().tets.iter().enumerate().filter_map(|(e,tet)|
        tet.iter().all(|&v|mesh.positions()[v as usize][0]<=end).then_some(e)).collect()
}
fn contains(bound:&fs_conduction::verification::region::RegionMeanFieldBound,truth:f64) {
    let interval=bound.bound.enclosure;
    assert!(interval.lo<=truth && truth<=interval.hi,"{interval:?} excludes {truth}");
}

#[test]
fn regional_native_dual_does_not_smear_load_across_shared_vertices() {
    with_cx(|cx,_| {
        let (complex,points)=box_grid([4,2,2],[1.0;3]);
        let mesh=ConductionMesh::new(complex,points).unwrap();
        let boundary=|t|ThermalBoundaryBuilder::new(&mesh)
            .region("ends",|f|on_box_face(f.centroid[0],0.0)||on_box_face(f.centroid[0],1.0),
                ThermalBc::dirichlet(t).unwrap()).unwrap().adiabatic_remainder().finish().unwrap();
        let primal_bc=boundary(300.0);
        let material=ConductivityModel::isotropic_declared(1.0).unwrap();
        let source=ScalarField::Nodal(mesh.positions().iter().map(|p|6.0*p[0]).collect());
        let problem=ConductionProblem {mesh:&mesh,boundary:&primal_bc,material:&material,
            element_materials:None,source:&source};
        let (primal,configuration)=config(mesh.vertex_count());
        let selected=cells(&mesh,0.5);
        let result=solve_with_region_mean_bound(cx,problem,&[],&selected,primal.clone(),configuration).unwrap();
        contains(&result.region,300.21875);
        let native=fs_conduction::solve(cx,problem,primal.clone()).unwrap();
        assert_eq!(result.primal.temperature,native.temperature);
        assert_eq!(result.region.dual_analysis.stability_iterations,0);
        let mut short=configuration;short.dual.max_iterations=1;
        let bounded=bound_temperature_region_mean(cx,problem,&[],&native.temperature,&selected,short).unwrap();
        contains(&bounded,300.21875);
        assert!(bounded.dual_analysis.dual_iterations<=1);
        assert!(result.region.dual_analysis.dual_iterations<=configuration.dual.max_iterations);
        assert!(mesh.positions().iter().zip(&result.region.dual_temperature)
            .any(|(p,z)|p[0]>0.5 && p[0]<1.0 && *z>0.0),"unselected domain still transmits heat");
        // This tempting nodal substitute adds a ramp load in the adjacent slab.
        let wrong_source=ScalarField::Nodal(mesh.positions().iter()
            .map(|p|if p[0]<=0.5 {1.0}else{0.0}).collect());
        let dual_bc=boundary(0.0);let mut dual_config=primal;dual_config.initial=InitialGuess::Uniform(0.0);
        let wrong=fs_conduction::solve(cx,ConductionProblem {boundary:&dual_bc,source:&wrong_source,..problem},dual_config).unwrap();
        assert!(wrong.temperature.iter().zip(&result.region.dual_temperature)
            .any(|(a,b)|(a-b).abs()>1e-5));
        let mut reversed=selected;reversed.reverse();
        let replay=bound_temperature_region_mean(cx,problem,&[],&result.primal.temperature,&reversed,configuration).unwrap();
        assert_eq!(replay.bound.enclosure,result.region.bound.enclosure);
        assert_eq!(replay.dual_temperature,result.region.dual_temperature);
        let unfinished=vec![300.0;mesh.vertex_count()];
        let bound=bound_temperature_region_mean(cx,problem,&[],&unfinished,&reversed,configuration).unwrap();
        contains(&bound,300.21875);
        assert_eq!(unfinished,vec![300.0;mesh.vertex_count()]);
    });
}

#[test]
fn nonuniform_volume_and_original_admission_are_retained() {
    with_cx(|cx,gate| {
        let (complex,mut points)=box_grid([4,2,2],[1.0;3]);
        for p in &mut points {if p[0]==0.25 {p[0]=0.125;}}
        let mesh=ConductionMesh::new(complex,points).unwrap();
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("left",|f|on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap()
            .region("right",|f|on_box_face(f.centroid[0],1.0),ThermalBc::dirichlet(302.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material=ConductivityModel::isotropic_declared(2.0).unwrap();let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,element_materials:None,source:&source};
        let (primal,configuration)=config(mesh.vertex_count());let selected=cells(&mesh,0.5);
        let result=solve_with_region_mean_bound(cx,problem,&[],&selected,primal.clone(),configuration).unwrap();
        contains(&result.region,300.5);
        assert!(result.region.bound.enclosure.hi-result.region.bound.enclosure.lo<1e-4);
        for bad in [vec![],vec![0,0],vec![mesh.element_count()]] {
            assert!(solve_with_region_mean_bound(cx,problem,&[],&bad,primal.clone(),configuration).is_err());
        }
        let mut low=configuration;low.residual_limits.max_rows=0;
        assert!(matches!(solve_with_region_mean_bound(cx,problem,&[],&selected,primal.clone(),low),
            Err(ConductionBoundError::Verification(TetError::Budget))));
        for field in [vec![],vec![f64::NAN;mesh.vertex_count()],vec![300.0;mesh.vertex_count()]] {
            assert!(bound_temperature_region_mean(cx,problem,&[],&field,&selected,configuration).is_err());
        }
        gate.request();
        assert!(matches!(solve_with_region_mean_bound(cx,problem,&[],&selected,primal,configuration),
            Err(ConductionBoundError::Verification(TetError::Cancelled))));
    });
}

fn joint_surface(mesh:&ConductionMesh,r:f64)->Vec<fs_conduction::InterfaceSurface> {
    use fs_matdb::{ClaimSet,InterfaceSystemCard,InterpolationPolicy,MaterialStateId,
        PropertyClaim,PropertyKey,PropertyValue,Provenance,QueryPoint,SelectionPolicy,
        SurfaceSpec,SystemContext,UncertaintyModel};
    use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as DIMS,
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as PROPERTY,InterfaceResistance,InterfaceSurface,ThermalInterfaces};
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(PROPERTY,DIMS),
        value:PropertyValue::Scalar{value:r,dims:DIMS},validity:fs_evidence::ValidityDomain::unconstrained(),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:Vec::new(),provenance:Provenance{source:"regional contact fixture".into(),license:"internal-test-use".into(),artifact:None}}).unwrap();
    let face=|name:&str|SurfaceSpec{material:MaterialStateId{chemistry:name.into(),phase:"solid".into(),
        process:"as-fixtured".into(),revision:0},texture_frame:"contact-normal".into()};
    let card=InterfaceSystemCard::assemble(face("a"),face("b"),SystemContext{medium:"dry".into(),third_body:None,
        environment:"vacuum".into(),history:"unaged".into()},claims,Vec::new()).unwrap();
    let resistance=InterfaceResistance::from_card("joint",&card,&QueryPoint::new(),SelectionPolicy::SingleClaimOnly).unwrap();
    vec![InterfaceSurface::new("joint",ThermalInterfaces::coincident_face_pairs(mesh).unwrap(),resistance).unwrap()]
}
fn two_slabs()->(ConductionMesh,usize) {
    let mut points=Vec::new();let mut tets=Vec::new();let mut per=0;
    for side in 0..2 {
        let (complex,p)=box_grid([2;3],[1.0;3]);per=p.len();
        let offset=points.len() as u32;
        tets.extend(complex.tets.into_iter().map(|tet|tet.map(|v|v+offset)));
        points.extend(p.into_iter().map(|[x,y,z]|[x+side as f64,y,z]));
    }
    (ConductionMesh::new(fs_rep_mesh::TetComplex::from_tets(points.len(),tets),points).unwrap(),per)
}

#[test]
fn hot_component_is_not_diluted_by_cold_solid_and_contact_flux_is_preserved() {
    with_cx(|cx,_| {
        let (mesh,_)=two_slabs();
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("hot",|f|on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(400.0).unwrap()).unwrap()
            .region("cold",|f|on_box_face(f.centroid[0],2.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material=ConductivityModel::isotropic_declared(1.0).unwrap();let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem{mesh:&mesh,boundary:&boundary,material:&material,element_materials:None,source:&source};
        let (primal,configuration)=config(mesh.vertex_count());
        let selected=cells(&mesh,1.0);let surfaces=joint_surface(&mesh,2.0);let before=surfaces.clone();
        let bound=solve_with_region_mean_bound(cx,problem,&surfaces,&selected,primal.clone(),configuration).unwrap();
        // Two unit slabs, k=1 and R=2: q=25; left mean387.5, right312.5, whole350.
        contains(&bound.region,387.5);
        assert!(bound.region.bound.enclosure.lo>387.49);
        let interfaces=fs_conduction::ThermalInterfaces::new(&mesh,&boundary,surfaces.clone()).unwrap();
        let reference=fs_conduction::solve_with_interfaces(cx,problem,&interfaces,primal.clone()).unwrap();
        assert_eq!(reference.temperature,bound.primal.temperature);
        assert_eq!(surfaces,before);
        let right:Vec<_>=(0..mesh.element_count()).filter(|e|!selected.contains(e)).collect();
        let cold=bound_temperature_region_mean(cx,problem,&surfaces,&reference.temperature,&right,configuration).unwrap();
        contains(&cold,312.5);
        assert!(cold.bound.enclosure.hi<bound.region.bound.enclosure.lo);
        assert!(solve_with_region_mean_bound(cx,problem,&[],&selected,primal,configuration).is_err());
    });
}

#[test]
fn region_heating_can_drain_only_through_an_unselected_contacted_component() {
    with_cx(|cx,_| {
        let (mesh,per)=two_slabs();
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("cold",|f|on_box_face(f.centroid[0],2.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material=ConductivityModel::isotropic_declared(1.0).unwrap();
        let source=ScalarField::Nodal(mesh.positions().iter().enumerate().map(|(i,p)|if i<per {6.0*p[0]}else{0.0}).collect());
        let problem=ConductionProblem{mesh:&mesh,boundary:&boundary,material:&material,element_materials:None,source:&source};
        let surfaces=joint_surface(&mesh,2.0);let selected=cells(&mesh,1.0);
        let (primal,configuration)=config(mesh.vertex_count());
        let full=solve_with_region_mean_bound(cx,problem,&surfaces,&selected,primal,configuration).unwrap();
        // Left T=310-x^3; right T=300+3*(2-x). Regional mean is309.75, not305.625.
        contains(&full.region,309.75);
        assert!(full.region.dual_temperature[per..].iter().any(|v|*v>0.0));
        let reused=bound_temperature_region_mean(cx,problem,&surfaces,&full.primal.temperature,&selected,configuration).unwrap();
        assert_eq!(full.region.bound.enclosure,reused.bound.enclosure);
        assert_eq!(full.region.dual_temperature,reused.dual_temperature);
        let incomplete=vec![300.0;mesh.vertex_count()];
        let result=bound_temperature_region_mean(cx,problem,&surfaces,&incomplete,&selected,configuration).unwrap();
        contains(&result,309.75);
        assert!(result.bound.integral.residual_correction.hi>0.0);
        assert_eq!(incomplete,vec![300.0;mesh.vertex_count()]);
    });
}

#![cfg(feature = "thermoelastic-coupling")]
//! Actual source laws, constrained 3-D FEM, mass transport and rigid response.
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_evidence::ValidityDomain;
use fs_ga::{Se3,Vec3,Wrench};
use fs_matdb::{ClaimSet,MaterialCard,MaterialStateId,PropertyClaim,PropertyKey,PropertyValue,
    Provenance,InterpolationPolicy,UncertaintyModel,QueryPoint};
use fs_material::state_point::*;
use fs_mbd::articulated::{ArticulatedModel,BaseState,JointModel,Link,SpatialInertia,forward_dynamics};
use fs_mbd_material::{RegionBinding,ResolvedBody,VolumeMoments,resolve_body};
use fs_mbd_material::deformation::TetrahedralMass;
use fs_mbd_material::deformation::thermoelastic::{ThermalBodyError,ThermalCell,ThermalSolveConfig,ThermalStiffnessModel};
use fs_qty::{Density,Dims,Pressure,Temperature};
use fs_solid::linear3::{TetAssemblyBudget,TetElasticError,TetStaticSolveConfig};

const NODES:[[f64;3];5]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]];
const TETS:[[usize;4];2]=[[0,1,2,3],[0,2,1,4]];
const KEYS:[&str;2]=["a","b"];
const SUPPORT:[usize;6]=[0,1,2,4,5,8];
fn near(a:f64,b:f64) {assert!((a-b).abs()<2e-9*a.abs().max(b.abs()).max(1.),"{a} != {b}");}
fn with_cx<R>(f:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate,arena,StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic))
    })
}
fn card(name:&str,rho:f64,young:[f64;2],alpha:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (key,dims,lo,hi) in [
        (DENSITY_PROPERTY,Density::DIMS,rho,rho),
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,young[0],young[1]),
        (POISSON_RATIO_PROPERTY,Dims([0;6]),0.25,0.25),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS,alpha,alpha),
    ] {
        // Density has ONLY reference-temperature support. A new density query
        // at 400 K would refuse and is unnecessary for this closed inventory.
        let density=key==DENSITY_PROPERTY;
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),
            value:if density {PropertyValue::Scalar {value:lo,dims}} else {PropertyValue::Curve {
                abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots:vec![(300.,lo),(400.,hi)],dims}},
            validity:ValidityDomain::unconstrained().with("T",300.,if density {300.}else{400.})
                .with("pressure",100000.,100000.),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:if density {InterpolationPolicy::ConstantWithinValidity}else{InterpolationPolicy::LinearInside},
            observations:vec![],provenance:Provenance {source:format!("synthetic {name}/{key}; not measured material data"),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:name.into(),phase:"solid".into(),
        process:"synthetic fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn reference(cards:&[MaterialCard])->ResolvedBody {
    let p=QueryPoint::new().with("T",300.).unwrap().with("pressure",100000.).unwrap();
    let bindings:Vec<_>=cards.iter().enumerate().map(|(i,card)|RegionBinding {key:KEYS[i],material:card,point:&p,
        geometry:VolumeMoments::tetrahedron(TETS[i].map(|v|{let p=NODES[v];Vec3::new(p[0],p[1],p[2])})).unwrap(),
        selection:MaterialPropertySelection::SingleClaimOnly}).collect();
    resolve_body(&bindings,2,||false).unwrap()
}
fn inputs<'a>(cards:&'a [MaterialCard],temperatures:&[f64])->Vec<ThermalCell<'a>> {
    cards.iter().zip(temperatures).enumerate().map(|(i,(card,&temperature_k))|ThermalCell {
        key:KEYS[i],card,temperature_k,stiffness_selection:MaterialPropertySelection::SingleClaimOnly,
        stiffness_model:ThermalStiffnessModel::Isotropic,
        expansion_model:fs_mbd_material::deformation::thermoelastic::ThermalExpansionModel::Isotropic,
        expansion_selection:MaterialPropertySelection::SingleClaimOnly}).collect()
}
fn config()->ThermalSolveConfig {
    ThermalSolveConfig {assembly:TetAssemblyBudget {maximum_nodes:5,maximum_tetrahedra:2,
        maximum_free_dofs:15,minimum_scaled_jacobian:1e-10},
        solve:TetStaticSolveConfig {relative_residual_tolerance:1e-11,maximum_iterations:100},
        mesh_frame:ContentHash([91;32]),maximum_absolute_mandel_strain:0.05,maximum_nodal_displacement_m:0.1}
}
fn acceleration(inertia:SpatialInertia)->f64 {
    let model=ArticulatedModel::new(vec![
        Link::new("fixed",None,Se3::identity(),JointModel::FIXED,inertia),
        Link::new("rotor",Some(0),Se3::identity(),JointModel::revolute(Vec3::new(0.,0.,1.),None).unwrap(),inertia),
    ]).unwrap();
    forward_dynamics(&model,BaseState::stationary(Se3::identity()),&[0.],&[0.],&[1.],Vec3::default(),
        &[Wrench::default();2]).unwrap().generalized_acceleration[0]
}

#[test]
fn free_expansion_reaches_actual_nodes_and_dynamics_without_replacing_mass() {
    let cards=[card("a",6.,[1000.,700.],1e-4)];let original=reference(&cards);
    let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let result=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.]),&SUPPORT,config(),cx)).unwrap();
    for (p,q) in NODES[..4].iter().zip(result.body().nodes_m()) {
        for i in 0..3 {near(q[i],1.01*p[i]);}
    }
    near(result.stress().elastic_energy_j(),0.);
    assert_eq!(result.body().inertia().mass().to_bits(),original.inertia().mass().to_bits());
    near(acceleration(result.body().inertia())/acceleration(original.inertia()),1./1.01_f64.powi(2));
    near(result.sources()[0].stiffness().property(YOUNG_MODULUS_PROPERTY).unwrap().value_si(),700.);
    assert!(result.sources()[0].stiffness().property(DENSITY_PROPERTY).is_none());
    assert_eq!(result.geometry().nodes_m(),result.body().nodes_m());
    assert_eq!(result.stress().solution_identity(),result.displacement().identity());
}

#[test]
fn constrained_material_change_produces_reactions_not_fictitious_expansion() {
    for young in [1000.,2000.] {
        let cards=[card("a",6.,[young,young],1e-4)];let original=reference(&cards);
        let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
        let fixed:Vec<_>=(0..12).collect();
        let result=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.]),&fixed,config(),cx)).unwrap();
        assert_eq!(result.body().inertia(),original.inertia());
        assert_eq!(result.body().nodes_m(),&NODES[..4]);
        for &s in &result.stress().elements()[0].stress_mandel_pa()[..3] {near(s,-young*0.01/(1.-2.*0.25));}
        assert!(result.stress().elastic_energy_j()>0.);
        for axis in 0..3 {near(result.stress().nodal_internal_force_n().iter().map(|f|f[axis]).sum(),0.);}
        assert_eq!(result.displacement().iterations(),0);
    }
}

#[test]
fn heterogeneous_bonded_cells_follow_independent_compatible_solution() {
    for top_e in [1000.,2000.] {
        let cards=[card("a",6.,[top_e,top_e],1e-4),card("b",12.,[2000.,2000.],1e-4)];
        let original=reference(&cards);
        let mesh=TetrahedralMass::bind(&original,&NODES,&TETS,&KEYS,5,||false).unwrap();
        let result=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[400.,300.]),&SUPPORT,config(),cx)).unwrap();
        // Independent P1 patch equations, equal reference volumes and nu=1/4.
        let shared=top_e*0.01/(top_e+2000.);
        near(result.displacement().displacement_m()[1][0],shared);
        near(result.displacement().displacement_m()[2][1],shared);
        near(result.displacement().displacement_m()[3][2],0.01+2./3.*(0.01-shared));
        near(result.displacement().displacement_m()[4][2],2./3.*shared);
        assert!(shared>0. && shared<0.01);
        // The unheated cell moves and is stressed by its neighbour.
        assert!(result.stress().elements()[1].elastic_energy_j()>0.);
        near(result.stress().elements()[0].stress_mandel_pa()[0],top_e/(1.-0.25)*(shared-0.01));
        near(result.stress().elements()[1].stress_mandel_pa()[0],2000./(1.-0.25)*shared);
        for (a,b) in original.regions().iter().zip(result.body().regions()) {assert_eq!(a.inertia().mass(),b.inertia().mass());}
        let i=result.body().inertia();let c=i.center_of_mass();
        near(acceleration(i),1./(i.inertia_com().m[8]+i.mass()*(c.x*c.x+c.y*c.y)));
    }
}

#[test]
fn input_order_reset_and_original_source_context_are_preserved() {
    let cards=[card("a",6.,[1000.,900.],1e-4),card("b",12.,[2000.,1500.],1e-4)];let original=reference(&cards);
    let mesh=TetrahedralMass::bind(&original,&NODES,&[TETS[1],TETS[0]],&["b","a"],5,||false).unwrap();
    let mut samples=inputs(&cards,&[400.,350.]);
    let hot=with_cx(|cx|mesh.solve_thermal(&samples,&SUPPORT,config(),cx)).unwrap();samples.reverse();
    assert_eq!(hot,with_cx(|cx|mesh.solve_thermal(&samples,&SUPPORT,config(),cx)).unwrap());
    let reset=with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[300.,300.]),&SUPPORT,config(),cx)).unwrap();
    assert_eq!(reset.body().inertia(),original.inertia());
    for (source,original) in hot.sources().iter().zip(original.regions()) {
        assert_eq!(source.expansion().reference().query_point(),original.material().query_point());
        assert_eq!(source.stiffness().query_point().iter().find(|(a,_)|a=="pressure").unwrap().1,100000.);
    }
}

#[test]
fn bad_sources_support_limits_and_geometry_limits_publish_no_solution() {
    let cards=[card("a",6.,[1000.,900.],1e-4)];let original=reference(&cards);let before=original.clone();
    let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    for t in [0.,f64::NAN,401.] {assert!(with_cx(|cx|mesh.solve_thermal(&inputs(&cards,&[t]),&SUPPORT,config(),cx)).is_err());}
    let foreign=[card("other",6.,[1000.,900.],1e-4)];
    assert!(with_cx(|cx|mesh.solve_thermal(&inputs(&foreign,&[350.]),&SUPPORT,config(),cx)).is_err());
    let samples=inputs(&cards,&[400.]);
    assert!(with_cx(|cx|mesh.solve_thermal(&[],&SUPPORT,config(),cx)).is_err());
    assert!(with_cx(|cx|mesh.solve_thermal(&samples,&[],config(),cx)).is_err());
    let tiny=ThermalSolveConfig {maximum_nodal_displacement_m:1e-4,..config()};
    assert!(matches!(with_cx(|cx|mesh.solve_thermal(&samples,&SUPPORT,tiny,cx)),Err(ThermalBodyError::Solid(_))));
    let strain=ThermalSolveConfig {maximum_absolute_mandel_strain:1e-4,..config()};
    assert!(with_cx(|cx|mesh.solve_thermal(&samples,&SUPPORT,strain,cx)).is_err());
    assert_eq!(before,original);
}

#[test]
fn cancellation_before_assembly_preserves_the_original_inventory() {
    let cards=[card("a",6.,[1000.,900.],1e-4)];let original=reference(&cards);
    let mesh=TetrahedralMass::bind(&original,&NODES[..4],&TETS[..1],&KEYS[..1],4,||false).unwrap();
    let gate=CancelGate::new();gate.request();
    let result=ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:3,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        mesh.solve_thermal(&inputs(&cards,&[400.]),&SUPPORT,config(),&cx)
    });
    assert!(matches!(result,Err(ThermalBodyError::Solid(TetElasticError::Cancelled))));
    near(original.inertia().mass(),1.);
}

#[path = "thermoelastic_body/heating.rs"]
mod heating;

#[path = "thermoelastic_body/vibration.rs"]
mod vibration;

#[path = "thermoelastic_body/loads.rs"]
mod loads;

#[path = "thermoelastic_body/orthotropic.rs"]
mod orthotropic;

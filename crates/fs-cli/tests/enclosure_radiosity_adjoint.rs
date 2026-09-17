//! Direct tests of the same radiosity pullback consumed by coupled cooling.
#[path="../src/network_command/radiation/enclosure/radiosity_adjoint.rs"]
mod derivative;
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_conduction::{ConductionMesh,SurfaceEmissivity,SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS};
use fs_conduction::radiation::{GrayDiffuseEnclosure,RadiationSurface,ViewFactorMatrix,ViewFactorEvidence,ViewFactorTolerance};
use fs_conduction::fixtures::box_grid;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,Provenance,UncertaintyModel,
    InterpolationPolicy,MaterialCard,MaterialStateId,SelectionPolicy};

fn context<T>(gate:&CancelGate,f:impl FnOnce(&Cx<'_>)->T)->T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena|f(&Cx::new(gate,arena,
        StreamKey {seed:17,kernel_id:717,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn epsilon(name:&str,value:f64)->SurfaceEmissivity {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar {value,dims:EMISSIVITY_DIMS},validity:ValidityDomain::unconstrained(),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"declared test finish".into(),license:"test".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"fixture".into(),phase:"solid".into(),
        process:"declared".into(),revision:0},claims,vec![]).unwrap();
    SurfaceEmissivity::from_card(name,&card,300.0,SelectionPolicy::SingleClaimOnly).unwrap()
}
fn model(e:[f64;3])->GrayDiffuseEnclosure {
    let (complex,positions)=box_grid([1,1,1],[1.0,2.0,3.0]);
    let mesh=ConductionMesh::new(complex,positions).unwrap();
    let surfaces=(0..3).map(|axis|RadiationSurface::new(&mesh,format!("face-{axis}"),
        |face|face.centroid[axis]==0.0,epsilon(&format!("face-{axis}"),e[axis])).unwrap()).collect::<Vec<_>>();
    let area:Vec<_>=surfaces.iter().map(|s|s.area_m2()).collect();
    let exchanges=[[0.0,0.9,0.4],[0.9,0.0,0.2],[0.4,0.2,0.0]];
    let mut f=vec![vec![0.0;3];3];
    for i in 0..3 {for j in 0..3 {if i!=j {f[i][j]=exchanges[i][j]/area[i];}}
        f[i][i]=1.0-f[i].iter().sum::<f64>();}
    let factors=ViewFactorMatrix::admit(area,f,ViewFactorEvidence::Analytic {
        geometry:"synthetic reciprocal unequal-area test matrix; not geometric visibility".into()},
        ViewFactorTolerance::default()).unwrap();
    GrayDiffuseEnclosure::new(surfaces,factors).unwrap()
}
fn functional(cx:&Cx<'_>,e:[f64;3],t:[f64;3],w:[f64;3])->f64 {
    model(e).solve(cx,&t).unwrap().net_outward_flux_w_m2.iter().zip(w).map(|(a,b)|a*b).sum()
}
fn near(a:f64,b:f64,tol:f64){assert!((a-b).abs()<tol*b.abs().max(1.0),"{a} versus {b}");}

#[test]
fn nonsymmetric_reflection_transpose_matches_temperature_and_finish_perturbations() {
    context(&CancelGate::new_clock_free(),|cx| {
        let e=[0.2,0.6,0.9];let t=[320.0,290.0,310.0];let w=[0.7,-0.3,0.2];
        let linear=derivative::Linearization::new(cx,&model(e),&t,1e-8).unwrap();
        let g=linear.pullback(cx,&w,1e-11,100).unwrap();
        assert!(g.relative_residual<=1e-11);assert!(g.iterations>0);
        for i in 0..3 {
            let h=1e-3;let mut plus=t;let mut minus=t;plus[i]+=h;minus[i]-=h;
            near(g.temperatures[i],(functional(cx,e,plus,w)-functional(cx,e,minus,w))/(2.0*h),2e-7);
            let h=1e-5_f64;let mut plus=e;let mut minus=e;plus[i]*=h.exp();minus[i]*=(-h).exp();
            near(g.log_emissivities[i],(functional(cx,plus,t,w)-functional(cx,minus,t,w))/(2.0*h),2e-7);
        }
    });
}

#[test]
fn black_surfaces_and_closed_energy_objectives_do_not_create_singularities_or_heat() {
    context(&CancelGate::new_clock_free(),|cx| {
        let e=[1.0,0.4,0.7];let t=[320.0,290.0,310.0];let w=[1.0,0.0,0.0];
        let enclosure=model(e);let linear=derivative::Linearization::new(cx,&enclosure,&t,1e-8).unwrap();
        let g=linear.pullback(cx,&w,1e-11,100).unwrap();
        let h=1e-6_f64;let mut minus=e;minus[0]*=(-h).exp();
        near(g.log_emissivities[0],(functional(cx,e,t,w)-functional(cx,minus,t,w))/h,2e-5);
        let energy=linear.pullback(cx,enclosure.view_factors().areas_m2(),1e-11,100).unwrap();
        for value in energy.temperatures.iter().chain(&energy.log_emissivities) {assert!(value.abs()<1e-8);}
    });
}

#[test]
fn insufficient_linear_budget_invalid_weights_and_cancellation_return_no_partial_derivative() {
    context(&CancelGate::new_clock_free(),|cx| {
        let linear=derivative::Linearization::new(cx,&model([0.2,0.6,0.9]),&[320.0,290.0,310.0],1e-8).unwrap();
        assert!(linear.pullback(cx,&[0.7,-0.3,0.2],1e-12,1).is_err());
        assert!(linear.pullback(cx,&[1.0],1e-12,100).is_err());
        assert!(linear.pullback(cx,&[f64::NAN,0.0,0.0],1e-12,100).is_err());
        let gate=CancelGate::new_clock_free();gate.request();
        context(&gate,|stopped|assert!(linear.pullback(stopped,&[0.0;3],1e-12,100).is_err()));
        let zero=linear.pullback(cx,&[0.0;3],1e-12,100).unwrap();
        assert_eq!(zero.temperatures,vec![0.0;3]);assert_eq!(zero.iterations,0);
    });
}

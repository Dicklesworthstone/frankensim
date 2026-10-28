use super::*;
use super::super::tests::fixture;

fn geometry()->SurfaceDiffusionGeometry {
    SurfaceDiffusionGeometry { connecting_width_m:1.0,distance_m:1e-6,site_density_mol_m2:1.0 }
}
fn selections(model:&CompetitiveAdsorbent)->Vec<SurfaceDiffusionSelection> {
    model.channels.iter().map(|law|SurfaceDiffusionSelection { species:law.species.clone(),
        selection:MaterialPropertySelection::SingleClaimOnly }).collect()
}
fn link(solid:&CompetitiveSolid)->SurfaceDiffusionLink {
    SurfaceDiffusionLink::try_new(Arc::clone(&solid.model),solid.temperature_k,geometry(),&selections(&solid.model),||false).unwrap()
}
fn controls()->SurfaceDiffusionSolve {
    SurfaceDiffusionSolve { maximum_net_species_fraction:0.5,maximum_site_turnover:0.5,
        mass_tolerance_kg:1e-13,energy_tolerance_j:1e-9,entropy_tolerance_j_k:1e-10 }
}
fn near(a:f64,b:f64,tolerance:f64) { assert!((a-b).abs()<=tolerance,"{a:.17e} != {b:.17e}"); }
#[test]
fn finite_surface_transport_matches_independent_dense_equations_and_conserves_energy() {
    let (left,_,_)=fixture([0.03,0.005]); let (right,_,_)=fixture([0.005,0.02]);
    let saved=(left.clone(),right.clone()); let path=link(&left); let dt=0.1;
    let candidate=path.advance_candidate(&left,&right,dt,controls(),||false).unwrap();
    // Independently assemble both implicit equations from G*(theta_i theta_v).
    let c=0.05; let total=[0.035,0.025]; let v=0.04; let g=[0.01/(c*c),0.02/(c*c)];
    let a=1.0+dt*g[0]*(v+total[0]); let b=dt*g[0]*total[0];
    let d=1.0+dt*g[1]*(v+total[1]); let e=dt*g[1]*total[1];
    let rhs=[0.03+dt*g[0]*total[0]*c,0.005+dt*g[1]*total[1]*c]; let determinant=a*d-b*e;
    near(candidate.left().moles[0],(d*rhs[0]-b*rhs[1])/determinant,2e-14);
    near(candidate.left().moles[1],(a*rhs[1]-e*rhs[0])/determinant,2e-14);
    for i in 0..2 { near(candidate.left().moles[i]+candidate.right().moles[i],total[i],2e-16); }
    assert_eq!(candidate.left().temperature_k,300.0); assert_eq!(candidate.right().temperature_k,300.0);
    assert!(candidate.step().entropy_change_j_k>0.0);
    near(candidate.step().energy_residual_j,0.0,1e-10);
    assert!(candidate.step().species.iter().all(|row|row.kinetic_mass_residual_kg.abs()<1e-13));
    near(path.coefficients()[0].diffusivity_m2_s,1e-8,0.0);
    near(path.coefficients()[1].conductance_mol_s,0.02,1e-17);
    assert_eq!((left,right),saved);
}
#[test]
fn unequal_vacancies_can_drive_uphill_individual_surface_loading() {
    let (left,_,_)=fixture([0.01,0.0]); let (right,_,_)=fixture([0.005,0.04]);
    let result=link(&left).advance_candidate(&left,&right,0.001,controls(),||false).unwrap();
    // A is more concentrated on the left, but its activity theta_A/theta_v
    // is lower there. Independent concentration-difference Fick rows fail.
    assert!(left.moles[0]>right.moles[0]);
    assert!(result.step().species[0].moles_left_to_right<0.0);
    assert!(result.left().moles[0]>left.moles[0]);
    assert!(result.step().entropy_change_j_k>0.0);
}
#[test]
fn scalar_fick_limit_and_first_order_refinement() {
    let (fixture,_,_)=fixture([0.0,0.0]);
    let model=Arc::new(CompetitiveAdsorbent::try_new(vec![fixture.model.channels[0].clone()]).unwrap());
    let left=CompetitiveSolid::ordered(Arc::clone(&model),0.1,vec![0.03],300.0).unwrap();
    let right=CompetitiveSolid::ordered(model,0.1,vec![0.01],300.0).unwrap(); let path=link(&left);
    let exact=0.02+0.01*(-0.4_f64).exp();
    let run=|count:usize| {
        let (mut l,mut r)=(left.clone(),right.clone());
        for _ in 0..count { let (nl,nr,_)=path.advance_candidate(&l,&r,1.0/count as f64,controls(),||false).unwrap().into_parts(); l=nl; r=nr; }
        (l.moles[0]-exact).abs()
    };
    let ratio=run(32)/run(64); assert!(ratio>1.9 && ratio<2.1);
}
#[test]
fn swapped_owners_reverse_transport_and_different_capacities_remain_conservative() {
    let (left,_,_)=fixture([0.03,0.005]); let (right,_,_)=fixture([0.005,0.02]);
    let right=CompetitiveSolid::ordered(Arc::clone(&right.model),0.2,right.moles,300.0).unwrap();
    let path=link(&left); let a=path.advance_candidate(&left,&right,0.1,controls(),||false).unwrap();
    let b=path.advance_candidate(&right,&left,0.1,controls(),||false).unwrap();
    for i in 0..2 {
        near(a.left().moles[i],b.right().moles[i],2e-14);
        near(a.right().moles[i],b.left().moles[i],2e-14);
        near(a.step().species[i].moles_left_to_right,-b.step().species[i].moles_left_to_right,2e-14);
    }
    near(a.step().entropy_change_j_k,b.step().entropy_change_j_k,2e-13);
}
#[test]
fn no_vacancy_zero_width_and_zero_time_preserve_both_phases() {
    let (left,_,_)=fixture([0.03,0.02]); let (right,_,_)=fixture([0.01,0.04]);
    let path=link(&left); let a=path.advance_candidate(&left,&right,1.0,controls(),||false).unwrap();
    assert_eq!(a.left(),&left); assert_eq!(a.right(),&right); assert_eq!(a.step().site_turnover,0.0);
    let (left,_,_)=fixture([0.03,0.005]); let (right,_,_)=fixture([0.005,0.02]);
    let mut g=geometry(); g.connecting_width_m=0.0;
    let closed=SurfaceDiffusionLink::try_new(Arc::clone(&left.model),300.0,g,&selections(&left.model),||false).unwrap();
    let a=closed.advance_candidate(&left,&right,1.0,controls(),||false).unwrap();
    assert_eq!(a.left(),&left); assert_eq!(a.right(),&right);
    let a=link(&left).advance_candidate(&left,&right,0.0,controls(),||false).unwrap();
    assert_eq!(a.left(),&left); assert_eq!(a.right(),&right);
}
#[test]
fn missing_mobility_thermal_drift_and_late_failures_publish_no_endpoint() {
    let (left,_,_)=fixture([0.03,0.005]); let (right,_,_)=fixture([0.005,0.02]); let path=link(&left);
    let saved=(left.clone(),right.clone());
    assert!(SurfaceDiffusionLink::try_new(Arc::clone(&left.model),300.0,geometry(),&[],||false).is_err());
    let mut duplicate=selections(&left.model); duplicate[1].species=duplicate[0].species.clone();
    assert!(SurfaceDiffusionLink::try_new(Arc::clone(&left.model),300.0,geometry(),&duplicate,||false).is_err());
    let warm=CompetitiveSolid::ordered(Arc::clone(&right.model),right.dry_mass_kg,right.moles.clone(),301.0).unwrap();
    assert!(path.advance_candidate(&left,&warm,0.1,controls(),||false).is_err());
    let mut c=controls(); c.maximum_net_species_fraction=1e-12;
    assert!(path.advance_candidate(&left,&right,0.1,c,||false).is_err());
    let mut polls=0; assert!(matches!(path.advance_candidate(&left,&right,0.1,controls(),||{polls+=1;polls>9}),Err(AdsorptionError::Cancelled)));
    assert_eq!((left,right),saved);
}

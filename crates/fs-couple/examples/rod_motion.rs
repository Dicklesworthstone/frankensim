//! Manufactured free-flight rod; no measured material or vehicle qualification.
//! cargo run -p fs-couple --features rod-dynamics --example rod_motion
//! cargo run -p fs-couple --features rod-dynamics --example rod_motion -- --eccentric --refinement 2
use std::sync::Arc;
use fs_couple::rod_motion::{RodMotion,RodMotionControls,RodMotionInertia,RodNodalWrench,SpatialRodMotion};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection,ScalarAdmissibility,
    ScalarPropertyRequirement,resolve_material_state_point};
use fs_qty::{Density,Pressure};
use fs_solid::rod::material_field::{SpatialMaterialRod,motion::RodRootBoundary};
use fs_solid::rod::material_field::tension::{MaterialTensionLink,TensionLinkAttachment};

type Error=Box<dyn std::error::Error>;
fn invalid(message:&str)->Error {std::io::Error::new(std::io::ErrorKind::InvalidInput,message).into()}

fn main()->Result<(),Error> {
    let mut eccentric=false;let mut refinement=1usize;
    let mut arguments=std::env::args().skip(1);
    while let Some(argument)=arguments.next() {
        match argument.as_str() {
            "--eccentric"=>eccentric=true,
            "--refinement"=>{
                refinement=arguments.next().ok_or_else(||invalid("--refinement needs an integer in 1..=32"))?.parse()?;
                if !(1..=32).contains(&refinement) {return Err(invalid("--refinement must be in 1..=32"));}
            }
            _=>return Err(invalid("usage: rod_motion [--eccentric] [--refinement 1..=32]")),
        }
    }
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),("young_modulus",1e4,Pressure::DIMS),
        ("shear_modulus",4e3,Pressure::DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {source:"manufactured rod flight".into(),license:"example-only".into(),artifact:None}})?;
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive)?);
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"manufactured elastic solid".into(),
        phase:"solid".into(),process:"example".into(),revision:0},claims,vec![])?;
    let material=resolve_material_state_point(&card,&QueryPoint::new(),&requirements,MaterialPropertySelection::SingleClaimOnly)?;
    let section=ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {outer_radius_m:0.05,material}])?
        .with_shear_modulus("shear_modulus",0.8)?;
    let gate=CancelGate::new();
    let rod=SpatialMaterialRod::try_straight(0.8,vec![section.clone();4])?;
    let velocity=[[0.1,0.0,0.0];5];let omega=[[0.2,0.1,0.0];5];
    if eccentric {
        // A slightly pre-tensioned bundle connects root and tip through actual
        // material-frame offsets. Its original mass stays at those attachments.
        let link=MaterialTensionLink::try_new(section.elastic().clone(),[
            TensionLinkAttachment {node:0,offset_body_m:[0.03,0.04,-0.02]},
            TensionLinkAttachment {node:4,offset_body_m:[0.03,0.04,-0.02]},
        ],0.795,1,0.1)?;
        let model=Arc::new(rod.with_tension_links(vec![link],&gate)?
            .into_spatial_motion_model(RodRootBoundary::Free,&gate)?);
        let mut pose=model.initial_pose();pose.positions_m[2][1]=0.01;
        let state=SpatialRodMotion::new_spatial(model,pose,&velocity,&omega,0.0,1e-12,&gate)?;
        advance(state,refinement,"eccentric",&gate)
    } else {
        let model=Arc::new(rod.into_motion_model(RodRootBoundary::Free,&gate)?);
        let mut pose=model.initial_pose();pose.positions_m[2][1]=0.01;
        let state=RodMotion::new(model,pose,&velocity,&omega,0.0,&gate)?;
        advance(state,refinement,"centered",&gate)
    }
}

fn advance<I:RodMotionInertia>(mut state:RodMotion<I>,refinement:usize,mode:&str,gate:&CancelGate)->Result<(),Error> {
    let dt=1e-4/refinement as f64;let steps=200*refinement;
    let controls=RodMotionControls {maximum_step_s:dt,maximum_steps:steps,
        energy_defect_limit_j:1e-5,maximum_rotation_rad:0.25,balance_origin_m:[0.0;3]};
    let initial=state.observe(controls.balance_origin_m)?;
    let energy=initial.kinetic_energy_j+initial.potential_energy_j;
    let loads=vec![RodNodalWrench::default();state.pose().positions_m.len()];
    eprintln!("mode={mode}, dt_s={dt:e}, steps={steps}, inertia_tensor_residual={:e}; coefficients are manufactured",
        state.maximum_inertia_tensor_residual());
    println!("time_s,kinetic_j,potential_j,total_energy_change_j,step_energy_defect_j,max_linear_balance_ns,max_angular_balance_nms,tip_x_m,tip_y_m,tip_z_m");
    for _ in 0..steps {
        let (next,r)=state.step_candidate(dt,&loads,controls,gate)?;
        let lp=r.linear_balance_residual_n_s.into_iter().map(f64::abs).fold(0.0,f64::max);
        let ap=r.angular_balance_residual_nm_s.into_iter().map(f64::abs).fold(0.0,f64::max);
        let tip=next.pose().positions_m.last().ok_or_else(||invalid("admitted rod unexpectedly has no tip"))?;
        println!("{:.8},{:.12e},{:.12e},{:.6e},{:.6e},{:.6e},{:.6e},{:.12e},{:.12e},{:.12e}",next.time_s(),
            r.after.kinetic_energy_j,r.after.potential_energy_j,
            r.after.kinetic_energy_j+r.after.potential_energy_j-energy,r.integration_defect_j,lp,ap,tip[0],tip[1],tip[2]);
        state=next;
    }
    Ok(())
}

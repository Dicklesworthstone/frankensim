//! Manufactured free-flight rod. Coefficients are not a measured material pack.
//! cargo run -p fs-couple --features rod-dynamics --example rod_motion
use std::sync::Arc;
use fs_couple::rod_motion::{RodMotion,RodMotionControls,RodNodalWrench};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection,ScalarAdmissibility,
    ScalarPropertyRequirement,resolve_material_state_point};
use fs_qty::{Density,Pressure};
use fs_solid::rod::material_field::{SpatialMaterialRod,motion::RodRootBoundary};

fn main()->Result<(),Box<dyn std::error::Error>> {
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
    let model=Arc::new(SpatialMaterialRod::try_straight(0.8,vec![section;4])?
        .into_motion_model(RodRootBoundary::Free,&gate)?);
    let mut pose=model.initial_pose();pose.positions_m[2][1]=0.01;
    let mut state=RodMotion::new(model,pose,&[[0.1,0.0,0.0];5],&[[0.2,0.1,0.0];5],0.0,&gate)?;
    let controls=RodMotionControls {maximum_step_s:1e-4,maximum_steps:200,
        energy_defect_limit_j:1e-5,maximum_rotation_rad:0.25,balance_origin_m:[0.0;3]};
    println!("time_s,kinetic_j,potential_j,step_energy_defect_j,max_linear_balance_ns,max_angular_balance_nms");
    for _ in 0..200 {
        let (next,r)=state.step_candidate(1e-4,&[RodNodalWrench::default();5],controls,&gate)?;
        let lp=r.linear_balance_residual_n_s.into_iter().map(f64::abs).fold(0.0,f64::max);
        let ap=r.angular_balance_residual_nm_s.into_iter().map(f64::abs).fold(0.0,f64::max);
        println!("{:.8},{:.12e},{:.12e},{:.6e},{:.6e},{:.6e}",next.time_s(),
            r.after.kinetic_energy_j,r.after.potential_energy_j,r.integration_defect_j,lp,ap);
        state=next;
    }
    Ok(())
}

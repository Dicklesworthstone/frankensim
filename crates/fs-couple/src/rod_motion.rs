//! Native large-rotation rod dynamics by symmetric Hamiltonian splitting.
//!
//! fs-solid owns material geometry, analytical forces and reference inventory;
//! fs-mbd-material owns full-tensor principal lowering; fs-time owns free rotation.
//! This composition uses half momentum kicks, center-of-mass drift and the second
//! kick. Material nodes are not silently moved to their inertia centroids.
//! Smooth conservative fixed-step motion is second-order symplectic, NOT an exact
//! energy-momentum discrete-gradient method. Explicit wave stability and temporal
//! refinement remain necessary. Numerical energy defects are never physical heat.
//!
//! Eccentric link inventory is endpoint-lumped at original rotating offsets.
//! Free nodes and stationary frame clamps are supported. An eccentric free rotor
//! with constrained centerline translation needs a coupled constraint integrator
//! and refuses. No distributed cable inertia, moving supports, contact, material
//! axial/bending damping, damage or live solid-source replacement is inferred.
//! Explicit source-bound torsional viscosity uses its own passive loss flow.
//! Prescribed hydrostatic fluid refresh retains this trajectory and accounts for the
//! environmental parameter work at the original potential reference.

use std::sync::Arc;
use fs_exec::CancelGate;
use fs_solid::rod::material_field::SpatialRodError;
use fs_solid::rod::material_field::motion::{PreparedRodMotion, RodMotionPose,
    RodMotionResponse, RodNodalInertia, RodNodalSpatialInertia};
use fs_time::{LieStepError, So3, Vec3, rigid_body_momentum_step};

mod inertia {
//! Static centroid/principal-axis lowering, shared by both rod trajectories.
//! No new eigen or rigid-rotation solver: use the admitted multibody tensor,
//! existing principal-frame adapter and native fs-time free-rotation flow.
use fs_exec::CancelGate;
use fs_mbd::articulated::SpatialInertia;
use fs_mbd_material::principal::PrincipalFrame;
use fs_solid::rod::material_field::motion::{RodNodalInertia, RodNodalSpatialInertia};
use super::{RodMotionError, So3, a, v, rotate, finite, checkpoint, rigid_body_momentum_step};

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::RodNodalInertia {}
    impl Sealed for super::RodNodalSpatialInertia {}
}

/// Only the two source-owned reference-inertia views may enter this trajectory.
/// The trait does not permit arbitrary caller-authored replacement inventories.
pub trait RodMotionInertia: sealed::Sealed + Clone + core::fmt::Debug {
    #[doc(hidden)]
    fn motion_mass_kg(&self) -> f64;
    #[doc(hidden)]
    fn motion_centroid_body_m(&self) -> [f64; 3];
    #[doc(hidden)]
    fn motion_central_tensor_body_kg_m2(&self) -> [[f64; 3]; 3];
}
impl RodMotionInertia for RodNodalInertia {
    fn motion_mass_kg(&self) -> f64 { self.mass_kg() }
    fn motion_centroid_body_m(&self) -> [f64;3] { [0.0;3] }
    fn motion_central_tensor_body_kg_m2(&self) -> [[f64;3];3] {
        let [x,y,z]=self.principal_kg_m2();
        [[x,0.0,0.0],[0.0,y,0.0],[0.0,0.0,z]]
    }
}
impl RodMotionInertia for RodNodalSpatialInertia {
    fn motion_mass_kg(&self) -> f64 { self.mass_kg() }
    fn motion_centroid_body_m(&self) -> [f64;3] { self.center_of_mass_body_m() }
    fn motion_central_tensor_body_kg_m2(&self) -> [[f64;3];3] { self.inertia_com_body_kg_m2() }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct KineticFrame {
    pub(super) mass: f64,
    pub(super) center: [f64;3],
    pub(super) principal: [f64;3],
    pub(super) relative_tensor_residual: f64,
    body_from_principal: Option<So3>,
}
impl KineticFrame {
    pub(super) fn new(inertia: &impl RodMotionInertia, tolerance: f64, gate: &CancelGate)
        -> Result<Self,RodMotionError> {
        checkpoint(gate)?;
        if !tolerance.is_finite() || tolerance<=0.0 || tolerance>1e-8 {
            return Err(RodMotionError::Invalid("principal tensor residual tolerance must be in (0,1e-8]"));
        }
        let mass=inertia.motion_mass_kg();
        let center=inertia.motion_centroid_body_m();
        let tensor=inertia.motion_central_tensor_body_kg_m2();
        let mut matrix=So3::identity().matrix();
        matrix.m=std::array::from_fn(|k|tensor[k/3][k%3]);
        let spatial=SpatialInertia::new(mass,v(center),matrix).map_err(RodMotionError::SpatialInertia)?;
        // No basis selection or changed arithmetic for an already diagonal
        // tensor. This includes the original centered model exactly.
        let (principal,body_from_principal,relative_tensor_residual)=
            if (0..3).all(|i|(0..3).all(|j|i==j || tensor[i][j]==0.0)) {
                ([tensor[0][0],tensor[1][1],tensor[2][2]],None,0.0)
            } else {
                let lowered=PrincipalFrame::from_inertia(spatial,tolerance,||gate.is_requested())
                    .map_err(RodMotionError::PrincipalFrame)?;
                let p=lowered.properties().principal_inertia_body();
                ([p.x,p.y,p.z],Some(lowered.body_from_principal().rotation()),lowered.relative_tensor_residual())
            };
        checkpoint(gate)?;
        Ok(Self { mass,center,principal,body_from_principal,relative_tensor_residual })
    }
    pub(super) fn eccentric(self)->bool { self.center.iter().any(|c|*c!=0.0) }
    pub(super) fn omega(self,spin:[f64;3])->Result<[f64;3],RodMotionError> {
        let p=match self.body_from_principal { Some(q)=>rotate(q.inverse(),spin)?,None=>spin };
        let w=std::array::from_fn(|i|p[i]/self.principal[i]);
        let w=match self.body_from_principal { Some(q)=>rotate(q,w)?,None=>w };
        finite(&w)?;Ok(w)
    }
    pub(super) fn spin(self,omega:[f64;3])->Result<[f64;3],RodMotionError> {
        let w=match self.body_from_principal { Some(q)=>rotate(q.inverse(),omega)?,None=>omega };
        let p=std::array::from_fn(|i|w[i]*self.principal[i]);
        let p=match self.body_from_principal { Some(q)=>rotate(q,p)?,None=>p };
        finite(&p)?;Ok(p)
    }
    pub(super) fn drift(self,frame:So3,spin:[f64;3],dt:f64)
        ->Result<(So3,[f64;3]),RodMotionError> {
        match self.body_from_principal {
            None=>{
                let (r,p)=rigid_body_momentum_step(frame,v(spin),v(self.principal),dt)?;
                Ok((r,a(p)))
            }
            Some(q)=>{
                // An exact zero-spin pose must not drift through Q Q^-1 roundoff.
                if spin.iter().all(|p|*p==0.0) { return Ok((frame,spin)); }
                let r=frame.compose(q).map_err(fs_time::LieStepError::from)?;
                let p=rotate(q.inverse(),spin)?;
                let (r,p)=rigid_body_momentum_step(r,v(p),v(self.principal),dt)?;
                let frame=r.compose(q.inverse()).map_err(fs_time::LieStepError::from)?;
                Ok((frame,rotate(q,a(p))?))
            }
        }
    }
}

#[cfg(test)]
mod tests {
//! G0/G1/G3/G4: real material/rod owners, independent spatial mass equations.
use super::super::*;
use fs_evidence::ValidityDomain;
use fs_ga::{Mat3, So3Tangent};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink,TensionLinkAttachment};

fn section()->ConcentricShearSection {
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,value,dims) in [("density",10.0,Density::DIMS),("young_modulus",40.0,Pressure::DIMS),
        ("shear_modulus",16.0,Pressure::DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {source:"manufactured eccentric dynamics".into(),license:"test-only".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap();
    let material=resolve_material_state_point(&card,&QueryPoint::new(),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {outer_radius_m:0.2,material}])
        .unwrap().with_shear_modulus("shear_modulus",0.8).unwrap()
}
const OFFSET:[f64;3]=[0.05,0.08,-0.03];
fn specimen(offset:[f64;3])->SpatialMaterialRod {
    let source=section();
    let link=MaterialTensionLink::try_new(source.elastic().clone(),[
        TensionLinkAttachment {node:0,offset_body_m:offset},TensionLinkAttachment {node:1,offset_body_m:offset},
    ],1.1,1,0.2).unwrap();
    SpatialMaterialRod::try_straight(1.0,vec![source]).unwrap()
        .with_tension_links(vec![link],&CancelGate::new()).unwrap()
}
fn controls()->RodMotionControls {
    RodMotionControls {maximum_step_s:0.01,maximum_steps:10_000,energy_defect_limit_j:1e-6,
        maximum_rotation_rad:0.2,balance_origin_m:[0.0;3]}
}
fn near(a:f64,b:f64,t:f64) {assert!((a-b).abs()<=t*a.abs().max(b.abs()).max(1.0),"{a:e} != {b:e}");}
fn vector_near(a:[f64;3],b:[f64;3],t:f64) {for i in 0..3 {near(a[i],b[i],t);}}
fn norm(a:[f64;3])->f64 {a.into_iter().map(|x|x*x).sum::<f64>().sqrt()}
fn deformed()->SpatialRodMotion {
    let gate=CancelGate::new();
    let model=Arc::new(specimen(OFFSET).into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let mut pose=model.initial_pose();pose.positions_m[1][1]=0.01;
    pose.frames[1]=So3::exp(So3Tangent::new(v([0.03,-0.02,0.01]))).unwrap();
    SpatialRodMotion::new_spatial(model,pose,&[[-0.01,0.02,0.0],[0.03,0.0,-0.02]],
        &[[0.1,-0.05,0.02],[-0.02,0.04,0.08]],0.0,1e-12,&gate).unwrap()
}

#[test]
fn g0_initial_energy_and_momenta_equal_the_actual_constituent_inventory() {
    let gate=CancelGate::new();let source=specimen(OFFSET);
    let half=0.5*source.tension_links()[0].reference_mass_kg();
    let bare=SpatialMaterialRod::try_straight(1.0,vec![section()]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).unwrap();
    let model=Arc::new(source.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let pose=model.initial_pose();let positions=pose.positions_m.clone();
    let velocity=[[0.2,-0.1,0.3],[-0.2,0.1,0.05]];let omega=[[0.3,0.2,-0.1],[-0.1,0.2,0.15]];
    let motion=SpatialRodMotion::new_spatial(model,pose,&velocity,&omega,0.0,1e-12,&gate).unwrap();
    let mut energy=0.0;let mut total_p=[0.0;3];let mut total_l=[0.0;3];
    for node in 0..2 {
        let b=bare.inertia()[node];let point_velocity=add(velocity[node],cross(omega[node],OFFSET));
        let point_p=point_velocity.map(|x|half*x);let base_p=velocity[node].map(|x|b.mass_kg()*x);
        let base_spin=std::array::from_fn(|i|b.principal_kg_m2()[i]*omega[node][i]);
        let momentum=add(base_p,point_p);total_p=add(total_p,momentum);
        total_l=add(total_l,add(add(cross(positions[node],base_p),base_spin),cross(add(positions[node],OFFSET),point_p)));
        for c in 0..3 {energy+=0.5*b.mass_kg()*velocity[node][c]*velocity[node][c]
            +0.5*base_spin[c]*omega[node][c]+0.5*half*point_velocity[c]*point_velocity[c];}
        vector_near(motion.linear_momenta_kg_m_s()[node],momentum,1e-13);
        vector_near(motion.velocities_m_s().unwrap()[node],velocity[node],1e-13);
        vector_near(motion.angular_velocities_body_rad_s().unwrap()[node],omega[node],1e-13);
    }
    let observed=motion.observe([0.0;3]).unwrap();near(observed.kinetic_energy_j,energy,1e-13);
    vector_near(observed.linear_momentum_kg_m_s,total_p,1e-13);
    vector_near(observed.angular_momentum_kg_m2_s,total_l,1e-13);
    assert!(motion.maximum_inertia_tensor_residual()<=1e-12);
}

#[test]
fn g0_material_node_kick_matches_independent_six_by_six_spatial_mass_solve_and_work() {
    let gate=CancelGate::new();let model=Arc::new(specimen(OFFSET).into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let mut pose=model.initial_pose();let rotation=So3::exp(So3Tangent::new(v([0.4,-0.2,0.3]))).unwrap();
    for (x,r) in pose.positions_m.iter_mut().zip(&mut pose.frames) {*x=add(rotate(rotation,*x).unwrap(),[0.7,-0.2,0.1]);*r=rotation;}
    let mut motion=SpatialRodMotion::new_spatial(model,pose,&[[0.1,0.2,-0.1];2],&[[0.2,-0.1,0.3];2],0.0,1e-12,&gate).unwrap();
    let before=motion.observe([0.0;3]).unwrap();let old_v=motion.velocities_m_s().unwrap();let old_w=motion.angular_velocities_body_rad_s().unwrap();
    let loads=[RodNodalWrench {force_n:[0.3,-0.1,0.2],torque_nm:[0.02,0.01,-0.03]};2];let h=0.003;
    let mut expected=Vec::new();
    for node in 0..2 {
        let i=motion.model().inertia()[node];let tensor=i.inertia_com_body_kg_m2();
        let matrix=fs_mbd::articulated::SpatialInertia::new(i.mass_kg(),v(i.center_of_mass_body_m()),
            Mat3 {m:std::array::from_fn(|k|tensor[k/3][k%3])}).unwrap().matrix().unwrap();
        let f=rotate(rotation.inverse(),add(motion.response.force_n[node],loads[node].force_n)).unwrap();
        let t=add(motion.response.torque_body_nm[node],rotate(rotation.inverse(),loads[node].torque_nm).unwrap());
        let mut change=[t[0],t[1],t[2],f[0],f[1],f[2]];
        fs_la::factor::lu(&matrix.m,6).unwrap().solve(&mut change);
        let body_v=rotate(rotation.inverse(),old_v[node]).unwrap();
        expected.push((std::array::from_fn(|c|old_w[node][c]+h*change[c]),
            std::array::from_fn(|c|body_v[c]+h*change[c+3])));
    }
    let mut accounts=KickAccounts::new(2);motion.kick(h,StageLoads::World(&loads),[0.0;3],&mut accounts,&gate).unwrap();
    for (node,(w,bv)) in expected.into_iter().enumerate() {
        vector_near(motion.angular_velocities_body_rad_s().unwrap()[node],w,2e-12);
        vector_near(rotate(rotation.inverse(),motion.velocities_m_s().unwrap()[node]).unwrap(),bv,2e-12);
    }
    let after=motion.observe([0.0;3]).unwrap();near(after.kinetic_energy_j-before.kinetic_energy_j,accounts.work,2e-12);
    vector_near(sub(after.angular_momentum_kg_m2_s,before.angular_momentum_kg_m2_s),accounts.angular,2e-12);
}

#[test]
fn g0_fixed_frame_translation_has_the_required_eccentric_reaction_couple() {
    let gate=CancelGate::new();let rod=specimen(OFFSET).with_supports(&[RodSupport {node:1,
        position_m:[None;3],frame:Some(So3::identity())}],&gate).unwrap();
    let model=Arc::new(rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());let pose=model.initial_pose();
    let mut motion=SpatialRodMotion::new_spatial(model,pose,&[[0.0;3];2],&[[0.0;3];2],0.0,1e-12,&gate).unwrap();
    let c=motion.model().inertia()[1].center_of_mass_body_m();let mass=motion.model().inertia()[1].mass_kg();
    let load=RodNodalWrench {force_n:[0.2,-0.1,0.3],torque_nm:[0.01,0.02,-0.01]};
    let h=0.002;let mut accounts=KickAccounts::new(2);
    motion.kick(h,StageLoads::World(&[RodNodalWrench::default(),load]),[0.0;3],&mut accounts,&gate).unwrap();
    vector_near(accounts.supports[1].angular_nm_s,add(load.torque_nm.map(|x|-h*x),cross(c,load.force_n.map(|x|h*x))),1e-13);
    assert_eq!(motion.body_spins_kg_m2_s()[1],[0.0;3]);
    near(accounts.work,0.5*h*h*load.force_n.iter().map(|f|f*f).sum::<f64>()/mass,1e-13);
    vector_near(motion.observe([0.0;3]).unwrap().angular_momentum_kg_m2_s,accounts.angular,1e-13);
}

#[test]
fn g1_eccentric_free_fall_has_no_spurious_gravity_torque_or_double_counted_work() {
    let gate=CancelGate::new();let g=[0.0,-9.0,1.0];let rod=specimen(OFFSET).with_gravity(g,&gate).unwrap();
    let model=Arc::new(rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());let pose=model.initial_pose();let original=pose.clone();
    let mut motion=SpatialRodMotion::new_spatial(model,pose,&[[0.0;3];2],&[[0.0;3];2],0.0,1e-12,&gate).unwrap();
    for _ in 0..20 {
        let (next,report)=motion.step_candidate(0.001,&[RodNodalWrench::default();2],controls(),&gate).unwrap();
        assert_eq!(report.source_work_quadrature_j,0.0);
        assert!(norm(report.angular_balance_residual_nm_s)<1e-11);
        motion=next;
    }
    for node in 0..2 {
        vector_near(motion.pose().positions_m[node],add(original.positions_m[node],g.map(|a|0.5*a*0.02*0.02)),1e-11);
        vector_near(motion.angular_velocities_body_rad_s().unwrap()[node],[0.0;3],1e-10);
    }
}

#[test]
fn g1_deforming_eccentric_free_flight_refines_and_preserves_both_momenta() {
    let gate=CancelGate::new();let initial=deformed();let before=initial.observe([0.0;3]).unwrap();
    let energy=before.kinetic_energy_j+before.potential_energy_j;let mut errors=Vec::new();
    for count in [20,40,80] {
        let mut state=initial.clone();
        for _ in 0..count {
            let (next,receipt)=state.step_candidate(0.02/count as f64,&[RodNodalWrench::default();2],controls(),&gate).unwrap();
            assert!(norm(receipt.linear_balance_residual_n_s)<1e-11);
            assert!(norm(receipt.angular_balance_residual_nm_s)<1e-11);state=next;
        }
        let end=state.observe([0.0;3]).unwrap();vector_near(end.linear_momentum_kg_m_s,before.linear_momentum_kg_m_s,1e-10);
        vector_near(end.angular_momentum_kg_m2_s,before.angular_momentum_kg_m2_s,1e-10);
        errors.push((end.kinetic_energy_j+end.potential_energy_j-energy).abs());
    }
    assert!(errors[1]<0.35*errors[0] && errors[2]<0.35*errors[1],"{errors:?}");
}

#[test]
fn g3_centered_specialization_and_spatial_zero_offset_are_bitwise_identical() {
    let gate=CancelGate::new();let source=specimen([0.0;3]);
    let centered=Arc::new(source.clone().into_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let spatial=Arc::new(source.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let mut pose=centered.initial_pose();pose.positions_m[1][1]=0.01;
    let velocity=[[0.01,0.02,-0.01];2];let omega=[[0.1,-0.03,0.02];2];
    let mut a=RodMotion::new(centered,pose.clone(),&velocity,&omega,0.0,&gate).unwrap();
    let mut b=SpatialRodMotion::new_spatial(spatial,pose,&velocity,&omega,0.0,1e-12,&gate).unwrap();
    for _ in 0..10 {
        let (x,rx)=a.step_candidate(0.001,&[RodNodalWrench::default();2],controls(),&gate).unwrap();
        let (y,ry)=b.step_candidate(0.001,&[RodNodalWrench::default();2],controls(),&gate).unwrap();
        assert_eq!(x.pose(),y.pose());assert_eq!(x.linear_momenta_kg_m_s(),y.linear_momenta_kg_m_s());
        assert_eq!(x.body_spins_kg_m2_s(),y.body_spins_kg_m2_s());assert_eq!(rx.after,ry.after);a=x;b=y;
    }
}

#[test]
fn g4_unsupported_eccentric_pins_refuse_instead_of_projecting_the_wrong_inertia() {
    let gate=CancelGate::new();let rod=specimen(OFFSET).with_supports(&[RodSupport {node:1,
        position_m:[None,Some(0.0),None],frame:None}],&gate).unwrap();
    let model=Arc::new(rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap());let pose=model.initial_pose();
    assert!(matches!(SpatialRodMotion::new_spatial(model.clone(),pose.clone(),&[[0.0;3];2],&[[0.0;3];2],
        0.0,1e-12,&gate),Err(RodMotionError::EccentricPositionSupport {node:1})));
    assert_eq!(model.initial_pose(),pose);
    for tolerance in [0.0,-1.0,f64::NAN,1e-7] {
        assert!(SpatialRodMotion::new_spatial(model.clone(),pose.clone(),&[[0.0;3];2],&[[0.0;3];2],
            0.0,tolerance,&gate).is_err());
    }
}

#[test]
fn g4_rejection_cancellation_and_cloned_continuation_preserve_the_complete_motion() {
    let gate=CancelGate::new();let initial=deformed();let original=initial.observe([0.0;3]).unwrap();
    let loads=[RodNodalWrench::default();2];let mut tight=controls();tight.energy_defect_limit_j=0.0;
    assert!(matches!(initial.step_candidate(0.001,&loads,tight,&gate),Err(RodMotionError::Energy {..})));
    assert_eq!(initial.observe([0.0;3]).unwrap(),original);assert_eq!(initial.accepted_steps(),0);
    let (mut a,_)=initial.step_candidate(0.001,&loads,controls(),&gate).unwrap();let mut b=a.clone();
    for _ in 0..5 {a=a.step_candidate(0.001,&loads,controls(),&gate).unwrap().0;b=b.step_candidate(0.001,&loads,controls(),&gate).unwrap().0;}
    assert_eq!(a.pose(),b.pose());assert_eq!(a.linear_momenta_kg_m_s(),b.linear_momenta_kg_m_s());assert_eq!(a.body_spins_kg_m2_s(),b.body_spins_kg_m2_s());
    gate.request();assert!(matches!(initial.step_candidate(0.001,&loads,controls(),&gate),Err(RodMotionError::Solid(SpatialRodError::Cancelled))));
    assert_eq!(initial.observe([0.0;3]).unwrap(),original);
}

}

}
pub use inertia::RodMotionInertia;
use inertia::KineticFrame;

/// Continuing material-fluid snapshots on the same native trajectory.
pub mod environment;
pub use environment::RodFluidRefresh;

/// Prescribed section temperatures on the original continuing trajectory.
pub mod thermal;
pub use thermal::RodThermalRefresh;

/// Reciprocal source-bound section thermodynamics on the finite heat timeline.
pub mod enthalpy;

/// Held spatial forces and material-following couples at rotating rod nodes.
pub mod loads;
pub use loads::{PreparedRodLoads, RodPointForce};
use loads::StageLoads;

/// Source-bound intrinsic circular-shaft viscosity on the native motion owner.
pub mod torsional_loss;
pub use torsional_loss::{PreparedRodTorsionalLoss, RodTorsionalHeat};

/// Instantaneous terminal and attachment impulses on the same trajectory.
pub mod impulse;

/// A force AT the material node and a pure couple held in the WORLD frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RodNodalWrench {
    pub force_n: [f64; 3],
    pub torque_nm: [f64; 3],
}
#[derive(Clone, Copy, Debug)]
pub struct RodMotionControls {
    /// Explicit upper step bound, not an inferred stability limit.
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    /// Bound on |delta kinetic + delta potential + physical heat - source work| [J].
    pub energy_defect_limit_j: f64,
    /// Upper bound on each free rotor's swept angle, in (0, pi/2].
    pub maximum_rotation_rad: f64,
    /// Fixed origin for orbital angular momentum and applied-moment accounting.
    pub balance_origin_m: [f64; 3],
}
#[derive(Clone, Debug)]
pub enum RodMotionError {
    Solid(SpatialRodError),
    Rotation(LieStepError),
    SpatialInertia(fs_mbd::articulated::ArticulatedError),
    PrincipalFrame(fs_mbd_material::principal::PrincipalFrameError),
    SpatialLoads(fs_mbd::spatial_loads::SpatialLoadError),
    /// A noncentral pin/slider would couple the constrained translational drift
    /// to rotation. It cannot be implemented by independent momentum projection.
    EccentricPositionSupport { node: usize },
    Invalid(&'static str),
    Energy { defect_j: f64, limit_j: f64 },
    /// A requested impulse was not resolved by the stored momentum update.
    ImpulseBalance { linear_defect_n_s: f64, angular_defect_nm_s: f64, work_defect_j: f64 },
}
impl core::fmt::Display for RodMotionError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        match self {
            Self::Solid(e)=>write!(f,"rod motion: {e}"),
            Self::Rotation(e)=>write!(f,"rod motion: {e}"),
            Self::SpatialInertia(e)=>write!(f,"rod reference inertia: {e}"),
            Self::PrincipalFrame(e)=>write!(f,"rod principal inertia: {e}"),
            Self::SpatialLoads(e)=>write!(f,"rod applied load: {e}"),
            Self::EccentricPositionSupport {node}=>write!(f,
                "rod node {node}: eccentric free rotation with a position support requires coupled constraint dynamics"),
            Self::Invalid(s)=>write!(f,"rod motion: {s}"),
            Self::Energy {defect_j,limit_j}=>write!(f,
                "rod motion energy defect {defect_j:e} J exceeds {limit_j:e} J"),
            Self::ImpulseBalance {linear_defect_n_s,angular_defect_nm_s,work_defect_j}=>write!(f,
                "rod impulse transfer exceeds its tolerances: linear {linear_defect_n_s:e} N s, angular {angular_defect_nm_s:e} N m s, work {work_defect_j:e} J"),
        }
    }
}
impl std::error::Error for RodMotionError {}
impl From<SpatialRodError> for RodMotionError { fn from(e:SpatialRodError)->Self {Self::Solid(e)} }
impl From<LieStepError> for RodMotionError { fn from(e:LieStepError)->Self {Self::Rotation(e)} }

#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionObservables {
    pub kinetic_energy_j: f64,
    /// Native elastic/link storage plus gravity and hydrostatic potentials.
    pub potential_energy_j: f64,
    pub linear_momentum_kg_m_s: [f64; 3],
    /// Centroid orbital momentum about the requested origin PLUS central spin.
    pub angular_momentum_kg_m2_s: [f64; 3],
}
#[derive(Clone, Debug, PartialEq)]
pub struct RodSupportImpulse {
    /// Impulse ON the rod, world coordinates [N s].
    pub linear_n_s: [f64; 3],
    /// Pure support couple ABOUT THE MATERIAL NODE, world coordinates [N m s].
    pub angular_nm_s: [f64; 3],
}
#[derive(Clone, Debug)]
pub struct RodMotionStep {
    pub before: RodMotionObservables,
    pub after: RodMotionObservables,
    /// Symmetric kick quadrature of applied material-node wrench power.
    /// Gravity and buoyancy remain in potential energy, not counted twice as work.
    pub source_work_quadrature_j: f64,
    /// Delta mechanical energy plus physical heat minus applied-work quadrature.
    pub integration_defect_j: f64,
    /// Nonnegative intrinsic shear heat in original segment order, when enabled.
    /// The numerical defect accounts for this heat separately from source work.
    pub torsional_heat: Option<RodTorsionalHeat>,
    /// Applied loads, gravity/buoyancy wrenches and stationary-support impulses.
    pub external_linear_impulse_n_s: [f64; 3],
    pub external_angular_impulse_nm_s: [f64; 3],
    pub linear_balance_residual_n_s: [f64; 3],
    pub angular_balance_residual_nm_s: [f64; 3],
    pub support_impulses: Vec<RodSupportImpulse>,
}

/// Accepted trajectory. Both inertia views use this SAME time/accounting owner.
/// Cloning shares immutable sources and principal frames and copies motion only.
#[derive(Clone, Debug)]
pub struct RodMotion<I = RodNodalInertia> {
    model: Arc<PreparedRodMotion<I>>,
    kinetics: Arc<[KineticFrame]>,
    pose: RodMotionPose,
    linear: Vec<[f64; 3]>,
    /// Central angular momentum expressed in the ORIGINAL material frame.
    spin: Vec<[f64; 3]>,
    response: RodMotionResponse,
    time_s: f64,
    steps: usize,
}
/// Full eccentric-inertia specialization; material-coordinate API is unchanged.
pub type SpatialRodMotion = RodMotion<RodNodalSpatialInertia>;

impl RodMotion {
    /// Initialize the original centered lane; velocities refer to material nodes.
    pub fn new(model:Arc<PreparedRodMotion>,pose:RodMotionPose,
        velocity_m_s:&[[f64;3]],angular_velocity_body_rad_s:&[[f64;3]],time_s:f64,
        gate:&CancelGate)->Result<Self,RodMotionError> {
        Self::initialize(model,pose,velocity_m_s,angular_velocity_body_rad_s,time_s,1e-12,gate)
    }
}
impl SpatialRodMotion {
    /// Initialize eccentric dynamics from material-node velocities, NOT centroid
    /// velocities. Lower the full tensor once under an explicit residual budget.
    /// Original material frames never become arbitrary principal eigenvectors.
    #[allow(clippy::too_many_arguments)]
    pub fn new_spatial(model:Arc<PreparedRodMotion<RodNodalSpatialInertia>>,pose:RodMotionPose,
        velocity_m_s:&[[f64;3]],angular_velocity_body_rad_s:&[[f64;3]],time_s:f64,
        relative_tensor_tolerance:f64,gate:&CancelGate)->Result<Self,RodMotionError> {
        Self::initialize(model,pose,velocity_m_s,angular_velocity_body_rad_s,time_s,
            relative_tensor_tolerance,gate)
    }
}
impl<I:RodMotionInertia> RodMotion<I> {
    #[allow(clippy::too_many_arguments)]
    fn initialize(model:Arc<PreparedRodMotion<I>>,pose:RodMotionPose,
        velocity_m_s:&[[f64;3]],angular_velocity_body_rad_s:&[[f64;3]],time_s:f64,
        tolerance:f64,gate:&CancelGate)->Result<Self,RodMotionError> {
        checkpoint(gate)?;
        let n=model.inertia().len();
        if !time_s.is_finite() || velocity_m_s.len()!=n || angular_velocity_body_rad_s.len()!=n
            || velocity_m_s.iter().flatten().chain(angular_velocity_body_rad_s.iter().flatten()).any(|v|!v.is_finite()) {
            return Err(RodMotionError::Invalid("matching finite velocities and initial time required"));
        }
        let response=model.evaluate(&pose,gate)?;
        let mut kinetics=Vec::with_capacity(n);
        let mut linear=vec![[0.0;3];n];
        let mut spin=vec![[0.0;3];n];
        for node in 0..n {
            checkpoint(gate)?;
            let k=KineticFrame::new(&model.inertia()[node],tolerance,gate)?;
            if k.eccentric() && !model.fixed_frames()[node]
                && model.fixed_positions()[node].iter().any(|fixed|*fixed) {
                return Err(RodMotionError::EccentricPositionSupport {node});
            }
            for c in 0..3 {
                if (model.fixed_positions()[node][c] && velocity_m_s[node][c]!=0.0)
                    || (model.fixed_frames()[node] && angular_velocity_body_rad_s[node][c]!=0.0) {
                    return Err(RodMotionError::Invalid("stationary supports require zero supported velocity"));
                }
            }
            let omega=angular_velocity_body_rad_s[node];
            let velocity=if k.eccentric() {
                add(velocity_m_s[node],rotate(pose.frames[node],cross(omega,k.center))?)
            } else {velocity_m_s[node]};
            linear[node]=velocity.map(|value|k.mass*value);
            spin[node]=k.spin(omega)?;
            kinetics.push(k);
        }
        let state=Self {model,kinetics:Arc::from(kinetics),pose,linear,spin,response,time_s,steps:0};
        state.observe([0.0;3])?;
        checkpoint(gate)?;Ok(state)
    }
    pub fn model(&self)->&PreparedRodMotion<I> {&self.model}
    pub fn pose(&self)->&RodMotionPose {&self.pose}
    /// Total nodal linear momenta (mass times centroid velocity).
    pub fn linear_momenta_kg_m_s(&self)->&[[f64;3]] {&self.linear}
    /// Spin about each mass centroid in the original material axes.
    pub fn body_spins_kg_m2_s(&self)->&[[f64;3]] {&self.spin}
    pub const fn time_s(&self)->f64 {self.time_s}
    pub const fn accepted_steps(&self)->usize {self.steps}
    pub fn maximum_inertia_tensor_residual(&self)->f64 {
        self.kinetics.iter().map(|k|k.relative_tensor_residual).fold(0.0,f64::max)
    }
    pub fn angular_velocities_body_rad_s(&self)->Result<Vec<[f64;3]>,RodMotionError> {
        self.kinetics.iter().zip(&self.spin).map(|(k,p)|k.omega(*p)).collect()
    }
    /// Material-node velocity, including the rotation of its offset from COM.
    pub fn velocities_m_s(&self)->Result<Vec<[f64;3]>,RodMotionError> {
        self.kinetics.iter().enumerate().map(|(node,k)| {
            let velocity=self.linear[node].map(|p|p/k.mass);
            let velocity=if k.eccentric() {sub(velocity,rotate(self.pose.frames[node],cross(k.omega(self.spin[node])?,k.center))?)}
                else {velocity};
            finite(&velocity)?;Ok(velocity)
        }).collect()
    }
    pub fn observe(&self,origin_m:[f64;3])->Result<RodMotionObservables,RodMotionError> {
        finite(&origin_m)?;
        let mut kinetic=0.0;let mut linear=[0.0;3];let mut angular=[0.0;3];
        for (node,k) in self.kinetics.iter().enumerate() {
            let p=self.linear[node];let spin=self.spin[node];let omega=k.omega(spin)?;
            finite(&p)?;finite(&spin)?;
            for c in 0..3 {kinetic+=0.5*p[c]*(p[c]/k.mass)+0.5*spin[c]*omega[c];linear[c]+=p[c];}
            let mut arm=sub(self.pose.positions_m[node],origin_m);
            if k.eccentric() {arm=add(arm,rotate(self.pose.frames[node],k.center)?);}
            let orbital=cross(arm,p);let world=rotate(self.pose.frames[node],spin)?;
            for c in 0..3 {angular[c]+=orbital[c]+world[c];}
        }
        finite(&[kinetic,self.response.potential_energy_j])?;finite(&linear)?;finite(&angular)?;
        if kinetic<0.0 {return Err(RodMotionError::Invalid("negative nodal kinetic storage"));}
        Ok(RodMotionObservables {kinetic_energy_j:kinetic,potential_energy_j:self.response.potential_energy_j,
            linear_momentum_kg_m_s:linear,angular_momentum_kg_m2_s:angular})
    }
    /// One atomic candidate. One fresh analytical force assembly; no dense
    /// tangent or new eigendecomposition in the step. Sources remain immutable.
    pub fn step_candidate(&self,dt:f64,loads:&[RodNodalWrench],controls:RodMotionControls,
        gate:&CancelGate)->Result<(Self,RodMotionStep),RodMotionError> {
        self.step_with_loads(dt,StageLoads::World(loads),controls,gate)
    }
    /// Advance under world-held material-point forces and body/world free
    /// couples. Each kick evaluates the actual material frame, including the
    /// endpoint after free rotation. Reference points remain material nodes;
    /// the existing kick transports their moments to eccentric mass centroids.
    ///
    /// The loads are constant in their named frames over this interval. Native
    /// work remains a symmetric kick quadrature, not an exact path integral.
    pub fn step_spatial_candidate(&self,dt:f64,loads:&PreparedRodLoads,controls:RodMotionControls,
        gate:&CancelGate)->Result<(Self,RodMotionStep),RodMotionError> {
        self.step_with_loads(dt,StageLoads::Spatial(loads),controls,gate)
    }
    fn step_with_loads(&self,dt:f64,loads:StageLoads<'_>,controls:RodMotionControls,
        gate:&CancelGate)->Result<(Self,RodMotionStep),RodMotionError> {
        self.step_with_torsion_loads(dt,loads,None,controls,gate)
    }
    fn step_with_torsion_loads(&self,dt:f64,loads:StageLoads<'_>,loss:Option<&PreparedRodTorsionalLoss>,controls:RodMotionControls,
        gate:&CancelGate)->Result<(Self,RodMotionStep),RodMotionError> {
        checkpoint(gate)?;
        if !dt.is_finite() || dt<=0.0 || 0.5*dt==0.0
            || !controls.maximum_step_s.is_finite() || controls.maximum_step_s<=0.0 || dt>controls.maximum_step_s
            || !controls.energy_defect_limit_j.is_finite() || controls.energy_defect_limit_j<0.0
            || !controls.maximum_rotation_rad.is_finite() || controls.maximum_rotation_rad<=0.0
            || controls.maximum_rotation_rad>core::f64::consts::FRAC_PI_2
            || self.steps>=controls.maximum_steps || loads.len()!=self.kinetics.len() {
            return Err(RodMotionError::Invalid("invalid rod interval, work gate or step/rotation budget"));
        }
        loads.validate()?;
        if let Some(law)=loss {law.admit(&self.model,gate)?;}
        let next_time=self.time_s+dt;
        if !next_time.is_finite() || next_time<=self.time_s {
            return Err(RodMotionError::Invalid("rod clock cannot represent the requested advance"));
        }
        let before=self.observe(controls.balance_origin_m)?;
        let mut candidate=self.clone();
        let mut accounts=KickAccounts::new(loads.len());
        let mut torsional_heat=loss.map(PreparedRodTorsionalLoss::heat);
        if let (Some(law),Some(heat))=(loss,torsional_heat.as_mut()) {
            law.flow(&mut candidate,0.5*dt,&mut accounts,heat,gate)?;
        }
        candidate.kick(0.5*dt,loads,controls.balance_origin_m,&mut accounts,gate)?;
        for node in 0..loads.len() {
            checkpoint(gate)?;
            let k=candidate.kinetics[node];let old_frame=candidate.pose.frames[node];
            if !candidate.model.fixed_frames()[node] {
                let pi=candidate.spin[node];
                let swept=dt*(fs_math::det::hypot(fs_math::det::hypot(pi[0],pi[1]),pi[2])
                    /k.principal.into_iter().fold(f64::INFINITY,f64::min));
                if !swept.is_finite() || swept>controls.maximum_rotation_rad {
                    return Err(RodMotionError::Invalid("rotor swept-angle budget exceeded; reduce the interval"));
                }
                let (frame,spin)=k.drift(old_frame,pi,dt)?;
                candidate.pose.frames[node]=frame;candidate.spin[node]=spin;
            }
            let correction=if k.eccentric() && !candidate.model.fixed_frames()[node] {
                sub(rotate(old_frame,k.center)?,rotate(candidate.pose.frames[node],k.center)?)
            } else {[0.0;3]};
            for c in 0..3 {
                if !candidate.model.fixed_positions()[node][c] {
                    let drift=dt*(candidate.linear[node][c]/k.mass);
                    let shift=if k.eccentric() {drift+correction[c]} else {drift};
                    let old=candidate.pose.positions_m[node][c];let next=old+shift;
                    if !shift.is_finite() || !next.is_finite() || (shift!=0.0 && next==old) {
                        return Err(RodMotionError::Invalid("rod translation step is unrepresentable at this world position"));
                    }
                    candidate.pose.positions_m[node][c]=next;
                }
            }
        }
        candidate.response=candidate.model.evaluate(&candidate.pose,gate)?;
        candidate.kick(0.5*dt,loads,controls.balance_origin_m,&mut accounts,gate)?;
        if let (Some(law),Some(heat))=(loss,torsional_heat.as_mut()) {
            law.flow(&mut candidate,0.5*dt,&mut accounts,heat,gate)?;
        }
        let after=candidate.observe(controls.balance_origin_m)?;
        let mut defect=(after.kinetic_energy_j-before.kinetic_energy_j)
            +(after.potential_energy_j-before.potential_energy_j)-accounts.work;
        if let Some(heat)=&torsional_heat {defect+=heat.total_heat_j;}
        finite(&[accounts.work,defect])?;
        if defect.abs()>controls.energy_defect_limit_j {
            return Err(RodMotionError::Energy {defect_j:defect,limit_j:controls.energy_defect_limit_j});
        }
        finite(&accounts.linear)?;finite(&accounts.angular)?;
        let linear_residual=sub(sub(after.linear_momentum_kg_m_s,before.linear_momentum_kg_m_s),accounts.linear);
        let angular_residual=sub(sub(after.angular_momentum_kg_m2_s,before.angular_momentum_kg_m2_s),accounts.angular);
        finite(&linear_residual)?;finite(&angular_residual)?;
        candidate.time_s=next_time;
        candidate.steps=self.steps.checked_add(1).ok_or(RodMotionError::Invalid("rod step count overflow"))?;
        checkpoint(gate)?;
        Ok((candidate,RodMotionStep {before,after,source_work_quadrature_j:accounts.work,
            integration_defect_j:defect,torsional_heat,external_linear_impulse_n_s:accounts.linear,
            external_angular_impulse_nm_s:accounts.angular,linear_balance_residual_n_s:linear_residual,
            angular_balance_residual_nm_s:angular_residual,support_impulses:accounts.supports}))
    }
    fn kick(&mut self,h:f64,loads:StageLoads<'_>,origin:[f64;3],accounts:&mut KickAccounts,
        gate:&CancelGate)->Result<(),RodMotionError> {
        let body_loads=self.model.body_loads(&self.pose,gate)?;
        for node in 0..loads.len() {
            checkpoint(gate)?;
            let frame=self.pose.frames[node];let k=self.kinetics[node];
            let (load,applied_body)=loads.stage_at(node,frame)?;
            let total=add(self.response.force_n[node],load.force_n);
            let total_torque=add(self.response.torque_body_nm[node],applied_body);
            let central_torque=if k.eccentric() {sub(total_torque,cross(k.center,rotate(frame.inverse(),total)?))}
                else {total_torque};
            let applied_central=if k.eccentric() {sub(applied_body,cross(k.center,rotate(frame.inverse(),load.force_n)?))}
                else {applied_body};
            let old_p=self.linear[node];let old_omega=k.omega(self.spin[node])?;
            let mut support_f=[0.0;3];let mut support_t=[0.0;3];
            for c in 0..3 {
                if self.model.fixed_positions()[node][c] {support_f[c]=-h*total[c];}
                else {self.linear[node][c]=h.mul_add(total[c],old_p[c]);}
                if self.model.fixed_frames()[node] {support_t[c]=-h*central_torque[c];}
                else {self.spin[node][c]=h.mul_add(central_torque[c],self.spin[node][c]);}
            }
            // A fixed material frame may translate. Its balancing pure couple
            // includes the moment of both applied and constrained linear kicks.
            if k.eccentric() && self.model.fixed_frames()[node] {
                support_t=add(support_t,cross(k.center,rotate(frame.inverse(),support_f)?));
            }
            let omega=k.omega(self.spin[node])?;
            for c in 0..3 {
                accounts.work+=h*load.force_n[c]*(0.5*(old_p[c]/k.mass)+0.5*(self.linear[node][c]/k.mass));
                accounts.work+=h*applied_central[c]*(0.5*old_omega[c]+0.5*omega[c]);
            }
            let support_t=rotate(frame,support_t)?;
            let body_force=body_loads.force_n[node];
            let body_couple=rotate(frame,body_loads.torque_body_nm[node])?;
            let impulse=std::array::from_fn(|c|h*(load.force_n[c]+body_force[c])+support_f[c]);
            let orbital=cross(sub(self.pose.positions_m[node],origin),impulse);
            for c in 0..3 {
                accounts.linear[c]+=impulse[c];
                accounts.angular[c]+=orbital[c]+h*(load.torque_nm[c]+body_couple[c])+support_t[c];
                accounts.supports[node].linear_n_s[c]+=support_f[c];
                accounts.supports[node].angular_nm_s[c]+=support_t[c];
            }
            finite(&self.linear[node])?;finite(&self.spin[node])?;
            finite(&accounts.supports[node].linear_n_s)?;finite(&accounts.supports[node].angular_nm_s)?;
        }
        Ok(())
    }
}
struct KickAccounts {work:f64,linear:[f64;3],angular:[f64;3],supports:Vec<RodSupportImpulse>}
impl KickAccounts {
    fn new(n:usize)->Self {Self {work:0.0,linear:[0.0;3],angular:[0.0;3],
        supports:vec![RodSupportImpulse {linear_n_s:[0.0;3],angular_nm_s:[0.0;3]};n]}}
}
fn checkpoint(gate:&CancelGate)->Result<(),RodMotionError> {
    if gate.is_requested() {Err(SpatialRodError::Cancelled.into())} else {Ok(())}
}
fn finite(values:&[f64])->Result<(),RodMotionError> {
    if values.iter().any(|v|!v.is_finite()) {Err(RodMotionError::Invalid("nonfinite rod motion arithmetic"))} else {Ok(())}
}
fn v(x:[f64;3])->Vec3 {Vec3::new(x[0],x[1],x[2])}
fn a(v:Vec3)->[f64;3] {[v.x,v.y,v.z]}
fn rotate(r:So3,x:[f64;3])->Result<[f64;3],RodMotionError> {
    Ok(a(r.rotate(v(x)).map_err(LieStepError::from)?))
}
fn add(a:[f64;3],b:[f64;3])->[f64;3] {std::array::from_fn(|i|a[i]+b[i])}
fn sub(a:[f64;3],b:[f64;3])->[f64;3] {std::array::from_fn(|i|a[i]-b[i])}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]}

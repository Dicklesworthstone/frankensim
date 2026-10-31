//! G0/G1/G3/G4: real native rod sources and the existing free-rotation owner.
#![cfg(feature = "rod-dynamics")]
use std::sync::Arc;
use fs_couple::rod_motion::{RodMotion,RodMotionControls,RodMotionError,RodNodalWrench};
use fs_exec::CancelGate;
use fs_evidence::ValidityDomain;
use fs_ga::{So3,So3Tangent,Vec3};
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer,ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection,ScalarAdmissibility,
    ScalarPropertyRequirement,resolve_material_state_point};
use fs_qty::{Density,Pressure};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::{PreparedRodMotion,RodRootBoundary};
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink,TensionLinkAttachment};

#[path = "rod_motion/environment.rs"]
mod environment;

#[path = "rod_motion/thermal.rs"]
mod thermal;

#[path = "rod_motion/enthalpy.rs"]
mod enthalpy;

#[path = "rod_motion/impulse.rs"]
mod impulse;

fn section(scale:f64)->ConcentricShearSection {
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),
        ("young_modulus",1e4*scale,Pressure::DIMS),("shear_modulus",4e3*scale,Pressure::DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},validity:ValidityDomain::unconstrained(),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic rod dynamics".into(),
                license:"test-only".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap();
    let material=resolve_material_state_point(&card,&QueryPoint::new(),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {outer_radius_m:0.05,material}])
        .unwrap().with_shear_modulus("shear_modulus",0.8).unwrap()
}
fn controls(dt:f64)->RodMotionControls {RodMotionControls {maximum_step_s:dt,maximum_steps:10_000,
    energy_defect_limit_j:1e-4,maximum_rotation_rad:0.25,balance_origin_m:[0.0;3]}}
fn axial_model(link:bool)->Arc<PreparedRodMotion> {
    let gate=CancelGate::new();let source=section(1.0);
    let mut rod=SpatialMaterialRod::try_straight(0.8,vec![source.clone()]).unwrap();
    if link {
        let member=MaterialTensionLink::try_new(source.elastic().clone(),[
            TensionLinkAttachment {node:0,offset_body_m:[0.0;3]},
            TensionLinkAttachment {node:1,offset_body_m:[0.0;3]}],0.79,1,0.1).unwrap();
        rod=rod.with_tension_links(vec![member],&gate).unwrap();
    }
    rod=rod.with_supports(&[RodSupport {node:1,position_m:[None,Some(0.0),Some(0.0)],
        frame:Some(So3::identity())}],&gate).unwrap();
    Arc::new(rod.into_motion_model(RodRootBoundary::Clamped,&gate).unwrap())
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a:e} != {b:e}");}
fn evolve(mut state:RodMotion,dt:f64,n:usize)->RodMotion {
    let gate=CancelGate::new();let loads=vec![RodNodalWrench::default();state.model().inertia().len()];
    for _ in 0..n {state=state.step_candidate(dt,&loads,controls(dt),&gate).unwrap().0;}
    state
}

#[test]
fn g1_axial_motion_matches_oscillator_with_second_order_refinement() {
    let gate=CancelGate::new();let model=axial_model(false);
    let k=model.source().sections()[0].elastic().axial_stiffness_n()/0.8;
    let omega=(k/model.inertia()[1].mass_kg()).sqrt();let mut errors=Vec::new();
    for (dt,steps) in [(1e-3,40),(5e-4,80),(2.5e-4,160)] {
        let mut pose=model.initial_pose();pose.positions_m[1][0]+=0.008;
        let state=RodMotion::new(model.clone(),pose,&[[0.0;3];2],&[[0.0;3];2],0.0,&gate).unwrap();
        let after=evolve(state,dt,steps);
        let exact=0.8+0.008*fs_math::det::cos(omega*0.04);
        errors.push((after.pose().positions_m[1][0]-exact).abs());
        assert_eq!(after.pose().positions_m[0],[0.0;3]);
    }
    for i in 0..2 {assert!((3.8..4.2).contains(&(errors[i]/errors[i+1])),"{errors:?}");}
}

#[test]
fn g0_g1_deforming_free_flight_keeps_momentum_with_second_order_energy_error() {
    let gate=CancelGate::new();
    let model=Arc::new(SpatialMaterialRod::try_straight(0.8,vec![section(1.0);2]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let mut pose=model.initial_pose();pose.positions_m[1][1]=0.012;
    pose.frames[2]=So3::exp(So3Tangent::new(Vec3::new(0.18,0.04,0.03))).unwrap();
    let initial=RodMotion::new(model,pose,&[[0.12,0.03,0.0],[0.1,-0.02,0.03],[0.08,0.01,-0.02]],
        &[[0.4,0.2,0.1],[-0.3,0.1,0.2],[0.2,-0.1,0.3]],0.0,&gate).unwrap();
    let before=initial.observe([0.0;3]).unwrap();let mut drift=Vec::new();
    for (dt,n) in [(2e-4,100),(1e-4,200),(5e-5,400)] {
        let mut state=initial.clone();
        for _ in 0..n {
            let (next,report)=state.step_candidate(dt,&[RodNodalWrench::default();3],controls(dt),&gate).unwrap();
            for c in 0..3 {
                near(report.linear_balance_residual_n_s[c],0.0,5e-12);
                near(report.angular_balance_residual_nm_s[c],0.0,5e-12);
            }
            state=next;
        }
        let after=state.observe([0.0;3]).unwrap();
        for c in 0..3 {
            near(after.linear_momentum_kg_m_s[c],before.linear_momentum_kg_m_s[c],2e-11);
            near(after.angular_momentum_kg_m2_s[c],before.angular_momentum_kg_m2_s[c],2e-11);
        }
        drift.push(((after.kinetic_energy_j-before.kinetic_energy_j)
            +(after.potential_energy_j-before.potential_energy_j)).abs());
    }
    for i in 0..2 {assert!((3.8..4.2).contains(&(drift[i]/drift[i+1])),"{drift:?}");}
}

#[test]
fn g1_spinning_free_fall_includes_section_spin_and_reference_weight() {
    let gate=CancelGate::new();
    let source=SpatialMaterialRod::try_straight(0.8,vec![section(1.0);2]).unwrap()
        .with_gravity([0.0,0.0,-9.0],&gate).unwrap();
    let model=Arc::new(source.into_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let initial=model.initial_pose();
    let mut state=RodMotion::new(model,initial.clone(),&[[0.1,0.2,-0.15];3],&[[4.0,0.0,0.0];3],0.0,&gate).unwrap();
    for _ in 0..200 {
        let (next,report)=state.step_candidate(1e-4,&[RodNodalWrench::default();3],controls(1e-4),&gate).unwrap();
        near(report.integration_defect_j,0.0,2e-11);
        state=next;
    }
    let rotation=So3::exp(So3Tangent::new(Vec3::new(0.08,0.0,0.0))).unwrap();
    for node in 0..3 {
        for c in 0..3 {
            let expected=initial.positions_m[node][c]+0.02*[0.1,0.2,-0.15][c]
                +0.5*0.02_f64.powi(2)*[0.0,0.0,-9.0][c];
            near(state.pose().positions_m[node][c],expected,2e-11);
        }
        assert!(state.pose().frames[node].body_minus(rotation).unwrap().angular.norm()<2e-11);
    }
}

#[test]
fn g1_centerline_link_changes_actual_dynamic_stiffness_and_mass() {
    let gate=CancelGate::new();let model=axial_model(true);
    let ea=model.source().sections()[0].elastic().axial_stiffness_n();
    let k0=ea/0.8;let k1=ea/0.79;let equilibrium=(k0*0.8+k1*0.79)/(k0+k1);
    let omega=((k0+k1)/model.inertia()[1].mass_kg()).sqrt();
    let mut pose=model.initial_pose();pose.positions_m[1][0]=equilibrium+0.0005;
    let state=RodMotion::new(model,pose,&[[0.0;3];2],&[[0.0;3];2],0.0,&gate).unwrap();
    let after=evolve(state,5e-4,200);
    near(after.pose().positions_m[1][0],equilibrium+0.0005*fs_math::det::cos(0.1*omega),2e-9);
}

#[test]
fn g0_stationary_supports_receive_applied_force_couple_and_weight_impulses() {
    let gate=CancelGate::new();
    let source=SpatialMaterialRod::try_straight(0.8,vec![section(1.0)]).unwrap()
        .with_supports(&[RodSupport {node:1,position_m:[Some(0.8),Some(0.0),Some(0.0)],
            frame:Some(So3::identity())}],&gate).unwrap().with_gravity([0.0,-9.0,0.0],&gate).unwrap();
    let model=Arc::new(source.into_motion_model(RodRootBoundary::Clamped,&gate).unwrap());
    let pose=model.initial_pose();
    let state=RodMotion::new(model.clone(),pose.clone(),&[[0.0;3];2],&[[0.0;3];2],0.0,&gate).unwrap();
    let loads=[RodNodalWrench::default(),RodNodalWrench {force_n:[1.0,0.2,0.0],torque_nm:[0.0,0.0,0.3]}];
    let (next,report)=state.step_candidate(1e-3,&loads,controls(1e-3),&gate).unwrap();
    assert_eq!(next.pose(),&pose);assert_eq!(report.source_work_quadrature_j,0.0);
    for node in 0..2 {
        for c in 0..3 {
            near(report.support_impulses[node].linear_n_s[c],-1e-3*(loads[node].force_n[c]
                +model.source().gravity_nodal_forces_n().unwrap()[node][c]),1e-13);
            near(report.support_impulses[node].angular_nm_s[c],-1e-3*loads[node].torque_nm[c],1e-13);
            near(report.linear_balance_residual_n_s[c],0.0,1e-13);
            near(report.angular_balance_residual_nm_s[c],0.0,1e-13);
        }
    }
}

#[test]
fn g3_g4_checkpoint_replay_and_late_energy_failure_preserve_motion() {
    let gate=CancelGate::new();let model=axial_model(false);
    let mut pose=model.initial_pose();pose.positions_m[1][0]+=0.008;
    let state=RodMotion::new(model,pose,&[[0.0;3];2],&[[0.0;3];2],0.0,&gate).unwrap();
    let straight=evolve(state.clone(),1e-3,40);
    let resumed=evolve(evolve(state.clone(),1e-3,17),1e-3,23);
    assert_eq!(straight.pose(),resumed.pose());
    assert_eq!(straight.linear_momenta_kg_m_s(),resumed.linear_momenta_kg_m_s());
    assert_eq!(straight.body_spins_kg_m2_s(),resumed.body_spins_kg_m2_s());
    assert_eq!(straight.time_s().to_bits(),resumed.time_s().to_bits());
    let mut strict=controls(1e-3);strict.energy_defect_limit_j=0.0;
    assert!(matches!(state.step_candidate(1e-3,&[RodNodalWrench::default();2],strict,&gate),Err(RodMotionError::Energy {..})));
    let after=evolve(state.clone(),1e-3,40);assert_eq!(after.pose(),straight.pose());
    let mut exhausted=controls(1e-3);exhausted.maximum_steps=40;
    assert!(straight.step_candidate(1e-3,&[RodNodalWrench::default();2],exhausted,&gate).is_err());
    assert!(state.step_candidate(1e-3,&[],controls(1e-3),&gate).is_err());
    gate.request();
    assert!(state.step_candidate(1e-3,&[RodNodalWrench::default();2],controls(1e-3),&gate).is_err());
    assert_eq!(state.accepted_steps(),0);assert_eq!(state.time_s(),0.0);
}

#[test]
fn g1_held_force_and_couple_supply_exact_uniform_acceleration_work() {
    let gate=CancelGate::new();
    let model=Arc::new(SpatialMaterialRod::try_straight(0.8,vec![section(1.0),section(2.0)]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let pose=model.initial_pose();let initial=pose.clone();let acceleration=[0.2,-0.1,0.3];
    let loads:Vec<_>=model.inertia().iter().map(|i|RodNodalWrench {
        force_n:acceleration.map(|a|a*i.mass_kg()),torque_nm:[0.7*i.principal_kg_m2()[0],0.0,0.0],
    }).collect();
    let mut state=RodMotion::new(model,pose,&[[0.1,0.0,0.0];3],&[[4.0,0.0,0.0];3],0.0,&gate).unwrap();
    let before=state.observe([0.0;3]).unwrap();let mut work=0.0;
    let mut options=controls(1e-4);options.balance_origin_m=[1.0,-2.0,3.0];
    for _ in 0..100 {
        let (next,report)=state.step_candidate(1e-4,&loads,options,&gate).unwrap();
        near(report.integration_defect_j,0.0,2e-11);work+=report.source_work_quadrature_j;
        for value in report.angular_balance_residual_nm_s {near(value,0.0,2e-11);}
        state=next;
    }
    for node in 0..3 {
        for c in 0..3 {near(state.pose().positions_m[node][c],initial.positions_m[node][c]
            +0.01*[0.1,0.0,0.0][c]+0.5*0.01_f64.powi(2)*acceleration[c],2e-11);}
    }
    let after=state.observe([0.0;3]).unwrap();
    near(after.kinetic_energy_j-before.kinetic_energy_j,work,2e-11);
}

#[test]
fn g1_g4_frozen_thermal_eigenstrain_drives_motion_and_actual_domain_crossings_refuse() {
    use fs_material::section::thermal::ConcentricThermoelasticSection;
    use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS,integrate_isotropic_thermal_expansion};
    let gate=CancelGate::new();let mut claims=ClaimSet::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),("young_modulus",1e4,Pressure::DIMS),
        ("shear_modulus",4e3,Pressure::DIMS),("linear_thermal_expansion_coefficient",1e-5,INVERSE_TEMPERATURE_DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",300.0,400.0),uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {source:"synthetic thermal rod motion".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap();
    let reference=QueryPoint::new().with("T",300.0).unwrap();let point=QueryPoint::new().with("T",400.0).unwrap();
    let requirements=[("density",Density::DIMS),("young_modulus",Pressure::DIMS),("shear_modulus",Pressure::DIMS)]
        .map(|(name,dims)|ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive).unwrap());
    let material=resolve_material_state_point(&card,&point,&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let elastic=ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {outer_radius_m:0.05,material}]).unwrap();
    let expansion=integrate_isotropic_thermal_expansion(&card,&reference,&point,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let thermal=ConcentricThermoelasticSection::try_new(elastic.clone(),vec![expansion]).unwrap();
    let source=SpatialMaterialRod::try_straight(0.8,vec![elastic.with_shear_modulus("shear_modulus",0.8).unwrap()]).unwrap()
        .with_thermal_sections(vec![thermal],0.001001,&gate).unwrap()
        .with_supports(&[RodSupport {node:1,position_m:[None,Some(0.0),Some(0.0)],frame:Some(So3::identity())}],&gate).unwrap();
    let model=Arc::new(source.into_motion_model(RodRootBoundary::Clamped,&gate).unwrap());
    let omega=(model.source().sections()[0].elastic().axial_stiffness_n()/0.8/model.inertia()[1].mass_kg()).sqrt();
    let state=RodMotion::new(model.clone(),model.initial_pose(),&[[0.0;3];2],&[[0.0;3];2],0.0,&gate).unwrap();
    let after=evolve(state,5e-4,80);
    near(after.pose().positions_m[1][0],0.8+0.0008*(1.0-fs_math::det::cos(omega*0.04)),2e-9);
    let fast=RodMotion::new(model.clone(),model.initial_pose(),&[[0.0;3],[-1.0,0.0,0.0]],&[[0.0;3];2],0.0,&gate).unwrap();
    assert!(matches!(fast.step_candidate(1e-3,&[RodNodalWrench::default();2],controls(1e-3),&gate),Err(RodMotionError::Solid(_))));
    assert_eq!(fast.pose(),&model.initial_pose());assert_eq!(fast.accepted_steps(),0);
}

#[test]
fn g3_intrinsically_curved_heterogeneous_rod_translates_without_prestress() {
    use fs_solid::rod::material_field::reference::RodReferenceGeometry;
    let gate=CancelGate::new();
    let angles=[0.0,0.2,0.4];
    let nodes:Vec<_>=angles.iter().map(|&t|[2.0*fs_math::det::sin(t),2.0*(1.0-fs_math::det::cos(t)),0.0]).collect();
    let frames:Vec<_>=angles.iter().map(|&t|So3::exp(So3Tangent::new(Vec3::new(0.0,0.0,t))).unwrap()).collect();
    let reference=RodReferenceGeometry::try_new(0.8,nodes,frames).unwrap();
    let source=SpatialMaterialRod::try_from_reference(reference,vec![section(1.0),section(2.0)]).unwrap();
    let model=Arc::new(source.into_motion_model(RodRootBoundary::Free,&gate).unwrap());
    let initial=model.initial_pose();
    let state=RodMotion::new(model,initial.clone(),&[[0.1,0.05,0.02];3],&[[0.0;3];3],0.0,&gate).unwrap();
    let after=evolve(state,1e-4,100);
    for node in 0..3 {for c in 0..3 {near(after.pose().positions_m[node][c],
        initial.positions_m[node][c]+0.01*[0.1,0.05,0.02][c],2e-10);}}
    near(after.observe([0.0;3]).unwrap().potential_energy_j,0.0,1e-15);
}

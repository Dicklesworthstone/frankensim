//! Native impulses on mechanisms with simultaneous material and joint closures.
//! G0/G1/G3/G4: independent inertia reductions, work, frame changes and refusal.

use std::cell::Cell;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::MaterialPoint;
use fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialError;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralImpulseResponse, BilateralLoop, BilateralLoopError,
    DistanceLoop, PointCoincidenceLoop, joint,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionError, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits,
    JointModel, Link, SpatialInertia, apply_free_floating_impulses, forward_kinematics,
};

fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 3e-10*(1.+b.abs()), "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a-b).norm() < 2e-9, "{a:?} != {b:?}");
}
fn inertia(mass: f64, spin: f64) -> SpatialInertia {
    SpatialInertia::new(mass, Vec3::default(), Mat3 {
        m: [spin,0.,0.,0.,spin,0.,0.,0.,spin],
    }).unwrap()
}
fn model(limit: Option<f64>) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier", None, Se3::identity(), JointModel::FIXED, inertia(5.,4.),
    )];
    for (i,mass) in [2.,3.,4.].into_iter().enumerate() {
        links.push(Link::new(format!("slider {i}"),Some(0),Se3::identity(),
            JointModel::prismatic(Vec3::new(1.,0.,0.),
                limit.map(|v| JointLimits::new(-10.,10.,v,0.1).unwrap())).unwrap(),
            inertia(mass,mass)));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(Se3::identity(),
            Twist::new(Vec3::default(),Vec3::new(0.7,0.,0.))),
        joint_position: vec![1.,3.,1.5],
        joint_velocity: vec![0.4,0.4,0.2],
    }
}
fn loops() -> [BilateralLoop;1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint { link:1,point_body:Vec3::default() },
        b: MaterialPoint { link:2,point_body:Vec3::default() },
        length_m:2.,
    })]
}
fn row() -> JointConstraint {
    JointConstraint {
        coordinate:ConstraintCoordinate::Length, coefficients:vec![0.,1.,-2.],
        position:0.,position_tolerance:1e-9,velocity_tolerance:1e-9,
    }
}
fn source() -> ([f64;3],[Wrench;4]) {
    let mut body = [Wrench::default();4];
    body[0].force.x=0.42;
    ([0.6,-0.2,0.4],body)
}
fn work(result: &BilateralImpulseResponse) {
    close(result.impulse.impulse_work,
        result.source_work_j + result.constraint_work_j + result.work_partition_residual_j);
    close(result.impulse.kinetic_energy_change,
        result.impulse.impulse_work + result.impulse.kinetic_work_residual);
    close(result.impulse.kinetic_work_residual,0.);
    close(result.work_partition_residual_j,0.);
}
fn reaction_world(
    model: &ArticulatedModel, state: &FreeFloatingMotionState, result: &BilateralImpulseResponse,
) -> Wrench {
    let kin=forward_kinematics(model,
        BaseState::prescribed(state.base.world_from_base,state.base.twist_body,Twist::zero()),
        &state.joint_position,&state.joint_velocity).unwrap();
    result.constraint_body_impulses.iter().zip(&kin.world_from_link)
        .fold(Wrench::default(),|sum,(impulse,pose)| {
            let w=impulse.transform_by(pose);
            Wrench::new(sum.torque+w.torque,sum.force+w.force)
        })
}

#[test]
fn simultaneous_source_and_constraint_impulses_match_independent_geared_inertia() {
    let model=model(None);
    let start=initial();
    let (applied,body)=source();
    let rows=[row()];
    let result=joint::apply_impulses(
        &model,&start,&applied,&body,&loops(),&rows,Default::default(),||false,
    ).unwrap();
    work(&result);
    // qdot=t*u, t=(2,2,1). With M=diag(2,3,4)-m*m^T/14,
    // t^T M t=10 and base delta=J/14-delta_u.
    let du=(2.*applied[0]+2.*applied[1]+applied[2]-body[0].force.x)/10.;
    let db=body[0].force.x/14.-du;
    for (i,t) in [2.,2.,1.].into_iter().enumerate() {
        close(result.impulse.joint_velocity_after[i],t*(0.2+du));
    }
    close(result.impulse.base_after.twist_body.linear.x,0.7+db);
    assert_eq!(result.impulse.base_after.world_from_base,start.base.world_from_base);
    let mu=applied[0]-body[0].force.x/7.-2.*du;
    let lambda=0.5*(applied[2]-2.*body[0].force.x/7.);
    close(result.loop_multipliers_ns[0],mu);
    close(result.joint_multipliers[0],lambda);
    close(result.constraint_body_impulses[1].force.x,-mu);
    close(result.constraint_body_impulses[2].force.x,mu);
    close(result.constraint_joint_impulses[1],lambda);
    close(result.constraint_joint_impulses[2],-2.*lambda);
    // Independent normalized W from [[5/6,1/3],[1/3,23/15]].
    let cross=(2.0_f64/23.).sqrt();
    close(result.reciprocal_condition,(1.-cross)/(1.+cross));
    close(result.minimum_scaled_pivot,21./23.);
    close(result.constraint_work_j,0.);
    close(result.row_constraint_work_j,0.);
    let expected_work=body[0].force.x*(0.7+0.5*db)
        +(2.*applied[0]+2.*applied[1]+applied[2])*(0.2+0.5*du);
    close(result.source_work_j,expected_work);
    close(result.impulse.kinetic_energy_change,expected_work);
    vector(result.impulse.linear_momentum_change_world,Vec3::new(0.42,0.,0.));
    vector(reaction_world(&model,&start,&result).force,Vec3::default());
    vector(reaction_world(&model,&start,&result).torque,Vec3::default());
    for r in result.velocity_residual_after_mps.iter()
        .chain(&result.joint_velocity_residuals_after) { close(*r,0.); }

    // Equivalent all-joint rows exercise a genuinely empty body-geometry block.
    let joint_only=[JointConstraint {
        coefficients:vec![-1.,1.,0.],position:2.,..row()
    },row()];
    let pure=joint::apply_impulses(
        &model,&start,&applied,&body,&[],&joint_only,Default::default(),||false,
    ).unwrap();
    assert!(pure.rows.is_empty());
    assert!(pure.loop_multipliers_ns.is_empty());
    close(pure.joint_multipliers[0],mu);
    close(pure.joint_multipliers[1],lambda);
    for i in 0..3 {
        close(pure.impulse.joint_velocity_after[i],result.impulse.joint_velocity_after[i]);
    }
    vector(pure.impulse.base_after.twist_body.linear,result.impulse.base_after.twist_body.linear);
    close(pure.source_work_j,expected_work);
    work(&pure);

    // Both native sources are in material coordinates; a world transform and
    // common inertial boost retain the relative mechanism and its reactions.
    let rotation=So3::exp(So3Tangent::new(Vec3::new(0.3,-0.2,0.4))).unwrap();
    let translation=Vec3::new(1.2,-0.7,0.9);
    let boost=Vec3::new(-0.5,0.3,0.2);
    let mut moved=start.clone();
    moved.base.world_from_base=Se3::from_parts(rotation,translation).unwrap();
    moved.base.twist_body.linear=start.base.twist_body.linear
        +rotation.inverse().rotate(boost).unwrap();
    let shifted=joint::apply_impulses(
        &model,&moved,&applied,&body,&loops(),&rows,Default::default(),||false,
    ).unwrap();
    work(&shifted);
    for i in 0..3 { close(shifted.impulse.joint_velocity_after[i],result.impulse.joint_velocity_after[i]); }
    vector(shifted.impulse.base_after.twist_body.linear,
        result.impulse.base_after.twist_body.linear+rotation.inverse().rotate(boost).unwrap());
    close(shifted.loop_multipliers_ns[0],mu);
    close(shifted.joint_multipliers[0],lambda);
    close(shifted.constraint_work_j,0.);
    close(shifted.source_work_j,expected_work
        +rotation.rotate(body[0].force).unwrap().dot(boost));
    let force=rotation.rotate(body[0].force).unwrap();
    vector(shifted.impulse.linear_momentum_change_world,force);
    vector(shifted.impulse.angular_momentum_change_world,translation.cross(force));
    assert_eq!(start,initial());
}

#[test]
fn constrained_event_admits_valid_final_speed_without_applying_the_free_probe() {
    let model=model(Some(0.6));
    let start=initial();
    let (applied,body)=source();
    assert!(matches!(apply_free_floating_impulses(
        &model,start.base,&start.joint_position,&start.joint_velocity,&applied,&body,
    ),Err(ArticulatedError::JointVelocityLimitViolation { .. })));
    let valid=joint::apply_impulses(
        &model,&start,&applied,&body,&loops(),&[row()],Default::default(),||false,
    ).unwrap();
    close(valid.impulse.joint_velocity_after[0],0.556);
    close(valid.impulse.joint_velocity_after[1],0.556);
    close(valid.impulse.joint_velocity_after[2],0.278);
    let mut too_large=body;
    too_large[0].force.x*=10.;
    let over=applied.map(|v|10.*v);
    assert!(matches!(joint::apply_impulses(
        &model,&start,&over,&too_large,&loops(),&[row()],Default::default(),||false,
    ),Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
        FreeFloatingMotionError::Dynamics(ArticulatedError::JointVelocityLimitViolation { .. })
    )))));
}

#[test]
fn finite_tangent_work_is_explicit_and_zero_source_preserves_the_admitted_defect() {
    let model=model(None);
    let mut start=initial();
    let defect=1e-3;
    start.joint_velocity[2]-=0.5*defect;
    let rows=[JointConstraint { velocity_tolerance:2e-3,..row() }];
    let zero_body=[Wrench::default();4];
    // M*t*0.1 = (0.2,0.3,0), an independently tangent source increment.
    let result=joint::apply_impulses(
        &model,&start,&[0.2,0.3,0.],&zero_body,&loops(),&rows,
        Default::default(),||false,
    ).unwrap();
    work(&result);
    // Eliminate the body row: effective gear mobility is 7/5.
    close(result.loop_multipliers_ns[0],2.*defect/7.);
    close(result.joint_multipliers[0],-5.*defect/7.);
    close(result.constraint_work_j,-5.*defect*defect/14.);
    close(result.row_constraint_work_j,-5.*defect*defect/14.);
    close(result.source_work_j,0.25);
    close(result.impulse.kinetic_energy_change,0.25-5.*defect*defect/14.);
    close(result.joint_velocity_residuals_before[0],defect);
    close(result.joint_velocity_residuals_after[0],0.);
    let passive=joint::apply_impulses(
        &model,&start,&[0.;3],&zero_body,&loops(),&rows,Default::default(),||false,
    ).unwrap();
    assert_eq!(passive.impulse.base_after,start.base);
    assert_eq!(passive.impulse.joint_velocity_after,start.joint_velocity);
    assert_eq!(passive.constraint_joint_impulses,vec![0.;3]);
    assert_eq!(passive.constraint_body_impulses,vec![Wrench::default();4]);
    assert_eq!(passive.joint_velocity_residuals_after,passive.joint_velocity_residuals_before);
    close(passive.constraint_work_j,0.);
    close(passive.row_constraint_work_j,0.);
    close(passive.source_work_j,0.);
    assert!(passive.reciprocal_condition>0.);
}

fn combined(parts: &[SpatialInertia], offsets: &[Vec3]) -> SpatialInertia {
    let mass: f64=parts.iter().map(|p|p.mass()).sum();
    let center=parts.iter().zip(offsets).fold(Vec3::default(),|c,(p,r)|
        c+(p.center_of_mass()+*r).scale(p.mass()/mass));
    let mut matrix=[0.;9];
    for (p,r) in parts.iter().zip(offsets) {
        let d=p.center_of_mass()+*r-center;
        let v=[d.x,d.y,d.z];
        for i in 0..3 { for j in 0..3 {
            matrix[3*i+j]+=p.inertia_com().m[3*i+j]
                +p.mass()*((if i==j {d.dot(d)} else {0.})-v[i]*v[j]);
        } }
    }
    SpatialInertia::new(mass,center,Mat3 {m:matrix}).unwrap()
}

#[test]
fn spherical_and_angular_reactions_match_full_tensor_assembly_at_finite_admitted_gap() {
    let parts:Vec<_>=[5.,1.,1.,1.,2.].into_iter().enumerate().map(|(i,m)| {
        let f=i as f64;
        SpatialInertia::new(m,Vec3::new(0.03*f,-0.02*f,0.01*f),
            Mat3 {m:[m,0.1*m,-0.05*m,0.1*m,1.2*m,0.08*m,-0.05*m,0.08*m,1.5*m]}).unwrap()
    }).collect();
    let axes=[Vec3::new(1.,0.,0.),Vec3::new(0.,1.,0.),Vec3::new(0.,0.,1.)];
    let mut links=vec![Link::new("base",None,Se3::identity(),JointModel::FIXED,parts[0])];
    for i in 0..3 {
        links.push(Link::new(format!("slide {i}"),Some(i),Se3::identity(),
            JointModel::prismatic(axes[i],None).unwrap(),parts[i+1]));
    }
    links.push(Link::new("rotor",Some(3),Se3::identity(),
        JointModel::revolute(axes[2],None).unwrap(),parts[4]));
    let model=ArticulatedModel::new(links).unwrap();
    let q=vec![2e-4,-1e-4,1.5e-4,0.];
    let offsets=[Vec3::default(),Vec3::new(q[0],0.,0.),Vec3::new(q[0],q[1],0.),
        Vec3::new(q[0],q[1],q[2]),Vec3::new(q[0],q[1],q[2])];
    let reduced=ArticulatedModel::new(vec![Link::new(
        "assembled",None,Se3::identity(),JointModel::FIXED,combined(&parts,&offsets),
    )]).unwrap();
    let base=FreeFloatingBaseState::new(Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3,-0.2,0.1))).unwrap(),
        Vec3::new(0.5,-0.3,0.2),
    ).unwrap(),Twist::new(Vec3::new(0.2,-0.15,0.1),Vec3::new(0.3,-0.2,0.4)));
    let start=FreeFloatingMotionState {base,joint_position:q,joint_velocity:vec![0.;4]};
    let closures=[BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a:MaterialPoint {link:0,point_body:Vec3::default()},
        b:MaterialPoint {link:3,point_body:Vec3::default()},
    })];
    let lock=[JointConstraint {
        coordinate:ConstraintCoordinate::Angle,coefficients:vec![0.,0.,0.,1.],
        position:0.,position_tolerance:1e-9,velocity_tolerance:1e-9,
    }];
    let mut body=[Wrench::default();5];
    body[4]=Wrench::new(Vec3::new(0.4,0.3,-0.2),Vec3::new(2.,-3.,1.));
    let root_body=[Wrench::new(
        body[4].torque+offsets[4].cross(body[4].force),body[4].force,
    )];
    let config=BilateralImpulseConfig {
        position_tolerance_m:1e-3,velocity_tolerance_mps:1e-3,..Default::default()
    };
    let result=joint::apply_impulses(
        &model,&start,&[0.1,-0.2,0.3,0.4],&body,&closures,&lock,config,||false,
    ).unwrap();
    let reference=apply_free_floating_impulses(
        &reduced,base,&[],&[],&[],&root_body,
    ).unwrap();
    work(&result);
    vector(result.impulse.base_after.twist_body.linear,reference.base_after.twist_body.linear);
    vector(result.impulse.base_after.twist_body.angular,reference.base_after.twist_body.angular);
    for rate in &result.impulse.joint_velocity_after {close(*rate,0.);}
    close(result.impulse.kinetic_energy_change,reference.kinetic_energy_change);
    close(result.source_work_j,reference.impulse_work);
    close(result.constraint_work_j,0.);
    close(result.row_constraint_work_j,0.);
    assert_eq!(result.rows.len(),3);
    assert_eq!(result.joint_multipliers.len(),1);
    assert!(result.point_residual_after[0].separation_world_m.norm()>1e-4);
    assert!(result.point_residual_after[0].relative_velocity_world_mps.norm()>1e-6);
    vector(reaction_world(&model,&start,&result).force,Vec3::default());
    vector(reaction_world(&model,&start,&result).torque,Vec3::default());
    vector(result.impulse.external_linear_impulse_world,reference.external_linear_impulse_world);
    vector(result.impulse.external_angular_impulse_world,reference.external_angular_impulse_world);
}

#[test]
fn empty_constraints_delegate_and_zero_source_still_admits_rank_and_cancellation() {
    let model=model(None);
    let start=initial();
    let saved=start.clone();
    let (applied,body)=source();
    let native=apply_free_floating_impulses(
        &model,start.base,&start.joint_position,&start.joint_velocity,&applied,&body,
    ).unwrap();
    let empty=joint::apply_impulses(
        &model,&start,&applied,&body,&[],&[],Default::default(),||false,
    ).unwrap();
    assert_eq!(empty.impulse,native);
    assert_eq!(empty.source_work_j,native.impulse_work);
    assert_eq!(empty.constraint_work_j,0.);
    assert_eq!(empty.reciprocal_condition,1.);
    let duplicate=[JointConstraint {
        coefficients:vec![-1.,1.,0.],position:2.,..row()
    }];
    assert!(matches!(joint::apply_impulses(
        &model,&start,&[0.;3],&[Wrench::default();4],&loops(),&duplicate,
        Default::default(),||false,
    ),Err(BilateralLoopError::Mobility { .. })));
    let polls=Cell::new(0);
    let reference=joint::apply_impulses(
        &model,&start,&applied,&body,&loops(),&[row()],Default::default(),
        || {polls.set(polls.get()+1);false},
    ).unwrap();
    for threshold in [1,polls.get()/2,polls.get()] {
        let mut visited=0;
        let refused=joint::apply_impulses(
            &model,&start,&applied,&body,&loops(),&[row()],Default::default(),
            || {visited+=1;visited>=threshold},
        );
        assert!(matches!(refused,Err(BilateralLoopError::Motion(
            ImplicitSpatialError::Motion(FreeFloatingMotionError::Cancelled)
        ))));
    }
    assert_eq!(joint::apply_impulses(
        &model,&start,&applied,&body,&loops(),&[row()],Default::default(),||false,
    ).unwrap(),reference);
    assert_eq!(start,saved);
}

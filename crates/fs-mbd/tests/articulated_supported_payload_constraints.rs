//! G0/G1/G3/G4: supported finite payload capture through native ideal joint rows.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::inventory::supported::{
    SupportedPayloadControls, capture_payload as unconstrained_capture,
};
use fs_mbd::articulated::motion::inventory::supported::constraints::{
    ConstrainedSupportedPayloadCapture, SupportedPayloadConstraintControls, capture_payload,
};
use fs_mbd::articulated::motion::inventory::supported::release::{
    SupportedPayloadReleaseControls, release_payload,
};
use fs_mbd::articulated::motion::inventory::{PayloadError, PayloadState};
use fs_mbd::articulated::motion::prescribed::{PrescribedBaseMotionState, PrescribedBasePath};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitSpatialConfig,
    constraints::{SupportedConstraintControls, step_prescribed_base_spatial},
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia,
    forward_kinematics,
};
use std::cell::Cell;

fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-9 * (1. + a.abs().max(b.abs())), "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x,b.x); close(a.y,b.y); close(a.z,b.z);
}
fn pose(axis: Vec3, position: Vec3) -> Se3 {
    Se3::from_parts(So3::exp(So3Tangent::new(axis)).unwrap(), position).unwrap()
}
fn inertia(mass: f64, center: Vec3, diagonal: [f64; 3], axes: Vec3) -> SpatialInertia {
    let r = So3::exp(So3Tangent::new(axes)).unwrap().matrix();
    SpatialInertia::new(mass, center, r.compose(Mat3 {
        m: [diagonal[0],0.,0., 0.,diagonal[1],0., 0.,0.,diagonal[2]],
    }).compose(r.transpose())).unwrap()
}
fn sphere(mass: f64) -> SpatialInertia {
    inertia(mass, Vec3::default(), [1.,1.,1.], Vec3::default())
}
fn controls(tolerances: &[f64]) -> SupportedPayloadControls<'_> {
    SupportedPayloadControls {
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_removed_relative_energy_j: 1000.,
    }
}
fn row_controls() -> SupportedPayloadConstraintControls {
    SupportedPayloadConstraintControls {
        minimum_reciprocal_condition: 1e-10,
        maximum_reaction_work_j: 2e-8,
    }
}
fn row(coordinate: ConstraintCoordinate, coefficients: Vec<f64>, position: f64) -> JointConstraint {
    JointConstraint {
        coordinate, coefficients, position,
        position_tolerance: 2e-10, velocity_tolerance: 2e-10,
    }
}
fn geared_sliders(support_speed: f64, first_speed_limit: f64)
    -> (ArticulatedModel, PrescribedBaseMotionState, PayloadState, Vec<JointConstraint>)
{
    let slider = |speed| JointModel::prismatic(
        Vec3::new(1.,0.,0.), Some(JointLimits::new(-10.,10.,speed,0.01).unwrap()),
    ).unwrap();
    let model = ArticulatedModel::new(vec![
        Link::new("support frame", None, Se3::identity(), JointModel::FIXED, sphere(1.)),
        Link::new("catcher", Some(0), Se3::identity(), slider(first_speed_limit), sphere(3.)),
        Link::new("geared carriage", Some(0), Se3::identity(), slider(10.), sphere(5.)),
    ]).unwrap();
    let state = PrescribedBaseMotionState {
        time_s: 4.,
        base: BaseState::prescribed(
            Se3::identity(), Twist::new(Vec3::default(), Vec3::new(support_speed,0.,0.)),
            Twist::zero(),
        ),
        joint_position: vec![0.4,0.2],
        joint_velocity: vec![2.,1.],
    };
    let payload = PayloadState {
        inertia: sphere(2.),
        world_from_body: pose(Vec3::default(), Vec3::new(0.4,0.,0.)),
        twist_body: Twist::new(Vec3::default(), Vec3::new(support_speed+7.,0.,0.)),
    };
    (model,state,payload,vec![row(ConstraintCoordinate::Length,vec![1.,-2.],0.)])
}

struct Physical {
    linear: Vec3,
    angular: Vec3,
    energy: f64,
    joint: Vec<f64>,
}

// Only native FK is shared with the implementation. Momentum and kinetic
// energy use independent COM/spin formulas, and generalized momenta project
// those physical quantities on each ancestor joint's world screw.
fn physical(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    payload: Option<(usize, PayloadState)>,
) -> Physical {
    let kin = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let mut bodies: Vec<_> = model
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| (i, link.inertia(), kin.world_from_link[i], kin.body_twist[i]))
        .collect();
    if let Some((target, incoming)) = payload {
        bodies.push((
            target,
            incoming.inertia,
            incoming.world_from_body,
            incoming.twist_body,
        ));
    }
    let mut result = Physical {
        linear: Vec3::default(),
        angular: Vec3::default(),
        energy: 0.,
        joint: vec![0.; model.dof_count()],
    };
    for (target, mass, pose, twist) in bodies {
        let center = pose.transform_point(mass.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(mass.center_of_mass()))
            .unwrap();
        let linear = velocity.scale(mass.mass());
        let spin_body = mass.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        result.linear = result.linear + linear;
        result.angular = result.angular + spin + center.cross(linear);
        result.energy +=
            0.5 * mass.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
        let mut ancestor = Some(target);
        while let Some(index) = ancestor {
            if let Some(dof) = model.dof_index(index) {
                let screw = model.links()[index].joint().motion_subspace();
                let joint_pose = kin.world_from_link[index];
                let angular = joint_pose.rotation().rotate(screw.angular).unwrap();
                let linear_axis = joint_pose.rotation().rotate(screw.linear).unwrap();
                result.joint[dof] += angular
                    .dot(spin + (center - joint_pose.translation()).cross(linear))
                    + linear_axis.dot(linear);
            }
            ancestor = model.links()[index].parent();
        }
    }
    result
}


fn audit(
    model: &ArticulatedModel, state: &PrescribedBaseMotionState, target: usize,
    incoming: PayloadState, result: &ConstrainedSupportedPayloadCapture,
) {
    let capture = &result.capture;
    let before = physical(model,state,Some((target,incoming)));
    let after = physical(&capture.model_after,&capture.state_after,None);
    let balance = &capture.balance;
    vector(after.linear-before.linear,balance.support_impulse_world.force);
    vector(after.angular-before.angular,balance.support_impulse_world.torque);
    vector(balance.linear_momentum_change_world,after.linear-before.linear);
    vector(balance.angular_momentum_change_world,after.angular-before.angular);
    vector(balance.linear_momentum_residual_world,Vec3::default());
    vector(balance.angular_momentum_residual_world,Vec3::default());
    for (i,(&old,&new)) in before.joint.iter().zip(&after.joint).enumerate() {
        let applied = result.reaction.joint_impulse.get(i).copied().unwrap_or(0.);
        close(new-old,applied);
        close(balance.joint_momentum_residual[i],new-old-applied);
    }
    close(balance.kinetic_energy_change_j,after.energy-before.energy);
    close(balance.kinetic_energy_change_j + capture.removed_relative_energy_j,
          balance.support_work_j + result.reaction.final_rate_pairing_j);
    close(balance.support_work_j,
        balance.support_impulse_body.force.dot(state.base.twist_body.linear)
        + balance.support_impulse_body.torque.dot(state.base.twist_body.angular));
    close(balance.energy_residual_j,result.reaction.final_rate_pairing_j);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    assert_eq!(capture.state_after.base,state.base);
    assert_eq!(capture.state_after.joint_position,state.joint_position);
    assert_eq!(capture.state_after.time_s,state.time_s);
}

#[test]
fn geared_capture_has_finite_recoil_and_signed_support_work() {
    for support_speed in [0.,2.,-2.] {
        let (model,state,payload,rows) = geared_sliders(support_speed,3.8);
        // The free catcher would reach 4 m/s and violate its physical limit.
        // Its simultaneous geared endpoint is 3.6 m/s and must be accepted.
        assert!(matches!(
            unconstrained_capture(&model,&state,1,payload,controls(&[2e-9;2]),||false),
            Err(PayloadError::Dynamics(ArticulatedError::JointVelocityLimitViolation{..}))
        ));
        let result = capture_payload(&model,&state,1,payload,&rows,
            controls(&[2e-9;2]),row_controls(),||false).unwrap();
        close(result.capture.state_after.joint_velocity[0],3.6);
        close(result.capture.state_after.joint_velocity[1],1.8);
        close(result.reaction.multipliers[0],-2.);
        close(result.reaction.joint_impulse[0],-2.);
        close(result.reaction.joint_impulse[1],4.);
        vector(result.capture.balance.support_impulse_body.force,Vec3::new(2.,0.,0.));
        close(result.capture.balance.support_work_j,2.*support_speed);
        close(result.capture.balance.kinetic_energy_change_j,-17.+2.*support_speed);
        close(result.capture.removed_relative_energy_j,17.);
        close(result.reaction.midpoint_work_j,0.);
        close(result.reaction.final_rate_pairing_j,0.);
        close(result.reaction.velocity_residuals_after[0],0.);
        close(result.reaction.reciprocal_condition,1.);
        audit(&model,&state,1,payload,&result);
    }
}

#[test]
fn locked_scalar_root_retains_support_reaction_without_effort_or_trial_speed_clipping() {
    let model = ArticulatedModel::new(vec![Link::new(
        "locked catcher",None,Se3::identity(),
        JointModel::prismatic(Vec3::new(1.,0.,0.),
            Some(JointLimits::new(-1.,1.,0.1,0.01).unwrap())).unwrap(),
        sphere(3.),
    )]).unwrap();
    let state = PrescribedBaseMotionState {
        time_s: 1.,
        base: BaseState::prescribed(Se3::identity(),
            Twist::new(Vec3::default(),Vec3::new(2.,0.,0.)),Twist::zero()),
        joint_position: vec![0.4], joint_velocity: vec![0.],
    };
    let payload = PayloadState {
        inertia: sphere(2.),
        world_from_body: pose(Vec3::default(),Vec3::new(0.4,0.,0.)),
        twist_body: Twist::new(Vec3::default(),Vec3::new(7.,0.,0.)),
    };
    let rows = [row(ConstraintCoordinate::Length,vec![1.],0.4)];
    let result = capture_payload(&model,&state,0,payload,&rows,
        controls(&[2e-9]),row_controls(),||false).unwrap();
    close(result.capture.state_after.joint_velocity[0],0.);
    close(result.reaction.joint_impulse[0],-10.);
    vector(result.capture.balance.support_impulse_body.force,Vec3::new(-10.,0.,0.));
    close(result.capture.balance.support_work_j,-20.);
    close(result.capture.removed_relative_energy_j,25.);
    close(result.capture.balance.kinetic_energy_change_j,-45.);
    audit(&model,&state,0,payload,&result);
}

fn branched()
    -> (ArticulatedModel, PrescribedBaseMotionState, PayloadState, Vec<JointConstraint>)
{
    let model = ArticulatedModel::new(vec![
        Link::new("carrier",None,pose(Vec3::new(0.1,-0.2,0.3),Vec3::new(-0.2,0.3,0.1)),
            JointModel::FIXED,sphere(4.)),
        Link::new("capture arm",Some(0),pose(Vec3::new(0.2,0.1,-0.3),Vec3::new(0.5,-0.4,0.2)),
            JointModel::revolute(Vec3::new(0.,0.,1.),None).unwrap(),
            inertia(3.,Vec3::new(0.2,-0.1,0.4),[2.,3.,4.],Vec3::new(0.1,0.3,-0.2))),
        Link::new("helical sibling",Some(0),pose(Vec3::new(-0.1,0.4,0.2),Vec3::new(-0.5,0.2,-0.3)),
            JointModel::helical(Vec3::new(0.,1.,0.),0.3,None).unwrap(),
            inertia(2.,Vec3::new(-0.2,0.3,0.1),[0.8,1.1,1.4],Vec3::new(-0.4,0.2,0.1))),
    ]).unwrap();
    let state = PrescribedBaseMotionState {
        time_s: 3.,
        base: BaseState::prescribed(
            pose(Vec3::new(0.3,-0.4,0.2),Vec3::new(2.,-3.,1.)),
            Twist::new(Vec3::new(0.2,-0.3,0.4),Vec3::new(1.1,-0.6,0.2)),
            Twist::zero(),
        ),
        joint_position: vec![0.3,-0.15],joint_velocity: vec![0.6,-0.3],
    };
    let kin = forward_kinematics(&model,state.base,&state.joint_position,&state.joint_velocity).unwrap();
    let payload = PayloadState {
        inertia: inertia(1.5,Vec3::new(0.1,-0.2,0.3),[0.4,0.5,0.6],Vec3::new(0.2,-0.3,0.4)),
        world_from_body: kin.world_from_link[1]
            .compose(pose(Vec3::new(0.4,-0.2,0.1),Vec3::new(0.7,-0.2,0.3))).unwrap(),
        twist_body: Twist::new(Vec3::new(0.3,0.2,-0.4),Vec3::new(0.8,-1.1,0.5)),
    };
    (model,state,payload,vec![row(ConstraintCoordinate::Angle,vec![1.,2.],0.)])
}

fn attached(
    model: &ArticulatedModel, state: &PrescribedBaseMotionState,
    target: usize, incoming: PayloadState,
) -> PayloadState {
    let kin = forward_kinematics(model,state.base,&state.joint_position,&state.joint_velocity).unwrap();
    let body_from_link = incoming.world_from_body.inverse().unwrap()
        .compose(kin.world_from_link[target]).unwrap();
    PayloadState {twist_body: kin.body_twist[target].transform_by(&body_from_link),..incoming}
}

#[test]
fn full_tensor_branch_matches_independent_reduced_inertia_and_row_scaling() {
    let (model,state,payload,rows) = branched();
    let initial = physical(&model,&state,Some((1,payload)));
    let mut zero = state.clone();
    zero.joint_velocity = vec![0.,0.];
    let p0 = physical(&model,&zero,Some((1,attached(&model,&zero,1,payload))));
    let mut unit = state.clone();
    unit.joint_velocity = vec![2.,-1.];
    let p1 = physical(&model,&unit,Some((1,attached(&model,&unit,1,payload))));
    // Independent one-coordinate kinetic metric with tangent [2,-1].
    // This uses original separate body tensors and COM/spin sums, never ABA
    // or the production loaded tensor to construct the expected rate.
    let reduced_momentum = 2.*(initial.joint[0]-p0.joint[0])-(initial.joint[1]-p0.joint[1]);
    let reduced_mass = 2.*(p1.joint[0]-p0.joint[0])-(p1.joint[1]-p0.joint[1]);
    let expected = reduced_momentum/reduced_mass;
    let result = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    close(result.capture.state_after.joint_velocity[0],2.*expected);
    close(result.capture.state_after.joint_velocity[1],-expected);
    audit(&model,&state,1,payload,&result);
    assert_ne!(result.capture.state_after.joint_velocity[0],state.joint_velocity[0]);
    assert_ne!(result.capture.state_after.joint_velocity[1],state.joint_velocity[1]);
    let mut scaled = rows.clone();
    for coefficient in &mut scaled[0].coefficients { *coefficient *= 1000.; }
    scaled[0].position *= 1000.;
    scaled[0].position_tolerance *= 1000.;
    scaled[0].velocity_tolerance *= 1000.;
    let other = capture_payload(&model,&state,1,payload,&scaled,
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    for (&a,&b) in result.capture.state_after.joint_velocity.iter()
        .zip(&other.capture.state_after.joint_velocity) { close(a,b); }
    close(result.reaction.multipliers[0],1000.*other.reaction.multipliers[0]);
    close(result.capture.removed_relative_energy_j,other.capture.removed_relative_energy_j);
    audit(&model,&state,1,payload,&other);
}

#[test]
fn constrained_capture_is_covariant_under_world_transform_and_inertial_boost() {
    let (model,state,payload,rows) = branched();
    let result = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    let change = pose(Vec3::new(-0.4,0.6,0.3),Vec3::new(-5.,3.,7.));
    let boost = Vec3::new(17.,-9.,6.);
    let base_pose = change.compose(state.base.world_from_base).unwrap();
    let base_boost = base_pose.rotation().inverse().rotate(boost).unwrap();
    let shifted = PrescribedBaseMotionState {
        base: BaseState::prescribed(base_pose,
            state.base.twist_body.plus(Twist::new(Vec3::default(),base_boost)),
            Twist::new(state.base.acceleration_body.angular,
                state.base.acceleration_body.linear-state.base.twist_body.angular.cross(base_boost))),
        ..state.clone()
    };
    let payload_pose = change.compose(payload.world_from_body).unwrap();
    let incoming = PayloadState {
        world_from_body: payload_pose,
        twist_body: payload.twist_body.plus(Twist::new(Vec3::default(),
            payload_pose.rotation().inverse().rotate(boost).unwrap())),
        ..payload
    };
    let transformed = capture_payload(&model,&shifted,1,incoming,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    let old = &result.capture;
    let new = &transformed.capture;
    for i in 0..2 {
        close(old.state_after.joint_velocity[i],new.state_after.joint_velocity[i]);
        close(result.reaction.joint_impulse[i],transformed.reaction.joint_impulse[i]);
    }
    close(old.removed_relative_energy_j,new.removed_relative_energy_j);
    vector(old.balance.support_impulse_body.force,new.balance.support_impulse_body.force);
    vector(old.balance.support_impulse_body.torque,new.balance.support_impulse_body.torque);
    let force = change.rotation().rotate(old.balance.support_impulse_world.force).unwrap();
    let torque = change.rotation().rotate(old.balance.support_impulse_world.torque).unwrap()
        + change.translation().cross(force);
    vector(new.balance.support_impulse_world.force,force);
    vector(new.balance.support_impulse_world.torque,torque);
    close(new.balance.support_work_j,old.balance.support_work_j+boost.dot(force));
    close(new.balance.kinetic_energy_change_j,old.balance.kinetic_energy_change_j+boost.dot(force));
    close(new.balance.kinetic_before_j,old.balance.kinetic_before_j);
    close(new.balance.kinetic_after_j,old.balance.kinetic_after_j);
    close(new.balance.accounting_support_work_j,old.balance.accounting_support_work_j);
    close(transformed.reaction.midpoint_work_j,result.reaction.midpoint_work_j);
    audit(&model,&shifted,1,incoming,&transformed);
}

#[test]
fn finite_tangent_defect_has_separate_midpoint_work_without_double_counting_projection() {
    let model = ArticulatedModel::new(vec![Link::new(
        "almost locked",None,Se3::identity(),
        JointModel::prismatic(Vec3::new(1.,0.,0.),None).unwrap(),sphere(3.),
    )]).unwrap();
    let state = PrescribedBaseMotionState {
        time_s: 0.,base: BaseState::stationary(Se3::identity()),
        joint_position: vec![0.],joint_velocity: vec![0.01],
    };
    let payload = PayloadState {
        inertia: sphere(2.),world_from_body: Se3::identity(),
        twist_body: Twist::new(Vec3::default(),Vec3::new(0.01,0.,0.)),
    };
    let mut lock = row(ConstraintCoordinate::Length,vec![1.],0.);
    lock.velocity_tolerance = 0.02;
    let controls = SupportedPayloadControls { energy_tolerance_j: 1e-12,..controls(&[2e-9]) };
    let mut budget = SupportedPayloadConstraintControls {
        maximum_reaction_work_j: 0.0003,..row_controls()
    };
    let result = capture_payload(&model,&state,0,payload,&[lock.clone()],
        controls,budget,||false).unwrap();
    close(result.capture.removed_relative_energy_j,0.00025);
    close(result.reaction.midpoint_work_j,-0.00025);
    close(result.reaction.final_rate_pairing_j,0.);
    assert!(result.capture.balance.energy_residual_j.abs() < 1e-12);
    close(result.reaction.velocity_residuals_before[0],0.01);
    close(result.reaction.velocity_residuals_after[0],0.);
    budget.maximum_reaction_work_j = 0.0002;
    assert!(matches!(capture_payload(&model,&state,0,payload,&[lock],
        controls,budget,||false),Err(PayloadError::ConstraintWork{..})));
    audit(&model,&state,0,payload,&result);
}

#[test]
fn empty_rows_match_original_and_geared_capture_continues_through_native_motion_and_release() {
    let (model,state,payload,rows) = geared_sliders(0.,10.);
    let ordinary = unconstrained_capture(&model,&state,1,payload,controls(&[2e-9;2]),||false).unwrap();
    let empty = capture_payload(&model,&state,1,payload,&[],
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    assert_eq!(ordinary,empty.capture);
    assert!(empty.reaction.joint_impulse.is_empty());
    assert!(empty.reaction.multipliers.is_empty());
    let capture = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
    let next = step_prescribed_base_spatial(
        &capture.model_after,&capture.state_after,PrescribedBasePath::stationary(Se3::identity()),
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.,0.],
                external_body_wrench: &[Wrench::default();3],
                gravity_world: Vec3::default(),
            },world_points: &[],connections: &[],
        },
        &rows,0.01,ImplicitSpatialConfig::default(),
        SupportedConstraintControls {
            minimum_reciprocal_condition: 1e-10,maximum_reaction_work_j: 2e-8,
        },||false,
    ).unwrap();
    let endpoint = &next.implicit.spatial.motion.state_after;
    close(endpoint.time_s,4.01);
    close(endpoint.joint_position[0],0.436);
    close(endpoint.joint_position[1],0.218);
    close(endpoint.joint_velocity[0],3.6);
    close(endpoint.joint_velocity[1],1.8);
    let release = release_payload(&capture.model_after,endpoint,&capture.attachment,
        SupportedPayloadReleaseControls {
            angular_momentum_tolerance_nms: 2e-9,linear_momentum_tolerance_ns: 2e-9,
            joint_momentum_tolerances: &[2e-9;2],energy_tolerance_j: 2e-8,
        },||false).unwrap();
    assert_eq!(release.model_after,model);
    assert_eq!(release.state_after,*endpoint);
    close(release.payload_after.twist_body.linear.x,3.6);
    let ejected =
        fs_mbd::articulated::motion::inventory::supported::constraints::release_payload_with_impulse(
            &capture.model_after,endpoint,&capture.attachment,impulse(4.),&rows,
            separation_controls(&[2e-9;2]),row_controls(),||false,
        ).unwrap();
    close(ejected.release.state_after.time_s,4.01);
    close(ejected.release.state_after.joint_velocity[0],3.6-16./17.);
    close(ejected.release.payload_after.twist_body.linear.x,5.6);
    release_audit(&capture.model_after,endpoint,1,&ejected);
}

#[test]
fn native_speed_row_rank_state_and_budget_refusals_are_atomic_and_cancellable() {
    let (model,state,payload,rows) = geared_sliders(0.,3.);
    assert!(matches!(capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false),
        Err(PayloadError::Dynamics(ArticulatedError::JointVelocityLimitViolation{..}))));
    let (model,state,payload,rows) = geared_sliders(0.,10.);
    let saved = (model.clone(),state.clone(),payload,rows.clone());
    let mut bad_state = state.clone();
    bad_state.joint_velocity[0] += 0.1;
    assert!(matches!(capture_payload(&model,&bad_state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false),
        Err(PayloadError::Constraint(ConstraintError::Velocity{..}))));
    bad_state = state.clone();
    bad_state.joint_position[0] += 0.1;
    assert!(matches!(capture_payload(&model,&bad_state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false),
        Err(PayloadError::Constraint(ConstraintError::Position{..}))));
    let duplicate = vec![rows[0].clone(),rows[0].clone()];
    assert!(matches!(capture_payload(&model,&state,1,payload,&duplicate,
        controls(&[2e-9;2]),row_controls(),||false),Err(PayloadError::Constraint(_))));
    assert!(capture_payload(&model,&state,1,payload,&rows,controls(&[2e-9;2]),
        SupportedPayloadConstraintControls {minimum_reciprocal_condition: f64::NAN,..row_controls()},
        ||false).is_err());
    assert!(capture_payload(&model,&state,1,payload,&rows,controls(&[2e-9;2]),
        SupportedPayloadConstraintControls {maximum_reaction_work_j: -1.,..row_controls()},
        ||false).is_err());
    assert!(matches!(capture_payload(&model,&state,1,payload,&rows,
        SupportedPayloadControls {maximum_removed_relative_energy_j: 16.9,..controls(&[2e-9;2])},
        row_controls(),||false),Err(PayloadError::EnergyExchange{..})));
    let polls = Cell::new(0usize);
    let accepted = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||{polls.set(polls.get()+1);false}).unwrap();
    for stop in 1..=polls.get() {
        let seen = Cell::new(0usize);
        assert_eq!(capture_payload(&model,&state,1,payload,&rows,
            controls(&[2e-9;2]),row_controls(),||{
                seen.set(seen.get()+1);seen.get() >= stop
            }).unwrap_err(),PayloadError::Cancelled);
    }
    let retry = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap();
    assert_eq!(retry,accepted);
    assert_eq!((model,state,payload,rows),saved);
}


fn separation_controls(tolerances: &[f64])
    -> fs_mbd::articulated::motion::inventory::supported::release::SupportedPayloadSeparationControls<'_>
{
    use fs_mbd::articulated::motion::inventory::supported::release::SupportedPayloadSeparationControls;
    SupportedPayloadSeparationControls {
        balance: SupportedPayloadReleaseControls {
            angular_momentum_tolerance_nms: 2e-9,
            linear_momentum_tolerance_ns: 2e-9,
            joint_momentum_tolerances: tolerances,
            energy_tolerance_j: 2e-8,
        },
        maximum_separation_work_j: 1000.,
    }
}

fn impulse(value: f64)
    -> fs_mbd::articulated::motion::inventory::release::SeparationImpulse
{
    fs_mbd::articulated::motion::inventory::release::SeparationImpulse {
        point_link_m: Vec3::default(),
        impulse_on_payload_link_ns: Vec3::new(value,0.,0.),
        couple_on_payload_link_nms: Vec3::default(),
    }
}

fn release_audit(
    loaded: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target: usize,
    result: &fs_mbd::articulated::motion::inventory::supported::constraints::ConstrainedSupportedPayloadRelease,
) {
    let release = &result.release;
    let before = physical(loaded,state,None);
    let after = physical(&release.model_after,&release.state_after,Some((target,release.payload_after)));
    let balance = &release.balance;
    vector(after.linear-before.linear,balance.support_impulse_world.force);
    vector(after.angular-before.angular,balance.support_impulse_world.torque);
    for (i,(&old,&new)) in before.joint.iter().zip(&after.joint).enumerate() {
        let applied = result.reaction.joint_impulse.get(i).copied().unwrap_or(0.);
        close(new-old,applied);
        close(balance.joint_momentum_residual[i],new-old-applied);
    }
    close(balance.kinetic_energy_change_j,after.energy-before.energy);
    close(balance.kinetic_energy_change_j,
        balance.support_work_j+release.supplied_work_j+result.reaction.midpoint_work_j);
    close(balance.energy_residual_j,0.);
    assert_eq!(release.state_after.base,state.base);
    assert_eq!(release.state_after.joint_position,state.joint_position);
    assert_eq!(release.state_after.time_s,state.time_s);
}

#[test]
fn geared_ejection_retains_finite_recoil_support_work_and_native_endpoint_limits() {
    use fs_mbd::articulated::motion::inventory::supported::constraints::release_payload_with_impulse;
    use fs_mbd::articulated::motion::inventory::supported::release::release_payload_with_impulse
        as unconstrained_ejection;
    for support_speed in [0.,2.,-2.] {
        for applied in [4.,-4.] {
            let (model,state,payload,rows) = geared_sliders(support_speed,4.6);
            let capture = capture_payload(&model,&state,1,payload,&rows,
                controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
            if applied < 0. {
                // The unrestrained recoil reaches 4.933... m/s, while the
                // simultaneous gear reaction leaves 4.541... m/s admissible.
                assert!(matches!(unconstrained_ejection(
                    &capture.model_after,&capture.state_after,&capture.attachment,impulse(applied),
                    separation_controls(&[2e-9;2]),||false),
                    Err(PayloadError::Dynamics(ArticulatedError::JointVelocityLimitViolation{..}))));
            }
            let result = release_payload_with_impulse(
                &capture.model_after,&capture.state_after,&capture.attachment,impulse(applied),
                &rows,separation_controls(&[2e-9;2]),row_controls(),||false).unwrap();
            assert_eq!(result.release.model_after,model);
            close(result.release.state_after.joint_velocity[0],3.6-4.*applied/17.);
            close(result.release.state_after.joint_velocity[1],1.8-2.*applied/17.);
            close(result.release.payload_after.twist_body.linear.x,support_speed+3.6+applied/2.);
            close(result.reaction.joint_impulse[0],5.*applied/17.);
            close(result.reaction.joint_impulse[1],-10.*applied/17.);
            close(result.release.balance.support_impulse_body.force.x,-5.*applied/17.);
            close(result.release.balance.support_work_j,-5.*applied*support_speed/17.);
            close(result.release.supplied_work_j,100./17.);
            close(result.reaction.midpoint_work_j,0.);
            release_audit(&capture.model_after,&capture.state_after,1,&result);
        }
    }
}

// Independent three-component centroidal solve, not the production 6x6 map.
fn inverse3(matrix: Mat3, vector: Vec3) -> Vec3 {
    let a = matrix.m;
    let c00 = a[4]*a[8]-a[5]*a[7];
    let c01 = a[5]*a[6]-a[3]*a[8];
    let c02 = a[3]*a[7]-a[4]*a[6];
    let determinant = a[0]*c00+a[1]*c01+a[2]*c02;
    Vec3::new(
        c00*vector.x+(a[2]*a[7]-a[1]*a[8])*vector.y+(a[1]*a[5]-a[2]*a[4])*vector.z,
        c01*vector.x+(a[0]*a[8]-a[2]*a[6])*vector.y+(a[2]*a[3]-a[0]*a[5])*vector.z,
        c02*vector.x+(a[1]*a[6]-a[0]*a[7])*vector.y+(a[0]*a[4]-a[1]*a[3])*vector.z,
    ).scale(1./determinant)
}
fn point_velocity(pose: Se3, twist: Twist, point: Vec3) -> Vec3 {
    pose.rotation().rotate(twist.linear).unwrap()
        + pose.rotation().rotate(twist.angular).unwrap().cross(point-pose.translation())
}
fn oblique_impulse() -> fs_mbd::articulated::motion::inventory::release::SeparationImpulse {
    fs_mbd::articulated::motion::inventory::release::SeparationImpulse {
        point_link_m: Vec3::new(0.2,0.4,-0.1),
        impulse_on_payload_link_ns: Vec3::new(1.2,-0.6,0.8),
        couple_on_payload_link_nms: Vec3::new(0.3,-0.2,0.4),
    }
}

#[test]
fn full_tensor_geared_ejection_matches_reduced_recoil_and_point_couple_work() {
    use fs_mbd::articulated::motion::inventory::supported::constraints::release_payload_with_impulse;
    let (model,state,payload,rows) = branched();
    let capture = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
    let command = oblique_impulse();
    let result = release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,command,
        &rows,separation_controls(&[2e-9;2]),row_controls(),||false).unwrap();
    let mut zero = capture.state_after.clone();
    zero.joint_velocity = vec![0.,0.];
    let mut unit = zero.clone();
    unit.joint_velocity = vec![2.,-1.];
    let p0 = physical(&model,&zero,None);
    let p1 = physical(&model,&unit,None);
    let mass = 2.*(p1.joint[0]-p0.joint[0])-(p1.joint[1]-p0.joint[1]);
    let k0 = forward_kinematics(&model,zero.base,&zero.joint_position,&zero.joint_velocity).unwrap();
    let k1 = forward_kinematics(&model,unit.base,&unit.joint_position,&unit.joint_velocity).unwrap();
    let tangent = k1.body_twist[1].minus(k0.body_twist[1]);
    let j_payload_link = Wrench::new(
        command.point_link_m.cross(command.impulse_on_payload_link_ns)+command.couple_on_payload_link_nms,
        command.impulse_on_payload_link_ns,
    );
    let generalized = -j_payload_link.torque.dot(tangent.angular)-j_payload_link.force.dot(tangent.linear);
    let rate_jump = generalized/mass;
    close(result.release.state_after.joint_velocity[0],capture.state_after.joint_velocity[0]+2.*rate_jump);
    close(result.release.state_after.joint_velocity[1],capture.state_after.joint_velocity[1]-rate_jump);

    let j_payload = j_payload_link.transform_by(&capture.attachment.link_from_payload().inverse().unwrap());
    let inertia = payload.inertia;
    let spin_jump = inverse3(inertia.inertia_com(),j_payload.torque-inertia.center_of_mass().cross(j_payload.force));
    let linear_jump = j_payload.force.scale(1./inertia.mass())-spin_jump.cross(inertia.center_of_mass());
    vector(result.release.payload_after.twist_body.angular,
        capture.payload_after.twist_body.angular+spin_jump);
    vector(result.release.payload_after.twist_body.linear,
        capture.payload_after.twist_body.linear+linear_jump);

    let before = forward_kinematics(&capture.model_after,capture.state_after.base,
        &capture.state_after.joint_position,&capture.state_after.joint_velocity).unwrap();
    let after = forward_kinematics(&model,result.release.state_after.base,
        &result.release.state_after.joint_position,&result.release.state_after.joint_velocity).unwrap();
    let target_pose = before.world_from_link[1];
    let point = target_pose.transform_point(command.point_link_m).unwrap();
    let force = target_pose.rotation().rotate(command.impulse_on_payload_link_ns).unwrap();
    let couple = target_pose.rotation().rotate(command.couple_on_payload_link_nms).unwrap();
    let ret_point_jump = point_velocity(target_pose,after.body_twist[1],point)
        -point_velocity(target_pose,before.body_twist[1],point);
    let body_pose = capture.payload_after.world_from_body;
    let payload_point_jump = point_velocity(body_pose,result.release.payload_after.twist_body,point)
        -point_velocity(body_pose,capture.payload_after.twist_body,point);
    let ret_spin_jump = target_pose.rotation().rotate(after.body_twist[1].angular-before.body_twist[1].angular).unwrap();
    let body_spin_jump = body_pose.rotation().rotate(spin_jump).unwrap();
    close(result.release.supplied_work_j,
        0.5*(force.dot(payload_point_jump-ret_point_jump)+couple.dot(body_spin_jump-ret_spin_jump)));
    release_audit(&capture.model_after,&capture.state_after,1,&result);
}

#[test]
fn passive_rows_preserve_state_and_active_finite_reaction_work_is_separate() {
    use fs_mbd::articulated::motion::inventory::supported::constraints::{
        release_payload as constrained_passive,release_payload_with_impulse,
    };
    let (model,state,payload,mut rows) = geared_sliders(0.,10.);
    let capture = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
    let mut imperfect = capture.state_after.clone();
    imperfect.joint_velocity[0] += 0.001;
    rows[0].velocity_tolerance = 0.002;
    let control = separation_controls(&[2e-9;2]);
    let passive = release_payload(&capture.model_after,&imperfect,&capture.attachment,
        control.balance,||false).unwrap();
    let preserved = constrained_passive(&capture.model_after,&imperfect,&capture.attachment,
        &rows,control.balance,row_controls(),||false).unwrap();
    assert_eq!(preserved.release,passive);
    assert_eq!(preserved.release.state_after,imperfect);
    close(preserved.reaction.midpoint_work_j,0.);
    close(preserved.reaction.final_rate_pairing_j,0.);
    assert!(preserved.reaction.joint_impulse.iter().all(|x| *x == 0.));
    let zero = release_payload_with_impulse(&capture.model_after,&imperfect,&capture.attachment,
        impulse(0.),&rows,control,row_controls(),||false).unwrap();
    assert_eq!(zero.release.state_after,passive.state_after);
    assert_eq!(zero.release.payload_after,passive.payload_after);
    assert_eq!(zero.release.balance,passive.balance);
    let budget = SupportedPayloadConstraintControls {maximum_reaction_work_j: 0.001,..row_controls()};
    let driven = release_payload_with_impulse(&capture.model_after,&imperfect,&capture.attachment,
        impulse(4.),&rows,control,budget,||false).unwrap();
    close(driven.reaction.midpoint_work_j,(20.-0.015)/17.*0.0005);
    close(driven.release.supplied_work_j,100./17.+0.01/17.);
    release_audit(&capture.model_after,&imperfect,1,&driven);
    assert!(matches!(release_payload_with_impulse(
        &capture.model_after,&imperfect,&capture.attachment,impulse(4.),&rows,
        control,row_controls(),||false),Err(PayloadError::ConstraintWork{..})));
    // Existing source-only semantics refuse a negative pair work, which can
    // arise when a tiny command opposes a finite admitted row correction.
    assert!(matches!(release_payload_with_impulse(
        &capture.model_after,&imperfect,&capture.attachment,impulse(-0.0001),&rows,
        control,budget,||false),Err(PayloadError::Invalid(_))));
}

#[test]
fn constrained_ejection_preserves_rigid_frame_and_boost_covariance() {
    use fs_mbd::articulated::motion::inventory::supported::constraints::release_payload_with_impulse;
    let (model,state,payload,rows) = branched();
    let capture = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
    let command = oblique_impulse();
    let original = release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,command,&rows,
        separation_controls(&[2e-9;2]),row_controls(),||false).unwrap();
    let change = pose(Vec3::new(-0.4,0.6,0.3),Vec3::new(-5.,3.,7.));
    let boost = Vec3::new(17.,-9.,6.);
    let old_state = &capture.state_after;
    let new_pose = change.compose(old_state.base.world_from_base).unwrap();
    let local_boost = new_pose.rotation().inverse().rotate(boost).unwrap();
    let shifted = PrescribedBaseMotionState {
        base: BaseState::prescribed(new_pose,
            old_state.base.twist_body.plus(Twist::new(Vec3::default(),local_boost)),
            Twist::new(old_state.base.acceleration_body.angular,
                old_state.base.acceleration_body.linear-old_state.base.twist_body.angular.cross(local_boost))),
        ..old_state.clone()
    };
    let transformed = release_payload_with_impulse(
        &capture.model_after,&shifted,&capture.attachment,command,&rows,
        separation_controls(&[2e-9;2]),row_controls(),||false).unwrap();
    let a = &original.release;
    let b = &transformed.release;
    for i in 0..2 {
        close(a.state_after.joint_velocity[i],b.state_after.joint_velocity[i]);
        close(original.reaction.joint_impulse[i],transformed.reaction.joint_impulse[i]);
    }
    close(a.supplied_work_j,b.supplied_work_j);
    close(original.reaction.midpoint_work_j,transformed.reaction.midpoint_work_j);
    vector(a.balance.support_impulse_body.force,b.balance.support_impulse_body.force);
    vector(a.balance.support_impulse_body.torque,b.balance.support_impulse_body.torque);
    let force = change.rotation().rotate(a.balance.support_impulse_world.force).unwrap();
    let torque = change.rotation().rotate(a.balance.support_impulse_world.torque).unwrap()
        +change.translation().cross(force);
    vector(b.balance.support_impulse_world.force,force);
    vector(b.balance.support_impulse_world.torque,torque);
    close(b.balance.support_work_j,a.balance.support_work_j+boost.dot(force));
    close(b.balance.kinetic_energy_change_j,a.balance.kinetic_energy_change_j+boost.dot(force));
    close(b.balance.kinetic_before_j,a.balance.kinetic_before_j);
    close(b.balance.kinetic_after_j,a.balance.kinetic_after_j);
    release_audit(&capture.model_after,&shifted,1,&transformed);
}

#[test]
fn constrained_ejection_refuses_source_speed_rank_and_cancellation_atomically() {
    use fs_mbd::articulated::motion::inventory::supported::constraints::{
        release_payload as constrained_passive,release_payload_with_impulse,
    };
    let (model,state,payload,rows) = geared_sliders(0.,4.3);
    let capture = capture_payload(&model,&state,1,payload,&rows,
        controls(&[2e-9;2]),row_controls(),||false).unwrap().capture;
    let saved = capture.clone();
    assert!(matches!(release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,impulse(-4.),&rows,
        separation_controls(&[2e-9;2]),row_controls(),||false),
        Err(PayloadError::Dynamics(ArticulatedError::JointVelocityLimitViolation{..}))));
    let mut limited = separation_controls(&[2e-9;2]);
    limited.maximum_separation_work_j = 5.8;
    assert!(matches!(release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,impulse(4.),&rows,
        limited,row_controls(),||false),Err(PayloadError::EnergyExchange{..})));
    let duplicate = [rows[0].clone(),rows[0].clone()];
    assert!(matches!(constrained_passive(
        &capture.model_after,&capture.state_after,&capture.attachment,&duplicate,
        limited.balance,row_controls(),||false),Err(PayloadError::Constraint(_))));
    let polls = Cell::new(0usize);
    let good = release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,impulse(4.),&rows,
        separation_controls(&[2e-9;2]),row_controls(),||{polls.set(polls.get()+1);false}).unwrap();
    for stop in 1..=polls.get() {
        let seen = Cell::new(0usize);
        assert_eq!(release_payload_with_impulse(
            &capture.model_after,&capture.state_after,&capture.attachment,impulse(4.),&rows,
            separation_controls(&[2e-9;2]),row_controls(),||{
                seen.set(seen.get()+1);seen.get()>=stop
            }).unwrap_err(),PayloadError::Cancelled);
    }
    let retry = release_payload_with_impulse(
        &capture.model_after,&capture.state_after,&capture.attachment,impulse(4.),&rows,
        separation_controls(&[2e-9;2]),row_controls(),||false).unwrap();
    assert_eq!(good,retry);
    assert_eq!(capture,saved);
}

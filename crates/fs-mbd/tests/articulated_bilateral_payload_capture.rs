//! G0/G1/G3/G4: payload inventory and mixed rigid mechanism constraints share
//! the original loaded impulse response and independent physical accounts.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::{self, PayloadControls, PayloadError, PayloadState};
use fs_mbd::articulated::motion::inventory::bilateral::{
    BilateralPayloadCapture, capture_payload as capture,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop, BilateralLoopError, DistanceLoop,
    PointCoincidenceLoop, joint,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialError;
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-9 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.16e} != {expected:.16e}");
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn sphere(mass: f64) -> SpatialInertia {
    SpatialInertia::new(mass, Vec3::default(), Mat3 {
        m: [mass, 0., 0., 0., mass, 0., 0., 0., mass],
    }).unwrap()
}
fn translated(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0., 0.)).unwrap()
}
fn point(link: usize) -> MaterialPoint {
    MaterialPoint { link, point_body: Vec3::default() }
}
fn rows() -> [JointConstraint; 1] {
    [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }]
}
fn loops() -> [BilateralLoop; 1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: point(1), b: point(2), length_m: 2.,
    })]
}
fn model(maximum_speed: Option<f64>) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier", None, Se3::identity(), JointModel::FIXED, sphere(5.),
    )];
    for (i, mass) in [2., 3., 4.].into_iter().enumerate() {
        links.push(Link::new(
            format!("slider {i}"), Some(0), Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1., 0., 0.),
                maximum_speed.map(|speed| JointLimits::new(-10., 10., speed, 0.1).unwrap()),
            ).unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(), Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
    }
}
fn incoming() -> PayloadState {
    PayloadState {
        inertia: sphere(2.),
        world_from_body: translated(1.),
        twist_body: Twist::new(Vec3::default(), Vec3::new(2.3, 0., 0.)),
    }
}
fn controls(tolerances: &[f64]) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_energy_exchange_j: 20.,
        maximum_constraint_work_j: 1e-8,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn solve(model: &ArticulatedModel, state: &FreeFloatingMotionState, payload: PayloadState) -> BilateralPayloadCapture {
    capture(
        model, state, 1, payload, &loops(), &rows(), controls(&[2e-9; 3]),
        Default::default(), || false,
    ).unwrap()
}

#[test]
fn geared_closed_capture_transmits_payload_momentum_and_continues_loaded_motion() {
    let model = model(Some(0.8));
    let start = initial();
    let saved = start.clone();
    let payload = incoming();
    // The free catcher reaches relative speed1.0; the physically combined
    // endpoint below reaches0.75745 and therefore must be admitted.
    assert!(inventory::capture_payload(
        &model, &start, 1, payload, &rows(), controls(&[2e-9; 3]), || false,
    ).is_err());
    let result = solve(&model, &start, payload);
    let capture = &result.capture;
    let event = result.reaction.as_ref().unwrap();
    // Independent closed mechanism coordinates are carrier V and mode u,
    // qdot=(2,2,1)*u. Loaded mass matrix is [[16,18],[18,32]].
    let momentum = 17.2;
    let conjugate = 23.8;
    let determinant = 16.*32.-18.*18.;
    let velocity = (32.*momentum-18.*conjugate)/determinant;
    let mode = (16.*conjugate-18.*momentum)/determinant;
    close(capture.state_after.base.twist_body.linear.x, velocity);
    for (rate, factor) in capture.state_after.joint_velocity.iter().zip([2.,2.,1.]) {
        close(*rate, factor*mode);
    }
    assert_eq!(capture.state_after.joint_position, start.joint_position);
    assert_eq!(capture.state_after.base.world_from_base, start.base.world_from_base);
    close(capture.model_after.links()[1].inertia().mass(), 4.);
    close(capture.payload_after.twist_body.linear.x, velocity+2.*mode);
    close(capture.payload_after.inertia.mass(), 2.);
    assert_eq!(capture.attachment.target_link(), 1);
    assert_eq!(capture.attachment.retained_inertia(), model.links()[1].inertia());
    let body_multiplier = 6.8-4.*(velocity+2.*mode);
    let joint_multiplier = (3.6-4.*(velocity+mode))/2.;
    close(event.loop_multipliers_ns[0], body_multiplier);
    close(event.joint_multipliers[0], joint_multiplier);
    close(capture.balance.generalized_constraint_impulse[6], -body_multiplier);
    close(capture.balance.generalized_constraint_impulse[7], body_multiplier+joint_multiplier);
    close(capture.balance.generalized_constraint_impulse[8], -2.*joint_multiplier);
    let kinetic_before = 11.16;
    let kinetic_after = 0.5*(momentum*velocity+conjugate*mode);
    close(capture.balance.momentum_before.kinetic_energy, kinetic_before);
    close(capture.balance.momentum_after.kinetic_energy, kinetic_after);
    close(capture.removed_relative_energy_j, kinetic_before-kinetic_after);
    close(capture.balance.reaction_work_j, 0.);
    close(capture.balance.energy_residual_j, 0.);
    vector(capture.balance.momentum_after.linear_momentum_world, Vec3::new(momentum,0.,0.));
    vector(capture.balance.momentum_after.angular_momentum_world, Vec3::default());
    assert!(capture.balance.maximum_normalized_momentum_residual <= 1.);
    // The embedded event starts from a loaded, rigidly carried baseline;
    // its source work is not the incoming body's complete capture energy.
    close(event.impulse.momentum_before.kinetic_energy, 7.08);
    close(event.source_work_j, 2.4*0.5*(1.1+velocity+2.*mode));
    close(event.constraint_work_j, 0.);
    let correlation = (1./3.) / ((7.0_f64/12.)*(23./15.)).sqrt();
    close(event.reciprocal_condition, (1.-correlation)/(1.+correlation));

    // Continue the actual loaded inventory through the existing mixed-loop
    // motion owner. The retained effective mode inertia is47/4, not the old10.
    let wrenches = [Wrench::default();4];
    let motion = joint::step_with_bushings(
        &capture.model_after, &capture.state_after,
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.05,0.,0.],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::default(),
            },
            connections: &[], world_points: &[],
        },
        &loops(), &rows(), &[], 0.1, Default::default(), || false,
    ).unwrap();
    let after = &motion.bilateral.state_after;
    let acceleration = 0.4/47.;
    close(after.joint_velocity[2], mode+0.1*acceleration);
    close(after.base.twist_body.linear.x, velocity-0.1*(18./16.)*acceleration);
    close(after.joint_position[0], start.joint_position[0]+0.2*(mode+0.05*acceleration));
    vector(motion.bilateral.linear_momentum_residual_world, Vec3::default());
    close(motion.bilateral.mechanical_work_residual, 0.);
    assert_eq!(start, saved);
}

#[test]
fn closed_capture_preserves_frame_covariance_and_small_energy_under_a_large_boost() {
    let model = model(None);
    let reference = solve(&model, &initial(), incoming());
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3,-0.2,0.4))).unwrap();
    let frame = Se3::from_parts(rotation, Vec3::new(1.2,-0.7,0.9)).unwrap();
    let boost = Vec3::new(-0.5,0.3,0.2);
    let mut state = initial();
    state.base.world_from_base = frame;
    state.base.twist_body.linear = state.base.twist_body.linear + rotation.inverse().rotate(boost).unwrap();
    let mut payload = incoming();
    payload.world_from_body = frame.compose(payload.world_from_body).unwrap();
    payload.twist_body.linear = payload.twist_body.linear + rotation.inverse().rotate(boost).unwrap();
    let moved = solve(&model, &state, payload);
    for (a,b) in moved.capture.state_after.joint_velocity.iter().zip(&reference.capture.state_after.joint_velocity) {
        close(*a,*b);
    }
    vector(
        moved.capture.state_after.base.twist_body.linear,
        reference.capture.state_after.base.twist_body.linear + rotation.inverse().rotate(boost).unwrap(),
    );
    close(moved.capture.removed_relative_energy_j, reference.capture.removed_relative_energy_j);
    close(moved.capture.balance.reaction_work_j, reference.capture.balance.reaction_work_j);
    vector(moved.capture.balance.momentum_after.linear_momentum_world,
        rotation.rotate(reference.capture.balance.momentum_after.linear_momentum_world).unwrap()
        + boost.scale(16.));

    let mut state = initial();
    let mut payload = incoming();
    state.base.twist_body.linear.x += 1e9;
    payload.twist_body.linear.x += 1e9;
    let tolerance = [2e-7;3];
    let controls = PayloadControls {
        energy_tolerance_j: 3e-6,
        maximum_constraint_work_j: 2e-6,
        ..controls(&tolerance)
    };
    let result = capture(
        &model, &state, 1, payload, &loops(), &rows(), controls,
        Default::default(), || false,
    ).unwrap();
    // Build the independent solve from the representable incoming velocities
    // in the same translating frame; avoid subtracting billion-scale energies.
    let relative = payload.twist_body.linear.x-state.base.twist_body.linear.x;
    let p = 14.*0.2+2.*relative;
    let h = 24.*0.2+4.*relative;
    let dv = (32.*p-18.*h)/188.;
    let mode = (16.*h-18.*p)/188.;
    let energy_before = 0.5*24.*0.2_f64.powi(2)+relative*relative;
    let energy_after = 0.5*(p*dv+h*mode);
    near(result.capture.state_after.base.twist_body.linear.x-state.base.twist_body.linear.x, dv, 3e-7);
    near(result.capture.state_after.joint_velocity[2], mode, 3e-7);
    near(result.capture.removed_relative_energy_j, energy_before-energy_after, 3e-6);
    assert!(result.capture.removed_relative_energy_j > 1.);
    assert!(result.capture.balance.energy_residual_j.abs() <= controls.energy_tolerance_j);
}

fn assembled(parts: &[SpatialInertia]) -> SpatialInertia {
    let mass: f64 = parts.iter().map(|p| p.mass()).sum();
    let center = parts.iter().fold(Vec3::default(), |sum,p|
        sum+p.center_of_mass().scale(p.mass()/mass));
    let mut tensor = [0.;9];
    for part in parts {
        let d = part.center_of_mass()-center;
        let x = [d.x,d.y,d.z];
        for i in 0..3 {
            for j in 0..3 {
                tensor[3*i+j] += part.inertia_com().m[3*i+j]
                    + part.mass()*((if i==j { d.dot(d) } else { 0. })-x[i]*x[j]);
            }
        }
    }
    SpatialInertia::new(mass,center,Mat3 { m: tensor }).unwrap()
}

#[test]
fn spherical_capture_and_angular_lock_match_independent_full_tensor_assembly() {
    let parts: Vec<_> = [5.,1.,1.,1.,2.].into_iter().enumerate().map(|(i,mass)| {
        let s = i as f64;
        SpatialInertia::new(mass,Vec3::new(0.03*s,-0.02*s,0.01*s),Mat3 {
            m: [mass,0.1*mass,-0.05*mass,0.1*mass,1.2*mass,0.08*mass,
                -0.05*mass,0.08*mass,1.5*mass],
        }).unwrap()
    }).collect();
    let axes = [Vec3::new(1.,0.,0.),Vec3::new(0.,1.,0.),Vec3::new(0.,0.,1.)];
    let mut links = vec![Link::new("base",None,Se3::identity(),JointModel::FIXED,parts[0])];
    for i in 0..3 {
        links.push(Link::new(
            format!("slider {i}"),Some(i),Se3::identity(),
            JointModel::prismatic(axes[i],None).unwrap(),parts[i+1],
        ));
    }
    links.push(Link::new(
        "locked rotor",Some(3),Se3::identity(),
        JointModel::revolute(axes[2],None).unwrap(),parts[4],
    ));
    let model = ArticulatedModel::new(links).unwrap();
    let reduced = ArticulatedModel::new(vec![Link::new(
        "assembled",None,Se3::identity(),JointModel::FIXED,assembled(&parts),
    )]).unwrap();
    let base = FreeFloatingBaseState::new(
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.3,-0.2,0.1))).unwrap(),
            Vec3::new(0.5,-0.3,0.2),
        ).unwrap(),
        Twist::new(Vec3::new(0.2,-0.15,0.1),Vec3::new(0.3,-0.2,0.4)),
    );
    let state = FreeFloatingMotionState {
        base,joint_position:vec![0.;4],joint_velocity:vec![0.;4],
    };
    let rigid = FreeFloatingMotionState { base,joint_position:vec![],joint_velocity:vec![] };
    let attachment = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.2,0.1,-0.3))).unwrap(),
        Vec3::new(0.4,-0.2,0.1),
    ).unwrap();
    let payload = PayloadState {
        inertia: SpatialInertia::new(
            1.3,Vec3::new(0.1,-0.03,0.05),
            Mat3 { m:[0.7,0.09,-0.02,0.09,0.8,0.04,-0.02,0.04,1.1] },
        ).unwrap(),
        world_from_body: base.world_from_base.compose(attachment).unwrap(),
        twist_body: Twist::new(Vec3::new(0.4,0.1,-0.3),Vec3::new(1.2,-0.3,0.7)),
    };
    let loops = [BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a:point(0),b:point(3),
    })];
    let rows = [JointConstraint {
        coordinate:ConstraintCoordinate::Angle,coefficients:vec![0.,0.,0.,1.],
        position:0.,position_tolerance:1e-9,velocity_tolerance:1e-9,
    }];
    let actual = capture(
        &model,&state,4,payload,&loops,&rows,controls(&[2e-8;4]),
        Default::default(),|| false,
    ).unwrap();
    let expected = inventory::capture_payload(
        &reduced,&rigid,0,payload,&[],controls(&[]),|| false,
    ).unwrap();
    vector(actual.capture.state_after.base.twist_body.angular,expected.state_after.base.twist_body.angular);
    vector(actual.capture.state_after.base.twist_body.linear,expected.state_after.base.twist_body.linear);
    for value in &actual.capture.state_after.joint_velocity { close(*value,0.); }
    close(actual.capture.removed_relative_energy_j,expected.removed_relative_energy_j);
    // Each model chooses its own root-COM accounting boost. Compare their
    // common lab-frame kinetic energy, not those different frame coordinates.
    close(actual.capture.balance.momentum_after.kinetic_energy,
        expected.balance.momentum_after.kinetic_energy);
    vector(actual.capture.balance.momentum_after.linear_momentum_world,expected.balance.momentum_after.linear_momentum_world);
    vector(actual.capture.balance.momentum_after.angular_momentum_world,expected.balance.momentum_after.angular_momentum_world);
    vector(actual.capture.payload_after.twist_body.angular,expected.payload_after.twist_body.angular);
    vector(actual.capture.payload_after.twist_body.linear,expected.payload_after.twist_body.linear);
    let reaction = actual.reaction.as_ref().unwrap();
    assert_eq!(reaction.loop_multipliers_ns.len(),3);
    assert_eq!(reaction.joint_multipliers.len(),1);
    assert!(actual.capture.balance.generalized_constraint_impulse[6..].iter().any(|j| j.abs()>1e-6));
    close(actual.capture.balance.energy_residual_j,0.);
}

#[test]
fn finite_tangency_work_is_retained_once_and_exact_carried_capture_stays_passive() {
    let model = model(None);
    let mut state = initial();
    let defect = 1e-3;
    state.joint_velocity[2] -= 0.5*defect;
    let rows = [JointConstraint { velocity_tolerance:2e-3,..rows()[0].clone() }];
    let controls = PayloadControls {
        maximum_constraint_work_j:1e-3,..controls(&[2e-9;3])
    };
    let result = capture(
        &model,&state,1,incoming(),&loops(),&rows,controls,
        Default::default(),|| false,
    ).unwrap();
    let lambda = -12./47.-35.*defect/47.;
    let work = 0.5*lambda*defect;
    close(result.reaction.as_ref().unwrap().joint_multipliers[0],lambda);
    close(result.capture.balance.reaction_work_j,work);
    assert!(work < -1e-4);
    close(result.capture.balance.kinetic_after_j-result.capture.balance.kinetic_before_j,
        -result.capture.removed_relative_energy_j);
    close(result.capture.balance.energy_residual_j,0.);
    assert!(matches!(
        capture(
            &model,&state,1,incoming(),&loops(),&rows,
            PayloadControls { maximum_constraint_work_j:work.abs()/2.,..controls },
            Default::default(),|| false,
        ),
        Err(PayloadError::ConstraintWork { .. })
    ));

    let carried = PayloadState {
        twist_body:Twist::new(Vec3::default(),Vec3::new(
            state.base.twist_body.linear.x+state.joint_velocity[0],0.,0.,
        )),
        ..incoming()
    };
    let passive = capture(
        &model,&state,1,carried,&loops(),&rows,controls,
        Default::default(),|| false,
    ).unwrap();
    assert_eq!(passive.capture.state_after,state);
    assert_eq!(passive.capture.removed_relative_energy_j,0.);
    assert_eq!(passive.capture.balance.reaction_work_j,0.);
    let reaction = passive.reaction.unwrap();
    assert!(reaction.loop_multipliers_ns.iter().all(|j| *j==0.));
    assert!(reaction.joint_multipliers.iter().all(|j| *j==0.));
    close(passive.capture.state_after.joint_velocity[1]-2.*passive.capture.state_after.joint_velocity[2],defect);
}

#[test]
fn exact_empty_loop_delegation_and_native_inventory_refusals_are_atomic() {
    let limited = model(Some(0.7));
    let model = model(None);
    let state = initial();
    let payload = incoming();
    let config = BilateralImpulseConfig::default();
    let tolerances = [2e-9;3];
    let controls = controls(&tolerances);
    let legacy = inventory::capture_payload(
        &model,&state,1,payload,&rows(),controls,|| false,
    ).unwrap();
    let delegated = capture(
        &model,&state,1,payload,&[],&rows(),controls,config,|| false,
    ).unwrap();
    assert_eq!(delegated.capture,legacy);
    assert!(delegated.reaction.is_none());
    let reference = solve(&model,&state,payload);
    let condition = reference.reaction.as_ref().unwrap().reciprocal_condition;
    assert!(condition>0. && condition<1.);
    assert!(matches!(
        capture(
            &model,&state,1,payload,&loops(),&rows(),
            PayloadControls { minimum_reciprocal_condition:0.5*(1.+condition),..controls },
            config,|| false,
        ),
        Err(PayloadError::Constraint(ConstraintError::IllConditioned { .. }))
    ));
    let redundant = [JointConstraint {
        coordinate:ConstraintCoordinate::Length,coefficients:vec![-1.,1.,0.],
        position:2.,position_tolerance:1e-9,velocity_tolerance:1e-9,
    }];
    assert!(matches!(
        capture(&model,&state,1,payload,&loops(),&redundant,controls,config,|| false),
        Err(PayloadError::Bilateral(BilateralLoopError::Mobility { .. }))
    ));
    assert!(matches!(
        capture(
            &model,&state,1,payload,&loops(),&rows(),
            PayloadControls { maximum_energy_exchange_j:0.5,..controls },
            config,|| false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert!(matches!(
        capture(&limited,&state,1,payload,&loops(),&rows(),controls,config,|| false),
        Err(PayloadError::Bilateral(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::JointVelocityLimitViolation { .. })
        ))))
    ));
    let polls = Cell::new(0);
    let run = capture(
        &model,&state,1,payload,&loops(),&rows(),controls,config,
        || { polls.set(polls.get()+1); false },
    ).unwrap();
    let total = polls.get();
    assert!(total>4);
    for threshold in [1,total] {
        let calls = Cell::new(0);
        assert_eq!(
            capture(
                &model,&state,1,payload,&loops(),&rows(),controls,config,
                || { calls.set(calls.get()+1); calls.get()>=threshold },
            ),
            Err(PayloadError::Cancelled),
        );
    }
    assert_eq!(state,initial());
    assert_eq!(solve(&model,&state,payload),run);
}

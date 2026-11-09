//! Supported helical motion: independent electromechanical/contact work oracles.
use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_flyer::partitioned::articulated::actuation::contact::step_with_contacts;
use fs_flyer::partitioned::articulated::contact::{
    ArticulatedFrictionalContact, ContactPath, ContactStepConfig,
};
use fs_mbd::articulated::contact_geometry::{ContactCounterpart, LinkContactPoint};

const MB: f64 = 4.0;
const ML: f64 = 2.0;
const IB: f64 = 4.0;
const IL: f64 = 1.0;
const LEAD: f64 = 0.4;
const G: f64 = 9.81;
const GAP: f64 = 1e-6;

fn mechanism() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new("carrier", None, Se3::identity(), JointModel::FIXED,
            inertia(MB, [0.0;3], [3.0,3.0,IB])),
        Link::new("supported screw", Some(0), Se3::identity(),
            JointModel::helical(vector([0.0,0.0,1.0]), LEAD, None).unwrap(),
            inertia(ML, [0.0;3], [0.8,0.8,IL])),
    ]).unwrap()
}
fn initial(current: f64, rate: f64) -> VoltageDrivenState {
    let mut state = single_state(current, rate);
    state.mechanical.motion.base.world_from_base =
        pose([0.0;3], [0.0,0.0,GAP-LEAD*state.mechanical.motion.joint_position[0]]);
    state.mechanical.motion.base.twist_body.linear.z = -LEAD*rate;
    state
}
fn contact(time: f64, conveyor: f64, friction: f64) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 17,
        point: LinkContactPoint { system_index: 0, link_index: 1, point_body: Vec3::default() },
        other: ContactCounterpart::PrescribedPlane {
            point_world: vector([conveyor*time,0.0,0.0]),
            velocity_world: vector([conveyor,0.0,0.0]),
        },
        normal_world: vector([0.0,0.0,1.0]),
        tangent_world: vector([1.0,0.0,0.0]),
        activation_tolerance: 2.0*GAP,
        friction_coefficient: friction,
        normal_law: NormalVelocityLaw::Target { speed: 0.0 },
    }
}
fn contact_config() -> ContactStepConfig<'static> {
    ContactStepConfig { smooth: one_config(), contact: FrictionalContactConfig::default() }
}
fn geometry(rho: f64) -> Vec<(usize, Strip)> {
    if rho == 0.0 { Vec::new() } else {
        vec![(0, plate(1.0,2.0,[0.0;3],[0.0,0.0,1.0],[0.0,1.0,0.0]))]
    }
}
fn supported_loads(map: &GeneralizedKinematics, strips: &[LinkedStrip], rho: f64,
    external_x: f64) -> Result<ArticulatedLoads, Refusal> {
    let mut result = loads(map, strips, if rho == 0.0 {1.0} else {rho})?;
    result.gravity_world.z = -G;
    result.world_wrenches[1].force.x = external_x;
    Ok(result)
}
fn refused(code: &'static str) -> Refusal {
    Refusal { code, message: code.into(), ranked_repairs: Vec::new() }
}

// This fixture's material point is x(t)=x_base(t)+R_base(t)*LEAD*q(t)*ez.
// Differentiate that original constant-twist/linear-q path twice. The body-frame
// acceleration is affine before its rotation, giving a complete interval bound.
// No point sampling is used to admit coverage or nonpenetration.
fn admit_path(path: &ContactPath<'_>) -> Result<f64, Refusal> {
    let pose = path.initial.motion.base.world_from_base;
    let u = path.midpoint.motion.base.twist_body;
    let z = vector([0.0,0.0,1.0]);
    let r = z.scale(LEAD*path.initial.motion.joint_position[0]);
    let slide = z.scale(LEAD*path.midpoint.motion.joint_velocity[0]);
    let point = pose.transform_point(r).unwrap();
    let speed = pose.rotation().rotate(u.linear + u.angular.cross(r) + slide).unwrap();
    let constant = u.angular.cross(u.linear)
        + u.angular.cross(u.angular.cross(r)) + u.angular.cross(slide).scale(2.0);
    let linear = u.angular.cross(u.angular.cross(slide));
    let acceleration = constant.norm() + path.duration_s*linear.norm();
    let c = path.contacts[0];
    assert_eq!(c.point.link_index,1);
    assert_eq!(c.point.point_body,Vec3::default());
    let ContactCounterpart::PrescribedPlane { point_world, velocity_world } = c.other else {
        panic!("declared plane fixture")
    };
    let gap = c.normal_world.dot(point-point_world);
    let linear_end = gap+path.duration_s*c.normal_world.dot(speed-velocity_world);
    let bound = gap.min(linear_end)-0.5*acceleration*path.duration_s*path.duration_s;
    if bound < 0.0 { return Err(refused("voltage-contact-path")); }
    Ok(bound)
}

struct Oracle { current: f64, rate: f64, endpoint_rate: f64, normal: f64 }
fn oracle(state: &VoltageDrivenState, motor: JointVoltageDrive, passive: JointImpedance,
    rho: f64, dt: f64) -> Oracle {
    // The child COM is supported: carrier translation is -LEAD*q_dot. Angular
    // momentum gives carrier spin=-IL/(IB+IL)*q_dot. No native mass/response
    // query supplies this independent constrained inertia or circuit equation.
    let carrier_mass = MB+rho*std::f64::consts::PI/2.0;
    let inertia = IB*IL/(IB+IL)+carrier_mass*LEAD*LEAD;
    let q = state.mechanical.motion.joint_position[0];
    let w = state.mechanical.motion.joint_velocity[0];
    let h = dt/2.0;
    let k = motor.motor.drive.torque_constant_nm_per_a;
    let l = motor.motor.inductance_h;
    let r = motor.motor.drive.winding_resistance_ohm;
    let a = inertia+h*passive.damping+h*h*passive.stiffness;
    let free = (inertia*w+h*(MB*G*LEAD-passive.stiffness*(q-passive.rest_position)))/a;
    let current = (l*state.currents_a[0]+h*motor.command_voltage_v-h*k*free)
        /(l+h*r+h*h*k*k/a);
    let rate = free+h*k*current/a;
    let endpoint_rate = 2.0*rate-w;
    let normal = (MB+ML)*G-carrier_mass*LEAD*(endpoint_rate-w)/dt;
    assert!(normal > 0.0);
    Oracle { current, rate, endpoint_rate, normal }
}
fn storage(model: &ArticulatedModel, state: &ArticulatedState, plates: &[(usize,Strip)],
    rho: f64, passive: JointImpedance) -> f64 {
    let kinetic = complete_invariants(model,state,plates,rho).0;
    let m = &state.motion;
    let kin = forward_kinematics(model,
        BaseState::prescribed(m.base.world_from_base,m.base.twist_body,Twist::zero()),
        &m.joint_position,&m.joint_velocity).unwrap();
    kinetic+G*(MB*kin.world_from_link[0].translation().z+ML*kin.world_from_link[1].translation().z)
        +0.5*passive.stiffness*(m.joint_position[0]-passive.rest_position).powi(2)
}

#[test]
fn supported_voltage_motion_uses_constrained_dry_fluid_inertia_and_complete_energy() {
    let model = mechanism();
    let passive = JointImpedance { damping: 0.7, stiffness: 3.0, rest_position: 0.1 };
    let dt = 0.02;
    let mut powered_first_rate = Vec::new();
    for rho in [0.0,1.2] {
        let plates = geometry(rho);
        let strips: Vec<_> = plates.iter().map(|(i,p)|LinkedStrip::new(*i,p.clone()).unwrap()).collect();
        // Powered descent and gravity-driven regeneration through the same real
        // support; winding state is continued across every accepted contact tick.
        for (current,rate,voltage) in [(0.2,0.0,1.0),(-0.3,2.0,0.2)] {
            let mut state = initial(current,rate);
            let motor = drive(0,voltage,0.2,0.4);
            for tick in 0..8 {
                let expected = oracle(&state,motor,passive,rho,dt);
                let original = state.clone();
                let before = storage(&model,&state.mechanical,&plates,rho,passive);
                let momentum_before = complete_invariants(&model,&state.mechanical,&plates,rho);
                let calls = Cell::new(0);
                let out = step_with_contacts(&model,&[passive],&state,&[],&[motor],
                    tick as f64*dt,dt,&[contact(tick as f64*dt,0.0,0.3)],contact_config(),
                    |_,_,map| { calls.set(calls.get()+1); supported_loads(map,&strips,rho,0.8) },
                    admit_path,||false).unwrap();
                assert_eq!(calls.get(),out.mechanical.step.iterations);
                let end = &out.mechanical.step.state;
                let receipt = &out.drives[0];
                close(end.motion.joint_velocity[0],expected.endpoint_rate,2e-9);
                close(end.motion.base.twist_body.linear.z,-LEAD*expected.endpoint_rate,2e-9);
                close(end.motion.base.twist_body.angular.z,-IL/(IB+IL)*expected.endpoint_rate,2e-9);
                close(end.motion.joint_position[0],
                    state.mechanical.motion.joint_position[0]+dt*expected.rate,2e-9);
                close(receipt.midpoint_current_a,expected.current,2e-9);
                close(receipt.current_after_a,2.0*expected.current-state.currents_a[0],2e-9);
                close(receipt.midpoint_rate_rad_s,expected.rate,2e-9);
                close(out.mechanical.reaction_solve.impulses[0][0],dt*expected.normal,2e-9);
                close(out.mechanical.reaction_solve.impulses[0][1],-dt*0.8,2e-9);
                close(out.mechanical.midpoint_relative_velocities[0][0],0.0,2e-9);
                close(out.mechanical.endpoint_relative_velocities[0][0],0.0,2e-9);
                assert!(out.mechanical.path_admission >= 0.0);
                let momentum_after = complete_invariants(&model,end,&plates,rho);
                close(momentum_after.1.z-momentum_before.1.z,
                    dt*(expected.normal-(MB+ML)*G),2e-9);
                close((momentum_after.2-momentum_before.2).norm(),0.0,2e-9);
                let passive_heat = dt*passive.damping*receipt.midpoint_rate_rad_s.powi(2);
                let external_work = dt*0.8*out.mechanical.midpoint_relative_velocities[0][1];
                let balance = storage(&model,end,&plates,rho,passive)-before
                    +passive_heat+receipt.magnetic_energy_change_j+receipt.winding_heat_j
                    +receipt.bridge_heat_j-receipt.electrical_draw_j+receipt.electrical_return_j
                    -out.mechanical.reaction_work_j-external_work;
                close(balance,0.0,3e-9);
                if rate > 0.0 { assert!(receipt.electrical_return_j > 0.0); }
                if tick == 0 && rate == 0.0 { powered_first_rate.push(end.motion.joint_velocity[0]); }
                assert_eq!(state,original);
                state = out.into_state();
            }
        }
    }
    assert!(powered_first_rate[1] < powered_first_rate[0]-0.01);
}

#[test]
fn moving_support_friction_work_stays_separate_from_winding_and_source_work() {
    let model = mechanism();
    let rho = 1.2;
    let plates = geometry(rho);
    let strips: Vec<_> = plates.iter().map(|(i,p)|LinkedStrip::new(*i,p.clone()).unwrap()).collect();
    let state = initial(0.2,0.0);
    let motor = drive(0,1.0,0.2,0.4);
    let (dt,mu,conveyor) = (0.02,0.3,0.6);
    let expected = oracle(&state,motor,ZERO_IMPEDANCE,rho,dt);
    // A declared external load balances the conveyor's tangential force at the
    // SAME child origin. This isolates support/friction power without imposing
    // an extra kinematic restraint or altering the coupled axial motor solve.
    let external = -mu*expected.normal;
    let out = step_with_contacts(&model,&[ZERO_IMPEDANCE],&state,&[],&[motor],
        0.0,dt,&[contact(0.0,conveyor,mu)],contact_config(),
        |_,_,map|supported_loads(map,&strips,rho,external),admit_path,||false).unwrap();
    let tangential = mu*expected.normal*dt;
    close(out.mechanical.reaction_solve.impulses[0][0],expected.normal*dt,2e-9);
    close(out.mechanical.reaction_solve.impulses[0][1],tangential,2e-9);
    close(out.mechanical.friction_work_j,-tangential*conveyor,2e-9);
    close(out.mechanical.support_work_j,tangential*conveyor,2e-9);
    close(out.mechanical.reaction_work_j,0.0,2e-9);
    close(out.drives[0].midpoint_current_a,expected.current,2e-9);
    close(out.drives[0].electrical_balance_residual_j,0.0,2e-9);
    assert!(out.mechanical.friction_work_j < 0.0 && out.mechanical.support_work_j > 0.0);
    let end = &out.mechanical.step.state;
    let external_work = dt*external*(out.mechanical.midpoint_relative_velocities[0][1]+conveyor);
    close(storage(&model,end,&plates,rho,ZERO_IMPEDANCE)
        -storage(&model,&state.mechanical,&plates,rho,ZERO_IMPEDANCE)
        -out.drives[0].joint_work_j-out.mechanical.reaction_work_j-external_work,0.0,3e-9);
}

#[test]
fn no_contact_reduces_to_the_original_voltage_owner_without_extra_load_evaluation() {
    let model = single_rotor();
    let state = single_state(0.2,0.6);
    let motor = drive(0,0.8,0.3,0.4);
    let smooth = step_with_voltage_drives(&model,&[ZERO_IMPEDANCE],&state,&[],&[motor],
        0.0,0.02,one_config(),|_,_,map|loads(map,&[],1.0),||false).unwrap();
    let calls = Cell::new(0);
    let constrained = step_with_contacts(&model,&[ZERO_IMPEDANCE],&state,&[],&[motor],
        0.0,0.02,&[],contact_config(),
        |_,_,map| { calls.set(calls.get()+1); loads(map,&[],1.0) }, |_|Ok(()),||false).unwrap();
    assert_eq!(constrained.mechanical.step.state,smooth.mechanical.state);
    assert_eq!(constrained.drives,smooth.drives);
    assert_eq!(calls.get(),constrained.mechanical.step.iterations);
}

#[test]
fn capability_geometry_incoming_impact_and_late_cancel_preserve_both_states() {
    let model = mechanism();
    let state = initial(0.0,0.0);
    let original = state.clone();
    let motor = drive(0,1.0,0.2,0.4);
    let c = contact(0.0,0.0,0.3);
    for source in [false,true] {
        let mut bounded = motor;
        if source { bounded.control.source.maximum_source_input_w = 0.01; }
        else { bounded.motor.drive.maximum_current_a = 0.01; }
        let calls = Cell::new(0);
        let geometry = Cell::new(false);
        let error = step_with_contacts(&model,&[ZERO_IMPEDANCE],&state,&[],&[bounded],
            0.0,0.02,&[c],contact_config(),
            |_,_,map| { calls.set(calls.get()+1); supported_loads(map,&[],0.0,0.0) },
            |path| { geometry.set(true); admit_path(path) },||false).unwrap_err();
        assert_eq!(error.code,if source {"articulated-drive-power-limit"} else {"articulated-drive-current-limit"});
        assert!(calls.get() > 0 && !geometry.get());
        assert_eq!(state,original);
    }
    let error = step_with_contacts(&model,&[ZERO_IMPEDANCE],&state,&[],&[motor],
        0.0,0.02,&[c],contact_config(),|_,_,map|supported_loads(map,&[],0.0,0.0),
        |_|Err::<(),_>(refused("voltage-contact-geometry")),||false).unwrap_err();
    assert_eq!(error.code,"voltage-contact-geometry");
    let cancel = Cell::new(false);
    let error = step_with_contacts(&model,&[ZERO_IMPEDANCE],&state,&[],&[motor],
        0.0,0.02,&[c],contact_config(),|_,_,map|supported_loads(map,&[],0.0,0.0),
        |path| { let lower = admit_path(path)?; cancel.set(true); Ok(lower) },||cancel.get()).unwrap_err();
    assert_eq!(error.code,"articulated-step-cancelled");
    assert_eq!(state,original);
    let mut incoming = state.clone();
    incoming.mechanical.motion.base.twist_body.linear.z = -0.1;
    let saved = incoming.clone();
    let error = step_with_contacts(&model,&[ZERO_IMPEDANCE],&incoming,&[],&[motor],
        0.0,0.02,&[c],contact_config(),|_,_,map|supported_loads(map,&[],0.0,0.0),
        admit_path,||false).unwrap_err();
    assert_eq!(error.code,"articulated-contact-impact-required");
    assert_eq!(incoming,saved);
}

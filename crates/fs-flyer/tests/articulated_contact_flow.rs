//! G1: contact-coupled carriers in an accelerating uniform current.
//!
//! Independently close dry kinetic/gravity plus relative fluid storage against
//! ambient source work. Internal contact must transfer momentum without work.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_flyer::Refusal;
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{
    LinkedStrip, UniformFlow, assemble_articulated_strips_in_uniform_flow,
};
use fs_flyer::partitioned::articulated::contact::ArticulatedFrictionalContact;
use fs_flyer::partitioned::articulated::contact::systems::{
    ContactSystem, ContactSystemsPath, step_with_contacts,
};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{ContactCounterpart, LinkContactPoint};
use fs_mbd::articulated::generalized::free_floating_terms;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

const DRY: [f64; 2] = [3.0, 1.0];
const ADDED: [f64; 2] = [0.4, 2.3];
const RHO: f64 = 1.2;
const G: f64 = 9.81;
const FLOW_ACCEL: f64 = 2.4;
const FLOW_INITIAL: f64 = -0.4;
const INITIAL_SPEED: f64 = 0.3;
const GAP: f64 = 1e-7;

fn z(value: f64) -> Vec3 {
    Vec3::new(0.0, 0.0, value)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn body(mass: f64) -> ArticulatedModel {
    let inertia = SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}
fn plate(mass: f64) -> LinkedStrip {
    LinkedStrip::new(
        0,
        Strip {
            name: "normal fluid inertia",
            chord_m: 1.0,
            span_m: 4.0 * mass / (RHO * std::f64::consts::PI),
            position_m: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            span_axis: [0.0, 1.0, 0.0],
            control_coord: None,
            control_gain: 0.0,
            control_angular_gain: [0.0; 3],
        },
    )
    .unwrap()
}
fn initial(height: f64) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), Vec3::new(1.0, -0.5, height)).unwrap(),
                Twist::new(Vec3::default(), z(INITIAL_SPEED)),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        },
        memory: vec![],
    }
}
fn flow(time: f64) -> UniformFlow {
    UniformFlow {
        velocity_world_mps: z(FLOW_INITIAL + FLOW_ACCEL * time),
        acceleration_world_mps2: z(FLOW_ACCEL),
    }
}
fn energy(states: &[ArticulatedState], time: f64) -> f64 {
    states
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let v = s.motion.base.twist_body.linear.z;
            let relative = v - flow(time).velocity_world_mps.z;
            0.5 * DRY[i] * v * v
                + 0.5 * ADDED[i] * relative * relative
                + DRY[i] * G * s.motion.base.world_from_base.translation().z
        })
        .sum()
}

// Every body has exactly zero angular twist and no joints. Its declared material
// point follows an affine path, so the complete relative gap attains its minimum
// at an endpoint. This is an analytic admission for this fixture's complete
// feature set, rather than a general CCD claim from two samples.
fn admit_affine(path: &ContactSystemsPath<'_>) -> Result<f64, Refusal> {
    assert!(
        path.midpoints
            .iter()
            .all(|s| s.motion.base.twist_body.angular == Vec3::default()
                && s.motion.joint_position.is_empty())
    );
    let start = path.sample(0.0, || false)?;
    let end = path.sample(1.0, || false)?;
    let gap = |k: &[fs_mbd::articulated::Kinematics]| {
        k[0].world_from_link[0].translation().z - k[1].world_from_link[0].translation().z
    };
    let lower = gap(&start).min(gap(&end));
    assert!(lower >= 0.0, "negative complete affine gap {lower:e}");
    Ok(lower)
}

#[test]
fn accelerating_current_transfers_contact_momentum_and_closes_ambient_source_work() {
    let models = [body(DRY[0]), body(DRY[1])];
    let plates = [plate(ADDED[0]), plate(ADDED[1])];
    let mut states = vec![initial(2.0 + GAP), initial(2.0)];
    let contact = ArticulatedFrictionalContact {
        id: 71,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::default(),
        }),
        normal_world: z(1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 2.0 * GAP,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::Target { speed: 0.0 },
    };
    let smooth = ArticulatedStepConfig {
        max_iterations: 64,
        linear_tolerance_mps: 1e-12,
        angular_tolerance_rad_s: 1e-12,
        joint_velocity_tolerances: &[],
    };
    let dt = 0.01;
    let dry_mass = DRY.iter().sum::<f64>();
    let fluid_mass = ADDED.iter().sum::<f64>();
    let total_mass = dry_mass + fluid_mass;
    let acceleration = (fluid_mass * FLOW_ACCEL - dry_mass * G) / total_mass;
    let impulse = dt * ((DRY[0] + ADDED[0]) * acceleration - (ADDED[0] * FLOW_ACCEL - DRY[0] * G));
    assert!(impulse > 0.0);
    let initial_energy = energy(&states, 0.0);
    let initial_momentum = total_mass * INITIAL_SPEED - fluid_mass * FLOW_INITIAL;
    let mut ambient_work = 0.0;

    for tick in 0..16 {
        let time = tick as f64 * dt;
        let before = states.clone();
        let systems: Vec<_> = (0..2)
            .map(|i| ContactSystem {
                model: &models[i],
                state: &states[i],
                impedance: &[],
                memory_inputs: &[],
                smooth,
            })
            .collect();
        let accepted = step_with_contacts(
            &systems,
            time,
            dt,
            &[contact],
            FrictionalContactConfig::default(),
            |stage_time, stages| {
                stages
                    .iter()
                    .enumerate()
                    .map(|(i, stage)| {
                        Ok(ArticulatedLoads {
                            fluid: assemble_articulated_strips_in_uniform_flow(
                                RHO,
                                std::slice::from_ref(&plates[i]),
                                stage.kinematics,
                                flow(stage_time),
                            )?
                            .into_loads(),
                            world_wrenches: vec![Wrench::default()],
                            gravity_world: z(-G),
                            joint_effort: vec![],
                        })
                    })
                    .collect()
            },
            admit_affine,
            || false,
        )
        .unwrap();
        assert_eq!(states, before, "borrowed states changed during the solve");
        close(accepted.reaction_solve.impulses[0][0], impulse);
        close(accepted.reaction_work_j, 0.0);
        close(accepted.external_linear_impulse_world.norm(), 0.0);
        close(accepted.external_angular_impulse_world.norm(), 0.0);
        states = accepted.steps.into_iter().map(|s| s.state).collect();

        let endpoint_time = time + dt;
        let midpoint_time = time + 0.5 * dt;
        let mut momentum = 0.0;
        for i in 0..2 {
            let old = &before[i].motion.base;
            let new = &states[i].motion.base;
            close(
                new.twist_body.linear.z,
                INITIAL_SPEED + endpoint_time * acceleration,
            );
            close(new.twist_body.angular.norm(), 0.0);
            let midpoint = FreeFloatingBaseState::new(
                Se3::from_parts(
                    So3::identity(),
                    (old.world_from_base.translation() + new.world_from_base.translation())
                        .scale(0.5),
                )
                .unwrap(),
                Twist::new(
                    Vec3::default(),
                    (old.twist_body.linear + new.twist_body.linear).scale(0.5),
                ),
            );
            let map = free_floating_terms(&models[i], midpoint, &[], &[], || false)
                .unwrap()
                .kinematics;
            let fluid = assemble_articulated_strips_in_uniform_flow(
                RHO,
                std::slice::from_ref(&plates[i]),
                &map,
                flow(midpoint_time),
            )
            .unwrap();
            let a = (new.twist_body.linear.z - old.twist_body.linear.z) / dt;
            let reaction = fluid.reaction(&[0.0, 0.0, a, 0.0, 0.0, 0.0]).unwrap();
            let expected_force = ADDED[i] * (FLOW_ACCEL - acceleration);
            close(reaction.force_world.z, expected_force);
            close(
                reaction.body_power_w,
                expected_force * midpoint.twist_body.linear.z,
            );
            close(
                reaction.ambient_power_w,
                expected_force * flow(midpoint_time).velocity_world_mps.z,
            );
            close(reaction.balance_residual_w, 0.0);
            ambient_work += dt * reaction.ambient_power_w;

            let map = free_floating_terms(&models[i], *new, &[], &[], || false)
                .unwrap()
                .kinematics;
            let fluid = assemble_articulated_strips_in_uniform_flow(
                RHO,
                std::slice::from_ref(&plates[i]),
                &map,
                flow(endpoint_time),
            )
            .unwrap();
            let relative = new.twist_body.linear.z - flow(endpoint_time).velocity_world_mps.z;
            close(
                fluid.disturbance_energy_j(),
                0.5 * ADDED[i] * relative * relative,
            );
            momentum += DRY[i] * new.twist_body.linear.z + fluid.generalized_fluid_momentum()[2];
        }
        close(momentum, initial_momentum - dry_mass * G * endpoint_time);
        close(
            energy(&states, endpoint_time) - initial_energy,
            ambient_work,
        );
    }
}

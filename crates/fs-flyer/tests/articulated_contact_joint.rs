//! Sustained contact of a moving prismatic linkage with complete fluid inertia.
//! All geometry is linear: the child's fixed world height gives an independent
//! scalar carrier equation, including the implicit joint spring and damper.

use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::articulated::contact::{ContactStepConfig, step_with_contacts};
use fs_flyer::partitioned::articulated::impact::friction::ArticulatedFrictionalContact;
use fs_flyer::partitioned::articulated::impact::{ContactCounterpart, LinkContactPoint};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

fn v(z: f64) -> Vec3 {
    Vec3::new(0., 0., z)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10 * (1. + b.abs()), "{a:e} != {b:e}");
}
fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [1., 0., 0., 0., 1., 0., 0., 0., 1.],
        },
    )
    .unwrap()
}

#[test]
fn supported_slider_uses_complete_implicit_fluid_recoil_and_closes_physical_work() {
    let (mb, ml, q0, v0, g, k, c, tau, dt) = (2., 1., 0.05, -0.1, 9.81, 100., 3., 1.3, 0.02);
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(mb),
        ),
        Link::new(
            "supported slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(v(1.), None).unwrap(),
            inertia(ml),
        ),
    ])
    .unwrap();
    let initial = ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), v(-q0)).unwrap(),
                Twist::new(Vec3::default(), v(v0)),
            ),
            joint_position: vec![q0],
            joint_velocity: vec![-v0],
        },
        memory: Vec::new(),
    };
    let retained = initial.clone();
    let contact = ArticulatedFrictionalContact {
        id: 17,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 1,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        normal_world: v(1.),
        tangent_world: Vec3::new(1., 0., 0.),
        activation_tolerance: 1e-9,
        friction_coefficient: 0.,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0. },
    };
    let impedance = [JointImpedance {
        stiffness: k,
        damping: c,
        rest_position: 0.,
    }];
    let config = ContactStepConfig {
        smooth: ArticulatedStepConfig {
            max_iterations: 32,
            linear_tolerance_mps: 1e-13,
            angular_tolerance_rad_s: 1e-13,
            joint_velocity_tolerances: &[1e-13],
        },
        contact: FrictionalContactConfig::default(),
    };
    let strip = LinkedStrip::new(
        0,
        Strip {
            name: "carrier attached fluid",
            chord_m: 1.,
            span_m: 2.,
            position_m: [0.; 3],
            normal: [0., 0., 1.],
            span_axis: [0., 1., 0.],
            control_coord: None,
            control_gain: 0.,
            control_angular_gain: [0.; 3],
        },
    )
    .unwrap();
    let mut endpoint_speeds = Vec::new();
    for rho in [0., 1.2] {
        let result = step_with_contacts(
            &model,
            &impedance,
            &initial,
            &[],
            0.2,
            dt,
            &[contact],
            config,
            |_, _, map| {
                Ok(ArticulatedLoads {
                    fluid: assemble_articulated_strips(
                        if rho == 0. { 1. } else { rho },
                        if rho == 0. {
                            &[]
                        } else {
                            core::slice::from_ref(&strip)
                        },
                        map,
                    )?,
                    world_wrenches: vec![Wrench::new(Vec3::default(), Vec3::default()); 2],
                    gravity_world: v(-g),
                    joint_effort: vec![tau],
                })
            },
            |path| {
                // For this declared nonrotating z-slider the exact path height
                // is affine. Its two endpoint heights establish the whole-path
                // lower bound; a midpoint tangent alone would not do so for a
                // curved or rotating contact path.
                assert_eq!(path.initial.motion.base.twist_body.angular, Vec3::default());
                assert_eq!(
                    path.midpoint.motion.base.twist_body.angular,
                    Vec3::default()
                );
                assert_eq!(
                    path.endpoint.motion.base.twist_body.angular,
                    Vec3::default()
                );
                let height = |s: &ArticulatedState| {
                    s.motion.base.world_from_base.translation().z + s.motion.joint_position[0]
                };
                let bounds = (height(path.initial), height(path.endpoint));
                assert!(bounds.0 >= 0. && bounds.1 >= 0.);
                Ok(bounds)
            },
            || false,
        )
        .unwrap();
        let end = &result.step.state;
        let added_mass = rho * std::f64::consts::PI / 2.;
        let b = mb + added_mass;
        let h = 0.5 * dt;
        let a = c + h * k;
        let expected_v = ((b - h * a) * v0 + dt * (-mb * g - tau + k * q0)) / (b + h * a);
        let vm = 0.5 * (v0 + expected_v);
        let expected_q = q0 - dt * vm;
        let expected_p = dt * (ml * g - tau - c * vm + k * (q0 - h * vm));
        assert!(expected_p > 0.);
        close(end.motion.base.twist_body.linear.z, expected_v);
        close(end.motion.joint_velocity[0], -expected_v);
        close(end.motion.joint_position[0], expected_q);
        close(end.motion.base.world_from_base.translation().z, -expected_q);
        close(result.reaction_solve.impulses[0][0], expected_p);
        close(result.reaction_solve.impulses[0][1], 0.);
        close(result.reaction_solve.impulses[0][2], 0.);
        assert_eq!(result.contact_ids, vec![17]);
        close(result.midpoint_relative_velocities[0][0], 0.);
        close(result.endpoint_relative_velocities[0][0], 0.);
        close(result.path_admission.0, 0.);
        close(result.path_admission.1, 0.);
        close(result.reaction_work_j, 0.);
        close(result.normal_work_j, 0.);
        close(result.friction_work_j, 0.);
        close(result.support_work_j, 0.);
        // Child velocity and world height remain zero. Gravity acts on dry
        // mass only; the attached fluid contributes kinetic storage. Midpoint
        // quadratics give the exact spring, gravity and damping work identity.
        let kinetic_change = 0.5 * b * (expected_v * expected_v - v0 * v0);
        let spring_change = 0.5 * k * (expected_q * expected_q - q0 * q0);
        let gravity_change = -mb * g * (expected_q - q0);
        let damping_loss = dt * c * vm * vm;
        let actuator_work = tau * (expected_q - q0);
        close(
            kinetic_change + spring_change + gravity_change + damping_loss - actuator_work,
            0.,
        );
        let actual_kinetic_change =
            0.5 * b * (end.motion.base.twist_body.linear.z.powi(2) - v0 * v0)
                + 0.5
                    * ml
                    * (end.motion.base.twist_body.linear.z + end.motion.joint_velocity[0]).powi(2);
        let actual_q = end.motion.joint_position[0];
        let actual_z = end.motion.base.world_from_base.translation().z;
        let actual_gravity_change = mb * g * (actual_z + q0) + ml * g * (actual_z + actual_q);
        let actual_damping_loss = dt
            * c
            * (0.5 * (initial.motion.joint_velocity[0] + end.motion.joint_velocity[0])).powi(2);
        close(
            actual_kinetic_change
                + 0.5 * k * (actual_q * actual_q - q0 * q0)
                + actual_gravity_change
                + actual_damping_loss
                - tau * (actual_q - q0),
            0.,
        );
        endpoint_speeds.push(end.motion.base.twist_body.linear.z);
        assert_eq!(initial, retained);
    }
    assert!(
        endpoint_speeds[1] > endpoint_speeds[0] + 0.05,
        "attached fluid must reduce the carrier's downward acceleration"
    );
}

//! Actual native motion under opposing loads with a resolvable net force/couple.
//! G0/G3: independent binary-arithmetic trajectories, cross-owner assembly,
//! and material-point moments; no high-precision substitute for the solver.

use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{PrescribedBaseMotionState, PrescribedBasePath};
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MaterialFrame, SpatialBushing, evaluate_bushing_forces,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::{
    step_prescribed_base_bushings, step_prescribed_base_spatial,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad, step_free_floating_spatial,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

const LARGE: f64 = 18_014_398_509_481_984.0;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        },
    )
    .unwrap()
}

fn single(joint: JointModel) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        inertia(2.0),
    )])
    .unwrap()
}

fn point(link: usize, point_body: Vec3) -> MaterialPoint {
    MaterialPoint { link, point_body }
}

fn force(point: MaterialPoint, x: f64) -> WorldPointLoad {
    WorldPointLoad {
        point,
        force_world: Vec3::new(x, 0.0, 0.0),
        couple_world: Vec3::default(),
    }
}

fn supported(speed: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(Se3::identity()),
        joint_position: vec![0.0],
        joint_velocity: vec![speed],
    }
}

#[test]
fn held_and_world_forces_keep_the_net_supported_acceleration_and_work() {
    let model = single(JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap());
    let state = supported(1.0);
    let held = [Wrench::new(Vec3::default(), Vec3::new(LARGE, 0.0, 0.0))];
    let p = point(0, Vec3::default());
    for world_points in [
        [force(p, 1.0), force(p, -LARGE)],
        [force(p, -LARGE), force(p, 1.0)],
    ] {
        let result = step_prescribed_base_spatial(
            &model,
            &state,
            PrescribedBasePath::stationary(Se3::identity()),
            SpatialMotionLoads {
                held: FreeFloatingMotionLoads {
                    generalized_force: &[0.0],
                    external_body_wrench: &held,
                    gravity_world: Vec3::default(),
                },
                world_points: &world_points,
                connections: &[],
            },
            0.125,
            || false,
        )
        .unwrap();
        // m=2, F=1, h=1/8: q=h+F h^2/(2m), v=1+F h/m.
        close(
            result.motion.state_after.joint_position[0],
            0.12890625,
            1.0e-15,
        );
        close(result.motion.state_after.joint_velocity[0], 1.0625, 1.0e-15);
        close(result.motion.external_work, 0.12890625, 1.0e-14);
        close(result.motion.kinetic_work_residual, 0.0, 1.0e-14);
    }
}

#[test]
fn free_flight_keeps_small_linear_and_angular_impulses_across_load_owners() {
    let model = single(JointModel::FIXED);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let p = point(0, Vec3::new(0.0, 1.0, 0.0));
    let held = [Wrench::new(
        Vec3::new(0.0, 0.0, -LARGE),
        Vec3::new(LARGE, 0.0, 0.0),
    )];
    let world_points = [force(p, 1.0), force(p, -LARGE)];
    let result = step_free_floating_spatial(
        &model,
        &state,
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &held,
                gravity_world: Vec3::default(),
            },
            world_points: &world_points,
            connections: &[],
        },
        0.125,
        || false,
    )
    .unwrap();
    close(
        result.motion.momentum_after.linear_momentum_world.x,
        2.125,
        1.0e-14,
    );
    close(
        result.motion.momentum_after.angular_momentum_world.z,
        -0.125,
        1.0e-14,
    );
    close(
        result.motion.external_linear_impulse_world.x,
        0.125,
        1.0e-15,
    );
    close(
        result.motion.external_angular_impulse_world.z,
        -0.125,
        1.0e-15,
    );
    close(
        result.motion.center_of_mass_world_after.x,
        0.12890625,
        1.0e-15,
    );
    close(result.motion.external_work, 0.1328125, 1.0e-14);
    close(result.motion.kinetic_work_residual, 0.0, 1.0e-13);
}

#[test]
fn off_center_force_retains_a_small_moment_between_large_products() {
    let model = single(JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap());
    let state = supported(0.0);
    let n = 134_217_728.0;
    // (n+1)(n-1) - n*n = -1 exactly for these stored binary inputs.
    let applied = WorldPointLoad {
        point: point(0, Vec3::new(n + 1.0, n, 0.0)),
        force_world: Vec3::new(n, n - 1.0, 0.0),
        couple_world: Vec3::default(),
    };
    let result = step_prescribed_base_spatial(
        &model,
        &state,
        PrescribedBasePath::stationary(Se3::identity()),
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.0],
                external_body_wrench: &[Wrench::default()],
                gravity_world: Vec3::default(),
            },
            world_points: &[applied],
            connections: &[],
        },
        0.125,
        || false,
    )
    .unwrap();
    close(
        result.motion.state_after.joint_position[0],
        -0.00390625,
        1.0e-15,
    );
    close(
        result.motion.state_after.joint_velocity[0],
        -0.0625,
        1.0e-15,
    );
    close(result.motion.external_work, 0.00390625, 1.0e-14);
}

fn slider() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0),
        ),
    ])
    .unwrap()
}

fn bushing(stiffness: f64, rest: f64) -> SpatialBushing {
    SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 1,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::new(rest, 0.0, 0.0),
        translation_stiffness: Vec3::new(stiffness, 0.0, 0.0),
        translation_damping: Vec3::default(),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
    }
}

#[test]
fn held_world_and_attachment_forces_are_rounded_only_after_complete_assembly() {
    let model = slider();
    let state = supported(0.0);
    let held = [
        Wrench::default(),
        Wrench::new(Vec3::default(), Vec3::new(LARGE, 0.0, 0.0)),
    ];
    let world = [force(point(1, Vec3::default()), 1.0)];
    let connection = AxialConnection {
        a: point(0, Vec3::default()),
        b: point(1, Vec3::default()),
        rest_length_m: 0.0,
        stiffness_n_per_m: LARGE,
        damping_ns_per_m: 0.0,
    };
    let h = 2.0_f64.powi(-30);
    for use_bushing in [false, true] {
        let connections = if use_bushing {
            &[][..]
        } else {
            std::slice::from_ref(&connection)
        };
        let elements = [bushing(LARGE, 0.0)];
        let bushings = if use_bushing { &elements[..] } else { &[][..] };
        let result = step_prescribed_base_bushings(
            &model,
            &state,
            PrescribedBasePath::stationary(Se3::identity()),
            SpatialMotionLoads {
                held: FreeFloatingMotionLoads {
                    generalized_force: &[0.0],
                    external_body_wrench: &held,
                    gravity_world: Vec3::default(),
                },
                world_points: &world,
                connections,
            },
            bushings,
            h,
            || false,
        )
        .unwrap();
        // The first midpoint has the original extension. Its assembled force
        // is LARGE + 1 - LARGE = 1, independently of attachment representation.
        close(
            result.motion.state_after.joint_velocity[0],
            0.5 * h,
            1.0e-24,
        );
        close(
            result.motion.state_after.joint_position[0],
            0.25 * h * h,
            1.0e-33,
        );
    }
}

#[test]
fn opposed_bushings_retain_the_resolvable_reciprocal_force() {
    let model = slider();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: vec![0.0],
        joint_velocity: vec![0.0],
    };
    let elements = [bushing(LARGE, 0.0), bushing(1.0, 0.0), bushing(LARGE, 2.0)];
    let result = evaluate_bushing_forces(&model, &state, &elements, || false).unwrap();
    close(result.body_wrench[0].force.x, 1.0, 0.0);
    close(result.body_wrench[1].force.x, -1.0, 0.0);
    close(result.body_wrench[0].torque.norm(), 0.0, 0.0);
    close(result.body_wrench[1].torque.norm(), 0.0, 0.0);
}

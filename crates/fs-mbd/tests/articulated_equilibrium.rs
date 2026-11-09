//! G1 static force/torque balance and G0 continuation/refusal on native owners.

#[path = "articulated_equilibrium/rotating.rs"]
mod rotating;

#[path = "articulated_equilibrium/linearized.rs"]
mod linearized;

use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::prescribed::equilibrium::{
    EquilibriumBushings, EquilibriumConfig, EquilibriumError, equilibrate, equilibrate_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::step_prescribed_base_spatial;
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}

fn inertia(mass: f64, center: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
    )
    .unwrap()
}

fn seed(position: Vec<f64>) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 17.0,
        base: BaseState::stationary(Se3::identity()),
        joint_velocity: vec![0.0; position.len()],
        joint_position: position,
    }
}

fn controls<'a>(coordinates: &'a [f64], efforts: &'a [f64]) -> EquilibriumConfig<'a> {
    EquilibriumConfig {
        coordinate_scales: coordinates,
        effort_scales: efforts,
        residual_tolerance: 1e-11,
        difference_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 1024,
        krylov_restart: 16,
    }
}

fn slider() -> (ArticulatedModel, AxialConnection) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default()),
        ),
    ])
    .unwrap();
    let connection = AxialConnection {
        a: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        rest_length_m: 1.0,
        stiffness_n_per_m: 100.0,
        damping_ns_per_m: 7.0,
    };
    (model, connection)
}

#[test]
fn loaded_axial_equilibrium_balances_weight_actuation_and_world_force_then_stays_at_rest() {
    let (model, connection) = slider();
    let initial = seed(vec![0.1]);
    let snapshot = initial.clone();
    let point = WorldPointLoad {
        point: connection.b,
        force_world: Vec3::new(3.0, 0.0, 0.0),
        couple_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[1.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(9.81, 0.0, 0.0),
        },
        world_points: &[point],
        connections: &[connection],
    };
    let solved = equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[0.1], &[20.0]),
        || false,
    )
    .unwrap();
    let exact = (2.0 * 9.81 + 3.0 + 1.0) / connection.stiffness_n_per_m;
    close(solved.state.joint_position[0], exact, 1e-12);
    close(solved.required_effort[0], 1.0, 2e-10);
    close(solved.joint_effort_residual[0], 0.0, 2e-10);
    close(
        solved.support_wrench_world.force.x,
        -(5.0 * 9.81 + 3.0),
        1e-12,
    );
    close(
        solved.attachment_energy_j,
        0.5 * connection.stiffness_n_per_m * exact * exact,
        1e-11,
    );
    assert_eq!(solved.state.time_s, initial.time_s);
    assert_eq!(initial, snapshot);
    assert_eq!(solved.state.joint_velocity, vec![0.0]);
    let mut state = solved.state;
    for _ in 0..20 {
        state = step_prescribed_base_spatial(
            &model,
            &state,
            PrescribedBasePath::stationary(Se3::identity()),
            loads,
            0.002,
            || false,
        )
        .unwrap()
        .motion
        .state_after;
    }
    close(state.joint_position[0], exact, 2e-12);
    close(state.joint_velocity[0], 0.0, 2e-12);
}

#[test]
fn finite_angle_gravity_and_the_original_torsional_bushing_find_the_same_balance() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "pendulum",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.7, 0.0, 0.0)),
        ),
    ])
    .unwrap();
    let bushings = [SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 1,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::default(),
        translation_damping: Vec3::default(),
        rotation_stiffness: Vec3::new(0.0, 0.0, 20.0),
        rotation_damping: Vec3::new(0.0, 0.0, 3.0),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(0.0, -9.81, 0.0),
        },
        world_points: &[],
        connections: &[],
    };
    let solved = equilibrate(
        &model,
        &seed(vec![-0.1]),
        loads,
        EquilibriumBushings::Diagonal(&bushings),
        controls(&[0.5], &[20.0]),
        || false,
    )
    .unwrap();
    let q = solved.state.joint_position[0];
    assert!(q < -0.4 && q > -0.7);
    // Independent moment balance about the hinge: k*q + m*g*l*cos(q) = 0.
    close(20.0 * q + 2.0 * 9.81 * 0.7 * q.cos(), 0.0, 3e-10);
    close(solved.support_wrench_world.force.y, 5.0 * 9.81, 1e-12);
    close(
        solved.support_wrench_world.torque.z,
        2.0 * 9.81 * 0.7 * q.cos(),
        2e-12,
    );
}

#[test]
fn coupled_translation_rotation_bushing_preserves_cross_terms_and_physical_scaling() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(1.0, Vec3::default()),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, Vec3::default()),
        ),
        Link::new(
            "rotor",
            Some(1),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0, Vec3::default()),
        ),
    ])
    .unwrap();
    let mut matrix = [[0.0; 6]; 6];
    matrix[0][0] = 100.0;
    matrix[5][5] = 20.0;
    matrix[0][5] = 4.0;
    matrix[5][0] = 4.0;
    let scales = BushingCoordinateScales {
        translation_m: 0.1,
        rotation_rad: 0.5,
    };
    let bushing = CoupledSpatialBushing::try_new(
        MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        Vec3::default(),
        PassiveBushingMatrix::try_new(matrix, scales).unwrap(),
        PassiveBushingMatrix::try_new([[0.0; 6]; 6], scales).unwrap(),
    )
    .unwrap();
    let load = WorldPointLoad {
        point: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        force_world: Vec3::new(6.0, 0.0, 0.0),
        couple_world: Vec3::new(0.0, 0.0, 2.0),
    };
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0; 2],
            external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::default(),
        },
        world_points: &[load],
        connections: &[],
    };
    let bushings = [bushing];
    let a = equilibrate(
        &model,
        &seed(vec![0.0; 2]),
        loads,
        EquilibriumBushings::Coupled(&bushings),
        controls(&[0.1, 0.5], &[10.0, 3.0]),
        || false,
    )
    .unwrap();
    let b = equilibrate(
        &model,
        &seed(vec![0.0; 2]),
        loads,
        EquilibriumBushings::Coupled(&bushings),
        controls(&[2.0, 0.03], &[0.7, 80.0]),
        || false,
    )
    .unwrap();
    let det = 100.0 * 20.0 - 4.0 * 4.0;
    for result in [a, b] {
        close(
            result.state.joint_position[0],
            (20.0 * 6.0 - 4.0 * 2.0) / det,
            3e-12,
        );
        close(
            result.state.joint_position[1],
            (100.0 * 2.0 - 4.0 * 6.0) / det,
            3e-12,
        );
        close(result.support_wrench_world.force.x, -6.0, 3e-12);
        close(result.support_wrench_world.torque.z, -2.0, 3e-12);
    }
}

#[test]
fn original_joint_effort_law_uses_stationary_trial_coordinates_and_unchanged_clock() {
    let (model, _) = slider();
    let initial = seed(vec![0.0]);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[4.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let result = equilibrate_with_effort(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[0.2], &[4.0]),
        |elapsed, state, effort| {
            assert_eq!(elapsed, 0.0);
            assert_eq!(state.time_s, initial.time_s);
            assert_eq!(state.joint_velocity, vec![0.0]);
            assert_eq!(effort, &[4.0]);
            effort[0] -= 20.0 * state.joint_position[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    close(result.state.joint_position[0], 0.2, 1e-12);
    close(result.applied_effort[0], 0.0, 4e-11);
    close(result.attachment_energy_j, 0.0, 0.0); // Callback storage stays with its owner.
}

#[test]
fn no_equilibrium_budget_and_cancellation_do_not_publish_a_candidate() {
    let (model, _) = slider();
    let initial = seed(vec![0.25]);
    let snapshot = initial.clone();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[1.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    assert!(matches!(
        equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls(&[1.0], &[1.0]),
            || false
        ),
        Err(EquilibriumError::Unresolved(_))
    ));
    let mut tiny = controls(&[1.0], &[1.0]);
    tiny.maximum_evaluations = 2;
    assert!(matches!(
        equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            tiny,
            || false
        ),
        Err(EquilibriumError::EvaluationBudget {
            used: 2,
            maximum: 2
        })
    ));
    assert!(matches!(
        equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls(&[1.0], &[1.0]),
            || true
        ),
        Err(EquilibriumError::Motion(
            PrescribedBaseMotionError::Cancelled
        ))
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn all_fixed_tree_recovers_support_and_rejects_moving_seeds_or_bad_scales() {
    let model = ArticulatedModel::new(vec![Link::new(
        "root",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, Vec3::new(0.5, 0.0, 0.0)),
    )])
    .unwrap();
    let initial = seed(vec![]);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::new(0.0, -10.0, 0.0),
        },
        world_points: &[],
        connections: &[],
    };
    let result = equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[], &[]),
        || false,
    )
    .unwrap();
    assert_eq!(result.state, initial);
    assert_eq!(result.residual_evaluations, 2);
    close(result.support_wrench_world.force.y, 20.0, 0.0);
    close(result.support_wrench_world.torque.z, 10.0, 0.0);
    let mut moving = initial.clone();
    moving.base.twist_body = Twist::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0));
    assert!(matches!(
        equilibrate(
            &model,
            &moving,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls(&[], &[]),
            || false
        ),
        Err(EquilibriumError::Invalid(_))
    ));
    let mut excessive = controls(&[], &[]);
    excessive.maximum_iterations = 65;
    assert!(matches!(
        equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            excessive,
            || false
        ),
        Err(EquilibriumError::Invalid(_))
    ));
    assert!(matches!(
        equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls(&[f64::NAN], &[]),
            || false
        ),
        Err(EquilibriumError::Invalid(_))
    ));
}

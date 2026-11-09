//! Independent lock/gear statics and continuation through native dynamics.
use super::*;
use fs_mbd::articulated::JointLimits;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::equilibrium::constrained::{
    self as owner, JointEquilibriumControls,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitSpatialConfig,
    constraints::{SupportedConstraintControls, step_prescribed_base_spatial as constrained_step},
};

fn row(coefficients: Vec<f64>, position: f64) -> JointConstraint {
    JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients,
        position,
        position_tolerance: 1e-10,
        velocity_tolerance: 1e-10,
    }
}
fn joint_controls<'a>(reaction: &'a [f64], position: &'a [f64]) -> JointEquilibriumControls<'a> {
    JointEquilibriumControls {
        reaction_scales: reaction,
        position_scales: position,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn two_sliders() -> (ArticulatedModel, [AxialConnection; 2]) {
    let mut links = vec![Link::new(
        "support",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(3.0, Vec3::default()),
    )];
    let connections = core::array::from_fn(|index| {
        links.push(Link::new(
            format!("slider-{index}"),
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(1.0, index as f64, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0 + index as f64, Vec3::default()),
        ));
        AxialConnection {
            a: MaterialPoint {
                link: 0,
                point_body: Vec3::new(0.0, index as f64, 0.0),
            },
            b: MaterialPoint {
                link: index + 1,
                point_body: Vec3::default(),
            },
            rest_length_m: 1.0,
            stiffness_n_per_m: 100.0 * (index + 1) as f64,
            damping_ns_per_m: 7.0,
        }
    });
    (ArticulatedModel::new(links).unwrap(), connections)
}

#[test]
fn g1_geared_loaded_root_matches_virtual_work_and_native_constrained_motion() {
    let (model, connections) = two_sliders();
    let initial = seed(vec![0.12, -0.04]); // Deliberately off the final manifold.
    let original = initial.clone();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0, 40.0],
            external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &connections,
    };
    let rows = [row(vec![-2.0, 1.0], 0.1)];
    let root = owner::equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        &rows,
        controls(&[0.1; 2], &[40.0; 2]),
        joint_controls(&[20.0], &[0.1]),
        || false,
    )
    .unwrap();
    assert_eq!(initial, original);
    // q1=2*q0+0.1, and virtual work requires (100*q0-10)+2*(200*q1-40)=0.
    close(root.state.joint_position[0], 1.0 / 18.0, 2e-12);
    close(root.state.joint_position[1], 19.0 / 90.0, 3e-12);
    close(root.multipliers[0], 20.0 / 9.0, 2e-10);
    close(root.reaction_effort[0], -40.0 / 9.0, 3e-10);
    close(root.reaction_effort[1], 20.0 / 9.0, 2e-10);
    close(
        root.reaction_effort[0] + 2.0 * root.reaction_effort[1],
        0.0,
        1e-14,
    );
    assert!(root.joint_effort_residual.iter().all(|r| r.abs() < 1e-9));
    assert_eq!(root.applied_effort, [10.0, 40.0]);
    assert!(root.position_residuals[0].abs() <= rows[0].position_tolerance);
    let path = PrescribedBasePath {
        midpoint: root.state.base,
        endpoint: root.state.base,
    };
    let step = constrained_step(
        &model,
        &root.state,
        path,
        loads,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        SupportedConstraintControls {
            minimum_reciprocal_condition: 1e-10,
            maximum_reaction_work_j: 1e-10,
        },
        || false,
    )
    .unwrap();
    for (&q, &expected) in step
        .implicit
        .spatial
        .motion
        .state_after
        .joint_position
        .iter()
        .zip(&root.state.joint_position)
    {
        close(q, expected, 2e-11);
    }
    for &v in &step.implicit.spatial.motion.state_after.joint_velocity {
        close(v, 0.0, 2e-10);
    }
    close(
        step.midpoint_reaction.multipliers[0],
        root.multipliers[0],
        1e-8,
    );
    close(step.implicit.spatial.motion.constraint_work, 0.0, 1e-12);
}

#[test]
fn g1_a_locked_loaded_joint_needs_no_stiffness_or_artificial_actuator_capacity() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "locked",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1.0, 0.0, 0.0),
                Some(JointLimits::new(-1.0, 1.0, 1.0, 0.1).unwrap()),
            )
            .unwrap(),
            inertia(2.0, Vec3::default()),
        ),
    ])
    .unwrap();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.05],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(9.81, 0.0, 0.0),
        },
        world_points: &[],
        connections: &[],
    };
    let root = owner::equilibrate(
        &model,
        &seed(vec![-0.2]),
        loads,
        EquilibriumBushings::Diagonal(&[]),
        &[row(vec![1.0], 0.3)],
        controls(&[0.1], &[20.0]),
        joint_controls(&[20.0], &[0.1]),
        || false,
    )
    .unwrap();
    close(root.state.joint_position[0], 0.3, 1e-12);
    close(root.required_effort[0], -19.62, 1e-12);
    close(root.applied_effort[0], 0.05, 1e-15);
    close(root.multipliers[0], -19.67, 3e-10);
    assert!(root.reaction_effort[0].abs() > 100.0 * 0.1);
    close(root.joint_effort_residual[0], 0.0, 3e-10);
}

#[test]
fn g3_rescaled_signed_rows_preserve_the_physical_root_and_reaction() {
    let (model, connections) = two_sliders();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0, 40.0],
            external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &connections,
    };
    for scale in [1e-6_f64, -1.0, 1e6] {
        let mut constraint = row(vec![-2.0 * scale, scale], 0.1 * scale);
        constraint.position_tolerance *= scale.abs();
        constraint.velocity_tolerance *= scale.abs();
        let root = owner::equilibrate(
            &model,
            &seed(vec![0.0; 2]),
            loads,
            EquilibriumBushings::Diagonal(&[]),
            &[constraint],
            controls(&[0.1; 2], &[40.0; 2]),
            joint_controls(&[20.0 / scale.abs()], &[0.1 * scale.abs()]),
            || false,
        )
        .unwrap();
        close(root.state.joint_position[0], 1.0 / 18.0, 3e-12);
        close(root.state.joint_position[1], 19.0 / 90.0, 3e-12);
        close(root.multipliers[0] * scale, 20.0 / 9.0, 3e-10);
        close(root.reaction_effort[0], -40.0 / 9.0, 5e-10);
    }
}

#[test]
fn g1_rotating_lock_and_repeatable_actuator_feedback_use_original_force_owner() {
    let (model, connection) = slider();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    let mut initial = seed(vec![0.0]);
    initial.base.twist_body = Twist::new(Vec3::new(0.0, 0.0, 3.0), Vec3::default());
    let rows = [row(vec![1.0], 0.2)];
    let root = owner::equilibrate_rotating(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        &rows,
        controls(&[0.1], &[100.0]),
        joint_controls(&[20.0], &[0.1]),
        || false,
    )
    .unwrap();
    // Outward centrifugal load uses the actual 1+q radius and native mass 2.
    close(root.multipliers[0], 100.0 * 0.2 - 2.0 * 9.0 * 1.2, 3e-10);
    assert_eq!(root.state.base.twist_body, initial.base.twist_body);
    initial.base.twist_body = Twist::zero();
    let root = owner::equilibrate_with_effort(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        &rows,
        controls(&[0.1], &[100.0]),
        joint_controls(&[20.0], &[0.1]),
        |elapsed, state, force| {
            assert_eq!(elapsed, 0.0);
            assert_eq!(state.time_s, 17.0);
            force[0] = 3.0 + 10.0 * state.joint_position[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    close(root.multipliers[0], 15.0, 3e-10);
    close(root.applied_effort[0], 5.0, 2e-11);
}

#[test]
fn g0_empty_rows_retain_original_result_and_invalid_rows_never_publish() {
    let (model, connection) = slider();
    let initial = seed(vec![0.1]);
    let snapshot = initial.clone();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    let config = controls(&[0.1], &[20.0]);
    let plain = equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        config,
        || false,
    )
    .unwrap();
    let empty = owner::equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        &[],
        config,
        joint_controls(&[], &[]),
        || false,
    )
    .unwrap();
    assert_eq!(empty.state, plain.state);
    assert_eq!(empty.required_effort, plain.required_effort);
    assert_eq!(empty.joint_effort_residual, plain.joint_effort_residual);
    assert_eq!(empty.nonlinear, plain.nonlinear);
    assert_eq!(empty.residual_evaluations, plain.residual_evaluations);
    let duplicate = [row(vec![1.0], 0.0), row(vec![2.0], 0.0)];
    assert!(matches!(
        owner::equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            &duplicate,
            config,
            joint_controls(&[1.0; 2], &[1.0; 2]),
            || false
        ),
        Err(EquilibriumError::Constraint(_))
    ));
    let rows = [row(vec![1.0], 0.3)];
    assert!(
        owner::equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            &rows,
            config,
            joint_controls(&[0.0], &[1.0]),
            || false
        )
        .is_err()
    );
    let mut short = config;
    short.maximum_evaluations = 2;
    assert!(matches!(
        owner::equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            &rows,
            short,
            joint_controls(&[20.0], &[0.1]),
            || false
        ),
        Err(EquilibriumError::EvaluationBudget { .. })
    ));
    assert!(
        owner::equilibrate(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            &rows,
            config,
            joint_controls(&[20.0], &[0.1]),
            || true
        )
        .is_err()
    );
    assert_eq!(initial, snapshot);
}

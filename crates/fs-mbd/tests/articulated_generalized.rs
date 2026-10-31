//! Physical and independent-recursion oracles for native coupled inertia.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::{
    GeneralizedError, MAX_GENERALIZED_JOINTS, MAX_GENERALIZED_LINKS, free_floating_terms,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, free_floating_forward_dynamics, free_floating_momentum, inverse_dynamics,
};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-11 * (1.0 + b.abs()), "{a} != {b}");
}
fn twist_close(a: Twist, b: Twist) {
    for (a, b) in a.to_array().into_iter().zip(b.to_array()) {
        close(a, b);
    }
}
fn inertia(mass: f64, center: Vec3, diagonal: [f64; 3]) -> SpatialInertia {
    let [x, y, z] = diagonal;
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
        },
    )
    .unwrap()
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn pose(v: Vec3, p: Vec3) -> Se3 {
    Se3::from_parts(rotation(v), p).unwrap()
}
fn product(matrix: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum())
        .collect()
}

#[test]
fn serial_rotors_retain_locked_base_and_full_joint_mass_without_unit_rate_probes() {
    let limits = JointLimits::new(-2.0, 2.0, 0.2, 1e-6).unwrap();
    let joint = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), Some(limits)).unwrap();
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, Vec3::default(), [1.0; 3]),
        ),
        Link::new(
            "first",
            Some(0),
            Se3::identity(),
            joint,
            inertia(2.0, Vec3::default(), [2.0; 3]),
        ),
        Link::new(
            "second",
            Some(1),
            Se3::identity(),
            joint,
            inertia(2.0, Vec3::default(), [3.0; 3]),
        ),
    ])
    .unwrap();
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 0.4), Vec3::default()),
    );
    let terms = free_floating_terms(&model, base, &[0.3, -0.2], &[0.1, -0.15], || false).unwrap();
    // Absolute spins are [omega, omega+qdot1, omega+qdot1+qdot2].
    let angular = [[6.0, 5.0, 3.0], [5.0, 5.0, 3.0], [3.0, 3.0, 3.0]];
    for (i, row) in [2, 6, 7].into_iter().enumerate() {
        for (j, column) in [2, 6, 7].into_iter().enumerate() {
            close(terms.mass[row][column], angular[i][j]);
        }
    }
    for axis in 3..6 {
        close(terms.mass[axis][axis], 6.0);
    }
    for value in terms.bias {
        close(value, 0.0);
    }
    assert_eq!(terms.kinematics.dof_count(), 2);
    assert_eq!(terms.kinematics.coordinate_count(), 8);
}

fn mixed_tree() -> ArticulatedModel {
    let tensor = SpatialInertia::new(
        1.7,
        Vec3::new(0.3, -0.1, 0.05),
        Mat3 {
            m: [0.2, 0.02, -0.01, 0.02, 0.3, 0.015, -0.01, 0.015, 0.4],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![
        Link::new(
            "root",
            None,
            pose(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.1, 0.2, -0.1)),
            JointModel::FIXED,
            inertia(3.0, Vec3::new(0.2, 0.1, -0.05), [0.7, 0.9, 1.0]),
        ),
        Link::new(
            "revolute",
            Some(0),
            pose(Vec3::new(0.2, 0.3, -0.1), Vec3::new(0.7, 0.1, 0.2)),
            JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            tensor,
        ),
        Link::new(
            "screw",
            Some(1),
            pose(Vec3::new(-0.2, 0.1, 0.4), Vec3::new(0.4, -0.2, 0.3)),
            JointModel::helical(Vec3::new(0.2, 0.3, 0.4), 0.1, None).unwrap(),
            inertia(1.2, Vec3::new(-0.1, 0.2, 0.15), [0.3, 0.4, 0.5]),
        ),
        Link::new(
            "fixed branch",
            Some(1),
            pose(Vec3::new(0.3, 0.0, -0.2), Vec3::new(-0.2, 0.4, 0.1)),
            JointModel::FIXED,
            inertia(0.6, Vec3::new(0.1, -0.1, 0.2), [0.1; 3]),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(Vec3::new(0.1, -0.3, 0.2), Vec3::new(-0.6, 0.3, 0.2)),
            JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(0.8, Vec3::new(0.2, 0.0, -0.1), [0.1, 0.12, 0.15]),
        ),
    ])
    .unwrap()
}

#[test]
fn mixed_tree_mass_bias_and_jacobians_match_native_energy_rnea_and_aba() {
    let model = mixed_tree();
    let base = FreeFloatingBaseState::new(
        pose(Vec3::new(-0.2, 0.4, 0.3), Vec3::new(1.0, -2.0, 0.8)),
        Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(0.2, -0.5, 0.7)),
    );
    let q = [0.35, -0.2, 0.12];
    let v = [0.2, -0.15, 0.09];
    let terms = free_floating_terms(&model, base, &q, &v, || false).unwrap();
    let mut u = base.twist_body.to_array().to_vec();
    u.extend(v);
    let momentum = product(&terms.mass, &u);
    close(
        0.5 * u.iter().zip(momentum).map(|(a, b)| a * b).sum::<f64>(),
        free_floating_momentum(&model, base, &q, &v)
            .unwrap()
            .kinetic_energy,
    );
    for i in 0..terms.mass.len() {
        for j in 0..terms.mass.len() {
            assert_eq!(terms.mass[i][j], terms.mass[j][i]);
        }
    }
    let base_acceleration = Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(-0.2, 0.4, 0.1));
    let qdd = [0.3, -0.5, 0.2];
    let inverse = inverse_dynamics(
        &model,
        BaseState::prescribed(base.world_from_base, base.twist_body, base_acceleration),
        &q,
        &v,
        &qdd,
        Vec3::default(),
        &[Wrench::default(); 5],
    )
    .unwrap();
    let mut acceleration = base_acceleration.to_array().to_vec();
    acceleration.extend(qdd);
    let base_force =
        inverse.body_wrench[0].transform_by(&model.links()[0].parent_from_child_zero());
    let mut expected = base_force.to_array().to_vec();
    expected.extend(&inverse.generalized_force);
    for ((mass_action, bias), expected) in product(&terms.mass, &acceleration)
        .iter()
        .zip(&terms.bias)
        .zip(expected)
    {
        close(mass_action + bias, expected);
    }
    for i in 0..model.link_count() {
        let apply = |values: &[f64]| {
            terms
                .kinematics
                .body_jacobian(i)
                .iter()
                .zip(values)
                .fold(Twist::zero(), |sum, (column, value)| {
                    sum.plus(column.scale(*value))
                })
        };
        twist_close(apply(&u), terms.kinematics.body_twist(i));
        twist_close(
            apply(&acceleration).plus(terms.kinematics.body_bias_acceleration(i)),
            inverse.body_acceleration[i],
        );
    }
    let mut loads = [Wrench::default(); 5];
    loads[2] = Wrench::new(Vec3::new(0.2, -0.1, 0.3), Vec3::new(-0.4, 0.5, 0.2));
    loads[4] = Wrench::new(Vec3::new(-0.1, 0.2, 0.1), Vec3::new(0.3, -0.2, 0.4));
    let effort = [0.4, -0.2, 0.1];
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let native =
        free_floating_forward_dynamics(&model, base, &q, &v, &effort, gravity, &loads).unwrap();
    let mut solved = native.base_spatial_acceleration_body.to_array().to_vec();
    solved.extend(native.generalized_acceleration);
    let projected = terms.kinematics.project_body_wrenches(&loads).unwrap();
    let gravitational = terms.kinematics.project_uniform_gravity(gravity).unwrap();
    for (i, mass_action) in product(&terms.mass, &solved).iter().enumerate() {
        let actuator = if i < 6 { 0.0 } else { effort[i - 6] };
        close(
            mass_action + terms.bias[i],
            projected[i] + gravitational[i] + actuator,
        );
    }
}

#[test]
fn moving_slider_retains_com_recoil_and_world_load_covariance() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, Vec3::default(), [1.0; 3]),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.3, 0.0, 0.0), [0.1; 3]),
        ),
    ])
    .unwrap();
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let terms = free_floating_terms(&model, base, &[0.2], &[0.0], || false).unwrap();
    close(terms.mass[3][3], 6.0);
    close(terms.mass[3][6], 2.0);
    close(terms.mass[6][6], 2.0);
    close(terms.mass[2][2], 1.0 + 0.1 + 2.0 * 0.5_f64.powi(2));
    close(terms.mass[2][4], 2.0 * 0.5);
    let moved = free_floating_terms(&model, base, &[0.7], &[0.0], || false).unwrap();
    close(moved.mass[2][2], 1.0 + 0.1 + 2.0);
    close(moved.mass[2][4], 2.0);
    // Internal actuator force cannot add total external momentum.
    let acceleration = [0.0, 0.0, 0.0, -0.75, 0.0, 0.0, 2.25];
    for (i, force) in product(&terms.mass, &acceleration).iter().enumerate() {
        close(*force, if i == 6 { 3.0 } else { 0.0 });
    }
    let gravity = terms
        .kinematics
        .project_uniform_gravity(Vec3::new(0.0, -9.81, 0.0))
        .unwrap();
    close(gravity[4], -6.0 * 9.81);
    close(gravity[2], -2.0 * 0.5 * 9.81);
    close(gravity[6], 0.0);

    let world_change = pose(Vec3::new(0.2, -0.4, 0.3), Vec3::new(30.0, -20.0, 10.0));
    let transformed = free_floating_terms(
        &model,
        FreeFloatingBaseState::stationary(world_change),
        &[0.2],
        &[0.0],
        || false,
    )
    .unwrap();
    assert_eq!(terms.mass, transformed.mass);
    let original = [
        Wrench::default(),
        Wrench::new(Vec3::new(0.1, 0.2, -0.1), Vec3::new(0.3, -0.2, 0.4)),
    ];
    let rotated = original.map(|w| {
        Wrench::new(
            world_change.rotation().rotate(w.torque).unwrap(),
            world_change.rotation().rotate(w.force).unwrap(),
        )
    });
    let before = terms.kinematics.project_world_wrenches(&original).unwrap();
    let after = transformed
        .kinematics
        .project_world_wrenches(&rotated)
        .unwrap();
    for (a, b) in before.into_iter().zip(after) {
        close(a, b);
    }
}

fn many_links(count: usize, moving: bool) -> ArticulatedModel {
    let joint = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap();
    ArticulatedModel::new(
        (0..count)
            .map(|i| {
                Link::new(
                    format!("link {i}"),
                    if i == 0 { None } else { Some(0) },
                    Se3::identity(),
                    if i == 0 || !moving {
                        JointModel::FIXED
                    } else {
                        joint
                    },
                    inertia(1.0, Vec3::default(), [1.0; 3]),
                )
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn load_projection_preserves_power_and_keeps_actuator_limits_distinct() {
    let limits = JointLimits::new(-1.0, 1.0, 0.5, 0.2).unwrap();
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, Vec3::default(), [1.0; 3]),
        ),
        Link::new(
            "moving",
            Some(0),
            pose(Vec3::new(0.2, -0.1, 0.3), Vec3::new(0.4, 0.2, -0.1)),
            JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), Some(limits)).unwrap(),
            inertia(2.0, Vec3::new(0.3, 0.1, 0.0), [0.1; 3]),
        ),
    ])
    .unwrap();
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0.2, -0.1, 0.3), Vec3::new(-0.3, 0.4, 0.1)),
    );
    let terms = free_floating_terms(&model, base, &[0.3], &[0.4], || false).unwrap();
    let loads = [
        Wrench::default(),
        Wrench::new(Vec3::new(2.0, -3.0, 4.0), Vec3::new(100.0, -200.0, 50.0)),
    ];
    let projected = terms.kinematics.project_body_wrenches(&loads).unwrap();
    let mut u = base.twist_body.to_array().to_vec();
    u.push(0.4);
    let generalized_power: f64 = projected
        .iter()
        .zip(u)
        .map(|(force, speed)| force * speed)
        .sum();
    let physical_power: f64 = loads
        .iter()
        .enumerate()
        .map(|(i, force)| force.pairing(terms.kinematics.body_twist(i)))
        .sum();
    close(generalized_power, physical_power);
    assert!(model.validate_applied_effort(&[0.2]).is_ok());
    assert!(matches!(
        model.validate_applied_effort(&[0.2_f64.next_up()]),
        Err(ArticulatedError::JointEffortLimitViolation { .. })
    ));
    assert!(free_floating_terms(&model, base, &[0.3], &[0.5_f64.next_up()], || false).is_err());
}

#[test]
fn native_limits_projection_shapes_and_cancellation_are_preserved() {
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let model = mixed_tree();
    let q = [0.1, 0.2, 0.3];
    let v = [0.0; 3];
    let mut polls = 0;
    assert_eq!(
        free_floating_terms(&model, base, &q, &v, || {
            polls += 1;
            polls == 4
        })
        .unwrap_err(),
        GeneralizedError::Cancelled
    );
    assert_eq!(q, [0.1, 0.2, 0.3]);
    assert!(free_floating_terms(&model, base, &[0.0], &v, || false).is_err());
    assert!(free_floating_terms(&model, base, &q, &[f64::NAN, 0.0, 0.0], || false).is_err());
    let terms = free_floating_terms(&model, base, &q, &v, || false).unwrap();
    assert!(terms.kinematics.project_body_wrenches(&[]).is_err());
    assert!(
        terms
            .kinematics
            .project_uniform_gravity(Vec3::new(f64::INFINITY, 0.0, 0.0))
            .is_err()
    );
    let mut nonfixed = model.links().to_vec();
    nonfixed[0] = Link::new(
        "invalid root",
        None,
        Se3::identity(),
        JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
        inertia(1.0, Vec3::default(), [1.0; 3]),
    );
    assert!(matches!(
        free_floating_terms(
            &ArticulatedModel::new(nonfixed).unwrap(),
            base,
            &[0.0; 4],
            &[0.0; 4],
            || false
        ),
        Err(GeneralizedError::Dynamics(
            ArticulatedError::FloatingBaseRootJointNotFixed
        ))
    ));
    let cap_joints = many_links(MAX_GENERALIZED_JOINTS + 1, true);
    assert!(
        free_floating_terms(
            &cap_joints,
            base,
            &[0.0; MAX_GENERALIZED_JOINTS],
            &[0.0; MAX_GENERALIZED_JOINTS],
            || false
        )
        .is_ok()
    );
    let over_joints = many_links(MAX_GENERALIZED_JOINTS + 2, true);
    assert!(matches!(
        free_floating_terms(&over_joints, base, &[], &[], || false),
        Err(GeneralizedError::SizeLimit { .. })
    ));
    let cap_links = many_links(MAX_GENERALIZED_LINKS, false);
    assert!(free_floating_terms(&cap_links, base, &[], &[], || false).is_ok());
    let over_links = many_links(MAX_GENERALIZED_LINKS + 1, false);
    assert!(matches!(
        free_floating_terms(&over_links, base, &[], &[], || false),
        Err(GeneralizedError::SizeLimit { .. })
    ));
}

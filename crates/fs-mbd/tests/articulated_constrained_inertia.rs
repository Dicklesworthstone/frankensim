//! Gear and locked-joint inertia changes, checked against scalar mechanics and
//! individual link COM momentum sums rather than a second generalized solver.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstrainedInertiaReconfiguration, ConstraintControls, ConstraintCoordinate, ConstraintError,
    JointConstraint, transfer_inertia_constrained,
};
use fs_mbd::articulated::motion::reconfigure::{
    ReconfigurationControls, ReconfigurationError, transfer_inertia,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, forward_kinematics,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * (1. + actual.abs().max(expected.abs())),
        "{actual:.16e} != {expected:.16e}",
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}
fn inertia(mass: f64, center: Vec3, diagonal: Vec3, axes: Vec3) -> SpatialInertia {
    let r = rotation(axes).matrix();
    let principal = Mat3 {
        m: [diagonal.x, 0., 0., 0., diagonal.y, 0., 0., 0., diagonal.z],
    };
    SpatialInertia::new(mass, center, r.compose(principal).compose(r.transpose())).unwrap()
}
fn sphere(mass: f64, moment: f64) -> SpatialInertia {
    inertia(
        mass,
        Vec3::default(),
        Vec3::new(moment, moment, moment),
        Vec3::default(),
    )
}
fn with_inertias(model: &ArticulatedModel, values: &[SpatialInertia]) -> ArticulatedModel {
    ArticulatedModel::new(
        model
            .links()
            .iter()
            .zip(values)
            .map(|(link, value)| {
                Link::new(
                    link.name(),
                    link.parent(),
                    link.parent_from_child_zero(),
                    link.joint(),
                    *value,
                )
            })
            .collect(),
    )
    .unwrap()
}
fn controls(tolerances: &[f64], work: f64) -> ConstraintControls<'_> {
    ConstraintControls {
        reconfiguration: ReconfigurationControls {
            relative_momentum_tolerance: 2e-12,
            angular_momentum_tolerance_nms: 2e-10,
            linear_momentum_tolerance_ns: 2e-10,
            joint_momentum_tolerances: tolerances,
            maximum_parameter_work_j: work,
        },
        maximum_reaction_work_j: 1e-8,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn row(coordinate: ConstraintCoordinate, coefficients: Vec<f64>, position: f64) -> JointConstraint {
    JointConstraint {
        coordinate,
        coefficients,
        position,
        position_tolerance: 2e-10,
        velocity_tolerance: 2e-10,
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

struct Observables {
    linear: Vec3,
    angular: Vec3,
    kinetic: f64,
    generalized: Vec<f64>,
}

fn physical(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Observables {
    let kinematics = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let mut linear = Vec3::default();
    let mut angular = Vec3::default();
    let mut about_base = Vec3::default();
    let mut kinetic = 0.;
    let mut momenta = Vec::new();
    for (i, link) in model.links().iter().enumerate() {
        let value = link.inertia();
        let pose = kinematics.world_from_link[i];
        let twist = kinematics.body_twist[i];
        let x = pose.transform_point(value.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(value.center_of_mass()))
            .unwrap();
        let p = velocity.scale(value.mass());
        let spin_body = value.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        linear = linear + p;
        angular = angular + spin + x.cross(p);
        about_base = about_base + spin + (x - state.base.world_from_base.translation()).cross(p);
        kinetic += 0.5 * value.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
        momenta.push((x, p, spin));
    }
    let inverse = state.base.world_from_base.rotation().inverse();
    let local_angular = inverse.rotate(about_base).unwrap();
    let local_linear = inverse.rotate(linear).unwrap();
    let mut generalized = vec![
        local_angular.x,
        local_angular.y,
        local_angular.z,
        local_linear.x,
        local_linear.y,
        local_linear.z,
    ];
    generalized.resize(6 + model.dof_count(), 0.);
    for (joint, link) in model.links().iter().enumerate() {
        let Some(dof) = model.dof_index(joint) else {
            continue;
        };
        let pose = kinematics.world_from_link[joint];
        let screw = link.joint().motion_subspace();
        let angular_axis = pose.rotation().rotate(screw.angular).unwrap();
        let linear_axis = pose.rotation().rotate(screw.linear).unwrap();
        for (body, &(x, p, spin)) in momenta.iter().enumerate() {
            let mut ancestor = Some(body);
            while let Some(index) = ancestor {
                if index == joint {
                    generalized[6 + dof] += angular_axis
                        .dot(spin + (x - pose.translation()).cross(p))
                        + linear_axis.dot(p);
                    break;
                }
                ancestor = model.links()[index].parent();
            }
        }
    }
    Observables {
        linear,
        angular,
        kinetic,
        generalized,
    }
}

fn check_receipt(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    before: &FreeFloatingMotionState,
    rows: &[JointConstraint],
    result: &ConstrainedInertiaReconfiguration,
) {
    let old = physical(previous, before);
    let new = physical(next, &result.state_after);
    assert_eq!(
        result.state_after.base.world_from_base,
        before.base.world_from_base
    );
    assert_eq!(result.state_after.joint_position, before.joint_position);
    vector(new.linear, old.linear);
    vector(new.angular, old.angular);
    vector(result.momentum_before.linear_momentum_world, old.linear);
    vector(result.momentum_after.linear_momentum_world, new.linear);
    vector(result.momentum_before.angular_momentum_world, old.angular);
    vector(result.momentum_after.angular_momentum_world, new.angular);
    close(result.momentum_before.kinetic_energy, old.kinetic);
    close(result.momentum_after.kinetic_energy, new.kinetic);
    close(result.kinetic_energy_change_j, new.kinetic - old.kinetic);
    close(
        result.parameter_work_j + result.reaction_work_j,
        new.kinetic - old.kinetic,
    );
    assert!(result.maximum_normalized_momentum_residual <= 1.);
    assert!(result.reciprocal_condition >= 1e-10);
    assert_eq!(
        result.generalized_constraint_impulse.len(),
        old.generalized.len()
    );
    assert_eq!(result.multipliers.len(), rows.len());
    for i in 0..old.generalized.len() {
        let expected_impulse = if i < 6 {
            0.
        } else {
            rows.iter()
                .zip(&result.multipliers)
                .map(|(row, impulse)| row.coefficients[i - 6] * impulse)
                .sum()
        };
        close(result.generalized_constraint_impulse[i], expected_impulse);
        close(new.generalized[i] - old.generalized[i], expected_impulse);
        close(result.generalized_momentum_before[i], old.generalized[i]);
        close(result.generalized_momentum_after[i], new.generalized[i]);
    }
    let mut reaction_work = 0.;
    for (index, row) in rows.iter().enumerate() {
        let position = dot(&row.coefficients, &before.joint_position) - row.position;
        let old_rate = dot(&row.coefficients, &before.joint_velocity);
        let new_rate = dot(&row.coefficients, &result.state_after.joint_velocity);
        assert!(position.abs() <= row.position_tolerance);
        assert!(old_rate.abs() <= row.velocity_tolerance);
        assert!(new_rate.abs() <= row.velocity_tolerance);
        close(result.position_residuals[index], position);
        close(result.velocity_residuals_before[index], old_rate);
        close(result.velocity_residuals_after[index], new_rate);
        reaction_work += 0.5 * result.multipliers[index] * (old_rate + new_rate);
    }
    close(result.reaction_work_j, reaction_work);
}

fn geared_rotors() -> (ArticulatedModel, ArticulatedModel, FreeFloatingMotionState) {
    let previous = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., 2.),
        ),
        Link::new(
            "first rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            sphere(1., 1.),
        ),
        Link::new(
            "second rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            sphere(1., 3.),
        ),
    ])
    .unwrap();
    let next = with_inertias(&previous, &[sphere(2., 4.), sphere(1., 2.), sphere(1., 1.)]);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0., 0., 0.4), Vec3::default()),
        ),
        joint_position: vec![0.2, -0.4],
        joint_velocity: vec![1.2, -2.4],
    };
    (previous, next, state)
}

#[test]
fn gear_recoil_reactions_work_and_reverse_match_reduced_scalar_mechanics() {
    let (previous, next, state) = geared_rotors();
    let rows = [row(ConstraintCoordinate::Angle, vec![2., 1.], 0.)];
    let settings = controls(&[2e-10; 2], 20.);
    let result =
        transfer_inertia_constrained(&previous, &next, &state, &rows, settings, || false).unwrap();
    // With s=qdot1=-qdot2/2, [p_base,p_s]=[[A,B],[B,C]]*[omega,s].
    // Old (A,B,C)=(6,-5,13); new (7,0,6). The two conserved momenta
    // therefore give omega'=-3.6/7, s'=13.6/6 independently of the API.
    let omega = -3.6 / 7.;
    let rate = 13.6 / 6.;
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., omega),
    );
    close(result.state_after.joint_velocity[0], rate);
    close(result.state_after.joint_velocity[1], -2. * rate);
    let reaction = omega - 2. * rate - (-6.);
    close(result.multipliers[0], reaction);
    close(2. * (omega + rate) - 1.6, 2. * reaction);
    close(
        result.generalized_momentum_after[6] - 2. * result.generalized_momentum_after[7],
        13.6,
    );
    let old_energy = 0.5 * 2. * 0.4_f64.powi(2) + 0.5 * 1.6_f64.powi(2) + 0.5 * 3. * 2_f64.powi(2);
    let new_energy = 0.5 * 7. * omega * omega + 0.5 * 6. * rate * rate;
    close(result.parameter_work_j, new_energy - old_energy);
    close(result.reaction_work_j, 0.);
    check_receipt(&previous, &next, &state, &rows, &result);

    let reversed = transfer_inertia_constrained(
        &next,
        &previous,
        &result.state_after,
        &rows,
        settings,
        || false,
    )
    .unwrap();
    check_receipt(&next, &previous, &result.state_after, &rows, &reversed);
    vector(
        reversed.state_after.base.twist_body.angular,
        state.base.twist_body.angular,
    );
    vector(
        reversed.state_after.base.twist_body.linear,
        state.base.twist_body.linear,
    );
    for (actual, expected) in reversed
        .state_after
        .joint_velocity
        .iter()
        .zip(&state.joint_velocity)
    {
        close(*actual, *expected);
    }
    close(result.parameter_work_j + reversed.parameter_work_j, 0.);

    let scale = -1e6;
    let rescaled = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![2. * scale, scale],
        position: 0.,
        position_tolerance: rows[0].position_tolerance * scale.abs(),
        velocity_tolerance: rows[0].velocity_tolerance * scale.abs(),
    }];
    let equivalent =
        transfer_inertia_constrained(&previous, &next, &state, &rescaled, settings, || false)
            .unwrap();
    vector(
        equivalent.state_after.base.twist_body.angular,
        result.state_after.base.twist_body.angular,
    );
    for (actual, expected) in equivalent
        .state_after
        .joint_velocity
        .iter()
        .zip(&result.state_after.joint_velocity)
    {
        close(*actual, *expected);
    }
    close(equivalent.multipliers[0] * scale, result.multipliers[0]);
    close(equivalent.parameter_work_j, result.parameter_work_j);
    close(equivalent.reciprocal_condition, result.reciprocal_condition);
}

fn serial_tree() -> (ArticulatedModel, ArticulatedModel, FreeFloatingMotionState) {
    let centers = [
        Vec3::new(0.12, -0.21, 0.08),
        Vec3::new(-0.08, 0.17, 0.04),
        Vec3::new(0.11, 0.03, -0.09),
        Vec3::new(0.04, -0.12, 0.06),
    ];
    let masses = [2.3, 1.1, 0.8, 0.4];
    let old_diagonals = [
        Vec3::new(1.2, 1.6, 2.),
        Vec3::new(0.5, 0.7, 0.9),
        Vec3::new(0.3, 0.45, 0.6),
        Vec3::new(0.2, 0.25, 0.3),
    ];
    let new_diagonals = [
        Vec3::new(1.7, 2.2, 2.8),
        Vec3::new(0.8, 1.05, 1.25),
        Vec3::new(0.35, 0.55, 0.7),
        Vec3::new(0.45, 0.55, 0.7),
    ];
    let old: Vec<_> = (0..4)
        .map(|i| {
            inertia(
                masses[i],
                centers[i],
                old_diagonals[i],
                Vec3::new(0.2 + 0.1 * i as f64, -0.15, 0.3),
            )
        })
        .collect();
    let new: Vec<_> = (0..4)
        .map(|i| {
            inertia(
                masses[i],
                centers[i],
                new_diagonals[i],
                Vec3::new(-0.2, 0.35, -0.1 - 0.13 * i as f64),
            )
        })
        .collect();
    let pose = |r, x| Se3::from_parts(rotation(r), x).unwrap();
    let previous = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            pose(Vec3::new(0.2, -0.1, 0.15), Vec3::new(0.17, -0.04, 0.11)),
            JointModel::FIXED,
            old[0],
        ),
        Link::new(
            "rotor",
            Some(0),
            pose(Vec3::new(0.1, 0.25, -0.15), Vec3::new(0.6, 0.08, -0.11)),
            JointModel::revolute(Vec3::new(1., 2., -1.), None).unwrap(),
            old[1],
        ),
        Link::new(
            "slider",
            Some(1),
            pose(Vec3::new(-0.18, 0.11, 0.2), Vec3::new(0.1, -0.35, 0.27)),
            JointModel::prismatic(Vec3::new(-1., 1., 2.), None).unwrap(),
            old[2],
        ),
        Link::new(
            "screw",
            Some(2),
            pose(Vec3::new(0.11, -0.12, 0.18), Vec3::new(-0.2, 0.17, 0.33)),
            JointModel::helical(Vec3::new(2., -1., 1.), 0.13, None).unwrap(),
            old[3],
        ),
    ])
    .unwrap();
    let next = with_inertias(&previous, &new);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(Vec3::new(0.31, -0.27, 0.19), Vec3::new(2., -1., 0.6)),
            Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(-0.2, 0.6, 0.35)),
        ),
        joint_position: vec![0.2, 0.1, 0.15],
        joint_velocity: vec![0.4, 0.2, 0.3],
    };
    (previous, next, state)
}

#[test]
fn mixed_joint_full_tensor_constraints_preserve_physical_momentum_across_frames() {
    let (previous, next, state) = serial_tree();
    // A rack relation (0.5 m/rad) and a second angular ratio leave tangent
    // joint direction [1,0.5,0.75]. Their reactions cannot change its momentum.
    let rows = [
        row(ConstraintCoordinate::Length, vec![-0.5, 1., 0.], 0.),
        row(ConstraintCoordinate::Angle, vec![-0.75, 0., 1.], 0.),
    ];
    let settings = controls(&[2e-10; 3], 100.);
    let result =
        transfer_inertia_constrained(&previous, &next, &state, &rows, settings, || false).unwrap();
    check_receipt(&previous, &next, &state, &rows, &result);
    let old = physical(&previous, &state);
    let new = physical(&next, &result.state_after);
    close(
        dot(&[1., 0.5, 0.75], &new.generalized[6..]),
        dot(&[1., 0.5, 0.75], &old.generalized[6..]),
    );
    assert!(
        (result.state_after.base.twist_body.angular - state.base.twist_body.angular).norm() > 1e-3
    );
    let scene = Se3::from_parts(
        rotation(Vec3::new(-0.4, 0.2, 0.3)),
        Vec3::new(4., -2.5, 1.7),
    )
    .unwrap();
    for boost in [Vec3::default(), Vec3::new(3., -2., 4.)] {
        let mut transformed = state.clone();
        transformed.base.world_from_base = scene.compose(state.base.world_from_base).unwrap();
        let local_boost = transformed
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
        transformed.base.twist_body.linear = state.base.twist_body.linear + local_boost;
        let changed =
            transfer_inertia_constrained(&previous, &next, &transformed, &rows, settings, || false)
                .unwrap();
        check_receipt(&previous, &next, &transformed, &rows, &changed);
        vector(
            changed.state_after.base.twist_body.angular,
            result.state_after.base.twist_body.angular,
        );
        vector(
            changed.state_after.base.twist_body.linear,
            result.state_after.base.twist_body.linear + local_boost,
        );
        for (actual, expected) in changed
            .state_after
            .joint_velocity
            .iter()
            .zip(&result.state_after.joint_velocity)
        {
            close(*actual, *expected);
        }
        for (actual, expected) in changed.multipliers.iter().zip(&result.multipliers) {
            close(*actual, *expected);
        }
        close(changed.parameter_work_j, result.parameter_work_j);
        close(changed.reaction_work_j, result.reaction_work_j);
    }
}

#[test]
fn final_admission_avoids_free_intermediate_limits_and_keeps_failures_atomic() {
    let limits = Some(JointLimits::new(-2., 2., 1., 100.).unwrap());
    let previous = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., 2.),
        ),
        Link::new(
            "locked rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), limits).unwrap(),
            sphere(1., 1.),
        ),
    ])
    .unwrap();
    let next = with_inertias(&previous, &[sphere(2., 0.1), sphere(1., 100.)]);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0., 0., 1.), Vec3::default()),
        ),
        joint_position: vec![0.25],
        joint_velocity: vec![0.],
    };
    let retained = (previous.clone(), next.clone(), state.clone());
    let rows = [row(ConstraintCoordinate::Angle, vec![1.], 0.25)];
    let settings = controls(&[2e-10], 2.);
    assert!(matches!(
        transfer_inertia(&previous, &next, &state, settings.reconfiguration, || false),
        Err(ReconfigurationError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    // The free intermediate would have omega=20, qdot=-19.99 and gain
    // 18.505 J. The actual locked endpoint has omega=3/100.1, qdot=0 and
    // releases about 1.455 J: both physical endpoint limits are admissible.
    let result =
        transfer_inertia_constrained(&previous, &next, &state, &rows, settings, || false).unwrap();
    close(result.state_after.base.twist_body.angular.z, 3. / 100.1);
    close(result.state_after.joint_velocity[0], 0.);
    close(result.parameter_work_j, 4.5 / 100.1 - 1.5);
    check_receipt(&previous, &next, &state, &rows, &result);
    let identity = transfer_inertia_constrained(
        &previous,
        &previous,
        &state,
        &rows,
        controls(&[2e-10], 0.),
        || false,
    )
    .unwrap();
    assert_eq!(identity.state_after, state);
    close(identity.parameter_work_j, 0.);
    close(identity.reaction_work_j, 0.);

    let (free_previous, free_next, free_state) = geared_rotors();
    let free_settings = controls(&[2e-10; 2], 100.);
    let plain = transfer_inertia(
        &free_previous,
        &free_next,
        &free_state,
        free_settings.reconfiguration,
        || false,
    )
    .unwrap();
    let empty = transfer_inertia_constrained(
        &free_previous,
        &free_next,
        &free_state,
        &[],
        free_settings,
        || false,
    )
    .unwrap();
    assert_eq!(empty.state_after, plain.state_after);
    assert_eq!(empty.momentum_before, plain.momentum_before);
    assert_eq!(empty.momentum_after, plain.momentum_after);
    assert_eq!(
        empty.generalized_momentum_before,
        plain.generalized_momentum_before
    );
    assert_eq!(
        empty.generalized_momentum_after,
        plain.generalized_momentum_after
    );
    assert_eq!(empty.parameter_work_j, plain.parameter_work_j);

    for invalid in [
        vec![rows[0].clone(), rows[0].clone()],
        vec![row(ConstraintCoordinate::Angle, vec![0.], 0.)],
    ] {
        assert!(
            transfer_inertia_constrained(&previous, &next, &state, &invalid, settings, || false)
                .is_err()
        );
    }
    let mut wrong_position = state.clone();
    wrong_position.joint_position[0] += 0.01;
    let mut wrong_rate = state.clone();
    wrong_rate.joint_velocity[0] = 0.01;
    for bad in [&wrong_position, &wrong_rate] {
        assert!(
            transfer_inertia_constrained(&previous, &next, bad, &rows, settings, || false).is_err()
        );
    }
    assert!(
        transfer_inertia_constrained(
            &previous,
            &next,
            &state,
            &rows,
            controls(&[2e-10], 0.5),
            || false
        )
        .is_err()
    );

    // A deliberately admitted finite tangency residual gives REAL reported
    // midpoint reaction work; a small reaction-work budget must reject it.
    let mut residual_state = state.clone();
    residual_state.joint_velocity[0] = 0.001;
    let residual_rows = [JointConstraint {
        velocity_tolerance: 0.002,
        ..rows[0].clone()
    }];
    assert!(matches!(
        transfer_inertia_constrained(
            &previous,
            &next,
            &residual_state,
            &residual_rows,
            settings,
            || false
        ),
        Err(ConstraintError::ReactionWork { .. })
    ));
    let admitted = transfer_inertia_constrained(
        &previous,
        &next,
        &residual_state,
        &residual_rows,
        ConstraintControls {
            maximum_reaction_work_j: 0.01,
            ..settings
        },
        || false,
    )
    .unwrap();
    check_receipt(&previous, &next, &residual_state, &residual_rows, &admitted);
    assert!(admitted.reaction_work_j.abs() > 1e-5);

    let mut polls = 0;
    transfer_inertia_constrained(&previous, &next, &state, &rows, settings, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(
            transfer_inertia_constrained(&previous, &next, &state, &rows, settings, || {
                calls += 1;
                calls == stop
            })
            .is_err()
        );
    }
    assert_eq!((previous, next, state), retained);
}

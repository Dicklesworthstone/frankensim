//! Native payload release checked against scalar mechanics, independent COM
//! momentum sums, and a principal-axis impulse/work oracle in the offset case.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::release::{
    PayloadRelease, SeparationImpulse, release_payload,
};
use fs_mbd::articulated::motion::inventory::{
    PayloadCapture, PayloadControls, PayloadError, PayloadState, capture_payload,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Kinematics, Link,
    SpatialInertia, forward_kinematics,
};

#[path = "articulated_payload_release/offset.rs"]
mod offset;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * (1. + actual.abs().max(expected.abs())),
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
fn controls(tolerances: &[f64], maximum_j: f64) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_energy_exchange_j: maximum_j,
        maximum_constraint_work_j: 2e-9,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn kinematics(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Kinematics {
    forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap()
}
fn attached(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target: usize,
    link_from_payload: Se3,
    value: SpatialInertia,
) -> PayloadState {
    let kin = kinematics(model, state);
    PayloadState {
        inertia: value,
        world_from_body: kin.world_from_link[target]
            .compose(link_from_payload)
            .unwrap(),
        twist_body: kin.body_twist[target].transform_by(&link_from_payload.inverse().unwrap()),
    }
}
fn load(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target: usize,
    link_from_payload: Se3,
    value: SpatialInertia,
    rows: &[JointConstraint],
) -> PayloadCapture {
    let payload = attached(model, state, target, link_from_payload, value);
    let tolerances = vec![2e-9; model.dof_count()];
    let result = capture_payload(
        model,
        state,
        target,
        payload,
        rows,
        controls(&tolerances, 1e-6),
        || false,
    )
    .unwrap();
    close(result.removed_relative_energy_j, 0.);
    result
}
fn advance(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    duration: f64,
) -> FreeFloatingMotionState {
    let efforts = vec![0.; model.dof_count()];
    let loads = vec![Wrench::default(); model.link_count()];
    step_free_floating(
        model,
        state,
        FreeFloatingMotionLoads {
            generalized_force: &efforts,
            external_body_wrench: &loads,
            gravity_world: Vec3::default(),
        },
        duration,
        || false,
    )
    .unwrap()
    .state_after
}

struct Observables {
    linear: Vec3,
    angular: Vec3,
    kinetic: f64,
}
fn physical(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    payload: Option<&PayloadState>,
    boost: Vec3,
) -> Observables {
    let kin = kinematics(model, state);
    let mut result = Observables {
        linear: Vec3::default(),
        angular: Vec3::default(),
        kinetic: 0.,
    };
    let bodies = model
        .links()
        .iter()
        .enumerate()
        .map(|(index, link)| PayloadState {
            inertia: link.inertia(),
            world_from_body: kin.world_from_link[index],
            twist_body: kin.body_twist[index],
        });
    for body in bodies.chain(payload.copied()) {
        let value = body.inertia;
        let pose = body.world_from_body;
        let twist = body.twist_body;
        let center = pose.transform_point(value.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(value.center_of_mass()))
            .unwrap();
        let p = velocity.scale(value.mass());
        let spin_body = value.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        result.linear = result.linear + p;
        result.angular = result.angular + spin + center.cross(p);
        let relative = velocity - boost;
        result.kinetic +=
            0.5 * value.mass() * relative.dot(relative) + 0.5 * twist.angular.dot(spin_body);
    }
    result
}
fn check_balance(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    result: &PayloadRelease,
) {
    let boost = result.balance.accounting_boost_world;
    let before = physical(model, state, None, boost);
    let after = physical(
        &result.model_after,
        &result.state_after,
        Some(&result.payload_after),
        boost,
    );
    vector(after.linear, before.linear);
    vector(after.angular, before.angular);
    vector(
        result.balance.momentum_before.linear_momentum_world,
        before.linear,
    );
    vector(
        result.balance.momentum_before.angular_momentum_world,
        before.angular,
    );
    vector(
        result.balance.momentum_after.linear_momentum_world,
        after.linear,
    );
    vector(
        result.balance.momentum_after.angular_momentum_world,
        after.angular,
    );
    close(result.balance.kinetic_before_j, before.kinetic);
    close(result.balance.kinetic_after_j, after.kinetic);
    close(
        after.kinetic - before.kinetic,
        result.supplied_work_j + result.constraint_work_j,
    );
    close(result.balance.energy_residual_j, 0.);
    assert!(result.balance.maximum_normalized_momentum_residual <= 1.);
    for index in 0..result.balance.generalized_momentum_before.len() {
        close(
            result.balance.generalized_momentum_after[index]
                - result.balance.generalized_momentum_before[index],
            result.balance.generalized_constraint_impulse[index],
        );
    }
}
fn coax(
    spin: f64,
    rate: f64,
    speed_limit: Option<f64>,
    boost: Vec3,
) -> (ArticulatedModel, FreeFloatingMotionState) {
    let limits = speed_limit.map(|limit| JointLimits::new(-10., 10., limit, 100.).unwrap());
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., 2.),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), limits).unwrap(),
            sphere(1., 1.),
        ),
    ])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0., 0., spin), boost),
        ),
        joint_position: vec![0.],
        joint_velocity: vec![rate],
    };
    (model, state)
}

#[test]
fn zero_release_preserves_evolved_state_and_point_couple_recoil_retains_small_boosted_work() {
    let (original, initial) = coax(0.4, 1.2, None, Vec3::default());
    let captured = load(&original, &initial, 1, Se3::identity(), sphere(3., 3.), &[]);
    let evolved = advance(&captured.model_after, &captured.state_after, 0.2);
    let inherited = attached(
        &captured.model_after,
        &evolved,
        1,
        captured.attachment.link_from_payload(),
        sphere(3., 3.),
    );
    let zero = release_payload(
        &captured.model_after,
        &evolved,
        &captured.attachment,
        None,
        &[],
        controls(&[2e-9], 0.),
        || false,
    )
    .unwrap();
    assert_eq!(zero.model_after, original);
    assert_eq!(zero.state_after, evolved);
    assert_eq!(zero.payload_after, inherited);
    assert_eq!(zero.supplied_work_j, 0.);
    check_balance(&captured.model_after, &evolved, &zero);
    assert!(
        release_payload(
            &zero.model_after,
            &zero.state_after,
            &captured.attachment,
            None,
            &[],
            controls(&[2e-9], 0.),
            || false
        )
        .is_err()
    );
    advance(&zero.model_after, &zero.state_after, 0.01);
    // Both returned systems remain valid inputs to the original native motion.
    let detached_model = ArticulatedModel::new(vec![Link::new(
        "payload",
        None,
        Se3::identity(),
        JointModel::FIXED,
        zero.payload_after.inertia,
    )])
    .unwrap();
    let detached_state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            zero.payload_after.world_from_body,
            zero.payload_after.twist_body,
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    advance(&detached_model, &detached_state, 0.01);

    let pair = SeparationImpulse {
        point_link_m: Vec3::new(0., 0.2, 0.),
        impulse_on_payload_link_ns: Vec3::new(1.5, 0., 0.),
        couple_on_payload_link_nms: Vec3::new(0., 0., 0.9),
    };
    for boost in [
        Vec3::default(),
        Vec3::new(2_f64.powi(30), -2_f64.powi(29), 2_f64.powi(28)),
    ] {
        let (model, state) = coax(0.4, 1.2, None, boost);
        let loaded = load(&model, &state, 1, Se3::identity(), sphere(3., 3.), &[]);
        let result = release_payload(
            &loaded.model_after,
            &loaded.state_after,
            &loaded.attachment,
            Some(pair),
            &[],
            controls(&[2e-9], 1.),
            || false,
        )
        .unwrap();
        assert_eq!(result.model_after, model);
        vector(
            result.state_after.base.twist_body.angular,
            Vec3::new(0., 0., 0.4),
        );
        close(result.state_after.joint_velocity[0], 0.6);
        // Subtract the exact common boost before comparison: a relative
        // tolerance on billion-scale velocities would hide missing recoil.
        vector(
            result.state_after.base.twist_body.linear - boost,
            Vec3::new(-0.5, 0., 0.),
        );
        vector(
            result.payload_after.twist_body.linear - boost,
            Vec3::new(0.5, 0., 0.),
        );
        vector(
            result.payload_after.twist_body.angular,
            Vec3::new(0., 0., 1.8),
        );
        // r cross J + C = 0.6 z. Combined force/couple recoil gives
        // 0.5*1.5^2*(1/3+1/3) + 0.5*0.6^2*(1/3+1) = 0.99 J.
        close(result.supplied_work_j, 0.99);
        close(result.constraint_work_j, 0.);
        close(result.balance.generalized_momentum_before[2], 7.2);
        close(result.balance.generalized_momentum_before[6], 6.4);
        check_balance(&loaded.model_after, &loaded.state_after, &result);
        assert!(matches!(
            release_payload(
                &loaded.model_after,
                &loaded.state_after,
                &loaded.attachment,
                Some(pair),
                &[],
                controls(&[2e-9], 0.9),
                || false
            ),
            Err(PayloadError::EnergyExchange { .. })
        ));
    }
}

#[test]
fn geared_release_admits_only_the_actual_endpoint_and_keeps_refusals_atomic() {
    let (model, state) = coax(1., 0., Some(0.5), Vec3::default());
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![1.],
        position: 0.,
        position_tolerance: 2e-10,
        velocity_tolerance: 2e-10,
    }];
    let loaded = load(&model, &state, 1, Se3::identity(), sphere(3., 3.), &rows);
    let snapshot = loaded.clone();
    let pair = SeparationImpulse {
        couple_on_payload_link_nms: Vec3::new(0., 0., 0.8),
        ..SeparationImpulse::default()
    };
    // The free intermediate has qdot=-0.8 and exceeds the 0.5 speed limit.
    // Simultaneous lock recoil instead gives qdot=0 and omega=11/15.
    let result = release_payload(
        &loaded.model_after,
        &loaded.state_after,
        &loaded.attachment,
        Some(pair),
        &rows,
        controls(&[2e-9], 0.22),
        || false,
    )
    .unwrap();
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 11. / 15.),
    );
    close(result.state_after.joint_velocity[0], 0.);
    vector(
        result.payload_after.twist_body.angular,
        Vec3::new(0., 0., 19. / 15.),
    );
    close(result.balance.multipliers[0], 8. / 15.);
    close(result.balance.generalized_constraint_impulse[6], 8. / 15.);
    close(result.supplied_work_j, 16. / 75.);
    close(result.constraint_work_j, 0.);
    check_balance(&loaded.model_after, &loaded.state_after, &result);
    advance(&result.model_after, &result.state_after, 0.01);
    assert!(matches!(
        release_payload(
            &loaded.model_after,
            &loaded.state_after,
            &loaded.attachment,
            Some(pair),
            &rows,
            controls(&[2e-9], 0.2),
            || false
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert!(
        release_payload(
            &loaded.model_after,
            &loaded.state_after,
            &loaded.attachment,
            Some(pair),
            &[],
            controls(&[2e-9], 1.),
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    release_payload(
        &loaded.model_after,
        &loaded.state_after,
        &loaded.attachment,
        Some(pair),
        &rows,
        controls(&[2e-9], 1.),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        release_payload(
            &loaded.model_after,
            &loaded.state_after,
            &loaded.attachment,
            Some(pair),
            &rows,
            controls(&[2e-9], 1.),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(PayloadError::Cancelled)
    ));
    let invalid = SeparationImpulse {
        point_link_m: Vec3::new(f64::NAN, 0., 0.),
        ..pair
    };
    assert!(
        release_payload(
            &loaded.model_after,
            &loaded.state_after,
            &loaded.attachment,
            Some(invalid),
            &rows,
            controls(&[2e-9], 1.),
            || false
        )
        .is_err()
    );
    assert_eq!(loaded, snapshot);
}

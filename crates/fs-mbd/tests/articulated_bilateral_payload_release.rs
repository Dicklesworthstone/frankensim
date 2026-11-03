//! G0/G1/G3/G4: retained-inertia ejection through material loops and joint gears.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::bilateral::capture_payload;
use fs_mbd::articulated::motion::inventory::bilateral::release::release_payload;
use fs_mbd::articulated::motion::inventory::release::{
    SeparationImpulse, release_payload as release_original,
};
use fs_mbd::articulated::motion::inventory::{
    PayloadCapture, PayloadControls, PayloadError, PayloadState,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, DistanceLoop, PointCoincidenceLoop, joint,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, forward_kinematics,
};

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 5e-9 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}");
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
fn controls(tolerances: &[f64]) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-8,
        linear_momentum_tolerance_ns: 2e-8,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_energy_exchange_j: 100.,
        maximum_constraint_work_j: 1.,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn point(link: usize) -> MaterialPoint {
    MaterialPoint { link, point_body: Vec3::default() }
}
fn loops() -> [BilateralLoop; 1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: point(1), b: point(2), length_m: 2.,
    })]
}
fn gear() -> JointConstraint {
    JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }
}
fn sliders() -> ArticulatedModel {
    let mut parts = vec![Link::new(
        "carrier", None, Se3::identity(), JointModel::FIXED, sphere(5.),
    )];
    for (i, mass) in [2., 3., 4.].into_iter().enumerate() {
        parts.push(Link::new(
            format!("slider {i}"), Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.),
                Some(JointLimits::new(-10., 10., 1.5, 0.1).unwrap())).unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(parts).unwrap()
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
fn separation(j: f64) -> SeparationImpulse {
    SeparationImpulse {
        point_link_m: Vec3::default(),
        impulse_on_payload_link_ns: Vec3::new(j, 0., 0.),
        couple_on_payload_link_nms: Vec3::default(),
    }
}
#[allow(clippy::too_many_arguments)]
fn attach(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target: usize,
    local: Se3,
    inertia: SpatialInertia,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
) -> PayloadCapture {
    let kin = forward_kinematics(
        model, BaseState::prescribed(
            state.base.world_from_base, state.base.twist_body, Twist::zero(),
        ), &state.joint_position, &state.joint_velocity,
    ).unwrap();
    let payload = PayloadState {
        inertia,
        world_from_body: kin.world_from_link[target].compose(local).unwrap(),
        twist_body: kin.body_twist[target].transform_by(&local.inverse().unwrap()),
    };
    let tolerances = vec![2e-8; model.dof_count()];
    let result = capture_payload(
        model, state, target, payload, loops, rows, controls(&tolerances),
        Default::default(), || false,
    ).unwrap().capture;
    close(result.removed_relative_energy_j, 0.);
    result
}

#[test]
fn captured_geared_linkage_advances_and_ejects_using_retained_inertia() {
    let model = sliders();
    let rows = [gear()];
    let closures = loops();
    let loaded = attach(
        &model, &initial(), 1, Se3::identity(), sphere(2.), &closures, &rows,
    );
    let zero = [Wrench::default(); 4];
    let moved = joint::step_with_bushings(
        &loaded.model_after, &loaded.state_after,
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.; 3],
                external_body_wrench: &zero,
                gravity_world: Vec3::default(),
            }, world_points: &[], connections: &[],
        },
        &closures, &rows, &[], 0.05, Default::default(), || false,
    ).unwrap().bilateral.state_after;
    let saved = moved.clone();
    for j in [4.0, -4.0] {
        let result = release_payload(
            &loaded.model_after, &moved, &loaded.attachment, Some(separation(j)),
            &closures, &rows, controls(&[2e-8; 3]), Default::default(), || false,
        ).unwrap();
        let released = &result.release;
        let reaction = result.reaction.as_ref().unwrap();
        assert_eq!(released.model_after, model);
        assert_eq!(released.state_after.joint_position, moved.joint_position);
        assert_eq!(released.state_after.base.world_from_base, moved.base.world_from_base);
        // Independent retained reduction: t=(2,2,1), effective joint mass10.
        // A target impulse -J gives du=-J/10 and carrier delta=J/35.
        for (i, factor) in [2., 2., 1.].into_iter().enumerate() {
            close(released.state_after.joint_velocity[i],
                moved.joint_velocity[i] - factor * j / 10.);
        }
        close(released.state_after.base.twist_body.linear.x,
            moved.base.twist_body.linear.x + j / 35.);
        close(released.payload_after.twist_body.linear.x,
            moved.base.twist_body.linear.x + moved.joint_velocity[0] + j / 2.);
        close(reaction.loop_multipliers_ns[0], -23. * j / 35.);
        close(reaction.joint_multipliers[0], j / 7.);
        // Point mobility is6/35 in the retained closed mechanism, plus1/2
        // for the separate2kg payload. Initial source-pair work cancels.
        close(released.supplied_work_j, 47. * j * j / 140.);
        close(released.constraint_work_j, 0.);
        close(released.balance.energy_residual_j, 0.);
        vector(released.balance.momentum_after.linear_momentum_world,
            released.balance.momentum_before.linear_momentum_world);
        vector(released.balance.momentum_after.angular_momentum_world,
            released.balance.momentum_before.angular_momentum_world);
        close(released.balance.kinetic_after_j - released.balance.kinetic_before_j,
            released.supplied_work_j);
        for value in &reaction.velocity_residual_after_mps { close(*value, 0.); }
        for value in &reaction.joint_velocity_residuals_after { close(*value, 0.); }
        assert!(released.balance.maximum_normalized_momentum_residual <= 1.);
        // At J=-4, the free first-slider candidate would be2.4 (>1.5);
        // the combined physical result is1.2 and must be admitted.
        if j < 0. { close(released.state_after.joint_velocity[0], 1.2); }
    }
    assert_eq!(moved, saved);
}

#[test]
fn passive_release_preserves_finite_slip_and_empty_loops_keep_legacy_results() {
    let model = sliders();
    let mut rows = [gear()];
    rows[0].velocity_tolerance = 2e-3;
    let loaded = attach(
        &model, &initial(), 1, Se3::identity(), sphere(2.), &loops(), &rows,
    );
    let mut state = loaded.state_after.clone();
    state.joint_velocity[2] -= 0.0005;
    state.base.twist_body.angular.z = -0.0;
    for impulse in [None, Some(SeparationImpulse::default())] {
        let result = release_payload(
            &loaded.model_after, &state, &loaded.attachment, impulse,
            &loops(), &rows, controls(&[2e-8; 3]), Default::default(), || false,
        ).unwrap();
        assert_eq!(result.release.state_after, state);
        for (a, b) in result.release.state_after.base.twist_body.to_array().into_iter()
            .zip(state.base.twist_body.to_array())
        { assert_eq!(a.to_bits(), b.to_bits()); }
        assert_eq!(result.release.model_after, model);
        close(result.release.supplied_work_j, 0.);
        close(result.release.constraint_work_j, 0.);
        let reaction = result.reaction.unwrap();
        close(reaction.joint_velocity_residuals_after[0], 0.001);
        assert_eq!(reaction.joint_velocity_residuals_after,
            reaction.joint_velocity_residuals_before);
        assert_eq!(reaction.loop_multipliers_ns, vec![0.]);
        assert_eq!(reaction.joint_multipliers, vec![0.]);
    }
    let strict_condition = PayloadControls {
        minimum_reciprocal_condition: 0.99, ..controls(&[2e-8; 3])
    };
    assert!(matches!(release_payload(
        &loaded.model_after, &state, &loaded.attachment, None,
        &loops(), &rows, strict_condition, Default::default(), || false,
    ), Err(PayloadError::Constraint(ConstraintError::IllConditioned { .. }))));
    let old = release_original(
        &loaded.model_after, &loaded.state_after, &loaded.attachment,
        Some(separation(0.3)), &rows, controls(&[2e-8; 3]), || false,
    ).unwrap();
    let empty = release_payload(
        &loaded.model_after, &loaded.state_after, &loaded.attachment,
        Some(separation(0.3)), &[], &rows, controls(&[2e-8; 3]),
        Default::default(), || false,
    ).unwrap();
    assert!(empty.reaction.is_none());
    assert_eq!(empty.release.state_after, old.state_after);
    assert_eq!(empty.release.model_after, old.model_after);
    assert_eq!(empty.release.payload_after, old.payload_after);
    assert_eq!(empty.release.balance, old.balance);
    assert_eq!(empty.release.supplied_work_j, old.supplied_work_j);
}

fn assembled(parts: &[SpatialInertia]) -> SpatialInertia {
    let mass: f64 = parts.iter().map(|p| p.mass()).sum();
    let center = parts.iter().fold(Vec3::default(), |c, p|
        c + p.center_of_mass().scale(p.mass() / mass));
    let mut a = [0.; 9];
    for p in parts {
        let d = p.center_of_mass() - center;
        let v = [d.x, d.y, d.z];
        for i in 0..3 { for j in 0..3 {
            a[3 * i + j] += p.inertia_com().m[3 * i + j]
                + p.mass() * ((if i == j { d.dot(d) } else { 0. }) - v[i] * v[j]);
        } }
    }
    SpatialInertia::new(mass, center, Mat3 { m: a }).unwrap()
}
fn solve3(a: Mat3, b: Vec3) -> Vec3 {
    // Independent cofactor inverse, without the articulated or fs-la solver.
    let m = a.m;
    let determinant = m[0] * (m[4] * m[8] - m[5] * m[7])
        - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    Vec3::new(
        ((m[4] * m[8] - m[5] * m[7]) * b.x
            + (m[2] * m[7] - m[1] * m[8]) * b.y
            + (m[1] * m[5] - m[2] * m[4]) * b.z) / determinant,
        ((m[5] * m[6] - m[3] * m[8]) * b.x
            + (m[0] * m[8] - m[2] * m[6]) * b.y
            + (m[2] * m[3] - m[0] * m[5]) * b.z) / determinant,
        ((m[3] * m[7] - m[4] * m[6]) * b.x
            + (m[1] * m[6] - m[0] * m[7]) * b.y
            + (m[0] * m[4] - m[1] * m[3]) * b.z) / determinant,
    )
}
fn rigid_response(inertia: SpatialInertia, impulse: Wrench) -> Twist {
    let angular = solve3(inertia.inertia_com(),
        impulse.torque - inertia.center_of_mass().cross(impulse.force));
    Twist::new(angular,
        impulse.force.scale(1. / inertia.mass())
            - angular.cross(inertia.center_of_mass()))
}

#[test]
fn spherical_lock_ejection_matches_independent_full_tensor_responses_and_pair_work() {
    let parts: Vec<_> = [5., 1., 1., 1., 2.].into_iter().enumerate().map(|(i, m)| {
        let f = i as f64;
        SpatialInertia::new(m, Vec3::new(0.03 * f, -0.02 * f, 0.01 * f),
            Mat3 { m: [m, 0.1*m, -0.05*m, 0.1*m, 1.2*m, 0.08*m,
                -0.05*m, 0.08*m, 1.5*m] }).unwrap()
    }).collect();
    let axes = [Vec3::new(1., 0., 0.), Vec3::new(0., 1., 0.), Vec3::new(0., 0., 1.)];
    let mut links = vec![Link::new(
        "base", None, Se3::identity(), JointModel::FIXED, parts[0],
    )];
    for i in 0..3 {
        links.push(Link::new(format!("slide {i}"), Some(i), Se3::identity(),
            JointModel::prismatic(axes[i], None).unwrap(), parts[i + 1]));
    }
    links.push(Link::new("rotor", Some(3), Se3::identity(),
        JointModel::revolute(axes[2], None).unwrap(), parts[4]));
    let model = ArticulatedModel::new(links).unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: vec![0.; 4], joint_velocity: vec![0.; 4],
    };
    let closures = [BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a: point(0), b: point(3),
    })];
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle, coefficients: vec![0., 0., 0., 1.],
        position: 0., position_tolerance: 1e-9, velocity_tolerance: 1e-9,
    }];
    let payload_inertia = SpatialInertia::new(
        2., Vec3::new(-0.1, 0.07, 0.03),
        Mat3 { m: [0.9, 0.08, -0.04, 0.08, 1.1, 0.06, -0.04, 0.06, 1.4] },
    ).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let offset = Vec3::new(0.3, -0.2, 0.4);
    let loaded = attach(
        &model, &state, 4, Se3::from_parts(rotation, offset).unwrap(),
        payload_inertia, &closures, &rows,
    );
    let applied = SeparationImpulse {
        point_link_m: Vec3::new(0.2, -0.1, 0.25),
        impulse_on_payload_link_ns: Vec3::new(0.7, -0.4, 0.3),
        couple_on_payload_link_nms: Vec3::new(0.08, 0.03, -0.06),
    };
    let result = release_payload(
        &loaded.model_after, &loaded.state_after, &loaded.attachment, Some(applied),
        &closures, &rows, controls(&[2e-8; 4]), Default::default(), || false,
    ).unwrap();
    let retained = assembled(&parts);
    let at_link = Wrench::new(
        applied.point_link_m.cross(applied.impulse_on_payload_link_ns)
            + applied.couple_on_payload_link_nms,
        applied.impulse_on_payload_link_ns,
    );
    let retained_impulse = Wrench::new(at_link.torque.scale(-1.), at_link.force.scale(-1.));
    let payload_impulse = Wrench::new(
        rotation.inverse().rotate(
            (applied.point_link_m - offset).cross(applied.impulse_on_payload_link_ns)
                + applied.couple_on_payload_link_nms,
        ).unwrap(),
        rotation.inverse().rotate(applied.impulse_on_payload_link_ns).unwrap(),
    );
    let expected_base = rigid_response(retained, retained_impulse);
    let expected_payload = rigid_response(payload_inertia, payload_impulse);
    vector(result.release.state_after.base.twist_body.angular, expected_base.angular);
    vector(result.release.state_after.base.twist_body.linear, expected_base.linear);
    vector(result.release.payload_after.twist_body.angular, expected_payload.angular);
    vector(result.release.payload_after.twist_body.linear, expected_payload.linear);
    for rate in &result.release.state_after.joint_velocity { close(*rate, 0.); }
    close(result.release.supplied_work_j,
        0.5 * (retained_impulse.pairing(expected_base) + payload_impulse.pairing(expected_payload)));
    close(result.release.constraint_work_j, 0.);
    close(result.release.balance.energy_residual_j, 0.);
    vector(result.release.balance.momentum_after.linear_momentum_world, Vec3::default());
    vector(result.release.balance.momentum_after.angular_momentum_world, Vec3::default());
    assert_eq!(result.release.model_after, model);
    let reaction = result.reaction.unwrap();
    assert_eq!(reaction.rows.len(), 3);
    assert_eq!(reaction.joint_multipliers.len(), 1);
    assert!(reaction.joint_multipliers[0].abs() > 1e-6);
}

#[test]
fn finite_constraint_work_and_independent_budgets_keep_signed_release_physics() {
    let model = sliders();
    let mut rows = [gear()];
    rows[0].velocity_tolerance = 2e-3;
    let loaded = attach(
        &model, &initial(), 1, Se3::identity(), sphere(2.), &loops(), &rows,
    );
    let mut state = loaded.state_after.clone();
    let defect = 1e-3;
    state.joint_velocity[2] -= 0.5 * defect;
    let saved = state.clone();
    for j in [4., -0.0001] {
        let result = release_payload(
            &loaded.model_after, &state, &loaded.attachment, Some(separation(j)),
            &loops(), &rows, controls(&[2e-8; 3]), Default::default(), || false,
        ).unwrap();
        // Solve the independent2x2 row system with before rate[0,defect].
        let lambda = (j - 5. * defect) / 7.;
        let expected_constraint = 0.5 * lambda * defect;
        let expected_source = 47. * j * j / 140. + j * defect / 14.;
        let reaction = result.reaction.as_ref().unwrap();
        close(reaction.loop_multipliers_ns[0], -23. * j / 35. + 2. * defect / 7.);
        close(reaction.joint_multipliers[0], lambda);
        assert!((result.release.constraint_work_j - expected_constraint).abs() < 2e-11);
        assert!((result.release.supplied_work_j - expected_source).abs() < 2e-11);
        close(result.release.balance.kinetic_after_j - result.release.balance.kinetic_before_j,
            expected_source + expected_constraint);
        close(result.release.balance.energy_residual_j, 0.);
        if j < 0. {
            // Free-release policy bounds absolute work and permits extraction.
            assert!(result.release.supplied_work_j < 0.);
        }
    }
    let reaction_cap = PayloadControls {
        maximum_constraint_work_j: 1e-5, ..controls(&[2e-8; 3])
    };
    assert!(matches!(release_payload(
        &loaded.model_after, &state, &loaded.attachment, Some(separation(4.)),
        &loops(), &rows, reaction_cap, Default::default(), || false,
    ), Err(PayloadError::ConstraintWork { .. })));
    let source_cap = PayloadControls {
        maximum_energy_exchange_j: 1., ..controls(&[2e-8; 3])
    };
    assert!(matches!(release_payload(
        &loaded.model_after, &state, &loaded.attachment, Some(separation(4.)),
        &loops(), &rows, source_cap, Default::default(), || false,
    ), Err(PayloadError::EnergyExchange { .. })));
    assert!(matches!(release_payload(
        &loaded.model_after, &state, &loaded.attachment, Some(separation(4.)),
        &loops(), &rows, controls(&[2e-8; 3]), Default::default(), || true,
    ), Err(PayloadError::Cancelled)));
    assert_eq!(state, saved);
}

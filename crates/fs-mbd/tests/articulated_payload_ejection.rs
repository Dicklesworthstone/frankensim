//! Prescribed final relative motion, checked against independent mechanics.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::inventory::bilateral::release::{
    eject_payloads as eject_mixed, release_payloads as release_mixed,
};
use fs_mbd::articulated::motion::inventory::release::ejection::{
    MAX_PAYLOAD_EJECTION_TARGETS, PayloadEjection, PayloadEjectionControls, PayloadEjectionError,
    PayloadEjectionInput, eject_payloads,
};
use fs_mbd::articulated::motion::inventory::release::{
    PayloadBatchRelease, PayloadReleaseInput, release_payloads,
};
use fs_mbd::articulated::motion::inventory::{
    PayloadBatchCapture, PayloadCaptureInput, PayloadControls, PayloadError, PayloadState,
    capture_payloads,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::MaterialPoint;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop, DistanceLoop,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Kinematics, Link,
    SpatialInertia, forward_kinematics,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, allowance {tolerance}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn twist(a: Twist, b: Twist, tolerance: f64) {
    vector(a.angular, b.angular, tolerance);
    vector(a.linear, b.linear, tolerance);
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn inertia(mass: f64, center: Vec3, principal: Vec3, axes: Vec3) -> SpatialInertia {
    let r = rotation(axes).matrix();
    let diagonal = Mat3 {
        m: [
            principal.x,
            0.,
            0.,
            0.,
            principal.y,
            0.,
            0.,
            0.,
            principal.z,
        ],
    };
    SpatialInertia::new(mass, center, r.compose(diagonal).compose(r.transpose())).unwrap()
}
fn sphere(mass: f64, moment: f64) -> SpatialInertia {
    inertia(
        mass,
        Vec3::default(),
        Vec3::new(moment, moment, moment),
        Vec3::default(),
    )
}
fn limits(joints: &[f64], maximum_work: f64) -> PayloadEjectionControls<'_> {
    PayloadEjectionControls {
        payload: PayloadControls {
            relative_momentum_tolerance: 2e-12,
            angular_momentum_tolerance_nms: 2e-8,
            linear_momentum_tolerance_ns: 2e-8,
            joint_momentum_tolerances: joints,
            energy_tolerance_j: 2e-8,
            maximum_energy_exchange_j: maximum_work,
            maximum_constraint_work_j: 0.01,
            minimum_reciprocal_condition: 1e-10,
        },
        relative_velocity_tolerance_m_s: 2e-10,
        relative_angular_velocity_tolerance_rad_s: 2e-10,
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
// Carrier velocity evaluated at the PAYLOAD ORIGIN, expressed in payload axes.
fn carried(value: Twist, local: Se3) -> Twist {
    let inverse = local.rotation().inverse();
    Twist::new(
        inverse.rotate(value.angular).unwrap(),
        inverse
            .rotate(value.linear + value.angular.cross(local.translation()))
            .unwrap(),
    )
}
fn load(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    parts: &[(usize, Se3, SpatialInertia)],
) -> PayloadBatchCapture {
    let k = kinematics(model, state);
    let inputs: Vec<_> = parts
        .iter()
        .map(|&(target_link, local, value)| PayloadCaptureInput {
            target_link,
            payload: PayloadState {
                inertia: value,
                world_from_body: k.world_from_link[target_link].compose(local).unwrap(),
                twist_body: carried(k.body_twist[target_link], local),
            },
        })
        .collect();
    let allowances = vec![2e-8; model.dof_count()];
    // Build already-carried inventory without projecting its authored rates.
    // Ejection must encounter the finite initial constraint defects itself;
    // supplying rows here would erase those defects during capture setup.
    let loaded = capture_payloads(
        model,
        state,
        &inputs,
        &[],
        limits(&allowances, 1e-6).payload,
        || false,
    )
    .unwrap();
    near(loaded.removed_relative_energy_j, 0., 1e-17);
    assert_eq!(loaded.state_after.joint_position, state.joint_position);
    twist(
        loaded.state_after.base.twist_body,
        state.base.twist_body,
        2e-12,
    );
    for (&actual, &authored) in loaded
        .state_after
        .joint_velocity
        .iter()
        .zip(&state.joint_velocity)
    {
        near(actual, authored, 2e-12);
    }
    loaded
}
fn scalar(boost: Vec3) -> PayloadBatchCapture {
    let model = ArticulatedModel::new(vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(4., 4.),
    )])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0.25, -0.125, 0.5), boost),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    load(
        &model,
        &state,
        &[
            (0, Se3::identity(), sphere(1., 1.)),
            (0, Se3::identity(), sphere(2., 2.)),
        ],
    )
}
fn targets<'a>(
    loaded: &'a PayloadBatchCapture,
    requested: [Twist; 2],
) -> [PayloadEjectionInput<'a>; 2] {
    [
        PayloadEjectionInput {
            attachment: &loaded.attachments[1],
            relative_twist_payload: requested[1],
        },
        PayloadEjectionInput {
            attachment: &loaded.attachments[0],
            relative_twist_payload: requested[0],
        },
    ]
}
fn assert_relative(result: &PayloadEjection, inputs: &[PayloadEjectionInput<'_>], tolerance: f64) {
    let release = &result.release;
    let k = kinematics(&release.model_after, &release.state_after);
    for (part, input) in release.payloads_after.iter().zip(inputs) {
        let carrier = carried(
            k.body_twist[input.attachment.target_link()],
            input.attachment.link_from_payload(),
        );
        twist(
            part.payload_after.twist_body.minus(carrier),
            input.relative_twist_payload,
            tolerance,
        );
    }
}
fn kinetic(value: SpatialInertia, velocity: Twist) -> f64 {
    let com = velocity.linear + velocity.angular.cross(value.center_of_mass());
    0.5 * (value.mass() * com.dot(com)
        + velocity
            .angular
            .dot(value.inertia_com().apply(velocity.angular)))
}
fn physical(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    extras: &[PayloadState],
    boost: Vec3,
) -> (Vec3, Vec3, f64) {
    let k = kinematics(model, state);
    let bodies = model
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| PayloadState {
            inertia: link.inertia(),
            world_from_body: k.world_from_link[i],
            twist_body: k.body_twist[i],
        });
    let mut p = Vec3::default();
    let mut l = Vec3::default();
    let mut energy = 0.;
    for body in bodies.chain(extras.iter().copied()) {
        let value = body.inertia;
        let center = body
            .world_from_body
            .transform_point(value.center_of_mass())
            .unwrap();
        let velocity = body
            .world_from_body
            .rotation()
            .rotate(body.twist_body.linear + body.twist_body.angular.cross(value.center_of_mass()))
            .unwrap();
        let linear = velocity.scale(value.mass());
        let spin_body = value.inertia_com().apply(body.twist_body.angular);
        let spin = body.world_from_body.rotation().rotate(spin_body).unwrap();
        p = p + linear;
        l = l + spin + center.cross(linear);
        let relative = velocity - boost;
        energy +=
            0.5 * (value.mass() * relative.dot(relative) + body.twist_body.angular.dot(spin_body));
    }
    (p, l, energy)
}
fn check_balance(loaded: &PayloadBatchCapture, release: &PayloadBatchRelease, tolerance: f64) {
    let boost = loaded
        .state_after
        .base
        .world_from_base
        .rotation()
        .rotate(loaded.state_after.base.twist_body.linear)
        .unwrap();
    let before = physical(&loaded.model_after, &loaded.state_after, &[], boost);
    let extras: Vec<_> = release
        .payloads_after
        .iter()
        .map(|p| p.payload_after)
        .collect();
    let after = physical(&release.model_after, &release.state_after, &extras, boost);
    vector(after.0, before.0, tolerance);
    vector(after.1, before.1, tolerance);
    near(
        after.2 - before.2,
        release.supplied_work_j + release.constraint_work_j,
        tolerance,
    );
}

#[test]
fn g1_simultaneous_final_targets_share_recoil_spin_and_work_under_common_boosts() {
    let command = [
        Twist::new(Vec3::new(0., 0., 2.), Vec3::new(3., 0., 0.)),
        Twist::new(Vec3::new(0., 0., -0.5), Vec3::new(-1., 0., 0.)),
    ];
    for boost in [
        Vec3::default(),
        Vec3::new(2_f64.powi(30), -2_f64.powi(29), 2_f64.powi(28)),
    ] {
        let loaded = scalar(boost);
        let inputs = targets(&loaded, command);
        let mut controls = limits(&[], 8.);
        controls.relative_velocity_tolerance_m_s = 5e-7;
        controls.payload.energy_tolerance_j = 2e-6;
        let result = eject_payloads(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &[],
            controls,
            || false,
        )
        .unwrap();
        let initial = loaded.state_after.base.twist_body;
        vector(
            result.release.state_after.base.twist_body.linear - initial.linear,
            Vec3::new(-1. / 7., 0., 0.),
            2e-7,
        );
        vector(
            result.release.state_after.base.twist_body.angular - initial.angular,
            Vec3::new(0., 0., -1. / 7.),
            2e-12,
        );
        for (part, i) in result.release.payloads_after.iter().zip([1, 0]) {
            vector(
                part.impulse_on_payload_body.force,
                Vec3::new(if i == 0 { 20. / 7. } else { -16. / 7. }, 0., 0.),
                2e-11,
            );
            vector(
                part.impulse_on_payload_body.torque,
                Vec3::new(0., 0., if i == 0 { 13. / 7. } else { -9. / 7. }),
                2e-11,
            );
        }
        near(result.release.supplied_work_j, 213. / 28., 2e-6);
        near(result.release.constraint_work_j, 0., 0.);
        assert_relative(&result, &inputs, 5e-7);
        check_balance(&loaded, &result.release, 2e-6);
    }
    // A sequential ejection commands different events: the later carrier
    // recoil changes the earlier payload's relative motion. The simultaneous
    // analytical endpoint above cannot be recovered by independent calls.
    let loaded = scalar(Vec3::default());
    let inputs = targets(&loaded, command);
    let first = eject_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &inputs[..1],
        &[],
        limits(&[], 10.),
        || false,
    )
    .unwrap();
    let second = eject_payloads(
        &first.release.model_after,
        &first.release.state_after,
        &inputs[1..],
        &[],
        limits(&[], 10.),
        || false,
    )
    .unwrap();
    near(
        second.release.state_after.base.twist_body.linear.x,
        -11. / 35.,
        2e-11,
    );
    near(
        second.release.state_after.base.twist_body.angular.z
            - loaded.state_after.base.twist_body.angular.z,
        -9. / 35.,
        2e-11,
    );
    let early = first.release.payloads_after[0].payload_after.twist_body;
    near(
        early.linear.x - second.release.state_after.base.twist_body.linear.x,
        -0.4,
        2e-11,
    );
}

// The independent free-carrier oracle uses ordinary COM linear momentum and
// a dense 3x3 centroidal inertia, never the production spatial/ABA inverse.
fn determinant(a: [f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}
fn solve3(a: [f64; 9], rhs: Vec3) -> Vec3 {
    let v = [rhs.x, rhs.y, rhs.z];
    let determinant = determinant(a);
    assert!(determinant > 0.);
    let mut result = [0.; 3];
    for column in 0..3 {
        let mut replaced = a;
        for row in 0..3 {
            replaced[3 * row + column] = v[row];
        }
        result[column] = self::determinant(replaced) / determinant;
    }
    Vec3::new(result[0], result[1], result[2])
}
fn carrier_recoil(values: [SpatialInertia; 3], local: [Se3; 2], command: [Twist; 2]) -> Twist {
    let masses = values.map(|v| v.mass());
    let centers = [
        values[0].center_of_mass(),
        local[0]
            .transform_point(values[1].center_of_mass())
            .unwrap(),
        local[1]
            .transform_point(values[2].center_of_mass())
            .unwrap(),
    ];
    let total_mass: f64 = masses.iter().sum();
    let center = centers
        .into_iter()
        .zip(masses)
        .fold(Vec3::default(), |s, (r, m)| s + r.scale(m / total_mass));
    let mut inertia = [0.; 9];
    for i in 0..3 {
        let r = centers[i] - center;
        let xyz = [r.x, r.y, r.z];
        let axes = if i == 0 {
            So3::identity()
        } else {
            local[i - 1].rotation()
        };
        let rotation = axes.matrix();
        let rotated = rotation
            .compose(values[i].inertia_com())
            .compose(rotation.transpose());
        for row in 0..3 {
            for col in 0..3 {
                inertia[3 * row + col] += rotated.m[3 * row + col]
                    + masses[i] * (if row == col { r.dot(r) } else { 0. } - xyz[row] * xyz[col]);
            }
        }
    }
    let mut linear_bias = Vec3::default();
    let mut angular_bias = Vec3::default();
    for i in 0..2 {
        let value = values[i + 1];
        let dc = command[i].linear + command[i].angular.cross(value.center_of_mass());
        let p = local[i].rotation().rotate(dc.scale(value.mass())).unwrap();
        let l = local[i]
            .rotation()
            .rotate(value.inertia_com().apply(command[i].angular))
            .unwrap();
        linear_bias = linear_bias + p;
        angular_bias = angular_bias + l + (centers[i + 1] - center).cross(p);
    }
    let angular = solve3(inertia, angular_bias.scale(-1.));
    Twist::new(
        angular,
        linear_bias.scale(-1. / total_mass) - angular.cross(center),
    )
}

#[test]
fn g1_g3_offcenter_full_tensors_resolve_payload_origin_targets_and_scene_covariance() {
    let values = [
        inertia(
            3.5,
            Vec3::new(0.15, -0.2, 0.1),
            Vec3::new(1.2, 1.6, 2.),
            Vec3::new(0.2, -0.3, 0.1),
        ),
        inertia(
            1.25,
            Vec3::new(-0.1, 0.15, 0.05),
            Vec3::new(0.3, 0.4, 0.5),
            Vec3::new(-0.1, 0.4, 0.2),
        ),
        inertia(
            0.75,
            Vec3::new(0.1, -0.05, 0.2),
            Vec3::new(0.2, 0.25, 0.3),
            Vec3::new(0.3, 0.2, -0.4),
        ),
    ];
    let local = [
        Se3::from_parts(
            rotation(Vec3::new(0.2, 0.5, -0.3)),
            Vec3::new(-0.7, 0.2, 0.3),
        )
        .unwrap(),
        Se3::from_parts(
            rotation(Vec3::new(-0.4, 0.1, 0.5)),
            Vec3::new(0.5, -0.6, 0.4),
        )
        .unwrap(),
    ];
    let command = [
        Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(-0.6, 0.2, 0.5)),
        Twist::new(Vec3::new(-0.2, 0.5, -0.4), Vec3::new(0.3, -0.5, 0.1)),
    ];
    let recoil = carrier_recoil(values, local, command);
    let mut expected_work = kinetic(values[0], recoil);
    for i in 0..2 {
        expected_work += kinetic(values[i + 1], carried(recoil, local[i]).plus(command[i]));
    }
    assert!(expected_work > 0.1);
    for (scene, boost) in [
        (Se3::identity(), Vec3::default()),
        (
            Se3::from_parts(rotation(Vec3::new(0.4, -0.2, 0.3)), Vec3::new(4., -3., 2.)).unwrap(),
            Vec3::new(0.3, -0.4, 0.1),
        ),
    ] {
        let model = ArticulatedModel::new(vec![Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            values[0],
        )])
        .unwrap();
        let state = FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                scene,
                Twist::new(
                    Vec3::new(0.2, -0.3, 0.4),
                    Vec3::new(-0.1, 0.4, 0.2) + scene.rotation().inverse().rotate(boost).unwrap(),
                ),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        };
        let loaded = load(
            &model,
            &state,
            &[(0, local[0], values[1]), (0, local[1], values[2])],
        );
        let inputs = targets(&loaded, command);
        let result = eject_payloads(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &[],
            limits(&[], 100.),
            || false,
        )
        .unwrap();
        assert_eq!(result.release.model_after, model);
        twist(
            result.release.state_after.base.twist_body,
            loaded.state_after.base.twist_body.plus(recoil),
            2e-10,
        );
        for (part, i) in result.release.payloads_after.iter().zip([1, 0]) {
            let delta = carried(recoil, local[i]).plus(command[i]);
            let original = carried(loaded.state_after.base.twist_body, local[i]);
            twist(part.payload_after.twist_body, original.plus(delta), 3e-10);
            let value = values[i + 1];
            let j =
                (delta.linear + delta.angular.cross(value.center_of_mass())).scale(value.mass());
            let k = value.inertia_com().apply(delta.angular) + value.center_of_mass().cross(j);
            vector(part.impulse_on_payload_body.force, j, 2e-10);
            vector(part.impulse_on_payload_body.torque, k, 2e-10);
            // The COM relative velocity includes delta-omega cross COM.
            // Equating it with the origin command would solve another problem.
            let carrier = carried(result.release.state_after.base.twist_body, local[i]);
            let relative = part.payload_after.twist_body.minus(carrier);
            let com_relative = relative.linear + relative.angular.cross(value.center_of_mass());
            vector(
                com_relative,
                command[i].linear + command[i].angular.cross(value.center_of_mass()),
                3e-10,
            );
            assert!((com_relative - command[i].linear).norm() > 0.01);
        }
        assert_relative(&result, &inputs, 3e-10);
        near(result.release.supplied_work_j, expected_work, 3e-10);
        check_balance(&loaded, &result.release, 2e-9);
    }
}

fn linkage(
    defect: f64,
) -> (
    ArticulatedModel,
    FreeFloatingMotionState,
    [JointConstraint; 1],
    [BilateralLoop; 1],
) {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(5., 5.),
    )];
    for (index, mass) in [2., 3., 4.].into_iter().enumerate() {
        links.push(Link::new(
            format!("slider {index}"),
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            sphere(mass, mass),
        ));
    }
    let model = ArticulatedModel::new(links).unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2 - 0.5 * defect],
    };
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 2e-3,
    }];
    let loops = [BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        length_m: 2.,
    })];
    (model, state, rows, loops)
}

#[test]
fn g1_mixed_constraints_retain_the_original_tangent_defect_in_actual_final_targets() {
    for defect in [0., 1e-3] {
        let (model, state, rows, loops) = linkage(defect);
        let loaded = load(
            &model,
            &state,
            &[
                (1, Se3::identity(), sphere(2., 2.)),
                (3, Se3::identity(), sphere(1., 1.)),
            ],
        );
        let inputs = [
            PayloadEjectionInput {
                attachment: &loaded.attachments[0],
                relative_twist_payload: Twist::new(Vec3::default(), Vec3::new(3., 0., 0.)),
            },
            PayloadEjectionInput {
                attachment: &loaded.attachments[1],
                relative_twist_payload: Twist::new(Vec3::default(), Vec3::new(-1., 0., 0.)),
            },
        ];
        let config = BilateralImpulseConfig {
            velocity_tolerance_mps: 2e-3,
            ..Default::default()
        };
        let result = eject_mixed(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &loops,
            &rows,
            limits(&[2e-8; 3], 10.),
            config,
            || false,
        )
        .unwrap();
        let release = &result.ejection.release;
        let reaction = result.reaction.as_ref().unwrap();
        // Loaded reduced inertia in (carrier speed, third-slider rate) is
        // [[17,19],[19,33]]. Its RHS retains the ORIGINAL admitted gear defect.
        let v = 23. / 25. - 7. * defect / 40.;
        let u = -13. / 50. + defect / 40.;
        near(release.state_after.base.twist_body.linear.x, v, 3e-10);
        for (i, scale) in [2., 2., 1.].into_iter().enumerate() {
            near(release.state_after.joint_velocity[i], scale * u, 3e-10);
        }
        near(
            release.payloads_after[0].impulse_on_payload_body.force.x,
            23. / 5. - defect / 4.,
            3e-10,
        );
        near(
            release.payloads_after[1].impulse_on_payload_body.force.x,
            -31. / 25. + 7. * defect / 20.,
            3e-10,
        );
        near(
            reaction.loop_multipliers_ns[0],
            -16. / 5. + defect / 2.,
            3e-10,
        );
        let lambda = 11. / 10. - 7. * defect / 8.;
        near(reaction.joint_multipliers[0], lambda, 3e-10);
        near(
            release.supplied_work_j,
            188. / 25. - 11. * defect / 20.,
            4e-10,
        );
        near(release.constraint_work_j, 0.5 * defect * lambda, 4e-10);
        assert_relative(&result.ejection, &inputs, 4e-10);
        for &residual in &reaction.joint_velocity_residuals_after {
            near(residual, 0., 2e-10);
        }
        for &residual in &reaction.velocity_residual_after_mps {
            near(residual, 0., 2e-10);
        }
        check_balance(&loaded, release, 2e-9);

        // The gear-only branch has a different native response preparation.
        // Its loaded reduced matrix is [[17,4,11],[4,4,0],[11,0,17]].
        let free = eject_payloads(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &rows,
            limits(&[2e-8; 3], 10.),
            || false,
        )
        .unwrap();
        near(
            free.release.state_after.base.twist_body.linear.x,
            19. / 25. - 3. * defect / 20.,
            3e-10,
        );
        for (actual, expected) in free.release.state_after.joint_velocity.iter().zip([
            -29. / 25. + 3. * defect / 20.,
            11. / 25. - defect / 10.,
            11. / 50. - defect / 20.,
        ]) {
            near(*actual, expected, 3e-10);
        }
        near(
            free.release.payloads_after[0]
                .impulse_on_payload_body
                .force
                .x,
            3.,
            3e-10,
        );
        near(
            free.release.payloads_after[1]
                .impulse_on_payload_body
                .force
                .x,
            -23. / 25. + 3. * defect / 10.,
            3e-10,
        );
        near(
            free.release.supplied_work_j,
            124. / 25. - 3. * defect / 20.,
            4e-10,
        );
        near(
            free.release.constraint_work_j,
            0.5 * defect * (3. / 10. - 3. * defect / 4.),
            4e-10,
        );
        assert_relative(&free, &inputs, 4e-10);
        check_balance(&loaded, &free.release, 2e-9);
    }
}

#[test]
fn g3_zero_commands_preserve_passive_bits_and_cancelling_targets_preserve_the_carrier() {
    let (model, state, rows, loops) = linkage(1e-3);
    let mut loaded = load(
        &model,
        &state,
        &[
            (1, Se3::identity(), sphere(2., 2.)),
            (3, Se3::identity(), sphere(1., 1.)),
        ],
    );
    loaded.state_after.base.twist_body.angular.x = -0.;
    let inputs = [
        PayloadEjectionInput {
            attachment: &loaded.attachments[0],
            relative_twist_payload: Twist::zero(),
        },
        PayloadEjectionInput {
            attachment: &loaded.attachments[1],
            relative_twist_payload: Twist::zero(),
        },
    ];
    let passive = [
        PayloadReleaseInput {
            attachment: &loaded.attachments[0],
            impulse: None,
        },
        PayloadReleaseInput {
            attachment: &loaded.attachments[1],
            impulse: None,
        },
    ];
    let config = BilateralImpulseConfig {
        velocity_tolerance_mps: 2e-3,
        ..Default::default()
    };
    let result = eject_mixed(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &loops,
        &rows,
        limits(&[2e-8; 3], 0.),
        config,
        || false,
    )
    .unwrap();
    let original = release_mixed(
        &loaded.model_after,
        &loaded.state_after,
        &passive,
        &loops,
        &rows,
        limits(&[2e-8; 3], 0.).payload,
        config,
        || false,
    )
    .unwrap();
    assert_eq!(result.ejection.release.state_after, loaded.state_after);
    assert_eq!(
        result.ejection.release.state_after,
        original.release.state_after
    );
    for (a, b) in result
        .ejection
        .release
        .state_after
        .base
        .twist_body
        .to_array()
        .into_iter()
        .zip(loaded.state_after.base.twist_body.to_array())
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for (a, b) in result
        .ejection
        .release
        .payloads_after
        .iter()
        .zip(&original.release.payloads_after)
    {
        assert_eq!(a.payload_after, b.payload_after);
        near(a.supplied_work_j, 0., 0.);
    }
    near(result.ejection.release.supplied_work_j, 0., 0.);
    near(result.ejection.release.constraint_work_j, 0., 0.);
    assert_relative(&result.ejection, &inputs, 2e-12);
    let free = eject_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &rows,
        limits(&[2e-8; 3], 0.),
        || false,
    )
    .unwrap();
    assert_eq!(free.release.state_after, loaded.state_after);
    near(free.release.supplied_work_j, 0., 0.);

    // Opposing nonzero relative targets exactly cancel the loaded momentum
    // bias. Detached bodies move; no carrier projection or stopping is added.
    let loaded = scalar(Vec3::default());
    let inputs = targets(
        &loaded,
        [
            Twist::new(Vec3::new(0., 0., 1.), Vec3::new(2., 0., 0.)),
            Twist::new(Vec3::new(0., 0., -0.5), Vec3::new(-1., 0., 0.)),
        ],
    );
    let result = eject_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &[],
        limits(&[], 4.),
        || false,
    )
    .unwrap();
    assert_eq!(result.release.state_after, loaded.state_after);
    near(result.release.supplied_work_j, 15. / 4., 2e-12);
    assert_relative(&result, &inputs, 2e-12);
    check_balance(&loaded, &result.release, 2e-10);
}

#[test]
fn g4_count_inventory_work_and_actual_rate_refusal_are_atomic_through_final_cancellation() {
    let loaded = scalar(Vec3::default());
    let saved = loaded.clone();
    let command = [
        Twist::new(Vec3::new(0., 0., 2.), Vec3::new(3., 0., 0.)),
        Twist::new(Vec3::new(0., 0., -0.5), Vec3::new(-1., 0., 0.)),
    ];
    let inputs = targets(&loaded, command);
    let run = |inputs: &[PayloadEjectionInput<'_>], control, stop: &mut dyn FnMut() -> bool| {
        eject_payloads(
            &loaded.model_after,
            &loaded.state_after,
            inputs,
            &[],
            control,
            stop,
        )
    };
    assert!(run(&[], limits(&[], 10.), &mut || false).is_err());
    let oversized = vec![inputs[0]; MAX_PAYLOAD_EJECTION_TARGETS + 1];
    assert!(matches!(
        run(&oversized, limits(&[], 10.), &mut || false),
        Err(PayloadEjectionError::Payload(PayloadError::Invalid(_)))
    ));
    let mut wrong_order = inputs;
    wrong_order.swap(0, 1);
    assert!(run(&wrong_order, limits(&[], 10.), &mut || false).is_err());
    let mut malformed = inputs;
    malformed[1].relative_twist_payload.linear.x = f64::NAN;
    assert!(run(&malformed, limits(&[], 10.), &mut || false).is_err());
    let mut bad = limits(&[], 10.);
    bad.relative_angular_velocity_tolerance_rad_s = -1.;
    assert!(run(&inputs, bad, &mut || false).is_err());
    assert!(matches!(
        run(&inputs, limits(&[], 0.), &mut || false),
        Err(PayloadEjectionError::Payload(
            PayloadError::EnergyExchange { .. }
        ))
    ));

    let mut polls = 0;
    let expected = run(&inputs, limits(&[], 10.), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        run(&inputs, limits(&[], 10.), &mut || {
            count += 1;
            count == polls
        }),
        Err(PayloadEjectionError::Payload(PayloadError::Cancelled))
    ));
    assert_eq!(count, polls);
    assert_eq!(loaded, saved);
    let retried = run(&inputs, limits(&[], 10.), &mut || false).unwrap();
    assert_eq!(retried.release.model_after, expected.release.model_after);
    assert_eq!(retried.release.state_after, expected.release.state_after);
    assert_eq!(
        retried.release.payloads_after,
        expected.release.payloads_after
    );
    assert_eq!(retried.release.balance, expected.release.balance);
    assert_eq!(
        retried.maximum_relative_velocity_residual_m_s,
        expected.maximum_relative_velocity_residual_m_s
    );
    assert_eq!(
        retried.maximum_relative_angular_velocity_residual_rad_s,
        expected.maximum_relative_angular_velocity_residual_rad_s
    );

    // At this boost both output velocities are multiples of an ULP far
    // coarser than 0.1. Exact relative equality is unrepresentable: admission
    // must inspect the returned states, not merely the computed source.
    let boosted = scalar(Vec3::new(2_f64.powi(30), 0., 0.));
    let input = targets(
        &boosted,
        [
            Twist::new(Vec3::default(), Vec3::new(0.1, 0., 0.)),
            Twist::zero(),
        ],
    );
    let mut actual = limits(&[], 1.);
    actual.relative_velocity_tolerance_m_s = 1e-5;
    actual.payload.energy_tolerance_j = 2e-6;
    let rounded = eject_payloads(
        &boosted.model_after,
        &boosted.state_after,
        &input,
        &[],
        actual,
        || false,
    )
    .unwrap();
    assert!(rounded.maximum_relative_velocity_residual_m_s > 0.);
    actual.relative_velocity_tolerance_m_s = 0.;
    assert!(matches!(
        eject_payloads(
            &boosted.model_after,
            &boosted.state_after,
            &input,
            &[],
            actual,
            || false
        ),
        Err(PayloadEjectionError::RelativeVelocity { .. })
    ));
}

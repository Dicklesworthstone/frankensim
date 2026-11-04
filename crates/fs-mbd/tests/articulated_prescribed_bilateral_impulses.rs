//! G0/G1/G3/G4: native support momentum and coupled material/joint impulse events.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionError;
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::MaterialPoint;
use fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialError;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::prescribed::apply_impulses;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop, BilateralLoopError, DistanceLoop, PointCoincidenceLoop,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link,
    PrescribedBaseImpulseResponse, SpatialInertia, apply_prescribed_base_impulses,
};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * (1. + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    for (a, b) in [actual.x, actual.y, actual.z]
        .into_iter()
        .zip([expected.x, expected.y, expected.z])
    {
        near(a, b);
    }
}
fn wrench(actual: Wrench, expected: Wrench) {
    vector(actual.force, expected.force);
    vector(actual.torque, expected.torque);
}
fn accounts(receipt: &PrescribedBaseImpulseResponse) {
    vector(receipt.linear_momentum_residual_world, Vec3::default());
    vector(receipt.angular_momentum_residual_world, Vec3::default());
    near(receipt.kinetic_work_residual, 0.);
    near(
        receipt.impulse_work,
        receipt.applied_impulse_work + receipt.support_work,
    );
}
fn sphere(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [mass, 0., 0., 0., mass, 0., 0., 0., mass],
        },
    )
    .unwrap()
}
fn sliders(maximum_speed: Option<f64>) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(5.),
    )];
    for (i, mass) in [2., 3., 4.].into_iter().enumerate() {
        links.push(Link::new(
            format!("slider {i}"),
            Some(0),
            Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1., 0., 0.),
                maximum_speed.map(|v| JointLimits::new(-10., 10., v, 0.01).unwrap()),
            )
            .unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 7.,
        base: BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, -0.2, 0.1)),
            Twist::new(Vec3::new(0.2, 0.3, -0.4), Vec3::new(-1., 2., 3.)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
    }
}
fn rows() -> [JointConstraint; 1] {
    [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }]
}
fn loops() -> [BilateralLoop; 1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        length_m: 2.,
    })]
}

#[test]
fn fixed_and_driven_supports_share_the_coupled_response_and_full_work_partition() {
    let model = sliders(Some(0.55));
    let mut start = initial();
    let source_joint = [0.6, -0.2, 0.4];
    let mut source_body = vec![Wrench::default(); 4];
    source_body[0] = Wrench::new(Vec3::new(0.2, -0.1, 0.05), Vec3::new(0.3, -0.4, 0.2));
    // The source alone takes the first slider to 0.7; only the actual coupled
    // endpoint at 0.5 is subject to its 0.55 speed limit. Impulses are not forces.
    assert!(
        apply_prescribed_base_impulses(
            &model,
            start.base,
            &start.joint_position,
            &start.joint_velocity,
            &source_joint,
            &source_body,
        )
        .is_err()
    );
    let support = Wrench::new(Vec3::new(-0.2, 0.1, -0.05), Vec3::new(0.4, 0.4, -0.2));
    let driven = apply_impulses(
        &model,
        &start,
        &source_joint,
        &source_body,
        &loops(),
        &rows(),
        Default::default(),
        || false,
    )
    .unwrap();
    // Independent joint mass diag(2,3,4), rows [-1,1,0] and [0,1,-2]:
    // W=[[5/6,1/3],[1/3,4/3]]. Its coupled solution is (0.4,0.1).
    for (actual, expected) in driven
        .impulse
        .joint_velocity_after
        .iter()
        .zip([0.5, 0.5, 0.25])
    {
        near(*actual, expected);
    }
    near(driven.loop_multipliers_ns[0], 0.4);
    near(driven.joint_multipliers[0], 0.1);
    near(driven.minimum_scaled_pivot, 0.9);
    let correlation = 1. / 10.0_f64.sqrt();
    near(
        driven.reciprocal_condition,
        (1. - correlation) / (1. + correlation),
    );
    wrench(driven.impulse.support_impulse_body, support);
    wrench(driven.impulse.external_impulse_world, source_body[0]);
    near(driven.source_work_j, 0.58);
    near(driven.constraint_work_j, 0.);
    near(driven.row_constraint_work_j, 0.);
    near(driven.work_partition_residual_j, 0.);
    near(driven.impulse.support_work, 0.18);
    near(driven.impulse.kinetic_energy_change, 0.76);
    assert_eq!(driven.impulse.base_after, start.base);
    accounts(&driven.impulse);

    let pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
        Vec3::new(3., -2., 5.),
    )
    .unwrap();
    start.base.world_from_base = pose;
    let transformed = apply_impulses(
        &model,
        &start,
        &source_joint,
        &source_body,
        &loops(),
        &rows(),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        transformed.impulse.joint_velocity_after,
        driven.impulse.joint_velocity_after
    );
    wrench(transformed.impulse.support_impulse_body, support);
    let force = pose.rotation().rotate(support.force).unwrap();
    wrench(
        transformed.impulse.support_impulse_world,
        Wrench::new(
            pose.rotation().rotate(support.torque).unwrap() + pose.translation().cross(force),
            force,
        ),
    );
    near(transformed.source_work_j, driven.source_work_j);
    near(
        transformed.impulse.support_work,
        driven.impulse.support_work,
    );
    accounts(&transformed.impulse);

    start.base = BaseState::stationary(Se3::identity());
    let fixed = apply_impulses(
        &model,
        &start,
        &source_joint,
        &source_body,
        &loops(),
        &rows(),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        fixed.impulse.joint_velocity_after,
        driven.impulse.joint_velocity_after
    );
    wrench(fixed.impulse.support_impulse_body, support);
    near(fixed.source_work_j, 0.27);
    near(fixed.impulse.support_work, 0.);
    near(fixed.impulse.kinetic_energy_change, 0.27);
    accounts(&fixed.impulse);
}

#[test]
fn mobile_root_joint_uses_native_serial_inertia_and_empty_constraints_delegate_exactly() {
    let model = ArticulatedModel::new(
        [2., 3., 4.]
            .into_iter()
            .enumerate()
            .map(|(i, mass)| {
                Link::new(
                    format!("serial {i}"),
                    i.checked_sub(1),
                    Se3::identity(),
                    JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
                    sphere(mass),
                )
            })
            .collect(),
    )
    .unwrap();
    let start = PrescribedBaseMotionState {
        time_s: 0.,
        base: BaseState::stationary(Se3::identity()),
        joint_position: vec![1., 1., 1.],
        joint_velocity: vec![0.2, 0.2, -0.2],
    };
    let closure = [BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        length_m: 2.,
    })];
    let gear = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![1., -1., 0.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }];
    let source = [0.9, -0.2, 0.3];
    let body = [Wrench::default(); 3];
    let result = apply_impulses(
        &model,
        &start,
        &source,
        &body,
        &closure,
        &gear,
        Default::default(),
        || false,
    )
    .unwrap();
    // M=[[9,7,4],[7,7,4],[4,4,4]], tangent t=(1,1,-1), t'Mt=18.
    // t'J=0.4 gives du=1/45, with a nonzero root actuator support reaction.
    for (actual, expected) in
        result
            .impulse
            .joint_velocity_after
            .iter()
            .zip([2. / 9., 2. / 9., -2. / 9.])
    {
        near(*actual, expected);
    }
    near(result.loop_multipliers_ns[0], -19. / 90.);
    near(result.joint_multipliers[0], -19. / 30.);
    wrench(
        result.impulse.support_impulse_body,
        Wrench::new(Vec3::default(), Vec3::new(4. / 15., 0., 0.)),
    );
    near(result.source_work_j, 19. / 225.);
    near(result.constraint_work_j, 0.);
    near(result.impulse.support_work, 0.);
    accounts(&result.impulse);

    let native = apply_prescribed_base_impulses(
        &model,
        start.base,
        &start.joint_position,
        &start.joint_velocity,
        &source,
        &body,
    )
    .unwrap();
    let empty = apply_impulses(
        &model,
        &start,
        &source,
        &body,
        &[],
        &[],
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(empty.impulse, native);
    assert_eq!(empty.source_work_j, native.applied_impulse_work);
    assert_eq!(empty.constraint_work_j, 0.);
    assert_eq!(empty.minimum_scaled_pivot, 1.);
    assert_eq!(empty.reciprocal_condition, 1.);
    assert!(empty.rows.is_empty());
    assert!(empty.joint_multipliers.is_empty());
}

#[test]
fn full_tensor_spherical_closure_retains_geometric_couples_and_all_support_components() {
    let masses = [5., 1., 1., 1., 2.];
    let centers = [0., 1., 2., 3., 4.].map(|i| Vec3::new(0.03 * i, -0.02 * i, 0.01 * i));
    let axes = [
        Vec3::new(0., 0., 1.),
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
        Vec3::new(0., 0., 1.),
    ];
    let model = ArticulatedModel::new(
        (0..5)
            .map(|i| {
                let mass = masses[i];
                let joint = if i == 0 || i == 4 {
                    JointModel::revolute(axes[i], None).unwrap()
                } else {
                    JointModel::prismatic(axes[i], None).unwrap()
                };
                Link::new(
                    format!("tensor {i}"),
                    i.checked_sub(1),
                    Se3::identity(),
                    joint,
                    SpatialInertia::new(
                        mass,
                        centers[i],
                        Mat3 {
                            m: [
                                mass,
                                0.1 * mass,
                                -0.05 * mass,
                                0.1 * mass,
                                1.2 * mass,
                                0.08 * mass,
                                -0.05 * mass,
                                0.08 * mass,
                                1.5 * mass,
                            ],
                        },
                    )
                    .unwrap(),
                )
            })
            .collect(),
    )
    .unwrap();
    let gap = Vec3::new(2e-4, -1e-4, 1.5e-4);
    let origins = [
        Vec3::default(),
        Vec3::new(gap.x, 0., 0.),
        Vec3::new(gap.x, gap.y, 0.),
        gap,
        gap,
    ];
    let base_twist = Twist::new(Vec3::new(0.1, -0.2, 0.15), Vec3::new(0.3, -0.4, 0.1));
    let start = PrescribedBaseMotionState {
        time_s: 0.,
        base: BaseState::prescribed(Se3::identity(), base_twist, Twist::zero()),
        joint_position: vec![0., gap.x, gap.y, gap.z, 0.],
        joint_velocity: vec![0.25, 0., 0., 0., 0.],
    };
    let closure = [BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 3,
            point_body: Vec3::default(),
        },
    })];
    let lock = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![0., 0., 0., 0., 1.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }];
    let controls = BilateralImpulseConfig {
        position_tolerance_m: 1e-3,
        velocity_tolerance_mps: 1e-3,
        ..Default::default()
    };
    let force = Vec3::new(0.7, -0.4, 0.3);
    let torque = Vec3::new(0.2, -0.1, 0.25).cross(force) + Vec3::new(0.08, 0.03, -0.06);
    let mut source_body = [Wrench::default(); 5];
    source_body[4] = Wrench::new(torque, force);
    let result = apply_impulses(
        &model,
        &start,
        &[0.; 5],
        &source_body,
        &closure,
        &lock,
        controls,
        || false,
    )
    .unwrap();
    // With the relative translations and rotor locked, the sole mode is root
    // rotation. Independently sum the original COM tensors and parallel axes.
    let total_torque = torque + gap.cross(force);
    let iz: f64 = (0..5)
        .map(|i| {
            let c = origins[i] + centers[i];
            masses[i] * (1.5 + c.x * c.x + c.y * c.y)
        })
        .sum();
    let du = total_torque.z / iz;
    near(result.impulse.joint_velocity_after[0], 0.25 + du);
    for rate in &result.impulse.joint_velocity_after[1..] {
        near(*rate, 0.);
    }
    let mut dp = Vec3::default();
    let mut dl = Vec3::default();
    for i in 0..5 {
        let c = origins[i] + centers[i];
        let linear = Vec3::new(0., 0., du).cross(c).scale(masses[i]);
        dp = dp + linear;
        dl = dl + Vec3::new(-0.05, 0.08, 1.5).scale(masses[i] * du) + c.cross(linear);
    }
    let support = Wrench::new(dl - total_torque, dp - force);
    wrench(result.impulse.support_impulse_body, support);
    vector(result.impulse.linear_momentum_change_world, dp);
    vector(result.impulse.angular_momentum_change_world, dl);
    let omega = base_twist.angular + Vec3::new(0., 0., 0.25);
    let source_work =
        total_torque.dot(omega) + force.dot(base_twist.linear) + 0.5 * total_torque.z * du;
    near(result.source_work_j, source_work);
    near(result.impulse.support_work, support.pairing(base_twist));
    near(result.constraint_work_j, 0.);
    near(result.row_constraint_work_j, 0.);
    near(result.work_partition_residual_j, 0.);
    near(
        result.impulse.kinetic_energy_change,
        source_work + support.pairing(base_twist),
    );
    vector(result.point_residual_before[0].separation_a_m, gap);
    vector(
        result.point_residual_before[0].relative_velocity_world_mps,
        omega.cross(gap),
    );
    vector(
        result.point_residual_after[0].relative_velocity_world_mps,
        (omega + Vec3::new(0., 0., du)).cross(gap),
    );
    accounts(&result.impulse);

    let stricter_speed = BilateralImpulseConfig {
        velocity_tolerance_mps: 1e-6,
        ..controls
    };
    assert!(matches!(
        apply_impulses(
            &model,
            &start,
            &[0.; 5],
            &source_body,
            &closure,
            &lock,
            stricter_speed,
            || false,
        ),
        Err(BilateralLoopError::PointVelocity {
            stage: "initial",
            ..
        })
    ));
}

#[test]
fn admitted_slip_has_measured_reaction_work_and_zero_source_stays_passive() {
    let model = sliders(None);
    let mut start = initial();
    start.base = BaseState::stationary(Se3::identity());
    let slip = 1e-3;
    start.joint_velocity[2] -= slip / 2.;
    let mut gear = rows();
    gear[0].velocity_tolerance = 2e-3;
    let source = [0.6, -0.2, 0.4];
    let body = [Wrench::default(); 4];
    let result = apply_impulses(
        &model,
        &start,
        &source,
        &body,
        &loops(),
        &gear,
        Default::default(),
        || false,
    )
    .unwrap();
    let mu = 0.4 + slip / 3.;
    let lambda = 0.1 - 5. * slip / 6.;
    near(result.loop_multipliers_ns[0], mu);
    near(result.joint_multipliers[0], lambda);
    for (actual, expected) in result.impulse.joint_velocity_after.iter().zip([
        0.5 - slip / 6.,
        0.5 - slip / 6.,
        0.25 - slip / 12.,
    ]) {
        near(*actual, expected);
    }
    near(result.source_work_j, 0.27 - 3. * slip / 20.);
    near(result.constraint_work_j, lambda * slip / 2.);
    near(result.row_constraint_work_j, lambda * slip / 2.);
    near(
        result.impulse.support_impulse_body.force.x,
        0.7 + 5. * slip / 6.,
    );
    assert!(result.constraint_work_j > 0.);
    near(result.work_partition_residual_j, 0.);
    accounts(&result.impulse);

    let passive = apply_impulses(
        &model,
        &start,
        &[0.; 3],
        &body,
        &loops(),
        &gear,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(passive.impulse.joint_velocity_after, start.joint_velocity);
    assert_eq!(passive.impulse.base_after, start.base);
    assert_eq!(passive.loop_multipliers_ns, vec![0.]);
    assert_eq!(passive.joint_multipliers, vec![0.]);
    assert_eq!(passive.constraint_body_impulses, vec![Wrench::default(); 4]);
    assert_eq!(passive.constraint_joint_impulses, vec![0.; 3]);
    near(passive.joint_velocity_residuals_before[0], slip);
    near(passive.joint_velocity_residuals_after[0], slip);
    assert_eq!(passive.source_work_j, 0.);
    assert_eq!(passive.constraint_work_j, 0.);
    accounts(&passive.impulse);
}

#[test]
fn combined_rank_limits_and_cancellation_refuse_atomically() {
    let model = sliders(Some(0.55));
    let start = initial();
    let snapshot = start.clone();
    let source = [0.6, -0.2, 0.4];
    let body = [Wrench::default(); 4];
    let duplicate = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![-1., 1., 0.],
        position: 2.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }];
    // Zero source does not waive admission of a redundant material/joint system.
    assert!(matches!(
        apply_impulses(
            &model,
            &start,
            &[0.; 3],
            &body,
            &loops(),
            &duplicate,
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Mobility { .. })
    ));
    assert!(matches!(
        apply_impulses(
            &model,
            &start,
            &[6., -2., 4.],
            &body,
            &loops(),
            &rows(),
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::JointVelocityLimitViolation { .. })
        )))
    ));
    assert!(matches!(
        apply_impulses(
            &model,
            &start,
            &source,
            &body,
            &vec![loops()[0]; 33],
            &[],
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Invalid(_))
    ));
    let mut bad_time = start.clone();
    bad_time.time_s = f64::NAN;
    assert!(matches!(
        apply_impulses(
            &model,
            &bad_time,
            &source,
            &body,
            &[],
            &[],
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Invalid(_))
    ));

    let polls = Cell::new(0);
    let expected = apply_impulses(
        &model,
        &start,
        &source,
        &body,
        &loops(),
        &rows(),
        Default::default(),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    assert!(polls.get() > 10);
    for stop in [1, polls.get() / 2, polls.get()] {
        let calls = Cell::new(0);
        assert!(matches!(
            apply_impulses(
                &model,
                &start,
                &source,
                &body,
                &loops(),
                &rows(),
                Default::default(),
                || {
                    calls.set(calls.get() + 1);
                    calls.get() == stop
                },
            ),
            Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::Cancelled
            )))
        ));
        assert_eq!(start, snapshot);
    }
    let retry = apply_impulses(
        &model,
        &start,
        &source,
        &body,
        &loops(),
        &rows(),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(start, snapshot);
}

//! G0/G1/G3/G4: joint gears and material closures share one physical mobility.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::MaterialFrame;
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, DistanceLoop, MixedBilateralLoopStep,
    PointCoincidenceLoop, joint, step_free_floating_bilateral_loops,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialError, step_free_floating_spatial_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::{
    MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-9 * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a - b).norm() <= 2e-8, "{a:?} != {b:?}");
}
fn inertia(m: f64, spin: f64) -> SpatialInertia {
    SpatialInertia::new(
        m,
        Vec3::default(),
        Mat3 {
            m: [spin, 0., 0., 0., spin, 0., 0., 0., spin],
        },
    )
    .unwrap()
}
fn point(link: usize) -> MaterialPoint {
    MaterialPoint {
        link,
        point_body: Vec3::default(),
    }
}
fn loads<'a>(u: &'a [f64], w: &'a [Wrench]) -> SpatialMotionLoads<'a> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: u,
            external_body_wrench: w,
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    }
}
fn row(a: Vec<f64>, coordinate: ConstraintCoordinate) -> JointConstraint {
    JointConstraint {
        coefficients: a,
        coordinate,
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }
}
fn distance() -> [BilateralLoop; 1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: point(1),
        b: point(2),
        length_m: 2.,
    })]
}
fn sliders(limited: bool) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(5., 4.),
    )];
    for (i, mass) in [2., 3., 4.].into_iter().enumerate() {
        links.push(Link::new(
            format!("slider {i}"),
            Some(0),
            Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1., 0., 0.),
                limited.then(|| JointLimits::new(-10., 10., 10., 0.1).unwrap()),
            )
            .unwrap(),
            inertia(mass, mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
    }
}
fn account(s: &MixedBilateralLoopStep, joint_row_scale: f64) {
    let b = &s.bilateral;
    vector(b.linear_momentum_residual_world, Vec3::default());
    vector(b.angular_momentum_residual_world, Vec3::default());
    vector(b.reaction_linear_impulse_world, Vec3::default());
    vector(b.reaction_angular_impulse_world, Vec3::default());
    close(
        b.mechanical_work_residual,
        b.constraint_impulse_work
            + b.midpoint.mechanical_work_residual
            + b.initial_reaction.kinetic_work_residual
            + b.final_reaction.kinetic_work_residual,
    );
    // Joint receipts retain the caller's original row units and scaling.
    for r in &s.joint_position_residuals_after {
        close(*r / joint_row_scale, 0.);
    }
    for r in &s.joint_velocity_residuals_after {
        close(*r / joint_row_scale, 0.);
    }
}

#[test]
fn geared_sliders_share_cross_mobility_carrier_recoil_and_actual_event_work() {
    let model = sliders(false);
    let start = initial();
    let held = [3., -1., 2.];
    let mut w = [Wrench::default(); 4];
    w[0].force.x = 2.1;
    let rows = [row(vec![0., 1., -2.], ConstraintCoordinate::Length)];
    let h = 0.2;
    let step = joint::step_with_bushings(
        &model,
        &start,
        loads(&held, &w),
        &distance(),
        &rows,
        &[],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&step, 1.0);
    // Independent elimination: qdot=t*udot, t=(2,2,1).
    // Mjoint=diag(2,3,4)-m*m^T/14; t^T Mjoint t=10.
    // Carrier acceleration is F/14-uddot, including the finite root mass.
    let a = (2. * held[0] + 2. * held[1] + held[2] - w[0].force.x) / 10.;
    let carrier_a = w[0].force.x / 14. - a;
    let vm = 0.2 + 0.5 * h * a;
    for (i, t) in [2., 2., 1.].into_iter().enumerate() {
        close(
            step.bilateral.state_after.joint_position[i],
            start.joint_position[i] + h * t * vm,
        );
        close(
            step.bilateral.state_after.joint_velocity[i],
            t * (0.2 + h * a),
        );
    }
    close(
        step.bilateral.state_after.base.twist_body.linear.x,
        0.7 + h * carrier_a,
    );
    close(
        step.bilateral
            .state_after
            .base
            .world_from_base
            .translation()
            .x,
        0.7 * h + 0.5 * h * h * carrier_a,
    );
    let mu = 0.5 * h * (held[0] - w[0].force.x / 7. - 2. * a);
    let lambda = 0.25 * h * (held[2] - 2. * w[0].force.x / 7.);
    close(step.bilateral.initial_impulse_ns[0], mu);
    close(step.bilateral.final_impulse_ns[0], mu);
    close(step.initial_joint_multipliers[0], lambda);
    close(step.final_joint_multipliers[0], lambda);
    // W=[[5/6,1/3],[1/3,23/15]]: its cross block is essential.
    let joint_work = 0.5 * lambda * (mu / 3. + 23. * lambda / 15.);
    let total_work =
        0.5 * (5. * mu * mu / 6. + 2. * mu * lambda / 3. + 23. * lambda * lambda / 15.);
    close(step.initial_joint_reaction_work_j, joint_work);
    close(step.final_joint_reaction_work_j, -joint_work);
    close(step.bilateral.initial_reaction.impulse_work, total_work);
    close(step.bilateral.final_reaction.impulse_work, -total_work);
    assert!(joint_work > 0.005);
    close(step.bilateral.constraint_impulse_work, 0.);
    close(step.bilateral.mechanical_work_residual, 0.);
    close(
        step.bilateral.midpoint.motion.actuator_work,
        h * (2. * held[0] + 2. * held[1] + held[2]) * vm,
    );
    close(
        step.bilateral.midpoint.motion.external_work,
        w[0].force.x * (0.7 * h + 0.5 * h * h * carrier_a),
    );

    // Rigid world transformation and a common inertial boost preserve the
    // relative trajectory, internal event work and actual carrier recoil.
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let translation = Vec3::new(1.2, -0.7, 0.9);
    let boost = Vec3::new(-0.5, 0.3, 0.2);
    let mut moved = start.clone();
    moved.base.world_from_base = Se3::from_parts(rotation, translation).unwrap();
    moved.base.twist_body.linear =
        start.base.twist_body.linear + rotation.inverse().rotate(boost).unwrap();
    let transformed = joint::step_with_bushings(
        &model,
        &moved,
        loads(&held, &w),
        &distance(),
        &rows,
        &[],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&transformed, 1.0);
    for i in 0..3 {
        close(
            transformed.bilateral.state_after.joint_position[i],
            step.bilateral.state_after.joint_position[i],
        );
        close(
            transformed.bilateral.state_after.joint_velocity[i],
            step.bilateral.state_after.joint_velocity[i],
        );
    }
    vector(
        transformed
            .bilateral
            .state_after
            .base
            .world_from_base
            .translation(),
        translation
            + rotation
                .rotate(
                    step.bilateral
                        .state_after
                        .base
                        .world_from_base
                        .translation(),
                )
                .unwrap()
            + boost.scale(h),
    );
    vector(
        transformed.bilateral.state_after.base.twist_body.linear,
        step.bilateral.state_after.base.twist_body.linear
            + rotation.inverse().rotate(boost).unwrap(),
    );
    close(transformed.initial_joint_reaction_work_j, joint_work);

    // The same physical constraints may contain no body-loop rows at all.
    let joint_only = [
        JointConstraint {
            position: 2.,
            ..row(vec![-1., 1., 0.], ConstraintCoordinate::Length)
        },
        rows[0].clone(),
    ];
    let pure = joint::step_with_bushings(
        &model,
        &start,
        loads(&held, &w),
        &[],
        &joint_only,
        &[],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&pure, 1.0);
    assert!(pure.bilateral.rows.is_empty());
    assert!(pure.bilateral.initial_impulse_ns.is_empty());
    for i in 0..3 {
        close(
            pure.bilateral.state_after.joint_position[i],
            step.bilateral.state_after.joint_position[i],
        );
        close(
            pure.bilateral.state_after.joint_velocity[i],
            step.bilateral.state_after.joint_velocity[i],
        );
    }
    vector(
        pure.bilateral
            .state_after
            .base
            .world_from_base
            .translation(),
        step.bilateral
            .state_after
            .base
            .world_from_base
            .translation(),
    );
    vector(
        pure.bilateral.state_after.base.twist_body.linear,
        step.bilateral.state_after.base.twist_body.linear,
    );
    close(pure.initial_joint_multipliers[0], mu);
    close(pure.initial_joint_multipliers[1], lambda);
    close(pure.final_joint_multipliers[0], mu);
    close(pure.final_joint_multipliers[1], lambda);
    close(pure.initial_joint_reaction_work_j, total_work);
    close(
        pure.bilateral.constraint_impulse_work,
        step.bilateral.constraint_impulse_work,
    );
    assert_eq!(start, initial());
}

fn coaxial() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(5., 4.),
        ),
        Link::new(
            "first slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2., 2.),
        ),
        Link::new(
            "second slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(3., 3.),
        ),
        Link::new(
            "rotor",
            Some(2),
            Se3::identity(),
            JointModel::revolute(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(4., 4.),
        ),
    ])
    .unwrap()
}
fn full_bushing() -> CoupledSpatialBushing {
    let prepare = |a: [[f64; 2]; 2]| {
        let mut m = [[0.; 6]; 6];
        for i in 0..2 {
            for j in 0..2 {
                m[3 * i][3 * j] = a[i][j];
            }
        }
        PassiveBushingMatrix::try_new(
            m,
            BushingCoordinateScales {
                translation_m: 0.4,
                rotation_rad: 0.2,
            },
        )
        .unwrap()
    };
    CoupledSpatialBushing::try_new(
        MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        MaterialFrame {
            link: 3,
            link_from_attachment: Se3::identity(),
        },
        Vec3::default(),
        prepare([[12., 4.], [4., 8.]]),
        prepare([[0.6, 0.1], [0.1, 0.4]]),
    )
    .unwrap()
}

#[test]
fn full_matrices_feedback_and_length_angle_gears_share_one_reduced_mode() {
    let model = coaxial();
    let q0 = 0.15;
    let v0 = -0.1;
    let start = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(
                Vec3::new(-4. * v0 / 13., 0., 0.),
                Vec3::new(-9. * v0 / 7., 0., 0.),
            ),
        ),
        joint_position: vec![2. * q0 - 2., 2. * q0, q0],
        joint_velocity: vec![2. * v0, 2. * v0, v0],
    };
    let held = [0.3, -0.2, 0.5];
    let w = [Wrench::default(); 4];
    let h = 0.6;
    // Translation and spin decouple before the gear:
    // t^T M_translation t=90/7, J_rotor=36/13.
    // [x,theta]=[2,1]*q gives K=72,C=3.2. Rotor feedback adds 2,.3.
    let mass = 90. / 7. + 36. / 13.;
    let vm = (v0 + h * (0.7 - 74. * q0) / (2. * mass))
        / (1. + h * 3.5 / (2. * mass) + h * h * 74. / (4. * mass));
    let qm = q0 + 0.5 * h * vm;
    let q1 = q0 + h * vm;
    let v1 = 2. * vm - v0;
    let config = BilateralLoopConfig {
        nonlinear: fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialConfig {
            angular_velocity_scale: 0.4,
            linear_velocity_scale: 2.,
            residual_tolerance: 1e-12,
            ..Default::default()
        },
        ..Default::default()
    };
    // Include a sign change, metre/radian reinterpretation, and extreme row
    // scaling with matching absolute tolerances. Physical impulses stay equal.
    let mut reference = None;
    for (factor, coordinate) in [
        (1.0_f64, ConstraintCoordinate::Length),
        (-0.5, ConstraintCoordinate::Angle),
        (1e6, ConstraintCoordinate::Length),
        (1e-6, ConstraintCoordinate::Length),
    ] {
        let rows = [JointConstraint {
            position_tolerance: 1e-9 * factor.abs(),
            velocity_tolerance: 1e-9 * factor.abs(),
            ..row(vec![0., factor, -2. * factor], coordinate)
        }];
        let mut last = None;
        let step = joint::step_with_coupled_bushings_and_effort(
            &model,
            &start,
            loads(&held, &w),
            &distance(),
            &rows,
            &[full_bushing()],
            h,
            config,
            |elapsed, stage, out| {
                close(elapsed, h / 2.);
                assert_eq!(out, &held);
                out[2] -= 2. * stage.joint_position[2] + 0.3 * stage.joint_velocity[2];
                last = Some(stage.clone());
                Ok(())
            },
            || false,
        )
        .unwrap();
        account(&step, factor);
        let last = last.unwrap();
        close(last.joint_position[2], qm);
        close(last.joint_velocity[2], vm);
        close(last.base.twist_body.angular.x, -4. * vm / 13.);
        close(last.base.twist_body.linear.x, -9. * vm / 7.);
        let s = &step.bilateral;
        close(s.state_after.joint_position[0], 2. * q1 - 2.);
        close(s.state_after.joint_position[1], 2. * q1);
        close(s.state_after.joint_position[2], q1);
        close(s.state_after.joint_velocity[0], 2. * v1);
        close(s.state_after.joint_velocity[1], 2. * v1);
        close(s.state_after.joint_velocity[2], v1);
        close(s.state_after.base.twist_body.angular.x, -4. * v1 / 13.);
        close(s.state_after.base.twist_body.linear.x, -9. * v1 / 7.);
        close(
            s.final_reaction.momentum_after.kinetic_energy,
            0.5 * mass * v1 * v1,
        );
        close(s.midpoint.stored_energy_before, 36. * q0 * q0);
        close(s.midpoint.stored_energy_after, 36. * q1 * q1);
        close(s.midpoint.dissipated_energy, h * 3.2 * vm * vm);
        close(
            s.midpoint.motion.actuator_work,
            h * (0.7 - 2. * qm - 0.3 * vm) * vm,
        );
        close(s.mechanical_work_residual, 0.);
        close(s.constraint_impulse_work, 0.);
        if let Some(value) = reference {
            close(factor * step.initial_joint_multipliers[0], value);
        } else {
            reference = Some(step.initial_joint_multipliers[0]);
        }
    }
}

#[test]
fn admitted_tangent_error_has_real_signed_joint_work_without_double_counting() {
    let model = sliders(false);
    let mut start = initial();
    let defect = 1e-3;
    start.joint_velocity[2] -= defect / 2.;
    let rows = [JointConstraint {
        velocity_tolerance: 2e-3,
        ..row(vec![0., 1., -2.], ConstraintCoordinate::Length)
    }];
    let w = [Wrench::default(); 4];
    let step = joint::step_with_bushings(
        &model,
        &start,
        loads(&[0.; 3], &w),
        &distance(),
        &rows,
        &[],
        0.1,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&step, 1.0);
    // Eliminate the distance row from W: W_AA-W_AB^2/W_BB=7/5.
    // The initial finite slip is removed by a genuine combined impulse.
    close(step.joint_velocity_residuals_before[0], defect);
    close(step.initial_joint_multipliers[0], -5. * defect / 7.);
    close(step.bilateral.initial_impulse_ns[0], 2. * defect / 7.);
    close(
        step.initial_joint_reaction_work_j,
        -5. * defect * defect / 14.,
    );
    close(
        step.bilateral.constraint_impulse_work,
        -5. * defect * defect / 14.,
    );
    close(
        step.bilateral.mechanical_work_residual,
        -5. * defect * defect / 14.,
    );
    close(step.bilateral.midpoint.dissipated_energy, 0.);
    close(step.final_joint_reaction_work_j, 0.);
    close(step.final_joint_multipliers[0], 0.);
    let strict = [JointConstraint {
        velocity_tolerance: defect / 2.,
        ..rows[0].clone()
    }];
    assert!(matches!(
        joint::step_with_bushings(
            &model,
            &start,
            loads(&[0.; 3], &w),
            &distance(),
            &strict,
            &[],
            0.1,
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Motion(
            ImplicitSpatialError::Constraint(ConstraintError::Velocity { .. })
        ))
    ));
}

fn combined_inertia(parts: &[SpatialInertia]) -> SpatialInertia {
    let mass: f64 = parts.iter().copied().map(SpatialInertia::mass).sum();
    let center = parts.iter().fold(Vec3::default(), |c, p| {
        c + p.center_of_mass().scale(p.mass() / mass)
    });
    let mut matrix = [0.; 9];
    for p in parts {
        let d = p.center_of_mass() - center;
        let x = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in 0..3 {
                matrix[3 * i + j] += p.inertia_com().m[3 * i + j]
                    + p.mass() * ((if i == j { d.dot(d) } else { 0. }) - x[i] * x[j]);
            }
        }
    }
    SpatialInertia::new(mass, center, Mat3 { m: matrix }).unwrap()
}

#[test]
fn spherical_closure_and_angular_lock_match_independent_full_tensor_assembly() {
    let parts: Vec<_> = [5., 1., 1., 1., 2.]
        .into_iter()
        .enumerate()
        .map(|(i, m)| {
            let f = i as f64;
            SpatialInertia::new(
                m,
                Vec3::new(0.03 * f, -0.02 * f, 0.01 * f),
                Mat3 {
                    m: [
                        m,
                        0.1 * m,
                        -0.05 * m,
                        0.1 * m,
                        1.2 * m,
                        0.08 * m,
                        -0.05 * m,
                        0.08 * m,
                        1.5 * m,
                    ],
                },
            )
            .unwrap()
        })
        .collect();
    let axes = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ];
    let mut links = vec![Link::new(
        "base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        parts[0],
    )];
    for i in 0..3 {
        links.push(Link::new(
            format!("slide {i}"),
            Some(i),
            Se3::identity(),
            JointModel::prismatic(axes[i], None).unwrap(),
            parts[i + 1],
        ));
    }
    links.push(Link::new(
        "rotor",
        Some(3),
        Se3::identity(),
        JointModel::revolute(axes[2], None).unwrap(),
        parts[4],
    ));
    let model = ArticulatedModel::new(links).unwrap();
    let reduced = ArticulatedModel::new(vec![Link::new(
        "assembled body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        combined_inertia(&parts),
    )])
    .unwrap();
    let base = FreeFloatingBaseState::new(
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.1))).unwrap(),
            Vec3::new(0.5, -0.3, 0.2),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.2, -0.15, 0.1), Vec3::new(0.3, -0.2, 0.4)),
    );
    let start = FreeFloatingMotionState {
        base,
        joint_position: vec![0.; 4],
        joint_velocity: vec![0.; 4],
    };
    let reference = FreeFloatingMotionState {
        base,
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let world = [WorldPointLoad {
        point: point(4),
        force_world: Vec3::new(0.3, -0.2, 0.1),
        couple_world: Vec3::new(0.05, 0.03, -0.04),
    }];
    let reduced_world = [WorldPointLoad {
        point: point(0),
        ..world[0]
    }];
    let w = [Wrench::default(); 5];
    let wr = [Wrench::default(); 1];
    let loops = [BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a: point(0),
        b: point(3),
    })];
    let rows = [row(vec![0., 0., 0., 1.], ConstraintCoordinate::Angle)];
    let h = 0.025;
    let actual = joint::step_with_bushings(
        &model,
        &start,
        SpatialMotionLoads {
            world_points: &world,
            ..loads(&[0.; 4], &w)
        },
        &loops,
        &rows,
        &[],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    let expected = step_free_floating_spatial_implicit(
        &reduced,
        &reference,
        SpatialMotionLoads {
            world_points: &reduced_world,
            ..loads(&[], &wr)
        },
        h,
        Default::default(),
        || false,
    )
    .unwrap()
    .spatial;
    account(&actual, 1.0);
    assert_eq!(actual.bilateral.rows.len(), 3);
    assert_eq!(actual.initial_joint_multipliers.len(), 1);
    for q in &actual.bilateral.state_after.joint_position {
        close(*q, 0.);
    }
    for v in &actual.bilateral.state_after.joint_velocity {
        close(*v, 0.);
    }
    let a = actual.bilateral.state_after.base;
    let b = expected.motion.state_after.base;
    vector(
        a.world_from_base.translation(),
        b.world_from_base.translation(),
    );
    vector(a.twist_body.linear, b.twist_body.linear);
    vector(a.twist_body.angular, b.twist_body.angular);
    for (x, y) in a
        .world_from_base
        .rotation()
        .matrix()
        .m
        .iter()
        .zip(b.world_from_base.rotation().matrix().m)
    {
        close(*x, y);
    }
    close(
        actual
            .bilateral
            .final_reaction
            .momentum_after
            .kinetic_energy,
        expected.motion.momentum_after.kinetic_energy,
    );
    close(
        actual.bilateral.midpoint.motion.external_work,
        expected.motion.external_work,
    );
    assert!(actual.initial_joint_multipliers[0].abs() > 1e-7);
}

#[test]
fn empty_joint_rows_are_exact_and_cross_family_redundancy_and_cancellation_refuse() {
    let model = sliders(false);
    let start = initial();
    let held = [0.2, 0.1, -0.3];
    let w = [Wrench::default(); 4];
    let load = loads(&held, &w);
    let loops = distance();
    let native = step_free_floating_bilateral_loops(
        &model,
        &start,
        load,
        &loops,
        0.1,
        Default::default(),
        || false,
    )
    .unwrap();
    let empty = joint::step_with_bushings(
        &model,
        &start,
        load,
        &loops,
        &[],
        &[],
        0.1,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(empty.bilateral, native);
    assert!(empty.initial_joint_multipliers.is_empty());
    assert_eq!(empty.initial_joint_impulse, vec![0.; 3]);
    let duplicate = [JointConstraint {
        position: 2.,
        ..row(vec![-1., 1., 0.], ConstraintCoordinate::Length)
    }];
    assert!(matches!(
        joint::step_with_bushings(
            &model,
            &start,
            load,
            &loops,
            &duplicate,
            &[],
            0.1,
            Default::default(),
            || false,
        ),
        Err(BilateralLoopError::Mobility { .. })
    ));
    let rows = [row(vec![0., 1., -2.], ConstraintCoordinate::Length)];
    // Exercise cancellation before, during and after shared row admission.
    // Existing body-only APIs retain their original cancellation variant.
    for threshold in 1..=32 {
        let mut polls = 0;
        let error = joint::step_with_bushings(
            &model,
            &start,
            load,
            &loops,
            &rows,
            &[],
            0.1,
            Default::default(),
            || {
                polls += 1;
                polls >= threshold
            },
        )
        .unwrap_err();
        assert!(
            matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Cancelled
                ))
            ),
            "{error:?}"
        );
        let mut polls = 0;
        let error = step_free_floating_bilateral_loops(
            &model,
            &start,
            load,
            &loops,
            0.1,
            Default::default(),
            || {
                polls += 1;
                polls >= threshold
            },
        )
        .unwrap_err();
        assert!(
            matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Cancelled
                ))
            ),
            "{error:?}"
        );
    }
    assert_eq!(start, initial());
}

#[test]
fn passive_joint_impulses_do_not_consume_motor_effort_limits() {
    let model = sliders(true);
    let start = initial();
    let mut w = [Wrench::default(); 4];
    w[0].force.x = 14.;
    let rows = [row(vec![0., 1., -2.], ConstraintCoordinate::Length)];
    let h = 0.05;
    let step = joint::step_with_bushings(
        &model,
        &start,
        loads(&[0.; 3], &w),
        &distance(),
        &rows,
        &[],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&step, 1.0);
    // A passive force of 4 N on joint 3 is transmitted through a 0.1 N motor.
    close(step.initial_joint_impulse[2], 0.5 * h * 4.);
    close(step.bilateral.midpoint.motion.actuator_work, 0.);
    let error = joint::step_with_bushings_and_effort(
        &model,
        &start,
        loads(&[0.; 3], &w),
        &distance(),
        &rows,
        &[],
        h,
        Default::default(),
        |_, _, effort| {
            effort[2] = 0.2;
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::JointEffortLimitViolation { .. })
        ))
    ));
}

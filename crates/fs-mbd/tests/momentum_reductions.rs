//! G0/G3: small physical transfers amid large, cancelling momenta and work.
use fs_ga::{Mat3, Se3, Twist, Vec3 as GaVec, Wrench as GaWrench};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    apply_free_floating_impulses, apply_prescribed_base_impulses, free_floating_momentum,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, RigidKineticEnergyChange, UnitQuaternion, Vec3,
    Wrench, apply_equal_and_opposite_impulse_at_body_points,
};

const LARGE: f64 = (1u64 << 55) as f64;

#[test]
fn native_impulse_storage_exposes_nondyadic_velocity_rounding() {
    let p = 2.0_f64.powi(54);
    let property = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap();
    let before = state(Vec3::ZERO, Vec3::new(p, 0.0, 0.0));
    let impulse = Vec3::new(4.0, 0.0, 0.0);
    let resultant = before.apply_momentum_impulse(property, impulse, Vec3::ZERO).unwrap();
    let point = before.apply_impulse_at_body_point(property, Vec3::ZERO, impulse).unwrap();
    // Exact storage is ((2^54+4)^2-(2^54)^2)/6 = 24019198012642648.
    // Native COM velocities pair to 24019198012642646; their rounded work
    // equals rounded storage, but the independent 2 J defect must survive.
    for receipt in [resultant.work, point.work] {
        assert_eq!(receipt.kinetic_energy_change_j, 24_019_198_012_642_648.0);
        assert_eq!(receipt.kinetic_energy_change_j, receipt.impulse_work_j);
        assert_eq!(receipt.work_energy_residual_j, 2.0);
    }
    let spin = RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(p, 0.0, 0.0)).unwrap();
    let result = spin.apply_momentum_impulse(property, Vec3::ZERO, Vec3::new(8.0, 0.0, 0.0)).unwrap();
    // Exact rotational change is 144115188075855904/3; native angular
    // velocities give work 48038396025285300, leaving a 4/3 J defect.
    assert_eq!(result.work.kinetic_energy_change_j, 48_038_396_025_285_304.0);
    assert!((result.work.work_energy_residual_j - 4.0 / 3.0).abs() < 2e-14);
}

#[test]
fn public_body_reduction_and_action_reaction_share_the_same_physical_storage() {
    let p = 2.0_f64.powi(54);
    let before = state(Vec3::ZERO, Vec3::new(p, 0.0, 0.0));
    for mass in [2.0, 3.0, 7.0] {
        let property = MassProperties::new(mass, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap();
        let pair = apply_equal_and_opposite_impulse_at_body_points(
            before, property, Vec3::ZERO, before, property, Vec3::ZERO,
            Vec3::new(4.0, 0.0, 0.0),
        ).unwrap();
        let mut change = RigidKineticEnergyChange::new();
        change.add_body(property, before, pair.first.state_after).unwrap();
        change.add_body(property, before, pair.second.state_after).unwrap();
        let expected = 16.0 / mass;
        assert!((pair.kinetic_energy_change_j - expected).abs() <= 8.0 * f64::EPSILON * expected);
        assert_eq!(change.finish().unwrap(), pair.kinetic_energy_change_j);
        assert!((pair.kinetic_energy_change_j - pair.impulse_work_j - pair.work_energy_residual_j).abs()
            <= 16.0 * f64::EPSILON * expected);
    }
    assert_eq!(RigidKineticEnergyChange::new().finish().unwrap(), 0.0);
}

fn properties() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}

#[test]
fn gravity_and_applied_force_keep_their_actual_residual_momentum() {
    // Exact integer oracle: (n+1)*(n-1) - n*n = -1, across binary unit scales.
    let n = 1_i128 << 27;
    let mass = (n + 1) as f64;
    let property =
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let initial = state(Vec3::ZERO, Vec3::ZERO);
    for g_scale in [2.0_f64.powi(-24), 1.0, 2.0_f64.powi(20)] {
        for axis in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ] {
            let gravity = Gravity::new(axis.scale((n - 1) as f64 * g_scale)).unwrap();
            let engine = RigidBodyIntegrator::new(gravity);
            let wrench = Wrench {
                force_world: axis.scale(-(n * n) as f64 * g_scale),
                torque_body: Vec3::ZERO,
            };
            let dt = 2.0;
            let expected_momentum = axis.scale(-g_scale * dt);
            let expected_position = axis.scale(-g_scale * dt * dt / (2.0 * mass));
            for receipt in [
                engine.step(initial, property, wrench, dt).unwrap(),
                engine.step_cayley(initial, property, wrench, dt).unwrap(),
                engine.step_symplectic(initial, property, wrench, dt).unwrap(),
                engine
                    .begin_verlet(initial, property, wrench, dt)
                    .unwrap()
                    .finish(wrench)
                    .unwrap(),
            ] {
                assert_eq!(receipt.state_after.linear_momentum_world(), expected_momentum);
                assert_eq!(receipt.state_after.pose().position_world(), expected_position);
                let energy = g_scale * g_scale * dt * dt / (2.0 * mass);
                assert!(
                    (receipt.work.kinetic_energy_change_j - energy).abs()
                        <= 4.0 * f64::EPSILON * energy
                );
                assert!(receipt.work.work_energy_residual_j.abs() <= 16.0 * f64::EPSILON * energy);
            }
            let kicked = engine
                .step_kick_drift(initial, property, wrench, dt)
                .unwrap();
            assert_eq!(kicked.state_after.linear_momentum_world(), expected_momentum);
            assert_eq!(
                kicked.state_after.pose().position_world(),
                expected_position.scale(2.0)
            );
            let motion = engine
                .prepare_point_motion(initial, property, wrench)
                .unwrap()
                .point_motion(Vec3::ZERO)
                .unwrap();
            assert_eq!(
                motion.center_of_mass_acceleration_world,
                axis.scale(-g_scale / mass)
            );

            // Different endpoint forces must retain their own gravitational
            // residuals; a held-load shortcut cannot replace their trapezoid.
            let last = Wrench {
                force_world: axis.scale((-(n * n) + 4) as f64 * g_scale),
                torque_body: Vec3::ZERO,
            };
            let varied = engine
                .begin_verlet(initial, property, wrench, dt)
                .unwrap()
                .finish(last)
                .unwrap();
            assert_eq!(
                varied.state_after.linear_momentum_world(),
                axis.scale(g_scale * dt)
            );
            assert_eq!(varied.state_after.pose().position_world(), expected_position);
        }
    }
}

#[test]
fn instantaneous_acceleration_does_not_require_a_representable_weight() {
    let property =
        MassProperties::new(1e308, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let engine = RigidBodyIntegrator::new(Gravity::new(Vec3::new(2.0, 0.0, 0.0)).unwrap());
    let motion = engine
        .prepare_point_motion(
            state(Vec3::ZERO, Vec3::ZERO),
            property,
            Wrench {
                force_world: Vec3::new(1e308, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
        )
        .unwrap()
        .point_motion(Vec3::ZERO)
        .unwrap();
    assert_eq!(motion.center_of_mass_acceleration_world, Vec3::new(3.0, 0.0, 0.0));
}

fn state(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

fn inertia() -> SpatialInertia {
    SpatialInertia::new(1.0, GaVec::default(), Mat3::identity()).unwrap()
}

#[test]
fn cross_products_preserve_centrality_and_antisymmetry() {
    // Binary scaling keeps the supplied vectors exactly parallel. A single
    // FMA against an already rounded second product leaves a false torque.
    for scale in [2.0_f64.powi(-100), 1.0, 2.0_f64.powi(100)] {
        let arm = Vec3::new(0.1, -0.3, 0.7).scale(scale);
        assert_eq!(arm.cross(arm), Vec3::ZERO);
        assert_eq!(arm.cross(arm.scale(32.0)), Vec3::ZERO);
        let other = Vec3::new(-0.2, 0.4, 0.9).scale(scale);
        assert_eq!(arm.cross(other), other.cross(arm).scale(-1.0));
    }
}

#[test]
fn native_central_impulses_and_axial_motion_do_not_invent_spin() {
    let arm = Vec3::new(0.1, -0.3, 0.7);
    let impulse = arm.scale(32.0);
    let before = state(arm, Vec3::ZERO);
    let receipt = before
        .apply_impulse_at_body_point(properties(), arm, impulse)
        .unwrap();
    assert_eq!(receipt.angular_impulse_body, Vec3::ZERO);
    assert_eq!(receipt.state_after.angular_momentum_body(), Vec3::ZERO);
    assert_eq!(receipt.state_after.linear_momentum_world(), impulse);
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    assert_eq!(
        integrator
            .diagnostics(receipt.state_after, properties())
            .unwrap()
            .angular_momentum_world,
        Vec3::ZERO
    );

    // Spherical inertia has no gyroscopic acceleration, and a material point
    // on its rotation axis has neither rotational velocity nor centripetal
    // acceleration. Exercise the actual native derivative/kinematics owners.
    let spinning = RigidBodyState::new(Pose::identity(), Vec3::ZERO, arm).unwrap();
    let point = integrator
        .prepare_point_motion(spinning, properties(), Wrench::ZERO)
        .unwrap()
        .point_motion(arm.scale(64.0))
        .unwrap();
    assert_eq!(point.angular_acceleration_body, Vec3::ZERO);
    assert_eq!(point.kinematics.point_velocity_world, Vec3::ZERO);
    assert_eq!(point.point_acceleration_world, Vec3::ZERO);
}

#[test]
fn nearly_central_impulses_retain_small_exact_angular_momentum() {
    // For integer n = 2^27, n*n - (n+1)*(n-1) = 1 exactly. Both
    // independently rounded products are 2^54, so an uncompensated second
    // product loses the complete physical torque. i128 is an independent
    // exact oracle, including sign and all three axis permutations.
    let n = 1_i128 << 27;
    let arms = [[n, n + 1, 0], [0, n, n + 1], [n + 1, 0, n]];
    let impulses = [[n - 1, n, 0], [0, n - 1, n], [n, 0, n - 1]];
    for (r, j) in arms.into_iter().zip(impulses) {
        let exact = [
            r[1] * j[2] - r[2] * j[1],
            r[2] * j[0] - r[0] * j[2],
            r[0] * j[1] - r[1] * j[0],
        ];
        let arm = Vec3::new(r[0] as f64, r[1] as f64, r[2] as f64);
        let impulse = Vec3::new(j[0] as f64, j[1] as f64, j[2] as f64);
        for sign in [-1.0, 1.0] {
            let expected = Vec3::new(exact[0] as f64, exact[1] as f64, exact[2] as f64).scale(sign);
            let receipt = state(Vec3::ZERO, Vec3::ZERO)
                .apply_impulse_at_body_point(properties(), arm, impulse.scale(sign))
                .unwrap();
            assert_eq!(receipt.angular_impulse_body, expected);
            assert_eq!(receipt.state_after.angular_momentum_body(), expected);
        }
    }
}

#[test]
fn resultant_and_point_work_retain_the_half_ulp_midpoint() {
    let old = state(Vec3::ZERO, Vec3::new(LARGE, 0.0, LARGE));
    let impulse = Vec3::new(8.0, 1.0, -8.0);
    // All momentum increments are representable. The exact energy increment
    // is p dot J + |J|²/2 = (8L - 8L) + (64 + 1 + 64)/2.
    let result = old
        .apply_momentum_impulse(properties(), impulse, Vec3::ZERO)
        .unwrap();
    assert_eq!(result.linear_increment_world, impulse);
    assert_eq!(result.work.kinetic_energy_change_j, 64.5);
    assert_eq!(result.work.impulse_work_j, 64.5);
    assert_eq!(result.work.work_energy_residual_j, 0.0);
    let point = old
        .apply_impulse_at_body_point(properties(), Vec3::ZERO, impulse)
        .unwrap();
    assert_eq!(point.work.kinetic_energy_change_j, 64.5);
    assert_eq!(point.work.impulse_work_j, 64.5);
    assert_eq!(point.work.work_energy_residual_j, 0.0);
}

#[test]
fn action_reaction_work_is_reduced_before_rounding_individual_body_work() {
    let old = state(Vec3::ZERO, Vec3::new(LARGE, 0.0, 0.0));
    let result = apply_equal_and_opposite_impulse_at_body_points(
        old,
        properties(),
        Vec3::ZERO,
        old,
        properties(),
        Vec3::ZERO,
        Vec3::new(8.0, 0.0, 0.0),
    )
    .unwrap();
    // Each body gains J²/2, while the two +/- L*J terms cancel.
    assert_eq!(result.kinetic_energy_change_j, 64.0);
    assert_eq!(result.impulse_work_j, 64.0);
    assert_eq!(result.work_energy_residual_j, 0.0);
    assert_eq!(result.linear_momentum_change_world, Vec3::ZERO);
    assert_eq!(result.angular_momentum_change_world, Vec3::ZERO);
}

#[test]
fn rounded_pair_total_cannot_hide_a_lost_momentum_transfer() {
    let first = state(Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0));
    let second = state(Vec3::ZERO, Vec3::new(LARGE, 0.0, 0.0));
    let result = apply_equal_and_opposite_impulse_at_body_points(
        first,
        properties(),
        Vec3::ZERO,
        second,
        properties(),
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    assert_eq!(result.impulse_balance_world, Vec3::ZERO);
    assert_eq!(result.first.state_after.linear_momentum_world(), Vec3::ZERO);
    assert_eq!(result.second.state_after, second);
    // Summing the large absolute endpoints first incorrectly reports zero.
    assert_eq!(
        result.linear_momentum_change_world,
        Vec3::new(-1.0, 0.0, 0.0)
    );
}

#[test]
fn pair_angular_balance_includes_local_spin_and_orbital_changes() {
    for translation in [0.0, 1e14] {
        let first = state(Vec3::new(translation - 1.0, 0.0, 0.0), Vec3::ZERO);
        let second = state(Vec3::new(translation + 1.0, 0.0, 0.0), Vec3::ZERO);
        let common = apply_equal_and_opposite_impulse_at_body_points(
            first,
            properties(),
            Vec3::new(1.0, 0.0, 0.0),
            second,
            properties(),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
        )
        .unwrap();
        assert_eq!(common.angular_impulse_world, Vec3::ZERO);
        assert_eq!(common.angular_momentum_change_world, Vec3::ZERO);
        assert_eq!(common.angular_momentum_residual_world, Vec3::ZERO);
        let separated = apply_equal_and_opposite_impulse_at_body_points(
            first,
            properties(),
            Vec3::ZERO,
            second,
            properties(),
            Vec3::ZERO,
            Vec3::new(0.0, 2.0, 0.0),
        )
        .unwrap();
        assert_eq!(separated.angular_impulse_world, Vec3::new(0.0, 0.0, 4.0));
        assert_eq!(
            separated.angular_momentum_change_world,
            separated.angular_impulse_world
        );
        assert_eq!(separated.angular_momentum_residual_world, Vec3::ZERO);
    }
}

#[test]
fn articulated_countermotion_retains_small_linear_and_angular_momentum() {
    for rates in [
        [LARGE, 1.0, -LARGE],
        [LARGE, -LARGE, 1.0],
        [1.0, LARGE, -LARGE],
    ] {
        for rotary in [false, true] {
            let mut links = vec![Link::new(
                "root",
                None,
                Se3::identity(),
                JointModel::FIXED,
                inertia(),
            )];
            for i in 0..3 {
                let joint = if rotary {
                    JointModel::revolute(GaVec::new(0.0, 0.0, 1.0), None).unwrap()
                } else {
                    JointModel::prismatic(GaVec::new(1.0, 0.0, 0.0), None).unwrap()
                };
                links.push(Link::new(
                    format!("child{i}"),
                    Some(0),
                    Se3::identity(),
                    joint,
                    inertia(),
                ));
            }
            let model = ArticulatedModel::new(links).unwrap();
            let result = free_floating_momentum(
                &model,
                FreeFloatingBaseState::new(Se3::identity(), Twist::zero()),
                &[0.0; 3],
                &rates,
            )
            .unwrap();
            if rotary {
                assert_eq!(result.angular_momentum_world, GaVec::new(0.0, 0.0, 1.0));
            } else {
                assert_eq!(result.linear_momentum_world, GaVec::new(1.0, 0.0, 0.0));
            }
        }
    }
}

#[test]
fn free_articulated_impulse_retains_small_signed_work_on_a_large_translation() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(),
    )])
    .unwrap();
    let result = apply_free_floating_impulses(
        &model,
        FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(GaVec::default(), GaVec::new(LARGE, 0.0, LARGE)),
        ),
        &[],
        &[],
        &[],
        &[GaWrench::new(GaVec::default(), GaVec::new(8.0, 1.0, -8.0))],
    )
    .unwrap();
    assert_eq!(
        result.linear_momentum_change_world,
        GaVec::new(8.0, 1.0, -8.0)
    );
    assert_eq!(result.kinetic_energy_change, 64.5);
    assert_eq!(result.impulse_work, 64.5);
    assert_eq!(result.kinetic_work_residual, 0.0);
}

#[test]
fn prescribed_joints_and_body_loads_retain_the_same_small_impulse_work() {
    let axes = [
        GaVec::new(1.0, 0.0, 0.0),
        GaVec::new(0.0, 1.0, 0.0),
        GaVec::new(0.0, 0.0, 1.0),
    ];
    let mut links = vec![Link::new(
        "root",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(),
    )];
    for (index, axis) in axes.iter().enumerate() {
        links.push(Link::new(
            format!("slider{index}"),
            Some(0),
            Se3::identity(),
            JointModel::prismatic(*axis, None).unwrap(),
            inertia(),
        ));
    }
    let model = ArticulatedModel::new(links).unwrap();
    let jumps = [8.0, 1.0, -8.0];
    for body_loads in [false, true] {
        let mut loads = [GaWrench::default(); 4];
        if body_loads {
            for i in 0..3 {
                loads[i + 1] = GaWrench::new(GaVec::default(), axes[i].scale(jumps[i]));
            }
        }
        let result = apply_prescribed_base_impulses(
            &model,
            BaseState::prescribed(Se3::identity(), Twist::zero(), Twist::zero()),
            &[0.0; 3],
            &[LARGE, 0.0, LARGE],
            if body_loads { &[0.0; 3] } else { &jumps },
            &loads,
        )
        .unwrap();
        assert_eq!(result.kinetic_energy_change, 64.5);
        assert_eq!(result.applied_impulse_work, 64.5);
        assert_eq!(result.support_work, 0.0);
        assert_eq!(result.kinetic_work_residual, 0.0);
        assert_eq!(result.linear_momentum_residual_world, GaVec::default());
    }
}

#[test]
fn prescribed_support_balances_a_small_couple_at_a_distant_world_origin() {
    let model = ArticulatedModel::new(vec![
        Link::new("root", None, Se3::identity(), JointModel::FIXED, inertia()),
        Link::new(
            "child",
            Some(0),
            Se3::identity(),
            JointModel::FIXED,
            inertia(),
        ),
    ])
    .unwrap();
    let base = BaseState::prescribed(
        Se3::from_parts(fs_ga::So3::identity(), GaVec::new(LARGE, 0.0, 0.0)).unwrap(),
        Twist::zero(),
        Twist::zero(),
    );
    let result = apply_prescribed_base_impulses(
        &model,
        base,
        &[],
        &[],
        &[],
        &[
            GaWrench::new(GaVec::new(0.0, 0.0, 1.0), GaVec::new(0.0, 1.0, 0.0)),
            GaWrench::new(GaVec::default(), GaVec::new(0.0, -1.0, 0.0)),
        ],
    )
    .unwrap();
    assert_eq!(result.external_impulse_world.force, GaVec::default());
    assert_eq!(
        result.external_impulse_world.torque,
        GaVec::new(0.0, 0.0, 1.0)
    );
    assert_eq!(result.support_impulse_world.force, GaVec::default());
    assert_eq!(
        result.support_impulse_world.torque,
        GaVec::new(0.0, 0.0, -1.0)
    );
    assert_eq!(result.angular_momentum_residual_world, GaVec::default());
    assert_eq!(result.kinetic_energy_change, 0.0);
    assert_eq!(result.impulse_work, 0.0);
}

#[test]
fn large_impulses_report_absorbed_existing_linear_and_spin_momentum() {
    for sign in [-1.0, 1.0] {
        let old = RigidBodyState::new(
            Pose::identity(),
            Vec3::new(sign, 0.0, 0.0),
            Vec3::new(0.0, sign, 0.0),
        )
        .unwrap();
        let linear = Vec3::new(sign * 1e20, 0.0, 0.0);
        let angular = Vec3::new(0.0, sign * 1e20, 0.0);
        let result = old
            .apply_momentum_impulse(properties(), linear, angular)
            .unwrap();
        assert_eq!(result.state_after.linear_momentum_world(), linear);
        assert_eq!(result.state_after.angular_momentum_body(), angular);
        // Each f64 endpoint difference rounds back to the requested impulse.
        // The original one-unit momentum was nevertheless lost in each frame.
        assert_eq!(result.linear_increment_world, linear);
        assert_eq!(result.angular_increment_body, angular);
        assert_eq!(result.linear_residual_world, Vec3::new(-sign, 0.0, 0.0));
        assert_eq!(result.angular_residual_body, Vec3::new(0.0, -sign, 0.0));
        assert_eq!(result.work.work_energy_residual_j, -1e20);

        let point = state(Vec3::ZERO, Vec3::new(sign, 0.0, 0.0))
            .apply_impulse_at_body_point(properties(), Vec3::ZERO, linear)
            .unwrap();
        // Rounded total work and kinetic change coincide. Their physical
        // defect must come from the compensated momentum error, not subtraction.
        assert_eq!(
            point.work.kinetic_energy_change_j,
            point.work.impulse_work_j
        );
        assert_eq!(point.work.work_energy_residual_j, -5e19);
    }
}

#[test]
fn momentum_exchange_retains_the_small_original_kinetic_energy() {
    for large in [LARGE, 1e20] {
        for spin in [false, true] {
            let before = Vec3::new(1.0, 0.0, large);
            let applied = Vec3::new(large, 0.0, -large);
            let old = RigidBodyState::new(
                Pose::identity(),
                if spin { Vec3::ZERO } else { before },
                if spin { before } else { Vec3::ZERO },
            )
            .unwrap();
            let result = old
                .apply_momentum_impulse(
                    properties(),
                    if spin { Vec3::ZERO } else { applied },
                    if spin { applied } else { Vec3::ZERO },
                )
                .unwrap();
            // Unit mass/inertia: final squared momentum is large^2, while
            // initial squared momentum is large^2+1. No large K subtraction.
            assert_eq!(result.work.kinetic_energy_change_j, -0.5);
            assert_eq!(result.work.impulse_work_j, 0.5 * large);
            assert_eq!(result.work.work_energy_residual_j, -0.5 * large);
        }
    }
}

#[test]
fn equal_and_opposite_large_impulses_expose_swallowed_prior_pair_momentum() {
    for large in [LARGE, 1e20] {
        let position = Vec3::new(0.0, 2.0, 0.0);
        let first = state(position, Vec3::new(large, 0.0, 0.0));
        let second = state(position, Vec3::new(1.0, 0.0, 0.0));
        let result = apply_equal_and_opposite_impulse_at_body_points(
            first,
            properties(),
            Vec3::ZERO,
            second,
            properties(),
            Vec3::ZERO,
            Vec3::new(large, 0.0, 0.0),
        )
        .unwrap();
        assert_eq!(result.impulse_balance_world, Vec3::ZERO);
        assert_eq!(result.first.state_after.linear_momentum_world(), Vec3::ZERO);
        assert_eq!(
            result.second.state_after.linear_momentum_world(),
            Vec3::new(large, 0.0, 0.0)
        );
        assert_eq!(
            result.linear_momentum_change_world,
            Vec3::new(-1.0, 0.0, 0.0)
        );
        assert_eq!(result.angular_impulse_world, Vec3::ZERO);
        assert_eq!(
            result.angular_momentum_change_world,
            Vec3::new(0.0, 0.0, 2.0)
        );
        assert_eq!(
            result.angular_momentum_residual_world,
            Vec3::new(0.0, 0.0, 2.0)
        );
        assert_eq!(result.kinetic_energy_change_j, -0.5);
        assert_eq!(result.impulse_work_j, 0.5 * large);
        assert_eq!(result.work_energy_residual_j, -0.5 * large);
    }
}

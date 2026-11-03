//! G0/G3: small physical transfers amid large, cancelling momenta and work.
use fs_ga::{Mat3, Se3, Twist, Vec3 as GaVec, Wrench as GaWrench};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    apply_free_floating_impulses, apply_prescribed_base_impulses, free_floating_momentum,
};
use fs_mbd::{
    MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3,
    apply_equal_and_opposite_impulse_at_body_points,
};

const LARGE: f64 = (1u64 << 55) as f64;

fn properties() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
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

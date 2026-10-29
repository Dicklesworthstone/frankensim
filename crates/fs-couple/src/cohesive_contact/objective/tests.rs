use super::*;
use fs_material::cohesive::{BilinearCohesive, CohesiveViscosity, mixed_mode::MixedModeParameters};

fn law(viscous: bool) -> BilinearCohesive {
    let law = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap();
    if viscous {
        law.with_viscosity(CohesiveViscosity::new(2e5, 3e5, 5.0).unwrap())
    } else {
        law
    }
}

fn config() -> ObjectiveCohesiveConfig {
    ObjectiveCohesiveConfig {
        maximum_cayley_increment: 8.0,
        kinematic_tolerance_m: 2e-14,
        material_energy_tolerance_j_m2: 1e-8,
        spatial_work_tolerance_j: 1e-10,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    for (x, y) in array(a).into_iter().zip(array(b)) {
        close(x, y, tolerance);
    }
}

#[test]
fn stressed_viscous_bond_has_no_fictitious_damage_or_heat_under_rigid_transport() {
    let law = law(true);
    let original = law
        .bonded_point()
        .trial_over_time([3e-4, 2.5e-4, -1e-4], 0.02, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    assert!(original.damage() > 0.0 && original.damage() < 1.0);
    let mut point = ObjectiveCohesivePoint::new(
        original,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, -0.001, 0.002),
        0.002,
    )
    .unwrap();
    for angle in [2.1, -2.7, core::f64::consts::PI, 0.6, -1.8] {
        let frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), angle).unwrap();
        point = point.transport_frame(frame).unwrap();
        assert_eq!(point.material(), original);
        assert!(core::ptr::eq(point.material().law(), &law));
        assert_eq!(
            point.material().mixed_mode_history(),
            original.mixed_mode_history()
        );
        assert_eq!(
            point.material().stored_energy_j_m2().to_bits(),
            original.stored_energy_j_m2().to_bits()
        );
        close_vec(
            point.separation_world_m(),
            frame.rotate_body_to_world(point.reference_offset_m().add(vector(original.jump_m()))),
            1e-16,
        );
    }
    // An actual finite Cayley work step with declared zero physical jump rate
    // also keeps stressed dashpots and irreversible history unchanged.
    let frame = point.frame().right_exp(Vec3::new(0.5, -0.7, 0.4)).unwrap();
    let trial = point
        .trial_local_jump(frame, original.jump_m(), 0.01, config())
        .unwrap();
    assert_eq!(trial.point.material(), original);
    assert_eq!(trial.material.viscous_dissipated_j_m2, 0.0);
    assert_eq!(trial.material.equilibrium.energy.work, 0.0);
    assert_eq!(trial.material.equilibrium.energy.dissipated, 0.0);
    close(trial.spatial_work_j, 0.0, 2e-15);
    assert!(trial.discrete_wrench.frame_couple_n_m.norm_squared() > 1e-4);
}

#[test]
fn rotating_mixed_mode_rupture_and_reclosure_conserve_actual_force_torque_and_work() {
    for viscous in [false, true] {
        let law = law(viscous);
        let original = law
            .bonded_point()
            .trial_over_time([5e-5, 2e-5, -3e-5], 0.02, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        let mut point = ObjectiveCohesivePoint::new(
            original,
            UnitQuaternion::IDENTITY,
            Vec3::new(0.003, -0.002, 0.004),
            0.002,
        )
        .unwrap();
        let initial_storage = point.area_m2() * original.stored_energy_j_m2();
        let (mut work, mut fracture, mut heat) = (0.0, 0.0, 0.0);
        for (i, jump) in [
            [5e-4, 7e-4, -3e-4],
            [0.003, 0.006, -0.002],
            [-2e-4, 0.001, 0.0],
        ]
        .into_iter()
        .enumerate()
        {
            let frame = point.frame().right_exp(Vec3::new(0.3, -0.5, 0.4)).unwrap();
            let step = point.trial_local_jump(frame, jump, 0.02, config()).unwrap();
            let direct = point
                .material()
                .trial_over_time_with_mean_traction(jump, 0.02, 1e-8)
                .unwrap();
            assert_eq!(step.material, direct.0);
            assert_eq!(step.mean_traction_pa, direct.1);
            let wrench = step.discrete_wrench;
            close_vec(
                wrench.force_on_a_n.add(wrench.force_on_b_n),
                Vec3::ZERO,
                1e-16,
            );
            // Actual world attachment locations, with arbitrary common motion.
            let a0 = Vec3::new(0.6, -0.8, 1.1);
            let a1 = a0.add(Vec3::new(0.02, -0.01, 0.03));
            let b0 = a0.add(point.separation_world_m());
            let b1 = a1.add(step.point.separation_world_m());
            let amid = a0.scale(0.5).add(a1.scale(0.5));
            let bmid = b0.scale(0.5).add(b1.scale(0.5));
            let total_moment = amid
                .cross(wrench.force_on_a_n)
                .add(bmid.cross(wrench.force_on_b_n))
                .add(wrench.frame_couple_n_m);
            close_vec(total_moment, Vec3::ZERO, 3e-13);
            let body_work = wrench.force_on_a_n.dot(a1.sub(a0))
                + wrench.force_on_b_n.dot(b1.sub(b0))
                + wrench.frame_couple_n_m.dot(step.cayley_increment_world);
            close(
                body_work,
                -point.area_m2() * step.material.total_work_j_m2(),
                3e-13,
            );
            work += step.spatial_work_j;
            fracture += point.area_m2() * step.material.equilibrium.energy.dissipated;
            heat += point.area_m2() * step.material.viscous_dissipated_j_m2;
            point = step.point;
            if i == 0 {
                assert!(point.material().damage() > original.damage());
            }
            if i >= 1 {
                assert!(point.material().failed());
            }
        }
        assert!(fracture > 0.0);
        assert_eq!(heat > 0.0, viscous);
        close(
            work,
            point.area_m2() * point.material().stored_energy_j_m2() - initial_storage
                + fracture
                + heat,
            1e-12,
        );
        assert!(point.material().traction_pa()[0] < 0.0);
        assert_eq!(point.material().traction_pa()[1], 0.0);
    }
}

#[test]
fn observer_rotation_covaries_both_work_efforts_and_keeps_material_path_identical() {
    let law = law(true);
    let original = law
        .bonded_point()
        .trial_over_time([3e-4, 2e-4, -1e-4], 0.02, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    let point = ObjectiveCohesivePoint::new(
        original,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, -0.002, 0.004),
        0.002,
    )
    .unwrap();
    let observer = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -3.0, 2.0), 2.4).unwrap();
    let rotated = point.transport_frame(observer).unwrap();
    let phi = Vec3::new(0.7, -0.4, 0.2);
    let jump = [8e-5, -2.5e-4, 3e-4];
    let a = point
        .trial_local_jump(
            UnitQuaternion::IDENTITY.right_exp(phi).unwrap(),
            jump,
            0.02,
            config(),
        )
        .unwrap();
    let b = rotated
        .trial_local_jump(observer.right_exp(phi).unwrap(), jump, 0.02, config())
        .unwrap();
    assert_eq!(a.material, b.material);
    close_vec(
        b.cayley_increment_world,
        observer.rotate_body_to_world(a.cayley_increment_world),
        2e-14,
    );
    for (left, right) in [
        (a.discrete_wrench, b.discrete_wrench),
        (a.endpoint_wrench, b.endpoint_wrench),
    ] {
        close_vec(
            right.force_on_a_n,
            observer.rotate_body_to_world(left.force_on_a_n),
            2e-14,
        );
        close_vec(
            right.frame_couple_n_m,
            observer.rotate_body_to_world(left.frame_couple_n_m),
            2e-14,
        );
    }
    close(a.spatial_work_j, b.spatial_work_j, 2e-14);
    let actual = rotated
        .trial(
            b.point.frame(),
            b.point.separation_world_m(),
            0.02,
            config(),
        )
        .unwrap();
    close_vec(
        vector(actual.point.material().jump_m()),
        vector(jump),
        2e-16,
    );
    close_vec(
        actual.discrete_wrench.force_on_a_n,
        b.discrete_wrench.force_on_a_n,
        2e-12,
    );
}

#[test]
fn endpoint_force_and_balancing_couple_are_derivatives_of_the_original_elastic_energy() {
    let law = law(false);
    let point = ObjectiveCohesivePoint::new(
        law.bonded_point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, -0.002, 0.004),
        0.002,
    )
    .unwrap();
    let end = point
        .trial_local_jump(
            UnitQuaternion::IDENTITY,
            [4e-5, -3e-5, 2e-5],
            0.02,
            config(),
        )
        .unwrap();
    let d = end.point.separation_world_m();
    let energy = |r: UnitQuaternion, gap: Vec3| {
        let local = array(r.rotate_world_to_body(gap).sub(point.reference_offset_m()));
        point.area_m2()
            * point
                .material()
                .trial(local, 1e-8)
                .unwrap()
                .point
                .stored_energy_j_m2()
    };
    for (i, axis) in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .into_iter()
    .enumerate()
    {
        let h = 1e-9;
        let force = (energy(UnitQuaternion::IDENTITY, d.add(axis.scale(h)))
            - energy(UnitQuaternion::IDENTITY, d.sub(axis.scale(h))))
            / (2.0 * h);
        close(force, array(end.endpoint_wrench.force_on_a_n)[i], 2e-8);
        let h = 1e-7;
        let couple = -(energy(UnitQuaternion::from_axis_angle(axis, h).unwrap(), d)
            - energy(UnitQuaternion::from_axis_angle(axis, -h).unwrap(), d))
            / (2.0 * h);
        close(couple, array(end.endpoint_wrench.frame_couple_n_m)[i], 2e-8);
    }
    assert_ne!(end.discrete_wrench, end.endpoint_wrench);
}

#[test]
fn half_turn_range_and_invalid_physical_trials_preserve_history_and_replay() {
    let law = law(true);
    let point = ObjectiveCohesivePoint::new(
        law.bonded_point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, 0.0, 0.0),
        0.002,
    )
    .unwrap();
    let snapshot = point;
    let half_turn =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), core::f64::consts::PI).unwrap();
    assert!(matches!(
        point.trial_local_jump(half_turn, [1e-5; 3], 0.02, config()),
        Err(ObjectiveCohesiveError::RotationLimit { .. })
    ));
    assert!(
        point
            .trial_local_jump(UnitQuaternion::IDENTITY, [0.03, 0.0, 0.0], 0.02, config())
            .is_err()
    );
    assert!(
        point
            .trial_local_jump(UnitQuaternion::IDENTITY, [1e-4; 3], 1e-9, config())
            .is_err()
    );
    assert!(
        point
            .trial(
                UnitQuaternion::IDENTITY,
                Vec3::new(f64::NAN, 0.0, 0.0),
                0.02,
                config()
            )
            .is_err()
    );
    assert!(
        point
            .trial_local_jump(UnitQuaternion::IDENTITY, [0.0; 3], 0.0, config())
            .is_err()
    );
    assert_eq!(point, snapshot);
    let frame = UnitQuaternion::IDENTITY
        .right_exp(Vec3::new(0.4, -0.2, 0.3))
        .unwrap();
    assert_eq!(
        point
            .trial_local_jump(frame, [4e-4, 2e-4, 0.0], 0.02, config())
            .unwrap(),
        snapshot
            .trial_local_jump(frame, [4e-4, 2e-4, 0.0], 0.02, config())
            .unwrap()
    );
    assert_eq!(
        point.transport_frame(half_turn).unwrap().material(),
        point.material()
    );
}

#[test]
fn finite_rotation_path_subdivision_preserves_mixed_fracture_and_physical_viscous_work() {
    let law = law(true);
    let original = ObjectiveCohesivePoint::new(
        law.bonded_point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, -0.002, 0.004),
        0.002,
    )
    .unwrap();
    let end = [0.002, 0.005, -0.002];
    let phi = Vec3::new(0.3, -0.2, 0.7);
    let complete = original
        .trial_local_jump(
            UnitQuaternion::IDENTITY.right_exp(phi).unwrap(),
            end,
            0.02,
            config(),
        )
        .unwrap();
    let mut point = original;
    let (mut work, mut fracture, mut heat) = (0.0, 0.0, 0.0);
    for i in 1..=32 {
        let s = f64::from(i) / 32.0;
        let step = point
            .trial_local_jump(
                UnitQuaternion::IDENTITY.right_exp(phi.scale(s)).unwrap(),
                end.map(|v| v * s),
                0.02 / 32.0,
                config(),
            )
            .unwrap();
        work += step.spatial_work_j;
        fracture += step.material.equilibrium.energy.dissipated;
        heat += step.material.viscous_dissipated_j_m2;
        point = step.point;
    }
    assert!(point.material().failed());
    close(work, complete.spatial_work_j, 2e-12);
    close(
        fracture,
        complete.material.equilibrium.energy.dissipated,
        2e-12,
    );
    close(heat, complete.material.viscous_dissipated_j_m2, 2e-12);
    close(
        point.material().dissipated_j_m2(),
        complete.point.material().dissipated_j_m2(),
        2e-12,
    );
}

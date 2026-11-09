//! Opposing planes share one original liquid collision/stream/conversion step.

use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::{
    Cell3, ContactModel3, FreeSurface3, FreeSurfaceError3, PlanarRigidWallMotion3,
    PlanarRigidWallPairReport3, PlanarWall3,
};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 6e-12 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn tank() -> FreeSurface3 {
    FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap()
}

fn motions(axis: usize) -> [PlanarRigidWallMotion3; 2] {
    let faces = [
        [PlanarWall3::XMin, PlanarWall3::XMax],
        [PlanarWall3::YMin, PlanarWall3::YMax],
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
    ][axis];
    core::array::from_fn(|end| {
        let mut origin = [4.0; 3];
        origin[axis] = 8.0 * end as f64;
        let mut velocity_at_origin = [0.0; 3];
        velocity_at_origin[(axis + 1) % 3] = if end == 0 { 1.0 / 256.0 } else { -1.0 / 512.0 };
        velocity_at_origin[(axis + 2) % 3] = if end == 0 { -1.0 / 1024.0 } else { 1.0 / 256.0 };
        PlanarRigidWallMotion3 {
            wall: faces[end],
            origin,
            velocity_at_origin,
            angular_speed: if end == 0 {
                1.0 / 1024.0
            } else {
                -1.0 / 2048.0
            },
        }
    })
}

fn inventory(surface: &FreeSurface3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    let [nx, ny, nz] = surface.domain_cells();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let Some(sample) = surface.sample_at(x, y, z) else {
                    continue;
                };
                let mass = sample.density * sample.fill_fraction;
                let momentum = sample.velocity.map(|u| u * mass);
                let p = [x, y, z].map(|n| n as f64 + 0.5);
                let moment = cross(core::array::from_fn(|d| p[d] - origin[d]), momentum);
                for d in 0..3 {
                    linear[d] += momentum[d];
                    angular[d] += moment[d];
                }
            }
        }
    }
    (linear, angular)
}

fn check(report: PlanarRigidWallPairReport3) {
    near(report.mass_after, report.mass_before);
    near(
        report.fluid_work,
        report.walls.iter().map(|w| w.fluid_work).sum(),
    );
    assert_eq!(report.wall_work, -report.fluid_work);
    for d in 0..3 {
        near(
            report.fluid_impulse[d],
            report.walls.iter().map(|w| w.fluid_impulse[d]).sum(),
        );
        assert_eq!(report.wall_impulse[d], -report.fluid_impulse[d]);
    }
    for wall in report.walls {
        assert_eq!(wall.mass_before, report.mass_before);
        assert_eq!(wall.mass_after, report.mass_after);
        assert_eq!(wall.relative_mass_defect, report.relative_mass_defect);
    }
}

#[test]
fn g0_stationary_pair_is_the_original_single_liquid_tick() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        let mut initial = tank();
        initial
            .rebind_body_acceleration([1e-5, -2e-5, 1e-5], &gate)
            .unwrap();
        let mut ordinary = initial.clone();
        ordinary.step_with_gate(&gate).unwrap();
        let stationary = motions(axis).map(|m| PlanarRigidWallMotion3 {
            velocity_at_origin: [0.0; 3],
            angular_speed: 0.0,
            ..m
        });
        let (pair, report) = initial
            .planar_rigid_wall_pair_candidate(stationary, 1e-12, &gate)
            .unwrap();
        check(report);
        assert_eq!(pair.grid().state_bits(), ordinary.grid().state_bits());
        assert_eq!(pair.ledger_mass(), ordinary.ledger_mass());
        assert_eq!(pair.conversions(), ordinary.conversions());
        assert_eq!(report.mass_before, initial.ledger_mass());
        assert_eq!(report.fluid_work, 0.0);
    }
}

#[test]
fn g0_one_moving_plane_retains_the_original_native_population_path() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for selected in 0..2 {
            let mut initial = tank();
            let mut motion = motions(axis);
            initial
                .step_with_planar_rigid_wall_motion(motion[selected], 1e-12, &gate)
                .unwrap();
            motion[1 - selected].velocity_at_origin = [0.0; 3];
            motion[1 - selected].angular_speed = 0.0;
            let (single, old) = initial
                .planar_rigid_wall_motion_candidate(motion[selected], 1e-12, &gate)
                .unwrap();
            let (pair, new) = initial
                .planar_rigid_wall_pair_candidate(motion, 1e-12, &gate)
                .unwrap();
            check(new);
            assert_eq!(pair.grid().state_bits(), single.grid().state_bits());
            assert_eq!(pair.ledger_mass(), single.ledger_mass());
            assert_eq!(pair.conversions(), single.conversions());
            assert_eq!(new.walls[selected], old);
        }
    }
}

#[test]
fn g1_both_native_planes_match_independent_impedance_and_liquid_momenta() {
    let gate = CancelGate::new_clock_free();
    let origin = [1.0, 3.0, -2.0];
    for axis in 0..3 {
        let initial = tank();
        let motion = motions(axis);
        let (candidate, report) = initial
            .planar_rigid_wall_pair_candidate(motion, 1e-12, &gate)
            .unwrap();
        check(report);
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        // Independent sum at cell centres: rho/3 for each tangential mode.
        let mut rotation = 0.0;
        for j in 0..8 {
            for i in 0..8 {
                rotation += ((i as f64 + 0.5 - 4.0).powi(2) + (j as f64 + 0.5 - 4.0).powi(2)) / 3.0;
            }
        }
        for wall in report.walls {
            near(wall.fluid_response[0][0], 64.0 / 3.0);
            near(wall.fluid_response[1][1], 64.0 / 3.0);
            near(wall.fluid_response[2][2], rotation);
            near(
                wall.fluid_impulse[a],
                (64.0 / 3.0) * wall.motion.velocity_at_origin[a],
            );
            near(
                wall.fluid_impulse[b],
                (64.0 / 3.0) * wall.motion.velocity_at_origin[b],
            );
            near(
                wall.fluid_angular_impulse[axis],
                rotation * wall.motion.angular_speed,
            );
            near(
                wall.fluid_work,
                wall.motion.velocity_at_origin[a] * wall.fluid_impulse[a]
                    + wall.motion.velocity_at_origin[b] * wall.fluid_impulse[b]
                    + wall.motion.angular_speed * wall.fluid_angular_impulse[axis],
            );
        }
        let before = inventory(&initial, origin);
        let after = inventory(&candidate, origin);
        let angular = report.fluid_angular_impulse_at(origin).unwrap();
        assert_eq!(
            report.wall_angular_impulse_at(origin).unwrap(),
            angular.map(|x| -x)
        );
        for d in 0..3 {
            near(after.0[d] - before.0[d], report.fluid_impulse[d]);
            near(after.1[d] - before.1[d], angular[d]);
        }
    }
}

#[test]
fn g3_both_affine_responses_share_one_input_and_no_instantaneous_cross_term() {
    let gate = CancelGate::new_clock_free();
    let mut initial = tank();
    let original = motions(2);
    for _ in 0..3 {
        initial
            .step_with_planar_rigid_wall_pair(original, 1e-12, &gate)
            .unwrap();
    }
    let (_, pilot) = initial
        .planar_rigid_wall_pair_candidate(original, 1e-12, &gate)
        .unwrap();
    let changed = original.map(|m| PlanarRigidWallMotion3 {
        velocity_at_origin: m.velocity_at_origin.map(|v| -0.5 * v),
        angular_speed: -2.0 * m.angular_speed,
        ..m
    });
    let (_, replay) = initial
        .planar_rigid_wall_pair_candidate(changed, 1e-12, &gate)
        .unwrap();
    for end in 0..2 {
        let old = pilot.walls[end];
        let new = replay.walls[end];
        let delta = [
            changed[end].velocity_at_origin[0] - original[end].velocity_at_origin[0],
            changed[end].velocity_at_origin[1] - original[end].velocity_at_origin[1],
            changed[end].angular_speed - original[end].angular_speed,
        ];
        let old_load = [
            old.fluid_impulse[0],
            old.fluid_impulse[1],
            old.fluid_angular_impulse[2],
        ];
        let new_load = [
            new.fluid_impulse[0],
            new.fluid_impulse[1],
            new.fluid_angular_impulse[2],
        ];
        for d in 0..3 {
            near(
                new_load[d],
                old_load[d]
                    + (0..3)
                        .map(|j| old.fluid_response[d][j] * delta[j])
                        .sum::<f64>(),
            );
        }
        let mut only_other = original;
        only_other[1 - end] = changed[1 - end];
        let (_, independent) = initial
            .planar_rigid_wall_pair_candidate(only_other, 1e-12, &gate)
            .unwrap();
        assert_eq!(independent.walls[end].fluid_impulse, old.fluid_impulse);
        assert_eq!(
            independent.walls[end].fluid_angular_impulse,
            old.fluid_angular_impulse
        );
        assert_eq!(independent.walls[end].fluid_response, old.fluid_response);
        assert_eq!(independent.walls[end].fluid_work, old.fluid_work);
    }
}

#[test]
fn g3_caller_order_preserves_each_wrench_and_one_physical_candidate() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        let mut initial = tank();
        let motion = motions(axis);
        initial
            .step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
            .unwrap();
        let original_bits = initial.grid().state_bits();
        let (a, forward) = initial
            .planar_rigid_wall_pair_candidate(motion, 1e-12, &gate)
            .unwrap();
        let (b, reversed) = initial
            .planar_rigid_wall_pair_candidate([motion[1], motion[0]], 1e-12, &gate)
            .unwrap();
        assert_eq!(initial.grid().state_bits(), original_bits);
        assert_eq!(a.grid().state_bits(), b.grid().state_bits());
        assert_eq!(a.conversions(), b.conversions());
        assert_eq!(forward.walls, [reversed.walls[1], reversed.walls[0]]);
        assert_eq!(forward.fluid_impulse, reversed.fluid_impulse);
        assert_eq!(forward.fluid_work, reversed.fluid_work);
        let mut committed = initial.clone();
        let accepted = committed
            .step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
            .unwrap();
        assert_eq!(accepted, forward);
        assert_eq!(committed.grid().state_bits(), a.grid().state_bits());
    }
}

#[test]
fn g3_independent_origin_changes_preserve_motion_work_and_common_moment() {
    let gate = CancelGate::new_clock_free();
    let initial = tank();
    let original = motions(2);
    let shifts = [[1.0, -2.0, 3.0], [-2.0, 1.0, -1.0]];
    let changed = core::array::from_fn(|end| {
        let mut m = original[end];
        m.origin = core::array::from_fn(|d| m.origin[d] + shifts[end][d]);
        m.velocity_at_origin[0] -= m.angular_speed * shifts[end][1];
        m.velocity_at_origin[1] += m.angular_speed * shifts[end][0];
        m
    });
    let (a, before) = initial
        .planar_rigid_wall_pair_candidate(original, 1e-12, &gate)
        .unwrap();
    let (b, after) = initial
        .planar_rigid_wall_pair_candidate(changed, 1e-12, &gate)
        .unwrap();
    assert_eq!(a.grid().state_bits(), b.grid().state_bits());
    assert_eq!(before.fluid_impulse, after.fluid_impulse);
    near(before.fluid_work, after.fluid_work);
    for origin in [[0.0; 3], [1.0, 2.0, 3.0]] {
        let left = before.fluid_angular_impulse_at(origin).unwrap();
        let right = after.fluid_angular_impulse_at(origin).unwrap();
        for d in 0..3 {
            near(left[d], right[d]);
        }
    }
    assert!(
        after
            .fluid_angular_impulse_at([f64::NAN, 0.0, 0.0])
            .is_err()
    );
}

#[test]
fn g4_invalid_pair_unflooded_end_and_local_mach_refuse_without_publication() {
    let gate = CancelGate::new_clock_free();
    let mut original = tank();
    let motion = motions(2);
    let before = original.grid().state_bits();
    let mass = original.ledger_mass();
    let mut adjacent = motion;
    adjacent[1] = motions(0)[0];
    let mut fast = motion;
    fast[1].angular_speed = 1.0;
    let mut normal = motion;
    normal[1].velocity_at_origin[2] = 0.01;
    for invalid in [[motion[0]; 2], adjacent, fast, normal] {
        assert!(
            original
                .step_with_planar_rigid_wall_pair(invalid, 1e-12, &gate)
                .is_err()
        );
        assert_eq!(original.grid().state_bits(), before);
        assert_eq!(original.ledger_mass(), mass);
    }
    let mut half = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, z| z < 4,
    )
    .unwrap();
    assert_eq!(half.sample_at(0, 0, 0).unwrap().cell, Cell3::Fluid);
    let before = half.grid().state_bits();
    assert!(
        half.step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
            .is_err()
    );
    assert_eq!(half.grid().state_bits(), before);
}

#[test]
fn g4_requested_cancellation_retains_liquid_and_replay_is_deterministic() {
    let gate = CancelGate::new_clock_free();
    let mut initial = tank();
    let motion = motions(1);
    initial
        .step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
        .unwrap();
    let before = initial.grid().state_bits();
    let mass = initial.ledger_mass();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        initial.step_with_planar_rigid_wall_pair(motion, 1e-12, &stopped),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert_eq!(initial.grid().state_bits(), before);
    assert_eq!(initial.ledger_mass(), mass);
    let mut replay = initial.clone();
    let a = initial
        .step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
        .unwrap();
    let b = replay
        .step_with_planar_rigid_wall_pair(motion, 1e-12, &gate)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(initial.grid().state_bits(), replay.grid().state_bits());
}

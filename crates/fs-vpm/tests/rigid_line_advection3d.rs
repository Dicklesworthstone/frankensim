//! Finite rigid line geometry through the original near/far material wake.
use std::cell::Cell;

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{
    FilamentWake, Refusal, RigidLineMotion, RigidWakeBoundary, SelfInducedShed,
    WakeBoundary, WakeRateCertificate, WakeTimeSpec, WakeTreeSpec,
};
use fs_vpm::images3d::ImagePlane;

const DT: f64 = 0.25;
const CORE: f64 = 0.25;

fn ring(gamma: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate { shed_hz: 1.0 / DT, n_stations: 1, max_rows: 64 },
        vec![[0.0, -0.7, 1.0], [0.0, 0.7, 1.0]],
    ).unwrap();
    wake.shed(&[gamma], [0.0; 3]).unwrap();
    for point in &mut wake.rows[0].nodes {
        point[0] = -1.0;
    }
    wake.rows[0].core2_m2 = CORE * CORE;
    wake
}

fn settings(substeps: usize) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [-0.1, 0.0, 0.0],
        emitted_core_radius_m: CORE,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    }
}

fn time() -> WakeTimeSpec {
    WakeTimeSpec {
        absolute_position_tolerance_m: 1e-10,
        relative_position_tolerance: 0.0,
        min_step_s: 1e-8,
        max_trials: 4096,
    }
}

fn rotating(angle: f64) -> RigidWakeBoundary {
    RigidWakeBoundary {
        line_motion: RigidLineMotion {
            pivot_start_m: [0.0, 0.0, 1.0],
            pivot_velocity_mps: [0.05, -0.025, 0.0],
            angular_velocity_world_rad_s: [0.0, 0.0, angle / DT],
        },
        image_plane: None,
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(distance(a, b) <= tolerance, "{a:?} != {b:?}");
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
}

// Independent Rodrigues reference, without the production SO(3) or sampler.
fn reference(motion: RigidLineMotion, point: [f64; 3], elapsed: f64) -> ([f64; 3], [f64; 3]) {
    let omega = motion.angular_velocity_world_rad_s;
    let rate = omega[0].hypot(omega[1]).hypot(omega[2]);
    let arm = std::array::from_fn::<_, 3, _>(|i| point[i] - motion.pivot_start_m[i]);
    let rotated = if rate == 0.0 {
        arm
    } else {
        let axis = omega.map(|v| v / rate);
        let angle = rate * elapsed;
        let dot = axis[0]*arm[0] + axis[1]*arm[1] + axis[2]*arm[2];
        let tangent = cross(axis, arm);
        std::array::from_fn(|i| {
            arm[i]*angle.cos() + tangent[i]*angle.sin() + axis[i]*dot*(1.0-angle.cos())
        })
    };
    let spin = cross(omega, rotated);
    (
        std::array::from_fn(|i| motion.pivot_start_m[i] + motion.pivot_velocity_mps[i]*elapsed + rotated[i]),
        std::array::from_fn(|i| motion.pivot_velocity_mps[i] + spin[i]),
    )
}

fn incident(p: [f64; 3], elapsed: f64) -> Result<[f64; 3], Refusal> {
    // Smooth divergence-free, irrotational strain plus a time-dependent gust.
    Ok([0.03*p[1], 0.03*p[0], 0.02*elapsed])
}

#[test]
fn finite_rotation_preserves_span_and_has_the_analytic_material_velocity() {
    let motion = RigidLineMotion {
        pivot_start_m: [2.0, -1.0, 0.25],
        pivot_velocity_mps: [0.1, 0.2, -0.1],
        angular_velocity_world_rad_s: [0.3, -0.4, 0.5],
    };
    let first = [2.2, -0.3, 0.4];
    let second = [1.7, -1.4, 0.8];
    for elapsed in [0.0, 0.4, 1.7, 40.0] {
        let (point, velocity) = motion.point_and_velocity(first, elapsed).unwrap();
        let expected = reference(motion, first, elapsed);
        close(point, expected.0, 2e-14);
        close(velocity, expected.1, 2e-14);
        let other = motion.point_and_velocity(second, elapsed).unwrap().0;
        assert!((distance(point, other) - distance(first, second)).abs() < 3e-14);
        let pivot = motion.point_and_velocity(motion.pivot_start_m, elapsed).unwrap();
        close(pivot.0, reference(motion, motion.pivot_start_m, elapsed).0, 1e-15);
        assert_eq!(pivot.1, motion.pivot_velocity_mps);
    }
    let epsilon = 1e-6;
    let before = motion.point_and_velocity(first, 0.4-epsilon).unwrap().0;
    let after = motion.point_and_velocity(first, 0.4+epsilon).unwrap().0;
    let derivative = std::array::from_fn(|i| (after[i]-before[i])/(2.0*epsilon));
    close(derivative, motion.point_and_velocity(first, 0.4).unwrap().1, 5e-10);
}

#[test]
fn zero_motion_keeps_legacy_bits_and_translation_keeps_the_original_line_law() {
    let original = ring(0.4);
    for tree in [None, Some(WakeTreeSpec { absolute_velocity_error_mps: 0.0, leaf_capacity: 2 })] {
        let mut legacy = original.clone();
        let mut rigid = original.clone();
        if let Some(tree) = tree {
            legacy.shed_self_induced_with_tree(
                &[0.6], settings(16), WakeBoundary::default(), tree, incident, || false,
            ).unwrap();
        } else {
            legacy.shed_self_induced_with_environment(
                &[0.6], settings(16), WakeBoundary::default(), incident, || false,
            ).unwrap();
        }
        rigid.shed_self_induced_with_rigid_line(
            &[0.6], settings(16), RigidWakeBoundary::default(), tree, incident, || false,
        ).unwrap();
        assert_eq!(rigid, legacy);
    }
    let mut boundary = rotating(0.0);
    boundary.line_motion.pivot_start_m = [1e100, -1e100, 1e100];
    let endpoint: Vec<_> = original.line_nodes.iter().map(|&p| {
        std::array::from_fn(|i| p[i] + boundary.line_motion.pivot_velocity_mps[i]*DT)
    }).collect();
    let mut legacy = original.clone();
    let mut rigid = original;
    legacy.shed_self_induced_with_environment(
        &[0.6], settings(16), WakeBoundary { line_endpoint_m: Some(&endpoint), image_plane: None },
        incident, || false,
    ).unwrap();
    rigid.shed_self_induced_with_rigid_line(
        &[0.6], settings(16), boundary, None, incident, || false,
    ).unwrap();
    for (a, b) in legacy.rows.iter().zip(&rigid.rows) {
        for (&p, &q) in a.nodes.iter().zip(&b.nodes) { close(p, q, 2e-14); }
    }
    assert_eq!(rigid.line_nodes, endpoint);
}

#[test]
fn native_midpoint_stages_use_rigid_source_geometry_and_world_fluid_velocity() {
    let original = ring(0.4);
    let boundary = rotating(1.2);
    let config = settings(8);
    let mut expected = original.clone();
    let h = DT / config.substeps as f64;
    for step in 0..config.substeps {
        let elapsed = step as f64*h;
        expected.line_nodes = original.line_nodes.iter()
            .map(|&p| reference(boundary.line_motion, p, elapsed).0).collect();
        let start = expected.rows[0].nodes.clone();
        let first: Vec<_> = start.iter().map(|&p| {
            let induced = expected.induced_velocity(p);
            let external = incident(p, elapsed).unwrap();
            std::array::from_fn::<_, 3, _>(|i| induced[i]+config.ambient_velocity_mps[i]+external[i])
        }).collect();
        for ((point, &p), velocity) in expected.rows[0].nodes.iter_mut().zip(&start).zip(first) {
            *point = std::array::from_fn(|i| p[i]+0.5*h*velocity[i]);
        }
        expected.line_nodes = original.line_nodes.iter()
            .map(|&p| reference(boundary.line_motion, p, elapsed+0.5*h).0).collect();
        let second: Vec<_> = expected.rows[0].nodes.iter().map(|&p| {
            let induced = expected.induced_velocity(p);
            let external = incident(p, elapsed+0.5*h).unwrap();
            std::array::from_fn::<_, 3, _>(|i| induced[i]+config.ambient_velocity_mps[i]+external[i])
        }).collect();
        for ((point, &p), velocity) in expected.rows[0].nodes.iter_mut().zip(&start).zip(second) {
            *point = std::array::from_fn(|i| p[i]+h*velocity[i]);
        }
    }
    let mut actual = original.clone();
    actual.shed_self_induced_with_rigid_line(
        &[7.0], config, boundary, None, incident, || false,
    ).unwrap();
    for (&p, &q) in actual.rows[0].nodes.iter().zip(&expected.rows[0].nodes) { close(p, q, 3e-14); }
    for (&p, &start) in actual.line_nodes.iter().zip(&original.line_nodes) {
        close(p, reference(boundary.line_motion, start, DT).0, 3e-15);
    }
    assert_eq!(actual.rows[0].gamma, original.rows[0].gamma);
    assert_eq!(actual.rows[0].age_ticks, original.rows[0].age_ticks+1.0);
    assert_eq!(actual.rows[1].gamma, [7.0]);
    assert_eq!(actual.rows[1].nodes, actual.line_nodes);
    assert_eq!(actual.rows[1].age_ticks, 0.0);
    for row in 0..=actual.rows.len() {
        for edge in 0..=actual.cert.n_stations {
            assert!(actual.node_circulation_residual(row, edge).unwrap().abs() < 2e-15);
        }
    }
}

#[test]
fn actual_rotating_wake_refines_and_differs_from_endpoint_chord_interpolation() {
    let original = ring(0.4);
    let boundary = rotating(1.2);
    let solve = |substeps| {
        let mut wake = original.clone();
        wake.shed_self_induced_with_rigid_line(
            &[0.6], settings(substeps), boundary, None, incident, || false,
        ).unwrap();
        wake
    };
    let fine = solve(512);
    let errors: Vec<_> = [8, 16, 32].into_iter().map(|n| {
        distance(solve(n).rows[0].nodes[0], fine.rows[0].nodes[0])
    }).collect();
    for pair in errors.windows(2) {
        assert!((3.0..5.5).contains(&(pair[0]/pair[1])), "{errors:?}");
    }
    let endpoint = fine.line_nodes.clone();
    let mut linear = original;
    linear.shed_self_induced_with_environment(
        &[0.6], settings(512), WakeBoundary { line_endpoint_m: Some(&endpoint), image_plane: None },
        incident, || false,
    ).unwrap();
    assert!(distance(linear.rows[0].nodes[0], fine.rows[0].nodes[0]) > 1e-5);
    // Both paths have the same endpoint: the changed wake comes from the
    // physically different intermediate lifting-line geometry.
    assert_eq!(linear.line_nodes, fine.line_nodes);
}

#[test]
fn full_turn_uses_arc_travel_and_adaptive_resolution_can_refine_it() {
    let original = ring(0.2);
    let mut boundary = rotating(std::f64::consts::TAU);
    boundary.line_motion.pivot_velocity_mps = [0.0; 3];
    let mut coarse = original.clone();
    assert_eq!(coarse.shed_self_induced_with_rigid_line(
        &[0.2], settings(1), boundary, None, |_, _| Ok([0.0; 3]), || false,
    ).unwrap_err().code, "wake-advection-displacement");
    assert_eq!(coarse, original);
    let mut adaptive = original.clone();
    let report = adaptive.shed_self_induced_adaptive_with_rigid_line(
        &[0.2], settings(1), time(), boundary, None, |_, _| Ok([0.0; 3]), || false,
    ).unwrap();
    assert!(report.rejected_trials > 0);
    assert!(report.accepted_substeps > 1);
    assert!(report.max_accepted_step_s <= CORE/(0.7*std::f64::consts::TAU/DT));
    for (&p, &q) in adaptive.line_nodes.iter().zip(&original.line_nodes) { close(p, q, 2e-15); }
    assert!(distance(adaptive.rows[0].nodes[0], original.rows[0].nodes[0]) > 1e-3);
    assert_eq!(adaptive.rows.len(), original.rows.len()+1);
}

fn sheet() -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate { shed_hz: 1.0/DT, n_stations: 3, max_rows: 64 },
        [-0.75, -0.25, 0.25, 0.75].map(|y| [0.0, y, 1.0]).to_vec(),
    ).unwrap();
    for row in 0..4 {
        let r = row as f64;
        wake.shed(&[0.2+0.03*r, 0.4+0.02*r, 0.1+0.01*r], [0.0; 3]).unwrap();
    }
    for (row, value) in wake.rows.iter_mut().enumerate() {
        for point in &mut value.nodes { point[0] = -0.4*(4-row) as f64; }
        value.core2_m2 = CORE*CORE;
    }
    wake
}

fn partition(full: &FilamentWake) -> (FilamentWake, FarField) {
    let mut near = full.clone();
    let far = FarField::aggregate(&mut near, 2, 1, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    (near, far)
}

fn same_system(near: &FilamentWake, far: &FarField, full: &FilamentWake) {
    let (expected, expected_far) = partition(full);
    assert_eq!(near.line_nodes, expected.line_nodes);
    assert_eq!(near.rows.len(), expected.rows.len());
    for (a, b) in near.rows.iter().zip(&expected.rows) {
        assert_eq!(a.gamma, b.gamma);
        assert_eq!(a.core2_m2, b.core2_m2);
        assert_eq!(a.age_ticks, b.age_ticks);
        for (&p, &q) in a.nodes.iter().zip(&b.nodes) { close(p, q, 3e-11); }
    }
    assert_eq!(far.cells.len(), expected_far.cells.len());
    for (a, b) in far.cells.iter().zip(&expected_far.cells) {
        assert_eq!(a.source_core2_m2, b.source_core2_m2);
        assert_eq!(a.segs.len(), b.segs.len());
        for ((a0, a1, gamma), (b0, b1, expected_gamma)) in a.segs.iter().zip(&b.segs) {
            close(*a0, *b0, 3e-11);
            close(*a1, *b1, 3e-11);
            assert_eq!(gamma, expected_gamma);
        }
        close(a.exact([2.0, 1.5, 3.0]), b.exact([2.0, 1.5, 3.0]), 3e-11);
    }
}

#[test]
fn fixed_and_adaptive_rigid_motion_share_the_connected_near_far_sources() {
    for adaptive in [false, true] {
        for tree in [None, Some(WakeTreeSpec { absolute_velocity_error_mps: 0.0, leaf_capacity: 2 })] {
            let mut full = sheet();
            let (mut near, mut far) = partition(&full);
            let boundary = rotating(0.5);
            if adaptive {
                let a = full.shed_self_induced_adaptive_with_rigid_line(
                    &[0.3, 0.4, 0.2], settings(8), time(), boundary, tree, incident, || false,
                ).unwrap();
                let b = far.shed_self_induced_adaptive_with_rigid_line(
                    &mut near, &[0.3, 0.4, 0.2], settings(8), time(), boundary, tree, incident, || false,
                ).unwrap();
                assert_eq!(a.accepted_substeps, b.accepted_substeps);
            } else {
                full.shed_self_induced_with_rigid_line(
                    &[0.3, 0.4, 0.2], settings(32), boundary, tree, incident, || false,
                ).unwrap();
                far.shed_self_induced_with_rigid_line(
                    &mut near, &[0.3, 0.4, 0.2], settings(32), boundary, tree, incident, || false,
                ).unwrap();
            }
            same_system(&near, &far, &full);
        }
    }
}

#[test]
fn invalid_motion_sampled_wall_crossing_and_late_cancellation_publish_nothing() {
    let original = ring(0.2);
    let mut invalid = rotating(0.5);
    invalid.line_motion.angular_velocity_world_rad_s[0] = f64::NAN;
    let mut wake = original.clone();
    assert_eq!(wake.shed_self_induced_with_rigid_line(
        &[0.3], settings(32), invalid, None, incident, || false,
    ).unwrap_err().code, "wake-advection-line-invalid");
    assert_eq!(wake, original);

    // Start and endpoint clear the wall, but an intermediate rigid rotation
    // crosses it. These are sampled checks, not a swept-volume certificate.
    let crossing = RigidWakeBoundary {
        line_motion: RigidLineMotion {
            pivot_start_m: [0.0, 0.0, 0.4],
            pivot_velocity_mps: [0.0; 3],
            angular_velocity_world_rad_s: [std::f64::consts::TAU/DT, 0.0, 0.0],
        },
        image_plane: Some(ImagePlane { z0_m: 0.0 }),
    };
    assert_eq!(wake.shed_self_induced_with_rigid_line(
        &[0.3], settings(64), crossing, None, |_, _| Ok([0.0; 3]), || false,
    ).unwrap_err().code, "wake-advection-line-invalid");
    assert_eq!(wake, original);

    let (mut near, mut far) = partition(&sheet());
    let before_near = near.clone();
    let before_far = far.clone();
    let reached = Cell::new(0.0);
    assert_eq!(far.shed_self_induced_adaptive_with_rigid_line(
        &mut near, &[0.3, 0.4, 0.2], settings(16), time(), rotating(0.5), None,
        |p, t| { reached.set(t); incident(p, t) }, || reached.get() > 0.5*DT,
    ).unwrap_err().code, "wake-advection-cancelled");
    assert!(reached.get() > 0.5*DT);
    assert_eq!(near, before_near);
    assert_eq!(far, before_far);
}

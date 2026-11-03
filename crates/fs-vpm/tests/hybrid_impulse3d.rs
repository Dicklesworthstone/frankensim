//! G2 independent closed-source integrals, G3 origin/partition laws, G4 refusal.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs};
use fs_vpm::filament3d::{FilamentWake, SelfInducedShed, WakeRateCertificate};
use fs_vpm::images3d::{ImagePlane, ImageSystem};

const CAP: usize = 100_000;

fn empty_near() -> FilamentWake {
    FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 20.0,
            n_stations: 1,
            max_rows: 64,
        },
        vec![[0.0, -0.5, 2.0], [0.0, 0.5, 2.0]],
    )
    .unwrap()
}

fn empty_far() -> FarField {
    FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        shed_dt_s: 0.05,
        cells: Vec::new(),
    }
}

fn sheet() -> FilamentWake {
    let mut near = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 20.0,
            n_stations: 2,
            max_rows: 64,
        },
        vec![[0.0, -1.0, 2.0], [0.1, 0.0, 2.1], [0.0, 1.0, 1.9]],
    )
    .unwrap();
    for r in 0..4 {
        near.shed(
            &[0.3 + r as f64 * 0.1, -0.2 + r as f64 * 0.05],
            [-0.25, 0.02, 0.015],
        )
        .unwrap();
        near.rows.last_mut().unwrap().core2_m2 = 0.01 + 0.03 * r as f64;
    }
    for (r, row) in near.rows.iter_mut().enumerate() {
        for (s, p) in row.nodes.iter_mut().enumerate() {
            p[2] += 0.04 * (r + 1) as f64 * (s as f64 - 1.0);
        }
    }
    near
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(
        (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]) <= tolerance,
        "{a:?} != {b:?}"
    );
}

fn polygon(vertices: &[[f64; 3]], gamma: f64) -> FarField {
    FarField {
        cells: vec![cell_from_segs(
            (0..vertices.len())
                .map(|i| (vertices[i], vertices[(i + 1) % vertices.len()], gamma))
                .collect(),
            0.04,
        )],
        ..empty_far()
    }
}

fn square(center: [f64; 3], gamma: f64) -> FarField {
    polygon(
        &[
            add(center, [-0.5, -0.5, 0.0]),
            add(center, [0.5, -0.5, 0.0]),
            add(center, [0.5, 0.5, 0.0]),
            add(center, [-0.5, 0.5, 0.0]),
        ],
        gamma,
    )
}

// Independent Simpson integration: x×dx is linear along a straight edge,
// and |x-origin|² dx is quadratic, so these three nodes integrate exactly.
fn integrate_loop(vertices: &[[f64; 3]], gamma: f64, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for edge in 0..vertices.len() {
        let a = subtract(vertices[edge], origin);
        let b = subtract(vertices[(edge + 1) % vertices.len()], origin);
        let dl = subtract(b, a);
        for (point, weight) in [(a, 1.0), (add(a, b).map(|v| 0.5 * v), 4.0), (b, 1.0)] {
            let area = cross(point, dl);
            let radius2 = point.iter().map(|v| v * v).sum::<f64>();
            for i in 0..3 {
                linear[i] += gamma * weight * area[i] / 12.0;
                angular[i] -= gamma * weight * radius2 * dl[i] / 12.0;
            }
        }
    }
    (linear, angular)
}

fn full_sheet_oracle(near: &FilamentWake, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (r, row) in near.rows.iter().enumerate() {
        let next = near.rows.get(r + 1).map_or(&near.line_nodes, |r| &r.nodes);
        for s in 0..near.cert.n_stations {
            let (i, j) = integrate_loop(
                &[next[s], next[s + 1], row.nodes[s + 1], row.nodes[s]],
                row.gamma[s],
                origin,
            );
            linear = add(linear, i);
            angular = add(angular, j);
        }
    }
    (linear, angular)
}

#[test]
fn complete_near_and_every_far_partition_match_independent_closed_quad_integrals() {
    let full = sheet();
    let origin = [0.25, -0.5, 1.0];
    let expected = full_sheet_oracle(&full, origin);
    assert_close(
        full.hydrodynamic_impulse_per_density().unwrap(),
        expected.0,
        2e-15,
    );
    assert_close(
        full.hydrodynamic_angular_impulse_per_density(origin)
            .unwrap(),
        expected.1,
        2e-15,
    );
    let unsplit = empty_far()
        .hydrodynamic_impulse_per_density(&full, origin, CAP, || false)
        .unwrap();
    for removed in 1..4 {
        let mut near = full.clone();
        let far = FarField::aggregate(
            &mut near,
            removed,
            1,
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.02 },
            0.05,
        )
        .unwrap();
        let actual = far
            .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
            .unwrap();
        assert_eq!(actual, unsplit, "partition at row {removed}");
        assert_close(actual.linear_per_density_m4ps(), expected.0, 2e-15);
        assert_close(actual.angular_per_density_m5ps(), expected.1, 2e-15);
        assert_eq!(
            near.hydrodynamic_impulse_per_density().unwrap_err().code,
            "impulse-interface-attached"
        );
    }
}

#[test]
fn arbitrary_origins_proper_rotation_and_large_translation_obey_transport_laws() {
    let vertices = [
        [-1.0, -0.5, 0.25],
        [1.0, -0.5, -0.25],
        [0.75, 0.75, 0.5],
        [-0.75, 0.75, -0.5],
    ];
    let near = empty_near();
    let far = polygon(&vertices, 1.25);
    let origin = [0.25, -0.75, 0.5];
    let measured = far
        .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
        .unwrap();
    let direct = integrate_loop(&vertices, 1.25, origin);
    assert_close(measured.linear_per_density_m4ps(), direct.0, 2e-15);
    assert_close(measured.angular_per_density_m5ps(), direct.1, 2e-15);
    let other = [-0.5, 0.25, -1.0];
    let transported = measured.about_origin(other).unwrap();
    assert_eq!(
        transported,
        far.hydrodynamic_impulse_per_density(&near, other, CAP, || false)
            .unwrap()
    );
    assert_close(
        transported.angular_per_density_m5ps(),
        subtract(
            measured.angular_per_density_m5ps(),
            cross(subtract(other, origin), measured.linear_per_density_m4ps()),
        ),
        2e-15,
    );
    let shift = [2.0_f64.powi(40), -2.0_f64.powi(41), 2.0_f64.powi(42)];
    let shifted = polygon(&vertices.map(|p| add(p, shift)), 1.25)
        .hydrodynamic_impulse_per_density(&near, add(origin, shift), CAP, || false)
        .unwrap();
    assert_eq!(
        shifted.linear_per_density_m4ps(),
        measured.linear_per_density_m4ps()
    );
    assert_eq!(
        shifted.angular_per_density_m5ps(),
        measured.angular_per_density_m5ps()
    );
    let turn = |p: [f64; 3]| [p[2], p[0], p[1]];
    let rotated = polygon(&vertices.map(turn), 1.25)
        .hydrodynamic_impulse_per_density(&near, turn(origin), CAP, || false)
        .unwrap();
    assert_close(
        rotated.linear_per_density_m4ps(),
        turn(measured.linear_per_density_m4ps()),
        2e-15,
    );
    assert_close(
        rotated.angular_per_density_m5ps(),
        turn(measured.angular_per_density_m5ps()),
        2e-15,
    );
}

#[test]
fn opposite_rings_retain_pure_angular_impulse_when_linear_impulse_vanishes() {
    let mut far = square([-2.0, 0.0, 0.0], 2.0);
    far.cells.extend(square([2.0, 0.0, 0.0], -2.0).cells);
    let moment = far
        .hydrodynamic_impulse_per_density(&empty_near(), [0.0; 3], CAP, || false)
        .unwrap();
    assert_eq!(moment.linear_per_density_m4ps(), [0.0; 3]);
    assert_close(moment.angular_per_density_m5ps(), [0.0, 8.0, 0.0], 2e-15);
    assert_eq!(
        moment
            .about_origin([1e12, -2e12, 3e12])
            .unwrap()
            .angular_per_density_m5ps(),
        moment.angular_per_density_m5ps()
    );
}

#[test]
fn interval_changes_retain_small_moments_hidden_by_a_distant_reference() {
    let near = empty_near();
    let mut far = square([0.0; 3], 2.0);
    let origin = [1e20, 0.0, 0.0];
    let before = far
        .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
        .unwrap();
    far.translate([0.125, 0.0, 0.0]).unwrap();
    let after = far
        .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
        .unwrap();
    // The individual rounded angular readouts lose the 1/4-unit change.
    assert_eq!(
        before.angular_per_density_m5ps(),
        after.angular_per_density_m5ps()
    );
    let change = after.exchange_since(&before, 0.25, 1.6).unwrap();
    assert_eq!(change.wake_impulse_change_ns, [0.0; 3]);
    assert_close(
        change.wake_angular_impulse_change_nms,
        [0.0, -0.4, 0.0],
        1e-15,
    );
    assert_close(change.mean_reaction_moment_nm, [0.0, 1.6, 0.0], 4e-15);
    let only_reference_changed = before.about_origin([-1e20, 4e20, 7e20]).unwrap();
    let unchanged = only_reference_changed
        .exchange_since(&before, 0.25, 1.6)
        .unwrap();
    assert_eq!(unchanged.mean_reaction_force_n, [0.0; 3]);
    assert_eq!(unchanged.mean_reaction_moment_nm, [0.0; 3]);
}

#[test]
fn density_time_sign_and_common_fixed_origin_define_the_wrench_equivalent() {
    let near = empty_near();
    let before = empty_far()
        .hydrodynamic_impulse_per_density(&near, [0.0; 3], 0, || false)
        .unwrap();
    let after = square([0.5, 0.25, 1.0], 2.0)
        .hydrodynamic_impulse_per_density(&near, [0.0; 3], CAP, || false)
        .unwrap();
    let change = after.exchange_since(&before, 0.25, 0.5).unwrap();
    assert_eq!(change.wake_impulse_change_ns, [0.0, 0.0, 1.0]);
    assert_close(
        change.wake_angular_impulse_change_nms,
        [0.25, -0.5, 0.0],
        2e-15,
    );
    assert_eq!(change.mean_reaction_force_n, [0.0, 0.0, -4.0]);
    assert_close(change.mean_reaction_moment_nm, [-1.0, 2.0, 0.0], 2e-15);
    let origin = [1.0, 2.0, 3.0];
    let shifted = after
        .about_origin(origin)
        .unwrap()
        .exchange_since(&before.about_origin([-2.0, 1.0, -4.0]).unwrap(), 0.25, 0.5)
        .unwrap();
    assert_eq!(shifted.origin_m, origin);
    assert_eq!(shifted.mean_reaction_force_n, change.mean_reaction_force_n);
    assert_close(
        shifted.mean_reaction_moment_nm,
        subtract(
            change.mean_reaction_moment_nm,
            cross(origin, change.mean_reaction_force_n),
        ),
        2e-15,
    );
    for (dt, rho) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
        (1.0, f64::INFINITY),
    ] {
        assert_eq!(
            after.exchange_since(&before, dt, rho).unwrap_err().code,
            "impulse-exchange-invalid"
        );
    }
    assert_eq!(
        after
            .exchange_since(&before, 1.0, f64::MAX)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    assert_eq!(
        after
            .exchange_since(&before, f64::MIN_POSITIVE / 8.0, 1.0)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
}

#[test]
fn actual_supported_wake_motion_has_the_analytic_impulse_change_without_forcing_conservation() {
    let mut near = empty_near();
    near.cert.shed_hz = 4.0;
    near.line_nodes = vec![[0.0, -0.7, 0.0], [0.0, 0.7, 0.0]];
    near.shed(&[4.0], [0.0; 3]).unwrap();
    for p in &mut near.rows[0].nodes {
        p[0] = -1.0;
    }
    near.rows[0].core2_m2 = 0.04;
    let far = empty_far();
    let origin = [0.3, -0.2, 0.4];
    let before = far
        .hydrodynamic_impulse_per_density(&near, [-0.2, 0.4, 0.1], CAP, || false)
        .unwrap();
    near.shed_self_induced(
        &[9.0],
        SelfInducedShed {
            ambient_velocity_mps: [0.0; 3],
            emitted_core_radius_m: 0.2,
            core_mode: WakeCoreEvolutionMode::Frozen,
            substeps: 128,
            max_displacement_fraction: 1.0,
            max_work: 1_000_000,
        },
        || false,
    )
    .unwrap();
    let after = far
        .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
        .unwrap();
    let change = after.exchange_since(&before, 0.25, 1.2).unwrap();
    let angle = 4.0 * 1.4 * 0.25
        / (4.0 * std::f64::consts::PI * (1.0_f64 + 1.4 * 1.4 + 0.04).sqrt())
        * (1.0 / (1.0 + 0.04) + 1.0 / (1.4 * 1.4 + 0.04));
    let expected_force = [
        4.0 * 1.4 * angle.sin(),
        0.0,
        4.0 * 1.4 * (angle.cos() - 1.0),
    ]
    .map(|v| -1.2 * v / 0.25);
    assert_close(change.mean_reaction_force_n, expected_force, 2e-6);
    assert_close(
        change.mean_reaction_moment_nm,
        cross(origin, expected_force).map(|v| -v),
        2e-6,
    );
    assert!(change.wake_impulse_change_ns.iter().any(|v| v.abs() > 0.1));
    assert_eq!(near.rows.last().unwrap().gamma, [9.0]);
}

#[test]
fn core_spreading_multipole_caches_and_image_devices_do_not_enter_physical_impulse() {
    let mut near = sheet();
    let mut far = FarField::aggregate(
        &mut near,
        2,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.03 },
        0.05,
    )
    .unwrap();
    let origin = [0.5, 0.1, 0.0];
    let physical = far
        .hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
        .unwrap();
    let image = ImageSystem::build(&near, &far, ImagePlane { z0_m: 0.0 }).unwrap();
    let image_far = FarField {
        cells: image.far_image_cells().to_vec(),
        ..empty_far()
    };
    let image_moment = image_far
        .hydrodynamic_impulse_per_density(image.near_image_wake(), origin, CAP, || false)
        .unwrap();
    assert!(
        image_moment
            .linear_per_density_m4ps()
            .iter()
            .any(|v| v.abs() > 0.01)
    );
    assert_eq!(
        physical,
        far.hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
            .unwrap()
    );
    far.advance_core(0.75).unwrap();
    for row in &mut near.rows {
        row.age_ticks += 15.0;
        row.core2_m2 += 0.1;
    }
    for cell in &mut far.cells {
        cell.w_sum = [1e10; 3];
        cell.dip = [[1e10; 3]; 3];
    }
    assert_eq!(
        physical,
        far.hydrodynamic_impulse_per_density(&near, origin, CAP, || false)
            .unwrap()
    );
}

#[test]
fn invalid_or_incomplete_sources_caps_and_cancellation_return_no_partial_diagnostic() {
    let original = sheet();
    let mut near = original.clone();
    let far = FarField::aggregate(&mut near, 2, 1, WakeCoreEvolutionMode::Frozen, 0.05).unwrap();
    for case in 0..5 {
        let mut bad_near = near.clone();
        let mut bad_far = far.clone();
        let expected = match case {
            0 => {
                bad_far.cells[0].segs[0].0[0] -= 0.01;
                "impulse-topology-invalid"
            }
            1 => {
                bad_far.cells.remove(0);
                "impulse-topology-invalid"
            }
            2 => {
                bad_near.downstream_gamma = Some(vec![1.0]);
                "impulse-topology-invalid"
            }
            3 => {
                bad_far.cells[0].segs[0].2 = f64::NAN;
                "impulse-geometry-invalid"
            }
            _ => {
                bad_far.cells[0].segs[0].0[2] = f64::INFINITY;
                "impulse-geometry-invalid"
            }
        };
        let digests = (bad_near.digest(), bad_far.digest());
        assert_eq!(
            bad_far
                .hydrodynamic_impulse_per_density(&bad_near, [0.0; 3], CAP, || false)
                .unwrap_err()
                .code,
            expected
        );
        assert_eq!((bad_near.digest(), bad_far.digest()), digests);
    }
    assert_eq!(
        far.hydrodynamic_impulse_per_density(&near, [f64::NAN, 0.0, 0.0], CAP, || false)
            .unwrap_err()
            .code,
        "impulse-origin-invalid"
    );
    let mut polls = 0;
    far.hydrodynamic_impulse_per_density(&near, [0.0; 3], CAP, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert_eq!(
            far.hydrodynamic_impulse_per_density(&near, [0.0; 3], CAP, || {
                seen += 1;
                seen >= stop
            })
            .unwrap_err()
            .code,
            "impulse-cancelled"
        );
    }
    let closed = square([0.0; 3], 1.0);
    assert_eq!(
        closed
            .hydrodynamic_impulse_per_density(&empty_near(), [0.0; 3], 3, || false)
            .unwrap_err()
            .code,
        "impulse-work-exhausted"
    );
    closed
        .hydrodynamic_impulse_per_density(&empty_near(), [0.0; 3], 4, || false)
        .unwrap();
    let empty = empty_far()
        .hydrodynamic_impulse_per_density(&empty_near(), [8.0, 9.0, 10.0], 0, || false)
        .unwrap();
    assert_eq!(empty.linear_per_density_m4ps(), [0.0; 3]);
    assert_eq!(empty.angular_per_density_m5ps(), [0.0; 3]);
    let too_large = polygon(
        &[
            [-1e200, -1e200, 0.0],
            [1e200, -1e200, 0.0],
            [1e200, 1e200, 0.0],
            [-1e200, 1e200, 0.0],
        ],
        1.0,
    );
    assert_eq!(
        too_large
            .hydrodynamic_impulse_per_density(&empty_near(), [0.0; 3], CAP, || false)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
}

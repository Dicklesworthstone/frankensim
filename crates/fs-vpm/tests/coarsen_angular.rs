//! Independent physical-edge angular impulse and bounded conservative coarsening.

use fs_vpm::coarsen3d::coarsen_oldest;
use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{FilamentWake, ShedRow, WakeRateCertificate, segment_velocity};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    assert!(
        norm(sub(actual, expected)) <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

// Two-point Gauss integration of the exact quadratic edge integrand; this
// oracle does not use the production triangle-centroid moment implementation.
fn edge_moments(a: [f64; 3], b: [f64; 3], gamma: f64, origin: [f64; 3]) -> [f64; 6] {
    let a = sub(a, origin);
    let b = sub(b, origin);
    let delta = sub(b, a);
    let midpoint = add(a, b).map(|x| 0.5 * x);
    let offset = delta.map(|x| x / (2.0 * 3.0_f64.sqrt()));
    let q0 = sub(midpoint, offset);
    let q1 = add(midpoint, offset);
    let squared_integral: f64 = (0..3).map(|i| 0.5 * (q0[i] * q0[i] + q1[i] * q1[i])).sum();
    let linear = cross(a, b).map(|x| 0.5 * gamma * x);
    let angular = delta.map(|x| -0.5 * gamma * x * squared_integral);
    [
        linear[0], linear[1], linear[2], angular[0], angular[1], angular[2],
    ]
}
fn moments(wake: &FilamentWake, origin: [f64; 3]) -> [f64; 6] {
    let mut total = [0.0; 6];
    wake.visit_segments_with_core_growth(0.0, |edge| {
        let value = edge_moments(edge.start, edge.end, edge.circulation, origin);
        for i in 0..6 {
            total[i] += value[i];
        }
    });
    total
}
fn strip_moments(wake: &FilamentWake, station: usize, origin: [f64; 3]) -> [f64; 6] {
    let mut strip = wake.clone();
    for row in &mut strip.rows {
        for (s, gamma) in row.gamma.iter_mut().enumerate() {
            if s != station {
                *gamma = 0.0;
            }
        }
    }
    moments(&strip, origin)
}
fn close_moments(actual: [f64; 6], expected: [f64; 6], tolerance: f64) {
    close(
        actual[..3].try_into().unwrap(),
        expected[..3].try_into().unwrap(),
        tolerance,
    );
    close(
        actual[3..].try_into().unwrap(),
        expected[3..].try_into().unwrap(),
        tolerance,
    );
}
fn rectangle(rows: usize) -> FilamentWake {
    let nodes = |x| vec![[x, -1.0, 0.0], [x, 1.0, 0.0]];
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 1,
            max_rows: 32,
        },
        nodes(2.0),
    )
    .unwrap();
    wake.rows = (0..rows)
        .map(|r| ShedRow {
            nodes: nodes(2.0 * r as f64 / rows as f64),
            gamma: vec![1.0],
            core2_m2: 0.0,
            age_ticks: (rows - r) as f64,
        })
        .collect();
    wake
}

#[test]
fn straight_pair_retains_angular_impulse_instead_of_the_old_twenty_five_percent_error() {
    let original = rectangle(2);
    close(
        original.hydrodynamic_impulse_per_density().unwrap(),
        [0.0, 0.0, 4.0],
        1e-14,
    );
    close(
        original
            .hydrodynamic_angular_impulse_per_density([0.0; 3])
            .unwrap(),
        [0.0, -4.0, 0.0],
        1e-14,
    );
    // Exact rectangle centroid × impulse: the previous midpoint construction
    // gave the same I but J_y=-5 instead of -4, despite this elementary geometry.
    let mut old_midpoint = original.clone();
    old_midpoint.rows.truncate(1);
    for node in &mut old_midpoint.rows[0].nodes {
        node[0] = 0.5;
    }
    old_midpoint.rows[0].gamma[0] = 4.0 / 3.0;
    close(
        moments(&old_midpoint, [0.0; 3])[..3].try_into().unwrap(),
        [0.0, 0.0, 4.0],
        1e-14,
    );
    close(
        moments(&old_midpoint, [0.0; 3])[3..].try_into().unwrap(),
        [0.0, -5.0, 0.0],
        1e-14,
    );
    let mut coarse = original.clone();
    assert!(
        !coarsen_oldest(&mut coarse, 1)
            .unwrap()
            .used_stationwise_correction
    );
    assert_eq!(coarse.rows.len(), 1);
    assert_eq!(coarse.rows[0].nodes, original.rows[0].nodes);
    assert!((coarse.rows[0].gamma[0] - 1.0).abs() < 1e-14);
    close_moments(
        moments(&coarse, [0.0; 3]),
        [0.0, 0.0, 4.0, 0.0, -4.0, 0.0],
        1e-14,
    );
}

fn twisted() -> FilamentWake {
    // Dyadic, genuinely nonplanar geometry with six independent retained
    // moment directions, rather than a planar rank-two-only solver oracle.
    let nodes = |r: usize| {
        (0..=2)
            .map(|s: usize| {
                [
                    0.5 * r as f64 + ((7 * r + 3 * s) % 11) as f64 / 64.0,
                    s as f64 - 1.0 + ((5 * r + 7 * s) % 13) as f64 / 64.0,
                    ((11 * r + 17 * s) % 19) as f64 / 32.0,
                ]
            })
            .collect()
    };
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 2,
            max_rows: 32,
        },
        nodes(14),
    )
    .unwrap();
    wake.rows = (0..14)
        .map(|r| ShedRow {
            nodes: nodes(r),
            gamma: vec![
                1.0 + (r % 7) as f64 / 8.0,
                -0.5 - ((3 * r) % 5) as f64 / 8.0,
            ],
            core2_m2: 1.0 / 1024.0,
            age_ticks: (14 - r) as f64,
        })
        .collect();
    wake
}
fn map_geometry(wake: &mut FilamentWake, map: impl Fn([f64; 3]) -> [f64; 3]) {
    for node in wake
        .line_nodes
        .iter_mut()
        .chain(wake.rows.iter_mut().flat_map(|r| &mut r.nodes))
    {
        *node = map(*node);
    }
}

#[test]
fn six_constraint_correction_preserves_each_twisted_station_and_frame_covariance() {
    let original = twisted();
    let origin = [0.75, -0.5, 0.25];
    let mut reference = original.clone();
    let expected: Vec<_> = (0..2)
        .map(|s| strip_moments(&original, s, origin))
        .collect();
    assert!(
        coarsen_oldest(&mut reference, 6)
            .unwrap()
            .used_stationwise_correction
    );
    assert_eq!(reference.rows.len(), 8);
    for s in 0..2 {
        close_moments(strip_moments(&reference, s, origin), expected[s], 2e-11);
    }
    close(
        reference
            .hydrodynamic_angular_impulse_per_density(origin)
            .unwrap(),
        moments(&original, origin)[3..].try_into().unwrap(),
        2e-11,
    );
    let rotate = |p: [f64; 3]| [p[2], p[0], p[1]];
    let translation = [16.0, -8.0, 4.0];
    let transform = |p| add(rotate(p), translation);
    let mut moved = original.clone();
    map_geometry(&mut moved, transform);
    assert!(
        coarsen_oldest(&mut moved, 6)
            .unwrap()
            .used_stationwise_correction
    );
    close(
        moved.hydrodynamic_impulse_per_density().unwrap(),
        rotate(reference.hydrodynamic_impulse_per_density().unwrap()),
        2e-11,
    );
    close(
        moved
            .hydrodynamic_angular_impulse_per_density(transform(origin))
            .unwrap(),
        rotate(
            reference
                .hydrodynamic_angular_impulse_per_density(origin)
                .unwrap(),
        ),
        2e-11,
    );
    for (actual, expected) in moved.rows.iter().zip(&reference.rows) {
        for (a, b) in actual.gamma.iter().zip(&expected.gamma) {
            assert!((a - b).abs() < 2e-10);
        }
    }
}

#[test]
fn insufficient_candidate_geometry_refuses_even_for_zero_net_linear_impulse() {
    for gamma in [[1.0, 2.0], [1.0, -1.0]] {
        let mut wake = rectangle(2);
        for (row, value) in wake.rows.iter_mut().zip(gamma) {
            row.gamma[0] = value;
        }
        if gamma[1] == -1.0 {
            close(
                wake.hydrodynamic_impulse_per_density().unwrap(),
                [0.0; 3],
                1e-14,
            );
            close(
                wake.hydrodynamic_angular_impulse_per_density([0.0; 3])
                    .unwrap(),
                [0.0, 2.0, 0.0],
                1e-14,
            );
        }
        let before = wake.clone();
        assert_eq!(
            coarsen_oldest(&mut wake, 1).unwrap_err().code,
            "coarsen-impulse-unrepresentable"
        );
        assert_eq!(
            wake, before,
            "failure of the bounded constructions must retain all original rows"
        );
    }
}

#[test]
fn angular_origin_transport_remains_accurate_far_from_world_zero() {
    let original = rectangle(2);
    let impulse = original.hydrodynamic_impulse_per_density().unwrap();
    let angular = original
        .hydrodynamic_angular_impulse_per_density([0.0; 3])
        .unwrap();
    let shift = [2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38)];
    let mut moved = original.clone();
    map_geometry(&mut moved, |p| add(p, shift));
    close(
        moved
            .hydrodynamic_angular_impulse_per_density(shift)
            .unwrap(),
        angular,
        1e-14,
    );
    let world_expected = add(angular, cross(shift, impulse));
    close(
        moved
            .hydrodynamic_angular_impulse_per_density([0.0; 3])
            .unwrap(),
        world_expected,
        2e-15 * norm(world_expected),
    );
    coarsen_oldest(&mut moved, 1).unwrap();
    close(
        moved
            .hydrodynamic_angular_impulse_per_density(shift)
            .unwrap(),
        angular,
        1e-14,
    );
    let offset = [0.5, -0.25, 0.125];
    close(
        moved
            .hydrodynamic_angular_impulse_per_density(add(shift, offset))
            .unwrap(),
        sub(angular, cross(offset, impulse)),
        1e-14,
    );
}

fn geometric_velocity(wake: &FilamentWake, p: [f64; 3]) -> [f64; 3] {
    let mut velocity = [0.0; 3];
    wake.visit_segments_with_core_growth(0.0, |edge| {
        let value = segment_velocity(p, edge.start, edge.end);
        for i in 0..3 {
            velocity[i] += edge.circulation * value[i];
        }
    });
    velocity
}

#[test]
fn straight_refinement_and_near_far_partition_retain_the_complete_geometric_moments() {
    let single = rectangle(1);
    let origin = [0.5, -0.25, 0.125];
    let expected = moments(&single, origin);
    for subdivisions in [2, 4, 8, 16] {
        let mut wake = rectangle(subdivisions);
        coarsen_oldest(&mut wake, subdivisions / 2).unwrap();
        close_moments(moments(&wake, origin), expected, 2e-14);
        // Geometric filaments are identical under straight subdivision. The
        // deliberately enlarged coarse core remains a separate field model;
        // coarsening reports its actual regularized-field error honestly.
        for probe in [[0.25, 0.125, 0.75], [3.0, -2.0, 4.0], [16.0, 8.0, -4.0]] {
            close(
                geometric_velocity(&wake, probe),
                geometric_velocity(&single, probe),
                2e-14,
            );
        }
        if wake.rows.len() < 2 {
            continue;
        }
        let far =
            FarField::aggregate(&mut wake, 1, 1, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
        assert_eq!(
            wake.hydrodynamic_angular_impulse_per_density(origin)
                .unwrap_err()
                .code,
            "impulse-interface-attached"
        );
        let mut total = moments(&wake, origin);
        for cell in &far.cells {
            for &(a, b, gamma) in &cell.segs {
                let value = edge_moments(a, b, gamma, origin);
                for i in 0..6 {
                    total[i] += value[i];
                }
            }
        }
        close_moments(total, expected, 2e-14);
    }
}

#[test]
fn angular_diagnostic_rejects_invalid_origin_and_geometry() {
    let mut wake = rectangle(2);
    assert_eq!(
        wake.hydrodynamic_angular_impulse_per_density([f64::NAN, 0.0, 0.0])
            .unwrap_err()
            .code,
        "impulse-origin-invalid"
    );
    wake.rows[0].gamma.clear();
    assert_eq!(
        wake.hydrodynamic_angular_impulse_per_density([0.0; 3])
            .unwrap_err()
            .code,
        "filament-geometry-invalid"
    );
}

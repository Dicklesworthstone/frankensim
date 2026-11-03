//! G0/G1/G3: geometric vector impulse, closed-ring field and complete domains.
//! These oracles do not use the algebraic cell-circulation identity.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{FilamentWake, ShedRow, WakeRateCertificate};

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|x| s * x)
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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
fn wake(line: Vec<[f64; 3]>, rows: Vec<ShedRow>) -> FilamentWake {
    let mut result = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: line.len() - 1,
            max_rows: 64,
        },
        line,
    )
    .unwrap();
    result.rows = rows;
    result
}
fn rectangle(a: f64, b: f64, gamma: f64, core2: f64, subdivisions: usize) -> FilamentWake {
    wake(
        vec![[a, -b, 0.0], [a, b, 0.0]],
        (0..subdivisions)
            .map(|r| {
                let x = -a + 2.0 * a * r as f64 / subdivisions as f64;
                ShedRow {
                    nodes: vec![[x, -b, 0.0], [x, b, 0.0]],
                    gamma: vec![gamma],
                    core2_m2: core2,
                    age_ticks: (subdivisions - r) as f64,
                }
            })
            .collect(),
    )
}

#[test]
fn g1_closed_rectangle_matches_independent_axis_field_and_vector_impulse() {
    for (a, b) in [(1.5, 0.75), (0.25, 2.0)] {
        for core2 in [0.0, 0.04, 0.7] {
            let gamma = 1.2;
            let w = rectangle(a, b, gamma, core2, 1);
            for z in [0.0, 0.3, 2.5] {
                let q = z * z + core2;
                let expected = gamma * a * b / (std::f64::consts::PI * (a * a + b * b + q).sqrt())
                    * (1.0 / (a * a + q) + 1.0 / (b * b + q));
                close(
                    w.induced_velocity([0.0, 0.0, z]),
                    [0.0, 0.0, expected],
                    2e-14,
                );
            }
            close(
                w.hydrodynamic_impulse_per_density().unwrap(),
                [0.0, 0.0, 4.0 * a * b * gamma],
                2e-14,
            );
        }
    }
}

#[test]
fn g0_zero_area_fresh_row_cannot_create_velocity_or_impulse() {
    let mut w = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 2,
            max_rows: 8,
        },
        vec![[0.0, -1.0, 0.0], [0.0, 0.0, 0.25], [0.0, 1.0, 0.0]],
    )
    .unwrap();
    w.shed(&[1.0, -0.5], [0.0; 3]).unwrap();
    close(w.induced_velocity([0.5, 0.1, 0.7]), [0.0; 3], 1e-16);
    assert_eq!(w.hydrodynamic_impulse_per_density().unwrap(), [0.0; 3]);
}

#[test]
fn g3_constant_ring_subdivision_changes_neither_field_nor_impulse() {
    for core2 in [0.0, 0.125] {
        let original = rectangle(2.0, 1.0, 0.75, core2, 1);
        for n in [2, 4, 8, 16] {
            let refined = rectangle(2.0, 1.0, 0.75, core2, n);
            for p in [[0.2, 0.1, 0.3], [2.4, -1.3, 0.75], [-5.0, 3.0, 2.0]] {
                close(
                    refined.induced_velocity(p),
                    original.induced_velocity(p),
                    2e-14,
                );
            }
            close(
                refined.hydrodynamic_impulse_per_density().unwrap(),
                original.hydrodynamic_impulse_per_density().unwrap(),
                2e-14,
            );
        }
    }
}

fn warped() -> FilamentWake {
    let nodes = |r: usize| {
        let x = r as f64 * 0.5;
        vec![
            [x, -1.0, 0.125 * x],
            [x + 0.125, 0.0, 0.25 + 0.125 * x * x],
            [x, 1.0, -0.125 * x],
        ]
    };
    wake(
        nodes(6),
        (0..6)
            .map(|r| ShedRow {
                nodes: nodes(r),
                gamma: vec![1.0 + 0.25 * r as f64, -0.5 + 0.125 * r as f64],
                core2_m2: 0.03125 * (r + 1) as f64,
                age_ticks: (6 - r) as f64,
            })
            .collect(),
    )
}

#[test]
fn g0_complete_domain_impulse_matches_retained_physical_edge_integrals() {
    let full = warped();
    for rows_per_cell in [1, 2, 4] {
        let mut near = full.clone();
        let far = FarField::aggregate(
            &mut near,
            4,
            rows_per_cell,
            WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.25 },
            0.015625,
        )
        .unwrap();
        assert_eq!(
            near.hydrodynamic_impulse_per_density().unwrap_err().code,
            "impulse-interface-attached"
        );
        // Only the combined source domain is closed. All contributions must
        // use the same origin, including open near and far subsets.
        for origin in [[0.0; 3], [16.0, -8.0, 4.0]] {
            let mut impulse = [0.0; 3];
            let mut integrate = |a, b, gamma| {
                impulse = add(
                    impulse,
                    scale(cross(sub(a, origin), sub(b, origin)), 0.5 * gamma),
                );
            };
            near.visit_segments_with_core_growth(0.0, |s| integrate(s.start, s.end, s.circulation));
            for cell in &far.cells {
                for &(a, b, gamma) in &cell.segs {
                    integrate(a, b, gamma);
                }
            }
            close(
                impulse,
                full.hydrodynamic_impulse_per_density().unwrap(),
                3e-14,
            );
        }
    }
}

#[test]
fn g3_impulse_and_velocity_transform_with_the_physical_frame() {
    let w = warped();
    // Exact proper rotation cyclically permuting axes; offset is exactly
    // representable for the fixture, so this isolates frame handling.
    let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
    let translation = [16.0, -8.0, 4.0];
    let transform = |x| add(rotate(x), translation);
    let mut moved = w.clone();
    for node in &mut moved.line_nodes {
        *node = transform(*node);
    }
    for row in &mut moved.rows {
        for node in &mut row.nodes {
            *node = transform(*node);
        }
    }
    close(
        moved.hydrodynamic_impulse_per_density().unwrap(),
        rotate(w.hydrodynamic_impulse_per_density().unwrap()),
        2e-14,
    );
    let p = [1.0, 0.125, 1.25];
    close(
        moved.induced_velocity(transform(p)),
        rotate(w.induced_velocity(p)),
        2e-14,
    );
}

#[test]
fn g0_invalid_geometry_refuses_instead_of_producing_partial_impulse() {
    let mut w = warped();
    w.rows[0].gamma.pop();
    assert_eq!(
        w.hydrodynamic_impulse_per_density().unwrap_err().code,
        "filament-geometry-invalid"
    );
    let mut w = rectangle(1.0, 1.0, f64::MAX, 0.0, 1);
    assert_eq!(
        w.hydrodynamic_impulse_per_density().unwrap_err().code,
        "impulse-unrepresentable"
    );
    w.rows.clear();
    assert_eq!(w.hydrodynamic_impulse_per_density().unwrap(), [0.0; 3]);
}

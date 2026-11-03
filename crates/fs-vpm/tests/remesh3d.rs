//! Physical edge integrals and finite-core field invariants under refinement.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, hybrid_velocity};
use fs_vpm::filament3d::{FilamentWake, ShedRow, WakeRateCertificate};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
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

// Simpson integration is exact for the quadratic angular integrand. This
// integrates the physical edges, independently of production quad moments.
fn edge_moments(a: [f64; 3], b: [f64; 3], gamma: f64, origin: [f64; 3]) -> [f64; 6] {
    let a = sub(a, origin);
    let b = sub(b, origin);
    let squared = |p: [f64; 3]| p.iter().map(|x| x * x).sum::<f64>();
    let mid = std::array::from_fn(|i| 0.5 * (a[i] + b[i]));
    let integral = (squared(a) + 4.0 * squared(mid) + squared(b)) / 6.0;
    let linear = cross(a, b).map(|v| 0.5 * gamma * v);
    let angular = sub(b, a).map(|v| -0.5 * gamma * v * integral);
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
fn close_moments(actual: [f64; 6], expected: [f64; 6]) {
    close(
        actual[..3].try_into().unwrap(),
        expected[..3].try_into().unwrap(),
        2e-12,
    );
    close(
        actual[3..].try_into().unwrap(),
        expected[3..].try_into().unwrap(),
        2e-12,
    );
}
fn warped() -> FilamentWake {
    let nodes = |x: f64| {
        vec![
            [x, -1.0 + 0.125 * x, 0.25 * x],
            [x + 0.125, 0.125 * x * x, 0.5 + 0.125 * x * x],
            [x, 1.0 + 0.125 * x, -0.125 * x],
        ]
    };
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 120.0,
            n_stations: 2,
            max_rows: 64,
        },
        nodes(0.0),
    )
    .unwrap();
    wake.rows = [-3.0, -1.0, -0.25]
        .into_iter()
        .enumerate()
        .map(|(r, x)| ShedRow {
            nodes: nodes(x),
            gamma: vec![1.25 - 0.75 * r as f64, -0.625 + 0.375 * r as f64],
            core2_m2: 0.03125 * (r + 1) as f64,
            age_ticks: [5.75, 2.25, 0.5][r],
        })
        .collect();
    wake.ticks = 17;
    wake
}
fn straight(length: f64, max_rows: usize) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 120.0,
            n_stations: 1,
            max_rows,
        },
        vec![[length, -1.0, 0.0], [length, 1.0, 0.0]],
    )
    .unwrap();
    wake.rows.push(ShedRow {
        nodes: vec![[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
        gamma: vec![0.75],
        core2_m2: 0.125,
        age_ticks: 2.5,
    });
    wake
}

#[test]
fn warped_sheet_refinement_preserves_all_six_impulses_per_station_and_kelvin_balance() {
    let original = warped();
    let origin = [0.375, -0.2, 0.4];
    let mut refined = original.clone();
    let report = refined.refine_streamwise(0.3, 40, || false).unwrap();
    assert!(report.rows_after > report.rows_before);
    for station in 0..=original.cert.n_stations {
        let mut before = original.clone();
        let mut after = refined.clone();
        if station < original.cert.n_stations {
            for row in before.rows.iter_mut().chain(&mut after.rows) {
                for (s, gamma) in row.gamma.iter_mut().enumerate() {
                    if s != station {
                        *gamma = 0.0;
                    }
                }
            }
        }
        let expected = moments(&before, origin);
        close_moments(moments(&after, origin), expected);
        close(
            after.hydrodynamic_impulse_per_density().unwrap(),
            expected[..3].try_into().unwrap(),
            2e-12,
        );
        close(
            after
                .hydrodynamic_angular_impulse_per_density(origin)
                .unwrap(),
            expected[3..].try_into().unwrap(),
            2e-12,
        );
    }
    for r in 0..=refined.rows.len() {
        for edge in 0..=refined.cert.n_stations {
            assert_eq!(refined.node_circulation_residual(r, edge), Some(0.0));
        }
    }
}

#[test]
fn finite_core_velocity_gradient_and_vorticity_survive_without_resetting_material_age() {
    let mut refined = warped();
    let original = refined.clone();
    let buffers: Vec<_> = refined
        .rows
        .iter()
        .map(|r| (r.nodes.as_ptr(), r.gamma.as_ptr()))
        .collect();
    refined.refine_streamwise(0.4, 40, || false).unwrap();
    assert_eq!(refined.cert, original.cert);
    assert_eq!(refined.line_nodes, original.line_nodes);
    assert_eq!(refined.ticks, original.ticks);
    for row in &refined.rows {
        let parent = original.rows.iter().find(|r| r.gamma == row.gamma).unwrap();
        assert_eq!(row.core2_m2.to_bits(), parent.core2_m2.to_bits());
        assert_eq!(row.age_ticks.to_bits(), parent.age_ticks.to_bits());
    }
    for (old, pointers) in original.rows.iter().zip(buffers) {
        let retained = refined.rows.iter().find(|r| r.nodes == old.nodes).unwrap();
        assert_eq!(retained, old);
        assert_eq!((retained.nodes.as_ptr(), retained.gamma.as_ptr()), pointers);
    }
    for growth in [0.0, 0.013, 0.125] {
        for probe in [[-2.0, 0.2, 0.7], [-0.25, -1.0, -0.2], [4.0, -3.0, 2.0]] {
            let before = original
                .induced_kinematics_with_core_growth(probe, growth)
                .unwrap();
            let after = refined
                .induced_kinematics_with_core_growth(probe, growth)
                .unwrap();
            close(after.velocity, before.velocity, 2e-13);
            close(after.vorticity(), before.vorticity(), 2e-12);
            for i in 0..3 {
                close(after.gradient[i], before.gradient[i], 2e-12);
            }
        }
    }
}

#[test]
fn target_and_both_row_budgets_are_enforced_and_refinement_is_idempotent() {
    let original = straight(4.0, 4);
    for budget in [0, 2] {
        let mut wake = original.clone();
        assert_eq!(
            wake.refine_streamwise(1.0, budget, || false)
                .unwrap_err()
                .code,
            "wake-refinement-budget"
        );
        assert_eq!(wake, original);
    }
    let mut wake = original.clone();
    let report = wake.refine_streamwise(1.0, 3, || false).unwrap();
    assert_eq!((report.rows_before, report.rows_after), (1, 4));
    assert_eq!(
        (report.longest_edge_before_m, report.longest_edge_after_m),
        (4.0, 1.0)
    );
    let refined = wake.clone();
    wake.refine_streamwise(1.0, 0, || false).unwrap();
    assert_eq!(wake, refined);
    let mut capped = original.clone();
    capped.cert.max_rows = 3;
    let snapshot = capped.clone();
    assert_eq!(
        capped
            .refine_streamwise(1.0, usize::MAX, || false)
            .unwrap_err()
            .code,
        "wake-refinement-budget"
    );
    assert_eq!(capped, snapshot);

    let mut fractional = straight(1.0, 32);
    fractional.refine_streamwise(0.1, 20, || false).unwrap();
    let snapshot = fractional.clone();
    fractional.refine_streamwise(0.1, 0, || false).unwrap();
    assert_eq!(fractional, snapshot);
}

#[test]
fn empty_and_zero_length_rows_are_stable_but_zero_strength_does_not_skip_resolution() {
    let mut wake = straight(4.0, 16);
    wake.rows[0].gamma[0] = 0.0;
    wake.shed(&[0.0], [0.0; 3]).unwrap();
    let newest = wake.rows.last().unwrap().clone();
    let report = wake.refine_streamwise(1.0, 3, || false).unwrap();
    assert_eq!(report.rows_after, 5);
    assert_eq!(wake.rows.last().unwrap(), &newest);
    assert_eq!(wake.induced_velocity([0.3, 0.4, 0.5]), [0.0; 3]);
    wake.rows.clear();
    let empty = wake.clone();
    let report = wake.refine_streamwise(1.0, 0, || false).unwrap();
    assert_eq!(report.longest_edge_after_m, 0.0);
    assert_eq!(wake, empty);
}

#[test]
fn refinement_retains_an_attached_far_boundary_and_complete_hybrid_moments() {
    let mut near = warped();
    let far =
        FarField::aggregate(&mut near, 1, 1, WakeCoreEvolutionMode::Frozen, 1.0 / 120.0).unwrap();
    let before = near.clone();
    let seam = near.rows[0].nodes.clone();
    let origin = [0.375, -0.2, 0.4];
    let combined = |wake: &FilamentWake| {
        let mut total = moments(wake, origin);
        for cell in &far.cells {
            for &(a, b, gamma) in &cell.segs {
                let value = edge_moments(a, b, gamma, origin);
                for i in 0..6 {
                    total[i] += value[i];
                }
            }
        }
        total
    };
    let expected = combined(&before);
    near.refine_streamwise(0.2, 40, || false).unwrap();
    assert_eq!(near.rows[0].nodes, seam);
    assert_eq!(near.downstream_gamma, before.downstream_gamma);
    close_moments(combined(&near), expected);
    for p in [[-1.1, 0.2, 0.7], [4.0, -3.0, 2.0]] {
        close(
            hybrid_velocity(&near, &far, p),
            hybrid_velocity(&before, &far, p),
            2e-13,
        );
    }
}

#[test]
fn cancellation_at_every_checkpoint_preserves_every_buffer_and_state() {
    let original = warped();
    let mut completed = original.clone();
    let mut polls = 0;
    completed
        .refine_streamwise(0.3, 40, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut wake = original.clone();
        let pointer = wake.rows[0].nodes.as_ptr();
        let mut poll = 0;
        assert_eq!(
            wake.refine_streamwise(0.3, 40, || {
                poll += 1;
                poll == stop
            })
            .unwrap_err()
            .code,
            "wake-refinement-cancelled"
        );
        assert_eq!(wake, original);
        assert_eq!(wake.rows[0].nodes.as_ptr(), pointer);
    }
}

#[test]
fn invalid_state_and_late_coordinate_resolution_failure_are_atomic() {
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut wake = warped();
        let original = wake.clone();
        assert_eq!(
            wake.refine_streamwise(target, 40, || false)
                .unwrap_err()
                .code,
            "wake-refinement-invalid"
        );
        assert_eq!(wake, original);
    }
    let mut bad = warped();
    bad.rows[1].age_ticks = -1.0;
    let original = bad.clone();
    assert_eq!(
        bad.refine_streamwise(0.3, 40, || false).unwrap_err().code,
        "wake-refinement-invalid"
    );
    assert_eq!(bad, original);
    let mut coarse_coordinates = straight(4.0, 16);
    let offset = 2.0_f64.powi(54);
    for p in coarse_coordinates
        .rows
        .iter_mut()
        .flat_map(|r| &mut r.nodes)
        .chain(&mut coarse_coordinates.line_nodes)
    {
        p[0] += offset;
    }
    let original = coarse_coordinates.clone();
    assert_eq!(
        coarse_coordinates
            .refine_streamwise(1.0, 10, || false)
            .unwrap_err()
            .code,
        "wake-refinement-unrepresentable"
    );
    assert_eq!(coarse_coordinates, original);
}

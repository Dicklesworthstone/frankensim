//! End-to-end error, work and material-state checks for accelerated wake motion.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{
    FilamentWake, SelfInducedShed, WakeBoundary, WakeRateCertificate, WakeTreeSpec,
};
use fs_vpm::images3d::ImagePlane;

fn sheet(rows: usize) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1000.0,
            n_stations: 4,
            max_rows: 512,
        },
        (0..=4).map(|i| [0.0, -1.0 + 0.5 * i as f64, 1.2]).collect(),
    )
    .unwrap();
    for r in 0..rows {
        let t = r as f64;
        let gamma: Vec<_> = (0..4)
            .map(|s| 0.12 + 0.03 * (0.17 * t + 0.6 * s as f64).sin())
            .collect();
        wake.shed(&gamma, [-0.15, 0.001, 0.0003]).unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = 0.04 + 0.00002 * t;
    }
    wake
}

fn config() -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [-0.4, 0.01, 0.0],
        emitted_core_radius_m: 0.2,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps: 2,
        max_displacement_fraction: 1.0,
        max_work: 100_000_000,
    }
}
fn spec(tolerance: f64) -> WakeTreeSpec {
    WakeTreeSpec {
        absolute_velocity_error_mps: tolerance,
        leaf_capacity: 8,
    }
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}
fn error(a: &FilamentWake, b: &FilamentWake) -> f64 {
    a.rows
        .iter()
        .zip(&b.rows)
        .flat_map(|(a, b)| a.nodes.iter().zip(&b.nodes))
        .map(|(&a, &b)| distance(a, b))
        .fold(0.0, f64::max)
}

#[test]
fn zero_tolerance_keeps_exact_source_order_and_material_state() {
    let mut exact = sheet(24);
    let mut tree = exact.clone();
    let boundary = WakeBoundary {
        line_endpoint_m: None,
        image_plane: Some(ImagePlane { z0_m: 0.0 }),
    };
    let incident = |p: [f64; 3], t: f64| Ok([0.01 * p[1], 0.02 * t, 0.0]);
    exact
        .shed_self_induced_with_environment(&[0.2; 4], config(), boundary, incident, || false)
        .unwrap();
    let report = tree
        .shed_self_induced_with_tree(&[0.2; 4], config(), boundary, spec(0.0), incident, || false)
        .unwrap();
    assert_eq!(tree, exact);
    assert_eq!(report.accepted_cells, 0);
    assert_eq!(report.max_velocity_error_bound_mps, 0.0);
    assert!(report.direct_interactions > 0);
    assert_eq!(report.peak_tree_nodes, 0);
}

#[test]
fn tree_reduces_exact_work_and_tight_accuracy_recovers_the_direct_trajectory() {
    let original = sheet(128);
    let mut exact = original.clone();
    exact
        .shed_self_induced(&[0.2; 4], config(), || false)
        .unwrap();
    let mut errors = Vec::new();
    for tolerance in [1e-3, 1e-5] {
        let mut wake = original.clone();
        let report = wake
            .shed_self_induced_with_tree(
                &[0.2; 4],
                config(),
                WakeBoundary::default(),
                spec(tolerance),
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap();
        let actual = error(&wake, &exact);
        assert!(
            actual <= 1.2 * tolerance / 1000.0,
            "local step error {actual:e}, tolerance {tolerance:e}, {report:?}"
        );
        assert!(
            report.accepted_cells > 0,
            "fixture must exercise actual acceleration: {report:?}"
        );
        assert!(report.max_velocity_error_bound_mps <= tolerance);
        let mut source_count = 0;
        original.visit_segments_with_core_growth(0.0, |_| source_count += 1);
        let direct_pairs = 4 * (original.rows.len() * 5) as u64 * source_count;
        assert!(
            report.direct_interactions < direct_pairs,
            "{report:?} >= {direct_pairs}"
        );
        for (after, before) in wake.rows.iter().zip(&original.rows) {
            assert_eq!(after.gamma, before.gamma);
            assert_eq!(after.core2_m2, before.core2_m2);
            assert_eq!(after.age_ticks, before.age_ticks + 1.0);
        }
        eprintln!(
            "tolerance={tolerance:e} position_error={actual:e} direct_pairs={direct_pairs} report={report:?}"
        );
        errors.push(actual);
    }
    assert!(errors[1] < errors[0], "{errors:?}");
}

#[test]
fn tree_moments_and_shared_vertices_commute_with_physical_partitioning() {
    let mut full = sheet(80);
    let mut near = full.clone();
    let mut far =
        FarField::aggregate(&mut near, 40, 8, WakeCoreEvolutionMode::Frozen, 0.001).unwrap();
    let endpoint: Vec<_> = full
        .line_nodes
        .iter()
        .map(|p| [p[0] + 0.0001, p[1], p[2] + 0.0001])
        .collect();
    let boundary = WakeBoundary {
        line_endpoint_m: Some(&endpoint),
        image_plane: Some(ImagePlane { z0_m: 0.0 }),
    };
    let incident = |p: [f64; 3], t: f64| Ok([0.01 * p[1], 0.02 * t, 0.0]);
    let a = full
        .shed_self_induced_with_tree(&[0.2; 4], config(), boundary, spec(1e-4), incident, || {
            false
        })
        .unwrap();
    let b = far
        .shed_self_induced_with_tree(
            &mut near,
            &[0.2; 4],
            config(),
            boundary,
            spec(1e-4),
            incident,
            || false,
        )
        .unwrap();
    let expected_far =
        FarField::aggregate(&mut full, 40, 8, WakeCoreEvolutionMode::Frozen, 0.001).unwrap();
    assert_eq!(near, full);
    assert_eq!(far, expected_far);
    assert_eq!(a.accepted_cells, b.accepted_cells);
    assert_eq!(
        a.max_velocity_error_bound_mps,
        b.max_velocity_error_bound_mps
    );
    assert!(a.accepted_cells > 0);
}

#[test]
fn growing_cores_and_images_remain_within_the_declared_short_step_accuracy() {
    let original = sheet(64);
    let mut cfg = config();
    cfg.core_mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.03 };
    let boundary = WakeBoundary {
        line_endpoint_m: None,
        image_plane: Some(ImagePlane { z0_m: 0.0 }),
    };
    let mut direct = original.clone();
    let mut tree = original;
    direct
        .shed_self_induced_with_environment(&[0.2; 4], cfg, boundary, |_, _| Ok([0.0; 3]), || false)
        .unwrap();
    let report = tree
        .shed_self_induced_with_tree(
            &[0.2; 4],
            cfg,
            boundary,
            spec(1e-4),
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    assert!(
        error(&direct, &tree) <= 1.2e-7,
        "{:e} {report:?}",
        error(&direct, &tree)
    );
    assert!(report.accepted_cells > 0);
    for (a, b) in tree.rows.iter().zip(&direct.rows) {
        assert_eq!(a.age_ticks, b.age_ticks);
        assert_eq!(a.core2_m2, b.core2_m2);
    }
}

#[test]
fn measured_work_cap_and_mid_tree_cancellation_roll_back_the_entire_step() {
    let original = sheet(48);
    let mut count = 0;
    let mut reference = original.clone();
    let report = reference
        .shed_self_induced_with_tree(
            &[0.2; 4],
            config(),
            WakeBoundary::default(),
            spec(1e-3),
            |_, _| Ok([0.0; 3]),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    assert!(report.work_used > 1000);
    let mut cfg = config();
    cfg.max_work = report.work_used;
    let mut at_cap = original.clone();
    assert_eq!(
        at_cap
            .shed_self_induced_with_tree(
                &[0.2; 4],
                cfg,
                WakeBoundary::default(),
                spec(1e-3),
                |_, _| Ok([0.0; 3]),
                || false
            )
            .unwrap(),
        report
    );
    assert_eq!(at_cap, reference);
    cfg.max_work -= 1;
    let mut failed = original.clone();
    assert_eq!(
        failed
            .shed_self_induced_with_tree(
                &[0.2; 4],
                cfg,
                WakeBoundary::default(),
                spec(1e-3),
                |_, _| Ok([0.0; 3]),
                || false
            )
            .unwrap_err()
            .code,
        "wake-advection-work-exhausted"
    );
    assert_eq!(failed, original);
    let mut cancelled = original.clone();
    let mut polls = 0;
    assert_eq!(
        cancelled
            .shed_self_induced_with_tree(
                &[0.2; 4],
                config(),
                WakeBoundary::default(),
                spec(1e-3),
                |_, _| Ok([0.0; 3]),
                || {
                    polls += 1;
                    polls >= count / 2
                }
            )
            .unwrap_err()
            .code,
        "wake-advection-cancelled"
    );
    assert_eq!(cancelled, original);
}

#[test]
fn invalid_accuracy_or_leaf_capacity_refuses_before_publication() {
    let original = sheet(8);
    for tree in [
        spec(f64::NAN),
        spec(-1.0),
        WakeTreeSpec {
            leaf_capacity: 0,
            ..spec(1e-4)
        },
        WakeTreeSpec {
            leaf_capacity: 257,
            ..spec(1e-4)
        },
    ] {
        let mut wake = original.clone();
        assert_eq!(
            wake.shed_self_induced_with_tree(
                &[0.2; 4],
                config(),
                WakeBoundary::default(),
                tree,
                |_, _| Ok([0.0; 3]),
                || false
            )
            .unwrap_err()
            .code,
            "wake-tree-invalid"
        );
        assert_eq!(wake, original);
    }
}

//! Adaptive motion of actual retained far sources and attached near interfaces.
use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs};
use fs_vpm::filament3d::{
    FilamentWake, SelfInducedShed, WakeBoundary, WakeRateCertificate, WakeTimeSpec, WakeTreeSpec,
};
use fs_vpm::images3d::ImagePlane;

const CORE: f64 = 0.2;

fn settings(mode: WakeCoreEvolutionMode) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [-0.3, 0.01, 0.02],
        emitted_core_radius_m: CORE,
        core_mode: mode,
        substeps: 1,
        max_displacement_fraction: 1.0,
        max_work: 100_000_000,
    }
}

fn time(tolerance: f64) -> WakeTimeSpec {
    WakeTimeSpec {
        absolute_position_tolerance_m: tolerance,
        relative_position_tolerance: 0.0,
        min_step_s: 1e-10,
        max_trials: 8192,
    }
}

fn sheet(rows: usize, dt: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1.0 / dt,
            n_stations: 3,
            max_rows: 256,
        },
        [-0.9, -0.3, 0.3, 0.9]
            .map(|y| [0.0, y, 1.0 + 0.04 * y * y])
            .to_vec(),
    )
    .unwrap();
    for r in 0..rows {
        let t = r as f64;
        wake.shed(
            &[0.25 + 0.04 * t, 0.45 + 0.02 * t, 0.18 - 0.01 * t],
            [-0.3, 0.0, 0.0],
        )
        .unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = CORE * CORE + 0.003 * t;
    }
    for (r, row) in wake.rows.iter_mut().take(rows - 1).enumerate() {
        for (s, p) in row.nodes.iter_mut().enumerate() {
            p[2] += 0.001 * (r + 1) as f64 * (s as f64 - 1.5);
        }
    }
    wake
}

fn partition(full: &FilamentWake, mode: WakeCoreEvolutionMode) -> (FilamentWake, FarField) {
    let mut near = full.clone();
    let dt = 1.0 / near.cert.shed_hz;
    let far = FarField::aggregate(&mut near, 3, 1, mode, dt).unwrap();
    (near, far)
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    let distance = (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
    assert!(distance <= tolerance, "{a:?} != {b:?}, error {distance:e}");
}

fn compare_geometry(
    near: &FilamentWake,
    far: &FarField,
    expected_near: &FilamentWake,
    expected_far: &FarField,
    tolerance: f64,
) {
    assert_eq!(near.rows.len(), expected_near.rows.len());
    assert_eq!(far.cells.len(), expected_far.cells.len());
    for (row, expected) in near.rows.iter().zip(&expected_near.rows) {
        assert_eq!(row.gamma, expected.gamma);
        assert_eq!(row.age_ticks, expected.age_ticks);
        assert_eq!(row.core2_m2, expected.core2_m2);
        for (&p, &q) in row.nodes.iter().zip(&expected.nodes) {
            assert_close(p, q, tolerance);
        }
    }
    for (cell, expected) in far.cells.iter().zip(&expected_far.cells) {
        assert_eq!(cell.segs.len(), expected.segs.len());
        for ((a, b, gamma), (ea, eb, egamma)) in cell.segs.iter().zip(&expected.segs) {
            assert_close(*a, *ea, tolerance);
            assert_close(*b, *eb, tolerance);
            assert_eq!(gamma, egamma);
        }
        assert_close(cell.centroid, expected.centroid, tolerance);
        assert_close(cell.w_sum, expected.w_sum, 20.0 * tolerance);
    }
}

fn free_square(mode: WakeCoreEvolutionMode) -> (FilamentWake, FarField) {
    let near = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 5.0,
            n_stations: 1,
            max_rows: 8,
        },
        vec![[3.0, -0.5, 1.0], [3.0, 0.5, 1.0]],
    )
    .unwrap();
    let vertices = [
        [-0.5, -0.5, 1.0],
        [0.5, -0.5, 1.0],
        [0.5, 0.5, 1.0],
        [-0.5, 0.5, 1.0],
    ];
    let cell = cell_from_segs(
        (0..4)
            .map(|i| (vertices[i], vertices[(i + 1) % 4], 2.0))
            .collect(),
        CORE * CORE,
    );
    (
        near,
        FarField {
            mode,
            shed_dt_s: 0.2,
            cells: vec![cell],
        },
    )
}

#[test]
fn free_ring_translation_matches_finite_core_speed_and_its_diffusion_integral() {
    let frozen = WakeCoreEvolutionMode::Frozen;
    let (mut near, mut far) = free_square(frozen);
    let original = far.clone();
    let mut config = settings(frozen);
    config.ambient_velocity_mps = [0.0; 3];
    let report = far
        .shed_self_induced_adaptive(
            &mut near,
            &[0.0],
            config,
            time(1e-13),
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    // At a corner only the two opposite unit edges contribute. Their exact
    // integrated speed is 2*Gamma/[4*pi*(1+a²)*sqrt(2+a²)].
    let speed =
        4.0 / (4.0 * std::f64::consts::PI * (1.0 + CORE * CORE) * (2.0 + CORE * CORE).sqrt());
    assert_eq!(report.accepted_substeps, 1);
    assert_eq!(report.rejected_trials, 0);
    for ((a, b, gamma), (old_a, old_b, old_gamma)) in
        far.cells[0].segs.iter().zip(&original.cells[0].segs)
    {
        assert_close(*a, [old_a[0], old_a[1], old_a[2] + 0.2 * speed], 2e-15);
        assert_close(*b, [old_b[0], old_b[1], old_b[2] + 0.2 * speed], 2e-15);
        assert_eq!(gamma, old_gamma);
    }
    assert_eq!(far.cells[0].w_sum, [0.0; 3]);
    assert_eq!(near.rows.len(), 1);
    assert_eq!(near.rows[0].age_ticks, 0.0);

    let nu = 0.1;
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: nu };
    // Integrate the changing speed with a²(t)=a²(0)+4*nu*t.
    let primitive = |core2: f64| {
        let u = (2.0 + core2).sqrt();
        2.0 / (8.0 * std::f64::consts::PI * nu) * ((u - 1.0) / (u + 1.0)).ln()
    };
    let expected = primitive(CORE * CORE + 4.0 * nu * 0.2) - primitive(CORE * CORE);
    let mut errors = Vec::new();
    for tolerance in [1e-6, 1e-10, 1e-13] {
        let (mut near, mut far) = free_square(mode);
        let mut config = settings(mode);
        config.ambient_velocity_mps = [0.0; 3];
        let report = far
            .shed_self_induced_adaptive(
                &mut near,
                &[0.0],
                config,
                time(tolerance),
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap();
        errors.push((far.cells[0].segs[0].0[2] - 1.0 - expected).abs());
        assert_eq!(far.cells[0].source_core2_m2, vec![CORE * CORE; 4]);
        assert_eq!(far.cells[0].diffusion_offset_m2, 4.0 * nu * 0.2);
        assert!(report.max_accepted_error_ratio <= 1.0);
        for &(_, _, gamma) in &far.cells[0].segs {
            assert_eq!(gamma, 2.0);
        }
    }
    assert!(
        errors[1] < errors[0] / 100.0 && errors[2] < errors[1] / 20.0,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-11, "{errors:?}");
}

#[test]
fn near_far_partition_shares_adaptive_stages_line_images_and_elapsed_cores() {
    for mode in [
        WakeCoreEvolutionMode::Frozen,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.02 },
    ] {
        let mut full = sheet(5, 0.05);
        let (mut near, mut far) = partition(&full, mode);
        let old_far = far.clone();
        let endpoint: Vec<_> = full
            .line_nodes
            .iter()
            .map(|p| [p[0] + 0.01, p[1], p[2] + 0.005])
            .collect();
        let boundary = WakeBoundary {
            line_endpoint_m: Some(&endpoint),
            image_plane: Some(ImagePlane { z0_m: 0.0 }),
        };
        let variation =
            |p: [f64; 3], t: f64| Ok([0.02 * p[1], 0.04 * p[2], -0.01 * p[0] + 0.03 * t]);
        let control = time(1e-10);
        let full_report = full
            .shed_self_induced_adaptive(
                &[0.5, 0.3, 0.2],
                settings(mode),
                control,
                boundary,
                None,
                variation,
                || false,
            )
            .unwrap();
        let split_report = far
            .shed_self_induced_adaptive(
                &mut near,
                &[0.5, 0.3, 0.2],
                settings(mode),
                control,
                boundary,
                None,
                variation,
                || false,
            )
            .unwrap();
        assert!(full_report.accepted_substeps > 1);
        assert_eq!(
            full_report.accepted_substeps,
            split_report.accepted_substeps
        );
        assert_eq!(full_report.rejected_trials, split_report.rejected_trials);
        let (expected_near, expected_far) = partition(&full, mode);
        compare_geometry(&near, &far, &expected_near, &expected_far, 3e-11);
        assert_eq!(near.line_nodes, endpoint);
        assert_eq!(near.rows.last().unwrap().nodes, endpoint);
        assert_eq!(near.rows.last().unwrap().age_ticks, 0.0);
        for ((cell, old), expected_cell) in far
            .cells
            .iter()
            .zip(&old_far.cells)
            .zip(&expected_far.cells)
        {
            assert_eq!(cell.source_core2_m2, old.source_core2_m2);
            let growth = match mode {
                WakeCoreEvolutionMode::Frozen => 0.0,
                WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => 4.0 * nu_eff_m2ps * 0.05,
            };
            assert_eq!(cell.diffusion_offset_m2, old.diffusion_offset_m2 + growth);
            for p in [[-0.7, 0.2, 1.3], [1.0, 2.0, 4.0]] {
                assert_close(cell.exact(p), expected_cell.exact(p), 5e-11);
            }
        }
    }
}

#[test]
fn adaptive_tree_uses_actual_far_sources_and_accounts_for_all_work() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let original = sheet(48, 0.001);
    let (initial_near, initial_far) = partition(&original, mode);
    let control = time(1e-7);
    let mut direct = original;
    let direct_report = direct
        .shed_self_induced_adaptive(
            &[0.5, 0.3, 0.2],
            settings(mode),
            control,
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    let tree = WakeTreeSpec {
        absolute_velocity_error_mps: 1e-3,
        leaf_capacity: 8,
    };
    let mut near = initial_near.clone();
    let mut far = initial_far.clone();
    let report = far
        .shed_self_induced_adaptive(
            &mut near,
            &[0.5, 0.3, 0.2],
            settings(mode),
            control,
            WakeBoundary::default(),
            Some(tree),
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    assert!(report.advection.accepted_cells > 0, "{report:?}");
    assert!(
        report.advection.direct_interactions < direct_report.advection.direct_interactions,
        "{report:?}"
    );
    assert!(report.advection.max_velocity_error_bound_mps <= tree.absolute_velocity_error_mps);
    let (expected_near, expected_far) = partition(&direct, mode);
    compare_geometry(&near, &far, &expected_near, &expected_far, 2.1e-6);
    let mut config = settings(mode);
    config.max_work = report.advection.work_used;
    let mut same_near = initial_near.clone();
    let mut same_far = initial_far.clone();
    let repeated = same_far
        .shed_self_induced_adaptive(
            &mut same_near,
            &[0.5, 0.3, 0.2],
            config,
            control,
            WakeBoundary::default(),
            Some(tree),
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    assert_eq!(repeated, report);
    assert_eq!(same_near, near);
    assert_eq!(same_far, far);
    config.max_work -= 1;
    let mut failed_near = initial_near.clone();
    let mut failed_far = initial_far.clone();
    assert_eq!(
        failed_far
            .shed_self_induced_adaptive(
                &mut failed_near,
                &[0.5, 0.3, 0.2],
                config,
                control,
                WakeBoundary::default(),
                Some(tree),
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err()
            .code,
        "wake-advection-work-exhausted"
    );
    assert_eq!(failed_near, initial_near);
    assert_eq!(failed_far, initial_far);
}

#[test]
fn cancellation_after_completed_stages_and_far_rebuild_rolls_back_both_partitions() {
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.02 };
    let original = sheet(5, 0.05);
    let (initial_near, initial_far) = partition(&original, mode);
    let mut near = initial_near.clone();
    let mut far = initial_far.clone();
    let mut full_checks = 0;
    let report = far
        .shed_self_induced_adaptive(
            &mut near,
            &[0.5, 0.3, 0.2],
            settings(mode),
            time(1e-10),
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || {
                full_checks += 1;
                false
            },
        )
        .unwrap();
    assert!(report.accepted_substeps > 1 && report.rejected_trials > 0);
    for stop in [full_checks / 2, full_checks] {
        let mut near = initial_near.clone();
        let mut far = initial_far.clone();
        let mut checks = 0;
        assert_eq!(
            far.shed_self_induced_adaptive(
                &mut near,
                &[0.5, 0.3, 0.2],
                settings(mode),
                time(1e-10),
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || {
                    checks += 1;
                    checks == stop
                },
            )
            .unwrap_err()
            .code,
            "wake-advection-cancelled"
        );
        assert_eq!(near, initial_near);
        assert_eq!(far, initial_far);
    }
}

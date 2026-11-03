//! G2 analytic far-ring motion, G3 partition equivalence, G4 atomic failures.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs};
use fs_vpm::filament3d::{FilamentWake, SelfInducedShed, WakeBoundary, WakeRateCertificate};
use fs_vpm::images3d::ImagePlane;

const DT: f64 = 0.05;
const CORE: f64 = 0.2;

fn sheet() -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1.0 / DT,
            n_stations: 3,
            max_rows: 64,
        },
        [-0.9, -0.3, 0.3, 0.9]
            .map(|y| [0.0, y, 1.0 + 0.04 * y * y])
            .to_vec(),
    )
    .unwrap();
    for r in 0..5 {
        let t = r as f64;
        wake.shed(
            &[0.25 + 0.04 * t, 0.45 + 0.02 * t, 0.18 - 0.01 * t],
            [-0.3, 0.0, 0.0],
        )
        .unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = CORE * CORE + 0.003 * t;
    }
    for (r, row) in wake.rows.iter_mut().take(4).enumerate() {
        for (s, p) in row.nodes.iter_mut().enumerate() {
            p[2] += 0.01 * (r + 1) as f64 * (s as f64 - 1.5);
        }
    }
    wake
}

fn settings(mode: WakeCoreEvolutionMode) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [-0.3, 0.01, 0.02],
        emitted_core_radius_m: CORE,
        core_mode: mode,
        substeps: 8,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    }
}

fn partition(full: &FilamentWake, mode: WakeCoreEvolutionMode) -> (FilamentWake, FarField) {
    let mut near = full.clone();
    let far = FarField::aggregate(&mut near, 3, 1, mode, DT).unwrap();
    (near, far)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(distance(a, b) <= tolerance, "{a:?} != {b:?}");
}

fn exact_hybrid(near: &FilamentWake, far: &FarField, p: [f64; 3]) -> [f64; 3] {
    let growth = match far.mode {
        WakeCoreEvolutionMode::Frozen => 0.0,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps } => 4.0 * nu_eff_m2ps * DT,
    };
    let mut v = near.induced_velocity_with_core_growth(p, growth);
    for cell in &far.cells {
        let far_v = cell.exact(p);
        for c in 0..3 {
            v[c] += far_v[c];
        }
    }
    v
}

#[test]
fn frozen_partition_preserves_the_complete_motion_and_rebuilds_all_moments() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let mut full = sheet();
    let (mut near, mut far) = partition(&full, mode);
    let original = far.clone();
    for interval in 0..3 {
        let gamma = [0.5 + 0.01 * interval as f64, 0.3, 0.2];
        full.shed_self_induced(&gamma, settings(mode), || false)
            .unwrap();
        far.shed_self_induced(&mut near, &gamma, settings(mode), || false)
            .unwrap();
        let (expected_near, expected_far) = partition(&full, mode);
        // Exact equality includes the two sides of every interface, all
        // moments/radii and all source cores, rather than only probe fields.
        assert_eq!(near, expected_near);
        assert_eq!(far, expected_far);
        for (before, after) in original.cells.iter().zip(&far.cells) {
            assert_eq!(before.source_core2_m2, after.source_core2_m2);
            assert_eq!(
                before.segs.iter().map(|s| s.2).collect::<Vec<_>>(),
                after.segs.iter().map(|s| s.2).collect::<Vec<_>>()
            );
        }
    }
    assert_ne!(far.cells[0].dip, original.cells[0].dip);
    assert_ne!(far.cells[0].quad, original.cells[0].quad);
    assert_ne!(far.cells[0].radius_m, original.cells[0].radius_m);
    assert_eq!(near.rows.last().unwrap().nodes, near.line_nodes);
    assert_eq!(near.rows.last().unwrap().age_ticks, 0.0);
}

#[test]
fn spreading_keeps_source_cores_and_common_elapsed_time_across_rebuilds() {
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.02 };
    let growth = 4.0 * 0.02 * DT;
    let mut full = sheet();
    let (mut near, mut far) = partition(&full, mode);
    let original = far.clone();
    for interval in 1..=3 {
        full.shed_self_induced(&[0.4, 0.2, 0.1], settings(mode), || false)
            .unwrap();
        far.shed_self_induced(&mut near, &[0.4, 0.2, 0.1], settings(mode), || false)
            .unwrap();
        let (expected_near, expected_far) = partition(&full, mode);
        for (a, b) in near.rows.iter().zip(&expected_near.rows) {
            assert_eq!(a.gamma, b.gamma);
            assert_eq!(a.age_ticks, b.age_ticks);
            assert_eq!(a.core2_m2, b.core2_m2);
            for (&a, &b) in a.nodes.iter().zip(&b.nodes) {
                assert_close(a, b, 3e-14);
            }
        }
        for ((initial, actual), expected) in original
            .cells
            .iter()
            .zip(&far.cells)
            .zip(&expected_far.cells)
        {
            assert_eq!(actual.source_core2_m2, initial.source_core2_m2);
            assert!((actual.diffusion_offset_m2 - interval as f64 * growth).abs() < 1e-17);
            for ((a, b, gamma), (ea, eb, egamma)) in actual.segs.iter().zip(&expected.segs) {
                assert_close(*a, *ea, 3e-14);
                assert_close(*b, *eb, 3e-14);
                assert_eq!(gamma, egamma);
            }
            assert!((actual.core2_m2 - expected.core2_m2).abs() < 3e-16);
            assert_close(actual.core_moment, expected.core_moment, 3e-16);
            for p in [[-0.7, 0.2, 1.3], [1.0, 2.0, 4.0]] {
                assert_close(actual.eval(p), expected.eval(p), 3e-14);
            }
        }
        for p in [[-0.7, 0.2, 1.3], [1.0, 2.0, 4.0]] {
            assert_close(
                exact_hybrid(&near, &far, p),
                full.induced_velocity_with_core_growth(p, growth),
                3e-14,
            );
        }
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
    let far = FarField {
        mode,
        shed_dt_s: 0.2,
        cells: vec![cell],
    };
    (near, far)
}

#[test]
fn retained_free_square_ring_translates_at_the_analytic_finite_core_speed() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let (mut near, mut far) = free_square(mode);
    let original = far.clone();
    let mut config = settings(mode);
    config.ambient_velocity_mps = [0.0; 3];
    // At a corner, incident straight edges contribute zero. Integrating
    // each opposite unit edge gives Γ / [4π(1+a²)sqrt(2+a²)].
    let speed =
        2.0 * 2.0 / (4.0 * std::f64::consts::PI * (1.0 + CORE * CORE) * (2.0 + CORE * CORE).sqrt());
    far.shed_self_induced(&mut near, &[0.0], config, || false)
        .unwrap();
    for ((a, b, gamma), (old_a, old_b, old_gamma)) in
        far.cells[0].segs.iter().zip(&original.cells[0].segs)
    {
        assert_close(*a, [old_a[0], old_a[1], old_a[2] + 0.2 * speed], 2e-15);
        assert_close(*b, [old_b[0], old_b[1], old_b[2] + 0.2 * speed], 2e-15);
        assert_eq!(gamma, old_gamma);
    }
    assert_eq!(far.cells[0].w_sum, [0.0; 3]);
    assert_eq!(far.cells[0].radius_m, original.cells[0].radius_m);
    // This symmetric free ring transports its vector impulse without a
    // support; no general anchored-wake conservation claim follows.
    let mut impulse = [0.0; 3];
    for &(a, b, g) in &far.cells[0].segs {
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        for c in 0..3 {
            impulse[c] += 0.5 * g * cross[c];
        }
    }
    assert_close(impulse, [0.0, 0.0, 2.0], 2e-15);
}

#[test]
fn far_core_growth_uses_midpoint_age_with_second_order_time_convergence() {
    let nu = 0.1;
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: nu };
    // Integrate the square's exact speed with a²(t)=a²(0)+4νt.
    // Substitution u=sqrt(2+a²) gives Γ/(8πν) ln((u-1)/(u+1)).
    let primitive = |core2: f64| {
        let u = (2.0 + core2).sqrt();
        2.0 / (8.0 * std::f64::consts::PI * nu) * ((u - 1.0) / (u + 1.0)).ln()
    };
    let exact = primitive(CORE * CORE + 4.0 * nu * 0.2) - primitive(CORE * CORE);
    let errors: Vec<_> = [2, 4, 8]
        .into_iter()
        .map(|substeps| {
            let (mut near, mut far) = free_square(mode);
            let mut config = settings(mode);
            config.ambient_velocity_mps = [0.0; 3];
            config.substeps = substeps;
            far.shed_self_induced(&mut near, &[0.0], config, || false)
                .unwrap();
            assert_eq!(far.cells[0].source_core2_m2, vec![CORE * CORE; 4]);
            assert_eq!(far.cells[0].diffusion_offset_m2, 4.0 * nu * 0.2);
            (far.cells[0].segs[0].0[2] - 1.0 - exact).abs()
        })
        .collect();
    for pair in errors.windows(2) {
        assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn moving_line_incident_shear_and_images_commute_with_partition() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let mut full = sheet();
    let (mut near, mut far) = partition(&full, mode);
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
        |p: [f64; 3], elapsed: f64| Ok([0.02 * p[1], 0.04 * p[2], -0.01 * p[0] + 0.03 * elapsed]);
    full.shed_self_induced_with_environment(
        &[0.5, 0.3, 0.2],
        settings(mode),
        boundary,
        variation,
        || false,
    )
    .unwrap();
    far.shed_self_induced_with_environment(
        &mut near,
        &[0.5, 0.3, 0.2],
        settings(mode),
        boundary,
        variation,
        || false,
    )
    .unwrap();
    let (expected_near, expected_far) = partition(&full, mode);
    assert_eq!(near, expected_near);
    assert_eq!(far, expected_far);
    assert_eq!(near.line_nodes, endpoint);
}

#[test]
fn mismatched_schedule_open_sources_bad_cores_and_work_caps_refuse_atomically() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let (near, far) = partition(&sheet(), mode);
    for case in 0..6 {
        let mut n = near.clone();
        let mut f = far.clone();
        let mut config = settings(mode);
        let expected = match case {
            0 => {
                f.shed_dt_s *= 2.0;
                "wake-advection-schedule-mismatch"
            }
            1 => {
                config.core_mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.01 };
                "wake-advection-schedule-mismatch"
            }
            2 => {
                f.cells[0].segs[0].0[0] -= 0.01;
                "wake-advection-topology-invalid"
            }
            3 => {
                f.cells.remove(0);
                "wake-advection-topology-invalid"
            }
            4 => {
                f.cells[0].source_core2_m2[0] = -1.0;
                "wake-advection-core-invalid"
            }
            _ => {
                config.max_work = 1_000;
                "wake-advection-work-exhausted"
            }
        };
        let n_before = n.clone();
        let f_before = f.clone();
        let error = f
            .shed_self_induced(&mut n, &[0.5, 0.3, 0.2], config, || false)
            .unwrap_err();
        assert_eq!(error.code, expected, "case {case}");
        assert_eq!(n, n_before);
        assert_eq!(f, f_before);
    }
}

#[test]
fn early_midstep_and_final_cancellation_publish_neither_partition() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let (near, far) = partition(&sheet(), mode);
    let mut polls = 0;
    far.clone()
        .shed_self_induced(&mut near.clone(), &[0.5, 0.3, 0.2], settings(mode), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 30);
    for stop in [1, polls / 2, polls] {
        let mut n = near.clone();
        let mut f = far.clone();
        let mut observed = 0;
        let error = f
            .shed_self_induced(&mut n, &[0.5, 0.3, 0.2], settings(mode), || {
                observed += 1;
                observed >= stop
            })
            .unwrap_err();
        assert_eq!(error.code, "wake-advection-cancelled");
        assert_eq!(n, near);
        assert_eq!(f, far);
    }
}

#[test]
fn image_plane_admission_includes_far_endpoints() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let (mut near, mut far) = free_square(mode);
    far.translate([0.0, 0.0, -1.1]).unwrap();
    let n_before = near.clone();
    let f_before = far.clone();
    let mut config = settings(mode);
    config.ambient_velocity_mps = [0.0; 3];
    let error = far
        .shed_self_induced_with_environment(
            &mut near,
            &[0.0],
            config,
            WakeBoundary {
                line_endpoint_m: None,
                image_plane: Some(ImagePlane { z0_m: 0.0 }),
            },
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap_err();
    assert_eq!(error.code, "image-plane-invalid");
    assert_eq!(near, n_before);
    assert_eq!(far, f_before);
}

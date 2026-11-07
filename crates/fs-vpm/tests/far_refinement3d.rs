//! G0/G2/G3: exact source invariants and material bending after subdivision.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs, hybrid_kinematics};
use fs_vpm::filament3d::{FilamentWake, SelfInducedShed, WakeBoundary, WakeRateCertificate};

fn close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    let error = (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
    assert!(error <= tolerance, "{a:?} != {b:?}: error {error}");
}

fn sheet(gamma: f64, spacing: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 1,
            max_rows: 32,
        },
        vec![[0.0, -0.5, 2.0], [0.0, 0.5, 2.0]],
    )
    .unwrap();
    for _ in 0..6 {
        wake.shed(&[gamma], [-spacing, 0.0, 0.0]).unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = 0.04;
    }
    wake
}

fn single(a: [f64; 3], b: [f64; 3], core: f64) -> FarField {
    FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        shed_dt_s: 0.01,
        cells: vec![cell_from_segs(vec![(a, b, 0.75)], core)],
    }
}

#[test]
fn exact_velocity_curl_and_both_impulses_survive_with_diffusion_history() {
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.02 };
    let mut near = sheet(0.5, 0.8);
    for (i, row) in near.rows.iter_mut().enumerate() {
        row.gamma[0] = [0.6, -0.2, 0.75, 0.25, -0.4, 0.5][i];
        row.core2_m2 += 0.01 * i as f64;
        for point in &mut row.nodes {
            point[2] += 0.03 * point[0] * point[1];
        }
    }
    let mut far = FarField::aggregate(&mut near, 4, 2, mode, 0.01).unwrap();
    far.advance_core(0.4).unwrap();
    let before = far.clone();
    let near_before = near.clone();
    let impulse = before
        .hydrodynamic_impulse_per_density(&near, [0.4, -0.2, 0.3], 10_000, || false)
        .unwrap();
    let report = far.refine_segments(0.15, 10_000, || false).unwrap();
    assert_eq!(report.cells_refined, 2);
    assert!(report.segments_after > report.segments_before);
    assert!(report.longest_segment_after_m <= 0.15 * (1.0 + 64.0 * f64::EPSILON));
    assert_eq!(near, near_before);
    assert_eq!(far.mode, before.mode);
    assert_eq!(far.shed_dt_s, before.shed_dt_s);
    for (old, new) in before.cells.iter().zip(&far.cells) {
        assert_eq!(
            old.diffusion_offset_m2.to_bits(),
            new.diffusion_offset_m2.to_bits()
        );
        // The original endpoints and their exact coordinate identities survive.
        for &(a, b, _) in &old.segs {
            assert!(new.segs.iter().any(|&(start, _, _)| start == a));
            assert!(new.segs.iter().any(|&(_, end, _)| end == b));
        }
        assert!(
            new.source_core2_m2
                .iter()
                .all(|core| old.source_core2_m2.contains(core))
        );
    }
    let actual = far
        .hydrodynamic_impulse_per_density(&near, [0.4, -0.2, 0.3], 10_000, || false)
        .unwrap();
    close(
        actual.linear_per_density_m4ps(),
        impulse.linear_per_density_m4ps(),
        2e-12,
    );
    close(
        actual.angular_per_density_m5ps(),
        impulse.angular_per_density_m5ps(),
        2e-11,
    );
    for point in [[-2.1, 0.2, 2.4], [0.1, -0.6, 2.0], [30.0, -40.0, 20.0]] {
        let expected = hybrid_kinematics(&near, &before, point).unwrap();
        let actual = hybrid_kinematics(&near, &far, point).unwrap();
        close(actual.velocity, expected.velocity, 3e-13);
        close(actual.vorticity(), expected.vorticity(), 3e-12);
        for axis in 0..3 {
            close(actual.gradient[axis], expected.gradient[axis], 3e-12);
        }
        for cell in &far.cells {
            let error = cell.expansion_error_bound(point);
            if error.is_finite() {
                close(cell.eval(point), cell.exact(point), error + 1e-14);
            }
        }
    }
    // Continued diffusion must use the same physical source cores and age.
    let mut coarse = before;
    coarse.advance_core(0.3).unwrap();
    far.advance_core(0.3).unwrap();
    close(
        far.exact_kinematics([-2.0, 0.3, 2.2]).unwrap().vorticity(),
        coarse
            .exact_kinematics([-2.0, 0.3, 2.2])
            .unwrap()
            .vorticity(),
        3e-12,
    );
}

#[test]
fn inserted_vertices_follow_native_coupled_motion_and_resolve_subsequent_bending() {
    // The divergence-free incident field u_z=0.05*x² has exact trajectories
    // z(t)=z(0)+0.05*x²*t. Tiny NONZERO circulation keeps induction negligible
    // while exercising the original vortex source and Kelvin-closure paths.
    let mut near = sheet(1e-10, 4.0);
    let mut far =
        FarField::aggregate(&mut near, 4, 2, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    let mut coarse_near = near.clone();
    let mut coarse_far = far.clone();
    far.refine_segments(1.0, 1_000, || false).unwrap();
    let settings = SelfInducedShed {
        ambient_velocity_mps: [0.0; 3],
        emitted_core_radius_m: 0.2,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps: 16,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    };
    for (near, far) in [(&mut near, &mut far), (&mut coarse_near, &mut coarse_far)] {
        far.shed_self_induced_with_environment(
            near,
            &[1e-10],
            settings,
            WakeBoundary::default(),
            |point, _| Ok([0.0, 0.0, 0.05 * point[0] * point[0]]),
            || false,
        )
        .unwrap();
    }
    let target_x = -18.0;
    let refined = far
        .cells
        .iter()
        .flat_map(|c| &c.segs)
        .flat_map(|&(a, b, _)| [a, b])
        .find(|p| (p[0] - target_x).abs() < 1e-8 && (p[1] + 0.5).abs() < 1e-8)
        .expect("refinement inserted this material point");
    let expected_z = 2.0 + 0.05 * target_x * target_x * 0.01;
    assert!((refined[2] - expected_z).abs() < 1e-8);
    let &(a, b, _) = coarse_far
        .cells
        .iter()
        .flat_map(|c| &c.segs)
        .find(|&&(a, b, _)| a[0] < target_x && target_x < b[0] && (a[1] + 0.5).abs() < 1e-8)
        .unwrap();
    let fraction = (target_x - a[0]) / (b[0] - a[0]);
    let chord_z = a[2] + fraction * (b[2] - a[2]);
    assert!((chord_z - expected_z).abs() > 0.001);
    far.hydrodynamic_impulse_per_density(&near, [0.0; 3], 10_000, || false)
        .unwrap();
    assert_eq!(near.rows.len(), coarse_near.rows.len());
    assert_eq!(near.ticks, coarse_near.ticks);
}

#[test]
fn complete_capacity_noop_and_every_cancellation_boundary_are_transactional() {
    let original = single([0.0; 3], [4.0, 0.0, 0.0], 0.04);
    let mut capped = original.clone();
    assert_eq!(
        capped.refine_segments(1.0, 3, || false).unwrap_err().code,
        "far-refinement-budget"
    );
    assert_eq!(capped, original);
    let mut refined = original.clone();
    let mut calls = 0;
    let report = refined
        .refine_segments(1.0, 4, || {
            calls += 1;
            false
        })
        .unwrap();
    assert_eq!((report.segments_before, report.segments_after), (1, 4));
    assert_eq!(report.longest_segment_after_m, 1.0);
    for stop in 1..=calls {
        let mut candidate = original.clone();
        let mut count = 0;
        assert_eq!(
            candidate
                .refine_segments(1.0, 4, || {
                    count += 1;
                    count == stop
                })
                .unwrap_err()
                .code,
            "far-refinement-cancelled"
        );
        assert_eq!(candidate, original, "cancellation at checkpoint {stop}");
    }
    let saved = refined.clone();
    let pointers = (
        refined.cells[0].segs.as_ptr(),
        refined.cells[0].source_core2_m2.as_ptr(),
    );
    assert_eq!(
        refined
            .refine_segments(1.0, 4, || false)
            .unwrap()
            .cells_refined,
        0
    );
    assert_eq!(refined, saved);
    assert_eq!(
        (
            refined.cells[0].segs.as_ptr(),
            refined.cells[0].source_core2_m2.as_ptr()
        ),
        pointers
    );
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut candidate = original.clone();
        assert!(candidate.refine_segments(target, 4, || false).is_err());
        assert_eq!(candidate, original);
    }
    let mut malformed = original.clone();
    malformed.cells[0].source_core2_m2.clear();
    let snapshot = malformed.clone();
    assert_eq!(
        malformed
            .refine_segments(1.0, 4, || false)
            .unwrap_err()
            .code,
        "far-refinement-invalid"
    );
    assert_eq!(malformed, snapshot);
}

#[test]
fn dimensional_scaling_and_coordinate_resolution_are_respected() {
    for scale in [1e-9, 1.0, 1e9] {
        let mut field = single([0.0; 3], [4.0 * scale, 0.0, 0.0], 0.04 * scale * scale);
        let original = field.clone();
        let report = field.refine_segments(scale, 4, || false).unwrap();
        assert_eq!(report.segments_after, 4);
        let point = [2.0 * scale, 0.3 * scale, -0.4 * scale];
        close(
            field
                .exact_kinematics(point)
                .unwrap()
                .velocity
                .map(|v| v * scale),
            original
                .exact_kinematics(point)
                .unwrap()
                .velocity
                .map(|v| v * scale),
            1e-13,
        );
    }
    let original = single([1e16, 0.0, 0.0], [1e16 + 4.0, 0.0, 0.0], 0.04);
    let mut candidate = original.clone();
    assert_eq!(
        candidate
            .refine_segments(1.0, 4, || false)
            .unwrap_err()
            .code,
        "far-refinement-unrepresentable"
    );
    assert_eq!(candidate, original);
    let mut candidate = single([0.0; 3], [4.0, 0.0, 0.0], 0.04);
    assert_eq!(
        candidate
            .refine_segments(1e-300, usize::MAX, || false)
            .unwrap_err()
            .code,
        "far-refinement-budget"
    );
}

#[test]
fn rounded_diagonal_zigzags_refuse_in_both_partitions_but_collinear_rounding_admits() {
    // All child lengths fit, but these rounded coordinates change the line.
    let diagonal = single([1e16, 0.0, 0.0], [1e16 + 4.0, 4.0, 0.0], 0.04);
    let mut candidate = diagonal.clone();
    assert_eq!(
        candidate
            .refine_segments(2.5, 3, || false)
            .unwrap_err()
            .code,
        "far-refinement-unrepresentable"
    );
    assert_eq!(candidate, diagonal);
    // Rounded lengths 4, 2, 4 are valid; equal parameter spacing is unnecessary.
    let mut horizontal = single([1e16, 0.0, 0.0], [1e16 + 10.0, 0.0, 0.0], 0.04);
    assert_eq!(
        horizontal
            .refine_segments(4.5, 3, || false)
            .unwrap()
            .segments_after,
        3
    );
    let mut near = sheet(0.5, 4.0);
    near.rows.truncate(1);
    near.rows[0].nodes = vec![[1e16, 0.0, 2.0], [1e16, 0.0, 3.0]];
    near.line_nodes = vec![[1e16 + 4.0, 4.0, 2.0], [1e16 + 4.0, 4.0, 3.0]];
    let snapshot = near.clone();
    assert_eq!(
        near.refine_streamwise(2.5, 2, || false).unwrap_err().code,
        "wake-refinement-unrepresentable"
    );
    assert_eq!(near, snapshot);
}

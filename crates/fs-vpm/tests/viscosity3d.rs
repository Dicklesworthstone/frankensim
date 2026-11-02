//! G0/G2/G3 regressions for the finite-core 3-D wake.
//! Independent quadrature validates the reduced regularized-filament model.
use fs_vpm::farfield3d::{
    FarField, WakeCoreEvolutionMode, cell_from_segs, cell_from_segs_with_cores,
};
use fs_vpm::filament3d::{
    FilamentWake, WakeRateCertificate, segment_velocity, segment_velocity_with_core,
};
use fs_vpm::images3d::{ImagePlane, ImageSystem, grounded_hybrid_velocity, mirror_cell};
use std::f64::consts::PI;

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn difference(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
}
fn wake() -> FilamentWake {
    let cert = WakeRateCertificate {
        shed_hz: 100.0,
        n_stations: 1,
        max_rows: 16,
    };
    let mut wake = FilamentWake::new(cert, vec![[0.0, -0.5, 1.0], [0.0, 0.5, 1.0]]).unwrap();
    for _ in 0..8 {
        wake.shed(&[1.0], [-0.1, 0.0, 0.0]).unwrap();
    }
    wake
}

#[test]
fn g2_finite_segment_matches_independent_regularized_line_quadrature() {
    let a = [-0.6, -0.7, 0.1];
    let b = [0.4, 0.9, -0.2];
    let p = [0.2, 0.1, 0.4];
    let core2 = 0.0625;
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let n = 2048;
    let mut oracle = [0.0; 3];
    // Composite Simpson integration, independent of the endpoint formula.
    for i in 0..=n {
        let t = f64::from(i) / f64::from(n);
        let r = [
            p[0] - (a[0] + t * ab[0]),
            p[1] - (a[1] + t * ab[1]),
            p[2] - (a[2] + t * ab[2]),
        ];
        let cross = [
            ab[1] * r[2] - ab[2] * r[1],
            ab[2] * r[0] - ab[0] * r[2],
            ab[0] * r[1] - ab[1] * r[0],
        ];
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + core2;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for j in 0..3 {
            oracle[j] += weight * cross[j] / (4.0 * PI * r2 * r2.sqrt());
        }
    }
    for v in &mut oracle {
        *v /= 3.0 * f64::from(n);
    }
    let actual = segment_velocity_with_core(p, a, b, core2);
    assert!(difference(actual, oracle) < 1e-10 * norm(oracle));
}

#[test]
fn g2_long_filament_has_finite_algebraic_core_and_softens_swirl() {
    let a = [0.0, 0.0, -1e6];
    let b = [0.0, 0.0, 1e6];
    let core2 = 0.04;
    for radius in [1e-12, 0.03, 0.2, 1.0] {
        let actual = segment_velocity_with_core([radius, 0.0, 0.0], a, b, core2);
        let expected = radius / (2.0 * PI * (radius * radius + core2));
        assert!(actual[0] == 0.0 && actual[2] == 0.0);
        assert!((actual[1] / expected - 1.0).abs() < 1e-10);
        assert!(actual[1] > 0.0, "no artificial zero-velocity core");
        assert!(actual[1] < 1.0 / (2.0 * PI * radius));
    }
    assert_eq!(segment_velocity_with_core([0.0; 3], a, b, core2), [0.0; 3]);
    assert_eq!(segment_velocity_with_core([1.0; 3], a, a, core2), [0.0; 3]);
}

#[test]
fn g3_core_kernel_scales_with_length_and_recovers_zero_core() {
    let a = [0.0, 0.0, -1.0];
    let b = [0.0, 0.0, 1.0];
    let p = [0.3, 0.0, 0.0];
    let inviscid = segment_velocity(p, a, b);
    let expected = 1.0 / (2.0 * PI * 0.3 * (1.0_f64 + 0.3 * 0.3).sqrt());
    assert!((inviscid[1] / expected - 1.0).abs() < 1e-14);
    assert_eq!(inviscid, segment_velocity_with_core(p, a, b, 0.0));
    for core2 in [0.0, 0.04] {
        let reference = segment_velocity_with_core(p, a, b, core2);
        for scale in [1e-9, 1e-3, 1e3, 1e9] {
            let actual = segment_velocity_with_core(
                p.map(|v| v * scale),
                a.map(|v| v * scale),
                b.map(|v| v * scale),
                core2 * scale * scale,
            );
            assert!(
                difference(actual.map(|v| v * scale), reference) < 1e-12 * norm(reference),
                "length scale {scale}, core² {core2}"
            );
        }
    }
}

#[test]
fn g0_near_wake_and_retained_exact_segments_use_the_same_core() {
    let cert = WakeRateCertificate {
        shed_hz: 100.0,
        n_stations: 1,
        max_rows: 2,
    };
    let a = [0.0, -0.5, 1.0];
    let b = [0.0, 0.5, 1.0];
    let mut near = FilamentWake::new(cert, vec![a, b]).unwrap();
    near.shed(&[2.0], [0.0; 3]).unwrap();
    let p = [0.1, 0.0, 1.1];
    let inviscid = near.induced_velocity(p);
    near.rows[0].core2_m2 = 0.09;
    let cell = cell_from_segs(vec![(a, b, 2.0)], 0.09);
    assert_eq!(near.induced_velocity(p), cell.exact(p));
    assert!(norm(near.induced_velocity(p)) < 0.5 * norm(inviscid));
}

#[test]
fn g0_spreading_changes_flow_and_zero_viscosity_recovers_frozen_flow() {
    let mut a = wake();
    let mut b = a.clone();
    let mut c = a.clone();
    let frozen = FarField::aggregate(&mut a, 4, 1, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    let zero = FarField::aggregate(
        &mut b,
        4,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.0 },
        0.01,
    )
    .unwrap();
    let viscous = FarField::aggregate(
        &mut c,
        4,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap();
    let p = [-0.7, 0.0, 1.2];
    assert_eq!(frozen.eval(p), zero.eval(p));
    assert!(difference(frozen.eval(p), viscous.eval(p)) > 1e-3);
    for (frozen, (zero, viscous)) in frozen
        .cells
        .iter()
        .zip(zero.cells.iter().zip(&viscous.cells))
    {
        assert_eq!(frozen.exact(p), zero.exact(p));
        assert!(difference(frozen.exact(p), viscous.exact(p)) > 1e-3);
        assert_eq!(
            frozen.w_sum, viscous.w_sum,
            "diffusion preserves circulation moments"
        );
    }
}

#[test]
fn g3_viscous_multipoles_converge_to_the_same_regularized_exact_field() {
    let segs = vec![
        ([-0.2, -0.5, 0.1], [0.1, 0.6, -0.1], 1.2),
        ([0.3, -0.2, -0.2], [0.5, 0.3, 0.2], -0.4),
    ];
    for core2 in [0.0, 0.04, 4.0] {
        let cell = cell_from_segs(segs.clone(), core2);
        let mut errors = Vec::new();
        for distance in [8.0, 16.0, 32.0] {
            let p = [
                cell.centroid[0] + 0.3 * distance,
                cell.centroid[1] + 0.6 * distance,
                cell.centroid[2] + 0.74 * distance,
            ];
            let exact = cell.exact(p);
            errors.push(difference(cell.eval(p), exact) / norm(exact));
        }
        assert!(errors[2] < 1e-3, "core² {core2}: {errors:?}");
        assert!(errors[0] > 3.0 * errors[1] && errors[1] > 3.0 * errors[2]);
    }
}

#[test]
fn g0_viscous_exact_field_stays_below_pruning_bound() {
    let cell = cell_from_segs(vec![([0.0, -0.5, 1.0], [0.0, 0.5, 1.0], 2.0)], 0.25);
    for p in [[0.2, 0.0, 1.1], [3.0, 0.0, 1.0], [10.0, 4.0, 3.0]] {
        assert!(norm(cell.exact(p)) <= cell.contribution_bound(p));
    }
    let probe = [10.0, 4.0, 3.0];
    let tolerance = 2.0 * cell.contribution_bound(probe);
    let mut far = FarField {
        mode: WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.05 },
        shed_dt_s: 0.01,
        cells: vec![cell],
    };
    let receipt = far.prune_audited(&[probe], tolerance).unwrap();
    assert_eq!(receipt.pruned.len(), 1);
    assert!(receipt.pruned[0].exact_worst <= receipt.pruned[0].bound_worst);
}

#[test]
fn g3_images_preserve_core_and_cancel_on_plane_normal_flow() {
    let mut near = wake();
    for row in &mut near.rows {
        row.core2_m2 = 0.16;
    }
    let far = FarField::aggregate(&mut near, 4, 2, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    let plane = ImagePlane { z0_m: 0.0 };
    let images = ImageSystem::build(&near, &far, plane).unwrap();
    let cell = &far.cells[0];
    let mirror = mirror_cell(cell, plane);
    let p = [-0.6, 0.2, 1.1];
    let expected = cell.exact(p);
    let actual = mirror.exact(plane.mirror(p));
    assert!(difference(actual, [expected[0], expected[1], -expected[2]]) < 1e-13);
    for p in [[-0.5, 0.2, 0.0], [0.2, -0.3, 0.0], [1.0, 0.0, 0.0]] {
        let v = grounded_hybrid_velocity(&near, &far, &images, p);
        assert!(v[2].abs() < 1e-12 * (1.0 + norm(v)));
    }
}

#[test]
fn g0_cell_uses_exact_near_field_and_softened_expansion_inside_a_wide_core() {
    let cell = cell_from_segs(vec![([0.0, -0.005, 0.0], [0.0, 0.005, 0.0], 1.0)], 1.0);
    let near = [0.001, 0.0, 0.0];
    assert_eq!(cell.eval(near), cell.exact(near));
    // Outside the source geometry, inside the viscous core.
    let overlap = [0.05, 0.0, 0.0];
    let exact = cell.exact(overlap);
    assert!(difference(cell.eval(overlap), exact) < 2e-5 * norm(exact));
    let inviscid = cell_from_segs(cell.segs.clone(), 0.0);
    assert!(norm(cell.eval(overlap)) < 1e-3 * norm(inviscid.eval(overlap)));
}

#[test]
fn g0_unrepresentable_spreading_refuses_without_removing_wake_rows() {
    let mut near = wake();
    let original = near.clone();
    let error = FarField::aggregate(
        &mut near,
        4,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 1.0 },
        f64::MAX,
    )
    .unwrap_err();
    assert_eq!(error.code, "farfield-core-invalid");
    assert_eq!(near, original);
}

#[test]
fn g2_thin_filament_retains_the_finite_field_beyond_both_endpoints() {
    let a = [0.0, 0.0, 0.0];
    let b = [0.0, 0.0, 1.0];
    let radius = 1e-12;
    // In this thin-filament limit the line integral is
    // radius * integral_1^2 z^-3 dz / (4 pi) = 3 radius / (32 pi).
    let expected = 3.0 * radius / (32.0 * PI);
    for z in [-1.0, 2.0] {
        for core2 in [0.0, 1e-24] {
            let v = segment_velocity_with_core([radius, 0.0, z], a, b, core2);
            assert!((v[1] / expected - 1.0).abs() < 1e-14);
            let reversed = segment_velocity_with_core([radius, 0.0, z], b, a, core2);
            assert!(difference(v, reversed.map(|x| -x)) < 1e-14 * norm(v));
        }
    }
}

#[test]
fn g0_mixed_source_cores_survive_exact_evaluation_and_images() {
    let segs = vec![
        ([-0.2, -0.4, 1.0], [-0.2, 0.4, 1.0], 1.2),
        ([0.3, -0.1, 1.2], [0.7, 0.3, 1.2], -0.5),
    ];
    let cores = vec![0.01, 0.9];
    let cell = cell_from_segs_with_cores(segs.clone(), cores.clone());
    let p = [0.1, 0.1, 1.05];
    let mut expected = [0.0; 3];
    for ((a, b, g), core2) in segs.iter().zip(&cores) {
        let v = segment_velocity_with_core(p, *a, *b, *core2);
        for k in 0..3 {
            expected[k] += g * v[k];
        }
    }
    assert_eq!(cell.exact(p), expected);
    let plane = ImagePlane { z0_m: 0.0 };
    let image = mirror_cell(&cell, plane);
    assert_eq!(image.source_core2_m2, cores);
    let mirrored = image.exact(plane.mirror(p));
    assert!(difference(mirrored, [expected[0], expected[1], -expected[2]]) < 1e-12);
}

#[test]
fn g0_large_core_variation_uses_original_sources_outside_geometry() {
    let cell = cell_from_segs_with_cores(
        vec![
            ([0.0, -0.5, 0.0], [0.0, 0.5, 0.0], 1.0),
            ([0.1, -0.5, 0.0], [0.1, 0.5, 0.0], 1.0),
        ],
        vec![0.0, 1e6],
    );
    let p = [3.0, 0.0, 1.0];
    let d2 = (p[0] - cell.centroid[0]).powi(2)
        + (p[1] - cell.centroid[1]).powi(2)
        + (p[2] - cell.centroid[2]).powi(2);
    assert!(d2 > cell.radius_m * cell.radius_m);
    assert!(cell.max_abs_core_delta_m2 > (d2 + cell.core2_m2) / 16.0);
    assert_eq!(cell.eval(p), cell.exact(p));
}

#[test]
fn g3_mixed_core_far_expansion_converges_and_scales_without_blanket_fallback() {
    let segs = vec![
        ([-0.2, -0.5, 0.1], [0.1, 0.6, -0.1], 1.2),
        ([0.3, -0.2, -0.2], [0.5, 0.3, 0.2], -0.4),
    ];
    let cores = vec![0.01, 4.0];
    let cell = cell_from_segs_with_cores(segs.clone(), cores.clone());
    let mut errors = Vec::new();
    for distance in [16.0, 32.0, 64.0] {
        let p = [
            cell.centroid[0] + 0.3 * distance,
            cell.centroid[1] + 0.6 * distance,
            cell.centroid[2] + 0.74 * distance,
        ];
        let d2 = (p[0] - cell.centroid[0]).powi(2)
            + (p[1] - cell.centroid[1]).powi(2)
            + (p[2] - cell.centroid[2]).powi(2);
        assert!(cell.max_abs_core_delta_m2 <= (d2 + cell.core2_m2) / 16.0);
        let exact = cell.exact(p);
        let approx = cell.eval(p);
        let error = difference(approx, exact) / norm(exact);
        assert!(
            error > 0.0,
            "a separated mixed-core cell still uses its expansion"
        );
        errors.push(error);
        for scale in [1e-3, 1e3] {
            let scaled = cell_from_segs_with_cores(
                segs.iter()
                    .map(|(a, b, g)| (a.map(|x| x * scale), b.map(|x| x * scale), *g))
                    .collect(),
                cores.iter().map(|c| c * scale * scale).collect(),
            );
            let value = scaled.eval(p.map(|x| x * scale)).map(|x| x * scale);
            assert!(difference(value, approx) < 1e-11 * norm(approx));
        }
    }
    assert!(errors[2] < 2e-4, "{errors:?}");
    assert!(
        errors[0] > 4.0 * errors[1] && errors[1] > 4.0 * errors[2],
        "{errors:?}"
    );
}

#[test]
fn g0_digest_identifies_which_source_carries_each_core() {
    let segs = vec![
        ([-1.0, -0.5, 0.0], [-1.0, 0.5, 0.0], 1.0),
        ([1.0, -0.5, 0.0], [1.0, 0.5, 0.0], 1.0),
    ];
    let a = cell_from_segs_with_cores(segs.clone(), vec![0.01, 1.0]);
    let b = cell_from_segs_with_cores(segs, vec![1.0, 0.01]);
    assert_eq!(a.core2_m2, b.core2_m2);
    assert_eq!(a.w_sum, b.w_sum);
    assert_eq!(a.dip, b.dip);
    assert_eq!(a.quad, b.quad);
    assert_eq!(a.core_moment, b.core_moment);
    assert!(difference(a.exact([-1.0, 0.0, 0.1]), b.exact([-1.0, 0.0, 0.1])) > 0.1);
    let left = FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        shed_dt_s: 0.01,
        cells: vec![a],
    };
    let right = FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        shed_dt_s: 0.01,
        cells: vec![b],
    };
    assert_ne!(left.digest(), right.digest());
}

#[test]
fn g0_rows_start_at_zero_age_and_conversion_does_not_reset_remaining_ages() {
    let mut near = wake();
    assert_eq!(near.rows.last().unwrap().age_ticks, 0.0);
    assert_eq!(near.rows[0].age_ticks, 7.0);
    near.shed(&[1.0], [-0.1, 0.0, 0.0]).unwrap();
    let ages: Vec<_> = near.rows.iter().map(|r| r.age_ticks).collect();
    assert_eq!(ages, vec![8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]);
    let before = near.clone();
    assert!(near.shed(&[f64::NAN], [0.0; 3]).is_err());
    assert_eq!(near, before, "rejected shedding must not age existing rows");
    let _far = FarField::aggregate(&mut near, 4, 1, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    assert_eq!(
        near.rows.iter().map(|r| r.age_ticks).collect::<Vec<_>>(),
        ages[4..]
    );
}

#[test]
fn g0_coarsening_retains_weighted_age_and_aggregation_uses_it() {
    let mut near = wake();
    near.rows[0].gamma[0] = 1.0;
    near.rows[1].gamma[0] = 3.0;
    fs_vpm::coarsen3d::coarsen_oldest(&mut near, 1).unwrap();
    assert!((near.rows[0].age_ticks - 6.25).abs() < 1e-14);
    assert_eq!(near.rows[1].age_ticks, 5.0);
    let expected = near.rows[0].core2_m2 + 0.02 * near.rows[0].age_ticks;
    let far = FarField::aggregate(
        &mut near,
        2,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap();
    assert!(
        far.cells[0]
            .source_core2_m2
            .iter()
            .all(|c| (*c - expected).abs() < 1e-14)
    );
    assert_eq!(
        near.rows[0].age_ticks, 4.0,
        "row removal is not elapsed time"
    );
}

#[test]
fn g3_coarsening_is_invariant_to_circulation_unit_scaling() {
    let mut normal = wake();
    for (i, row) in normal.rows.iter_mut().enumerate() {
        row.gamma[0] = (i + 1) as f64;
        row.core2_m2 = 0.002 * (i + 1) as f64;
    }
    let mut small = normal.clone();
    let scale = 1e-18;
    for row in &mut small.rows {
        for gamma in &mut row.gamma {
            *gamma *= scale;
        }
    }
    fs_vpm::coarsen3d::coarsen_oldest(&mut normal, 2).unwrap();
    fs_vpm::coarsen3d::coarsen_oldest(&mut small, 2).unwrap();
    for (a, b) in normal.rows.iter().zip(&small.rows) {
        assert!((a.age_ticks - b.age_ticks).abs() < 1e-13);
        assert!((a.core2_m2 - b.core2_m2).abs() < 1e-13);
        for (x, y) in a.nodes.iter().zip(&b.nodes) {
            assert!(difference(*x, *y) < 1e-13);
        }
        for (ga, gb) in a.gamma.iter().zip(&b.gamma) {
            assert!((ga - gb / scale).abs() < 1e-12 * ga.abs().max(1.0));
        }
    }
    let normal_impulse = fs_vpm::coarsen3d::station_impulse(&normal, 0);
    let scaled_impulse = fs_vpm::coarsen3d::station_impulse(&small, 0) / scale;
    assert!((normal_impulse - scaled_impulse).abs() < 1e-12 * normal_impulse.abs());
}

#[test]
fn g3_aggregation_keeps_one_viscous_field_across_near_and_far() {
    let full = wake();
    let growth_per_tick = 0.02;
    let p = [-0.25, 0.1, 1.12];
    let expected = full.induced_velocity_with_core_growth(p, growth_per_tick);
    let mut near = full.clone();
    let far = FarField::aggregate(
        &mut near,
        4,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap();
    let aged_near = near.induced_velocity_with_core_growth(p, growth_per_tick);
    let mut exact = aged_near;
    for cell in &far.cells {
        let v = cell.exact(p);
        for k in 0..3 {
            exact[k] += v[k];
        }
    }
    assert!(difference(exact, expected) < 1e-12 * (1.0 + norm(expected)));
    // These probes lie inside all removed cells' bounding spheres, so the
    // public hybrid selects exact sources too: this tests coupling, not truncation.
    let actual = fs_vpm::farfield3d::hybrid_velocity(&near, &far, p);
    assert!(difference(actual, expected) < 1e-12 * (1.0 + norm(expected)));
    assert!(
        difference(aged_near, near.induced_velocity(p)) > 1e-4,
        "the remaining near wake must also feel viscosity"
    );
}

#[test]
fn g3_images_cancel_viscous_near_and_continually_diffused_far_normal_flow() {
    let mut near = wake();
    let mut far = FarField::aggregate(
        &mut near,
        4,
        2,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.4 },
        0.01,
    )
    .unwrap();
    far.advance_core(0.05).unwrap();
    let plane = ImagePlane { z0_m: 0.0 };
    let images = ImageSystem::build(&near, &far, plane).unwrap();
    for (source, image) in far.cells.iter().zip(images.far_image_cells()) {
        assert_eq!(source.diffusion_offset_m2, image.diffusion_offset_m2);
        assert_eq!(source.source_core2_m2, image.source_core2_m2);
    }
    for p in [[-0.5, 0.2, 0.0], [0.2, -0.3, 0.0], [1.0, 0.0, 0.0]] {
        let v = grounded_hybrid_velocity(&near, &far, &images, p);
        assert!(
            v[2].abs() < 1e-12 * (1.0 + norm(v)),
            "normal residual {v:?}"
        );
    }
}

#[test]
fn g0_continued_diffusion_matches_each_original_source_and_repeated_steps() {
    let segs = vec![
        ([-0.2, -0.4, 1.0], [-0.2, 0.4, 1.0], 1.2),
        ([0.3, -0.1, 1.2], [0.7, 0.3, 1.2], -0.5),
    ];
    let cores = vec![0.01, 0.2];
    let cell = cell_from_segs_with_cores(segs.clone(), cores.clone());
    let original = FarField {
        mode: WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.25 },
        shed_dt_s: 0.01,
        cells: vec![cell],
    };
    let mut once = original.clone();
    let mut repeated = original.clone();
    once.advance_core(0.5).unwrap();
    for _ in 0..4 {
        repeated.advance_core(0.125).unwrap();
    }
    assert_eq!(once.cells[0].diffusion_offset_m2, 0.5);
    assert_eq!(once.cells[0].source_core2_m2, cores);
    assert_eq!(
        once.cells[0].segs, segs,
        "diffusion does not pretend to convect sources"
    );
    assert_eq!(once.cells[0].core_moment, original.cells[0].core_moment);
    let rebuilt = cell_from_segs_with_cores(segs, vec![0.51, 0.7]);
    for p in [[0.1, 0.1, 1.05], [10.0, 3.0, 4.0]] {
        let expected = rebuilt.exact(p);
        assert!(difference(once.cells[0].exact(p), expected) < 1e-12 * (1.0 + norm(expected)));
        assert!(difference(repeated.cells[0].exact(p), expected) < 1e-12 * (1.0 + norm(expected)));
        let expected_far = rebuilt.eval(p);
        assert!(difference(once.eval(p), expected_far) < 1e-12 * (1.0 + norm(expected_far)));
    }
    assert_ne!(once.digest(), original.digest());
}

#[test]
fn g0_diffusion_identity_and_overflow_refusal_are_transactional() {
    let segment = vec![([0.0, -0.5, 1.0], [0.0, 0.5, 1.0], 1.0)];
    let original = FarField {
        mode: WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.25 },
        shed_dt_s: 0.01,
        cells: vec![
            cell_from_segs(segment.clone(), 0.1),
            cell_from_segs(segment, f64::MAX),
        ],
    };
    let mut far = original.clone();
    let initial_digest = far.digest();
    far.advance_core(0.0).unwrap();
    assert_eq!(far, original);
    assert_eq!(far.digest(), initial_digest);
    // The first cell could advance, but the second would overflow. Nothing commits.
    assert_eq!(
        far.advance_core(f64::MAX).unwrap_err().code,
        "core-advance-overflow"
    );
    assert_eq!(far, original);
    assert_eq!(far.digest(), initial_digest);
    for bad in [-0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            far.advance_core(bad).unwrap_err().code,
            "core-advance-invalid"
        );
        assert_eq!(far, original);
    }
    for mode in [
        WakeCoreEvolutionMode::Frozen,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.0 },
    ] {
        let mut identity = original.clone();
        identity.mode = mode;
        let before = identity.clone();
        let digest = identity.digest();
        identity.advance_core(f64::MAX).unwrap();
        assert_eq!(identity, before);
        assert_eq!(identity.digest(), digest);
    }
}

#[test]
fn g0_aggregation_validates_ages_of_rows_that_remain_near() {
    let mut near = wake();
    near.rows.last_mut().unwrap().age_ticks = f64::NAN;
    let digest = near.digest();
    let count = near.rows.len();
    let error = FarField::aggregate(
        &mut near,
        4,
        1,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap_err();
    assert_eq!(error.code, "farfield-core-invalid");
    assert_eq!(near.rows.len(), count);
    assert_eq!(near.digest(), digest);
}

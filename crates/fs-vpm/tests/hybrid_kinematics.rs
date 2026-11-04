//! G0/G1/G3: actual near/far/image wake curl survives physical partitioning,
//! continued core spreading, reflection and independent velocity differentiation.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs, hybrid_kinematics};
use fs_vpm::filament3d::{FilamentWake, WakeKinematics, WakeRateCertificate};
use fs_vpm::images3d::{ImagePlane, ImageSystem, grounded_hybrid_kinematics};

const DT: f64 = 0.01;

fn wake() -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 3,
            max_rows: 12,
        },
        vec![
            [0.0, -0.8, 1.4],
            [0.0, -0.2, 1.5],
            [0.0, 0.4, 1.6],
            [0.0, 0.9, 1.5],
        ],
    )
    .unwrap();
    for k in 0..8 {
        let t = k as f64;
        wake.shed(
            &[1.0 + 0.1 * t, 0.4 - 0.05 * t, -0.3 + 0.02 * t],
            [-0.15, 0.02, 0.01],
        )
        .unwrap();
    }
    for (i, row) in wake.rows.iter_mut().enumerate() {
        row.core2_m2 = 0.008 + 0.003 * i as f64;
    }
    wake
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "actual={a:e}, expected={b:e}, relative/absolute tolerance={tolerance:e}"
    );
}

fn same_field(a: WakeKinematics, b: WakeKinematics, tolerance: f64) {
    for i in 0..3 {
        near(a.velocity[i], b.velocity[i], tolerance);
        near(a.vorticity()[i], b.vorticity()[i], tolerance);
        for j in 0..3 {
            near(a.gradient[i][j], b.gradient[i][j], tolerance);
        }
    }
    near(a.divergence(), b.divergence(), tolerance);
}

fn exact_velocity(near: &FilamentWake, far: &FarField, p: [f64; 3]) -> [f64; 3] {
    let mut u = near.induced_velocity_with_core_growth(p, far.near_core_growth_m2().unwrap());
    for cell in &far.cells {
        let v = cell.exact(p);
        for i in 0..3 {
            u[i] += v[i];
        }
    }
    u
}

#[test]
fn mixed_core_partition_and_continued_diffusion_preserve_every_derivative() {
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.07 };
    for transferred in [2, 4, 6] {
        let mut full = wake();
        let mut near_wake = full.clone();
        let mut far = FarField::aggregate(&mut near_wake, transferred, 2, mode, DT).unwrap();
        for age in [0.0, 0.1] {
            if age > 0.0 {
                far.advance_core(age).unwrap();
                for row in full.rows.iter_mut().chain(&mut near_wake.rows) {
                    row.age_ticks += age / DT;
                }
            }
            for p in [[-0.55, -0.1, 1.45], [1.0, 0.3, 1.9], [-3.0, 2.2, 5.0]] {
                let reference = full
                    .induced_kinematics_with_core_growth(p, far.near_core_growth_m2().unwrap())
                    .unwrap();
                let actual = hybrid_kinematics(&near_wake, &far, p).unwrap();
                same_field(actual, reference, 5e-13);
                let velocity = exact_velocity(&near_wake, &far, p);
                for i in 0..3 {
                    near(actual.velocity[i], velocity[i], 2e-14);
                }
            }
        }
    }
}

#[test]
fn exact_hybrid_jacobian_matches_independent_whole_velocity_differences() {
    let mut near_wake = wake();
    let far = FarField::aggregate(
        &mut near_wake,
        4,
        2,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.03 },
        DT,
    )
    .unwrap();
    for p in [[-0.55, -0.1, 1.45], [1.0, 0.3, 1.9], [-3.0, 2.2, 5.0]] {
        let sample = hybrid_kinematics(&near_wake, &far, p).unwrap();
        for axis in 0..3 {
            let mut plus = p;
            let mut minus = p;
            let h = 1e-5;
            plus[axis] += h;
            minus[axis] -= h;
            let up = exact_velocity(&near_wake, &far, plus);
            let um = exact_velocity(&near_wake, &far, minus);
            for component in 0..3 {
                near(
                    sample.gradient[component][axis],
                    (up[component] - um[component]) / (2.0 * h),
                    3e-8,
                );
            }
        }
        near(sample.divergence(), 0.0, 3e-13);
    }
}

#[test]
fn images_reflect_velocity_gradient_and_axial_vorticity_without_changing_sources() {
    let mut near_wake = wake();
    let far = FarField::aggregate(
        &mut near_wake,
        4,
        2,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.03 },
        DT,
    )
    .unwrap();
    let source_digests = (near_wake.digest(), far.digest());
    let plane = ImagePlane { z0_m: 0.25 };
    let images = ImageSystem::build(&near_wake, &far, plane).unwrap();
    let reflection = [1.0, 1.0, -1.0];
    for p in [[-0.55, -0.1, 1.45], [1.0, 0.3, 1.9], [-3.0, 2.2, 0.25]] {
        let source = hybrid_kinematics(&near_wake, &far, plane.mirror(p)).unwrap();
        let image = images.exact_kinematics(p).unwrap();
        for i in 0..3 {
            near(image.velocity[i], reflection[i] * source.velocity[i], 3e-14);
            near(
                image.vorticity()[i],
                -reflection[i] * source.vorticity()[i],
                3e-13,
            );
            for j in 0..3 {
                near(
                    image.gradient[i][j],
                    reflection[i] * reflection[j] * source.gradient[i][j],
                    3e-13,
                );
            }
        }
    }
    for p in [[-0.5, 0.0, 0.25], [0.3, 0.7, 0.25]] {
        let field = grounded_hybrid_kinematics(&near_wake, &far, &images, p).unwrap();
        near(field.velocity[2], 0.0, 3e-14);
        near(field.gradient[2][0], 0.0, 3e-14);
        near(field.gradient[2][1], 0.0, 3e-14);
        near(field.vorticity()[0], 0.0, 3e-14);
        near(field.vorticity()[1], 0.0, 3e-14);
    }
    assert_eq!((near_wake.digest(), far.digest()), source_digests);
}

#[test]
fn exact_far_sampling_uses_retained_cores_not_expansion_caches() {
    let mut cell = cell_from_segs(vec![([0.0, -1.0, 1.0], [0.0, 1.0, 1.0], 2.0)], 0.04);
    cell.diffusion_offset_m2 = 0.01;
    let p = [0.0, 0.0, 1.0];
    let sample = cell.exact_kinematics(p).unwrap();
    assert_eq!(sample.velocity, [0.0; 3]);
    let expected = 2.0 / (core::f64::consts::PI * 0.05 * 1.05_f64.sqrt());
    near(sample.vorticity()[1], expected, 2e-14);
    cell.core2_m2 = f64::NAN;
    cell.centroid = [f64::NAN; 3];
    cell.w_sum = [f64::NAN; 3];
    assert_eq!(cell.exact_kinematics(p).unwrap(), sample);
    // Inactive zero-core geometry has no singular physical source.
    let mut zero = cell_from_segs(vec![([0.0, -1.0, 1.0], [0.0, 1.0, 1.0], 0.0)], 0.0);
    assert_eq!(zero.exact_kinematics(p), Some(WakeKinematics::default()));
    zero.segs[0].2 = 1.0;
    assert!(zero.exact_kinematics(p).is_none());
    zero.source_core2_m2[0] = 0.05;
    assert!(zero.exact_kinematics(p).is_some());
}

#[test]
fn malformed_sources_refuse_before_image_construction_or_sampling() {
    let mut near_wake = wake();
    let far = FarField::aggregate(&mut near_wake, 4, 2, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let p = [-0.2, 0.1, 1.4];
    let plane = ImagePlane { z0_m: 0.0 };
    for defect in 0..7 {
        let mut invalid_near = near_wake.clone();
        let mut invalid_far = far.clone();
        match defect {
            0 => {
                invalid_near.rows[0].nodes.pop();
            }
            1 => {
                invalid_near.downstream_gamma.as_mut().unwrap().pop();
            }
            2 => invalid_near.rows[0].age_ticks = f64::NAN,
            3 => {
                invalid_far.cells[0].source_core2_m2.pop();
            }
            4 => invalid_far.cells[0].source_core2_m2[0] = -1.0,
            5 => invalid_far.cells[0].segs[0].2 = f64::NAN,
            6 => invalid_far.shed_dt_s *= 2.0,
            _ => unreachable!(),
        }
        assert!(
            invalid_far
                .validate_kinematics_sources(&invalid_near)
                .is_err()
        );
        assert!(hybrid_kinematics(&invalid_near, &invalid_far, p).is_none());
        assert!(ImageSystem::build(&invalid_near, &invalid_far, plane).is_err());
    }
    assert!(far.exact_kinematics([f64::NAN, 0.0, 0.0]).is_none());
    let mut invalid_far = far.clone();
    invalid_far.cells[0].diffusion_offset_m2 = f64::INFINITY;
    assert!(invalid_far.exact_kinematics(p).is_none());
}

#[test]
fn images_reject_far_sources_that_cross_the_certified_plane() {
    let mut near_wake = wake();
    let mut far =
        FarField::aggregate(&mut near_wake, 4, 2, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    // Near clearance alone cannot certify already-aggregated physical sources.
    far.cells[0].segs[0].0[2] = -0.01;
    assert_eq!(
        ImageSystem::build(&near_wake, &far, ImagePlane { z0_m: 0.0 })
            .unwrap_err()
            .code,
        "image-plane-invalid"
    );
}

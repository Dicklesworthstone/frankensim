//! G1/G3: physical lens inverse maps and explicit local failure boundaries.

use super::*;
use fs_render::camera::{SequentialProjectionError, SequentialProjectionOptions};

fn camera(lens: SequentialLens, eye: Point3, direction: Vec3) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        eye,
        eye.offset(direction),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.012).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

fn target_on_ray(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    tangent: [f64; 2],
    pupil: LensSample,
    wavelength: f64,
) -> Point3 {
    let [sample] = camera
        .generate_sequential_spectral_rays_from_tangent_offsets(
            cx,
            tangent[0],
            tangent[1],
            pupil,
            &[wavelength],
        )
        .unwrap();
    let sample = sample.expect("authored unvignetted native path");
    sample.ray.origin.offset(sample.ray.dir.scale(0.7))
}

fn options(seed: [f64; 2]) -> SequentialProjectionOptions {
    SequentialProjectionOptions {
        initial_tangent_xy: seed,
        position_tolerance_m: 1e-10,
        ..SequentialProjectionOptions::default()
    }
}

#[test]
fn inverse_projection_retains_native_spectral_pupil_ghost_and_world_pose() {
    with_cx(|_, cx| {
        let material = glass(CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(), 3.0);
        for ghosts in [false, true] {
            let lens = if ghosts {
                biconvex(material.clone())
                    .with_two_reflection_ghosts(0.7)
                    .unwrap()
            } else {
                biconvex(material.clone())
            };
            let camera = camera(lens, Point3::new(2.0, -1.0, 4.0), Vec3::new(0.2, -0.1, 1.0));
            // The same original u selects a fixed direct or ghost interval;
            // the native proposal alone remaps the conditional pupil.
            let pupils = if ghosts {
                [
                    LensSample::try_new(0.35, 0.5).unwrap(),
                    LensSample::try_new(0.85, 0.5).unwrap(),
                ]
            } else {
                [LensSample::CENTER, LensSample::try_new(0.55, 0.45).unwrap()]
            };
            for pupil in pupils {
                for wavelength in [450.0, 650.0] {
                    let expected = [0.004, -0.003];
                    let point = target_on_ray(&camera, cx, expected, pupil, wavelength);
                    let result = camera
                        .project_sequential_point(
                            cx,
                            point,
                            1.5,
                            wavelength,
                            pupil,
                            options([expected[0] + 0.002, expected[1] - 0.001]),
                        )
                        .unwrap();
                    for axis in 0..2 {
                        close(result.tangent_xy()[axis], expected[axis], 2e-8);
                    }
                    close(
                        result.ndc_xy()[0],
                        expected[0] / 1.5 / camera.raster_vertical_half_tan(),
                        3e-7,
                    );
                    close(
                        result.ndc_xy()[1],
                        expected[1] / camera.raster_vertical_half_tan(),
                        3e-7,
                    );
                    assert!(result.residual_m() <= 1e-10);
                    assert!(result.trace_evaluations() <= 128);
                    close(
                        result.depth_m(),
                        point.delta_from(camera.eye()).dot(camera.forward()),
                        1e-14,
                    );
                    let tangent = result.tangent_xy();
                    let [replayed] = camera
                        .generate_sequential_spectral_rays_from_tangent_offsets(
                            cx,
                            tangent[0],
                            tangent[1],
                            pupil,
                            &[wavelength],
                        )
                        .unwrap();
                    let ray = replayed.unwrap().ray;
                    // Independent closest-point check on the full world ray,
                    // rather than reproducing the solver's axial-plane residual.
                    let delta = point.delta_from(ray.origin);
                    let distance = delta.dot(ray.dir);
                    assert!(distance > 0.0);
                    let miss = Vec3::new(
                        delta.x - ray.dir.x * distance,
                        delta.y - ray.dir.y * distance,
                        delta.z - ray.dir.z * distance,
                    );
                    assert!(miss.dot(miss).sqrt() < 2e-9);
                }
            }
        }
    });
}

#[test]
fn inverse_projection_keeps_chromatic_geometry_and_steady_coating_eligibility() {
    with_cx(|_, cx| {
        let material = glass(CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(), 2.0);
        let bare_lens = biconvex(material);
        let film = coating(&[(1.3, 0.03, 110e-9)]);
        let coated_lens = lens(&[
            bare_lens.surfaces()[0].with_coating(film),
            bare_lens.surfaces()[1].with_coating(film),
        ]);
        let bare = camera(
            bare_lens,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let coated = camera(
            coated_lens,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let pupil = LensSample::try_new(0.55, 0.45).unwrap();
        let point = target_on_ray(&bare, cx, [0.008, -0.004], pupil, 450.0);
        let mut projections = [[0.0; 2]; 2];
        for (index, wavelength) in [450.0, 650.0].into_iter().enumerate() {
            let a = bare
                .project_sequential_point(
                    cx,
                    point,
                    1.0,
                    wavelength,
                    pupil,
                    options([0.01, -0.005]),
                )
                .unwrap();
            let b = coated
                .project_sequential_point(
                    cx,
                    point,
                    1.0,
                    wavelength,
                    pupil,
                    options([0.01, -0.005]),
                )
                .unwrap();
            assert_eq!(a, b);
            projections[index] = a.tangent_xy();
        }
        assert!((projections[0][0] - projections[1][0]).abs() > 1e-9);
        assert!(
            coated
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    pupil,
                    &[550.0],
                )
                .is_err()
        );
        assert_eq!(
            coated.project_sequential_point(cx, point, 1.0, 900.0, pupil, options([0.01, -0.005]),),
            Err(SequentialProjectionError::Camera(
                CameraError::SequentialLensOutsideModel
            ))
        );
    });
}

#[test]
fn inverse_projection_never_publishes_clipped_unresolved_or_cancelled_coordinates() {
    with_cx(|gate, cx| {
        let optical = biconvex(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0));
        let open = camera(
            optical.clone(),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let pupil = LensSample::try_new(0.55, 0.45).unwrap();
        let point = target_on_ray(&open, cx, [0.004, 0.002], pupil, 550.0);
        let mut limited = options([0.004, 0.002]);
        limited.maximum_trace_evaluations = 1;
        assert_eq!(
            open.project_sequential_point(cx, point, 1.0, 550.0, pupil, limited,)
                .unwrap()
                .trace_evaluations(),
            1
        );
        limited.initial_tangent_xy = [0.006, 0.004];
        assert_eq!(
            open.project_sequential_point(cx, point, 1.0, 550.0, pupil, limited,),
            Err(SequentialProjectionError::WorkLimit {
                trace_evaluations: 1
            })
        );
        assert_eq!(
            open.project_sequential_point(
                cx,
                Point3::new(0.0, 0.0, 0.04),
                1.0,
                550.0,
                pupil,
                options([0.0; 2]),
            ),
            Err(SequentialProjectionError::TargetOutsideModel)
        );
        let closed_lens = lens(&[
            optical.surfaces()[0],
            optical.surfaces()[1],
            surface(0.12, 0.0, 1e-12, None),
        ]);
        let closed = camera(
            closed_lens,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        assert_eq!(
            closed.project_sequential_point(cx, point, 1.0, 550.0, pupil, options([0.004, 0.002]),),
            Err(SequentialProjectionError::SeedClipped)
        );
        let mut invalid = options([0.0; 2]);
        invalid.position_tolerance_m = 0.0;
        assert_eq!(
            open.project_sequential_point(cx, point, 1.0, 550.0, pupil, invalid,),
            Err(SequentialProjectionError::InvalidOptions)
        );
        gate.request();
        assert_eq!(
            open.project_sequential_point(cx, point, 1.0, 550.0, pupil, options([0.0; 2]),),
            Err(SequentialProjectionError::Camera(CameraError::Cancelled))
        );
    });
}

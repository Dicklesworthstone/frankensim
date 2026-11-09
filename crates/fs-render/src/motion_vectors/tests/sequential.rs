//! Native wavelength/pupil motion retains the original mesh correspondence.
use super::*;
use crate::camera::{
    LensSurface, SequentialLens, SequentialProjectionError, SequentialProjectionOptions,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};

fn spectral_camera() -> PhysicalCamera {
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.5, 0.004, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.0).unwrap(),
        GlassProvenance::Custom,
    );
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        5.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(
        SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.012, Some(glass)).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}

fn native_miss(
    cx: &Cx<'_>,
    camera: &PhysicalCamera,
    projection: RasterProjection,
    raster: RasterSize,
    wavelength: f64,
    pupil: LensSample,
    target: Point3,
) -> f64 {
    let ndc = projection.ndc_xy();
    let half = camera.raster_vertical_half_tan();
    let ray = camera
        .generate_sequential_spectral_rays_from_tangent_offsets(
            cx,
            ndc[0] * raster.aspect_ratio() * half,
            ndc[1] * half,
            pupil,
            &[wavelength],
        )
        .unwrap()[0]
        .unwrap()
        .ray;
    let delta = target.delta_from(ray.origin);
    let distance = delta.dot(ray.dir) / ray.dir.dot(ray.dir);
    assert!(distance > 0.0);
    let closest = ray.origin.offset(ray.dir.scale(distance));
    target.delta_from(closest).norm()
}

#[test]
fn native_spectral_motion_follows_one_material_point_and_exact_raster_signs() {
    with_cx(|cx| {
        let camera = spectral_camera();
        let raster = RasterSize::try_new(640, 360).unwrap();
        let pupil = LensSample::try_new(0.73, 0.41).unwrap();
        let previous_pose = transform_z(0.0, [-0.99, 0.005, -5.0]);
        let current_pose = transform_z(0.0, [-1.0, 0.0, -5.0]);
        let next_pose = transform_z(0.0, [-1.02, -0.004, -5.0]);
        let primary = hit_sample(cx, &instance(current_pose));
        let previous = frame(0.0, 7, camera.clone(), previous_pose);
        let current = frame(1.0, 7, camera.clone(), current_pose);
        let next = frame(2.0, 7, camera.clone(), next_pose);
        let mut blue_displacement: Option<f64> = None;
        for wavelength in [450.0, 650.0] {
            let result = available(
                &compute_sequential_motion_vectors(
                    cx,
                    primary,
                    &previous,
                    &current,
                    &next,
                    raster,
                    wavelength,
                    pupil,
                    SequentialProjectionOptions::default(),
                )
                .unwrap(),
            );
            assert_eq!(result.identity, primary.identity());
            let (before, before_ndc, before_pixels) = projected(result.previous);
            let (after, after_ndc, after_pixels) = projected(result.next);
            for (projection, pose) in [
                (before, previous_pose),
                (result.current, current_pose),
                (after, next_pose),
            ] {
                let target = pose.transform_point(Point3::new(1.0, 0.0, 0.0));
                assert!(
                    native_miss(cx, &camera, projection, raster, wavelength, pupil, target) < 2e-8
                );
                assert_close(projection.depth_m(), 5.0);
            }
            for (target, ndc_delta, pixel_delta) in [
                (before, before_ndc, before_pixels),
                (after, after_ndc, after_pixels),
            ] {
                assert_close(
                    ndc_delta[0],
                    target.ndc_xy()[0] - result.current.ndc_xy()[0],
                );
                assert_close(
                    ndc_delta[1],
                    target.ndc_xy()[1] - result.current.ndc_xy()[1],
                );
                assert_close(pixel_delta[0], ndc_delta[0] * 320.0);
                assert_close(pixel_delta[1], -ndc_delta[1] * 180.0);
            }
            assert!(before_pixels[0].abs() > 0.1);
            if let Some(blue) = blue_displacement {
                assert!((before_pixels[0] - blue).abs() > 1e-6);
            } else {
                blue_displacement = Some(before_pixels[0]);
            }
        }
        assert!(matches!(
            compute_motion_vectors(primary, &previous, &current, &next, raster),
            Err(MotionVectorError::Camera(
                CameraError::SequentialLensOutsideModel
            ))
        ));
    });
}

#[test]
fn native_motion_preserves_cuts_behind_camera_identity_and_bounded_refusal() {
    with_cx(|cx| {
        let camera = spectral_camera();
        let raster = RasterSize::try_new(640, 360).unwrap();
        let pupil = LensSample::try_new(0.73, 0.41).unwrap();
        let pose = transform_z(0.0, [-1.0, 0.0, -5.0]);
        let primary = hit_sample(cx, &instance(pose));
        let current = frame(1.0, 7, camera.clone(), pose);
        let seed = camera
            .project_sequential_point(
                cx,
                Point3::new(0.0, 0.0, -5.0),
                raster.aspect_ratio(),
                550.0,
                pupil,
                SequentialProjectionOptions::default(),
            )
            .unwrap()
            .tangent_xy();
        let options = SequentialProjectionOptions {
            initial_tangent_xy: seed,
            maximum_trace_evaluations: 1,
            ..SequentialProjectionOptions::default()
        };
        // A cut is categorical even when the other shot has no native lens.
        let cut = frame(0.0, 8, camera_at(0.0, 0.0), pose);
        let result = available(
            &compute_sequential_motion_vectors(
                cx, primary, &cut, &current, &current, raster, 550.0, pupil, options,
            )
            .unwrap(),
        );
        assert_eq!(
            result.previous,
            MotionEndpoint::CameraCut {
                current_shot_id: 7,
                target_shot_id: 8,
            }
        );
        assert_eq!(projected(result.next).2.map(f64::to_bits), [0; 2]);

        let behind = frame(0.0, 7, camera.clone(), transform_z(0.0, [-1.0, 0.0, 5.0]));
        let result = available(
            &compute_sequential_motion_vectors(
                cx, primary, &behind, &current, &current, raster, 550.0, pupil, options,
            )
            .unwrap(),
        );
        assert_eq!(
            result.previous,
            MotionEndpoint::BehindCamera {
                signed_depth_m: -5.0
            }
        );

        let displaced = frame(0.0, 7, camera.clone(), transform_z(0.0, [-0.98, 0.0, -5.0]));
        assert!(matches!(
            compute_sequential_motion_vectors(
                cx, primary, &displaced, &current, &current, raster, 550.0, pupil, options,
            ),
            Err(MotionVectorError::SequentialProjection(
                SequentialProjectionError::WorkLimit {
                    trace_evaluations: 1
                }
            ))
        ));
        let foreign_pose = frame(1.0, 7, camera, transform_z(0.0, [-0.99, 0.0, -5.0]));
        assert!(matches!(
            compute_sequential_motion_vectors(
                cx,
                primary,
                &cut,
                &foreign_pose,
                &foreign_pose,
                raster,
                550.0,
                pupil,
                options,
            ),
            Err(MotionVectorError::ForeignCurrentPose)
        ));
    });
}

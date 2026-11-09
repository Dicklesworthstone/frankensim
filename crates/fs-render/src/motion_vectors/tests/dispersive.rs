//! Original material correspondence under a fixed dispersive pupil.
use super::*;
use crate::camera::{DispersiveThinLens, ThinLensElement};
use crate::dielectric::CauchyIor;

fn dispersive_camera() -> PhysicalCamera {
    let lens = DispersiveThinLens::try_new(
        550.0,
        &[ThinLensElement::try_new(
            CauchyIor::try_new(1.5, 0.004, 0.0).unwrap().into(),
            20.0,
            -20.0,
        )
        .unwrap()],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.005).unwrap(),
        5.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_geometric_dispersive_lens(lens)
    .unwrap()
}

#[test]
fn g1_geometric_lens_motion_replays_exact_correspondence_pupil_and_reference_frames() {
    with_cx(|cx| {
        let camera = dispersive_camera();
        let pupil = LensSample::try_new(0.83, 0.27).unwrap();
        let raster = RasterSize::try_new(640, 360).unwrap();
        let pose = transform_z(0.0, [-1.0, 0.0, -5.0]);
        let primary = hit_sample(cx, &instance(pose));
        let current = frame(1.0, 7, camera.clone(), pose);
        let previous = frame(
            0.0,
            7,
            camera.clone(),
            transform_z(0.0, [-0.98, 0.01, -6.0]),
        );
        for wavelength in [450.0, 550.0, 650.0] {
            let motion = available(
                &compute_geometric_dispersive_motion_vectors(
                    cx, primary, &previous, &current, &current, raster, wavelength, pupil, 1e-9,
                )
                .unwrap(),
            );
            let (target, ndc, pixels) = projected(motion.previous);
            for (projection, world) in [
                (motion.current, Point3::new(0.0, 0.0, -5.0)),
                (target, Point3::new(0.02, 0.01, -6.0)),
            ] {
                let n = projection.ndc_xy();
                let half = camera.raster_vertical_half_tan();
                let [ray] = camera
                    .generate_geometric_spectral_rays_from_tangent_offsets(
                        cx,
                        n[0] * raster.aspect_ratio() * half,
                        n[1] * half,
                        pupil,
                        &[wavelength],
                    )
                    .unwrap();
                let delta = world.delta_from(ray.origin);
                let distance = delta.dot(ray.dir) / ray.dir.dot(ray.dir);
                assert!(
                    world
                        .delta_from(ray.origin.offset(ray.dir.scale(distance)))
                        .norm()
                        < 2e-9
                );
            }
            assert_close(pixels[0], ndc[0] * 320.0);
            assert_close(pixels[1], -ndc[1] * 180.0);
            assert_eq!(projected(motion.next).2.map(f64::to_bits), [0; 2]);
        }
        let cut = frame(0.0, 8, camera_at(0.0, 0.0), pose);
        let result = available(
            &compute_geometric_dispersive_motion_vectors(
                cx, primary, &cut, &current, &current, raster, 550.0, pupil, 1e-9,
            )
            .unwrap(),
        );
        assert!(matches!(result.previous, MotionEndpoint::CameraCut { .. }));
        assert!(matches!(
            compute_geometric_dispersive_motion_vectors(
                cx, primary, &previous, &current, &current, raster, 550.0, pupil, 0.0,
            ),
            Err(MotionVectorError::Camera(CameraError::InvalidProjection))
        ));
    });
}

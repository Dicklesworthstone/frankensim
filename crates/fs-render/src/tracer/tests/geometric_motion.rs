//! Full cinematic consumption of native material-dispersive geometric rays.
use super::sequential_motion::{previous, provenance, scene, settings, shutter, trace};
use super::*;
use crate::camera::{
    CameraKeyframe, CameraProjection, CameraShot, DispersiveThinLens, ThinLensElement,
};
use crate::dielectric::DielectricIor;

fn physical(eye_x: f64, exterior: Option<DielectricIor>) -> PhysicalCamera {
    let element = ThinLensElement::try_new(
        CauchyIor::try_new(1.7, 0.004, 0.0).unwrap().into(),
        20.0,
        -20.0,
    )
    .unwrap();
    let lens = match exterior {
        Some(exterior) => DispersiveThinLens::try_new_in_medium(550.0, &[element], exterior),
        None => DispersiveThinLens::try_new(550.0, &[element]),
    }
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(eye_x, 0.0, 3.0),
        Point3::new(eye_x, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
        3.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_geometric_dispersive_lens(lens)
    .unwrap()
}

fn animated(exterior: Option<DielectricIor>) -> AnimatedCamera {
    AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            7,
            0.0,
            1.0,
            vec![
                CameraKeyframe::try_new(0.0, physical(0.0, exterior)).unwrap(),
                CameraKeyframe::try_new(1.0, physical(0.01, exterior)).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap()
}

#[test]
fn every_native_wavelength_hits_its_own_emitter_and_retains_the_actual_hero_witness() {
    with_test_cx(|cx| {
        let camera = animated(None);
        let exposure = camera.admit_shutter(cx, shutter(), CutSide::After).unwrap();
        let time = shutter().time_at(shutter().sample_for_stream(settings().seed, 17, 29));
        let physical = camera.evaluate_exposure(cx, exposure, time).unwrap();
        let half = physical.raster_vertical_half_tan();
        let tangent = [
            (2.0 * (3.0 + 0.42) / 7.0 - 1.0) * settings().camera_aspect() * half,
            (1.0 - 2.0 * (2.0 + 0.58) / 5.0) * half,
        ];
        let wavelengths = hero_wavelength_packet::<PACKET>(
            LAMBDA_MIN + 0.37 * (LAMBDA_MAX - LAMBDA_MIN),
            LAMBDA_MIN,
            LAMBDA_MAX,
        );
        let pupil = camera_lens_sample([0x89ab_cdef, 0x0123_4567], 17, 29).unwrap();
        let rays = physical
            .generate_geometric_spectral_rays_from_tangent_offsets(
                cx,
                tangent[0],
                tangent[1],
                pupil,
                &wavelengths,
            )
            .unwrap();
        let targets = rays.map(|ray| ray.origin.offset(ray.dir.scale(-ray.origin.z / ray.dir.z)));
        let mut separation = f64::INFINITY;
        for i in 0..PACKET {
            for j in 0..i {
                separation = separation.min(targets[i].delta_from(targets[j]).norm());
            }
        }
        assert!(
            separation > 1e-7,
            "fixture must resolve the actual chromatic directions"
        );
        let half_extent = separation / 8.0;
        let white = lift_rgb([1.0; 3]);
        let black = lift_rgb([0.0; 3]);
        let strengths = [1.0, 2.0, 4.0, 8.0];
        let mut scene = scene();
        scene.primitives.clear();
        scene.lights.clear();
        for (lane, point) in targets.into_iter().enumerate() {
            scene.primitives.push(Primitive {
                shape: instanced_test_mesh(
                    horizontal_quad_at(0.0, point.x, point.y, half_extent),
                    60 + lane as u64,
                ),
                material: Material::Lambertian { reflectance: black },
                emission: Some((white, strengths[lane])),
            });
            scene.lights.push(RectLight {
                corner: Point3::new(point.x - half_extent, point.y - half_extent, 0.0),
                edge_u: Vec3::new(2.0 * half_extent, 0.0, 0.0),
                edge_v: Vec3::new(0.0, 2.0 * half_extent, 0.0),
                prim: lane,
                emission: (white, strengths[lane]),
            });
        }
        let actual = trace(&scene, &camera, cx, true);
        let uncaptured = trace(&scene, &camera, cx, false);
        assert_eq!(
            actual.xyz.map(f64::to_bits),
            uncaptured.xyz.map(f64::to_bits)
        );
        assert!(uncaptured.spectral_motion.is_none());
        let mut expected = [0.0; 3];
        for lane in 0..PACKET {
            let wavelength = wavelengths[lane];
            let weight = white.eval(wavelength) * strengths[lane] * (LAMBDA_MAX - LAMBDA_MIN)
                / PACKET as f64
                / y_integral();
            expected[0] += weight * cie_x(wavelength);
            expected[1] += weight * cie_y(wavelength);
            expected[2] += weight * cie_z(wavelength);
        }
        for channel in 0..3 {
            assert_near(
                actual.xyz[channel],
                expected[channel],
                3e-12,
                "native four-ray emission integral",
            );
        }
        assert_eq!(actual.primary.as_ref().unwrap().primitive_index, 0);
        assert!(
            actual
                .primary
                .as_ref()
                .unwrap()
                .hit
                .point
                .delta_from(targets[0])
                .norm()
                < 1e-12
        );
        let witness = actual.spectral_motion.unwrap();
        assert_eq!(witness.wavelength_nm.to_bits(), wavelengths[0].to_bits());
        assert_eq!(witness.lens_sample, pupil);
        assert_eq!(
            witness.initial_tangent_xy.map(f64::to_bits),
            tangent.map(f64::to_bits)
        );
        assert!(matches!(
            previous(cx, &camera, &scene, &actual, None),
            Err(CinematicAovError::SampleAlignmentMismatch)
        ));
        let displacement = previous(cx, &camera, &scene, &actual, Some(witness))
            .unwrap()
            .unwrap();
        // Original translated camera, fixed world point and common raster
        // magnification: pupil chromatic offset cancels in this translation.
        let expected_x = (0.01 * 0.25) / 3.0 / (settings().camera_aspect() * half) * 3.5;
        assert_near(
            displacement[0],
            expected_x,
            2e-10,
            "fixed-pupil previous motion",
        );
        assert_near(displacement[1], 0.0, 2e-10, "unchanged vertical position");
    });
}

#[test]
fn public_geometric_aov_film_keeps_beauty_and_rejects_wrong_exterior_or_field() {
    with_test_cx(|cx| {
        let scene = scene();
        let camera = animated(None);
        let s = Settings {
            width: 3,
            height: 3,
            spp: 2,
            ..settings()
        };
        let config = CinematicAovConfig::new(
            CinematicAovProfile::FinalDiagnostic,
            provenance(),
            CinematicAovLimits::default(),
        );
        let film =
            render_cinematic_with_aovs(&scene, &camera, CutSide::After, cx, &s, shutter(), config)
                .unwrap();
        let beauty = render_cinematic(&scene, &camera, CutSide::After, cx, &s, shutter()).unwrap();
        for (actual, expected) in film.beauty().xyz.iter().zip(&beauty.xyz) {
            assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        }
        let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
        let validity = &decoded
            .channels
            .iter()
            .find(|c| c.name == "diagnostic.validity")
            .unwrap()
            .data;
        let guides = film.denoise_guides().unwrap();
        let mut valid = 0;
        for (index, &mask) in validity.iter().enumerate() {
            if (mask as u32) & crate::aov::validity::PREVIOUS_MOTION != 0 {
                valid += 1;
                assert!(guides.motion_prev_x()[index] > 1e-5);
                assert!((guides.axial_depth_m()[index] - 3.0).abs() < 1e-5);
            }
        }
        assert!(valid > 0);
        let wrong_medium = animated(Some(CauchyIor::try_constant(1.33).unwrap().into()));
        assert!(matches!(
            render_cinematic(&scene, &wrong_medium, CutSide::After, cx, &s, shutter()),
            Err(TracerError::Camera(CameraError::DispersiveLensOutsideModel))
        ));
        let wide = Settings {
            width: 100,
            height: 1,
            spp: 1,
            ..s
        };
        assert!(matches!(
            render_cinematic(&scene, &camera, CutSide::After, cx, &wide, shutter()),
            Err(TracerError::Camera(CameraError::DispersiveLensOutsideModel))
        ));
        // Existing shared-ray API still refuses distinct wavelength geometry.
        assert_eq!(
            physical(0.0, None)
                .generate_spectral_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    LensSample::CENTER,
                    &[550.0],
                )
                .unwrap_err(),
            CameraError::DispersionRequiresSpectralRays
        );
    });
}

//! Actual beauty-ray witnesses drive source-bound sequential-lens motion AOVs.
use super::*;
use crate::camera::{CameraKeyframe, CameraProjection, CameraShot, LensSurface, SequentialLens};
use crate::dielectric::GlassProvenance;

fn physical(eye_x: f64, sequential: bool, front_radius: f64) -> PhysicalCamera {
    let camera = PhysicalCamera::try_look_at(
        Point3::new(eye_x, 0.0, 3.0),
        Point3::new(eye_x, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        3.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap();
    if !sequential {
        return camera;
    }
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.5, 0.004, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.0).unwrap(),
        GlassProvenance::Custom,
    );
    camera
        .with_sequential_lens(
            SequentialLens::try_new(
                550.0,
                &[
                    LensSurface::try_new(0.08, 20.0, 0.012, Some(glass)).unwrap(),
                    LensSurface::try_new(0.09, -20.0, front_radius, None).unwrap(),
                ],
            )
            .unwrap(),
        )
        .unwrap()
}

fn animated(sequential: bool, front_radius: f64) -> AnimatedCamera {
    AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            7,
            0.0,
            1.0,
            vec![
                CameraKeyframe::try_new(0.0, physical(0.0, sequential, front_radius)).unwrap(),
                CameraKeyframe::try_new(1.0, physical(0.01, sequential, front_radius)).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap()
}

pub(super) fn scene() -> Scene {
    let white = lift_rgb([1.0; 3]);
    Scene {
        primitives: vec![Primitive {
            shape: instanced_test_mesh(horizontal_quad(0.0, 4.0), 41),
            material: Material::Lambertian { reflectance: white },
            emission: Some((white, 2.0)),
        }],
        lights: vec![RectLight {
            corner: Point3::new(-4.0, -4.0, 0.0),
            edge_u: Vec3::new(8.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 8.0, 0.0),
            prim: 0,
            emission: (white, 2.0),
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 3.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.24,
        },
    }
}

pub(super) fn settings() -> Settings {
    Settings {
        width: 7,
        height: 5,
        spp: 32,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x0123_4567_89ab_cdef,
    }
}

pub(super) fn shutter() -> ShutterInterval {
    ShutterInterval::resolve(
        0.55,
        0.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 8 },
        ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
    )
    .unwrap()
}

pub(super) fn provenance() -> CinematicAovProvenance {
    let identity = hash_domain(
        "org.frankensim.render.sequential-motion-aov-test",
        b"native",
    );
    CinematicAovProvenance::try_new(1, 0.5, 0.25, 0.75, identity, identity, identity).unwrap()
}

pub(super) fn trace(scene: &Scene, camera: &AnimatedCamera, cx: &Cx<'_>, capture: bool) -> PathTraceSample {
    let interval = shutter();
    let time = PathTime {
        interval,
        normalized: interval.sample_for_stream(settings().seed, 17, 29),
        cached_animated: std::array::from_fn(|_| None),
    };
    let exposure = camera.admit_shutter(cx, interval, CutSide::After).unwrap();
    let lighting = SceneLighting::from(AdmittedLighting::try_new(&scene.lights, None).unwrap());
    trace_path(
        scene,
        &lighting,
        cx,
        &settings(),
        1.0 / y_integral(),
        17,
        29,
        0.42,
        0.58,
        0.37,
        Some(&time),
        CameraPath::Cinematic { camera, exposure },
        capture,
        None,
        None,
    )
    .unwrap()
}

pub(super) fn previous(
    cx: &Cx<'_>,
    camera: &AnimatedCamera,
    scene: &Scene,
    traced: &PathTraceSample,
    witness: Option<SpectralMotionWitness>,
) -> Result<Option<[f64; 2]>, CinematicAovError> {
    let exposure = camera.admit_shutter(cx, shutter(), CutSide::After).unwrap();
    let time = traced.absolute_time_s.unwrap();
    aligned_previous_motion(
        cx,
        camera,
        exposure,
        CutSide::After,
        &scene.primitives[0].shape,
        traced.primary.as_ref().unwrap().surface.unwrap(),
        provenance(),
        &settings(),
        time,
        camera.evaluate_exposure(cx, exposure, time).unwrap(),
        witness,
    )
}

#[test]
fn original_hero_pupil_and_sensor_tangent_drive_aligned_motion_without_changing_beauty() {
    with_test_cx(|cx| {
        let scene = scene();
        let camera = animated(true, 0.012);
        let traced = trace(&scene, &camera, cx, true);
        let uncaptured = trace(&scene, &camera, cx, false);
        assert_eq!(
            traced.xyz.map(f64::to_bits),
            uncaptured.xyz.map(f64::to_bits)
        );
        assert!(uncaptured.spectral_motion.is_none());
        let witness = traced.spectral_motion.unwrap();
        let primary = traced.primary.as_ref().unwrap();
        let time = traced.absolute_time_s.unwrap();
        let exposure = camera.admit_shutter(cx, shutter(), CutSide::After).unwrap();
        let physical = camera.evaluate_exposure(cx, exposure, time).unwrap();
        let half = physical.raster_vertical_half_tan();
        let expected_tangent = [
            (2.0 * (3.0 + 0.42) / 7.0 - 1.0) * settings().camera_aspect() * half,
            (1.0 - 2.0 * (2.0 + 0.58) / 5.0) * half,
        ];
        assert_eq!(
            witness.initial_tangent_xy.map(f64::to_bits),
            expected_tangent.map(f64::to_bits)
        );
        assert_eq!(
            witness.wavelength_nm.to_bits(),
            hero_wavelength_packet::<PACKET>(
                LAMBDA_MIN + 0.37 * (LAMBDA_MAX - LAMBDA_MIN),
                LAMBDA_MIN,
                LAMBDA_MAX,
            )[0]
            .to_bits()
        );
        assert_eq!(
            witness.lens_sample,
            camera_lens_sample([0x89ab_cdef, 0x0123_4567], 17, 29).unwrap()
        );
        let ray = physical
            .generate_sequential_spectral_rays_from_tangent_offsets(
                cx,
                expected_tangent[0],
                expected_tangent[1],
                witness.lens_sample,
                &[witness.wavelength_nm],
            )
            .unwrap()[0]
            .unwrap()
            .ray;
        let distance = -ray.origin.z / ray.dir.z;
        let expected_hit = ray.origin.offset(ray.dir.scale(distance));
        assert!(primary.hit.point.delta_from(expected_hit).norm() < 1e-12);

        // The reference keeps the accepted 0.55 shutter phase: 0.55 - 0.25,
        // rather than silently projecting at the previous presentation time.
        let before = motion_frame_for_shape(
            cx,
            &camera,
            CutSide::After,
            &scene.primitives[0].shape,
            time - 0.25,
            None,
        )
        .unwrap()
        .unwrap();
        let current = motion_frame_for_shape(
            cx,
            &camera,
            CutSide::After,
            &scene.primitives[0].shape,
            time,
            Some((exposure.shot_id(), physical.clone())),
        )
        .unwrap()
        .unwrap();
        let MotionVectorComputation::Available(reference) = compute_sequential_motion_vectors(
            cx,
            primary.surface.unwrap(),
            &before,
            &current,
            &current,
            RasterSize::try_new(7, 5).unwrap(),
            witness.wavelength_nm,
            witness.lens_sample,
            SequentialProjectionOptions {
                initial_tangent_xy: expected_tangent,
                ..Default::default()
            },
        )
        .unwrap() else {
            panic!("native correspondence unavailable")
        };
        let MotionEndpoint::Projected {
            displacement_pixels: expected,
            ..
        } = reference.previous
        else {
            panic!("native previous projection unavailable")
        };
        let actual = previous(cx, &camera, &scene, &traced, Some(witness))
            .unwrap()
            .unwrap();
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert!(actual[0].abs() > 1e-5);

        let palette =
            CinematicAovPalette::try_from_scene(&scene, CinematicAovLimits::default(), true, cx)
                .unwrap();
        let cache = AovAlbedoCache::try_new(scene.primitives.len(), true).unwrap();
        let aligned = prepare_aligned_aov_primary(
            &scene,
            &camera,
            exposure,
            CutSide::After,
            cx,
            &settings(),
            provenance(),
            &palette,
            &cache,
            true,
            time,
            Some(primary),
            Some(witness),
        )
        .unwrap()
        .unwrap();
        assert_eq!(aligned.previous_motion_pixels, Some(actual));
        assert_near(aligned.depth_m, 3.0, 1e-12, "depth uses the camera eye");
    });
}

#[test]
fn native_clipping_missing_witness_local_inverse_and_camera_cut_have_distinct_outcomes() {
    with_test_cx(|cx| {
        let scene = scene();
        let camera = animated(true, 0.012);
        let traced = trace(&scene, &camera, cx, true);
        let witness = traced.spectral_motion.unwrap();
        assert!(matches!(
            previous(cx, &camera, &scene, &traced, None),
            Err(CinematicAovError::SampleAlignmentMismatch)
        ));
        let clipped_seed = SpectralMotionWitness {
            initial_tangent_xy: [1e6, 1e6],
            ..witness
        };
        assert_eq!(
            previous(cx, &camera, &scene, &traced, Some(clipped_seed)).unwrap(),
            None
        );
        let invalid = SpectralMotionWitness {
            initial_tangent_xy: [f64::NAN, 0.0],
            ..witness
        };
        assert!(matches!(
            previous(cx, &camera, &scene, &traced, Some(invalid)),
            Err(CinematicAovError::Tracer(TracerError::MotionVector(
                MotionVectorError::SequentialProjection(SequentialProjectionError::InvalidOptions)
            )))
        ));
        let cut = AnimatedCamera::try_new(vec![
            CameraShot::try_new(
                8,
                0.0,
                0.4,
                vec![CameraKeyframe::try_new(0.0, physical(0.0, false, 0.012)).unwrap()],
            )
            .unwrap(),
            CameraShot::try_new(
                7,
                0.4,
                1.0,
                vec![CameraKeyframe::try_new(0.4, physical(0.0055, true, 0.012)).unwrap()],
            )
            .unwrap(),
        ])
        .unwrap();
        let at_cut = trace(&scene, &cut, cx, true);
        assert_eq!(
            previous(cx, &cut, &scene, &at_cut, at_cut.spectral_motion).unwrap(),
            None
        );

        // A physically tiny front stop clips the actual packet. AOV geometry
        // stays absent; another wavelength must never supply a hero witness.
        let stopped = trace(&scene, &animated(true, 1e-9), cx, true);
        assert!(stopped.primary.is_none());
        assert!(stopped.spectral_motion.is_none());
        assert_eq!(stopped.xyz, [0.0; 3]);
        assert!(matches!(
            TracerError::from(MotionVectorError::from(SequentialProjectionError::Camera(
                CameraError::Cancelled
            ))),
            TracerError::Cancelled
        ));
    });
}

#[test]
fn ordinary_camera_previous_motion_retains_the_original_projection_bits() {
    with_test_cx(|cx| {
        let scene = scene();
        let camera = animated(false, 0.012);
        let traced = trace(&scene, &camera, cx, true);
        assert!(traced.spectral_motion.is_none());
        let time = traced.absolute_time_s.unwrap();
        let before = motion_frame_for_shape(
            cx,
            &camera,
            CutSide::After,
            &scene.primitives[0].shape,
            time - 0.25,
            None,
        )
        .unwrap()
        .unwrap();
        let current = motion_frame_for_shape(
            cx,
            &camera,
            CutSide::After,
            &scene.primitives[0].shape,
            time,
            None,
        )
        .unwrap()
        .unwrap();
        let MotionVectorComputation::Available(reference) = compute_motion_vectors(
            traced.primary.as_ref().unwrap().surface.unwrap(),
            &before,
            &current,
            &current,
            RasterSize::try_new(7, 5).unwrap(),
        )
        .unwrap() else {
            panic!("ordinary correspondence unavailable")
        };
        let MotionEndpoint::Projected {
            displacement_pixels,
            ..
        } = reference.previous
        else {
            panic!("ordinary previous projection unavailable")
        };
        assert_eq!(
            previous(cx, &camera, &scene, &traced, None)
                .unwrap()
                .unwrap()
                .map(f64::to_bits),
            displacement_pixels.map(f64::to_bits)
        );
    });
}

#[test]
fn public_sequential_aov_film_exports_real_motion_and_keeps_beauty_bits() {
    with_test_cx(|cx| {
        let scene = scene();
        let camera = animated(true, 0.012);
        let settings = Settings {
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
        let film = render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            config,
        )
        .unwrap();
        let beauty =
            render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter()).unwrap();
        for (actual, expected) in film.beauty().xyz.iter().zip(&beauty.xyz) {
            assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        }
        let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
        let validity = &decoded
            .channels
            .iter()
            .find(|channel| channel.name == "diagnostic.validity")
            .unwrap()
            .data;
        let motion = film.denoise_guides().unwrap();
        assert!(
            validity
                .iter()
                .any(|value| (*value as u32) & crate::aov::validity::PREVIOUS_MOTION != 0)
        );
        assert!(
            motion
                .motion_prev_x()
                .iter()
                .any(|value| value.abs() > 1e-5)
        );
        for (index, value) in validity.iter().enumerate() {
            if (*value as u32) & crate::aov::validity::PREVIOUS_MOTION != 0 {
                assert_eq!(
                    (*value as u32) & crate::aov::validity::PRIMARY,
                    crate::aov::validity::PRIMARY
                );
                assert!(motion.motion_prev_x()[index].is_finite());
                assert!((motion.axial_depth_m()[index] - 3.0).abs() < 1e-5);
            }
        }
    });
}

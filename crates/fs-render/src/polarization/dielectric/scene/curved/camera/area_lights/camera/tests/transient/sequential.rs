//! Real-lens group clocks on both native finite-area MIS strategies.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;

fn physical_camera() -> PhysicalCamera {
    let material = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(
        SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(material)).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
                LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}

fn optical(camera: &PhysicalCamera, cx: &Cx<'_>, wavelength: f64) -> SequentialCameraBeam {
    camera
        .generate_sequential_spectral_timed_beams_from_tangent_offsets(
            cx,
            0.0,
            0.0,
            LensSample::CENTER,
            &[wavelength],
        )
        .unwrap()[0]
        .unwrap()
}

fn binding<'a, 's, 'm>(
    source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    seeds: &'a [BundleImageSeed<'a>],
) -> AreaLightBinding<'a, 's, 'm> {
    AreaLightBinding {
        source,
        receiver_object_id: 1,
        source_side: BoundarySide::Positive,
        initial_medium: 0,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        seeds,
        image_options: image_options(1),
        position_tolerance_m: 2e-7,
    }
}

#[test]
fn resolved_lens_prefix_preserves_both_mis_strategies_and_sensor_epoch() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let vacuum = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [binding(&light, &seeds)];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let history = History::new(|_| Ok(1.0));
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let physical = physical_camera();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(20e-9);
        assert_eq!(clock.to_f64(), epoch.to_f64());

        for wavelength in [450.0, 650.0] {
            let beam = optical(&physical, cx, wavelength);
            let baseline = lights
                .sample_with_modulation(
                    cx,
                    || camera.launch_sequential(&beam, 0),
                    0,
                    settings(),
                    &mut workspace,
                    |_| [0.0; 3],
                    |_, _| AREA_DRAW,
                    |_| Ok(1.0),
                )
                .unwrap();
            let mut surface_draws = Vec::new();
            let mut area_draws = Vec::new();
            let result = observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &beam,
                    0,
                    clock,
                    settings(),
                    &mut workspace,
                    |event| {
                        surface_draws.push(event);
                        [0.0; 3]
                    },
                    |event, binding| {
                        area_draws.push((event, binding));
                        AREA_DRAW
                    },
                )
                .unwrap();
            assert_eq!(result.area_sample(), &baseline);
            assert!(baseline.area_radiance_w_m2_sr_nm() > 0.0);
            assert!(
                baseline.radiance_sample().radiance_w_m2_sr_nm()
                    > baseline.area_radiance_w_m2_sr_nm()
            );
            assert_eq!(surface_draws, [0]);
            assert_eq!(area_draws, [(0, 0)]);
            assert_eq!(result.observation_time_s(), clock);

            let queries = history.queries.lock().unwrap();
            let pair = &queries[queries.len() - 2..];
            // Independent axial path: 9 cm air + 1 cm Cauchy glass in
            // the lens, then 0.9 m to the diffuse panel + 2 m to the source.
            let ng = 1.48 + 0.018 / (wavelength * 1e-3).powi(2);
            let expected_delay = (2.99 + 0.01 * ng) / C;
            for query in pair {
                assert_eq!(query.observation_time_s, clock);
                assert_eq!(query.wavelength_nm, wavelength);
                near(query.group_delay_s.to_f64(), expected_delay, 2e-8);
                near(
                    (query.emission_time_s - epoch).to_f64(),
                    20e-9 - expected_delay,
                    2e-8,
                );
            }
            let camera_hit = pair[1];
            let expected_hit_delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(
                    baseline
                        .radiance_sample()
                        .path()
                        .group_optical_path_length_m()
                        .unwrap(),
                ) / Dd::from_f64(C);
            assert_eq!(camera_hit.group_delay_s, expected_hit_delay);
            assert_eq!(
                pair[0].position_m,
                light.sample_area(AREA_DRAW).unwrap().position_m()
            );
        }
        assert_eq!(history.queries.lock().unwrap().len(), 4);
    });
}

#[test]
fn resolved_lens_shutter_retards_the_complete_pulse_once_and_checks_raw_beams() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let vacuum = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [binding(&light, &seeds)];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let physical = physical_camera();
        let beam = optical(&physical, cx, 550.0);
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = lights
            .sample_with_modulation(
                cx,
                || camera.launch_sequential(&beam, 0),
                0,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
                |_| Ok(1.0),
            )
            .unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let pulse = PiecewiseLinearHistory::try_new(
            cx,
            epoch,
            &[
                (-1e-6, 0.0),
                (0.0, 0.0),
                (1e-9, 2.0),
                (2e-9, 0.0),
                (1e-6, 0.0),
            ],
        )
        .unwrap();
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&pulse)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let ng = 1.48 + 0.018 / 0.55_f64.powi(2);
        let arrival = epoch + Dd::from_f64((2.99 + 0.01 * ng) / C);
        let window = ExposureWindow::try_new(
            arrival - Dd::from_f64(0.5e-9),
            arrival + Dd::from_f64(2.5e-9),
        )
        .unwrap();
        let result = observer
            .sample_sequential_ray_with_samples_during_exposure(
                cx,
                &beam,
                0,
                window,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
            )
            .unwrap();
        near(
            result.mean_radiance_sample().radiance_w_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * 2.0 / 3.0,
            2e-12,
        );
        near(
            result.mean_area_sample().area_radiance_w_m2_sr_nm(),
            baseline.area_radiance_w_m2_sr_nm() * 2.0 / 3.0,
            2e-12,
        );
        near(
            result.integrated_radiance_j_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * 2e-9,
            2e-12,
        );
        assert_eq!(result.mean_area_sample().images(), baseline.images());
        assert_eq!(
            result.mean_area_sample().trace_evaluations(),
            baseline.trace_evaluations()
        );
        assert_eq!(
            result.mean_radiance_sample().path(),
            baseline.radiance_sample().path()
        );

        let untimed = physical
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(
            observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &untimed,
                    0,
                    arrival,
                    settings(),
                    &mut workspace,
                    |_| panic!("untimed beam must be refused before transport"),
                    |_, _| panic!("untimed beam must be refused before NEE"),
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::MissingLensGroupDelay)
        );
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &[])
            .unwrap()
            .with_analyzer(CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap());
        let analyzer_observer =
            TransientCurvedAreaLightCameraScene::try_new(cx, &analyzed, &lights).unwrap();
        assert_eq!(
            analyzer_observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &beam,
                    0,
                    arrival,
                    settings(),
                    &mut workspace,
                    |_| panic!("analyzer must be installed at the sensor"),
                    |_, _| AREA_DRAW,
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::AnalyzerMustPrecedeLens)
        );

        let failure = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let failed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&failure)]).unwrap();
        let failed_observer =
            TransientCurvedAreaLightCameraScene::try_new(cx, &failed, &lights).unwrap();
        assert_eq!(
            failed_observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &beam,
                    0,
                    arrival,
                    settings(),
                    &mut workspace,
                    |_| [0.0; 3],
                    |_, _| AREA_DRAW,
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::OutsideTimeDomain,
            })
        );
        assert_eq!(failure.queries.lock().unwrap().len(), 1);
    });
}

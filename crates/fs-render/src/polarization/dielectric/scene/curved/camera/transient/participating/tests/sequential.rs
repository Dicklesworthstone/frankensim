use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;

mod exposure;

fn lens(exterior: DielectricGlass) -> PhysicalCamera {
    let interior = DielectricGlass::new(
        CauchyIor::try_new(1.9, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.7).unwrap(),
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
        SequentialLens::try_new_with_sensor_medium(
            550.0,
            None,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(interior)).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, Some(exterior)).unwrap(),
                LensSurface::try_new(0.1, 0.0, 0.02, Some(exterior)).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn timed_immersed_lens_keeps_two_molecular_turns_and_chromatic_sensor_clock() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let physical = lens(glass);
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let unity = History::new(|_| Ok(1.0));
        let static_camera =
            TransientCurvedCameraScene::try_new(cx, &camera, &[unity.binding()]).unwrap();
        let static_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &static_camera, &scattering).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let pulse = History::new(|query: SurfaceEmissionQuery| {
            Ok((1.0 - ((query.emission_time_s - epoch).to_f64() / 1e-10).abs()).max(0.0))
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[pulse.binding()]).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let blue_group = 1.4 + 0.06 / 0.45_f64.powi(2);
        let blue_lens =
            (0.08 + 0.01 * (1.9 + 0.018 / 0.45_f64.powi(2)) + 0.01 * blue_group) / LIGHT_M_S;
        let clock = epoch + Dd::from_f64(blue_lens + 1.4 * blue_group / LIGHT_M_S);
        for wavelength in [450.0, 650.0] {
            let optical = physical
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[wavelength],
                )
                .unwrap()[0]
                .unwrap();
            let expected = steady
                .sample_sequential_ray_with_samples(
                    cx,
                    &optical,
                    0,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            let static_sample = static_scene
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &optical,
                    0,
                    clock,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            assert_eq!(static_sample.radiance_sample(), expected.radiance_sample());
            assert_eq!(static_sample.volume_events(), 2);
            assert_eq!(static_sample.free_flight_draws(), 3);
            // From the 0.1 m front stop: 0.25 m forward, 0.5 m sideways,
            // then 0.65 m to the source. The lens prefix is not in this path.
            near(expected.radiance_sample().path().traveled_m(), 1.4, 2e-7);
            let observed = scene
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &optical,
                    0,
                    clock,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            let expected_radiance = if wavelength == 450.0 {
                expected.radiance_sample().radiance_w_m2_sr_nm()
            } else {
                0.0
            };
            near(
                observed.radiance_sample().radiance_w_m2_sr_nm(),
                expected_radiance,
                4e-5,
            );
            assert_eq!(observed.observation_time_s(), clock);
            let query = *pulse.queries.lock().unwrap().last().unwrap();
            let group = 1.4 + 0.06 / (wavelength * 1e-3).powi(2);
            let lens_delay =
                (0.08 + 0.01 * (1.9 + 0.018 / (wavelength * 1e-3).powi(2)) + 0.01 * group)
                    / LIGHT_M_S;
            near(optical.group_delay_s().unwrap().to_f64(), lens_delay, 1e-23);
            let total = optical.group_delay_s().unwrap()
                + Dd::from_f64(
                    expected
                        .radiance_sample()
                        .path()
                        .group_optical_path_length_m()
                        .unwrap(),
                ) / Dd::from_f64(LIGHT_M_S);
            assert_eq!(query.observation_time_s, clock);
            assert_eq!(query.group_delay_s, total);
            assert_eq!(query.emission_time_s, clock - total);
            near(
                query.group_delay_s.to_f64(),
                lens_delay + 1.4 * group / LIGHT_M_S,
                3e-16,
            );
        }
        assert_eq!(pulse.queries.lock().unwrap().len(), 2);
    });
}

#[test]
fn timed_lens_null_candidates_preserve_curved_delay_and_zero_scattering_identity() {
    with_cx(|cx, gate| {
        let exterior = DielectricGlass::new(
            CauchyIor::try_constant(1.2).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let physical = lens(exterior);
        let optical = physical
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        let geometry = geometry();
        let media = [SceneMedium::GradientIndex {
            field: &AffineIndex,
            identity: ContentHash([74; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Ok(0.75));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history.binding()]).unwrap();
        let clock = Dd::from_f64(1e-8);
        let volumes = [None];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let clear = transient
            .sample_sequential_ray_at_time(cx, &optical, 0, clock, settings().surfaces, 31, 42)
            .unwrap();
        let zero = scene
            .sample_sequential_ray_at_time(cx, &optical, 0, clock, settings(), 31, 42)
            .unwrap();
        assert_eq!(zero.transient_sample(), &clear);
        assert_eq!(zero.volume_events(), 0);
        assert_eq!(zero.free_flight_draws(), 0);
        let ordinary = physical
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
            scene
                .sample_sequential_ray_at_time(cx, &ordinary, 0, clock, settings(), 31, 42)
                .unwrap_err(),
            TransientCameraError::MissingLensGroupDelay
        );
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &[history.binding()])
            .unwrap()
            .with_analyzer(analyzer);
        let analyzed_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &analyzed, &scattering).unwrap();
        assert_eq!(
            analyzed_scene
                .sample_sequential_ray_at_time(cx, &optical, 0, clock, settings(), 31, 42)
                .unwrap_err(),
            TransientCameraError::AnalyzerMustPrecedeLens
        );
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &sigma,
                &VacuumDensity,
                1.0,
                ContentHash([79; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let before = history.queries.lock().unwrap().len();
        let observed = scene
            .sample_sequential_ray_with_tracking_samples_at_time(
                cx,
                &optical,
                0,
                clock,
                settings(),
                |_| panic!("opaque"),
                |flight| match flight {
                    0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.5],
                    1 => [0.99, 0.5, 0.0, 0.5],
                    _ => panic!("source precedes next candidate"),
                },
            )
            .unwrap();
        assert_eq!(observed.volume_events(), 0);
        assert_eq!(observed.free_flight_draws(), 2);
        assert_eq!(history.queries.lock().unwrap().len() - before, 1);
        // The curved leg spans z=0.1..1.0. Native n=1.2+0.2x gives
        // n=1.2 cosh((z-0.1)/6), ds/dz=cosh((z-0.1)/6).
        let arc = 6.0 * (0.9_f64 / 6.0).sinh();
        let group_m = 0.6 * (0.9 + 3.0 * (0.9_f64 / 3.0).sinh());
        near(observed.radiance_sample().path().traveled_m(), arc, 2e-7);
        let query = *history.queries.lock().unwrap().last().unwrap();
        near(
            query.group_delay_s.to_f64(),
            optical.group_delay_s().unwrap().to_f64() + group_m / LIGHT_M_S,
            3e-16,
        );
        assert_eq!(query.observation_time_s, clock);
        let failure = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let failed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[failure.binding()]).unwrap();
        let failed_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &failed, &scattering).unwrap();
        assert!(matches!(
            failed_scene.sample_sequential_ray_with_tracking_samples_at_time(
                cx,
                &optical,
                0,
                clock,
                settings(),
                |_| panic!("opaque"),
                |_| [0.99, 0.5, 0.0, 0.5]
            ),
            Err(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::OutsideTimeDomain
            })
        ));
        gate.request();
        assert!(matches!(
            scene.sample_sequential_ray_at_time(cx, &optical, 0, clock, settings(), 31, 42),
            Err(TransientCameraError::Camera(CameraRadianceError::Scene(
                CurvedSceneError::Scene(SceneError::Cancelled)
            )))
        ));
    });
}

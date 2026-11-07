use super::*;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};

fn counts<P>(
    sample: &CurvedSequentialRayleighPixelSample<P>,
    rays: &[Option<Ray>; CAMERA_WAVELENGTHS],
    coefficient: f64,
) {
    assert_eq!(sample.volume_events(), &[0; CAMERA_WAVELENGTHS]);
    assert_eq!(
        sample.free_flight_draws(),
        &rays.map(|ray| u32::from(ray.is_some() && coefficient > 0.0))
    );
}

#[test]
fn molecular_pixels_preserve_zero_scattering_bits_and_clipped_lane_counters() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = Field {
            curved: false,
            calls: AtomicUsize::new(0),
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([71; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [emitter(380.0, 780.0)];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let wide = camera(visible_glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &surface, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let radii = stop_radii(&optical_packet(&wide, cx, full.wavelengths_nm()));
        let physical = camera(visible_glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let analyzer = CameraAnalyzer::try_circular(Vec3::new(1.0, 0.0, 0.0), true).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(380.0, [0.0, 1.0]), (780.0, [1.0, 0.0])])
                .unwrap();
        let visible = raster
            .sample_pixel_sequential(cx, &surface, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let analyzed = raster
            .sample_pixel_sequential_analyzed(
                cx,
                &surface,
                0,
                analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let qe = raster
            .sample_pixel_sensor_sequential(cx, &surface, &sensor, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let qe_analyzed = raster
            .sample_pixel_sensor_sequential_analyzed(
                cx,
                &surface,
                &sensor,
                0,
                analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(visible.camera_rays().iter().flatten().count(), 2);
        assert_eq!(qe.camera_rays().iter().flatten().count(), 2);
        let trace_settings = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 8,
        };
        let phase = PhaseFunction::Rayleigh;
        for coefficient in [0.0, 1e-12] {
            let spectrum = MassScatteringSpectrum::constant(coefficient).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
            // At 1e-12/m every representable Philox flight exceeds this
            // scene. The surviving conditional signal is identical, without
            // an extra scattering extinction factor on its Beer response.
            let actual = raster
                .sample_pixel_sequential_rayleigh(cx, &scene, 0, trace_settings, SEED, 0, SAMPLE)
                .unwrap();
            assert_eq!(actual.pixel_sample(), &visible);
            counts(&actual, visible.camera_rays(), coefficient);
            assert_eq!(
                actual,
                raster
                    .sample_pixel_sequential_rayleigh(
                        cx,
                        &scene,
                        0,
                        trace_settings,
                        SEED,
                        0,
                        SAMPLE
                    )
                    .unwrap()
            );
            let actual = raster
                .sample_pixel_sequential_rayleigh_analyzed(
                    cx,
                    &scene,
                    0,
                    analyzer,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(actual.pixel_sample(), &analyzed);
            counts(&actual, analyzed.camera_rays(), coefficient);
            let actual = raster
                .sample_pixel_sensor_sequential_rayleigh(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(actual.pixel_sample(), &qe);
            counts(&actual, qe.camera_rays(), coefficient);
            let actual = raster
                .sample_pixel_sensor_sequential_rayleigh_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    analyzer,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(actual.pixel_sample(), &qe_analyzed);
            counts(&actual, qe_analyzed.camera_rays(), coefficient);
        }
    });
}

#[test]
fn real_volume_limit_aborts_all_sequential_observers() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = Field {
            curved: false,
            calls: AtomicUsize::new(0),
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([72; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [emitter(380.0, 780.0)];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(100.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let physical = camera(visible_glass(), 550.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(380.0, [1.0]), (780.0, [1.0])]).unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let trace_settings = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 1,
        };
        let expected = CurvedRasterError::Transport(CameraRadianceError::Unresolved(
            PathTermination::VolumeLimit,
        ));
        assert_eq!(
            raster
                .sample_pixel_sequential_rayleigh(cx, &scene, 0, trace_settings, SEED, 0, SAMPLE)
                .unwrap_err(),
            expected
        );
        assert_eq!(
            raster
                .sample_pixel_sequential_rayleigh_analyzed(
                    cx,
                    &scene,
                    0,
                    analyzer,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap_err(),
            expected
        );
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_rayleigh(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap_err(),
            expected
        );
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_rayleigh_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    analyzer,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap_err(),
            expected
        );
    });
}

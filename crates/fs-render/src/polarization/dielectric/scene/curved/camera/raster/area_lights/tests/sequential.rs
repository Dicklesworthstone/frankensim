//! Steady coated-lens pixels preserve original detector and clipping semantics.

use super::*;
use crate::camera::{Aperture, CameraProjection, LensSurface, PhysicalCamera, SequentialLens};
use crate::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::DielectricCoating;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, SpectralCameraAnalyzer,
};
use crate::spectral::detector::QuantumEfficiencySensor;

fn coated_camera(stop_radius: f64) -> PhysicalCamera {
    let film = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.02, 110e-9)
                    .unwrap(),
            ],
            1,
            ConductorSource::try_new(
                hash_domain(
                    "org.frankensim.test.steady-area-pixel",
                    b"authored absorbing coating",
                ),
                ConductorDataStatus::Representative,
            )
            .unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(glass))
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.09, -20.0, 0.012, None)
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.12, 0.0, stop_radius, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

#[test]
fn steady_coated_area_pixels_keep_original_cie_qe_and_sensor_analyzers() {
    fixtures::with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![panel(9, 1.0, 10.0)]).unwrap();
        let vacuum = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let physical = coated_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let analyzer = CameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            Stokes {
                i: 0.55,
                q: 0.2,
                u: -0.1,
                v: 0.05,
            },
        )
        .unwrap();
        let spectral = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (
                    LAMBDA_MIN,
                    Stokes {
                        i: 0.45,
                        q: 0.1,
                        u: 0.0,
                        v: 0.05,
                    },
                ),
                (
                    LAMBDA_MAX,
                    Stokes {
                        i: 0.65,
                        q: -0.15,
                        u: 0.08,
                        v: -0.03,
                    },
                ),
            ],
        )
        .unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75]),
            (LAMBDA_MAX, [0.25, 0.75]),
        ])
        .unwrap();
        let (seed, pixel, sample) = (79, 0, 6);

        // With no area proposals, the same original optical seed must produce
        // the exact native emission walk, including nonzero Q/U/V and coating.
        let expected = raster
            .sample_pixel_sequential(cx, &camera, 0, settings().surfaces, seed, pixel, sample)
            .unwrap();
        assert!(expected.spectral_samples().iter().any(Option::is_some));
        let actual = raster
            .sample_pixel_sequential_area_lights(
                cx,
                &area,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);
        assert_eq!(actual.area_samples(), &[0; 4]);
        assert_eq!(actual.trace_evaluations(), &[0; 4]);

        let expected = raster
            .sample_pixel_sequential_analyzed(
                cx,
                &camera,
                0,
                analyzer,
                settings().surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sequential_area_lights_analyzed(
                cx,
                &area,
                0,
                analyzer,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);
        let expected = raster
            .sample_pixel_sequential_spectrally_analyzed(
                cx,
                &camera,
                0,
                &spectral,
                settings().surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sequential_area_lights_spectrally_analyzed(
                cx,
                &area,
                0,
                &spectral,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);

        let expected = raster
            .sample_pixel_sensor_sequential(
                cx,
                &camera,
                &sensor,
                0,
                settings().surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sensor_sequential_area_lights(
                cx,
                &area,
                &sensor,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);
        let expected = raster
            .sample_pixel_sensor_sequential_analyzed(
                cx,
                &camera,
                &sensor,
                0,
                analyzer,
                settings().surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sensor_sequential_area_lights_analyzed(
                cx,
                &area,
                &sensor,
                0,
                analyzer,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);
        let expected = raster
            .sample_pixel_sensor_sequential_spectrally_analyzed(
                cx,
                &camera,
                &sensor,
                0,
                &spectral,
                settings().surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sensor_sequential_area_lights_spectrally_analyzed(
                cx,
                &area,
                &sensor,
                0,
                &spectral,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), &expected);

        let closed_physical = coated_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_physical, 1, 1).unwrap();
        let dark = closed
            .sample_pixel_sequential_area_lights(
                cx,
                &area,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(dark.pixel_sample().xyz(), [0.0; 3]);
        assert_eq!(dark.area_samples(), &[0; 4]);
        let dark = closed
            .sample_pixel_sensor_sequential_area_lights(
                cx,
                &area,
                &sensor,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 2]
        );
        assert_eq!(dark.trace_evaluations(), &[0; 4]);
        // Analyzer support is admitted before a physically clipped packet.
        let incomplete = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (500.0, Stokes::unpolarized(0.5)),
                (600.0, Stokes::unpolarized(0.5)),
            ],
        )
        .unwrap();
        assert!(
            closed
                .sample_pixel_sequential_area_lights_spectrally_analyzed(
                    cx,
                    &area,
                    0,
                    &incomplete,
                    settings(),
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .is_err()
        );
        assert!(
            closed
                .sample_pixel_sensor_sequential_area_lights_spectrally_analyzed(
                    cx,
                    &area,
                    &sensor,
                    0,
                    &incomplete,
                    settings(),
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .is_err()
        );
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_area_lights(
                    cx,
                    &area,
                    &sensor,
                    0,
                    settings(),
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .is_err()
        );
    });
}

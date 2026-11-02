use super::super::tests::{bindings, camera, geometry, settings, with_cx};
use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};

#[test]
fn zero_scattering_is_bit_identical_and_thin_gas_retains_conditional_beer_weight() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.2).unwrap(),
            BeerLambertAbsorption::try_constant(0.1).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(2.0).unwrap(),
            EmissionSides::Both,
        )];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let physical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let expected = raster
            .sample_pixel(cx, &surface, 0, settings(), 73, 4, 9)
            .unwrap();
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
            let result = raster
                .sample_pixel_rayleigh(cx, &scene, 0, trace_settings, 73, 4, 9)
                .unwrap();
            assert_eq!(result.pixel_sample(), &expected);
            assert_eq!(result.volume_events(), &[0; CAMERA_WAVELENGTHS]);
            assert_eq!(
                result.free_flight_draws(),
                &[u32::from(coefficient > 0.0); CAMERA_WAVELENGTHS]
            );
            // At 1e-12/m every possible finite Philox flight exceeds this
            // 2 m scene. Its conditional weight is Beer absorption alone;
            // scattering survival belongs to flight probability, not a second
            // extinction multiplier. Different scheduling keeps every bit.
            let _ = raster
                .sample_pixel_rayleigh(cx, &scene, 0, trace_settings, 73, 2, 5)
                .unwrap();
            assert_eq!(
                result,
                raster
                    .sample_pixel_rayleigh(cx, &scene, 0, trace_settings, 73, 4, 9)
                    .unwrap()
            );
        }
    });
}

#[test]
fn volume_limit_propagates_through_the_whole_spectral_pixel() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(2.0).unwrap(),
            EmissionSides::Both,
        )];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(100.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let trace_settings = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 1,
        };
        assert!(matches!(
            raster.sample_pixel_rayleigh(cx, &scene, 0, trace_settings, 73, 0, 0),
            Err(CurvedRasterError::Transport(
                CameraRadianceError::Unresolved(PathTermination::VolumeLimit)
            ))
        ));
    });
}

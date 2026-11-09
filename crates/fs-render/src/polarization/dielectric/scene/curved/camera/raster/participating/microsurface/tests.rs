use super::super::super::tests::{camera, geometry, settings, with_cx};
use super::*;
use crate::camera::{Aperture, CameraProjection, LensSurface, PhysicalCamera, SequentialLens};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::lighting::EnvironmentMap;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithDielectricCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EnvironmentRadiance,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_geom::{Point3, Vec3};

fn glass(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn lens(stop: Option<f64>) -> PhysicalCamera {
    let mut surfaces = vec![
        LensSurface::try_new(0.08, 20.0, 0.01, Some(glass(1.5))).unwrap(),
        LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
    ];
    if let Some(radius) = stop {
        surfaces.push(LensSurface::try_new(0.12, 0.0, radius, None).unwrap());
    }
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
    .unwrap()
}

#[test]
fn full_glass_volume_pixels_preserve_native_zero_scattering_and_detector_measure() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let air = glass(1.0);
        let material = glass(1.5);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let bindings = vec![SurfaceBinding {
            object_id: 7,
            interaction: Interaction::Dielectric {
                outside: 1,
                inside: 0,
                surface: DielectricSurface::try_anisotropic_rough(
                    0.9,
                    0.6,
                    Vec3::new(1.0, 0.0, 0.0),
                )
                .unwrap(),
            },
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings, 100_000).unwrap();
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let surface = CurvedCameraScene::try_new(
            cx,
            &native,
            &[],
            Some(EnvironmentRadiance::try_new(&map, 3.0).unwrap()),
        )
        .unwrap();
        let scattering = CurvedRayleighScene::try_new(cx, &native, &[None, None]).unwrap();
        let volume = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let limits = SmithWalkSettings::try_new(512).unwrap();
        let smith = SmithDielectricRayleighCameraScene::new(&volume, limits);
        let plain = SmithDielectricCameraScene::new(&surface, limits);
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 2, 2).unwrap();
        let trace = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 8,
        };
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        for sample in 0..8 {
            let result = raster
                .sample_pixel_smith_dielectric_rayleigh(
                    cx,
                    &smith,
                    0,
                    Some(analyzer),
                    trace,
                    73,
                    1,
                    sample,
                )
                .unwrap();
            let expected = raster
                .sample_pixel_smith_dielectric(
                    cx,
                    &plain,
                    0,
                    Some(analyzer),
                    settings(),
                    73,
                    1,
                    sample,
                )
                .unwrap();
            assert_eq!(result.pixel_sample(), &expected);
            assert_eq!(result.volume_events(), &[0; CAMERA_WAVELENGTHS]);
            assert_eq!(result.free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        }
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(380.0, [0.0, 1.0]), (780.0, [1.0, 0.0])])
                .unwrap();
        let result = raster
            .sample_pixel_sensor_smith_dielectric_rayleigh(
                cx,
                &smith,
                &sensor,
                0,
                Some(analyzer),
                trace,
                73,
                1,
                4,
            )
            .unwrap();
        let pixel = result.pixel_sample();
        let mut energy = [0.0; 2];
        let mut photons = [0.0; 2];
        for (draw, radiance) in pixel
            .detector_wavelengths()
            .iter()
            .zip(pixel.spectral_samples())
        {
            let w = draw.wavelength_nm();
            for (channel, qe) in [(w - 380.0) / 400.0, (780.0 - w) / 400.0]
                .into_iter()
                .enumerate()
            {
                // Complementary channels have a uniform 1/400 nm proposal.
                let power = radiance.radiance_w_m2_sr_nm() * qe * 100.0;
                energy[channel] += power;
                photons[channel] += power * w * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for c in 0..2 {
            assert!((pixel.qe_weighted_radiance_w_m2_sr()[c] - energy[c]).abs() < 1e-10);
            assert!(
                (pixel.photoelectron_rate_per_m2_sr()[c] - photons[c]).abs() < 1e-12 * photons[c]
            );
        }
        assert_eq!(
            result,
            raster
                .sample_pixel_sensor_smith_dielectric_rayleigh(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    Some(analyzer),
                    trace,
                    73,
                    1,
                    4,
                )
                .unwrap()
        );

        for stop in [None, Some(1e-8)] {
            let physical = lens(stop);
            let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
            let actual = raster
                .sample_pixel_sequential_smith_dielectric_rayleigh(
                    cx,
                    &smith,
                    0,
                    Some(analyzer),
                    trace,
                    79,
                    0,
                    6,
                )
                .unwrap();
            let expected = raster
                .sample_pixel_sequential_smith_dielectric(
                    cx,
                    &plain,
                    0,
                    Some(analyzer),
                    settings(),
                    79,
                    0,
                    6,
                )
                .unwrap();
            assert_eq!(actual.pixel_sample(), &expected);
            assert_eq!(actual.volume_events(), &[0; CAMERA_WAVELENGTHS]);
            assert_eq!(actual.free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
            let qe = raster
                .sample_pixel_sensor_sequential_smith_dielectric_rayleigh(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    Some(analyzer),
                    trace,
                    79,
                    0,
                    6,
                )
                .unwrap();
            let live = qe.pixel_sample().camera_rays().iter().flatten().count();
            if stop.is_some() {
                assert_eq!(live, 0);
                assert_eq!(qe.pixel_sample().qe_weighted_radiance_w_m2_sr(), &[0.0; 2]);
            } else {
                assert!(live > 0);
                assert!(
                    qe.pixel_sample()
                        .qe_weighted_radiance_w_m2_sr()
                        .iter()
                        .all(|x| *x > 0.0)
                );
            }
            assert_eq!(qe.free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        }
    });
}

#[test]
fn unresolved_volume_or_cancelled_pixel_never_publishes_a_partial_sample() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let air = glass(1.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 7,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let spectrum = MassScatteringSpectrum::constant(100.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let media = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &media).unwrap();
        let camera_scene = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let smith = SmithDielectricRayleighCameraScene::new(
            &camera_scene,
            SmithWalkSettings::try_new(512).unwrap(),
        );
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let trace = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 1,
        };
        assert!(matches!(
            raster.sample_pixel_smith_dielectric_rayleigh(cx, &smith, 0, None, trace, 73, 0, 0,),
            Err(CurvedRasterError::Transport(
                CameraRadianceError::Unresolved(PathTermination::VolumeLimit)
            ))
        ));
        gate.request();
        assert!(
            raster
                .sample_pixel_smith_dielectric_rayleigh(cx, &smith, 0, None, trace, 73, 0, 0,)
                .is_err()
        );
    });
}

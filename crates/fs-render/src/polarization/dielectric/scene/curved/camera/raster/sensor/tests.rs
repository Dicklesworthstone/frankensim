use super::super::tests::{bindings, camera, geometry, settings, with_cx};
use super::*;
use crate::dielectric::{
    BeerLambertAbsorption, DielectricGlass, DielectricIorSample, GlassProvenance, TabulatedIor,
};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

mod geometric_dispersion;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn infrared_glass() -> DielectricGlass {
    let ior = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.2).unwrap()
        }),
        750.0,
    )
    .unwrap();
    DielectricGlass::new(ior, BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}

fn infrared_emitter(value: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        7,
        RadianceSpectrum::try_new_in_band(vec![(8000.0, value), (14000.0, value)]).unwrap(),
        EmissionSides::Both,
    )
}

#[test]
fn g1_native_infrared_pixels_integrate_detector_energy_and_wavelength_correct_electrons() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = infrared_glass();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 1, 4)
            .unwrap();
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for (draw, path) in sample
            .detector_wavelengths()
            .iter()
            .zip(sample.spectral_samples())
        {
            let wavelength = draw.wavelength_nm();
            assert!((8000.0..=14000.0).contains(&wavelength));
            let eta = [
                (wavelength - 8000.0) / 6000.0,
                (14000.0 - wavelength) / 6000.0,
                0.0,
            ];
            near(path.radiance_w_m2_sr_nm(), 3.0, 2e-12);
            for c in 0..3 {
                // The complementary source curves have a uniform sum, so
                // p(lambda)=1/6000nm independently of the production sampler.
                let weighted_power = 3.0 * eta[c] * 6000.0 / 4.0;
                energy[c] += weighted_power;
                electrons[c] +=
                    weighted_power * (wavelength * 1e-9) / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for c in 0..3 {
            near(sample.qe_weighted_radiance_w_m2_sr()[c], energy[c], 2e-8);
            near(
                sample.photoelectron_rate_per_m2_sr()[c],
                electrons[c],
                2e-12 * electrons[c].max(1.0),
            );
        }
        near(energy[0] + energy[1], 18_000.0, 2e-10);
        assert_eq!(sample.photoelectron_rate_per_m2_sr()[2], 0.0);
        raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 0, 19)
            .unwrap();
        assert_eq!(
            sample,
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 1, 4)
                .unwrap()
        );
    });
}

struct InfraredVerticalIndex;
impl GradientIndexField for InfraredVerticalIndex {
    fn sample(
        &self,
        p: Point3,
        wavelength: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if !(8000.0..=14000.0).contains(&wavelength) {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        let c = (11000.0 / wavelength).powi(2);
        let n = 1.28 + 0.03 * p.z + (0.02 + 0.01 * p.z) * c;
        GradientIndexSample::try_new(
            n,
            Vec3::new(0.0, 0.0, 0.03 + 0.01 * c),
            0.12,
            Some(n + 2.0 * (0.02 + 0.01 * p.z) * c),
        )
    }
}

#[test]
fn g1_detector_wavelengths_keep_independent_curved_geometry_and_absorption() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = InfraredVerticalIndex;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([22; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [0.8]), (14000.0, [0.8])]).unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 79, 0, 6)
            .unwrap();
        let initial = sample.camera_ray();
        let transverse = initial.dir.x.powi(2) + initial.dir.y.powi(2);
        let mut expected_energy = 0.0;
        for (draw, path) in sample
            .detector_wavelengths()
            .iter()
            .zip(sample.spectral_samples())
        {
            let c = (11000.0 / draw.wavelength_nm()).powi(2);
            let n0 = 1.28 + 0.02 * c;
            let n1 = 1.34 + 0.04 * c;
            let p2 = n0 * n0 * transverse;
            let length = 2.0 * (n1 + n0) / ((n1 * n1 - p2).sqrt() + (n0 * n0 - p2).sqrt());
            let radiance = 3.0 * (n0 / n1).powi(2) * (-0.12 * length).exp();
            near(path.radiance_w_m2_sr_nm(), radiance, 2e-7);
            near(path.path().traveled_m(), length, 2e-7);
            near(
                path.path().beam().ray().dir.x,
                initial.dir.x * n0 / n1,
                2e-8,
            );
            expected_energy += radiance * 6000.0 * 0.8 / 4.0;
        }
        near(
            sample.qe_weighted_radiance_w_m2_sr()[0],
            expected_energy,
            1e-3,
        );
        assert_ne!(
            sample.spectral_samples()[0].path().beam().ray().dir,
            sample.spectral_samples()[1].path().beam().ray().dir
        );
    });
}

#[test]
fn g0_missing_source_support_and_cancelled_transport_never_publish_black_channels() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let glass = infrared_glass();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        let outside =
            QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (600.0, [1.0])]).unwrap();
        assert!(
            raster
                .sample_pixel_sensor(cx, &scene, &outside, 0, settings(), 2, 0, 0)
                .is_err()
        );
        let supported =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [1.0]), (14000.0, [1.0])]).unwrap();
        assert!(
            raster
                .sample_pixel_sensor(cx, &scene, &supported, 0, settings(), 2, 0, 0)
                .is_ok()
        );
        gate.request();
        assert!(
            raster
                .sample_pixel_sensor(cx, &scene, &supported, 0, settings(), 2, 0, 0)
                .is_err()
        );
    });
}

#[test]
fn g3_infrared_detector_preserves_native_scattering_weights_and_volume_refusals() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = infrared_glass();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [0.2, 0.9]), (14000.0, [0.8, 0.1])])
                .unwrap();
        let expected = raster
            .sample_pixel_sensor(cx, &surface, &sensor, 0, settings(), 73, 1, 4)
            .unwrap();
        let phase = PhaseFunction::Rayleigh;
        let trace_settings = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 1,
        };
        for coefficient in [0.0, 1e-12, 100.0] {
            let spectrum = MassScatteringSpectrum::source_table(
                &[(8000.0, coefficient), (14000.0, coefficient)],
                ContentHash([41; 32]),
            )
            .unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let scene = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
            let result = raster.sample_pixel_sensor_rayleigh(
                cx,
                &scene,
                &sensor,
                0,
                trace_settings,
                73,
                1,
                4,
            );
            if coefficient == 100.0 {
                assert!(matches!(
                    result,
                    Err(CurvedRasterError::Transport(
                        CameraRadianceError::Unresolved(PathTermination::VolumeLimit)
                    ))
                ));
            } else {
                let result = result.unwrap();
                assert_eq!(result.pixel_sample(), &expected);
                assert_eq!(result.volume_events(), &[0; CAMERA_WAVELENGTHS]);
                assert_eq!(
                    result.free_flight_draws(),
                    &[u32::from(coefficient > 0.0); CAMERA_WAVELENGTHS]
                );
            }
        }
    });
}

//! G1/G3: original infrared Planck/Kirchhoff paths and physical detector weights.

use super::super::raster::CurvedCameraRaster;
use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::thermal_surface::{SpectralDiffuseThermalSurface, SpectralEmissivity};

fn infrared_medium(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(8_000.0 + 750.0 * i as f64, index).unwrap()
            }),
            750.0,
        )
        .unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn source(id: u64, epsilon: f64, temperature: f64) -> DiffuseThermalEmitter {
    DiffuseThermalEmitter::new(
        id,
        SpectralDiffuseThermalSurface::try_new(
            temperature,
            SpectralEmissivity::constant_in_band(epsilon, (8_000.0, 14_000.0)).unwrap(),
            ContentHash([id as u8; 32]),
        )
        .unwrap(),
    )
}

fn blackbody(temperature: f64, wavelength_nm: f64) -> f64 {
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    let lambda = wavelength_nm * 1e-9;
    2.0 * h * c * c * 1e-9 / (lambda.powi(5) * (h * c / (lambda * k * temperature)).exp_m1())
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-30),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_native_reflection_and_emission_share_one_infrared_kirchhoff_law() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0), plane(9, -1.0, 0.0, 10.0)])
                .unwrap();
        let vacuum = infrared_medium(1.0);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let front = source(1, 0.37, 300.0);
        let back = source(9, 1.0, 300.0);
        let reflectances = [front.reflectance_spectrum().unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[front, back])
            .unwrap();
        for wavelength in [8_000.0, 9_250.0, 11_000.0, 14_000.0] {
            let sample = camera
                .sample_ray_with_samples(cx, ray(0.0, 1.0), wavelength, 0, settings(), |_| {
                    [0.0, 0.0, 0.0]
                })
                .unwrap();
            // Direct epsilon*B plus the native reflected (1-epsilon)*B is B.
            relative(sample.radiance_w_m2_sr_nm(), blackbody(300.0, wavelength));
            assert_eq!(sample.path().last_contact().unwrap().object_id, 9);
        }
        assert_eq!(front.surface().temperature_k(), 300.0);
        assert_eq!(back.surface().temperature_k(), 300.0);
    });
}

#[test]
fn g1_infrared_sensor_observes_planck_at_its_actual_wavelengths_without_tabulation() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let vacuum = infrared_medium(1.0);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitter = source(9, 1.0, 300.0);
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[emitter])
            .unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.1).unwrap(),
            1.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8_000.0, [0.0, 1.0]), (14_000.0, [1.0, 0.0])])
                .unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 0, 4)
            .unwrap();
        let mut power = [0.0; 2];
        let mut electrons = [0.0; 2];
        for (draw, lane) in sample
            .detector_wavelengths()
            .iter()
            .zip(sample.spectral_samples())
        {
            let wavelength = draw.wavelength_nm();
            let radiance = blackbody(300.0, wavelength);
            relative(lane.radiance_w_m2_sr_nm(), radiance);
            let fraction = (wavelength - 8_000.0) / 6_000.0;
            for (channel, qe) in [fraction, 1.0 - fraction].into_iter().enumerate() {
                let weighted = radiance * qe * 6_000.0 / 4.0;
                power[channel] += weighted;
                electrons[channel] +=
                    weighted * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for c in 0..2 {
            relative(sample.qe_weighted_radiance_w_m2_sr()[c], power[c]);
            relative(sample.photoelectron_rate_per_m2_sr()[c], electrons[c]);
        }
        assert_eq!(
            sample,
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 0, 4)
                .unwrap()
        );
    });
}

#[test]
fn g0_thermal_bindings_refuse_other_reflectance_double_sources_and_nonvacuum_emission() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let vacuum = infrared_medium(1.0);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let black = source(9, 1.0, 300.0);
        let camera = || CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        assert!(
            camera()
                .with_diffuse_thermal_emitters(cx, &[source(9, 0.5, 300.0)])
                .is_err()
        );
        assert!(
            camera()
                .with_diffuse_thermal_emitters(cx, &[source(8, 1.0, 300.0)])
                .is_err()
        );
        assert!(
            camera()
                .with_diffuse_thermal_emitters(cx, &[black, black])
                .is_err()
        );
        let calibrated = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new_in_band(vec![(8_000.0, 1.0), (14_000.0, 1.0)]).unwrap(),
            EmissionSides::Both,
        )];
        assert!(
            CurvedCameraScene::try_new(cx, &native, &calibrated, None)
                .unwrap()
                .with_diffuse_thermal_emitters(cx, &[black])
                .is_err()
        );
        let disabled = camera()
            .with_diffuse_thermal_emitters(cx, &[black])
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[])
            .unwrap();
        assert_eq!(
            disabled.sample_ray(cx, ray(0.0, 1.0), 10_000.0, 0, settings(), 31, 2),
            camera().sample_ray(cx, ray(0.0, 1.0), 10_000.0, 0, settings(), 31, 2)
        );
        let host = infrared_medium(1.5);
        let refractive_media = [SceneMedium::Homogeneous(&host)];
        let refractive_native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &refractive_media,
            vec![opaque(9)],
            100_000,
        )
        .unwrap();
        let refractive = CurvedCameraScene::try_new(cx, &refractive_native, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[black])
            .unwrap();
        assert!(matches!(
            refractive.sample_ray(cx, ray(0.0, 1.0), 10_000.0, 0, settings(), 31, 2),
            Err(CameraRadianceError::UnsupportedSurface { object_id: 9 })
        ));
        let pending = camera();
        gate.request();
        assert!(pending.with_diffuse_thermal_emitters(cx, &[black]).is_err());
    });
}

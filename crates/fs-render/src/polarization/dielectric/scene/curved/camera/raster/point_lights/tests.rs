use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::instances::RigidTransform;
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::CurvedCameraScene;
use crate::polarization::dielectric::scene::curved::camera::point_lights::PointLightBranch;
use crate::polarization::dielectric::scene::curved::camera::raster::tests::{
    camera, geometry, settings, with_cx,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::spectral::{SpectralReflectance, cie_x, cie_y, cie_z, y_integral};
use fs_geom::Vec3;

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.0).unwrap(),
        GlassProvenance::Custom,
    )
}

fn point_settings() -> PointLightCameraSettings {
    let mut connection = BundleConnectionOptions::default();
    connection.scene.trace = settings();
    PointLightCameraSettings {
        surfaces: settings(),
        connection,
        max_connections: 4,
        max_trace_evaluations: 96,
    }
}

fn bindings() -> Vec<SurfaceBinding> {
    vec![SurfaceBinding {
        object_id: 7,
        interaction: Interaction::DepolarizingDiffuse {
            medium: 0,
            reflectance: 0,
        },
    }]
}

#[test]
fn g3_empty_fixture_binding_preserves_physical_camera_spectral_pixel_bits() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.37).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(),
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let connected = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
        let camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&camera, 3, 2).unwrap();
        let original = raster
            .sample_pixel(cx, &surface, 0, settings(), 73, 4, 9)
            .unwrap();
        let result = raster
            .sample_pixel_point_lights(cx, &connected, 0, point_settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(*result.pixel_sample(), original);
        assert_eq!(
            result.fixture_radiance_w_m2_sr_nm(),
            &[0.0; CAMERA_WAVELENGTHS]
        );
        assert_eq!(result.connections(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.trace_evaluations(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.field_evaluations(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.query_work(), &[0; CAMERA_WAVELENGTHS]);
    });
}

#[test]
fn g1_each_pixel_lane_receives_calibrated_diffuse_point_light_with_atomic_refusals() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let rho = 0.37;
        let reflectances = [SpectralReflectance::constant(rho).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(),
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::identity(),
            1000.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins(core::array::from_fn(|i| 1.0 + i as f64 / 80.0)).unwrap(),
        )
        .unwrap();
        let lights = [PointLightBranch {
            source: &source,
            receiver_object_id: 7,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
            position_tolerance_m: 1.0e-6,
        }];
        let connected = CurvedPointLightCameraScene::try_new(cx, &surface, &lights).unwrap();
        let camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&camera, 3, 2).unwrap();
        let result = raster
            .sample_pixel_point_lights(cx, &connected, 0, point_settings(), 73, 4, 9)
            .unwrap();
        let repeat = raster
            .sample_pixel_point_lights(cx, &connected, 0, point_settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(result, repeat);
        let original = raster
            .sample_pixel(cx, &surface, 0, settings(), 73, 4, 9)
            .unwrap();
        let pixel = result.pixel_sample();
        assert_eq!(pixel.camera_ray(), original.camera_ray());
        assert_eq!(pixel.wavelengths_nm(), original.wavelengths_nm());
        let direction = pixel.camera_ray().dir;
        // Camera and fixture are at the origin; the actual receiving plane is
        // z=2. Its irradiance is I*cos(theta)/r². No light-path proposal is used.
        let distance = 2.0 / direction.z;
        let mut expected_xyz = [0.0; 3];
        for (lane, &wavelength) in pixel.wavelengths_nm().iter().enumerate() {
            let intensity = source
                .spectral_intensity_w_sr_nm(direction, wavelength)
                .unwrap();
            let expected =
                rho * intensity * direction.z / (core::f64::consts::PI * distance * distance);
            let actual = pixel.spectral_samples()[lane].radiance_w_m2_sr_nm();
            assert!((actual - expected).abs() <= 3.0e-6 * expected + 1.0e-14);
            assert_eq!(actual, result.fixture_radiance_w_m2_sr_nm()[lane]);
            assert_eq!(
                pixel.spectral_samples()[lane].path(),
                original.spectral_samples()[lane].path()
            );
            assert_eq!(result.connections()[lane], 1);
            assert!(result.trace_evaluations()[lane] >= 1);
            assert!(result.field_evaluations()[lane] > 0);
            assert!(result.query_work()[lane] > 0);
            for (sum, observer) in expected_xyz.iter_mut().zip([
                cie_x(wavelength),
                cie_y(wavelength),
                cie_z(wavelength),
            ]) {
                *sum += expected * observer * 100.0 / y_integral();
            }
        }
        for (actual, expected) in pixel.xyz().into_iter().zip(expected_xyz) {
            assert!((actual - expected).abs() <= 3.0e-6 * expected + 1.0e-14);
        }
        let mut exhausted = point_settings();
        exhausted.max_trace_evaluations = 0;
        assert!(matches!(
            raster.sample_pixel_point_lights(cx, &connected, 0, exhausted, 73, 4, 9),
            Err(PointLightRasterError::Transport { lane: 0, .. })
        ));
        assert!(matches!(
            raster.sample_pixel_point_lights(cx, &connected, 0, point_settings(), 73, 6, 9),
            Err(PointLightRasterError::Raster(
                CurvedRasterError::InvalidRaster
            ))
        ));
    });
}

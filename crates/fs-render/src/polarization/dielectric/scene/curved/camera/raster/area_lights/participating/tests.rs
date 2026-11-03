//! Native spectral streams and independently evaluated molecular source loss.

use super::*;
use super::super::super::tests as fixtures;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::BoundarySide;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleBranchPolicy;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::camera::area_lights::CurvedAreaEmitter;
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::{
    AreaLightBinding, AreaLightCameraSettings, AreaLightCameraWorkspace,
    CurvedAreaLightCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, SpectralReflectance};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::hash_domain;
use fs_geom::{Point3, Vec3};

fn material() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.3).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn settings() -> RayleighAreaLightCameraSettings {
    let mut area_lights = AreaLightCameraSettings::default();
    area_lights.surfaces = fixtures::settings();
    area_lights.surfaces.max_distance_m = 12.0;
    area_lights.connection.scene.trace = area_lights.surfaces;
    area_lights.connection.scene.position_derivative_tolerance_m = 1e-10;
    area_lights.connection.scene.direction_derivative_tolerance = 1e-11;
    area_lights.max_area_samples = 8;
    area_lights.max_trace_evaluations = 4096;
    RayleighAreaLightCameraSettings {
        area_lights,
        max_volume_events: 8,
        max_source_free_flight_draws: 16,
        max_source_candidates: 16,
    }
}

fn panel(id: u64, z: f64, scale: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("molecular-area-raster-panel", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-scale, -scale, z], [scale, -scale, z], [0.0, 2.0 * scale, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    ).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 5e-6 * expected.abs().max(1e-14),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn g0_empty_area_bindings_preserve_native_molecular_pixels_and_aperture() {
    fixtures::with_cx(|cx, _| {
        let geometry = fixtures::geometry();
        let glass = material();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native = CurvedPolarizedScene::try_new(
            cx, &geometry, &media, fixtures::bindings(), 100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(7,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let phase = PhaseFunction::Rayleigh;
        for coefficient in [0.0, 1e-12] {
            let spectrum = MassScatteringSpectrum::constant(coefficient).unwrap();
            let volumes = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let scene = CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &area).unwrap();
            let ordinary = raster.sample_pixel_rayleigh(cx, &molecular, 0,
                CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces,
                    max_volume_events: settings().max_volume_events }, 73, 4, 9).unwrap();
            let result = raster.sample_pixel_rayleigh_area_lights(
                cx, &scene, 0, settings(), &mut workspace, 73, 4, 9,
            ).unwrap();
            assert_eq!(result.pixel_sample(), ordinary.pixel_sample());
            assert_eq!(result.volume_events(), ordinary.volume_events());
            assert_eq!(result.free_flight_draws(), ordinary.free_flight_draws());
            assert_eq!(result.area_pixel_sample().area_radiance_w_m2_sr_nm(), &[0.0; 4]);
            assert_eq!(result.area_pixel_sample().area_samples(), &[0; 4]);
            assert_eq!(result.area_pixel_sample().trace_evaluations(), &[0; 4]);
            assert_eq!(result.source_free_flight_draws(), &[0; 4]);
            assert_eq!(result.source_candidates(), &[0; 4]);
            assert_eq!(result, raster.sample_pixel_rayleigh_area_lights(
                cx, &scene, 0, settings(), &mut workspace, 73, 4, 9,
            ).unwrap());
        }
    });
}

#[test]
fn g1_spectral_area_pixels_apply_source_survival_without_conditioned_camera_extinction() {
    fixtures::with_cx(|cx, _| {
        // The last black panel closes these sampled camera continuations;
        // positive scattering never relies on a vacuum escape assumption.
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0), panel(9, -1.0, 0.3), panel(10, -1.5, 50.0),
        ]).unwrap();
        let glass = material();
        let media = [SceneMedium::Homogeneous(&glass)];
        let rho = 0.6;
        let reflectances = [SpectralReflectance::constant(rho).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(cx, &geometry, &media,
            &[], &[], &reflectances, vec![
                SurfaceBinding { object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse { medium: 0, reflectance: 0 } },
                SurfaceBinding { object_id: 9, interaction: Interaction::Opaque },
                SurfaceBinding { object_id: 10, interaction: Interaction::Opaque },
            ], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed { initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir }];
        let bindings = [AreaLightBinding { source: &source, receiver_object_id: 1,
            source_side: BoundarySide::Positive, initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0), seeds: &seeds,
            image_options: BundleImageOptions { max_images: 1, ..BundleImageOptions::default() },
            position_tolerance_m: 2e-7 }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let mut surface_workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let key = [seed as u32, (seed >> 32) as u32];
        let phase = PhaseFunction::Rayleigh;
        for coefficient in [0.0, 0.01] {
            let spectrum = MassScatteringSpectrum::constant(coefficient).unwrap();
            let volumes = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let scene = CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &area).unwrap();
            let result = raster.sample_pixel_rayleigh_area_lights(
                cx, &scene, 0, settings(), &mut workspace, seed, pixel, sample,
            ).unwrap();
            let ordinary = raster.sample_pixel_rayleigh(cx, &molecular, 0,
                CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces,
                    max_volume_events: settings().max_volume_events }, seed, pixel, sample).unwrap();
            assert_eq!(result.volume_events(), &[0; CAMERA_WAVELENGTHS]);
            assert_eq!(result.free_flight_draws(), ordinary.free_flight_draws());
            assert_eq!(result.source_free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
            assert_eq!(result.source_candidates(), &[0; CAMERA_WAVELENGTHS]);
            let ray = result.pixel_sample().camera_ray();
            let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
            let receiver = Point3::new(ray.origin.x + camera_length * ray.dir.x,
                ray.origin.y + camera_length * ray.dir.y, 1.0);
            for lane in 0..CAMERA_WAVELENGTHS {
                assert_eq!(result.pixel_sample().spectral_samples()[lane].path(),
                    ordinary.pixel_sample().spectral_samples()[lane].path());
                let wavelength = result.pixel_sample().wavelengths_nm()[lane];
                let draw = philox4x32_10([pixel, sample, 0, 0x6372_6100 + lane as u32 * 64], key);
                let root = unit(draw[0]).sqrt();
                let v = unit(draw[1]);
                let point = Point3::new(-0.3 + 0.6 * root * (1.0 - v) + 0.3 * root * v,
                    -0.3 + 0.9 * root * v, -1.0);
                let displacement = receiver.delta_from(point);
                let r2 = displacement.dot(displacement);
                let distance = r2.sqrt();
                let cosine = 2.0 / distance;
                let le = 1.0 + 4.0 * (wavelength - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
                let irradiance = le * 0.27 * cosine * cosine / r2
                    * (-(0.2 + coefficient) * distance).exp();
                let light = r2 / (0.27 * cosine);
                let continuation = cosine / core::f64::consts::PI;
                let weight = light * light / (light * light + continuation * continuation);
                // The selected camera flights are uncollided. Their conditional
                // prefix contains Beer absorption, not a second survival factor.
                let expected = (-0.2 * camera_length).exp() * rho / core::f64::consts::PI
                    * irradiance * weight;
                near(result.area_pixel_sample().area_radiance_w_m2_sr_nm()[lane], expected);
                assert_eq!(result.area_pixel_sample().area_samples()[lane], 1);
                assert_eq!(result.area_pixel_sample().images()[lane], 1);
            }
            if coefficient == 0.0 {
                let surface = raster.sample_pixel_area_lights(cx, &area, 0,
                    settings().area_lights, &mut surface_workspace, seed, pixel, sample).unwrap();
                for lane in 0..CAMERA_WAVELENGTHS {
                    near(result.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                        surface.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm());
                }
            }
            let mut limited = settings();
            limited.area_lights.max_trace_evaluations = 1;
            assert!(matches!(raster.sample_pixel_rayleigh_area_lights(
                cx, &scene, 0, limited, &mut workspace, seed, pixel, sample,
            ), Err(AreaLightRasterError::Transport { lane: 0, .. })));
        }
    });
}

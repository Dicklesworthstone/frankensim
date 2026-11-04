//! Physical-camera pixels retain area-light radiometry in every spectral lane.

use super::*;
use super::super::{CAMERA_WAVELENGTHS, CurvedCameraRaster};
use super::super::tests as fixtures;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance, evaluate_rough_dielectric,
};
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
    AreaLightBinding, AreaLightCameraSettings, AreaLightCameraWorkspace, CurvedAreaLightCameraScene,
};
use crate::gradient_index::event::BoundarySide;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use crate::spectral::SpectralReflectance;
use fs_blake3::hash_domain;
use fs_geom::{Point3, Vec3};
use fs_rand::philox::philox4x32_10;

mod sensor;

fn material() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.3).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn settings() -> AreaLightCameraSettings {
    let mut settings = AreaLightCameraSettings::default();
    settings.surfaces = fixtures::settings();
    settings.connection.scene.trace = settings.surfaces;
    settings.connection.scene.position_derivative_tolerance_m = 1e-10;
    settings.connection.scene.direction_derivative_tolerance = 1e-11;
    settings.max_area_samples = 8;
    settings.max_trace_evaluations = 4096;
    settings
}

fn panel(id: u64, z: f64, scale: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("area-camera-raster-panel", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-scale, -scale, z], [scale, -scale, z], [0.0, 2.0 * scale, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn unit(word: u32) -> f64 {
    (f64::from(word) + 0.5) / 4_294_967_296.0
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-6 * expected.abs().max(1e-14),
        "{actual:.17e} != {expected:.17e}",
    );
}

#[test]
fn g0_area_pixels_with_empty_bindings_preserve_aperture_spectrum_and_sample_identity() {
    fixtures::with_cx(|cx, _| {
        let geometry = fixtures::geometry();
        let glass = material();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native = CurvedPolarizedScene::try_new(
            cx, &geometry, &media, fixtures::bindings(), 100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let ordinary = raster.sample_pixel(cx, &camera, 0, settings().surfaces, 73, 4, 9).unwrap();
        let result = raster.sample_pixel_area_lights(
            cx, &area, 0, settings(), &mut workspace, 73, 4, 9,
        ).unwrap();
        assert_eq!(result.pixel_sample(), &ordinary);
        assert_eq!(result.area_radiance_w_m2_sr_nm(), &[0.0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.area_samples(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.images(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.trace_evaluations(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.field_evaluations(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(result.query_work(), &[0; CAMERA_WAVELENGTHS]);
        let repeated = raster.sample_pixel_area_lights(
            cx, &area, 0, settings(), &mut workspace, 73, 4, 9,
        ).unwrap();
        assert_eq!(result.pixel_sample(), repeated.pixel_sample());
        let origin = result.pixel_sample().camera_ray().origin;
        assert_eq!(origin.z, 0.0);
        assert!(origin.x * origin.x + origin.y * origin.y <= 0.01_f64.powi(2));
    });
}

#[test]
fn g1_area_pixels_match_independent_spectral_irradiance_and_native_lane_streams() {
    fixtures::with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![panel(1, 1.0, 10.0), panel(9, -1.0, 0.3)]).unwrap();
        let glass = material();
        let media = [SceneMedium::Homogeneous(&glass)];
        let rho = 0.6;
        let reflectances = [SpectralReflectance::constant(rho).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![
                SurfaceBinding { object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse { medium: 0, reflectance: 0 } },
                SurfaceBinding { object_id: 9, interaction: Interaction::Opaque },
            ],
            100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions { max_images: 1, ..BundleImageOptions::default() },
            position_tolerance_m: 2e-7,
        }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let result = raster.sample_pixel_area_lights(
            cx, &area, 0, settings(), &mut workspace, seed, pixel, sample,
        ).unwrap();
        let repeated = raster.sample_pixel_area_lights(
            cx, &area, 0, settings(), &mut workspace, seed, pixel, sample,
        ).unwrap();
        assert_eq!(result.pixel_sample(), repeated.pixel_sample());
        assert_eq!(result.area_radiance_w_m2_sr_nm(), repeated.area_radiance_w_m2_sr_nm());
        assert_eq!(result.trace_evaluations(), repeated.trace_evaluations());
        assert_eq!(result.field_evaluations(), repeated.field_evaluations());
        assert_eq!(result.query_work(), repeated.query_work());
        let ray = result.pixel_sample().camera_ray();
        let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
        let receiver = Point3::new(
            ray.origin.x + camera_length * ray.dir.x,
            ray.origin.y + camera_length * ray.dir.y,
            1.0,
        );
        let key = [seed as u32, (seed >> 32) as u32];
        let render_cx = cx.with_stream_seed(seed);
        for lane in 0..CAMERA_WAVELENGTHS {
            let wavelength = result.pixel_sample().wavelengths_nm()[lane];
            let draw = philox4x32_10([pixel, sample, 0, 0x6372_6100 + lane as u32 * 64], key);
            let root = unit(draw[0]).sqrt();
            let v = unit(draw[1]);
            // Independent area coordinates on the authored triangle. No
            // prepared-source sampler, bundle Jacobian or MIS helper is used.
            let point = Point3::new(
                -0.3 + 0.6 * root * (1.0 - v) + 0.3 * root * v,
                -0.3 + 0.9 * root * v,
                -1.0,
            );
            let displacement = receiver.delta_from(point);
            let distance_squared = displacement.dot(displacement);
            let distance = distance_squared.sqrt();
            let cosine = 2.0 / distance;
            let radiance = 1.0 + 4.0 * (wavelength - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
            let irradiance = radiance * 0.27 * cosine * cosine / distance_squared
                * (-0.2 * distance).exp();
            let p_light = distance_squared / (0.27 * cosine);
            let p_bsdf = cosine / core::f64::consts::PI;
            let light_weight = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
            let expected = (-0.2 * camera_length).exp() * rho / core::f64::consts::PI
                * irradiance * light_weight;
            near(result.area_radiance_w_m2_sr_nm()[lane], expected);
            assert_eq!(result.area_samples()[lane], 1);
            assert_eq!(result.images()[lane], 1);

            // Replay each lane directly with the specified independent streams.
            // Exact equality also checks aperture/raster identity and full
            // continuation/emitter-hit transport beyond the analytic NEE term.
            let manual = area.sample_ray_with_samples(
                &render_cx, ray, wavelength, 0, settings(), &mut workspace,
                |event| {
                    let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                    core::array::from_fn(|i| unit(block[i]))
                },
                |event, binding| {
                    let block = philox4x32_10(
                        [pixel, sample, event, 0x6372_6100 + lane as u32 * 64 + binding as u32], key,
                    );
                    [unit(block[0]), unit(block[1])]
                },
            ).unwrap();
            assert_eq!(&result.pixel_sample().spectral_samples()[lane], manual.radiance_sample());
            assert_eq!(result.area_radiance_w_m2_sr_nm()[lane], manual.area_radiance_w_m2_sr_nm());
            assert_eq!(result.trace_evaluations()[lane], u64::from(manual.trace_evaluations()));
            assert_eq!(result.field_evaluations()[lane], manual.field_evaluations());
            assert_eq!(result.query_work()[lane], manual.query_work());
        }
        let mut limited = settings();
        limited.max_trace_evaluations = 1;
        assert!(matches!(
            raster.sample_pixel_area_lights(
                cx, &area, 0, limited, &mut workspace, seed, pixel, sample,
            ),
            Err(AreaLightRasterError::Transport { lane: 0, .. })
        ));
    });
}

#[test]
fn g1_rough_glass_pixels_match_dispersive_transmission_and_preserve_native_paths() {
    fixtures::with_cx(|cx, _| {
        // The original winding faces the light-side glass; the physical camera
        // approaches from the inside. The source emits back toward that interface.
        let emitter_geometry = GeometryInstance::try_new(
            9,
            hash_domain("rough-area-raster-source", b"9"),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-0.3, -0.3, 2.0], [0.0, 0.6, 2.0], [0.3, -0.3, 2.0]],
                vec![[0, 1, 2]],
            )),
            RigidTransform::identity(),
        ).unwrap();
        let geometry = InstanceScene::try_new(vec![panel(1, 1.0, 10.0), emitter_geometry]).unwrap();
        let camera_medium = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.1).unwrap(),
            GlassProvenance::Custom,
        );
        let source_medium = DielectricGlass::new(
            CauchyIor::try_new(1.45, 0.004, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.15).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [
            SceneMedium::Homogeneous(&camera_medium),
            SceneMedium::Homogeneous(&source_medium),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx, &geometry, &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::try_rough(0.25).unwrap(),
                    },
                },
                SurfaceBinding { object_id: 9, interaction: Interaction::Opaque },
            ],
            100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, -1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 1,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions { max_images: 1, ..BundleImageOptions::default() },
            position_tolerance_m: 2e-7,
        }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let ordinary = raster.sample_pixel(
            cx, &camera, 0, settings().surfaces, seed, pixel, sample,
        ).unwrap();
        let result = raster.sample_pixel_area_lights(
            cx, &area, 0, settings(), &mut workspace, seed, pixel, sample,
        ).unwrap();
        assert_eq!(result.pixel_sample().camera_ray(), ordinary.camera_ray());
        assert_eq!(result.pixel_sample().wavelengths_nm(), ordinary.wavelengths_nm());
        let ray = result.pixel_sample().camera_ray();
        let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
        let receiver = Point3::new(
            ray.origin.x + camera_length * ray.dir.x,
            ray.origin.y + camera_length * ray.dir.y,
            1.0,
        );
        let key = [seed as u32, (seed >> 32) as u32];
        for lane in 0..CAMERA_WAVELENGTHS {
            let wavelength = result.pixel_sample().wavelengths_nm()[lane];
            let draw = philox4x32_10([pixel, sample, 0, 0x6372_6100 + lane as u32 * 64], key);
            let root = unit(draw[0]).sqrt();
            let v = unit(draw[1]);
            let source_point = Point3::new(
                -0.3 + 0.3 * root * (1.0 - v) + 0.6 * root * v,
                -0.3 + 0.9 * root * (1.0 - v),
                2.0,
            );
            let displacement = source_point.delta_from(receiver);
            let distance_squared = displacement.dot(displacement);
            let distance = distance_squared.sqrt();
            let direction_to_source = displacement.scale(1.0 / distance);
            let cosine = 1.0 / distance;
            let wavelength_um = wavelength * 1e-3;
            let eta = 1.45 + 0.004 / (wavelength_um * wavelength_um);
            let bsdf = evaluate_rough_dielectric(
                Vec3::new(0.0, 0.0, -1.0), ray.dir.scale(-1.0), direction_to_source,
                1.0, eta, 0.25,
            ).unwrap();
            assert_eq!(bsdf.event, DielectricEvent::Transmission);
            assert!(bsdf.value > 0.0 && bsdf.pdf > 0.0);
            let radiance = 1.0 + 4.0 * (wavelength - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
            let irradiance = radiance * 0.27 * cosine * cosine / distance_squared
                * (-0.15 * distance).exp();
            let p_light = distance_squared / (0.27 * cosine);
            let light_weight = p_light * p_light / (p_light * p_light + bsdf.pdf * bsdf.pdf);
            // The native scalar BSDF already has the radiance eta convention.
            // No extra light-side n^2 factor belongs in this camera-side result.
            let expected = (-0.1 * camera_length).exp() * bsdf.value * irradiance * light_weight;
            near(result.area_radiance_w_m2_sr_nm()[lane], expected);
            assert_eq!(result.area_samples()[lane], 1);
            assert_eq!(result.images()[lane], 1);
            assert_eq!(
                result.pixel_sample().spectral_samples()[lane].path(),
                ordinary.spectral_samples()[lane].path(),
            );
        }
    });
}

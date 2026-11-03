//! Area-light camera radiometry and complementary emitter-hit regressions.

use super::*;
use crate::charts::{Ray, TriMesh};
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::{BoundarySide, GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionLaunch, BundleConnectionOptions, BundleConnectionTarget,
    BundleSurfaceSource,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::camera::area_lights::CurvedAreaEmitter;
use crate::spectral::SpectralReflectance;
use crate::surface_detail::{BrushingPattern, SurfaceDetail};
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Point3, Vec3};

const RHO: f64 = 0.6;
const LE: f64 = 4.0;
const AREA: f64 = 0.27;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey { seed: 61, kernel_id: 12, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn panel(id: u64, x: f64, z: f64, scale: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("native-area-camera-panel", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[x - scale, -scale, z], [x + scale, -scale, z], [x, 2.0 * scale, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn diffuse(id: u64, medium: usize) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::DepolarizingDiffuse { medium, reflectance: 0 },
    }
}

fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding { object_id: id, interaction: Interaction::Opaque }
}

fn emitter(id: u64) -> SurfaceEmitter {
    SurfaceEmitter::new(id, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)
}

fn unit(v: Vec3) -> Vec3 {
    let length = v.dot(v).sqrt();
    v.scale(1.0 / length)
}

fn seed(direction: Vec3) -> BundleImageSeed<'static> {
    BundleImageSeed {
        initial_direction: unit(direction),
        policy: BundleBranchPolicy::TransmitOrReflectTir,
    }
}

fn image_options(max_images: usize) -> BundleImageOptions {
    BundleImageOptions {
        max_images,
        min_separation_rad: 0.01,
        duplicate_tolerance_rad: 1e-5,
        max_trace_evaluations: 512,
    }
}

fn settings() -> AreaLightCameraSettings {
    let surfaces = CurvedTraceSettings {
        max_surface_events: 12,
        max_distance_m: 12.0,
        max_total_steps: 32_768,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-12,
                relative_tolerance: 1e-11,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    };
    let mut connection = BundleConnectionOptions::default();
    connection.scene.trace = surfaces;
    connection.scene.position_derivative_tolerance_m = 1e-10;
    connection.scene.direction_derivative_tolerance = 1e-11;
    AreaLightCameraSettings {
        surfaces,
        connection,
        max_area_samples: 8,
        max_trace_evaluations: 4096,
    }
}

fn camera_ray() -> Ray {
    Ray { origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-14),
        "{actual:.17e} != {expected:.17e}",
    );
}

fn metal() -> ConductorOptics {
    let source = ConductorSource::try_new(fs_blake3::ContentHash([91; 32]),
        ConductorDataStatus::Representative).unwrap();
    ConductorOptics::try_new(core::array::from_fn(|i| {
        ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
    }), source).unwrap()
}

fn oriented_panel(id: u64, center: Point3, normal: Vec3, scale: f64) -> GeometryInstance {
    let normal = unit(normal);
    let tangent = unit(Vec3::new(normal.z, 0.0, -normal.x));
    let bitangent = Vec3::new(normal.y * tangent.z - normal.z * tangent.y,
        normal.z * tangent.x - normal.x * tangent.z,
        normal.x * tangent.y - normal.y * tangent.x);
    let point = |u, v| [center.x + scale * (u * tangent.x + v * bitangent.x),
        center.y + scale * (u * tangent.y + v * bitangent.y),
        center.z + scale * (u * tangent.z + v * bitangent.z)];
    GeometryInstance::try_new(id, hash_domain("polarized-rough-area-panel", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vec![point(-1.0, -1.0), point(1.0, -1.0), point(0.0, 2.0)],
            vec![[0, 1, 2]])), RigidTransform::identity()).unwrap()
}

#[test]
fn g1_rough_area_connection_preserves_non_coplanar_polarized_camera_prefix() {
    with_cx(|cx| {
        let mirror_normal = unit(Vec3::new(-1.0, 0.0, 1.0));
        let receiver = Point3::new(1.0, 0.0, 1.0);
        let source = Point3::new(-1.0, 1.2, 1.8);
        let incoming = unit(receiver.delta_from(source));
        let outgoing = incoming.scale(-1.0);
        let geometry = InstanceScene::try_new(vec![
            oriented_panel(1, Point3::new(0.0, 0.0, 1.0), mirror_normal, 0.15),
            oriented_panel(2, receiver, Vec3::new(1.0, 0.0, 0.0), 0.5),
            oriented_panel(9, source, incoming, 0.3),
        ]).unwrap();
        let material = glass(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [metal()];
        let alpha = 0.3;
        let scene = CurvedPolarizedScene::try_new_with_materials(cx, &geometry, &media,
            &conductors, &[], &[], vec![
                SurfaceBinding { object_id: 1, interaction: Interaction::SmoothConductor { optics: 0, medium: 0 } },
                SurfaceBinding { object_id: 2, interaction: Interaction::RoughConductor {
                    optics: 0, medium: 0, surface: ConductorSurface::try_rough(alpha).unwrap(),
                }}, opaque(9),
            ], 100_000).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(incoming)];
        let bindings = [AreaLightBinding {
            source: &light, receiver_object_id: 2, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
            settings().surfaces, |_| [0.0; 3]).unwrap();
        let result = owner.sample_ray_with_samples(cx, camera_ray(), 550.0, 0, settings(), &mut workspace,
            |_| [0.0; 3], |event, binding| {
                assert_eq!((event, binding), (1, 0)); [4.0 / 9.0, 0.5]
            }).unwrap();
        // An independent native smooth-mirror prefix and facet Mueller
        // reflection, followed by analytic GGX density/masking and flat source
        // geometry. Neither the connection evaluator nor its MIS helper is
        // used by this reference.
        let prefix = camera.launch(camera_ray(), 550.0, 0).unwrap().advance(1.0).unwrap()
            .reflect_conductor(mirror_normal, &conductors[0]).unwrap().beam.advance(1.0).unwrap();
        let facet = unit(Vec3::new(-1.0 + outgoing.x, outgoing.y, outgoing.z));
        let reflected = prefix.reflect_conductor(facet, &conductors[0]).unwrap().beam;
        let cosine = -outgoing.x;
        let half_cosine = -facet.x;
        let denominator = half_cosine * half_cosine * (alpha * alpha - 1.0) + 1.0;
        let d = alpha * alpha / (core::f64::consts::PI * denominator * denominator);
        let masking = 2.0 * cosine /
            (cosine + (alpha * alpha + (1.0 - alpha * alpha) * cosine * cosine).sqrt());
        let distance2 = receiver.delta_from(source).dot(receiver.delta_from(source));
        let p_light = distance2 / AREA;
        let p_bsdf = d / 4.0;
        let weight = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
        let irradiance = LE * AREA * cosine / distance2 * (-0.1 * distance2.sqrt()).exp();
        let expected = reflected.power().i * d * masking / (4.0 * cosine) * irradiance * weight;
        near(result.area_radiance_w_m2_sr_nm(), expected, 1e-5);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), expected, 1e-5);
        assert_eq!(baseline.radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(result.radiance_sample().path(), baseline.path());
        let depolarized = PolarizedBeam::try_new(prefix.ray().origin, prefix.frame(), 550.0, prefix.medium(),
            Stokes::unpolarized(prefix.power().i)).unwrap()
            .reflect_conductor(facet, &conductors[0]).unwrap().beam;
        assert!((reflected.power().i - depolarized.power().i).abs() > 1e-5 * reflected.power().i,
            "depolarizing the camera prefix must change this non-coplanar metallic response");
        assert_eq!((result.area_samples(), result.images()), (1, 1));
    });
}

#[test]
fn g1_rough_and_brushed_receivers_share_native_continuation_and_complementary_mis() {
    with_cx(|cx| {
        // Offset the radial brush from its native pole so this exercises the
        // anisotropic finish rather than its documented isotropic fallback.
        let x = 0.2;
        let ray = Ray { origin: Point3::new(x, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) };
        let geometry = InstanceScene::try_new(vec![panel(1, x, 1.0, 1.0), panel(9, x, -1.0, 0.3)]).unwrap();
        let material = glass(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [metal()];
        let detail = SurfaceDetail::try_new(0.12, 3.0, BrushingPattern::Radial, 0.0, None).unwrap();
        let interactions = [
            (Interaction::RoughConductor { optics: 0, medium: 0,
                surface: ConductorSurface::try_rough(0.0001).unwrap() }, (0.0001, 0.0001)),
            (Interaction::RoughConductor { optics: 0, medium: 0,
                surface: ConductorSurface::try_rough(0.2).unwrap() }, (0.2, 0.2)),
            (Interaction::BrushedConductor { optics: 0, medium: 0, detail }, detail.alpha_pair()),
        ];
        for (interaction, (alpha_x, alpha_y)) in interactions {
            let scene = CurvedPolarizedScene::try_new_with_materials(cx, &geometry, &media,
                &conductors, &[], &[], vec![SurfaceBinding { object_id: 1, interaction }, opaque(9)],
                100_000).unwrap();
            let emitters = [emitter(9)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
            let bindings = [AreaLightBinding {
                source: &light, receiver_object_id: 1, source_side: BoundarySide::Positive,
                initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let baseline = camera.sample_ray_with_samples(cx, ray, 550.0, 0,
                settings().surfaces, |_| [0.0; 3]).unwrap();
            let mut surface_draws = 0;
            let mut area_draws = 0;
            let result = owner.sample_ray_with_samples(cx, ray, 550.0, 0, settings(), &mut workspace,
                |event| { assert_eq!(event, 0); surface_draws += 1; [0.0; 3] },
                |event, binding| {
                    assert_eq!((event, binding), (0, 0)); area_draws += 1; [4.0 / 9.0, 0.5]
                }).unwrap();
            // Independent normal-incidence conductor Fresnel, anisotropic GGX
            // D(0), flat receiving geometry, Beer law and power partition.
            let fresnel = ((1.0_f64 - 0.2).powi(2) + 9.0) / ((1.0_f64 + 0.2).powi(2) + 9.0);
            let p_bsdf = 1.0 / (4.0 * core::f64::consts::PI * alpha_x * alpha_y);
            let p_light = 4.0 / AREA;
            let light_weight = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
            let hit_weight = p_bsdf * p_bsdf / (p_light * p_light + p_bsdf * p_bsdf);
            let attenuation = (-0.3_f64).exp();
            let expected_nee = fresnel * p_bsdf * LE * AREA / 4.0 * attenuation * light_weight;
            let expected_hit = fresnel * LE * attenuation * hit_weight;
            near(baseline.radiance_w_m2_sr_nm(), fresnel * LE * attenuation, 2e-6);
            near(result.area_radiance_w_m2_sr_nm(), expected_nee, 3e-6);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_nee + expected_hit, 3e-6);
            assert_eq!(result.radiance_sample().path(), baseline.path());
            assert_eq!((surface_draws, area_draws), (1, 1));
            assert_eq!((result.area_samples(), result.images()), (1, 1));
            assert!(result.trace_evaluations() >= 3);
            let repeated = owner.sample_ray_with_samples(cx, ray, 550.0, 0, settings(), &mut workspace,
                |_| [0.0; 3], |_, _| [4.0 / 9.0, 0.5]).unwrap();
            assert_eq!(result.radiance_sample(), repeated.radiance_sample());
            assert_eq!(result.area_radiance_w_m2_sr_nm(), repeated.area_radiance_w_m2_sr_nm());
        }
    });
}

#[test]
fn g1_area_camera_power_mis_matches_analytic_plane_estimator() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![diffuse(1, 0), opaque(9)], 100_000,
        ).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings().surfaces, |_| [0.0; 3],
        ).unwrap();
        let mut surface_draws = 0;
        let mut area_draws = 0;
        let result = owner.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings(), &mut workspace,
            |event| { assert_eq!(event, 0); surface_draws += 1; [0.0; 3] },
            |event, binding| {
                assert_eq!((event, binding), (0, 0));
                area_draws += 1;
                [4.0 / 9.0, 0.5]
            },
        ).unwrap();
        let p_light = 4.0 / AREA;
        let p_bsdf = 1.0 / core::f64::consts::PI;
        let denominator = p_light * p_light + p_bsdf * p_bsdf;
        let raw_nee = RHO * LE * AREA / (4.0 * core::f64::consts::PI) * (-0.6_f64).exp();
        let raw_hit = RHO * LE * (-0.6_f64).exp();
        let expected_nee = raw_nee * p_light * p_light / denominator;
        let expected_hit = raw_hit * p_bsdf * p_bsdf / denominator;
        near(baseline.radiance_w_m2_sr_nm(), raw_hit, 2e-6);
        near(result.area_radiance_w_m2_sr_nm(), expected_nee, 3e-6);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_nee + expected_hit, 3e-6);
        assert_eq!(result.radiance_sample().path(), baseline.path());
        assert_eq!((surface_draws, area_draws), (1, 1));
        assert_eq!(result.area_samples(), 1);
        assert_eq!(result.images(), 1);
        assert!(result.trace_evaluations() >= 3, "NEE and emitter-hit membership both perform native work");
        let repeated = owner.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings(), &mut workspace,
            |_| [0.0; 3], |_, _| [4.0 / 9.0, 0.5],
        ).unwrap();
        assert_eq!(result.radiance_sample(), repeated.radiance_sample());
        assert_eq!(result.area_radiance_w_m2_sr_nm(), repeated.area_radiance_w_m2_sr_nm());
        assert_eq!(result.trace_evaluations(), repeated.trace_evaluations());
        assert_eq!(result.field_evaluations(), repeated.field_evaluations());
        assert_eq!(result.query_work(), repeated.query_work());
    });
}

#[test]
fn g0_area_camera_empty_bindings_and_direct_emitters_preserve_native_samples() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![diffuse(1, 0), opaque(9)], 100_000,
        ).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let empty = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings().surfaces, |_| [0.0; 3],
        ).unwrap();
        let result = empty.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings(), &mut workspace,
            |_| [0.0; 3], |_, _| panic!("empty bindings consume no area draws"),
        ).unwrap();
        assert_eq!(result.radiance_sample(), &baseline);
        assert_eq!(result.area_radiance_w_m2_sr_nm(), 0.0);
        assert_eq!((result.area_samples(), result.images(), result.trace_evaluations()), (0, 0, 0));

        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let direct = Ray { origin: camera_ray().origin, dir: Vec3::new(0.0, 0.0, -1.0) };
        let direct_baseline = camera.sample_ray_with_samples(
            cx, direct, 550.0, 0, settings().surfaces, |_| panic!("opaque emitter has no surface draw"),
        ).unwrap();
        let result = owner.sample_ray_with_samples(
            cx, direct, 550.0, 0, settings(), &mut workspace,
            |_| panic!("opaque emitter has no surface draw"),
            |_, _| panic!("direct emitter hit has no diffuse area-light vertex"),
        ).unwrap();
        assert_eq!(result.radiance_sample(), &direct_baseline);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), LE * (-0.2_f64).exp(), 2e-6);
        assert_eq!((result.area_samples(), result.images(), result.trace_evaluations()), (0, 0, 0));
    });
}

#[test]
fn g0_area_camera_rejects_foreign_sources_and_duplicate_bindings() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)]).unwrap();
        let material = glass(1.3, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![diffuse(1, 0), opaque(9)], 100_000,
        ).unwrap();
        let foreign_scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![diffuse(1, 0), opaque(9)], 100_000,
        ).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let foreign = CurvedAreaEmitter::try_new(cx, &foreign_scene, &emitters[0], 1).unwrap();
        let cloned_emitter = emitters[0].clone();
        let foreign_emitter = CurvedAreaEmitter::try_new(cx, &scene, &cloned_emitter, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        for source in [&foreign, &foreign_emitter] {
            let bindings = [AreaLightBinding {
                source, receiver_object_id: 1, source_side: BoundarySide::Positive,
                initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
            }];
            assert!(matches!(CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings),
                Err(AreaLightCameraError::InvalidBinding { index: 0 })));
        }
        let binding = || AreaLightBinding {
            source: &light, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        };
        let duplicate = [binding(), binding()];
        assert!(matches!(CurvedAreaLightCameraScene::try_new(cx, &camera, &duplicate),
            Err(AreaLightCameraError::DuplicateBinding { first: 0, second: 1 })));
    });
}

#[test]
fn g1_area_camera_unbound_emitters_keep_full_weight_and_area_budget_is_cumulative() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 0.0, 1.0, 1.0), panel(8, 2.0, -1.0, 0.3), panel(9, 0.0, -1.0, 0.3),
        ]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![diffuse(1, 0), opaque(8), opaque(9)], 100_000,
        ).unwrap();
        let emitters = [emitter(8), emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let off_axis = CurvedAreaEmitter::from_camera(cx, &camera, 8, 1).unwrap();
        let on_axis = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let off_seeds = [seed(Vec3::new(-1.0, 0.0, 1.0))];
        let on_seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let off_binding = || AreaLightBinding {
            source: &off_axis, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &off_seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        };
        let bindings = [off_binding()];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings().surfaces, |_| [0.0; 3],
        ).unwrap();
        let result = owner.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, settings(), &mut workspace,
            |_| [0.0; 3], |_, _| [4.0 / 9.0, 0.5],
        ).unwrap();
        near(result.radiance_sample().radiance_w_m2_sr_nm() - result.area_radiance_w_m2_sr_nm(),
            baseline.radiance_w_m2_sr_nm(), 2e-12);
        assert_eq!(result.radiance_sample().path(), baseline.path());
        assert_eq!(result.area_samples(), 1);
        assert!(result.area_radiance_w_m2_sr_nm() > 0.0);

        let both = [off_binding(), AreaLightBinding {
            source: &on_axis, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &on_seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let limited = CurvedAreaLightCameraScene::try_new(cx, &camera, &both).unwrap();
        let mut limits = settings();
        limits.max_area_samples = 1;
        let mut draws = 0;
        let error = limited.sample_ray_with_samples(
            cx, camera_ray(), 550.0, 0, limits, &mut workspace,
            |_| [0.0; 3], |_, _| { draws += 1; [4.0 / 9.0, 0.5] },
        );
        assert!(matches!(error, Err(AreaLightCameraError::AreaSampleLimit)));
        assert_eq!(draws, 1, "budget failure precedes a second source sample");
    });
}

fn advance(point: Point3, direction: Vec3, distance: f64) -> Point3 {
    Point3::new(
        point.x + distance * direction.x,
        point.y + distance * direction.y,
        point.z + distance * direction.z,
    )
}

fn refract(direction: Vec3, outward: Vec3, from: f64, to: f64) -> (Vec3, f64, f64) {
    let normal = if direction.dot(outward) > 0.0 { outward.scale(-1.0) } else { outward };
    let cosine = -normal.dot(direction);
    let ratio = from / to;
    let transmitted_cosine = (1.0 - ratio * ratio * (1.0 - cosine * cosine)).sqrt();
    assert!(transmitted_cosine.is_finite());
    let rs = (from * cosine - to * transmitted_cosine) / (from * cosine + to * transmitted_cosine);
    let rp = (to * cosine - from * transmitted_cosine) / (to * cosine + from * transmitted_cosine);
    let correction = ratio * cosine - transmitted_cosine;
    let outgoing = unit(Vec3::new(
        ratio * direction.x + correction * normal.x,
        ratio * direction.y + correction * normal.y,
        ratio * direction.z + correction * normal.z,
    ));
    (outgoing, 1.0 - rs * rs, 1.0 - rp * rp)
}

// Independent two-interface sphere optics: no native bundle, Jacobian,
// branch-probability accumulator or camera helper participates in this oracle.
fn ball_reference(origin: Point3, direction: Vec3, receiver_z: f64) -> (Point3, f64, f64, f64) {
    let offset = Vec3::new(origin.x, origin.y, origin.z);
    let b = offset.dot(direction);
    let entry_distance = -b - (b * b - offset.dot(offset) + 1.0).sqrt();
    assert!(entry_distance.is_finite() && entry_distance > 0.0);
    let entry = advance(origin, direction, entry_distance);
    let entry_normal = unit(Vec3::new(entry.x, entry.y, entry.z));
    let (inside, ts_in, tp_in) = refract(direction, entry_normal, 1.0, 1.5);
    let chord = -2.0 * Vec3::new(entry.x, entry.y, entry.z).dot(inside);
    let exit = advance(entry, inside, chord);
    let exit_normal = unit(Vec3::new(exit.x, exit.y, exit.z));
    let (outside, ts_out, tp_out) = refract(inside, exit_normal, 1.5, 1.0);
    let endpoint = advance(exit, outside, (receiver_z - exit.z) / outside.z);
    let proposal = 0.25 * (ts_in + tp_in) * (ts_out + tp_out);
    let physical_power = 0.5 * (ts_in * ts_out + tp_in * tp_out);
    (endpoint, chord, proposal, physical_power)
}

fn reverse_area_jacobian(origin: Point3, direction: Vec3, source_z: f64) -> f64 {
    let first = unit(Vec3::new(direction.z, 0.0, -direction.x));
    let second = Vec3::new(0.0, 1.0, 0.0);
    let h = 1e-6;
    let column = |axis: Vec3| {
        let ray = |sign: f64| unit(Vec3::new(
            direction.x + sign * h * axis.x,
            direction.y + sign * h * axis.y,
            direction.z + sign * h * axis.z,
        ));
        let plus = ball_reference(origin, ray(1.0), source_z).0;
        let minus = ball_reference(origin, ray(-1.0), source_z).0;
        Vec3::new((plus.x - minus.x) / (2.0 * h), (plus.y - minus.y) / (2.0 * h),
            (plus.z - minus.z) / (2.0 * h))
    };
    let a = column(first);
    let b = column(second);
    (a.x * b.y - a.y * b.x).abs()
}

#[test]
fn g1_area_camera_ball_lens_preserves_unseeded_images_and_reciprocal_density() {
    with_cx(|cx| {
        let receiver = Point3::new(0.03, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![
            panel(1, receiver.x, receiver.z, 1.0),
            GeometryInstance::try_new(
                2, hash_domain("native-area-camera-ball", b"unit"),
                SharedGeometry::chart(SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 }),
                RigidTransform::identity(),
            ).unwrap(),
            panel(9, 0.0, -3.0, 0.3),
        ]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&lens), SceneMedium::Homogeneous(&air)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![
                diffuse(1, 1),
                SurfaceBinding { object_id: 2, interaction: Interaction::Dielectric {
                    outside: 1, inside: 0, surface: DielectricSurface::SMOOTH,
                }},
                opaque(9),
            ], 100_000,
        ).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let area_sample = light.sample_area([4.0 / 9.0, 0.5]).unwrap();
        let seeds = [
            seed(Vec3::new((-0.18_f64).sin(), 0.0, (-0.18_f64).cos())),
            seed(Vec3::new(0.14_f64.sin(), 0.0, 0.14_f64.cos())),
        ];
        let target = BundleConnectionTarget::try_diffuse_receiver(
            1, receiver, Vec3::new(0.0, 0.0, -1.0), 2e-7,
        ).unwrap();
        let mut directions = [Vec3::new(0.0, 0.0, -1.0); 2];
        let mut nee = [0.0; 2];
        let mut hit = [0.0; 2];
        let mut complementary = [0.0; 2];
        for index in 0..2 {
            // The native solve locates each image; all density and power
            // expectations below use independent sphere optics.
            let connected = scene.connect_bundle_from_surface_with_policy(
                cx,
                BundleConnectionLaunch {
                    ray: Ray { origin: area_sample.position_m(), dir: seeds[index].initial_direction },
                    transverse_axis: Vec3::new(0.0, 1.0, 0.0), wavelength_nm: 550.0,
                    normalized_source_stokes: Stokes::unpolarized(1.0), initial_medium: 1,
                },
                BundleSurfaceSource::try_new(9, BoundarySide::Positive).unwrap(),
                target, settings().connection, seeds[index].policy,
            ).unwrap();
            let direction = connected.path().bundle().propagation().propagation().ray.dir.scale(-1.0);
            directions[index] = direction;
            let (source_point, chord, branch_probability, physical_power) = ball_reference(receiver, direction, -3.0);
            assert!(source_point.x.abs() < 2e-6 && source_point.y.abs() < 2e-6);
            let j_reverse = reverse_area_jacobian(receiver, direction, -3.0);
            let cosine = -direction.z;
            let p_light = j_reverse / AREA;
            let p_bsdf = cosine * branch_probability / core::f64::consts::PI;
            let denominator = p_light * p_light + p_bsdf * p_bsdf;
            let attenuation = (-0.2 * chord).exp();
            let raw_nee = RHO * LE * physical_power * attenuation * cosine
                / (core::f64::consts::PI * p_light);
            nee[index] = raw_nee * p_light * p_light / denominator;
            hit[index] = RHO * LE * physical_power * attenuation / branch_probability;
            complementary[index] = p_bsdf * p_bsdf / denominator;
        }
        assert!((directions[0].x - directions[1].x).abs() > 0.01);
        let outgoing = directions[1];
        // Native diffuse sampling around the facing -Z normal has tangent +Y,
        // bitangent +X. This inverts that cosine proposal without changing it.
        let diffuse_draw = [outgoing.x * outgoing.x + outgoing.y * outgoing.y,
            if outgoing.x >= 0.0 { 0.25 } else { 0.75 }, 0.5];
        let draws = |event| if event == 0 { diffuse_draw } else { [0.0, 0.0, 0.5] };
        let ray = Ray { origin: Point3::new(receiver.x, 0.0, 1.5), dir: Vec3::new(0.0, 0.0, 1.0) };
        let baseline = camera.sample_ray_with_samples(
            cx, ray, 550.0, 1, settings().surfaces, draws,
        ).unwrap();
        near(baseline.radiance_w_m2_sr_nm(), hit[1], 3e-5);
        let mut workspace = AreaLightCameraWorkspace::try_new(2, 2, 8).unwrap();
        for count in [1, 2] {
            let bindings = [AreaLightBinding {
                source: &light, receiver_object_id: 1, source_side: BoundarySide::Positive,
                initial_medium: 1, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds[..count], image_options: image_options(count), position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let result = owner.sample_ray_with_samples(
                cx, ray, 550.0, 1, settings(), &mut workspace,
                draws, |_, _| [4.0 / 9.0, 0.5],
            ).unwrap();
            let expected_nee: f64 = nee[..count].iter().sum();
            let hit_weight = if count == 1 { 1.0 } else { complementary[1] };
            near(result.area_radiance_w_m2_sr_nm(), expected_nee, 4e-4);
            near(result.radiance_sample().radiance_w_m2_sr_nm(),
                expected_nee + hit[1] * hit_weight, 4e-4);
            assert_eq!(result.radiance_sample().path(), baseline.path());
            assert_eq!(result.images(), count as u32);
            assert_eq!(result.area_samples(), 1);
        }
    });
}

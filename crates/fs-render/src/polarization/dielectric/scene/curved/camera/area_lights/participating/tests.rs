//! Native molecular NEE, matched emitter hits and independent source survival.

mod sphere;

use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, CurvedTraceSettings};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::AreaLightBinding;
use crate::polarization::dielectric::scene::curved::camera::area_lights::AreaEmitterLimits;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedDensityField, CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::MassScatteringSpectrum;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

const LE: f64 = 4.0;
const AREA: f64 = 0.27;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena,
            StreamKey { seed: 72, kernel_id: 12, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}

fn glass() -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(1.4).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(), GlassProvenance::Custom)
}

fn panel(axis: usize) -> GeometryInstance {
    let vertices = if axis == 0 {
        vec![[1.0, -0.3, 0.2], [1.0, 0.3, 0.2], [1.0, 0.0, 1.1]]
    } else {
        vec![[0.2, -0.3, 1.0], [0.8, -0.3, 1.0], [0.5, 0.6, 1.0]]
    };
    GeometryInstance::try_new(9, hash_domain("molecular-area-camera-panel", &[axis as u8]),
        SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 2, 1]])),
        RigidTransform::identity()).unwrap()
}

fn settings() -> RayleighAreaLightCameraSettings {
    let surfaces = CurvedTraceSettings {
        max_surface_events: 12, max_distance_m: 12.0, max_total_steps: 32_768,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11, direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11, optical_depth_tolerance: 1e-12,
                relative_tolerance: 1e-11, ..GradientIndexOptions::default()
            }, surface_tolerance_m: 1e-8, ..GradientEventOptions::default()
        },
    };
    let mut connection = BundleConnectionOptions::default();
    connection.scene.trace = surfaces;
    connection.scene.position_derivative_tolerance_m = 1e-10;
    connection.scene.direction_derivative_tolerance = 1e-11;
    RayleighAreaLightCameraSettings {
        area_lights: AreaLightCameraSettings {
            surfaces, connection, max_area_samples: 8, max_trace_evaluations: 4096,
        }, max_volume_events: 8, max_source_free_flight_draws: 32, max_source_candidates: 32,
    }
}

fn camera_ray() -> Ray {
    Ray { origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) }
}

fn seeds(direction: Vec3) -> [BundleImageSeed<'static>; 1] {
    [BundleImageSeed { initial_direction: direction, policy: BundleBranchPolicy::TransmitOrReflectTir }]
}

fn binding<'a, 's, 'm>(source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    seeds: &'a [BundleImageSeed<'a>]) -> MolecularAreaLightBinding<'a, 's, 'm> {
    MolecularAreaLightBinding {
        source, receiver_medium: 0, source_side: BoundarySide::Positive, initial_medium: 0,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0), seeds,
        image_options: BundleImageOptions {
            max_images: 1, min_separation_rad: 0.01,
            duplicate_tolerance_rad: 1e-5, max_trace_evaluations: 512,
        }, position_tolerance_m: 2e-7,
    }
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-5 * expected.abs().max(1e-14),
        "{actual:.17e} != {expected:.17e}");
}

fn partition(light: f64, phase: f64) -> (f64, f64) {
    let denominator = light * light + phase * phase;
    (light * light / denominator, phase * phase / denominator)
}

#[test]
fn g1_molecular_area_and_emitter_terms_use_complementary_angular_weights() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(0)]).unwrap();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let seed = seeds(Vec3::new(-1.0, 0.0, 0.0));
        let bindings = [binding(&source, &seed)];
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        for delta in [0.0, 0.03] {
            let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
            let media = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &media).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
                cx, &molecular, &surface, &bindings).unwrap();
            let isotropic = 3.0 * delta / (2.0 + delta);
            let right_angle = isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON;
            let flight = |draw| match draw {
                0 => [1.0 - (-0.35_f64).exp(), right_angle, 0.0, 0.0],
                1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
                _ => panic!("source must precede a second collision"),
            };
            let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
                CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
                |_| panic!("opaque emitter"), flight).unwrap();
            let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let mut area_draws = 0;
            let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace, |_| panic!("opaque emitter"), flight,
                |_, _| panic!("no diffuse binding"),
                |event, index| { assert_eq!((event, index), (0, 0)); area_draws += 1; [4.0 / 9.0, 0.5] },
                |_, _, _| panic!("uniform source survival is analytic"),
            ).unwrap();
            let phase_pdf = 1.5 / (4.0 * core::f64::consts::PI) * (1.0 + delta) / (2.0 + delta);
            let (light_weight, hit_weight) = partition(1.0 / AREA, phase_pdf);
            let expected_area = LE * AREA * phase_pdf * (-1.0_f64).exp() * light_weight;
            let expected_hit = LE * (-0.3_f64).exp() * hit_weight;
            near(result.area_radiance_w_m2_sr_nm(), expected_area);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + expected_hit);
            near(baseline.radiance_sample().radiance_w_m2_sr_nm(), LE * (-0.3_f64).exp());
            assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());
            assert_eq!((result.volume_events(), result.free_flight_draws()), (1, 2));
            assert_eq!((area_draws, result.area_samples(), result.images()), (1, 1, 1));
            assert_eq!((result.source_free_flight_draws(), result.source_candidates()), (0, 0));
            assert!(result.trace_evaluations() >= 4);
        }
    });
}

#[test]
fn g1_molecular_native_face_selections_keep_survival_and_hit_support_separate() {
    with_cx(|cx| {
        let source_geometry = GeometryInstance::try_new(
            9, hash_domain("molecular-area-disjoint-native-faces", b"original"),
            SharedGeometry::mesh(TriMesh::new(
                vec![
                    [1.0, -0.3, 0.2], [1.0, 0.3, 0.2], [1.0, 0.0, 1.1],
                    [1.0, 0.7, 0.2], [1.0, 1.3, 0.2], [1.0, 1.0, 1.1],
                ], vec![[0, 2, 1], [3, 5, 4]],
            )), RigidTransform::identity(),
        ).unwrap();
        let geometry = InstanceScene::try_new(vec![source_geometry]).unwrap();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let limits = AreaEmitterLimits {
            max_triangles: 2, max_preparation_work: 10_000, max_lookup_work: 10_000,
        };
        let on_axis = CurvedAreaEmitter::from_camera_faces(cx, &camera, 9, &[0], limits).unwrap();
        let off_axis = CurvedAreaEmitter::from_camera_faces(cx, &camera, 9, &[1], limits).unwrap();
        let whole = CurvedAreaEmitter::from_camera_faces(cx, &camera, 9, &[0, 1], limits).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let on_seeds = seeds(Vec3::new(-1.0, 0.0, 0.0));
        let inverse_root2 = 1.0 / 2.0_f64.sqrt();
        let off_seeds = seeds(Vec3::new(-inverse_root2, -inverse_root2, 0.0));
        let on_binding = binding(&on_axis, &on_seeds);
        let off_binding = binding(&off_axis, &off_seeds);
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let flight = |draw| match draw {
            0 => [1.0 - (-0.35_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("source must precede a second molecular collision"),
        };
        let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
            |_| panic!("opaque emitter"), flight).unwrap();
        let bindings = [off_binding, on_binding];
        let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &bindings).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let mut draws = [0; 2];
        let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            settings(), &mut workspace, |_| panic!("opaque emitter"), flight,
            |_, _| panic!("no surface receiver"), |event, index| {
                assert_eq!(event, 0);
                draws[index] += 1;
                [4.0 / 9.0, 0.5]
            }, |_, _, _| panic!("uniform source survival is analytic"),
        ).unwrap();
        let phase_pdf = 3.0 / (16.0 * core::f64::consts::PI);
        let mut expected_area = 0.0;
        let mut off_area = 0.0;
        let mut on_hit_weight = 0.0;
        for y in [0.0_f64, 1.0] {
            let distance2 = 1.0 + y * y;
            let distance = distance2.sqrt();
            let source_cosine = 1.0 / distance;
            let (light_weight, hit_weight) = partition(distance2 / (AREA * source_cosine), phase_pdf);
            // Both incoming directions are perpendicular to the original
            // camera direction. Absorption and scattering survival belong to
            // the source radiometry, not to either angular proposal density.
            let term = LE * AREA * source_cosine / distance2 * phase_pdf
                * (-0.1 - 0.9 * distance).exp() * light_weight;
            expected_area += term;
            if y == 0.0 { on_hit_weight = hit_weight; } else { off_area = term; }
        }
        let raw_hit = LE * (-0.3_f64).exp();
        near(result.area_radiance_w_m2_sr_nm(), expected_area);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + raw_hit * on_hit_weight);
        assert_eq!(draws, [1, 1]);
        assert_eq!((result.volume_events(), result.free_flight_draws()), (1, 2));
        assert_eq!((result.area_samples(), result.images()), (2, 2));
        assert_eq!((result.source_free_flight_draws(), result.source_candidates()), (0, 0));
        assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());

        let uncovered_bindings = [off_binding];
        let uncovered = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &uncovered_bindings).unwrap();
        let result = uncovered.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            settings(), &mut workspace, |_| panic!("opaque emitter"), flight,
            |_, _| panic!("no surface receiver"), |_, _| [4.0 / 9.0, 0.5],
            |_, _, _| panic!("uniform source survival is analytic"),
        ).unwrap();
        near(result.area_radiance_w_m2_sr_nm(), off_area);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), off_area + raw_hit);
        assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());

        let overlap = [on_binding, binding(&whole, &off_seeds)];
        assert!(matches!(CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &overlap),
            Err(AreaLightCameraError::DuplicateBinding { first: 0, second: 1 })));
    });
}

struct HalfDensity;
impl CurvedDensityField for HalfDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> { Ok(0.5) }
}

#[test]
fn g1_rough_surface_receivers_keep_noisy_survival_out_of_emitter_hit_partition() {
    with_cx(|cx| {
        let panel_at = |id: u64, z: f64, scale: f64| GeometryInstance::try_new(id,
            hash_domain("rough-participating-area-panel", &id.to_le_bytes()),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-scale, -scale, z], [scale, -scale, z], [0.0, 2.0 * scale, z]],
                vec![[0, 1, 2]])), RigidTransform::identity()).unwrap();
        let geometry = InstanceScene::try_new(vec![panel_at(1, 1.0, 1.0), panel_at(9, -1.0, 0.3)]).unwrap();
        let material = glass();
        let other = DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(), GlassProvenance::Custom);
        let media = [SceneMedium::Homogeneous(&material), SceneMedium::Homogeneous(&other)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 1, interaction: Interaction::Dielectric {
                outside: 1, inside: 0, surface: DielectricSurface::try_rough(0.2).unwrap(),
            }}, SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seed = seeds(Vec3::new(0.0, 0.0, 1.0));
        let proposal = binding(&source, &seed);
        let bindings = [AreaLightBinding {
            source: &source, receiver_object_id: 1, source_side: proposal.source_side,
            initial_medium: 0, source_transverse_axis: proposal.source_transverse_axis,
            seeds: &seed, image_options: proposal.image_options,
            position_tolerance_m: proposal.position_tolerance_m,
        }];
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new_heterogeneous(&spectrum, &HalfDensity,
            1.0, ContentHash([84; 32]), &phase).unwrap()), None];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &surface).unwrap();
        let flight = |_| [0.99, 0.5, 0.0, 0.9];
        let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
            |_| [0.0; 3], flight).unwrap();
        let fresnel = ((1.4_f64 - 1.0) / (1.4 + 1.0)).powi(2);
        let pdf = fresnel / (4.0 * core::f64::consts::PI * 0.2 * 0.2);
        let (light_weight, hit_weight) = partition(4.0 / AREA, pdf);
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let mut hit_terms = [0.0; 2];
        for (slot, candidate) in [false, true].into_iter().enumerate() {
            let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace, |_| [0.0; 3], flight,
                |event, index| { assert_eq!((event, index), (0, 0)); [4.0 / 9.0, 0.5] },
                |_, _| panic!("no real molecular vertex"),
                |area, image, draw| {
                    assert_eq!((area, image), (0, 0));
                    if candidate && draw == 0 { 1.0 - (-0.25_f64).exp() } else { 0.99 }
                }).unwrap();
            let ratio = if candidate { 0.5 } else { 1.0 };
            let expected_area = LE * AREA / 4.0 * pdf * (-0.6_f64).exp() * ratio * light_weight;
            let expected_hit = LE * (-0.6_f64).exp() * hit_weight;
            near(result.area_radiance_w_m2_sr_nm(), expected_area);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + expected_hit);
            assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());
            assert_eq!(result.volume_events(), 0);
            assert_eq!((result.area_samples(), result.images()), (1, 1));
            assert_eq!(result.source_candidates(), u32::from(candidate));
            hit_terms[slot] = result.radiance_sample().radiance_w_m2_sr_nm()
                - result.area_radiance_w_m2_sr_nm();
        }
        near(hit_terms[0], hit_terms[1]);
    });
}

#[test]
fn g1_null_candidates_preserve_molecular_anchor_and_noisy_survival_does_not_change_mis() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(0)]).unwrap();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let seed = seeds(Vec3::new(-1.0, 0.0, 0.0));
        let bindings = [binding(&source, &seed)];
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new_heterogeneous(&spectrum, &HalfDensity,
            1.0, ContentHash([83; 32]), &phase).unwrap())];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &bindings).unwrap();
        let flight = |draw| match draw {
            0 => [1.0 - (-0.5_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - (-0.25_f64).exp(), 0.1, 0.4, 0.9],
            2 => [0.99, 0.5, 0.0, 0.9],
            _ => panic!("unexpected camera flight"),
        };
        let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
            |_| panic!("opaque emitter"), flight).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let mut outcomes = Vec::new();
        for one_candidate in [false, true] {
            let mut source_draws = 0;
            let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace, |_| panic!("opaque emitter"), flight,
                |_, _| panic!("no diffuse binding"),
                |event, binding| { assert_eq!((event, binding), (0, 0)); [4.0 / 9.0, 0.5] },
                |area, image, draw| {
                    assert_eq!((area, image), (0, 0)); source_draws += 1;
                    if one_candidate && draw == 0 { 1.0 - (-0.25_f64).exp() } else { 0.99 }
                },
            ).unwrap();
            let phase_pdf = 3.0 / (16.0 * core::f64::consts::PI);
            let (light_weight, hit_weight) = partition(1.0 / AREA, phase_pdf);
            let ratio = if one_candidate { 0.5 } else { 1.0 };
            let expected_area = LE * AREA * phase_pdf * (-0.3_f64).exp() * ratio * light_weight;
            let expected_hit = LE * (-0.3_f64).exp() * hit_weight;
            near(result.area_radiance_w_m2_sr_nm(), expected_area);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + expected_hit);
            assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());
            assert_eq!((result.volume_events(), result.free_flight_draws()), (1, 3));
            assert_eq!((result.area_samples(), result.images()), (1, 1));
            assert_eq!(result.source_free_flight_draws(), source_draws);
            assert_eq!(result.source_candidates(), if one_candidate { 1 } else { 0 });
            outcomes.push(result);
        }
        near(outcomes[0].area_radiance_w_m2_sr_nm(), 2.0 * outcomes[1].area_radiance_w_m2_sr_nm());
        near(outcomes[0].radiance_sample().radiance_w_m2_sr_nm() - outcomes[0].area_radiance_w_m2_sr_nm(),
            outcomes[1].radiance_sample().radiance_w_m2_sr_nm() - outcomes[1].area_radiance_w_m2_sr_nm());
    });
}

#[test]
fn g1_second_molecular_vertex_retains_polarized_prefix_and_cumulative_area_limit() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(2)]).unwrap();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let seed = seeds(Vec3::new(0.0, 0.0, -1.0));
        let bindings = [binding(&source, &seed)];
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap())];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &bindings).unwrap();
        let flight = |draw| match draw {
            0 => [1.0 - (-0.25_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - (-0.5_f64).exp(), 0.5 + f64::EPSILON, 0.5, 0.0],
            2 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("unexpected camera flight"),
        };
        let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
            |_| panic!("opaque emitter"), flight).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            settings(), &mut workspace, |_| panic!("opaque emitter"), flight,
            |_, _| panic!("no diffuse binding"), |_, index| { assert_eq!(index, 0); [4.0 / 9.0, 0.5] },
            |_, _, _| panic!("uniform survival"),
        ).unwrap();
        let p90 = 3.0 / (16.0 * core::f64::consts::PI);
        let first_distance = (0.5_f64 * 0.5 + 0.75 * 0.75).sqrt();
        let first_cosine = 0.75 / first_distance;
        let first_phase = p90 * (1.0 + first_cosine * first_cosine);
        let first_light_pdf = first_distance * first_distance / (AREA * first_cosine);
        let (first_weight, _) = partition(first_light_pdf, first_phase);
        let first = LE * AREA * first_cosine / (first_distance * first_distance)
            * first_phase * (-0.05 - 1.2 * first_distance).exp() * first_weight;
        let (second_weight, hit_weight) = partition(0.75 * 0.75 / AREA, p90);
        // Both scattering planes are xz. The first 90-degree event completely
        // polarizes along y; the second response is twice the scalar phase.
        // Its proposal is still p90, so using 2*p90 in MIS would be incorrect.
        let second = LE * AREA / (0.75 * 0.75) * (2.0 * p90)
            * (-1.05_f64).exp() * second_weight;
        let hit = 2.0 * LE * (-0.3_f64).exp() * hit_weight;
        near(result.area_radiance_w_m2_sr_nm(), first + second);
        near(result.radiance_sample().radiance_w_m2_sr_nm(), first + second + hit);
        near(baseline.radiance_sample().radiance_w_m2_sr_nm(), 2.0 * LE * (-0.3_f64).exp());
        assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());
        assert_eq!((result.volume_events(), result.free_flight_draws(), result.area_samples()), (2, 3, 2));
        let mut limited = settings();
        limited.area_lights.max_area_samples = 1;
        assert!(matches!(owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            limited, &mut workspace, |_| panic!("opaque emitter"), flight,
            |_, _| panic!("no diffuse binding"), |_, _| [4.0 / 9.0, 0.5],
            |_, _, _| panic!("uniform survival")), Err(AreaLightCameraError::AreaSampleLimit)));
    });
}

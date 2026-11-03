use super::*;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance,
};
#[cfg(feature = "volumes")]
use crate::dielectric::DielectricSurface;
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::GradientIndexOptions;
use crate::instances::{InstanceScene, RigidTransform};
use crate::polarization::dielectric::scene::curved::{CurvedTraceSettings, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleBranchPolicy;
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
#[cfg(feature = "volumes")]
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighMedium;
#[cfg(feature = "volumes")]
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_geom::fixtures::SphereChart;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey { seed: 43, kernel_id: 8, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(), GlassProvenance::Custom)
}

fn mesh(id: u64, vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>, transform: RigidTransform)
    -> GeometryInstance {
    GeometryInstance::try_new(id, hash_domain("native-area-emitter-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices, triangles)), transform).unwrap()
}

fn source(z: f64) -> GeometryInstance {
    mesh(2, vec![[-0.3, -0.3, z], [0.3, -0.3, z], [0.0, 0.6, z]],
        vec![[0, 1, 2]], RigidTransform::identity())
}

fn detector(point: Point3) -> GeometryInstance {
    mesh(9, vec![[point.x - 0.5, -0.5, point.z], [point.x + 0.5, -0.5, point.z],
        [point.x, 0.5, point.z]], vec![[0, 1, 2]], RigidTransform::identity())
}

fn detector_binding() -> SurfaceBinding {
    SurfaceBinding { object_id: 9,
        interaction: Interaction::Detector { reference_axis: Vec3::new(0.0, 1.0, 0.0) } }
}

fn options() -> BundleConnectionOptions {
    let mut options = BundleConnectionOptions::default();
    options.scene.trace = CurvedTraceSettings {
        max_surface_events: 8, max_distance_m: 8.0,
        max_total_steps: 16_384, max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11, direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11, optical_depth_tolerance: 1e-12,
                relative_tolerance: 1e-11, ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8, ..GradientEventOptions::default()
        },
    };
    options.scene.position_derivative_tolerance_m = 1e-10;
    options.scene.direction_derivative_tolerance = 1e-11;
    options
}

fn seed(angle: f64) -> BundleImageSeed<'static> {
    BundleImageSeed { initial_direction: Vec3::new(angle.sin(), 0.0, angle.cos()),
        policy: BundleBranchPolicy::TransmitOrReflectTir }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs().max(1e-14),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn g0_original_panel_area_pdf_and_winding_survive_rigid_placement() {
    with_cx(|cx| {
        let half = 0.5_f64.sqrt();
        let placement = RigidTransform::try_new([0.0, half, 0.0, half], [1.0, 2.0, 3.0]).unwrap();
        let geometry = InstanceScene::try_new(vec![mesh(2,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -2.0, 0.0]],
            vec![[0, 1, 2], [1, 0, 3]], placement)]).unwrap();
        let material = glass(1.3, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front);
        let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 2).unwrap();
        near(panel.total_area_m2(), 3.0, 1e-15);
        for (u, triangle, root, a, b, c) in [
            (0.1, 0, 0.3_f64.sqrt(), [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            (0.8, 1, 0.7_f64.sqrt(), [2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, -2.0, 0.0]),
        ] {
            let sampled = panel.sample_area([u, 0.25]).unwrap();
            assert_eq!(sampled.triangle_index(), Some(triangle));
            let contextual = panel.sample_area_with_context(cx, [u, 0.25]).unwrap();
            assert_eq!(contextual.triangle_index(), sampled.triangle_index());
            assert_eq!(contextual.position_m(), sampled.position_m());
            assert_eq!(contextual.normal(), sampled.normal());
            assert_eq!(contextual.pdf_area_per_m2(), sampled.pdf_area_per_m2());
            assert_eq!(sampled.object_id(), 2);
            near(sampled.pdf_area_per_m2(), 1.0 / 3.0, 1e-14);
            let expected_local = Point3::new(
                a[0] * (1.0 - root) + b[0] * root * 0.75 + c[0] * root * 0.25,
                a[1] * (1.0 - root) + b[1] * root * 0.75 + c[1] * root * 0.25, 0.0);
            let expected = placement.transform_point(expected_local);
            near(sampled.position_m().x, expected.x, 1e-14);
            near(sampled.position_m().y, expected.y, 1e-14);
            near(sampled.position_m().z, expected.z, 1e-14);
            let normal = sampled.normal();
            assert!((normal.x - 1.0).abs() < 1e-13 && normal.y.abs() < 1e-13 && normal.z.abs() < 1e-13);
        }
        assert!(matches!(CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1),
            Err(AreaEmitterError::TriangleBudget { required: 2, limit: 1 })));
        assert!(matches!(panel.sample_area([0.0, 0.5]), Err(AreaEmitterError::InvalidSample)));
        let overlapping = TriMesh::new(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.5, 0.0]],
            vec![[0, 1, 2], [1, 0, 3]]);
        assert!(matches!(overlapping.validate_area_source_faces(cx, &[0, 1], 100),
            Err(MeshAreaError::OverlappingFaces { .. })));
    });
}

#[test]
fn g1_closed_native_mesh_and_conditional_faces_keep_original_area_and_transport() {
    with_cx(|cx| {
        // A tetrahedron with four noncoplanar original faces of areas
        // 1, 3, 1.5, 3.5 m². Each face retains its original outward winding.
        let tetrahedron = mesh(2,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 3.0]],
            vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], RigidTransform::identity());
        let point = Point3::new(2.0 / 3.0, 1.0 / 3.0, -4.0);
        let geometry = InstanceScene::try_new(vec![tetrahedron, detector(point)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }, detector_binding()],
            100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front);
        let whole = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 4).unwrap();
        near(whole.total_area_m2(), 9.0, 1e-14);
        let normals = [Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0), Vec3::new(3.0 / 7.0, 6.0 / 7.0, 2.0 / 7.0)];
        let mut lower = 0.0;
        for (face, area) in [1.0, 3.0, 1.5, 3.5].into_iter().enumerate() {
            let sample = whole.sample_area([(lower + area * 4.0 / 9.0) / 9.0, 0.5]).unwrap();
            assert_eq!(sample.triangle_index(), Some(face));
            near(sample.pdf_area_per_m2(), 1.0 / 9.0, 1e-14);
            assert!(sample.normal().dot(normals[face]) > 1.0 - 1e-14);
            lower += area;
        }
        let limits = AreaEmitterLimits { max_triangles: 4, max_preparation_work: 1000, max_lookup_work: 100 };
        let even = CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[0, 2], limits).unwrap();
        let odd = CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[1, 3], limits).unwrap();
        assert!(!even.support_overlaps(&odd));
        assert!(whole.support_overlaps(&odd));
        near(odd.total_area_m2(), 6.5, 1e-14);
        let selected = odd.sample_area([0.8, 0.5]).unwrap();
        assert_eq!(selected.triangle_index(), Some(3));
        near(selected.pdf_area_per_m2(), 1.0 / 6.5, 1e-14);
        assert!(matches!(CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[0],
            AreaEmitterLimits { max_triangles: 3, ..limits }),
            Err(AreaEmitterError::TriangleBudget { required: 4, limit: 3 })));
        for faces in [&[][..], &[1, 0][..], &[0, 0][..], &[4][..]] {
            assert!(matches!(CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, faces, limits),
                Err(AreaEmitterError::Mesh(MeshAreaError::InvalidSupport))));
        }
        let front = CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[0], limits).unwrap();
        let seeds = [BundleImageSeed { initial_direction: Vec3::new(0.0, 0.0, -1.0),
            policy: BundleBranchPolicy::Sequence(&[]) }];
        assert!(matches!(front.launch(whole.sample_area([4.0 / 81.0, 0.5]).unwrap(),
            BoundarySide::Positive, 0, Vec3::new(0.0, 1.0, 0.0), 550.0, &seeds),
            Err(AreaEmitterError::SampleOwnerMismatch)));
        let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
        let target = BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, 1.0), 2e-7).unwrap();
        let result = front.incident_irradiance(cx, front.sample_area([4.0 / 9.0, 0.5]).unwrap(),
            BoundarySide::Positive, 0, Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(), &seeds,
            BundleImageOptions { max_images: 1, ..BundleImageOptions::default() }, &mut workspace).unwrap();
        near(result.spectral_irradiance_w_m2_nm(), 4.0 / 16.0 * (-0.8_f64).exp(), 2e-6);
        assert_eq!(result.image_count(), 1);
    });
}

#[test]
fn g0_mesh_source_refuses_excluded_coplanar_overlap_and_preparation_exhaustion() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![mesh(2,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
                [0.25, 0.25, 0.0], [2.25, 0.25, 0.0], [0.25, 2.25, 0.0]],
            vec![[0, 1, 2], [3, 4, 5]], RigidTransform::identity())]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, super::super::RadianceSpectrum::try_constant(1.0).unwrap(),
            EmissionSides::Front);
        let limits = AreaEmitterLimits { max_triangles: 2, max_preparation_work: 100, max_lookup_work: 100 };
        assert!(matches!(CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[0], limits),
            Err(AreaEmitterError::Mesh(MeshAreaError::OverlappingFaces { first: 0, second: 1 }))));
        assert!(matches!(CurvedAreaEmitter::try_new_faces(cx, &scene, &emitter, &[0],
            AreaEmitterLimits { max_preparation_work: 1, ..limits }),
            Err(AreaEmitterError::Mesh(MeshAreaError::WorkLimit { limit: 1 }))));
    });
}

#[test]
fn g1_area_source_detector_and_transverse_volume_irradiance_keep_physical_units() {
    with_cx(|cx| {
        let point = Point3::new(0.0, 0.0, 4.0);
        let geometry = InstanceScene::try_new(vec![source(0.0), detector(point)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque },
                detector_binding()], 100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front);
        let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
        let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
        near(sample.pdf_area_per_m2(), 1.0 / 0.27, 1e-14);
        let seeds = [seed(0.0); 2];
        let images = BundleImageOptions { max_images: 1, ..BundleImageOptions::default() };
        let mut workspace = BundleImageWorkspace::try_new(1, 2).unwrap();
        let target = BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
        let detector = panel.incident_irradiance(cx, sample, BoundarySide::Positive, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(), &seeds, images, &mut workspace).unwrap();
        // Le * source cosine * area / distance² * Beer. The surrounding
        // refractive index 1.3 introduces no extra source-leg n² multiplier.
        near(detector.spectral_irradiance_w_m2_nm(), 4.0 * 0.27 / 16.0 * (-0.8_f64).exp(), 2e-6);
        assert_eq!(detector.image_count(), 1);
        assert_eq!(detector.trace_evaluations(), 2);
        let volume = BundleConnectionTarget::try_volume_receiver(0, Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.6, 0.0, -0.8), 2e-7).unwrap();
        let transverse = panel.incident_irradiance(cx, sample, BoundarySide::Positive, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, volume, options(), &seeds, images, &mut workspace).unwrap();
        near(transverse.spectral_irradiance_w_m2_nm(), 4.0 * 0.27 / 4.0 * (-0.4_f64).exp(), 2e-6);
        assert!(matches!(panel.incident_irradiance(cx, sample, BoundarySide::Negative, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, volume, options(), &seeds, images, &mut workspace),
            Err(AreaEmitterError::UnsupportedEmissionSide)));
        let both = SurfaceEmitter::new(2, emitter.spectrum().clone(), EmissionSides::Both);
        let both_panel = CurvedAreaEmitter::try_new(cx, &scene, &both, 1).unwrap();
        assert!(matches!(both_panel.incident_irradiance(cx, sample, BoundarySide::Positive, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, volume, options(), &seeds, images, &mut workspace),
            Err(AreaEmitterError::SampleOwnerMismatch)));
        let negative = BundleConnectionTarget::try_volume_receiver(0, Point3::new(0.0, 0.0, -2.0),
            Vec3::new(0.0, 0.0, 1.0), 2e-7).unwrap();
        let negative_seeds = [BundleImageSeed { initial_direction: Vec3::new(0.0, 0.0, -1.0),
            policy: BundleBranchPolicy::Sequence(&[]) }; 2];
        let backward = both_panel.incident_irradiance(cx,
            both_panel.sample_area([4.0 / 9.0, 0.5]).unwrap(), BoundarySide::Negative, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, negative, options(), &negative_seeds,
            images, &mut workspace).unwrap();
        near(backward.spectral_irradiance_w_m2_nm(),
            transverse.spectral_irradiance_w_m2_nm(), 2e-6);
    });
}

#[cfg(feature = "volumes")]
#[test]
fn g1_area_lens_images_sum_original_branch_power_with_independent_source_survival() {
    with_cx(|cx| {
        let point = Point3::new(0.03, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(1, hash_domain("native-area-emitter-test", b"unit-lens"),
                SharedGeometry::chart(SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 }),
                RigidTransform::identity()).unwrap(), source(-3.0), detector(point)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&lens), SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![
            SurfaceBinding { object_id: 1, interaction: Interaction::Dielectric {
                outside: 1, inside: 0, surface: DielectricSurface::SMOOTH } },
            SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }, detector_binding()],
            100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, super::super::RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Front);
        let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
        let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
        let target = BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
        let seeds = [seed(-0.18), seed(0.065), seed(0.14), seed(-0.18)];
        let spectrum = MassScatteringSpectrum::constant(0.3).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap()), None];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let mut expected = 0.0;
        let mut expected_traces = 0;
        for (index, seed) in seeds[..3].iter().enumerate() {
            let connected = scene.connect_bundle_from_surface_with_policy(cx,
                BundleConnectionLaunch { ray: Ray { origin: sample.position_m(), dir: seed.initial_direction },
                    transverse_axis: Vec3::new(0.0, 1.0, 0.0), wavelength_nm: 550.0,
                    normalized_source_stokes: Stokes::unpolarized(1.0), initial_medium: 1 },
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                target, options(), seed.policy).unwrap();
            let bundle = connected.path().bundle();
            let cosine = bundle.propagation().initial_ray().dir.z;
            let intensity = 3.0 * cosine / sample.pdf_area_per_m2();
            let clear = bundle.point_source_irradiance(target.normal(), options().projection, intensity).unwrap();
            // Only the lens absorbs/scatters. Its original Beer depth gives
            // its physical internal chord independently of survival replay.
            let chord = bundle.propagation().propagation().optical_depth / 0.2;
            expected += clear * (-0.3 * chord).exp();
            expected_traces += connected.trace_evaluations() + 1;
            if index == 0 { expected_traces += connected.trace_evaluations(); }
        }
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(3, 4).unwrap();
        let images = BundleImageOptions { max_images: 3, min_separation_rad: 0.01,
            duplicate_tolerance_rad: 1e-5, max_trace_evaluations: 512 };
        let result = panel.incident_irradiance_with_transmittance_samples(cx, sample,
            BoundarySide::Positive, 1, Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(),
            &seeds, images, &scattering, CurvedTransmittanceOptions {
                max_free_flight_draws: 0, max_tracking_candidates: 0 }, &mut workspace,
            |_, _| panic!("uniform source survival consumes no draws")).unwrap();
        near(result.spectral_irradiance_w_m2_nm(), expected, 3e-6);
        assert_eq!(result.image_count(), 3);
        assert_eq!(result.trace_evaluations(), expected_traces);
        assert_eq!(result.source_free_flight_draws(), 0);
        assert_eq!(result.source_candidates(), 0);
        let repeat = panel.incident_irradiance_with_transmittance_samples(cx, sample,
            BoundarySide::Positive, 1, Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(),
            &seeds, images, &scattering, CurvedTransmittanceOptions::default(), &mut workspace,
            |_, _| panic!("uniform source survival consumes no draws")).unwrap();
        assert_eq!(result, repeat);
    });
}

mod mis;
mod chart;


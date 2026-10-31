//! G1/G3: a launch point carries its actual host from the first physical metre.

use super::*;
use fs_geom::fixtures::SphereChart;

fn cube(half: f64, center_x: f64) -> TriMesh {
    let source = box_mesh_at(0.0, half, -half, false);
    let vertices = source
        .vertices
        .iter()
        .map(|v| [center_x + v[0] * half / 4.0, v[1] * half / 4.0, v[2]])
        .collect();
    TriMesh::new(vertices, source.triangles)
}

fn interior_scene(index: f64, absorption: f64) -> Scene {
    let mut scene = scene_with_boundaries(
        &[(
            1.0,
            -1.0,
            material(index, absorption, DielectricSurface::SMOOTH),
            false,
        )],
        3.0,
        20.0,
    );
    scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
    scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
    scene
}

#[test]
fn g1_inside_camera_exit_has_correct_absorption_fresnel_and_radiance_eta() {
    with_cx(|cx| {
        let policy = settings(DirectStrategy::BsdfOnly, 16, 2);
        let baseline = render(&interior_scene(1.0, 0.0), cx, &policy).unwrap();
        let absorbed = render(&interior_scene(1.0, 0.7), cx, &policy).unwrap();
        assert_channel_ratio(
            &absorbed,
            &baseline,
            fs_math::det::exp(-0.7),
            3.0e-14,
            "inside physical absorption length",
        );
        let glass = render(&interior_scene(1.5, 0.7), cx, &policy).unwrap();
        let fresnel = ((1.5_f64 - 1.0) / (1.5 + 1.0)).powi(2);
        let expected = fs_math::det::exp(-0.7) * (1.0 - fresnel) * 1.5 * 1.5;
        assert_channel_ratio(
            &glass,
            &baseline,
            expected,
            3.0e-14,
            "glass-to-air initial exit",
        );
    });
}

#[test]
fn g1_initial_total_internal_reflection_does_not_leak_to_exterior() {
    with_cx(|cx| {
        let mut scene = interior_scene(1.5, 0.0);
        scene.camera.forward = Vec3::new(0.8, 0.0, 0.6);
        let policy = settings(DirectStrategy::BsdfOnly, 8, 2);
        let trapped = render(&scene, cx, &policy).unwrap();
        assert_eq!(trapped.xyz[0], [0.0; 3]);
        scene.primitives[0].material = material(1.0, 0.0, DielectricSurface::SMOOTH);
        let escaping = render(&scene, cx, &policy).unwrap();
        assert!(escaping.xyz[0].iter().all(|v| *v > 0.0));
    });
}

#[test]
fn g1_nested_initial_media_use_geometric_order_and_restore_outer_absorption() {
    with_cx(|cx| {
        let mut scene = interior_scene(1.0, 0.0);
        // Deliberately put the innermost primitive first. Scene order cannot
        // decide the initial medium or the host restored after its exit.
        scene.primitives[0].shape = Shape::Mesh(cube(1.0, 0.0));
        scene.primitives[0].material = material(1.0, 0.2, DielectricSurface::SMOOTH);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(cube(2.0, 0.0)),
            material: material(1.0, 0.7, DielectricSurface::SMOOTH),
            emission: None,
        });
        let policy = settings(DirectStrategy::Mis, 16, 3);
        let absorbed = render(&scene, cx, &policy).unwrap();
        scene.primitives[0].material = material(1.0, 0.0, DielectricSurface::SMOOTH);
        scene.primitives[2].material = material(1.0, 0.0, DielectricSurface::SMOOTH);
        let clear = render(&scene, cx, &policy).unwrap();
        assert_channel_ratio(
            &absorbed,
            &clear,
            fs_math::det::exp(-0.9),
            3.0e-14,
            "one metre in each nested host",
        );
    });
}

#[test]
fn g1_certified_chart_origin_uses_the_same_medium_transport() {
    with_cx(|cx| {
        let mut scene = interior_scene(1.0, 0.6);
        scene.primitives[0].shape = Shape::Chart(Box::new(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }));
        let policy = settings(DirectStrategy::BsdfOnly, 8, 2);
        let absorbed = render(&scene, cx, &policy).unwrap();
        scene.primitives[0].material = material(1.0, 0.0, DielectricSurface::SMOOTH);
        let clear = render(&scene, cx, &policy).unwrap();
        assert_channel_ratio(
            &absorbed,
            &clear,
            fs_math::det::exp(-0.6),
            1.0e-6,
            "certified signed chart interior",
        );
    });
}

#[test]
fn g3_instanced_interior_is_classified_in_the_native_local_frame() {
    use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
    with_cx(|cx| {
        let glass = material(1.0, 0.4, DielectricSurface::SMOOTH);
        let mut scene = scene_with_boundaries_at(12.0, &[(1.0, -1.0, glass, false)], 3.0, 20.0);
        scene.camera.eye = Point3::new(12.0, 0.0, 0.0);
        scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
        let identity =
            fs_blake3::hash_domain("org.frankensim.render.initial-medium-test", b"local-cube");
        scene.primitives[0].shape = Shape::Instance(
            GeometryInstance::try_new(
                91,
                identity,
                SharedGeometry::mesh(cube(1.0, 0.0)),
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [12.0, 0.0, 0.0]).unwrap(),
            )
            .unwrap(),
        );
        let policy = settings(DirectStrategy::BsdfOnly, 8, 2);
        let absorbed = render(&scene, cx, &policy).unwrap();
        scene.primitives[0].material = material(1.0, 0.0, DielectricSurface::SMOOTH);
        let clear = render(&scene, cx, &policy).unwrap();
        assert_channel_ratio(
            &absorbed,
            &clear,
            fs_math::det::exp(-0.4),
            3.0e-14,
            "translated interior",
        );
    });
}

#[test]
fn g0_initial_overlaps_boundary_points_and_open_meshes_refuse_transactionally() {
    with_cx(|cx| {
        let mut overlap = interior_scene(1.0, 0.0);
        overlap.primitives[0].shape = Shape::Mesh(cube(1.0, -0.5));
        overlap.primitives.push(Primitive {
            shape: Shape::Mesh(cube(1.0, 0.5)),
            material: material(1.0, 0.0, DielectricSurface::SMOOTH),
            emission: None,
        });
        let mut boundary = interior_scene(1.0, 0.0);
        boundary.camera.eye.z = 1.0;
        let mut open = interior_scene(1.0, 0.0);
        let source = cube(1.0, 0.0);
        open.primitives[0].shape = Shape::Mesh(TriMesh::new(
            source.vertices,
            source.triangles[..10].to_vec(),
        ));
        let mut reversed = interior_scene(1.0, 0.0);
        reversed.primitives[0].shape = Shape::Mesh(box_mesh_at(0.0, 1.0, -1.0, true));
        let policy = settings(DirectStrategy::Mis, 4, 4);
        for scene in [overlap, boundary, open, reversed] {
            let mut film = Film::new(1, 1);
            let original = film.clone();
            assert!(matches!(
                render_range(&scene, cx, &policy, &mut film, 0, 4),
                Err(TracerError::InitialMedium { .. })
            ));
            assert_eq!(film, original);
        }
    });
}

#[test]
fn g3_inside_camera_progressive_ranges_preserve_sample_bits() {
    with_cx(|cx| {
        let scene = interior_scene(1.5, 0.7);
        let policy = settings(DirectStrategy::Mis, 16, 5);
        let complete = render(&scene, cx, &policy).unwrap();
        let mut progressive = Film::new(1, 1);
        render_range(&scene, cx, &policy, &mut progressive, 0, 5).unwrap();
        render_range(&scene, cx, &policy, &mut progressive, 5, 16).unwrap();
        assert_film_bits_eq(&complete, &progressive, "interior progressive replay");
    });
}

use super::*;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena,
            fs_exec::StreamKey { seed: 59, kernel_id: 8, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic));
    });
}

#[test]
fn g1_native_mesh_lookup_distinguishes_actual_planes_and_original_faces() {
    with_cx(|cx| {
        let mesh = TriMesh::new(vec![
            [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0],
            [4.0, 0.0, 0.0], [4.0, 2.0, 0.0], [4.0, 0.0, 2.0],
        ], vec![[0, 1, 2], [3, 4, 5], [6, 7, 8]]);
        mesh.validate_area_source_faces(cx, &[0, 1, 2], 100).unwrap();
        for (point, face) in [
            (Point3::new(0.25, 0.25, 0.0), 0),
            (Point3::new(0.25, 0.25, 2.0), 1),
            (Point3::new(4.0, 0.25, 0.25), 2),
        ] {
            assert_eq!(mesh.area_source_face_at(cx, point, 1e-9, 100).unwrap(), Some(face));
        }
        assert_eq!(mesh.area_source_face_at(cx, Point3::new(0.25, 0.25, 1.0), 1e-9, 100).unwrap(), None);
        assert_eq!(mesh.area_source_face_at(cx, Point3::new(0.25, 0.0, 0.0), 1e-9, 100),
            Err(MeshAreaError::AmbiguousPoint));
        assert_eq!(mesh.area_source_face_at(cx, Point3::new(0.25, 0.25, 5e-10), 1e-9, 100).unwrap(), Some(0));
    });
}

#[test]
fn g0_contact_and_inverse_uncertainty_cannot_choose_a_conditional_face_across_an_edge() {
    with_cx(|cx| {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let square = TriMesh::new(vertices.clone(), vec![[0, 1, 2], [0, 2, 3]]);
        for delta in [-1e-12, 0.0, 1e-12] {
            assert_eq!(square.area_source_face_at(cx, Point3::new(0.5, 0.5 + delta, 0.0), 1e-7, 100),
                Err(MeshAreaError::AmbiguousPoint));
        }
        // A definitely outside edge may exclude a face. An outside point whose
        // uncertainty still reaches that face cannot silently exclude it.
        let lower = TriMesh::new(vertices, vec![[0, 1, 2]]);
        assert_eq!(lower.area_source_face_at(cx, Point3::new(0.5, 0.5 + 1e-8, 0.0), 1e-7, 100),
            Err(MeshAreaError::AmbiguousPoint));
        assert_eq!(lower.area_source_face_at(cx, Point3::new(0.5, 0.5 + 1e-5, 0.0), 1e-7, 100).unwrap(), None);
        assert_eq!(square.area_source_face_at(cx, Point3::new(-5e-8, 0.25, 0.0), 1e-7, 100),
            Err(MeshAreaError::AmbiguousPoint));
        for (point, face) in [(Point3::new(0.75, 0.25, 0.0), 0), (Point3::new(0.25, 0.75, 0.0), 1)] {
            assert_eq!(square.area_source_face_at(cx, point, 1e-7, 100).unwrap(), Some(face));
        }

        // Identity rotation with a large proper-rigid translation. Ordinary
        // world placement and inverse subtraction move this point across the
        // shared diagonal: the conditional source PDF must not follow that
        // rounding accident. Use the same placement allowance as camera MIS.
        let original = Point3::new(0.50000026, 0.50000025, 0.0);
        let translation = [1_073_741_824.0, 2_147_483_648.0, 3.0];
        let world = Point3::new(original.x + translation[0], original.y + translation[1],
            original.z + translation[2]);
        let local = Point3::new(world.x - translation[0], world.y - translation[1], world.z - translation[2]);
        assert!(original.x > original.y && local.x < local.y);
        let scale = [world.x, world.y, world.z, local.x, local.y, local.z]
            .into_iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
        let tolerance = (64.0 * f64::EPSILON * scale).next_up();
        assert_eq!(square.area_source_face_at(cx, local, tolerance, 100), Err(MeshAreaError::AmbiguousPoint));
        assert_eq!(square.area_source_face_at(cx, Point3::new(0.75, 0.25, 0.0), tolerance, 100).unwrap(), Some(0));
    });
}

#[test]
fn g0_coplanar_overlap_refuses_selected_and_excluded_faces_without_rejecting_shared_edges() {
    with_cx(|cx| {
        let square = TriMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], vec![[0, 1, 2], [0, 2, 3]]);
        square.validate_area_source_faces(cx, &[0, 1], 100).unwrap();
        let overlap = TriMesh::new(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
            [0.25, 0.25, 0.0], [2.25, 0.25, 0.0], [0.25, 2.25, 0.0]], vec![[0, 1, 2], [3, 4, 5]]);
        for faces in [&[0][..], &[1][..], &[0, 1][..]] {
            assert!(matches!(overlap.validate_area_source_faces(cx, faces, 100),
                Err(MeshAreaError::OverlappingFaces { .. })));
        }
        let duplicate = TriMesh::new(square.vertices, vec![[0, 1, 2], [0, 1, 2]]);
        assert!(matches!(duplicate.validate_area_source_faces(cx, &[0], 100),
            Err(MeshAreaError::OverlappingFaces { first: 0, second: 1 })));
        // Noncoplanar faces meet only along a line; the native contact owner
        // retains responsibility for ambiguity at that zero-area locus.
        let hinge = TriMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]], vec![[0, 1, 2], [1, 0, 3]]);
        hinge.validate_area_source_faces(cx, &[0, 1], 100).unwrap();
    });
}

#[test]
fn g3_original_bvh_prunes_sparse_emitter_lookup_and_preparation_is_bounded() {
    with_cx(|cx| {
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        for i in 0..512_u32 {
            let x = 4.0 * f64::from(i);
            vertices.extend([[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]]);
            faces.push([3 * i, 3 * i + 1, 3 * i + 2]);
        }
        let mesh = TriMesh::new(vertices, faces);
        mesh.validate_area_source_faces(cx, &[0, 511], 100).unwrap();
        assert_eq!(mesh.area_source_face_at(cx, Point3::new(0.25, 0.25, 0.0), 1e-9, 100).unwrap(), Some(0));
        assert_eq!(mesh.validate_area_source_faces(cx, &[0], 1), Err(MeshAreaError::WorkLimit { limit: 1 }));
        assert_eq!(mesh.area_source_face_at(cx, Point3::new(0.25, 0.25, 0.0), 1e-9, 1),
            Err(MeshAreaError::WorkLimit { limit: 1 }));
        for faces in [&[][..], &[1, 0][..], &[1, 1][..], &[512][..]] {
            assert_eq!(mesh.validate_area_source_faces(cx, faces, 100), Err(MeshAreaError::InvalidSupport));
        }
    });
}

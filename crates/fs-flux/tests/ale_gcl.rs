//! flux-007: affine-triangle geometric-conservation-law evidence.

use fs_flux::TriMesh;
use fs_flux::ale::{
    AleGclArithmetic2, AleGclError2, AleGclMeshIssue2, AleMeshState2, audit_affine_triangle_gcl2,
};
use fs_solid::Mesh2;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn single_triangle() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        elems: vec![vec![0, 1, 2]],
        patches: Vec::new(),
    })
}

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 64.0 * f64::EPSILON * scale.max(1.0)
}

#[test]
fn flux_007_fixed_translation_expansion_and_shear() {
    let triangle = single_triangle();
    let fixed = audit_affine_triangle_gcl2(&triangle, &triangle.verts, 0.25).unwrap();
    let zero = 0.0f64.to_bits();
    let fixed_pass = fixed.total_area_change().to_bits() == zero
        && fixed.boundary_swept_area().to_bits() == zero
        && fixed.total_area_defect().to_bits() == zero
        && fixed
            .owner_edge_swept_areas()
            .iter()
            .all(|sweep| sweep.to_bits() == zero)
        && fixed.cells().iter().all(|cell| {
            cell.outward_swept_area().to_bits() == zero && cell.area_defect().to_bits() == zero
        });

    let translated: Vec<_> = triangle
        .verts
        .iter()
        .map(|point| [point[0] + 0.25, point[1] - 0.5])
        .collect();
    let translation = audit_affine_triangle_gcl2(&triangle, &translated, 0.25).unwrap();
    let translation_pass = close(translation.total_area_change(), 0.0, 1.0)
        && close(translation.boundary_swept_area(), 0.0, 1.0)
        && close(translation.total_area_defect(), 0.0, 1.0)
        && close(translation.worst_abs_cell_defect(), 0.0, 1.0);
    verdict(
        "flux-007-fixed-and-translation",
        fixed_pass && translation_pass,
        &format!(
            "fixed_defect={:.3e} translation_defect={:.3e}",
            fixed.total_area_defect(),
            translation.total_area_defect()
        ),
    );

    let expanded: Vec<_> = triangle
        .verts
        .iter()
        .map(|point| [2.0 * point[0], 2.0 * point[1]])
        .collect();
    let expansion = audit_affine_triangle_gcl2(&triangle, &expanded, 0.5).unwrap();
    let diagonal = triangle
        .edges
        .iter()
        .position(|edge| edge.verts == (1, 2))
        .unwrap();
    let expansion_pass = expansion.cells()[0].area_before().to_bits() == 0.5f64.to_bits()
        && expansion.cells()[0].area_after().to_bits() == 2.0f64.to_bits()
        && expansion.cells()[0].minimum_area().to_bits() == 0.5f64.to_bits()
        && expansion.total_area_change().to_bits() == 1.5f64.to_bits()
        && expansion.boundary_swept_area().to_bits() == 1.5f64.to_bits()
        && expansion.total_area_defect().to_bits() == zero
        && expansion.owner_edge_swept_areas()[diagonal].to_bits() == 1.5f64.to_bits();
    verdict(
        "flux-007-uniform-expansion",
        expansion_pass,
        &format!(
            "change={} boundary_sweep={} defect={:.3e}",
            expansion.total_area_change(),
            expansion.boundary_swept_area(),
            expansion.total_area_defect()
        ),
    );

    let square = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let sheared: Vec<_> = square
        .verts
        .iter()
        .map(|point| [point[0] + 0.5 * point[1], point[1]])
        .collect();
    let shear = audit_affine_triangle_gcl2(&square, &sheared, 1.0).unwrap();
    verdict(
        "flux-007-affine-gcl",
        close(shear.total_area_change(), 0.0, 1.0)
            && close(shear.boundary_swept_area(), 0.0, 1.0)
            && close(shear.total_area_defect(), 0.0, 1.0)
            && close(shear.worst_abs_cell_defect(), 0.0, 1.0),
        &format!(
            "shear change={:.3e} boundary_sweep={:.3e} defect={:.3e} worst_cell={:.3e}",
            shear.total_area_change(),
            shear.boundary_swept_area(),
            shear.total_area_defect(),
            shear.worst_abs_cell_defect()
        ),
    );
}

#[test]
fn flux_007_interior_cancellation_replay_and_translation_covariance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let mut next = mesh.verts.clone();
    next[3][0] += 0.25;
    let first = audit_affine_triangle_gcl2(&mesh, &next, 0.125).unwrap();
    let replay = audit_affine_triangle_gcl2(&mesh, &next, 0.125).unwrap();
    let interior = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();

    let translated_start = Mesh2 {
        nodes: mesh
            .verts
            .iter()
            .map(|point| [point[0] + 8.0, point[1] - 4.0])
            .collect(),
        elems: mesh.tris.iter().map(|tri| tri.to_vec()).collect(),
        patches: Vec::new(),
    };
    let translated_mesh = TriMesh::from_mesh2(&translated_start);
    let translated_next: Vec<_> = next
        .iter()
        .map(|point| [point[0] + 8.0, point[1] - 4.0])
        .collect();
    let shifted = audit_affine_triangle_gcl2(&translated_mesh, &translated_next, 0.125).unwrap();
    let exact_interior_balance = first.owner_edge_swept_areas()[interior].to_bits()
        == (-0.125f64).to_bits()
        && first.cells()[0].outward_swept_area().to_bits() == 0.0f64.to_bits()
        && first.cells()[1].outward_swept_area().to_bits() == 0.125f64.to_bits()
        && first.total_area_change().to_bits() == 0.125f64.to_bits()
        && first.boundary_swept_area().to_bits() == 0.125f64.to_bits();
    verdict(
        "flux-007-interior-cancellation",
        first == replay
            && first == shifted
            && first.interior_edges() == 1
            && first.boundary_edges() == 4
            && exact_interior_balance
            && first.worst_abs_interior_cancellation().to_bits() == 0.0f64.to_bits()
            && close(first.total_area_defect(), 0.0, 1.0)
            && close(first.worst_abs_cell_defect(), 0.0, 1.0),
        &format!(
            "interior_edges={} owner_sweep={} worst_pair_defect={:.3e}",
            first.interior_edges(),
            first.owner_edge_swept_areas()[interior],
            first.worst_abs_interior_cancellation()
        ),
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One ordered refusal matrix shares the same canonical mesh.
fn flux_007_refuses_bad_time_geometry_and_topology() {
    let mesh = single_triangle();
    let empty = TriMesh {
        verts: Vec::new(),
        tris: Vec::new(),
        edges: Vec::new(),
        tri_edges: Vec::new(),
        areas: Vec::new(),
        centroids: Vec::new(),
    };
    let empty_refused = matches!(
        audit_affine_triangle_gcl2(&empty, &[], 1.0),
        Err(AleGclError2::EmptyMesh)
    );
    let zero_time_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &mesh.verts, 0.0),
        Err(AleGclError2::InvalidTimeStep { .. })
    );
    let nan_time_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &mesh.verts, f64::NAN),
        Err(AleGclError2::InvalidTimeStep { .. })
    );
    let count_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &mesh.verts[..2], 1.0),
        Err(AleGclError2::VertexCountMismatch { .. })
    );

    let mut nonfinite = mesh.verts.clone();
    nonfinite[1][0] = f64::NAN;
    let nonfinite_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &nonfinite, 1.0),
        Err(AleGclError2::NonFiniteCoordinate {
            state: AleMeshState2::Final,
            vertex: 1,
            axis: 0
        })
    );

    let inverted = vec![[0.0, 0.0], [0.0, 1.0], [1.0, 0.0]];
    let inverted_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &inverted, 1.0),
        Err(AleGclError2::DegenerateCell {
            state: AleMeshState2::Final,
            cell: 0,
            ..
        })
    );

    let midstep_collapse = vec![[0.0, 0.0], [-1.0, 0.0], [0.0, -1.0]];
    let symmetric_collapse_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &midstep_collapse, 1.0),
        Err(AleGclError2::TrajectoryCollapse {
            cell: 0,
            normalized_time,
            signed_area
        }) if normalized_time.to_bits() == 0.5f64.to_bits()
            && signed_area.to_bits() == 0.0f64.to_bits()
    );

    let off_center_collapse = vec![[0.0, 0.0], [-1.0, 0.0], [0.0, -2.0]];
    let off_center_collapse_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &off_center_collapse, 1.0),
        Err(AleGclError2::TrajectoryCollapse {
            cell: 0,
            normalized_time,
            signed_area
        }) if close(normalized_time, 5.0 / 12.0, 1.0)
            && close(signed_area, -1.0 / 48.0, 1.0)
    );

    let quarter_turn = vec![[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]];
    let rotation_admitted = matches!(
        audit_affine_triangle_gcl2(&mesh, &quarter_turn, 1.0),
        Ok(receipt) if receipt.cells()[0].minimum_area().to_bits() == 0.25f64.to_bits()
            && receipt.total_area_change().to_bits() == 0.0f64.to_bits()
    );
    verdict(
        "flux-007-trajectory-guard",
        inverted_refused
            && symmetric_collapse_refused
            && off_center_collapse_refused
            && rotation_admitted,
        &format!(
            "endpoint={} symmetric={} off_center={} rotation={}",
            inverted_refused,
            symmetric_collapse_refused,
            off_center_collapse_refused,
            rotation_admitted
        ),
    );

    let mut malformed = single_triangle();
    malformed.tri_edges[0][0].1 = 0.0;
    let sign_refused = matches!(
        audit_affine_triangle_gcl2(&malformed, &malformed.verts, 1.0),
        Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::LocalEdgeSign,
            index: 0
        })
    );

    let mut stale_area = single_triangle();
    stale_area.areas[0] = 0.75;
    let stale_area_refused = matches!(
        audit_affine_triangle_gcl2(&stale_area, &stale_area.verts, 1.0),
        Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::StoredArea,
            index: 0
        })
    );

    let mut stale_centroid = single_triangle();
    stale_centroid.centroids[0][0] = f64::INFINITY;
    let stale_centroid_refused = matches!(
        audit_affine_triangle_gcl2(&stale_centroid, &stale_centroid.verts, 1.0),
        Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::StoredCentroid,
            index: 0
        })
    );

    let base = 9.0e307_f64;
    let adjacent = f64::from_bits(base.to_bits() + 1);
    let narrow_height = 1.0 / (adjacent - base);
    let overflowing_centroid = TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![[base, 0.0], [adjacent, 0.0], [base, narrow_height]],
        elems: vec![vec![0, 1, 2]],
        patches: Vec::new(),
    });
    let centroid_overflow_refused = matches!(
        audit_affine_triangle_gcl2(&overflowing_centroid, &overflowing_centroid.verts, 1.0),
        Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::CellCentroid,
            index: 0
        })
    );

    let mut stale_edge = single_triangle();
    stale_edge.edges[0].len *= 2.0;
    let stale_edge_refused = matches!(
        audit_affine_triangle_gcl2(&stale_edge, &stale_edge.verts, 1.0),
        Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::StoredEdgeGeometry,
            index: 0
        })
    );

    let same_direction = TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        elems: vec![vec![0, 1, 2], vec![0, 1, 3]],
        patches: Vec::new(),
    });
    let orientation_refused = matches!(
        audit_affine_triangle_gcl2(&same_direction, &same_direction.verts, 1.0),
        Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::InteriorEdgeOrientation,
            ..
        })
    );

    let overflow = vec![[0.0, 0.0], [f64::MAX, f64::MAX], [-f64::MAX, f64::MAX]];
    let overflow_refused = matches!(
        audit_affine_triangle_gcl2(&mesh, &overflow, 1.0),
        Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::CellArea,
            index: 0
        })
    );
    verdict(
        "flux-007-fail-closed",
        empty_refused
            && zero_time_refused
            && nan_time_refused
            && count_refused
            && nonfinite_refused
            && sign_refused
            && stale_area_refused
            && stale_centroid_refused
            && centroid_overflow_refused
            && stale_edge_refused
            && orientation_refused
            && overflow_refused,
        &format!(
            "empty={} time={} nan_time={} count={} finite={} sign={} area_cache={} centroid_cache={} centroid_overflow={} edge_cache={} orientation={} overflow={}",
            empty_refused,
            zero_time_refused,
            nan_time_refused,
            count_refused,
            nonfinite_refused,
            sign_refused,
            stale_area_refused,
            stale_centroid_refused,
            centroid_overflow_refused,
            stale_edge_refused,
            orientation_refused,
            overflow_refused
        ),
    );
}

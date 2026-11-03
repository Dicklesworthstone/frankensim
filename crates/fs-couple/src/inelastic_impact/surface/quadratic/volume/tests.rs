use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solid::linear3::{TetAssemblyBudget, tet10::promote_p1_to_p2};

const CORNERS: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];

fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 23,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn controls() -> QuadraticTetPathControls {
    QuadraticTetPathControls {
        maximum_cells: 65_536,
        maximum_active_cells: 4_096,
        fraction_tolerance: 1e-3,
        parameter_tolerance: 1e-3,
        minimum_scaled_determinant: 1e-12,
    }
}

#[test]
fn complete_material_paths_detect_a_fold_hidden_by_unchanged_corners() {
    run(|_, cx| {
        let (nodes, tets) = promote_p1_to_p2(&CORNERS[..4], &[[0, 1, 2, 3]], None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &tets,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        let prepared = surface.prepare_volume_paths(controls(), cx).unwrap();
        assert!(matches!(
            prepared.reference_admission().verdict,
            QuadraticSolidVolumePathVerdict::Positive { .. }
        ));
        // X + 1.5*t*X*(1-X) fixes both endpoints of every x edge. Its
        // determinant at X=1 is 1-1.5*t, reaching zero at t=2/3.
        let folded: Vec<_> = nodes
            .iter()
            .map(|p| [p[0] + 1.5 * p[0] * (1.0 - p[0]), p[1], p[2]])
            .collect();
        assert_eq!(&nodes[..4], &folded[..4]);
        let report = prepared
            .affine_path_at_positions(&nodes, &folded, controls(), cx)
            .unwrap();
        match report.verdict {
            QuadraticSolidVolumePathVerdict::Unresolved {
                tetrahedron,
                fraction_start,
                fraction_end,
                barycentric_bounds,
                ..
            } => {
                assert_eq!(tetrahedron, 0);
                assert!(fraction_start <= 2.0 / 3.0);
                assert!(fraction_end - fraction_start <= controls().fraction_tolerance);
                assert!(barycentric_bounds[1][1] > 0.9);
            }
            other => panic!("corner-only orientation hid a material fold: {other:?}"),
        }
        assert_eq!(report.tested_tetrahedra, 1);
        assert!(report.examined_cells <= controls().maximum_cells);
    });
}

#[test]
fn reference_admission_cannot_be_replaced_by_positive_face_and_centroid_samples() {
    run(|_, cx| {
        let (mut nodes, tets) = promote_p1_to_p2(&CORNERS[..4], &[[0, 1, 2, 3]], None, cx).unwrap();
        // The native face-centroid and volume-centroid probes stay positive,
        // but the complete reference map has determinant -1/8 at X=1.
        for p in &mut nodes {
            p[0] += 1.125 * p[0] * (1.0 - p[0]);
        }
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &tets,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        assert!(matches!(
            surface.prepare_volume_paths(controls(), cx),
            Err(QuadraticSolidVolumeError::ReferenceUnresolved { .. })
        ));
    });
}

#[test]
fn curved_isochoric_paths_reuse_reference_admission_and_share_mesh_work_limits() {
    run(|gate, cx| {
        let (nodes, tets) =
            promote_p1_to_p2(&CORNERS, &[[0, 1, 2, 3], [0, 2, 1, 4]], None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &tets,
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let prepared = surface.prepare_volume_paths(controls(), cx).unwrap();
        // det F=1 over BOTH complete tetrahedra and the complete time path.
        let end: Vec<_> = nodes
            .iter()
            .map(|p| [p[0] - 8.0 * p[1] * p[2] + 1.25, p[1] - 2.0, p[2] + 3.0])
            .collect();
        let report = prepared
            .affine_path_at_positions(&nodes, &end, controls(), cx)
            .unwrap();
        assert_eq!(report.tested_tetrahedra, 2);
        assert!(matches!(report.verdict,
            QuadraticSolidVolumePathVerdict::Positive { minimum_scaled_determinant }
                if minimum_scaled_determinant > controls().minimum_scaled_determinant));
        let limited = QuadraticTetPathControls {
            maximum_cells: 1,
            maximum_active_cells: 1,
            ..controls()
        };
        assert!(
            prepared
                .affine_path_at_positions(&nodes, &end, limited, cx)
                .is_err()
        );
        assert!(surface.prepare_volume_paths(limited, cx).is_err());
        gate.request();
        assert!(matches!(
            prepared.affine_path_at_positions(&nodes, &end, controls(), cx),
            Err(QuadraticSolidVolumeError::Cancelled)
        ));
    });
}

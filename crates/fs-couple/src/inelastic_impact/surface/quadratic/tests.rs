use super::*;
use crate::inelastic_impact::multiple::FiniteSolidContactCounterpart;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion};
use fs_solid::linear3::{TetAssemblyBudget, tet10::promote_p1_to_p2};

const CORNERS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn controls() -> QuadraticSolidSurfaceControls {
    QuadraticSolidSurfaceControls {
        surface: FiniteSolidSurfaceControls {
            activation_tolerance_m: 0.02,
            maximum_contacts: 8,
            ..Default::default()
        },
        query: QuadraticQueryControls {
            distance_tolerance: 1e-9,
            witness_tolerance: 1e-6,
            ..Default::default()
        },
        maximum_witness_span_m: 3e-6,
    }
}

fn sphere(id: u64, rigid_body: usize) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    }
}

fn body(center: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(center, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap()
}

#[test]
fn conforming_internal_faces_cancel_and_every_original_midside_binds_identity() {
    with_cx(|cx| {
        let mut corners = CORNERS.to_vec();
        corners.push([0.0, 0.0, -1.0]);
        let (nodes, elements) =
            promote_p1_to_p2(&corners, &[[0, 1, 2, 3], [0, 2, 1, 4]], None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &elements,
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        assert_eq!(surface.faces().len(), 6);
        for face in surface.faces() {
            let mut corners = [face[0], face[1], face[2]];
            corners.sort_unstable();
            assert_ne!(corners, [0, 1, 2]);
        }
        assert!(surface.matches_reference(&nodes, &elements));
        let mut changed = nodes.clone();
        changed[elements[0][4]][2] += 0.03125;
        assert!(!surface.matches_reference(&changed, &elements));
        assert!(matches!(
            surface.current_positions_from_reference(
                &changed,
                &elements,
                &vec![[0.0; 3]; changed.len()],
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::ReferenceMismatch
            ))
        ));
        let mut reordered = elements.clone();
        reordered.swap(0, 1);
        assert!(!surface.matches_reference(&nodes, &reordered));

        // Identical coordinates do not make two DIFFERENT shared-edge node IDs
        // conforming. Native topology refuses before cancelling the shared face.
        let mut nonconforming_nodes = nodes.clone();
        let mut nonconforming = elements.clone();
        nonconforming_nodes.push(nodes[elements[1][4]]);
        nonconforming[1][4] = nonconforming_nodes.len() - 1;
        assert!(matches!(
            QuadraticSolidContactSurface::from_reference(
                &nonconforming_nodes,
                &nonconforming,
                TetAssemblyBudget::standard(),
                8,
                cx,
            ),
            Err(QuadraticSolidSurfaceError::SolidGeometry(
                TetElasticError::MidsideMismatch { .. }
            ))
        ));
    });
}

#[test]
fn curved_material_edges_and_vertices_have_stable_deduplicated_contact_ids() {
    with_cx(|cx| {
        let (nodes, elements) = promote_p1_to_p2(&CORNERS, &[[0, 1, 2, 3]], None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &elements,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        let bodies = [
            body(Vec3::new(-0.08, -0.08, -0.08)),
            body(Vec3::new(0.25, -0.1, -0.1)),
        ];
        let spheres = [sphere(9, 0), sphere(3, 1)];
        let query = surface
            .sphere_separations_at_positions(&nodes, &bodies, &spheres, controls(), cx)
            .unwrap();
        assert_eq!(query.contacts.len(), 2);
        assert!(query.candidate_triangles >= 5);
        assert_eq!(
            query
                .closest
                .iter()
                .map(|w| w.sphere_id)
                .collect::<Vec<_>>(),
            [3, 9]
        );
        let contacts = query.admitted_contacts().unwrap();
        for witness in &query.contacts {
            let FiniteSolidContactCounterpart::Solid(point) = witness.contact.counterpart else {
                panic!("quadratic exterior must retain the original solid point");
            };
            assert_eq!(point.tetrahedron, 0);
            let expected = if witness.sphere_id == 3 {
                [0.75, 0.25, 0.0, 0.0]
            } else {
                [1.0, 0.0, 0.0, 0.0]
            };
            for (actual, expected) in point.barycentric.into_iter().zip(expected) {
                assert!((actual - expected).abs() < 2e-6);
            }
            assert!(witness.gap_lower_m > 0.0);
            assert_ne!(witness.rigid_point_world_m, witness.solid_point_world_m);
        }
        let reordered = surface
            .sphere_separations_at_positions(
                &nodes,
                &bodies,
                &[spheres[1], spheres[0]],
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(contacts, reordered.admitted_contacts().unwrap());
        assert!(matches!(
            surface.sphere_separations_at_positions(
                &nodes,
                &bodies,
                &spheres,
                QuadraticSolidSurfaceControls {
                    surface: FiniteSolidSurfaceControls {
                        maximum_contacts: 1,
                        ..controls().surface
                    },
                    ..controls()
                },
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::ContactBudgetExhausted { .. }
            ))
        ));
    });
}

#[test]
fn interior_and_uncertain_activation_refuse_without_repairing_raw_geometry() {
    with_cx(|cx| {
        let (nodes, elements) = promote_p1_to_p2(&CORNERS, &[[0, 1, 2, 3]], None, cx).unwrap();
        let surface = QuadraticSolidContactSurface::from_reference(
            &nodes,
            &elements,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        assert!(matches!(
            surface.sphere_separations_at_positions(
                &nodes,
                &[body(Vec3::new(0.1, 0.1, 0.1))],
                &[sphere(1, 0)],
                controls(),
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::UnresolvedExterior { .. }
            ))
        ));
        let mut admission = controls();
        // Exact binary gap equals the band boundary. Directed numerical bounds
        // straddle it; raw detection may return, but admission cannot omit it.
        admission.surface.activation_tolerance_m = 0.03125;
        let raw = surface
            .sphere_separations_at_positions(
                &nodes,
                &[body(Vec3::new(0.25, 0.25, -0.15625))],
                &[sphere(1, 0)],
                admission,
                cx,
            )
            .unwrap();
        assert!(raw.closest[0].gap_m > 0.0);
        assert!(raw.unresolved_activation.is_some());
        assert!(matches!(
            raw.admitted_contacts(),
            Err(QuadraticSolidSurfaceError::UnresolvedActivation { .. })
        ));
    });
}

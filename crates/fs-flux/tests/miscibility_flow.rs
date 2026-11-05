//! G1/G3/G4: native directional flow carries the original conserved mixture.

use fs_flux::miscibility::advection::BinaryMiscibilityAdvectionControls;
use fs_flux::miscibility::flow::BinaryMiscibilityFlowError;
use fs_flux::miscibility::{
    BinaryMiscibilityCell, BinaryMiscibilityError, BinaryMiscibilityFace, BinaryMiscibilityNetwork,
};
use fs_flux::{FluxParams, FluxSystem, StokesBoundary, TriMesh};
use fs_qty::{AmountConcentration, Temperature};
use fs_solid::Mesh2;
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::BinaryRegularSolutionV1;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;

fn geometry(
    mesh: &TriMesh,
    depth: f64,
) -> (Vec<BinaryMiscibilityCell>, Vec<BinaryMiscibilityFace>) {
    let cells = mesh
        .areas
        .iter()
        .map(|&area| BinaryMiscibilityCell {
            volume_m3: area * depth,
        })
        .collect();
    let faces = mesh
        .edges
        .iter()
        .filter(|edge| edge.tris.1 != usize::MAX)
        .map(|edge| BinaryMiscibilityFace {
            cells: [edge.tris.0, edge.tris.1],
            area_m2: edge.len * depth,
            // Transport-only fixture: no diffusion geometry is inferred from
            // these generally nonorthogonal triangle-centroid connections.
            distance_m: 1.0,
            diffusivity_m2_s: 0.0,
        })
        .collect();
    (cells, faces)
}

fn network(
    cells: Vec<BinaryMiscibilityCell>,
    faces: Vec<BinaryMiscibilityFace>,
) -> BinaryMiscibilityNetwork {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        0.0,
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1.0),
        0.0,
    )
    .unwrap();
    BinaryMiscibilityNetwork::new(model, cells, faces, 16 * 1024 * 1024, || false).unwrap()
}

fn native_vortex(system: &FluxSystem<'_>) -> Vec<f64> {
    let mut velocity = vec![0.0; system.n_u];
    // curl(psi), psi=x(1-x)y(1-y), is closed and divergence free. The
    // three-point rule integrates its exact BDM mean and linear moments.
    for (index, edge) in system.mesh.edges.iter().enumerate() {
        let a = system.mesh.verts[edge.verts.0];
        let b = system.mesh.verts[edge.verts.1];
        let offset = 15.0_f64.sqrt() / 10.0;
        for (s, weight) in [
            (0.5 - offset, 5.0 / 18.0),
            (0.5, 4.0 / 9.0),
            (0.5 + offset, 5.0 / 18.0),
        ] {
            let x = (1.0 - s) * a[0] + s * b[0];
            let y = (1.0 - s) * a[1] + s * b[1];
            let u = [
                x * (1.0 - x) * (1.0 - 2.0 * y),
                -(1.0 - 2.0 * x) * y * (1.0 - y),
            ];
            let normal = u[0] * edge.normal[0] + u[1] * edge.normal[1];
            velocity[2 * index] += weight * normal;
            velocity[2 * index + 1] += weight * normal * (s - 0.5);
        }
    }
    velocity
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-13,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_zero_net_native_counterflow_mixes_and_respects_depth() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut velocity = vec![0.0; system.n_u];
    velocity[2 * edge + 1] = 0.2;
    for depth in [0.2, 2.0] {
        let (cells, faces) = geometry(&mesh, depth);
        let network = network(cells, faces);
        let initial = network
            .state_from_mole_fractions(vec![[1.0, 0.0], [0.0, 1.0]], || false)
            .unwrap();
        let flow = network
            .prepare_advection_from_bdm(&system, &velocity, depth, 0.0, || false)
            .unwrap();
        // u.n=2.4(s-1/2); each triangular lobe integrates to 0.3.
        let rate = 0.3 * mesh.edges[edge].len * depth;
        near(flow.face_rates()[0].forward_m3_s, rate);
        near(flow.face_rates()[0].reverse_m3_s, rate);
        let next = network
            .advect_candidate(
                &initial,
                &flow,
                0.1,
                BinaryMiscibilityAdvectionControls::default(),
                || false,
            )
            .unwrap();
        assert_eq!(next.step().substeps, 1);
        let moved = 0.1 * rate;
        near(next.state().component_moles()[0][0], 0.5 * depth - moved);
        near(next.state().component_moles()[0][1], moved);
        near(next.state().component_moles()[1][0], moved);
        near(next.step().transferred_component_moles[0][0], moved);
        near(next.step().transferred_component_moles[0][1], -moved);
        near(next.state().mole_fractions()[1][0], moved / (0.5 * depth));
        assert_eq!(initial.mole_fractions(), &[[1.0, 0.0], [0.0, 1.0]]);
    }
}

#[test]
fn g3_mesh_orientation_and_input_face_order_do_not_change_transport() {
    let mut mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let depth = 0.4;
    let (cells, mut faces) = geometry(&mesh, depth);
    faces.reverse();
    for face in &mut faces {
        face.cells.swap(0, 1);
    }
    let network = network(cells, faces);
    let initial = network
        .state_from_mole_fractions(
            mesh.centroids.iter().map(|p| [p[0], 1.0 - p[0]]).collect(),
            || false,
        )
        .unwrap();
    let system = FluxSystem::new(&mesh);
    let velocity = native_vortex(&system);
    let flow = network
        .prepare_advection_from_bdm(&system, &velocity, depth, 1e-12, || false)
        .unwrap();
    let expected = network
        .advect_candidate(
            &initial,
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    assert_ne!(
        expected.state().component_moles(),
        initial.component_moles()
    );
    assert!(
        expected
            .step()
            .component_residual_mol
            .iter()
            .all(|x| x.abs() < 1e-13)
    );
    // Reverse the stored global orientation without changing physical flow.
    for edge in &mut mesh.edges {
        if edge.tris.1 != usize::MAX {
            edge.tris = (edge.tris.1, edge.tris.0);
            edge.normal = edge.normal.map(|n| -n);
        }
    }
    for cell in &mut mesh.tri_edges {
        for (edge, sign) in cell {
            if mesh.edges[*edge].tris.1 != usize::MAX {
                *sign = -*sign;
            }
        }
    }
    let system = FluxSystem::new(&mesh);
    let reversed_velocity = native_vortex(&system);
    let reversed = network
        .prepare_advection_from_bdm(&system, &reversed_velocity, depth, 1e-12, || false)
        .unwrap();
    assert_eq!(flow.face_rates(), reversed.face_rates());
    let actual = network
        .advect_candidate(
            &initial,
            &reversed,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn g1_solved_closed_stokes_field_stirs_the_conserved_mixture() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundary: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let force = |p: [f64; 2]| [p[1] * (1.0 - p[1]), -p[0] * (1.0 - p[0])];
    let solved = system
        .assemble_stokes_mixed(
            FluxParams {
                tol: 2e-13,
                ..FluxParams::default()
            },
            &force,
            &boundary,
        )
        .unwrap()
        .solve()
        .unwrap();
    let (cells, faces) = geometry(&mesh, 0.5);
    let network = network(cells, faces);
    let initial = network
        .state_from_mole_fractions(
            mesh.centroids.iter().map(|p| [p[0], 1.0 - p[0]]).collect(),
            || false,
        )
        .unwrap();
    let flow = network
        .prepare_advection_from_bdm(&system, &solved.field.x[..system.n_u], 0.5, 1e-12, || false)
        .unwrap();
    let next = network
        .advect_candidate(
            &initial,
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    assert!(flow.maximum_outgoing_rate_s() > 0.0);
    assert_ne!(next.state().component_moles(), initial.component_moles());
    assert!(
        next.step()
            .component_residual_mol
            .iter()
            .all(|n| n.abs() < 1e-13)
    );
    assert!(next.step().maximum_capacity_residual_mol < 1e-13);
}

#[test]
fn g4_zero_mean_exterior_counterflow_and_hidden_mean_divergence_refuse() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let (cells, faces) = geometry(&mesh, 1.0);
    let network = network(cells, faces);
    let mut velocity = vec![0.0; system.n_u];
    let boundary = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 == usize::MAX)
        .unwrap();
    velocity[2 * boundary + 1] = 1.0;
    assert!(
        matches!(network.prepare_advection_from_bdm(&system, &velocity, 1.0, 1e-12, || false), Err(BinaryMiscibilityFlowError::OpenBoundary { edge }) if edge == boundary)
    );
    velocity.fill(0.0);
    let interior = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    velocity[2 * interior] = 1e-30;
    velocity[2 * interior + 1] = 1.0;
    let native = system.edge_volume_transport(&velocity).unwrap();
    assert_eq!(
        native[interior].along_normal_m2_s(),
        native[interior].against_normal_m2_s()
    );
    assert!(matches!(
        network.prepare_advection_from_bdm(&system, &velocity, 1.0, 0.0, || false),
        Err(BinaryMiscibilityFlowError::VolumeDivergence { .. })
    ));
}

#[test]
fn g4_mismatched_graph_geometry_and_native_input_refuse() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = vec![0.0; system.n_u];
    for defect in 0..4 {
        let (mut cells, mut faces) = geometry(&mesh, 1.0);
        match defect {
            0 => cells[0].volume_m3 *= 1.01,
            1 => faces[0].area_m2 *= 1.01,
            2 => {
                faces.pop();
            }
            _ => faces.push(BinaryMiscibilityFace {
                cells: [0, cells.len() - 1],
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_m2_s: 0.0,
            }),
        }
        let network = network(cells, faces);
        assert!(matches!(
            network.prepare_advection_from_bdm(&system, &velocity, 1.0, 1e-12, || false),
            Err(BinaryMiscibilityFlowError::Geometry(_))
        ));
    }
    let (cells, faces) = geometry(&mesh, 1.0);
    let network = network(cells, faces);
    assert!(matches!(
        network.prepare_advection_from_bdm(
            &system,
            &velocity[..velocity.len() - 1],
            1.0,
            1e-12,
            || false
        ),
        Err(BinaryMiscibilityFlowError::Transport(_))
    ));
    for depth in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            network
                .prepare_advection_from_bdm(&system, &velocity, depth, 1e-12, || false)
                .is_err()
        );
    }
}

#[test]
fn g4_cancellation_never_publishes_a_partial_prepared_flow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = native_vortex(&system);
    let (cells, faces) = geometry(&mesh, 1.0);
    let network = network(cells, faces);
    let mut calls = 0;
    network
        .prepare_advection_from_bdm(&system, &velocity, 1.0, 1e-12, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, mesh.tris.len() + 3, calls] {
        let mut count = 0;
        assert!(matches!(
            network.prepare_advection_from_bdm(&system, &velocity, 1.0, 1e-12, || {
                count += 1;
                count == stop
            }),
            Err(BinaryMiscibilityFlowError::Miscibility(
                BinaryMiscibilityError::Cancelled
            ))
        ));
    }
}

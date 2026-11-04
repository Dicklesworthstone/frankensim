#![cfg(feature = "continuum")]
//! Physical directional rates, including opposing flow on a single BDM1 face.

use fs_flux::ns::transport::FluxTransportError;
use fs_flux::{FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 32.0 * f64::EPSILON * scale,
        "{actual:e} != {expected:e} (scale {scale:e})"
    );
}

#[test]
fn zero_net_face_retains_both_physical_streams() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let length = mesh.edges[edge].len;
    let mut velocity = vec![0.0; system.n_u];
    for moment in [0.5, -0.5] {
        velocity[2 * edge + 1] = moment;
        let transport = system.edge_volume_transport(&velocity).unwrap();
        // u.n = +/- 6(s-1/2): each triangular lobe has area 3/4.
        near(transport[edge].along_normal_m2_s(), 0.75 * length, length);
        near(transport[edge].against_normal_m2_s(), 0.75 * length, length);
        assert_eq!(transport[edge].net_m2_s(), 0.0);
        for (index, flow) in transport.iter().enumerate() {
            if index != edge {
                assert_eq!(flow.along_normal_m2_s(), 0.0);
                assert_eq!(flow.against_normal_m2_s(), 0.0);
            }
        }
    }
}

#[test]
fn shifted_counterflow_and_constant_reversal_match_exact_face_areas() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(3.0, 2.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    for (mean, moment, positive, negative) in [
        (1.0, 0.5, 4.0 / 3.0, 1.0 / 3.0),
        (-1.0, -0.5, 1.0 / 3.0, 4.0 / 3.0),
        (2.5, 0.0, 2.5, 0.0),
        (-2.5, 0.0, 0.0, 2.5),
        (3.0, 0.5, 3.0, 0.0),
    ] {
        let velocity: Vec<_> = (0..mesh.edges.len()).flat_map(|_| [mean, moment]).collect();
        let transport = system.edge_volume_transport(&velocity).unwrap();
        for (edge, flow) in mesh.edges.iter().zip(&transport) {
            near(
                flow.along_normal_m2_s(),
                positive * edge.len,
                3.0 * edge.len,
            );
            near(
                flow.against_normal_m2_s(),
                negative * edge.len,
                3.0 * edge.len,
            );
            near(flow.net_m2_s(), mean * edge.len, 3.0 * edge.len);
            near(
                flow.along_normal_m2_s() - flow.against_normal_m2_s(),
                flow.net_m2_s(),
                3.0 * edge.len,
            );
        }
    }
}

#[test]
fn closed_recirculation_preserves_cellwise_continuity() {
    for n in [2, 4, 8] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let mut velocity = vec![0.0; system.n_u];
        // Independent divergence-free polynomial from psi=x(1-x)y(1-y),
        // with zero outward flux on the entire square boundary. Three-point
        // Gauss exactly projects its cubic normal trace and linear moment.
        for (index, edge) in mesh.edges.iter().enumerate() {
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
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
        let transport = system.edge_volume_transport(&velocity).unwrap();
        for face in &mesh.tri_edges {
            let mut net = 0.0;
            let mut incoming = 0.0;
            let mut outgoing = 0.0;
            for &(edge, sign) in face {
                let flow = transport[edge];
                net += sign * flow.net_m2_s();
                let [out, into] = if sign > 0.0 {
                    [flow.along_normal_m2_s(), flow.against_normal_m2_s()]
                } else {
                    [flow.against_normal_m2_s(), flow.along_normal_m2_s()]
                };
                outgoing += out;
                incoming += into;
            }
            near(net, 0.0, incoming.max(outgoing));
            near(incoming, outgoing, incoming.max(outgoing));
        }
        for (index, flow) in transport.iter().enumerate() {
            if mesh.is_boundary(index) {
                assert_eq!(flow.along_normal_m2_s(), 0.0);
                assert_eq!(flow.against_normal_m2_s(), 0.0);
            }
        }
    }
}

#[test]
fn original_net_moment_survives_large_counterflow_and_invalid_fields_refuse() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let mut velocity = vec![0.0; system.n_u];
    velocity[0] = 1e-30;
    velocity[1] = 1.0;
    let flow = system.edge_volume_transport(&velocity).unwrap()[0];
    assert_eq!(flow.net_m2_s(), 1e-30 * mesh.edges[0].len);
    assert!(flow.along_normal_m2_s() > 0.0 && flow.against_normal_m2_s() > 0.0);
    assert!(matches!(
        system.edge_volume_transport(&velocity[..system.n_u - 1]),
        Err(FluxTransportError::VelocityLength { .. })
    ));
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        velocity[1] = invalid;
        assert_eq!(
            system.edge_volume_transport(&velocity),
            Err(FluxTransportError::NonFiniteVelocity { dof: 1 })
        );
    }
    velocity.fill(0.0);
    let long_edge = mesh.edges.iter().position(|edge| edge.len > 1.0).unwrap();
    velocity[2 * long_edge] = f64::MAX;
    assert_eq!(
        system.edge_volume_transport(&velocity),
        Err(FluxTransportError::UnrepresentableRate { edge: long_edge })
    );
}

#[test]
fn nonaffine_bdf1_endpoints_preserve_local_continuity_near_stagnation() {
    use fs_flux::{
        MaterialFluxNumerics, MaterialFluxSystem, MaterialNavierStokesBoundary,
        NavierStokesControls,
    };
    use fs_material::fluid_transport::FluidTransportState;

    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let reference = FluidTransportState::try_new(2.0, 5.0, 0.2, 1.0).unwrap();
    let fluid = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 2e-13,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let walls = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &walls,
        })
        .collect();
    let force = |p: [f64; 2]| [p[1] * (1.0 - p[1]), -p[0] * (1.0 - p[0])];
    let mut previous = vec![0.0; system.n_u];
    for step in 1..=2 {
        let next = fluid
            .bdf1_navier_stokes_mixed_rheology_step(
                &force,
                &boundaries,
                &previous,
                0.01,
                NavierStokesControls {
                    maximum_iterations: 40,
                    relative_tolerance: 2e-11,
                    ..NavierStokesControls::default()
                },
                |_, reference| Ok(*reference),
            )
            .unwrap();
        let velocity = next.step.solution.velocity_dofs();
        let faces = system.edge_volume_transport(velocity).unwrap();
        for (edge, face) in faces.iter().enumerate() {
            if mesh.is_boundary(edge) {
                // Exact prescribed zero moments must not create ghost inflow.
                assert_eq!(velocity[2 * edge], 0.0);
                assert_eq!(velocity[2 * edge + 1], 0.0);
                assert_eq!(face.along_normal_m2_s(), 0.0);
                assert_eq!(face.against_normal_m2_s(), 0.0);
            }
        }
        assert!(
            next.step
                .solution
                .viscous_dissipation()
                .total_power_w_per_m()
                > 0.0
        );
        for (cell, edges) in mesh.tri_edges.iter().enumerate() {
            let mut incoming = 0.0_f64;
            let mut outgoing = 0.0_f64;
            for &(edge, sign) in edges {
                let face = faces[edge];
                let [out, into] = if sign > 0.0 {
                    [face.along_normal_m2_s(), face.against_normal_m2_s()]
                } else {
                    [face.against_normal_m2_s(), face.along_normal_m2_s()]
                };
                outgoing += out;
                incoming += into;
            }
            // The fixed-inventory transport bound is relative to THIS cell's
            // actual streams, so a global residual cannot hide a stagnant cell.
            assert!(
                (incoming - outgoing).abs() <= 64.0 * f64::EPSILON * incoming.max(outgoing),
                "step {step}, cell {cell}: incoming={incoming:e}, outgoing={outgoing:e}"
            );
        }
        previous = velocity.to_vec();
    }
}

use super::*;
use crate::ns::transport::scalar::{ScalarBoundary, ScalarTransportControls};
use crate::{FluxSystem, TriMesh};

#[test]
fn inactive_limiting_recovers_original_consistent_mass_and_complete_sip_reservoir() {
    let mut geometry = fs_solid::Mesh2::triangles(1.0, 1.0, 2, 2);
    // Obtuse triangles expose the signed weak SIP load that a right-triangle
    // fixture can miss despite testing the same boundary API.
    for point in &mut geometry.nodes {
        point[0] += 2.0 * point[1];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    let boundary = |p: [f64; 2]| 300.0 + p[0] - 0.2 * p[1];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| ScalarBoundary::Value {
            edge,
            value: &boundary,
        })
        .collect();
    let native = system
        .prepare_scalar_transport(
            &vec![0.0; system.n_u],
            &(0..mesh.tris.len())
                .map(|i| 0.02 + 0.005 * i as f64)
                .collect::<Vec<_>>(),
            &boundaries,
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    assert!(
        native.boundary_fluxes.iter().any(|f| f
            .reservoir_coefficients
            .iter()
            .flatten()
            .any(|&b| b < 0.0))
    );
    let bounded = native
        .prepare_bounded_transport(
            ScalarBounds {
                lower: 0.0,
                upper: 1000.0,
            },
            || false,
        )
        .unwrap();
    let previous: Vec<_> = (0..3 * mesh.tris.len())
        .map(|i| 299.0 + 0.1 * i as f64)
        .collect();
    let dt = 1e-5;
    let euler = |old: &[f64]| {
        let mut derivative = vec![0.0; old.len()];
        native.operator.spmv(old, &mut derivative);
        for cell in 0..mesh.tris.len() {
            let matrix: Vec<_> = (0..3)
                .flat_map(|i| (0..3).map(move |j| if i == j { 2.0 } else { 1.0 }))
                .map(|v| mesh.areas[cell] / 12.0 * v)
                .collect();
            let mut rhs: Vec<_> = (0..3)
                .map(|i| native.boundary_load[3 * cell + i] - derivative[3 * cell + i])
                .collect();
            // Independent original dense LU oracle, not the production block
            // inverse or a lumped-mass target that could hide missing fluxes.
            fs_la::factor::lu(&matrix, 3).unwrap().solve(&mut rhs);
            for i in 0..3 {
                derivative[3 * cell + i] = old[3 * cell + i] + dt * rhs[i];
            }
        }
        derivative
    };
    let first = euler(&previous);
    let second = euler(&first);
    let result = bounded
        .ssprk2_step(&previous, dt, BoundedScalarControls::default(), || false)
        .unwrap();
    assert_eq!(result.substeps, 1);
    assert_eq!(result.limited_fluxes, 0);
    for i in 0..previous.len() {
        let expected = 0.5 * previous[i] + 0.5 * second[i];
        assert!((result.values[i] - expected).abs() < 3e-13);
    }
    assert!(result.boundaries.iter().all(|b| b.limiter_diffusion == 0.0));
    assert!(result.balance_residual.abs() < 3e-13);
    // The SAME signed boundary loads now heat an initially cold field at its
    // authored lower bound. Any unremoved reservoir graph flux must appear in
    // the accepted exchange; restoring the high-order endpoint by clipping
    // would lose this account.
    let bounded_heat = native
        .prepare_bounded_transport(
            ScalarBounds {
                lower: 299.0,
                upper: 303.0,
            },
            || false,
        )
        .unwrap();
    let heating = bounded_heat
        .ssprk2_step(
            &vec![299.0; previous.len()],
            0.02,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    assert!(heating.minimum >= 299.0 - 1e-11 && heating.maximum <= 303.0 + 1e-11);
    assert!(heating.integral_after > heating.integral_before);
    assert!(
        heating
            .boundaries
            .iter()
            .any(|b| b.limiter_diffusion.abs() > 1e-12)
    );
    assert!(heating.balance_residual.abs() < 2e-11);
}

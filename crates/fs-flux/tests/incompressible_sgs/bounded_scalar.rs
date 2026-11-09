//! Physical sharp-front transport, conservation, SGS coupling and refinement.

use super::*;
use fs_flux::ns::transport::scalar::bounded::{
    BoundedScalarControls, BoundedScalarError, ScalarBounds,
};
use fs_flux::ns::transport::scalar::{
    ScalarBoundary, ScalarSolveControls, ScalarTransportControls,
};

fn boundaries<'a>(
    mesh: &TriMesh,
    zero: &'a dyn Fn([f64; 2]) -> f64,
    donor: Option<&'a dyn Fn([f64; 2]) -> f64>,
) -> Vec<ScalarBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| ScalarBoundary::Flux {
            edge,
            outward: zero,
            inflow: donor,
        })
        .collect()
}
fn nodal(mesh: &TriMesh, f: &dyn Fn([f64; 2]) -> f64) -> Vec<f64> {
    mesh.tris
        .iter()
        .flat_map(|tri| tri.map(|v| f(mesh.verts[v])))
        .collect()
}
fn fraction_bounds() -> ScalarBounds {
    ScalarBounds {
        lower: 0.0,
        upper: 1.0,
    }
}

#[test]
fn sharp_high_peclet_front_is_bounded_without_the_inventory_error_of_clipping() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 8, 3));
    let system = FluxSystem::new(&mesh);
    let velocity = moments(&mesh, &|_| [1.0, 0.0]);
    let zero = |_: [f64; 2]| 0.0;
    let inlet = |_: [f64; 2]| 1.0;
    let native = system
        .prepare_scalar_transport(
            &velocity,
            &vec![1e-8; mesh.tris.len()],
            &boundaries(&mesh, &zero, Some(&inlet)),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let previous: Vec<_> = mesh
        .centroids
        .iter()
        .flat_map(|p| [if p[0] < 0.5 { 1.0 } else { 0.0 }; 3])
        .collect();
    let snapshot = previous.clone();
    let unrestricted = native
        .bdf1_step(
            &previous,
            0.005,
            &|_, _| 0.0,
            ScalarSolveControls {
                relative_tolerance: 1e-12,
                maximum_iterations: 1000,
            },
            || false,
        )
        .unwrap();
    let minimum = unrestricted
        .values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let maximum = unrestricted
        .values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        minimum < -1e-5 || maximum > 1.0 + 1e-5,
        "missing unrestricted negative control: {minimum:e}..{maximum:e}"
    );
    let clipping_inventory_change: f64 = unrestricted
        .values
        .iter()
        .enumerate()
        .map(|(i, &v)| mesh.areas[i / 3] / 3.0 * (v.clamp(0.0, 1.0) - v))
        .sum();
    assert!(clipping_inventory_change.abs() > 1e-8);
    let bounded = native
        .prepare_bounded_transport(fraction_bounds(), || false)
        .unwrap();
    let result = bounded
        .ssprk2_step(&previous, 0.15, BoundedScalarControls::default(), || false)
        .unwrap();
    assert!(result.substeps > 1 && result.limited_fluxes > 0);
    assert!(result.minimum >= -3e-14 && result.maximum <= 1.0 + 3e-14);
    close(result.balance_residual, 0.0, 3e-12);
    assert!(result.integral_after > result.integral_before + 0.1);
    assert_eq!(previous, snapshot);
}

#[test]
fn skew_mesh_diffusion_retains_closed_inventory_and_temperature_bounds() {
    let mut geometry = Mesh2::triangles(1.0, 1.0, 4, 3);
    for p in &mut geometry.nodes {
        p[0] += 0.4 * p[1];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| 0.0;
    let diffusivity: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| if p[1] < 0.5 { 0.002 } else { 0.2 })
        .collect();
    let native = system
        .prepare_scalar_transport(
            &vec![0.0; system.n_u],
            &diffusivity,
            &boundaries(&mesh, &zero, None),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let bounded = native
        .prepare_bounded_transport(
            ScalarBounds {
                lower: 280.0,
                upper: 340.0,
            },
            || false,
        )
        .unwrap();
    let previous: Vec<_> = mesh
        .centroids
        .iter()
        .flat_map(|p| [if p[1] < 0.5 { 280.0 } else { 340.0 }; 3])
        .collect();
    let result = bounded
        .ssprk2_step(&previous, 0.2, BoundedScalarControls::default(), || false)
        .unwrap();
    assert!(result.minimum >= 280.0 - 1e-11 && result.maximum <= 340.0 + 1e-11);
    close(result.integral_after, result.integral_before, 5e-10);
    close(result.balance_residual, 0.0, 5e-10);
    assert!(
        result
            .values
            .iter()
            .zip(&previous)
            .any(|(a, b)| (a - b).abs() > 1.0)
    );
    for boundary in &result.boundaries {
        close(boundary.advection, 0.0, 0.0);
        close(boundary.diffusion, 0.0, 0.0);
        close(boundary.limiter_diffusion, 0.0, 0.0);
    }
    let constant = bounded
        .ssprk2_step(
            &vec![300.0; previous.len()],
            0.2,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    for value in constant.values {
        close(value, 300.0, 5e-10);
    }
}

#[test]
fn shared_sgs_diffusivity_and_reversing_reservoirs_remain_bounded_and_conservative() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = moments(&mesh, &|p| [2.0 * p[1] - 0.7, 0.8 - p[0]]);
    let zero = |_: [f64; 2]| 0.0;
    let donor = |p: [f64; 2]| 0.7 + 0.1 * p[1];
    let bc = boundaries(&mesh, &zero, Some(&donor));
    let n = mesh.tris.len();
    let previous: Vec<_> = mesh
        .centroids
        .iter()
        .flat_map(|p| [if p[0] < 0.5 { 0.0 } else { 1.0 }; 3])
        .collect();
    for model in [
        SgsModel::Smagorinsky { cs: 0.3 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.1 },
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let closure = PlanarSgs::with_model(&system, model, &vec![0.2; n], &vec![0.01; n]).unwrap();
        let coupled = closure
            .prepare_scalar_transport(
                &velocity,
                &vec![0.001; n],
                &vec![0.7; n],
                &bc,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let native = system
            .prepare_scalar_transport(
                &velocity,
                coupled.transport.diffusivity_m2_s(),
                &bc,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let result = coupled
            .transport
            .prepare_bounded_transport(fraction_bounds(), || false)
            .unwrap()
            .ssprk2_step(&previous, 0.1, BoundedScalarControls::default(), || false)
            .unwrap();
        let expected = native
            .prepare_bounded_transport(fraction_bounds(), || false)
            .unwrap()
            .ssprk2_step(&previous, 0.1, BoundedScalarControls::default(), || false)
            .unwrap();
        assert_eq!(result.values, expected.values);
        assert_eq!(
            result.balance_residual.to_bits(),
            expected.balance_residual.to_bits()
        );
        assert!(result.minimum >= -3e-14 && result.maximum <= 1.0 + 3e-14);
        close(result.balance_residual, 0.0, 2e-12);
        assert!(result.boundaries.iter().any(|b| b.advection < 0.0));
        assert!(result.boundaries.iter().any(|b| b.advection > 0.0));
    }
}

#[test]
fn bounded_smooth_advection_refines_in_space() {
    let pulse = |x: f64| {
        if (x - 0.5).abs() < 0.25 {
            (2.0 * std::f64::consts::PI * (x - 0.5)).cos().powi(2)
        } else {
            0.0
        }
    };
    let mut errors = Vec::new();
    for resolution in [8, 16, 32] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, resolution, 3));
        let system = FluxSystem::new(&mesh);
        let velocity = moments(&mesh, &|_| [1.0, 0.0]);
        let zero = |_: [f64; 2]| 0.0;
        let native = system
            .prepare_scalar_transport(
                &velocity,
                &vec![0.0; mesh.tris.len()],
                &boundaries(&mesh, &zero, Some(&zero)),
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let bounded = native
            .prepare_bounded_transport(fraction_bounds(), || false)
            .unwrap();
        let result = bounded
            .ssprk2_step(
                &nodal(&mesh, &|p| pulse(p[0])),
                0.05,
                BoundedScalarControls::default(),
                || false,
            )
            .unwrap();
        let mut error = 0.0;
        for (cell, tri) in mesh.tris.iter().enumerate() {
            let points = tri.map(|v| mesh.verts[v]);
            for (p, weight) in fs_flux::bdm::tri_quad(points, mesh.areas[cell]) {
                let mut actual = 0.0;
                for i in 0..3 {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    let cross = (points[j][0] - p[0]) * (points[k][1] - p[1])
                        - (points[j][1] - p[1]) * (points[k][0] - p[0]);
                    actual += result.values[3 * cell + i] * cross / (2.0 * mesh.areas[cell]);
                }
                error += weight * (actual - pulse(p[0] - 0.05)).powi(2);
            }
        }
        errors.push(error.sqrt());
        close(result.balance_residual, 0.0, 3e-12);
    }
    assert!(
        errors[0] / errors[1] > 1.7 && errors[1] / errors[2] > 1.7,
        "smooth bounded L2 errors: {errors:?}"
    );
}

#[test]
fn bounded_work_and_late_cancellation_are_atomic_and_boundary_sources_are_explicit() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| 0.0;
    let one = |_: [f64; 2]| 1.0;
    let velocity = moments(&mesh, &|_| [0.0, 1.0]);
    let native = system
        .prepare_scalar_transport(
            &velocity,
            &vec![0.02; mesh.tris.len()],
            &boundaries(&mesh, &zero, Some(&one)),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let bounded = native
        .prepare_bounded_transport(fraction_bounds(), || false)
        .unwrap();
    let initial = nodal(&mesh, &|p| p[1]);
    let snapshot = initial.clone();
    let mut polls = 0;
    let good = bounded
        .ssprk2_step(&initial, 0.2, BoundedScalarControls::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut interrupted = 0;
    assert!(matches!(
        bounded.ssprk2_step(&initial, 0.2, BoundedScalarControls::default(), || {
            interrupted += 1;
            interrupted == polls / 2
        }),
        Err(BoundedScalarError::Native(
            fs_flux::ns::transport::scalar::ScalarTransportError::Cancelled
        ))
    ));
    let replay = bounded
        .ssprk2_step(&initial, 0.2, BoundedScalarControls::default(), || false)
        .unwrap();
    assert_eq!(good.values, replay.values);
    assert_eq!(initial, snapshot);
    assert!(matches!(
        bounded.ssprk2_step(
            &initial,
            0.2,
            BoundedScalarControls {
                maximum_substeps: 1,
                ..BoundedScalarControls::default()
            },
            || false
        ),
        Err(BoundedScalarError::SubstepLimit { .. })
    ));
    assert!(
        bounded
            .ssprk2_step(
                &vec![-0.1; initial.len()],
                0.1,
                BoundedScalarControls::default(),
                || false
            )
            .is_err()
    );
    let flux = system
        .prepare_scalar_transport(
            &velocity,
            &vec![0.02; mesh.tris.len()],
            &boundaries(&mesh, &one, Some(&one)),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    assert!(matches!(
        flux.prepare_bounded_transport(fraction_bounds(), || false),
        Err(BoundedScalarError::PrescribedFlux { .. })
    ));
    let invalid_trace = |p: [f64; 2]| 1.5 * p[0] - 0.25;
    let face_mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let face_system = FluxSystem::new(&face_mesh);
    let trace = face_system
        .prepare_scalar_transport(
            &moments(&face_mesh, &|_| [0.0, 1.0]),
            &[0.0; 2],
            &boundaries(&face_mesh, &zero, Some(&invalid_trace)),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    assert!(matches!(
        trace.prepare_bounded_transport(fraction_bounds(), || false),
        Err(BoundedScalarError::BoundaryBounds { .. })
    ));
}

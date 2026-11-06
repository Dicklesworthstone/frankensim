//! Native planar SGS stress: physical scaling, endpoint coupling and refusal.
#![cfg(all(feature = "continuum", feature = "experimental-sgs"))]

use fs_flux::bdm::edge_gauss_pub;
use fs_flux::ns::boundary::{NavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};
use fs_flux::ns::sgs::{PlanarSgsError, PlanarSmagorinsky};
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = Vec::with_capacity(2 * mesh.edges.len());
    for edge in &mesh.edges {
        let a = mesh.verts[edge.verts.0];
        let b = mesh.verts[edge.verts.1];
        let mut pair = [0.0; 2];
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = velocity(point);
            let normal = u[0] * edge.normal[0] + u[1] * edge.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0]) + (point[1] - a[1]) * (b[1] - a[1]))
                / edge.len.powi(2)
                - 0.5;
            pair[0] += weight * normal / edge.len;
            pair[1] += weight * normal * s / edge.len;
        }
        values.extend(pair);
    }
    values
}

fn velocity_boundaries<'a>(
    mesh: &TriMesh,
    value: &'a dyn Fn([f64; 2]) -> [f64; 2],
) -> Vec<NavierStokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| NavierStokesBoundary::Velocity { edge, value })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn affine_shear_uses_each_authored_width_and_area_with_physical_units() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(3.0, 0.7, 3, 2));
    let system = FluxSystem::new(&mesh);
    let widths: Vec<_> = (0..mesh.tris.len()).map(|cell| 0.1 + 0.01 * cell as f64).collect();
    let molecular: Vec<_> = (0..mesh.tris.len()).map(|cell| 0.02 + 0.001 * cell as f64).collect();
    let model = PlanarSmagorinsky::new(&system, 0.2, &widths, &molecular).unwrap();
    let doubled = PlanarSmagorinsky::new(
        &system,
        0.2,
        &widths.iter().map(|width| 2.0 * width).collect::<Vec<_>>(),
        &molecular,
    )
    .unwrap();
    let velocity = moments(&mesh, &|p| [0.31 + 2.5 * p[1], -0.21]);
    let field = model.field(&velocity).unwrap();
    let wider = doubled.field(&velocity).unwrap();
    let mut expected_integral = 0.0;
    for (cell, diagnostic) in field.cells.iter().enumerate() {
        // Independent affine derivative: D_xy = D_yx = 1.25,
        // so sqrt(2 D:D) = 2.5, regardless of cell shape or translation.
        let eddy = (0.2 * widths[cell]).powi(2) * 2.5;
        let transfer = eddy * 6.25;
        close(diagnostic.strain_rate_per_s, 2.5, 3e-12);
        close(diagnostic.eddy_viscosity_m2_s, eddy, 3e-14);
        close(diagnostic.effective_viscosity_m2_s, molecular[cell] + eddy, 3e-14);
        close(diagnostic.dissipation_m2_s3, transfer, 3e-12);
        close(wider.cells[cell].eddy_viscosity_m2_s, 4.0 * eddy, 1e-13);
        expected_integral += mesh.areas[cell] * transfer;
    }
    close(field.integrated_dissipation_m4_s3, expected_integral, 3e-12);
    close(wider.integrated_dissipation_m4_s3, 4.0 * expected_integral, 1e-11);
}

#[test]
fn rigid_rotation_and_galilean_translation_do_not_create_modeled_strain() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.25; mesh.tris.len()];
    let molecular = vec![0.01; mesh.tris.len()];
    let model = PlanarSmagorinsky::new(&system, 0.2, &widths, &molecular).unwrap();
    let rotation = model
        .field(&moments(&mesh, &|p| [0.3 - 0.25 * p[1], -0.2 + 0.25 * p[0]]))
        .unwrap();
    for cell in rotation.cells {
        assert!(cell.strain_rate_per_s < 1e-13);
        assert!(cell.eddy_viscosity_m2_s < 1e-15);
    }
    let shear = model.field(&moments(&mesh, &|p| [0.5 * p[1], 0.0])).unwrap();
    let translated = model
        .field(&moments(&mesh, &|p| [0.5 * p[1] + 0.7, -0.4]))
        .unwrap();
    for (a, b) in shear.cells.iter().zip(&translated.cells) {
        close(a.eddy_viscosity_m2_s, b.eddy_viscosity_m2_s, 1e-14);
        close(a.dissipation_m2_s3, b.dissipation_m2_s3, 1e-14);
    }
}

#[test]
fn bdf1_and_bdf2_close_sgs_and_convection_at_the_accepted_traction_flow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.7; mesh.tris.len()];
    let molecular = vec![0.04; mesh.tris.len()];
    let model = PlanarSmagorinsky::new(&system, 0.5, &widths, &molecular).unwrap();
    let params = FluxParams { nu: f64::NAN, tol: 1e-12, ..FluxParams::default() };
    let controls = NavierStokesControls {
        maximum_iterations: 80,
        relative_tolerance: 1e-10,
        relaxation: 0.8,
    };
    let dt = 0.02;
    let old = |p: [f64; 2]| [0.1 * p[0] + 0.1 * p[1] + 0.2, -0.1 * p[1] + 0.05];
    let older_field = |p: [f64; 2]| [-0.05 * p[0] + 0.15 * p[1] + 0.1, 0.05 * p[1] + 0.03];
    let target = |p: [f64; 2]| [0.3 * p[0] + 0.2 * p[1] + 0.4, -0.3 * p[1] + 0.1];
    let previous = moments(&mesh, &old);
    let older = moments(&mesh, &older_field);
    let previous_snapshot = previous.clone();
    let older_snapshot = older.clone();
    let rate = 0.4_f64.sqrt();
    let length_squared = 0.35_f64.powi(2);
    let eddy = length_squared * rate;
    let nu = 0.04 + eddy;
    let traction = |p: [f64; 2], n: [f64; 2]| {
        let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
        [
            (0.6 * nu - pressure) * n[0] + 0.2 * nu * n[1],
            0.2 * nu * n[0] + (-0.6 * nu - pressure) * n[1],
        ]
    };
    let boundaries: Vec<_> = mesh.edges.iter().enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[0].abs() > 0.5 {
                NavierStokesBoundary::Open {
                    edge, traction: &traction, backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                NavierStokesBoundary::Velocity { edge, value: &target }
            }
        })
        .collect();
    for second_order in [false, true] {
        // Constant target strain makes div(2 nu_eff D) zero. Convection and
        // pressure gradients below are independent derivatives of target.
        let force = |p: [f64; 2]| {
            let u = target(p);
            let v = old(p);
            let w = older_field(p);
            let derivative: [f64; 2] = core::array::from_fn(|axis| {
                if second_order {
                    (1.5 * u[axis] - 2.0 * v[axis] + 0.5 * w[axis]) / dt
                } else {
                    (u[axis] - v[axis]) / dt
                }
            });
            [
                derivative[0] + 0.09 * p[0] + 0.14 + 2.0,
                derivative[1] + 0.09 * p[1] - 0.03 - 3.0,
            ]
        };
        let result = model.bdf2_step(
            params, &force, &boundaries, &previous,
            second_order.then_some(older.as_slice()), dt, controls,
        ).unwrap();
        let (error, divergence) = system.velocity_error(&result.step.solution.field.x, &target);
        assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
        assert!(result.step.nonlinear_residual <= controls.relative_tolerance);
        for cell in &result.sgs.cells {
            close(cell.strain_rate_per_s, rate, 3e-8);
            close(cell.eddy_viscosity_m2_s, eddy, 4e-9);
            close(cell.dissipation_m2_s3, eddy * 0.4, 4e-9);
        }
        close(result.sgs.integrated_dissipation_m4_s3, eddy * 0.4, 4e-9);
        // Reassemble independently at the accepted velocity; a lagged strain
        // must fail even though the momentum advection remains at the endpoint.
        let mass_history: Vec<_> = previous.iter().zip(&older)
            .map(|(&recent, &old)| if second_order { recent + (recent - old) / 3.0 } else { recent })
            .collect();
        let effective_dt = if second_order { (2.0 / 3.0) * dt } else { dt };
        let initial_rates = system.shear_rates(&previous);
        for accepted in [true, false] {
            let operator = system.assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                params,
                &|cell, _| if accepted {
                    result.sgs.cells[cell].effective_viscosity_m2_s
                } else {
                    0.04 + length_squared * initial_rates[cell]
                },
                &force, &boundaries, &result.step.solution.field.x[..system.n_u],
                &mass_history, effective_dt,
            ).unwrap();
            let mut action = vec![0.0; system.n];
            operator.matrix().spmv(&result.step.solution.field.x, &mut action);
            let defect = action.iter().zip(operator.rhs())
                .fold(0.0_f64, |norm, (ax, rhs)| norm.hypot(ax - rhs));
            if accepted {
                assert!(defect < 3e-8, "accepted endpoint defect {defect}");
            } else {
                assert!(defect > 1e-5, "lagged SGS strain was not detected: {defect}");
            }
        }
    }
    assert_eq!(previous, previous_snapshot);
    assert_eq!(older, older_snapshot);
}

#[test]
fn disabled_closure_matches_the_native_molecular_solver_bitwise() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.5; mesh.tris.len()];
    let molecular = vec![0.1; mesh.tris.len()];
    let model = PlanarSmagorinsky::new(&system, 0.0, &widths, &molecular).unwrap();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let force = |p: [f64; 2]| [p[1] - 0.5, 0.5 - p[0]];
    let previous = moments(&mesh, &|p| [0.2 * p[1], 0.1 * p[0]]);
    let older = moments(&mesh, &|p| [0.1 * p[1], 0.05 * p[0]]);
    let params = FluxParams { nu: f64::NAN, ..FluxParams::default() };
    let controls = NavierStokesControls::default();
    for old in [None, Some(older.as_slice())] {
        let reference = system.bdf2_navier_stokes_mixed_rheology_step(
            params, &|cell, _, _| molecular[cell], &force, &boundaries,
            &previous, old, 0.01, controls,
        ).unwrap();
        let result = model.bdf2_step(
            params, &force, &boundaries, &previous, old, 0.01, controls,
        ).unwrap();
        assert_eq!(result.step.solution.field.x, reference.solution.field.x);
        assert_eq!(result.step.nonlinear_residual, reference.nonlinear_residual);
        assert!(result.sgs.cells.iter().all(|cell| {
            cell.eddy_viscosity_m2_s == 0.0 && cell.dissipation_m2_s3 == 0.0
        }));
    }
    let first = model.bdf1_step(
        params, &force, &boundaries, &previous, 0.01, controls,
    ).unwrap();
    let startup = model.bdf2_step(
        params, &force, &boundaries, &previous, None, 0.01, controls,
    ).unwrap();
    assert_eq!(first.step.solution.field.x, startup.step.solution.field.x);
}

#[test]
fn unforced_no_slip_bdf1_damps_resolved_energy_with_positive_sgs_transfer() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let params = FluxParams { nu: 0.02, sigma: 20.0, ..FluxParams::default() };
    let controls = NavierStokesControls {
        relative_tolerance: 1e-10, relaxation: 0.8, ..NavierStokesControls::default()
    };
    // Start from a native converged, divergence-free flow instead of assuming
    // arbitrary interpolated moments obey the discrete no-slip constraints.
    let driven = system.bdf1_navier_stokes_mixed_step(
        params, &|p| [p[1] - 0.5, 0.5 - p[0]], &boundaries,
        &vec![0.0; system.n_u], 0.1, controls,
    ).unwrap();
    let previous = &driven.solution.field.x[..system.n_u];
    let before = system.kinetic_energy(previous);
    let model = PlanarSmagorinsky::new(
        &system, 0.5, &vec![0.5; mesh.tris.len()], &vec![0.02; mesh.tris.len()],
    ).unwrap();
    let result = model.bdf1_step(
        params, &zero, &boundaries, previous, 0.1, controls,
    ).unwrap();
    let after = system.kinetic_energy(&result.step.solution.field.x);
    assert!(before > 0.0 && after >= 0.0 && after < before, "{before} -> {after}");
    assert!(result.sgs.integrated_dissipation_m4_s3 > 0.0);
    assert!(result.sgs.cells.iter().all(|cell| cell.dissipation_m2_s3 >= 0.0));
}

#[test]
fn invalid_fields_and_unrepresentable_transport_refuse_without_changing_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.5; mesh.tris.len()];
    let molecular = vec![0.01; mesh.tris.len()];
    for coefficient in [-0.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            PlanarSmagorinsky::new(&system, coefficient, &widths, &molecular),
            Err(PlanarSgsError::InvalidInput(_))
        ));
    }
    assert!(PlanarSmagorinsky::new(&system, 0.2, &[], &molecular).is_err());
    assert!(PlanarSmagorinsky::new(&system, 0.2, &vec![0.0; widths.len()], &molecular).is_err());
    assert!(PlanarSmagorinsky::new(&system, 0.2, &widths, &vec![f64::NAN; widths.len()]).is_err());
    let model = PlanarSmagorinsky::new(&system, 0.2, &widths, &molecular).unwrap();
    assert!(model.field(&vec![0.0; system.n_u - 1]).is_err());
    assert!(model.field(&vec![f64::NAN; system.n_u]).is_err());

    // A representable eddy viscosity must not silently publish zero transfer
    // when its strictly positive dissipation underflows the physical units.
    let underflowing = PlanarSmagorinsky::new(
        &system, 1.0, &vec![1e-100; mesh.tris.len()], &molecular,
    ).unwrap();
    let tiny_shear = moments(&mesh, &|p| [1e-50 * p[1], 0.0]);
    assert!(matches!(
        underflowing.field(&tiny_shear),
        Err(PlanarSgsError::UnrepresentableDissipation { .. })
    ));

    let overflowing = PlanarSmagorinsky::new(
        &system, 1.0, &vec![1e160; mesh.tris.len()], &molecular,
    ).unwrap();
    let shear = |p: [f64; 2]| [p[1], 0.0];
    let previous = moments(&mesh, &shear);
    let snapshot = previous.clone();
    let boundaries = velocity_boundaries(&mesh, &shear);
    assert!(matches!(
        overflowing.bdf1_step(
            FluxParams::default(), &|_| [0.0; 2], &boundaries, &previous,
            0.01, NavierStokesControls::default(),
        ),
        Err(PlanarSgsError::UnrepresentableTransport { .. })
    ));
    assert_eq!(previous, snapshot);
}

use super::*;
use crate::{FluxSystem, TriMesh};
use fs_material::fluid_transport::rheology::CarreauYasuda;
use fs_solid::Mesh2;

fn zero(_: [f64; 2]) -> [f64; 2] { [0.0; 2] }
fn reference() -> FluidTransportState {
    FluidTransportState::try_new(2.0, 4.0, 0.5, 0.8).unwrap()
}
fn controls() -> NavierStokesControls {
    NavierStokesControls { maximum_iterations: 100, relative_tolerance: 1e-10, relaxation: 0.8 }
}
fn closed(mesh: &TriMesh) -> Vec<MaterialNavierStokesBoundary<'_>> {
    (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e)).map(|edge| {
        MaterialNavierStokesBoundary::Velocity { edge, value: &zero }
    }).collect()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} versus {b}");
}

#[test]
fn cell_load_preserves_hydrostatic_pressure_and_original_constraint_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let source = reference();
    let material = MaterialFluxSystem::new(&system, &source, super::super::super::MaterialFluxNumerics {
        linear_tolerance: 1e-12, ..super::super::super::MaterialFluxNumerics::default()
    }).unwrap();
    let velocity = vec![0.0; system.n_u];
    let force = [0.25, -0.5];
    let boundaries = closed(&mesh);
    let (flow, forces) = material.bdf1_navier_stokes_mixed_rheology_with_cell_forces(
        &zero, &boundaries, &velocity, 0.2, controls(), |_, _| Ok((source, force)),
    ).unwrap();
    assert_eq!(forces, vec![force; mesh.tris.len()]);
    for &u in flow.step.solution.velocity_dofs() { near(u, 0.0, 2e-10); }
    // The exact gradient load does no solenoidal work. P0 pressure equals the
    // original affine hydrostatic pressure's cell average, in its mean gauge.
    for (cell, &p) in flow.step.solution.pressure_pa().iter().enumerate() {
        let x = mesh.centroids[cell];
        near(p, force[0]*(x[0]-0.5) + force[1]*(x[1]-0.5), 2e-9);
    }
    let edges: Vec<_> = (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e)).collect();
    let reaction = flow.step.solution.boundary_loads(&edges, [0.0; 2], &[Some(0.0)]).unwrap();
    for axis in 0..2 { near(reaction.stress.force_n_per_m[axis], force[axis], 2e-9); }
    assert!(flow.step.solution.viscous_dissipation().total_power_w_per_m() < 1e-16);
    assert_eq!(velocity, vec![0.0; system.n_u]);
}

#[test]
fn source_cell_identity_matches_original_spatial_loading_with_nonlinear_rheology() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let source = reference();
    let material = MaterialFluxSystem::new(&system, &source, super::super::super::MaterialFluxNumerics {
        linear_tolerance: 1e-12, ..super::super::super::MaterialFluxNumerics::default()
    }).unwrap();
    let boundaries = closed(&mesh);
    let old = vec![0.0; system.n_u];
    // The discontinuity follows original cell faces. Opposing vertical forces
    // have nonzero curl and drive an actual closed-cavity circulation.
    let spatial = |p: [f64; 2]| [0.0, if p[0] < 0.5 { 0.5 } else { -0.5 }];
    let curve = CarreauYasuda::new(0.8, 0.1, 2.0, 0.5, 2.0, [0.0, 1e4]).unwrap();
    let resolve = |point: MaterialFluxRheologyPoint, _: &FluidTransportState| {
        FluidTransportState::try_new(2.0, 4.0, 0.5, curve.viscosity_pa_s(point.shear_rate_per_s)?)
    };
    let original = material.bdf1_navier_stokes_mixed_rheology_step(
        &spatial, &boundaries, &old, 0.1, controls(), resolve,
    ).unwrap();
    let (coupled, forces) = material.bdf1_navier_stokes_mixed_rheology_with_cell_forces(
        &zero, &boundaries, &old, 0.1, controls(), |point, reference| {
            Ok((resolve(point, reference)?, spatial(point.position_m)))
        },
    ).unwrap();
    assert!(coupled.step.solution.velocity_dofs().iter().any(|u| u.abs() > 1e-6));
    for (&a, &b) in coupled.step.solution.velocity_dofs().iter().zip(original.step.solution.velocity_dofs()) {
        near(a, b, 2e-10);
    }
    for (cell, &force) in forces.iter().enumerate() {
        assert_eq!(force, spatial(mesh.centroids[cell]));
        near(coupled.samples()[cell].transport().dynamic_viscosity_pa_s(),
            original.samples()[cell].transport().dynamic_viscosity_pa_s(), 2e-10);
    }
    assert!(coupled.step.nonlinear_residual < controls().relative_tolerance);
}

#[test]
fn a_late_source_or_nonfinite_cell_force_never_returns_a_partial_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let source = reference();
    let material = MaterialFluxSystem::new(&system, &source, super::super::super::MaterialFluxNumerics::default()).unwrap();
    let boundaries = closed(&mesh);
    let old = vec![0.0; system.n_u];
    for refuse_source in [false, true] {
        let result = material.bdf1_navier_stokes_mixed_rheology_with_cell_forces(
            &zero, &boundaries, &old, 0.1, controls(), |point, _| {
                if point.cell == 1 {
                    if refuse_source { return Err(MaterialError::State { what: "original source domain".into() }); }
                    return Ok((source, [0.0, f64::NAN]));
                }
                Ok((source, [0.0, 0.5]))
            },
        );
        assert!(matches!(result, Err(MaterialFluxError::Material(_))));
        assert_eq!(old, vec![0.0; system.n_u]);
    }
}

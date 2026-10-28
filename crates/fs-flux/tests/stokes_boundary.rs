//! Mixed steady-Stokes boundaries: physical traction, pressure datum and flux.

use fs_flux::{
    FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, StokesPressureReference, TriMesh,
};
use fs_solid::Mesh2;
use std::cell::Cell;

fn velocity_boundaries<'a>(
    mesh: &TriMesh,
    value: &'a dyn Fn([f64; 2]) -> [f64; 2],
) -> Vec<StokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity { edge, value })
        .collect()
}

/// Independent fully developed solution: u=y(1-y), p=5.3-0.8x,
/// nu=0.4 on [0,2]x[0,1], hence -nu*u''+p'=0 and Q=1/6.
/// The end's TOTAL stress includes its tangential viscous traction.
#[test]
fn pressure_driven_channel_refines_and_conserves_mass() {
    let mut errors = Vec::new();
    for n in [3, 6] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2 * n, n));
        let system = FluxSystem::new(&mesh);
        let params = FluxParams { nu: 0.4, ..FluxParams::default() };
        let wall = |_: [f64; 2]| [0.0; 2];
        let traction = |x: [f64; 2], normal: [f64; 2]| {
            let pressure = 5.3 - 0.8 * x[0];
            let shear = 0.4 * (1.0 - 2.0 * x[1]);
            [
                -pressure * normal[0] + shear * normal[1],
                shear * normal[0] - pressure * normal[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    StokesBoundary::TotalTraction { edge, value: &traction }
                } else {
                    StokesBoundary::Velocity { edge, value: &wall }
                }
            })
            .collect();
        let assembly = system
            .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
            .expect("admit pressure-driven channel");
        let solution = assembly.solve().expect("solve pressure-driven channel");
        assert_eq!(solution.pressure_reference, StokesPressureReference::TractionDatum);
        let (error, divergence) = system.velocity_error(
            &solution.field.x,
            &|x| [x[1] * (1.0 - x[1]), 0.0],
        );
        assert!(divergence < 2e-10, "n={n}, divergence={divergence:e}");
        assert!(solution.field.rel_residual <= params.tol);
        let mut inward = 0.0;
        let mut outward = 0.0;
        for (edge, geometry) in mesh.edges.iter().enumerate() {
            if !mesh.is_boundary(edge) {
                continue;
            }
            let flux = geometry.len * solution.field.x[2 * edge];
            if geometry.normal[0] < -0.5 {
                inward -= flux;
            } else if geometry.normal[0] > 0.5 {
                outward += flux;
            }
        }
        assert!((inward - outward).abs() < 2e-11);
        if n == 6 {
            assert!(error < 0.02, "channel velocity L2={error:e}");
            assert!((outward - 1.0 / 6.0).abs() < 0.008, "outlet Q={outward:e}");
            let mean = solution.pressure.iter().zip(&mesh.areas)
                .map(|(pressure, area)| pressure * area).sum::<f64>()
                / mesh.areas.iter().sum::<f64>();
            assert!((mean - 4.5).abs() < 0.05, "pressure datum mean={mean}");
        }
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 2.5, "channel errors={errors:?}");
}

/// One velocity edge family anchors affine shear; every other edge carries
/// the exact total stress. A pressure shift must change traction by -dp*n
/// while leaving the velocity and per-cell continuity unchanged.
#[test]
fn mixed_affine_traction_patch_retains_pressure_datum_and_shift() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams { nu: 0.7, ..FluxParams::default() };
    let velocity = |x: [f64; 2]| [0.2 + 1.3 * x[1], 0.0];
    let mut previous: Option<Vec<f64>> = None;
    for pressure in [7.0, 32.0] {
        let traction = |_: [f64; 2], normal: [f64; 2]| {
            let shear = params.nu * 1.3;
            [
                -pressure * normal[0] + shear * normal[1],
                shear * normal[0] - pressure * normal[1],
            ]
        };
        let boundaries: Vec<_> = mesh.edges.iter().enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0] < -0.5 {
                    StokesBoundary::Velocity { edge, value: &velocity }
                } else {
                    StokesBoundary::TotalTraction { edge, value: &traction }
                }
            })
            .collect();
        let solution = system
            .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
            .expect("admit mixed patch")
            .solve()
            .expect("solve mixed patch");
        let (error, divergence) = system.velocity_error(&solution.field.x, &velocity);
        assert!(error < 5e-11 && divergence < 5e-11, "error={error:e}, div={divergence:e}");
        assert!(solution.pressure.iter().all(|value| (value - pressure).abs() < 5e-10));
        assert_eq!(solution.pressure_reference, StokesPressureReference::TractionDatum);
        if let Some(old) = previous.as_ref() {
            assert!(old.iter().zip(&solution.field.x[..system.n_u])
                .all(|(a, b)| (a - b).abs() < 5e-11));
        }
        previous = Some(solution.field.x[..system.n_u].to_vec());
    }
}

/// New all-velocity admission delegates to the same symmetric-stress kernel;
/// it cannot perturb the established matrix, force or accepted solution.
#[test]
fn all_velocity_partition_preserves_existing_operator_bits() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams { nu: 0.37, ..FluxParams::default() };
    let velocity = |x: [f64; 2]| [0.3 - x[1], x[0] - 0.2];
    let force = |x: [f64; 2]| [x[0] - 0.4, -2.0 * x[1]];
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let (matrix, rhs) = system.assemble_variable_viscosity(
        params, &|_, _| params.nu, &force, &velocity, None,
    );
    let mixed = system.assemble_stokes_mixed(params, &force, &boundaries)
        .expect("compatible prescribed velocities");
    assert_eq!(mixed.pressure_reference(), StokesPressureReference::MeanZero);
    for row in 0..matrix.nrows() {
        let (old_columns, old_values) = matrix.row(row);
        let (new_columns, new_values) = mixed.matrix().row(row);
        assert_eq!(old_columns, new_columns);
        assert_eq!(
            old_values.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            new_values.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        );
    }
    assert_eq!(
        rhs.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        mixed.rhs().iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
    );
    let old = system.solve_linear(&matrix, &rhs, params);
    let new = mixed.solve().expect("compatible all-velocity solve");
    assert_eq!(
        old.x.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        new.field.x.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
    );
    assert_eq!(new.pressure, system.pressure(&old.x));
}

#[test]
fn boundary_ownership_and_closed_domain_flux_are_checked_before_solve() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let calls = Cell::new(0_usize);
    let velocity = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let mut duplicate = boundaries.clone();
    duplicate.push(boundaries[0]);
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &duplicate),
        Err(StokesBoundaryError::DuplicateEdge { .. })
    ));
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries[1..]),
        Err(StokesBoundaryError::MissingEdge { .. })
    ));
    let interior = (0..mesh.edges.len()).find(|edge| !mesh.is_boundary(*edge)).unwrap();
    let mut wrong_owner = boundaries.clone();
    wrong_owner.push(StokesBoundary::Velocity { edge: interior, value: &velocity });
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &wrong_owner),
        Err(StokesBoundaryError::InteriorEdge { .. })
    ));
    assert_eq!(calls.get(), 0, "ownership must fail before callback evaluation");

    let outward = |x: [f64; 2]| x;
    let incompatible = velocity_boundaries(&mesh, &outward);
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &incompatible),
        Err(StokesBoundaryError::IncompatibleFlux { .. })
    ));
    let zero_traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let unanchored: Vec<_> = mesh.edges.iter().enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction { edge, value: &zero_traction })
        .collect();
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &unanchored),
        Err(StokesBoundaryError::NoVelocityBoundary)
    ));
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [f64::NAN, 0.0], &boundaries),
        Err(StokesBoundaryError::NonFiniteForce)
    ));
}

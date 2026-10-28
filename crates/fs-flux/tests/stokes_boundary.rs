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
        let params = FluxParams {
            nu: 0.4,
            ..FluxParams::default()
        };
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
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::Velocity { edge, value: &wall }
                }
            })
            .collect();
        let assembly = system
            .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
            .expect("admit pressure-driven channel");
        let solution = assembly.solve().expect("solve pressure-driven channel");
        assert_eq!(
            solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        let (error, divergence) =
            system.velocity_error(&solution.field.x, &|x| [x[1] * (1.0 - x[1]), 0.0]);
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
            let mean = solution
                .pressure
                .iter()
                .zip(&mesh.areas)
                .map(|(pressure, area)| pressure * area)
                .sum::<f64>()
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
    let params = FluxParams {
        nu: 0.7,
        ..FluxParams::default()
    };
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
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0] < -0.5 {
                    StokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let solution = system
            .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
            .expect("admit mixed patch")
            .solve()
            .expect("solve mixed patch");
        let (error, divergence) = system.velocity_error(&solution.field.x, &velocity);
        assert!(
            error < 5e-11 && divergence < 5e-11,
            "error={error:e}, div={divergence:e}"
        );
        assert!(
            solution
                .pressure
                .iter()
                .all(|value| (value - pressure).abs() < 5e-10)
        );
        assert_eq!(
            solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        if let Some(old) = previous.as_ref() {
            assert!(
                old.iter()
                    .zip(&solution.field.x[..system.n_u])
                    .all(|(a, b)| (a - b).abs() < 5e-11)
            );
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
    let params = FluxParams {
        nu: 0.37,
        ..FluxParams::default()
    };
    let velocity = |x: [f64; 2]| [0.3 - x[1], x[0] - 0.2];
    let force = |x: [f64; 2]| [x[0] - 0.4, -2.0 * x[1]];
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let (matrix, rhs) =
        system.assemble_variable_viscosity(params, &|_, _| params.nu, &force, &velocity, None);
    let mixed = system
        .assemble_stokes_mixed(params, &force, &boundaries)
        .expect("compatible prescribed velocities");
    assert_eq!(
        mixed.pressure_reference(),
        StokesPressureReference::MeanZero
    );
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
    let interior = (0..mesh.edges.len())
        .find(|edge| !mesh.is_boundary(*edge))
        .unwrap();
    let mut wrong_owner = boundaries.clone();
    wrong_owner.push(StokesBoundary::Velocity {
        edge: interior,
        value: &velocity,
    });
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &wrong_owner),
        Err(StokesBoundaryError::InteriorEdge { .. })
    ));
    assert_eq!(
        calls.get(),
        0,
        "ownership must fail before callback evaluation"
    );

    let outward = |x: [f64; 2]| x;
    let incompatible = velocity_boundaries(&mesh, &outward);
    assert!(matches!(
        system.assemble_stokes_mixed(params, &|_| [0.0; 2], &incompatible),
        Err(StokesBoundaryError::IncompatibleFlux { .. })
    ));
    let zero_traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let unanchored: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &zero_traction,
        })
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

/// A natural boundary removes the pressure pin. Above the dense cutoff the
/// block Schur solve must recover that physical datum, not a mean-zero shift.
#[test]
fn sparse_mixed_traction_retains_pressure_datum_and_short_budgets_refuse() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 16, 16));
    let system = FluxSystem::new(&mesh);
    assert!(system.n > 1500);
    let params = FluxParams {
        nu: 0.7,
        tol: 1e-13,
        max_iters: 4000,
        ..FluxParams::default()
    };
    let pressure = 7.0;
    let velocity = |x: [f64; 2]| [0.2 + 1.3 * x[1], 0.0];
    let traction = |_: [f64; 2], normal: [f64; 2]| {
        let shear = params.nu * 1.3;
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
            if geometry.normal[0] < -0.5 {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let assembly = system
        .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
        .unwrap();
    assert_eq!(
        assembly.pressure_reference(),
        StokesPressureReference::TractionDatum
    );
    let solution = assembly.solve().unwrap();
    let (error, divergence) = system.velocity_error(&solution.field.x, &velocity);
    let pressure_error = solution
        .pressure
        .iter()
        .fold(0.0_f64, |error, value| error.max((value - pressure).abs()));
    eprintln!(
        "sparse mixed Stokes: dofs={}, iterations={}, L2={error:e}, div={divergence:e}, pressure={pressure_error:e}, residual={:e}, scaled={:e}",
        system.n,
        solution.field.iters,
        solution.field.rel_residual,
        solution.field.scaled_rel_residual,
    );
    assert!(solution.field.iters > 0 && solution.field.iters <= params.max_iters);
    assert!(solution.field.converged(params.tol));
    assert!(error < 2e-10 && divergence < 2e-9);
    assert!(pressure_error < 2e-9);

    // A partial restart is admitted and bounded; it must not be rounded up
    // to an old 100-iteration cycle or published as a converged flow.
    let short = system
        .assemble_stokes_mixed(
            FluxParams {
                max_iters: 1,
                ..params
            },
            &|_| [0.0; 2],
            &boundaries,
        )
        .unwrap();
    assert!(matches!(
        short.solve(),
        Err(StokesBoundaryError::NotConverged {
            iterations,
            relative_residual,
            scaled_relative_residual,
        }) if iterations <= 1
            && (relative_residual > params.tol || scaled_relative_residual > params.tol)
    ));
}

/// A pressure gradient drives two aligned viscosity layers. Integrating
/// nu(y) u'(y) = C - G*y gives a continuous piecewise quadratic profile,
/// continuous shear stress, and an analytic volume flux. Both open ends
/// carry total stress; no inlet velocity supplies the answer.
#[test]
fn variable_viscosity_pressure_channel_refines_with_continuous_stress() {
    let height = 0.5;
    let lower_nu = 0.2;
    let upper_nu = 4.0;
    let gradient = 0.8;
    let stress_constant = gradient
        * (height * height / (2.0 * lower_nu) + (1.0 - height * height) / (2.0 * upper_nu))
        / (height / lower_nu + (1.0 - height) / upper_nu);
    let offset = (stress_constant * height - 0.5 * gradient * height * height)
        * (1.0 / lower_nu - 1.0 / upper_nu);
    let exact = |point: [f64; 2]| {
        let y = point[1];
        let primitive = stress_constant * y - 0.5 * gradient * y * y;
        [
            if y <= height {
                primitive / lower_nu
            } else {
                primitive / upper_nu + offset
            },
            0.0,
        ]
    };
    let exact_flux = (0.5 * stress_constant * height * height
        - gradient * height * height * height / 6.0)
        / lower_nu
        + (0.5 * stress_constant * (1.0 - height * height)
            - gradient * (1.0 - height * height * height) / 6.0)
            / upper_nu
        + offset * (1.0 - height);
    let mut errors = Vec::new();
    let mut pressure_errors = Vec::new();
    for n in [8, 16] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2 * n, n));
        let system = FluxSystem::new(&mesh);
        let params = FluxParams {
            nu: f64::NAN,
            max_iters: 6000,
            ..FluxParams::default()
        };
        let viscosity = |cell: usize, _: [f64; 2]| {
            let y = mesh.tris[cell]
                .iter()
                .map(|&v| mesh.verts[v][1])
                .sum::<f64>()
                / 3.0;
            if y < height { lower_nu } else { upper_nu }
        };
        let wall = |_: [f64; 2]| [0.0; 2];
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = 5.3 - gradient * point[0];
            let shear = stress_constant - gradient * point[1];
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
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::Velocity { edge, value: &wall }
                }
            })
            .collect();
        let solution = system
            .assemble_stokes_mixed_variable_viscosity(
                params,
                &viscosity,
                &|_| [0.0; 2],
                &boundaries,
            )
            .unwrap()
            .solve()
            .unwrap();
        let (error, divergence) = system.velocity_error(&solution.field.x, &exact);
        assert!(solution.field.converged(params.tol));
        assert_eq!(
            solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        assert!(divergence < 2e-10, "n={n}, div={divergence:e}");
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
        let pressure_mean = solution
            .pressure
            .iter()
            .zip(&mesh.areas)
            .map(|(p, area)| p * area)
            .sum::<f64>()
            / mesh.areas.iter().sum::<f64>();
        eprintln!(
            "layered pressure channel n={n}, dofs={}, iterations={}: L2={error:e}, div={divergence:e}, Q={outward:e}, exact Q={exact_flux:e}, mean p={pressure_mean:e}",
            system.n, solution.field.iters
        );
        assert!((inward - outward).abs() < 2e-11);
        pressure_errors.push((pressure_mean - 4.5).abs());
        if n == 16 {
            assert!(
                system.n > 1500,
                "fine channel exercises the sparse mixed solve"
            );
            assert!((pressure_mean - 4.5).abs() < 0.001);
            assert!(error < 0.003, "fine velocity error={error:e}");
            assert!(
                (outward / exact_flux - 1.0).abs() < 0.01,
                "fine flux={outward:e}"
            );
        } else {
            let frozen = system
                .assemble_stokes_mixed_variable_viscosity(
                    params,
                    &|_, _| lower_nu,
                    &|_| [0.0; 2],
                    &boundaries,
                )
                .unwrap()
                .solve()
                .unwrap();
            let frozen_error = system.velocity_error(&frozen.field.x, &exact).0;
            assert!(
                frozen_error > 5.0 * error,
                "a single reference viscosity must miss the layered flow"
            );
        }
        errors.push(error);
    }
    assert!(
        pressure_errors[1] < 0.5 * pressure_errors[0],
        "pressure datum error must decrease with spatial refinement: {pressure_errors:?}"
    );
    let refinement = errors[0] / errors[1];
    assert!(
        (3.48..4.60).contains(&refinement),
        "second-order layered channel errors={errors:?}, ratio={refinement}"
    );
}

/// The affine profile changes slope by 100:1 while its physical shear
/// traction is the same on both sides of the material interface. BDM1 can
/// represent this field exactly, exposing a wrong trace or averaging law.
#[test]
fn variable_viscosity_mixed_patch_preserves_one_sided_stress() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: f64::NAN,
        ..FluxParams::default()
    };
    let layer_nu: Vec<_> = mesh
        .tris
        .iter()
        .map(|tri| {
            let y = tri.iter().map(|&v| mesh.verts[v][1]).sum::<f64>() / 3.0;
            if y < 0.5 { 0.1 } else { 10.0 }
        })
        .collect();
    let lower_traces = Cell::new(0usize);
    let upper_traces = Cell::new(0usize);
    let viscosity = |cell: usize, point: [f64; 2]| {
        if (point[1] - 0.5).abs() < 1e-14 {
            let counter = if layer_nu[cell] < 1.0 {
                &lower_traces
            } else {
                &upper_traces
            };
            counter.set(counter.get() + 1);
        }
        layer_nu[cell]
    };
    let shear = 0.3;
    let pressure = 7.0;
    let velocity = |point: [f64; 2]| {
        let y = point[1];
        [
            0.2 + shear * (y.min(0.5) / 0.1 + (y - 0.5).max(0.0) / 10.0),
            0.0,
        ]
    };
    let traction = |_: [f64; 2], normal: [f64; 2]| {
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
            if geometry.normal[0] < -0.5 {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let solution = system
        .assemble_stokes_mixed_variable_viscosity(params, &viscosity, &|_| [0.0; 2], &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    assert!(lower_traces.get() > 0 && upper_traces.get() > 0);
    let (error, divergence) = system.velocity_error(&solution.field.x, &velocity);
    let mut stress_error = 0.0_f64;
    for (cell, nu) in layer_nu.iter().enumerate() {
        let mut shear_gradient = 0.0;
        for basis in 0..6 {
            let edge = mesh.tri_edges[cell][basis / 2].0;
            let dof = 2 * edge + basis % 2;
            let gradient = system.bases[cell].grad[basis];
            shear_gradient += solution.field.x[dof] * (gradient[0][1] + gradient[1][0]);
        }
        stress_error = stress_error.max((nu * shear_gradient - shear).abs());
    }
    let pressure_error = solution
        .pressure
        .iter()
        .fold(0.0_f64, |e, p| e.max((p - pressure).abs()));
    eprintln!(
        "layered mixed patch: L2={error:e}, div={divergence:e}, shear error={stress_error:e}, pressure error={pressure_error:e}"
    );
    assert!(error < 1e-10 && divergence < 1e-10);
    assert!(stress_error < 1e-9 && pressure_error < 1e-8);
}

#[test]
fn mixed_variable_constant_law_preserves_existing_operator_bits() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.37,
        ..FluxParams::default()
    };
    let velocity = |point: [f64; 2]| [0.3 - point[1], point[0] - 0.2];
    let force = |point: [f64; 2]| [point[0] - 0.4, -2.0 * point[1]];
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let constant = system
        .assemble_stokes_mixed(params, &force, &boundaries)
        .unwrap();
    let variable = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams {
                nu: f64::NAN,
                ..params
            },
            &|_, _| params.nu,
            &force,
            &boundaries,
        )
        .unwrap();
    assert_eq!(constant.pressure_reference(), variable.pressure_reference());
    for row in 0..system.n {
        let (columns, values) = constant.matrix().row(row);
        let (other_columns, other_values) = variable.matrix().row(row);
        assert_eq!(columns, other_columns);
        assert!(
            values
                .iter()
                .zip(other_values)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
    assert!(
        constant
            .rhs()
            .iter()
            .zip(variable.rhs())
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

#[test]
fn variable_viscosity_refuses_bad_samples_without_publishing_fallbacks() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: f64::NAN,
        ..FluxParams::default()
    };
    let velocity = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let calls = Cell::new(0usize);
    let counted = |_: usize, _: [f64; 2]| {
        calls.set(calls.get() + 1);
        0.1
    };
    assert!(matches!(
        system.assemble_stokes_mixed_variable_viscosity(
            params,
            &counted,
            &|_| [0.0; 2],
            &boundaries[1..]
        ),
        Err(StokesBoundaryError::MissingEdge { .. })
    ));
    assert_eq!(
        calls.get(),
        0,
        "boundary admission precedes constitutive queries"
    );
    for bad in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        calls.set(0);
        let invalid = |_: usize, _: [f64; 2]| {
            calls.set(calls.get() + 1);
            bad
        };
        let result = system.assemble_stokes_mixed_variable_viscosity(
            params,
            &invalid,
            &|_| [0.0; 2],
            &boundaries,
        );
        assert!(
            matches!(result, Err(StokesBoundaryError::InvalidViscosity { cell: 0, point, value })
            if point.iter().all(|x| x.is_finite()) && value.to_bits() == bad.to_bits())
        );
        assert_eq!(
            calls.get(),
            1,
            "no constitutive queries after the first refusal"
        );
    }
    // A source can admit every interior point yet reject its high-side
    // interface trace. That actual one-sided refusal must remain visible.
    let restricted = |cell: usize, point: [f64; 2]| {
        let y = mesh.tris[cell]
            .iter()
            .map(|&v| mesh.verts[v][1])
            .sum::<f64>()
            / 3.0;
        if y > 0.5 && (point[1] - 0.5).abs() < 1e-14 {
            -1.0
        } else {
            0.1
        }
    };
    assert!(matches!(
        system.assemble_stokes_mixed_variable_viscosity(params, &restricted, &|_| [0.0; 2], &boundaries),
        Err(StokesBoundaryError::InvalidViscosity { cell, point, value })
            if mesh.tris[cell].iter().map(|&v| mesh.verts[v][1]).sum::<f64>() / 3.0 > 0.5
                && (point[1] - 0.5).abs() < 1e-14 && value == -1.0
    ));
}

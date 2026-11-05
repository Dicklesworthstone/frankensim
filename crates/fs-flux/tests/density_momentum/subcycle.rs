//! G1/G3: density microsteps feed ONE original native macro momentum solve.

use super::{
    Boundary, Cell, DensityTransportError, Error, FluxParams, FluxSystem, Mesh2, StokesBoundary,
    TriMesh, controls, moments, near, numerics,
};
use fs_flux::ns::boundary::conservative::subcycle::DensityMomentumSubcycleControls;

fn options(reference: f64, depth: f64) -> DensityMomentumSubcycleControls {
    let mut momentum = controls(reference, depth);
    momentum.maximum_courant = 0.3;
    DensityMomentumSubcycleControls {
        momentum,
        maximum_substeps: 128,
    }
}

#[test]
fn macro_subcycling_preserves_high_contrast_free_stream_and_boundary_accounts() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let u = [0.4, -0.2];
    let uniform = |_: [f64; 2]| u;
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -7.0 * n);
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    let boundaries: Vec<_> = selected
        .iter()
        .map(|&edge| Boundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &uniform,
        })
        .collect();
    let previous = moments(&system, &uniform);
    let dt = 2.0;
    let depth = 2.0;
    for (low, high) in [(2.0, 8.0), (2e-8, 8e8)] {
        let cells: Vec<_> = mesh
            .centroids
            .iter()
            .map(|point| Cell {
                density_kg_m3: if point[0] < 0.5 { low } else { high },
                specific: [123.0],
            })
            .collect();
        let inflow: Vec<_> = (0..mesh.edges.len())
            .map(|edge| {
                mesh.is_boundary(edge).then_some(Cell {
                    density_kg_m3: high,
                    specific: [123.0],
                })
            })
            .collect();
        for reference in [high, 0.125 * high] {
            let calls = std::cell::Cell::new(0);
            let coupled = system
                .density_momentum_subcycled_bdf1_step(
                    numerics(),
                    &|_, _, _| 0.05 * high,
                    &|_, _, _| {
                        calls.set(calls.get() + 1);
                        [0.0; 2]
                    },
                    &boundaries,
                    &previous,
                    &previous,
                    &cells,
                    &inflow,
                    dt,
                    options(reference, depth),
                    || false,
                )
                .unwrap();
            assert!(coupled.subcycles.substeps > 1);
            assert!(coupled.subcycles.transport.maximum_courant > 1.0);
            assert!(coupled.subcycles.maximum_substep_courant <= 0.3);
            assert_eq!(calls.get(), 6 * cells.len());
            let step = coupled.momentum;
            assert_eq!(step.transport, coupled.subcycles.transport);
            let (error, divergence) = system.velocity_error(&step.flow.field.x, &uniform);
            assert!(error < 2e-10 && divergence < 2e-10);
            for &pressure in &step.pressure_pa {
                near(pressure, 7.0, 2e-9 + 4096.0 * f64::EPSILON * high / dt);
            }
            assert!(step.cells.iter().all(|cell| cell.specific == [123.0]));
            assert!(step.cells.iter().all(|cell| cell.density_kg_m3 > 0.0));
            assert!(step.cells.iter().zip(&cells).any(|(a, b)| a != b));
            let scale = step
                .transport
                .mass_after_kg
                .max(step.transport.mass_before_kg)
                .max(1.0);
            near(step.transport.relative_mass_balance_residual, 0.0, 3e-14);
            near(
                step.transport.relative_carried_balance_residual[0],
                0.0,
                3e-14,
            );
            let out = step
                .flow
                .reactions
                .transport(&selected, [0.0; 2], [0.0; 2])
                .unwrap();
            near(out.outward_volume_flux_m2_s, 0.0, 2e-14);
            for axis in 0..2 {
                near(
                    step.observations.momentum_after_kg_m_s[axis],
                    step.transport.mass_after_kg * u[axis],
                    3e-10 * scale,
                );
                near(
                    step.observations.momentum_after_kg_m_s[axis]
                        - step.observations.momentum_before_kg_m_s[axis]
                        + dt * reference * depth * out.outward_momentum.force[axis],
                    0.0,
                    3e-10 * scale,
                );
            }
            near(
                step.observations.kinetic_energy_after_j
                    - step.observations.kinetic_energy_before_j
                    + dt * reference * depth * out.outward_kinetic_energy_flux,
                0.0,
                3e-10 * scale,
            );
            assert!(step.observations.temporal_dissipation_j < 1e-18 * scale);
        }
    }
}

#[test]
fn closed_macro_identity_uses_old_inertia_and_one_averaged_native_flux() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let (a, b) = mesh.edges[internal].tris;
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let mut cells = [Cell {
        density_kg_m3: 2.0,
        specific: [100.0],
    }; 2];
    cells[b].density_kg_m3 = 8.0;
    let expected = [1277.0 / 256.0, 1283.0 / 256.0];
    let zero = |_: [f64; 2]| [0.0; 2];
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let boundaries: Vec<_> = selected
        .iter()
        .map(|&edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let native: Vec<_> = selected
        .iter()
        .map(|&edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let reference = 4.0;
    let depth = 1.7;
    let dt = 5.0;
    let mu = 0.03;
    let viscosity = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams {
                nu: mu / reference,
                sigma: numerics().penalty,
                ..FluxParams::default()
            },
            &|_, _| mu / reference,
            &zero,
            &native,
        )
        .unwrap();
    let force_calls = std::cell::Cell::new(0);
    let coupled = system
        .density_momentum_subcycled_bdf1_step(
            numerics(),
            &|cell, _, state| {
                near(state.density_kg_m3, expected[usize::from(cell == b)], 5e-14);
                mu
            },
            &|cell, _, state| {
                force_calls.set(force_calls.get() + 1);
                near(state.density_kg_m3, expected[usize::from(cell == b)], 5e-14);
                [0.0; 2]
            },
            &boundaries,
            &previous,
            &previous,
            &cells,
            &vec![None; mesh.edges.len()],
            dt,
            options(reference, depth),
            || false,
        )
        .unwrap();
    assert_eq!(coupled.subcycles.substeps, 8);
    assert_eq!(force_calls.get(), 6 * cells.len());
    let step = coupled.momentum;
    near(step.cells[a].density_kg_m3, expected[0], 5e-14);
    near(step.cells[b].density_kg_m3, expected[1], 5e-14);
    assert!(step.cells.iter().all(|state| state.specific == [100.0]));
    near(
        step.transport.mass_after_kg,
        step.transport.mass_before_kg,
        5e-14,
    );
    // Independent one-mode macro solve: each triangle has integral |phi|²=12,
    // so M=12*(rho_a+rho_b). The frozen upwind coefficient is M*9sqrt(2)*m_old
    // and the original symmetric-stress SIP stiffness is 96*mu*(4*sigma-3).
    let lambda = 96.0 * mu * (4.0 * numerics().penalty - 3.0) / 120.0;
    let old_m = previous[2 * internal + 1];
    let exact_m = old_m / (1.0 + dt * (lambda + 9.0 * 2.0f64.sqrt() * old_m));
    near(step.flow.field.x[2 * internal + 1], exact_m, 5e-13);
    let mut velocity = step.flow.field.x.clone();
    velocity[system.n_u..].fill(0.0);
    let mut action = vec![0.0; system.n];
    viscosity.matrix().spmv(&velocity, &mut action);
    let viscous_work = velocity
        .iter()
        .zip(action)
        .map(|(u, au)| u * au)
        .sum::<f64>()
        * dt
        * reference
        * depth;
    let obs = step.observations;
    assert!(obs.temporal_dissipation_j > 0.0 && obs.interior_upwind_dissipation_j > 0.0);
    assert!(viscous_work > 0.0);
    near(
        obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
            + obs.temporal_dissipation_j
            + obs.interior_upwind_dissipation_j
            + viscous_work,
        0.0,
        3e-12,
    );
    let loads = step
        .flow
        .reactions
        .resultant(&selected, [0.0; 2], &[Some(0.0)])
        .unwrap();
    for axis in 0..2 {
        near(
            obs.momentum_after_kg_m_s[axis] - obs.momentum_before_kg_m_s[axis]
                + dt * reference * depth * loads.stress.force[axis],
            0.0,
            3e-12,
        );
    }
}

#[test]
fn one_substep_matches_original_native_matrix_rhs_solution_and_observations_exactly() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let cells = [
        Cell {
            density_kg_m3: 2.0,
            specific: [3.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [-1.0],
        },
    ];
    let inflow = vec![None; mesh.edges.len()];
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let settings = options(4.0, 1.7);
    let original = system
        .assemble_density_momentum_bdf1(
            numerics(),
            &|_, _, _| 0.03,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            0.1,
            settings.momentum,
            || false,
        )
        .unwrap();
    let subcycled = system
        .assemble_density_momentum_subcycled_bdf1(
            numerics(),
            &|_, _, _| 0.03,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            0.1,
            settings,
            || false,
        )
        .unwrap();
    assert_eq!(subcycled.subcycles().substeps, 1);
    assert_eq!(original.transported_cells(), subcycled.transported_cells());
    assert_eq!(
        original.momentum().matrix().to_dense(),
        subcycled.momentum().matrix().to_dense()
    );
    assert_eq!(original.momentum().rhs(), subcycled.momentum().rhs());
    let expected = original.solve(|| false).unwrap();
    let actual = subcycled.solve(|| false).unwrap().momentum;
    assert_eq!(actual.cells, expected.cells);
    assert_eq!(actual.flow.field.x, expected.flow.field.x);
    assert_eq!(actual.pressure_pa, expected.pressure_pa);
    assert_eq!(actual.observations, expected.observations);
}

#[test]
fn macro_work_refusal_and_cancellation_leave_original_state_and_no_candidate() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let cells = [Cell::<0> {
        density_kg_m3: 2.0,
        specific: [],
    }; 2];
    let retained = previous.clone();
    let inflow = vec![None; mesh.edges.len()];
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let calls = std::cell::Cell::new(0);
    let force = |_, _, _: &Cell<0>| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let mut settings = options(4.0, 1.0);
    settings.maximum_substeps = 4;
    let refused = system.assemble_density_momentum_subcycled_bdf1(
        numerics(),
        &|_, _, _| 0.03,
        &force,
        &boundaries,
        &previous,
        &previous,
        &cells,
        &inflow,
        5.0,
        settings,
        || false,
    );
    assert!(matches!(
        refused,
        Err(Error::Transport(DensityTransportError::SubstepLimit {
            maximum: 4
        }))
    ));
    assert_eq!(calls.get(), 0);
    settings.maximum_substeps = 8;
    let mut polls = 0;
    let accepted = system
        .assemble_density_momentum_subcycled_bdf1(
            numerics(),
            &|_, _, _| 0.03,
            &force,
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            5.0,
            settings,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(matches!(accepted.solve(|| true), Err(Error::Cancelled)));
    let mut seen = 0;
    let cancelled = system.assemble_density_momentum_subcycled_bdf1(
        numerics(),
        &|_, _, _| 0.03,
        &force,
        &boundaries,
        &previous,
        &previous,
        &cells,
        &inflow,
        5.0,
        settings,
        || {
            seen += 1;
            seen == polls
        },
    );
    assert!(matches!(
        cancelled,
        Err(Error::Cancelled | Error::Transport(DensityTransportError::Cancelled))
    ));
    assert_eq!(previous, retained);
    assert_eq!(
        cells,
        [Cell {
            density_kg_m3: 2.0,
            specific: []
        }; 2]
    );
}

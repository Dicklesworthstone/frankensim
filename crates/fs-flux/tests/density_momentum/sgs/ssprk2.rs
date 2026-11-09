use super::*;
use fs_flux::ns::boundary::conservative::sgs::DensityMomentumSgsStep;
use fs_flux::ns::boundary::conservative::ssprk2::DensityMomentumSspRk2Controls;
use fs_flux::ns::transport::density::DensityReconstruction;

fn options(reference: f64, depth: f64) -> DensityMomentumSspRk2Controls {
    let mut momentum = controls(reference, depth);
    momentum.maximum_courant = 0.3;
    DensityMomentumSspRk2Controls {
        momentum,
        maximum_substeps: 128,
    }
}

fn viscosity(state: &Cell<3>, rate: f64) -> f64 {
    0.02 + 2e-5 * state.specific[0] + 0.01 * state.specific[1] + 0.005 * rate / (1.0 + rate)
}

fn assert_closed_kinetic_balance<const N: usize>(
    system: &FluxSystem<'_>,
    solved: &DensityMomentumSgsStep<N>,
    dt: f64,
    reference: f64,
    depth: f64,
) {
    let zero = |_: [f64; 2]| [0.0; 2];
    let native: Vec<_> = (0..system.mesh.edges.len())
        .filter(|&edge| system.mesh.is_boundary(edge))
        .map(|edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let viscous = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams {
                sigma: numerics().penalty,
                ..FluxParams::default()
            },
            &|cell, _| solved.cells[cell].effective_dynamic_viscosity_pa_s / reference,
            &zero,
            &native,
        )
        .unwrap();
    let mut velocity = solved.step.flow.field.x.clone();
    velocity[system.n_u..].fill(0.0);
    let mut action = vec![0.0; system.n];
    viscous.matrix().spmv(&velocity, &mut action);
    let viscous_work = velocity
        .iter()
        .zip(action)
        .map(|(u, au)| u * au)
        .sum::<f64>()
        * reference
        * depth
        * dt;
    let obs = solved.step.observations;
    assert!(viscous_work > 0.0 && obs.temporal_dissipation_j > 0.0);
    near(
        obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
            + obs.temporal_dissipation_j
            + obs.interior_upwind_dissipation_j
            + viscous_work,
        0.0,
        3e-10,
    );
}

#[test]
fn rk_sgs_uses_analytic_stage_endpoint_for_all_material_and_closure_samples() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let (a, b) = mesh.edges[internal].tris;
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let mut original = [Cell {
        density_kg_m3: 2.0,
        specific: [100.0, 0.1, 0.9],
    }; 2];
    original[b] = Cell {
        density_kg_m3: 8.0,
        specific: [500.0, 0.8, 0.2],
    };
    // The closed two-cell exchange has the same linear operator for density
    // and each extensive inventory. Eight SSPRK2 intervals multiply its
    // antisymmetric mode by (1 - 1/2 + (1/2)^2/2)^8 = (5/8)^8.
    let decay = (5.0_f64 / 8.0).powi(8);
    let mut expected = original;
    for (cell, sign) in [(a, -1.0), (b, 1.0)] {
        let density = 5.0 + sign * 3.0 * decay;
        expected[cell].density_kg_m3 = density;
        for slot in 0..3 {
            let qa = original[a].density_kg_m3 * original[a].specific[slot];
            let qb = original[b].density_kg_m3 * original[b].specific[slot];
            expected[cell].specific[slot] =
                (0.5 * (qa + qb) + sign * 0.5 * (qb - qa) * decay) / density;
        }
    }
    let inflow = vec![None; mesh.edges.len()];
    let widths = [0.5; 2];
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    for model in models() {
        let force_calls = Counter::new(0);
        let source_calls = Counter::new(0);
        let boundary_calls = Counter::new(0);
        let wall = |_: [f64; 2]| {
            boundary_calls.set(boundary_calls.get() + 1);
            [0.0; 2]
        };
        let boundaries: Vec<_> = selected
            .iter()
            .map(|&edge| Boundary::Velocity { edge, value: &wall })
            .collect();
        let samples = RefCell::new([(0.0, 0.0, Cell::default()); 2]);
        let check_state = |cell: usize, state: &Cell<3>| {
            near(state.density_kg_m3, expected[cell].density_kg_m3, 8e-14);
            for slot in 0..3 {
                near(state.specific[slot], expected[cell].specific[slot], 2e-12);
            }
        };
        let coupled = system
            .density_momentum_ssprk2_bdf1_sgs_step(
                numerics(),
                SgsConfig {
                    model,
                    filter_widths_m: &widths,
                },
                &|cell, point, state, rate| {
                    assert_eq!(point, mesh.centroids[cell]);
                    check_state(cell, state);
                    source_calls.set(source_calls.get() + 1);
                    let mu = viscosity(state, rate);
                    samples.borrow_mut()[cell] = (rate, mu, *state);
                    mu
                },
                &|cell, _, state| {
                    check_state(cell, state);
                    force_calls.set(force_calls.get() + 1);
                    [0.0; 2]
                },
                &boundaries,
                &previous,
                &previous,
                &original,
                &inflow,
                5.0,
                options(4.0, 1.7),
                nonlinear(),
                || false,
            )
            .unwrap();
        assert_eq!(coupled.density_stages.substeps, 8);
        assert_eq!(coupled.density_stages.stage_evaluations, 16);
        assert!(coupled.density_stages.maximum_stage_courant <= 0.3);
        let solved = coupled.momentum;
        assert_eq!(solved.step.transport, coupled.density_stages.transport);
        assert_eq!(force_calls.get(), 6 * original.len());
        assert_eq!(boundary_calls.get(), 2 * selected.len());
        assert!(source_calls.get() >= original.len() * (solved.iterations + 1));
        assert!(
            solved.iterations > 1 && solved.nonlinear_residual <= nonlinear().relative_tolerance
        );
        let mut powers = [0.0; 2];
        for (cell, response) in solved.cells.iter().enumerate() {
            let state = &solved.step.cells[cell];
            check_state(cell, state);
            near(state.specific[1] + state.specific[2], 1.0, 2e-15);
            let k = system
                .try_cell_kinematics(&solved.step.flow.field.x, cell)
                .unwrap();
            let retained = samples.borrow()[cell];
            assert_eq!(
                (
                    response.strain_rate_per_s,
                    response.molecular_dynamic_viscosity_pa_s,
                    *state
                ),
                retained
            );
            let nu_t = eddy(model, widths[cell], k.gradient());
            near(response.eddy_viscosity_m2_s, nu_t, 1e-14);
            near(
                response.modeled_dynamic_viscosity_pa_s,
                state.density_kg_m3 * nu_t,
                1e-14,
            );
            let rate2 = k.shear_rate().powi(2);
            near(response.molecular_heating_w_m3, retained.1 * rate2, 1e-14);
            near(
                response.modeled_transfer_w_m3,
                state.density_kg_m3 * nu_t * rate2,
                1e-14,
            );
            powers[0] += mesh.areas[cell] * 1.7 * response.molecular_heating_w_m3;
            powers[1] += mesh.areas[cell] * 1.7 * response.modeled_transfer_w_m3;
        }
        near(solved.molecular_heating_w, powers[0], 1e-14);
        near(solved.modeled_transfer_w, powers[1], 1e-14);
        assert!(solved.molecular_heating_w > 0.0 && solved.modeled_transfer_w > 0.0);
        near(
            solved.step.transport.relative_mass_balance_residual,
            0.0,
            3e-14,
        );
        for residual in solved.step.transport.relative_carried_balance_residual {
            near(residual, 0.0, 3e-14);
        }
        assert_closed_kinetic_balance(&system, &solved, 5.0, 4.0, 1.7);
    }
}

#[test]
fn reconstructed_rk_inventories_keep_native_flux_energy_and_exact_zero_sgs_limit() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh)
        .with_density_reconstruction(DensityReconstruction::BoundedLinearInventories);
    let previous = moments(&system, &|[x, y]| {
        [
            x * (1.0 - x) * (1.0 - 2.0 * y),
            -(1.0 - 2.0 * x) * y * (1.0 - y),
        ]
    });
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|[x, y]| Cell {
            density_kg_m3: 2.0 + x + 2.0 * y,
            specific: [500.0 + 200.0 * x + 100.0 * y, 0.1 + 0.6 * x, 0.9 - 0.6 * x],
        })
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let widths = vec![0.5; cells.len()];
    let settings = options(4.0, 1.7);
    let baseline = system
        .density_momentum_ssprk2_bdf1_step(
            numerics(),
            &|_, _, state| viscosity(state, 0.0),
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            0.7,
            settings,
            || false,
        )
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
        SgsModel::Wale { cw: 0.7 },
    ] {
        let active = model.coefficient() > 0.0;
        let coupled = system
            .density_momentum_ssprk2_bdf1_sgs_step(
                numerics(),
                SgsConfig {
                    model,
                    filter_widths_m: &widths,
                },
                &|_, _, state, _| viscosity(state, 0.0),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &previous,
                &previous,
                &cells,
                &inflow,
                0.7,
                settings,
                StokesRheologyControls {
                    relaxation: if active { nonlinear().relaxation } else { 1.0 },
                    ..nonlinear()
                },
                || false,
            )
            .unwrap();
        assert_eq!(coupled.density_stages, baseline.density_stages);
        assert!(coupled.density_stages.substeps > 1);
        let solved = coupled.momentum;
        assert_eq!(solved.step.cells, baseline.momentum.cells);
        assert_eq!(solved.step.transport, baseline.momentum.transport);
        for cell in &solved.step.cells {
            assert!(cell.density_kg_m3 > 0.0);
            assert!((500.0..=800.0).contains(&cell.specific[0]));
            assert!(
                (0.0..=1.0).contains(&cell.specific[1]) && (0.0..=1.0).contains(&cell.specific[2])
            );
            near(cell.specific[1] + cell.specific[2], 1.0, 4e-15);
        }
        if active {
            assert!(solved.modeled_transfer_w > 0.0);
            assert!(
                solved.step.observations.kinetic_energy_after_j
                    < baseline.momentum.observations.kinetic_energy_after_j
            );
        } else {
            assert_eq!(solved.step.flow.field.x, baseline.momentum.flow.field.x);
            assert_eq!(solved.step.pressure_pa, baseline.momentum.pressure_pa);
            assert_eq!(solved.modeled_transfer_w, 0.0);
        }
        assert_closed_kinetic_balance(&system, &solved, 0.7, 4.0, 1.7);
    }
}

#[test]
fn rk_sgs_stage_exhaustion_late_source_refusal_and_cancellation_are_atomic() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut old = vec![0.0; system.n_u];
    old[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let original = [
        Cell::<1> {
            density_kg_m3: 2.0,
            specific: [100.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [500.0],
        },
    ];
    let retained_cells = original;
    let retained_velocity = old.clone();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let widths = [0.5; 2];
    let config = SgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.7 },
        filter_widths_m: &widths,
    };
    let source_calls = Counter::new(0);
    let force_calls = Counter::new(0);
    let molecular = |_, _, _: &Cell<1>, _| {
        source_calls.set(source_calls.get() + 1);
        0.03
    };
    let force = |_, _, _: &Cell<1>| {
        force_calls.set(force_calls.get() + 1);
        [0.0; 2]
    };
    let mut settings = options(4.0, 1.0);
    settings.maximum_substeps = 4;
    let refused = system.density_momentum_ssprk2_bdf1_sgs_step(
        numerics(),
        config,
        &molecular,
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        5.0,
        settings,
        nonlinear(),
        || false,
    );
    assert!(matches!(
        refused,
        Err(SgsError::Momentum(Error::Transport(
            DensityTransportError::SubstepLimit { maximum: 4 }
        )))
    ));
    assert_eq!((source_calls.get(), force_calls.get()), (0, 0));
    settings.maximum_substeps = 8;
    let refused = system.density_momentum_ssprk2_bdf1_sgs_step(
        numerics(),
        config,
        &|_, _, _, _| {
            source_calls.set(source_calls.get() + 1);
            if source_calls.get() > original.len() {
                f64::NAN
            } else {
                0.03
            }
        },
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        5.0,
        settings,
        nonlinear(),
        || false,
    );
    assert!(matches!(refused, Err(SgsError::Momentum(Error::Flow(_)))));
    assert!(source_calls.get() > original.len());
    assert_eq!(force_calls.get(), 6 * original.len());
    source_calls.set(0);
    let cancelled = system.density_momentum_ssprk2_bdf1_sgs_step(
        numerics(),
        config,
        &molecular,
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        5.0,
        settings,
        nonlinear(),
        || source_calls.get() > 0,
    );
    assert!(matches!(
        cancelled,
        Err(SgsError::Momentum(Error::Cancelled))
    ));
    assert_eq!(source_calls.get(), 1);
    let exhausted = system.density_momentum_ssprk2_bdf1_sgs_step(
        numerics(),
        config,
        &molecular,
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        5.0,
        settings,
        StokesRheologyControls {
            maximum_iterations: 1,
            ..nonlinear()
        },
        || false,
    );
    assert!(matches!(exhausted, Err(SgsError::Momentum(Error::Flow(_)))));
    assert_eq!(original, retained_cells);
    assert_eq!(old, retained_velocity);
}

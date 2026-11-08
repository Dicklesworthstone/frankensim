use super::*;
use fs_flux::ns::boundary::StokesRheologyControls;
use fs_flux::ns::boundary::conservative::sgs::{
    DensityMomentumSgsConfig as SgsConfig, DensityMomentumSgsError as SgsError,
};
use fs_flux::ns::sgs::SgsModel;
use std::cell::{Cell as Counter, RefCell};

fn nonlinear() -> StokesRheologyControls {
    StokesRheologyControls {
        maximum_iterations: 120,
        relative_tolerance: 2e-10,
        relaxation: 0.65,
    }
}

fn models() -> [SgsModel; 3] {
    [
        SgsModel::Smagorinsky { cs: 0.7 },
        SgsModel::Wale { cw: 0.7 },
        SgsModel::Vreman { cv: 0.15 },
    ]
}

fn eddy(model: SgsModel, width: f64, gradient: [[f64; 2]; 2]) -> f64 {
    fs_lbm::sgs3::eddy_viscosity(
        model,
        width,
        [
            [gradient[0][0], gradient[0][1], 0.0],
            [gradient[1][0], gradient[1][1], 0.0],
            [0.0; 3],
        ],
    )
    .unwrap()
}

#[test]
fn variable_density_sgs_closes_native_kinetic_identity_with_one_transport_and_endpoint_source() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.05;
    let original = [
        Cell::<1> {
            density_kg_m3: 2.0,
            specific: [3.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [19.0],
        },
    ];
    let retained_history = previous.clone();
    let inflow = vec![None; mesh.edges.len()];
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let native: Vec<_> = selected
        .iter()
        .map(|&edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let dt = 0.08;
    let depth = 1.7;
    let widths = vec![0.5; original.len()];
    for model in models() {
        let mut reference_result = None;
        for reference in [1.0, 9.0] {
            let boundary_calls = Counter::new(0);
            let wall = |_: [f64; 2]| {
                boundary_calls.set(boundary_calls.get() + 1);
                [0.0; 2]
            };
            let boundaries: Vec<_> = selected
                .iter()
                .map(|&edge| Boundary::Velocity { edge, value: &wall })
                .collect();
            let force_calls = Counter::new(0);
            let source_calls = Counter::new(0);
            let samples = RefCell::new(vec![(0.0, 0.0, Cell::default()); original.len()]);
            let molecular = |cell, point, transported: &Cell<1>, rate: f64| {
                assert_eq!(point, mesh.centroids[cell]);
                source_calls.set(source_calls.get() + 1);
                let mu = 0.02 + 0.003 * rate / (1.0 + rate);
                samples.borrow_mut()[cell] = (rate, mu, *transported);
                mu
            };
            let solved = system
                .density_momentum_bdf1_sgs_step(
                    numerics(),
                    SgsConfig {
                        model,
                        filter_widths_m: &widths,
                    },
                    &molecular,
                    &|_, _, _| {
                        force_calls.set(force_calls.get() + 1);
                        [0.0; 2]
                    },
                    &boundaries,
                    &previous,
                    &previous,
                    &original,
                    &inflow,
                    dt,
                    controls(reference, depth),
                    nonlinear(),
                    || false,
                )
                .unwrap();
            assert_eq!(force_calls.get(), 6 * original.len());
            assert_eq!(boundary_calls.get(), 2 * selected.len());
            assert!(source_calls.get() >= original.len() * (solved.iterations + 1));
            assert!(
                solved.iterations > 1
                    && solved.nonlinear_residual <= nonlinear().relative_tolerance
            );
            let step = &solved.step;
            assert!(step.cells[0].density_kg_m3 > original[0].density_kg_m3);
            assert!(step.cells[1].density_kg_m3 < original[1].density_kg_m3);
            let baseline = system
                .density_momentum_bdf1_step(
                    numerics(),
                    &|_, _, _| 0.02,
                    &|_, _, _| [0.0; 2],
                    &selected
                        .iter()
                        .map(|&edge| Boundary::Velocity { edge, value: &zero })
                        .collect::<Vec<_>>(),
                    &previous,
                    &previous,
                    &original,
                    &inflow,
                    dt,
                    controls(reference, depth),
                    || false,
                )
                .unwrap();
            assert_eq!(step.cells, baseline.cells);
            assert_eq!(step.transport, baseline.transport);
            assert!(
                step.observations.kinetic_energy_after_j
                    < baseline.observations.kinetic_energy_after_j
            );
            let mut molecular_power = 0.0;
            let mut modeled_power = 0.0;
            for (cell, output) in solved.cells.iter().enumerate() {
                let k = system
                    .try_cell_kinematics(&step.flow.field.x, cell)
                    .unwrap();
                let rate = k.shear_rate();
                let nu_t = eddy(model, widths[cell], k.gradient());
                let retained = samples.borrow()[cell];
                near(output.strain_rate_per_s, rate, 1e-14);
                assert_eq!(output.strain_rate_per_s, retained.0);
                assert_eq!(output.molecular_dynamic_viscosity_pa_s, retained.1);
                assert_eq!(step.cells[cell], retained.2);
                near(output.eddy_viscosity_m2_s, nu_t, 1e-14);
                near(
                    output.modeled_dynamic_viscosity_pa_s,
                    step.cells[cell].density_kg_m3 * nu_t,
                    1e-14,
                );
                near(
                    output.molecular_heating_w_m3,
                    retained.1 * rate * rate,
                    1e-14,
                );
                near(
                    output.modeled_transfer_w_m3,
                    step.cells[cell].density_kg_m3 * nu_t * rate * rate,
                    1e-14,
                );
                molecular_power += mesh.areas[cell] * depth * output.molecular_heating_w_m3;
                modeled_power += mesh.areas[cell] * depth * output.modeled_transfer_w_m3;
            }
            near(solved.molecular_heating_w, molecular_power, 1e-14);
            near(solved.modeled_transfer_w, modeled_power, 1e-14);
            assert!(modeled_power > 0.0);
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
            let mut velocity = step.flow.field.x.clone();
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
            let obs = step.observations;
            near(
                obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
                    + obs.temporal_dissipation_j
                    + obs.interior_upwind_dissipation_j
                    + viscous_work,
                0.0,
                2e-9,
            );
            let load = step
                .flow
                .reactions
                .resultant(&selected, [0.0; 2], &[Some(0.0)])
                .unwrap();
            for axis in 0..2 {
                near(
                    obs.momentum_after_kg_m_s[axis] - obs.momentum_before_kg_m_s[axis]
                        + dt * reference * depth * load.stress.force[axis],
                    0.0,
                    2e-9,
                );
            }
            if let Some(expected) = &reference_result {
                for (&a, &b) in step.flow.field.x[..system.n_u].iter().zip(expected) {
                    near(a, b, 2e-10);
                }
            } else {
                reference_result = Some(step.flow.field.x[..system.n_u].to_vec());
            }
        }
    }
    assert_eq!(previous, retained_history);
}

#[test]
fn full_gradient_sgs_and_actual_molecular_rheology_close_physical_traction_at_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 7.0;
    let original = vec![
        Cell::<0> {
            density_kg_m3: rho,
            specific: []
        };
        mesh.tris.len()
    ];
    let inflow = vec![None; mesh.edges.len()];
    let advector = vec![0.0; system.n_u];
    let width = 0.4;
    let widths = vec![width; original.len()];
    let dt = 0.2;
    let depth = 1.3;
    let gradients = [[[0.3, 0.2], [0.0, -0.3]], [[0.3, 0.6], [-0.4, -0.3]]];
    for model in models() {
        let mut expected_eddies = Vec::new();
        for gradient in gradients {
            let exact = |p: [f64; 2]| {
                [
                    0.4 + gradient[0][0] * p[0] + gradient[0][1] * p[1],
                    gradient[1][0] * p[0] + gradient[1][1] * p[1],
                ]
            };
            let old = moments(&system, &|p| exact(p).map(|value| 0.7 * value));
            let rate = 0.4_f64.sqrt();
            let mu = 0.1 + 0.2 / (1.0 + rate);
            let nu_t = eddy(model, width, gradient);
            expected_eddies.push(nu_t);
            let total_mu = mu + rho * nu_t;
            let pressure = |p: [f64; 2]| 5.0 + 2.0 * p[0] - 3.0 * p[1];
            let traction = |p: [f64; 2], n: [f64; 2]| {
                [
                    total_mu
                        * (2.0 * gradient[0][0] * n[0] + (gradient[0][1] + gradient[1][0]) * n[1])
                        - pressure(p) * n[0],
                    total_mu
                        * ((gradient[0][1] + gradient[1][0]) * n[0] + 2.0 * gradient[1][1] * n[1])
                        - pressure(p) * n[1],
                ]
            };
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&e| mesh.is_boundary(e))
                .map(|edge| {
                    if mesh.edges[edge].normal[0].abs() > 0.5 {
                        Boundary::Open {
                            edge,
                            traction: &traction,
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    } else {
                        Boundary::Velocity {
                            edge,
                            value: &exact,
                        }
                    }
                })
                .collect();
            let solved = system
                .density_momentum_bdf1_sgs_step(
                    numerics(),
                    SgsConfig {
                        model,
                        filter_widths_m: &widths,
                    },
                    &|_, _, state, gamma| {
                        assert_eq!(state.density_kg_m3, rho);
                        0.1 + 0.2 / (1.0 + gamma)
                    },
                    &|_, point, state| {
                        let u = exact(point);
                        [
                            state.density_kg_m3 * 0.3 / dt * u[0] + 2.0,
                            state.density_kg_m3 * 0.3 / dt * u[1] - 3.0,
                        ]
                    },
                    &boundaries,
                    &advector,
                    &old,
                    &original,
                    &inflow,
                    dt,
                    controls(3.0, depth),
                    nonlinear(),
                    || false,
                )
                .unwrap();
            let (error, divergence) = system.velocity_error(&solved.step.flow.field.x, &exact);
            assert!(
                error < 3e-8 && divergence < 3e-9,
                "{model:?}: {error:e}, {divergence:e}"
            );
            assert_eq!(solved.step.cells, original);
            for (cell, output) in solved.cells.iter().enumerate() {
                near(output.strain_rate_per_s, rate, 2e-7);
                near(output.molecular_dynamic_viscosity_pa_s, mu, 2e-8);
                near(output.eddy_viscosity_m2_s, nu_t, 2e-8);
                near(
                    solved.step.pressure_pa[cell],
                    pressure(mesh.centroids[cell]),
                    2e-7,
                );
            }
            near(solved.molecular_heating_w, depth * mu * 0.4, 3e-8);
            near(solved.modeled_transfer_w, depth * rho * nu_t * 0.4, 3e-8);
        }
        if matches!(model, SgsModel::Smagorinsky { .. }) {
            near(expected_eddies[0], expected_eddies[1], 1e-15);
        } else {
            assert!((expected_eddies[0] - expected_eddies[1]).abs() > 1e-5);
        }
    }
}

#[test]
fn zero_sgs_recovers_native_density_step_and_refusals_preserve_both_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut old = vec![0.0; system.n_u];
    old[2 * internal + 1] = 0.04;
    let original = [
        Cell::<0> {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 9.0,
            specific: [],
        },
    ];
    let retained_velocity = old.clone();
    let retained_cells = original;
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let baseline = system
        .density_momentum_bdf1_step(
            numerics(),
            &|_, _, _| 0.07,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &old,
            &old,
            &original,
            &inflow,
            0.05,
            controls(3.0, 1.0),
            || false,
        )
        .unwrap();
    let widths = vec![0.5; original.len()];
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let solved = system
            .density_momentum_bdf1_sgs_step(
                numerics(),
                SgsConfig {
                    model,
                    filter_widths_m: &widths,
                },
                &|_, _, _, _| 0.07,
                &|_, _, _| [0.0; 2],
                &boundaries,
                &old,
                &old,
                &original,
                &inflow,
                0.05,
                controls(3.0, 1.0),
                StokesRheologyControls {
                    relaxation: 1.0,
                    ..nonlinear()
                },
                || false,
            )
            .unwrap();
        assert_eq!(solved.step.cells, baseline.cells);
        assert_eq!(solved.step.transport, baseline.transport);
        assert_eq!(solved.step.flow.field.x, baseline.flow.field.x);
        assert_eq!(solved.modeled_transfer_w, 0.0);
        assert!(solved.molecular_heating_w > 0.0);
    }
    let source_calls = Counter::new(0);
    let force_calls = Counter::new(0);
    let force = |_, _, _: &Cell<0>| {
        force_calls.set(force_calls.get() + 1);
        [0.0; 2]
    };
    let invalid = system.density_momentum_bdf1_sgs_step(
        numerics(),
        SgsConfig {
            model: SgsModel::Smagorinsky { cs: f64::NAN },
            filter_widths_m: &widths,
        },
        &|_, _, _, _| {
            source_calls.set(source_calls.get() + 1);
            0.07
        },
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        0.05,
        controls(3.0, 1.0),
        nonlinear(),
        || false,
    );
    assert!(matches!(invalid, Err(SgsError::Subgrid(_))));
    assert_eq!((source_calls.get(), force_calls.get()), (0, 0));
    let config = SgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.7 },
        filter_widths_m: &widths,
    };
    let invalid = system.density_momentum_bdf1_sgs_step(
        numerics(),
        config,
        &|cell, _, _, _| if cell == 1 { f64::NAN } else { 0.07 },
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        0.05,
        controls(3.0, 1.0),
        nonlinear(),
        || false,
    );
    assert!(matches!(invalid, Err(SgsError::Momentum(Error::Flow(_)))));
    let exhausted = system.density_momentum_bdf1_sgs_step(
        numerics(),
        config,
        &|_, _, _, _| 0.07,
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        0.05,
        controls(3.0, 1.0),
        StokesRheologyControls {
            maximum_iterations: 1,
            ..nonlinear()
        },
        || false,
    );
    assert!(matches!(exhausted, Err(SgsError::Momentum(Error::Flow(_)))));
    let cancelled = system.density_momentum_bdf1_sgs_step(
        numerics(),
        config,
        &|_, _, _, _| {
            source_calls.set(source_calls.get() + 1);
            0.07
        },
        &force,
        &boundaries,
        &old,
        &old,
        &original,
        &inflow,
        0.05,
        controls(3.0, 1.0),
        nonlinear(),
        || source_calls.get() > 0,
    );
    assert!(matches!(
        cancelled,
        Err(SgsError::Momentum(Error::Cancelled))
    ));
    assert_eq!(old, retained_velocity);
    assert_eq!(original, retained_cells);
}

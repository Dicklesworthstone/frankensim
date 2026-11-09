//! Native bounded density work at the accepted momentum/pressure/rheology endpoint.

use super::*;
use fs_flux::material::{
    MaterialStokesWindkesselHistory as History, MaterialStokesWindkesselPort as Port,
};
use fs_flux::ns::boundary::conservative::coupled::windkessel::DensityMomentumWindkesselError as RcrError;
use fs_flux::ns::boundary::conservative::coupled::{
    CoupledDensityMomentumControls as CoupledControls, CoupledDensityMomentumStep,
    CoupledDensityTransportReport as Work,
};
use fs_flux::ns::boundary::conservative::ssprk2::DensityMomentumSspRk2Controls;
use fs_flux::ns::boundary::conservative::subcycle::DensityMomentumSubcycleControls;
use fs_flux::ns::transport::density::DensityReconstruction;
use fs_flux::ns::transport::density::ssprk2::DensitySspRk2Controls;
use fs_flux::ns::transport::density::subcycle::DensitySubcycleControls;
use fs_material::fluid_transport::rheology::CarreauYasuda;
use std::cell::RefCell;

fn options(rk: bool, momentum: Controls, maximum_substeps: usize) -> CoupledControls {
    if rk {
        DensityMomentumSspRk2Controls {
            momentum,
            maximum_substeps,
        }
        .into()
    } else {
        DensityMomentumSubcycleControls {
            momentum,
            maximum_substeps,
        }
        .into()
    }
}

fn replay<const N: usize>(
    system: &FluxSystem<'_>,
    result: &CoupledDensityMomentumStep<N>,
    old: &[Cell<N>],
    inflow: &[Option<Cell<N>>],
    control: Controls,
    maximum_substeps: usize,
) -> usize {
    let step = &result.step;
    let plan = system
        .prepare_density_transport(
            &step.flow.field.x[..system.n_u],
            control.depth_m,
            control.relative_volume_divergence_tolerance,
            control.maximum_transport_storage_bytes,
            || false,
        )
        .unwrap();
    let mut cells = vec![Cell::default(); old.len()];
    let substeps = match result.density_work {
        Work::SubcycledEuler(actual) => {
            let mut work = plan
                .subcycle_workspace(control.maximum_transport_storage_bytes)
                .unwrap();
            let expected = work
                .advance_into(
                    old,
                    inflow,
                    step.transport.timestep_s,
                    DensitySubcycleControls {
                        maximum_courant: control.maximum_courant,
                        maximum_substeps,
                    },
                    &mut cells,
                    || false,
                )
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(actual.transport, step.transport);
            assert!(actual.maximum_substep_courant <= control.maximum_courant);
            actual.substeps
        }
        Work::SspRk2(actual) => {
            let mut work = plan
                .ssprk2_workspace(control.maximum_transport_storage_bytes)
                .unwrap();
            let expected = work
                .advance_into(
                    old,
                    inflow,
                    step.transport.timestep_s,
                    DensitySspRk2Controls {
                        maximum_courant: control.maximum_courant,
                        maximum_substeps,
                    },
                    &mut cells,
                    || false,
                )
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(actual.transport, step.transport);
            assert_eq!(actual.stage_evaluations, 2 * actual.substeps);
            assert!(actual.maximum_stage_courant <= control.maximum_courant);
            actual.substeps
        }
        Work::SingleStep(_) => panic!("bounded density work was requested"),
    };
    assert_eq!(cells, step.cells);
    near(step.transport.mass_balance_residual_kg, 0.0, 2e-10);
    for residual in step.transport.carried_balance_residual {
        near(residual, 0.0, 3e-8);
    }
    substeps
}

#[test]
fn macro_cfl_above_one_accelerates_and_transports_once_with_original_native_stages() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let endpoint = |_: [f64; 2]| [0.4, -0.2];
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -7.0 * n);
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &endpoint,
        })
        .collect();
    let old: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| Cell {
            density_kg_m3: if p[0] < 0.5 { 2.0 } else { 80.0 },
            specific: if p[0] < 0.5 {
                [100.0, 0.2]
            } else {
                [300.0, 0.8]
            },
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 8.0,
                specific: [220.0, 0.4],
            })
        })
        .collect();
    let dt = 2.0;
    for reconstruction in [
        DensityReconstruction::DonorCell,
        DensityReconstruction::BoundedLinear,
        DensityReconstruction::BoundedLinearInventories,
    ] {
        let system = FluxSystem::new(&mesh).with_density_reconstruction(reconstruction);
        let previous = vec![0.0; system.n_u];
        let force = |cell: usize, _: [f64; 2], _: &Cell<2>| {
            [
                old[cell].density_kg_m3 * 0.4 / dt,
                -old[cell].density_kg_m3 * 0.2 / dt,
            ]
        };
        for rk in [false, true] {
            for reference in [10.0, 80.0] {
                let control = Controls {
                    maximum_courant: 0.3,
                    ..controls(reference, 1.7)
                };
                let result = system
                    .density_momentum_coupled_bdf1_step(
                        numerics(),
                        &|_, _, _| 0.2,
                        &force,
                        &boundaries,
                        &previous,
                        &old,
                        &inflow,
                        dt,
                        options(rk, control, 128),
                        nonlinear(),
                        || false,
                    )
                    .unwrap();
                let (error, divergence) =
                    system.velocity_error(&result.step.flow.field.x, &endpoint);
                assert!(
                    error < 3e-9 && divergence < 3e-9,
                    "{error:e}, {divergence:e}"
                );
                assert!(result.step.transport.maximum_courant > 1.0);
                assert!(replay(&system, &result, &old, &inflow, control, 128) > 1);
                assert!(result.step.cells.iter().zip(&old).any(|(a, b)| a != b));
                assert!(result.step.cells.iter().all(|s| s.density_kg_m3 > 0.0
                    && (100.0 - 1e-10..=300.0 + 1e-10).contains(&s.specific[0])
                    && (0.2 - 1e-12..=0.8 + 1e-12).contains(&s.specific[1])));
                for p in &result.step.pressure_pa {
                    near(*p, 7.0, 2e-7);
                }
                // Ordinary controls retain the existing strict macro-CFL refusal.
                assert!(matches!(
                    system.density_momentum_coupled_bdf1_step(
                        numerics(),
                        &|_, _, _| 0.2,
                        &force,
                        &boundaries,
                        &previous,
                        &old,
                        &inflow,
                        dt,
                        control,
                        nonlinear(),
                        || false,
                    ),
                    Err(Error::Transport(DensityTransportError::Courant { .. }))
                ));
            }
        }
        assert!(previous.iter().all(|u| *u == 0.0));
    }
    assert_eq!(old[0].specific, [100.0, 0.2]);
}

#[test]
fn accumulated_flux_closes_macro_kinetic_work_at_the_density_and_rheology_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.05;
    let old = [
        Cell {
            density_kg_m3: 2.0,
            specific: [3.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [19.0],
        },
    ];
    let inflow = vec![None; mesh.edges.len()];
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let native: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let dt = 5.0;
    let control = Controls {
        maximum_courant: 0.05,
        ..controls(7.0, 1.7)
    };
    for rk in [false, true] {
        let samples = RefCell::new(vec![(Cell::default(), 0.0, 0.0); old.len()]);
        let source = |cell, _, state: &Cell<1>, rate| {
            let curve = CarreauYasuda::new(
                0.004 + 0.0002 * state.density_kg_m3,
                0.001,
                4.0,
                0.5,
                2.0,
                [0.0, 100.0],
            )?;
            let mu = curve.viscosity_pa_s(rate)?;
            samples.borrow_mut()[cell] = (*state, rate, mu);
            Ok(mu)
        };
        let solved = system
            .density_momentum_coupled_bdf1_rheology_step(
                numerics(),
                &source,
                &|_, _, _| [0.0; 2],
                &boundaries,
                &previous,
                &old,
                &inflow,
                dt,
                options(rk, control, 128),
                nonlinear(),
                || false,
            )
            .unwrap();
        assert!(replay(&system, &solved.step, &old, &inflow, control, 128) > 1);
        let step = &solved.step.step;
        assert!(solved.viscous_heating_w > 0.0);
        for (cell, response) in solved.cells.iter().enumerate() {
            assert_eq!(
                samples.borrow()[cell],
                (
                    step.cells[cell],
                    response.strain_rate_per_s,
                    response.dynamic_viscosity_pa_s
                )
            );
            near(
                response.strain_rate_per_s,
                system
                    .try_cell_kinematics(&step.flow.field.x, cell)
                    .unwrap()
                    .shear_rate(),
                1e-13,
            );
        }
        let viscous = system
            .assemble_stokes_mixed_variable_viscosity(
                FluxParams {
                    sigma: numerics().penalty,
                    ..FluxParams::default()
                },
                &|cell, _| {
                    solved.cells[cell].dynamic_viscosity_pa_s / control.reference_density_kg_m3
                },
                &zero,
                &native,
            )
            .unwrap();
        let mut velocity = step.flow.field.x.clone();
        velocity[system.n_u..].fill(0.0);
        let mut action = vec![0.0; system.n];
        viscous.matrix().spmv(&velocity, &mut action);
        let work = velocity
            .iter()
            .zip(action)
            .map(|(u, au)| u * au)
            .sum::<f64>()
            * dt
            * control.depth_m
            * control.reference_density_kg_m3;
        let obs = step.observations;
        assert!(
            obs.temporal_dissipation_j > 0.0
                && obs.interior_upwind_dissipation_j > 0.0
                && work > 0.0
        );
        near(
            obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
                + obs.temporal_dissipation_j
                + obs.interior_upwind_dissipation_j
                + work,
            0.0,
            2e-9,
        );
    }
    assert_eq!(previous[2 * internal + 1], 0.05);
}

fn port_curve(state: &Cell<1>, rate: f64) -> Result<f64, fs_material::MaterialError> {
    CarreauYasuda::new(
        0.6 + 0.02 * state.density_kg_m3 + 0.002 * state.specific[0],
        0.05,
        4.0,
        0.45,
        2.0,
        [0.0, 100.0],
    )?
    .viscosity_pa_s(rate)
}

#[test]
fn bounded_density_work_reaches_sheared_pressure_ports_and_preserves_work_refusal() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let west: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].mid[0] == 0.0)
        .collect();
    let east: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].mid[0] == 1.0)
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let old: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| {
            if p[0] < 0.5 {
                Cell {
                    density_kg_m3: 2.0,
                    specific: [3.0],
                }
            } else {
                Cell {
                    density_kg_m3: 8.0,
                    specific: [19.0],
                }
            }
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 6.0,
                specific: [9.0],
            })
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let control = Controls {
        maximum_courant: 0.02,
        ..controls(13.0, 1.7)
    };
    for rk in [false, true] {
        let run = |maximum_substeps| {
            system.density_momentum_coupled_bdf1_rheology_flow_control_step(
                numerics(),
                &|_, _, state, rate| port_curve(state, rate),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &[&west],
                &[-0.15],
                2e-11,
                &previous,
                &old,
                &inflow,
                0.25,
                options(rk, control, maximum_substeps),
                nonlinear(),
                || false,
            )
        };
        let solved = run(128).unwrap();
        assert!(replay(&system, &solved.step.step, &old, &inflow, control, 128) > 1);
        near(solved.outward_volume_rates_m2_s[0], -0.15, 2e-11);
        assert!(solved.maximum_flow_residual_m2_s <= 2e-11);
        assert!(solved.step.viscous_heating_w > 0.0);
        assert!(solved.step.cells.iter().any(|c| c.strain_rate_per_s > 0.05));
        let reaction = solved
            .step
            .step
            .step
            .flow
            .reactions
            .resultant(&west, [0.0; 2], &[None])
            .unwrap();
        near(
            reaction.imposed.force[0] * control.reference_density_kg_m3,
            -solved.pressure_offsets_pa[0],
            2e-8,
        );
        assert!(matches!(
            run(1),
            Err(Error::Transport(DensityTransportError::SubstepLimit {
                maximum: 1
            }))
        ));
    }
    assert!(previous.iter().all(|v| *v == 0.0));
    assert_eq!(old[0].specific, [3.0]);
}

#[test]
fn bounded_density_stages_do_not_substep_external_rcr_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let west: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].mid[0] == 0.0)
        .collect();
    let east: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].mid[0] == 1.0)
        .collect();
    let ports = [&west[..], &east[..]];
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let old: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| {
            if p[0] < 0.5 {
                Cell {
                    density_kg_m3: 2.0,
                    specific: [3.0],
                }
            } else {
                Cell {
                    density_kg_m3: 8.0,
                    specific: [19.0],
                }
            }
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 6.0,
                specific: [9.0],
            })
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let old_x = [1.0, 0.2];
    let loads = [
        Port {
            reservoir_pressure_pa: 12.0,
            proximal_resistance_pa_s_per_m2: 1.0,
            distal_resistance_pa_s_per_m2: 20.0,
            compliance_m2_per_pa: 0.03,
            inertance_pa_s2_per_m2: 0.2,
        },
        Port {
            reservoir_pressure_pa: 0.0,
            proximal_resistance_pa_s_per_m2: 2.0,
            distal_resistance_pa_s_per_m2: 30.0,
            compliance_m2_per_pa: 0.02,
            inertance_pa_s2_per_m2: 0.3,
        },
    ];
    let dt = 0.08;
    let control = Controls {
        maximum_courant: 0.01,
        ..controls(13.0, 1.7)
    };
    for rk in [false, true] {
        let run = |maximum_substeps| {
            system.density_momentum_coupled_bdf1_rheology_windkessel_step(
                numerics(),
                &|_, _, state, rate| port_curve(state, rate),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &ports,
                &loads,
                History {
                    velocity: &previous,
                    compliance_pressures_pa: &old_x,
                },
                &old,
                &inflow,
                dt,
                options(rk, control, maximum_substeps),
                nonlinear(),
                2e-9,
                2e-11,
                || false,
            )
        };
        let solved = run(128).unwrap();
        assert!(replay(&system, &solved.flow.step, &old, &inflow, control, 128) > 1);
        assert!(solved.flow.viscous_heating_w > 0.0);
        for i in 0..2 {
            let load = loads[i];
            let q = solved.outward_volume_rates_m2_s[i];
            let x = (load.compliance_m2_per_pa * old_x[i] / dt + q)
                / (load.compliance_m2_per_pa / dt + 1.0 / load.distal_resistance_pa_s_per_m2);
            let p = load.reservoir_pressure_pa
                + load.proximal_resistance_pa_s_per_m2 * q
                + load.inertance_pa_s2_per_m2 * q / dt
                + x;
            let balance = solved.balances[i];
            near(balance.compliance_pressure_pa, x, 2e-9);
            near(solved.pressure_offsets_pa[i], p, 2e-8);
            near(balance.port_power_w_per_m, p * q, 2e-8);
            near(balance.power_balance_residual_w_per_m, 0.0, 2e-9);
            near(balance.compliance_flow_residual_m2_s, 0.0, 2e-11);
            let defect = (load.compliance_m2_per_pa * (x - old_x[i]).powi(2)
                + load.inertance_pa_s2_per_m2 * q * q)
                / (2.0 * dt);
            near(balance.bdf_energy_defect_w_per_m, defect, 2e-10);
        }
        assert!(matches!(
            run(1),
            Err(RcrError::Momentum(Error::Transport(
                DensityTransportError::SubstepLimit { maximum: 1 }
            )))
        ));
    }
    assert!(previous.iter().all(|v| *v == 0.0));
    assert_eq!(old_x, [1.0, 0.2]);
}

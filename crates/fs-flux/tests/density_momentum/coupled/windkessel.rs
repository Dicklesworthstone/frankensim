//! G1/G3: finite external storage and changing-density native momentum share one endpoint.

use super::*;
use fs_flux::material::{
    MaterialFluxError, MaterialStokesWindkesselHistory as History,
    MaterialStokesWindkesselPort as Port,
};
use fs_flux::ns::boundary::conservative::coupled::windkessel::DensityMomentumWindkesselError as RcrError;

fn ports(mesh: &TriMesh) -> (Vec<usize>, Vec<usize>) {
    let side = |x: f64| {
        (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].mid[0] == x)
            .collect()
    };
    (side(0.0), side(1.0))
}

// Independently assemble the ORIGINAL simultaneous momentum and two storage
// equations for a unit-length/height uniform plug, in physical per-depth units.
// Unknowns are (U, x_west, x_east), not the reduced port response parameters.
fn plug_endpoint(rho: f64, loads: &[Port; 2], old_u: f64, old_x: [f64; 2], dt: f64) -> [f64; 3] {
    let [west, east] = loads;
    let inertia = rho + west.inertance_pa_s2_per_m2 + east.inertance_pa_s2_per_m2;
    let matrix = [
        inertia / dt + west.proximal_resistance_pa_s_per_m2 + east.proximal_resistance_pa_s_per_m2,
        -1.0,
        1.0,
        1.0,
        west.compliance_m2_per_pa / dt + 1.0 / west.distal_resistance_pa_s_per_m2,
        0.0,
        -1.0,
        0.0,
        east.compliance_m2_per_pa / dt + 1.0 / east.distal_resistance_pa_s_per_m2,
    ];
    let mut state = [
        west.reservoir_pressure_pa - east.reservoir_pressure_pa + inertia * old_u / dt,
        west.compliance_m2_per_pa * old_x[0] / dt,
        east.compliance_m2_per_pa * old_x[1] / dt,
    ];
    fs_la::factor::lu(&matrix, 3).unwrap().solve(&mut state);
    state
}

#[test]
fn two_loaded_ports_advance_heterogeneous_density_and_both_storage_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let (west, east) = ports(&mesh);
    let open: Vec<_> = west.iter().chain(&east).copied().collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -7.0 * n);
    let initial: Vec<Cell<2>> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[0] < 0.5 { 2.0 } else { 80.0 },
            specific: if point[0] < 0.5 {
                [100.0, 0.2]
            } else {
                [300.0, 0.8]
            },
        })
        .collect();
    let initial_velocity = moments(&system, &|_| [0.04, 0.0]);
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| {
            open.contains(&edge).then_some(Cell {
                density_kg_m3: 8.0,
                specific: [220.0, 0.4],
            })
        })
        .collect();
    let loads = [
        Port {
            reservoir_pressure_pa: 90.0,
            proximal_resistance_pa_s_per_m2: 3.0,
            distal_resistance_pa_s_per_m2: 40.0,
            compliance_m2_per_pa: 0.03,
            inertance_pa_s2_per_m2: 0.4,
        },
        Port {
            reservoir_pressure_pa: 10.0,
            proximal_resistance_pa_s_per_m2: 5.0,
            distal_resistance_pa_s_per_m2: 60.0,
            compliance_m2_per_pa: 0.02,
            inertance_pa_s2_per_m2: 0.7,
        },
    ];
    for (reference, depth) in [(1.0, 1.7), (80.0, 1.7), (80.0, 0.5)] {
        let mut old_cells = initial.clone();
        let mut old_velocity = initial_velocity.clone();
        let mut old_x = [12.0, 4.0];
        let mut old_u = 0.04;
        for dt in [0.05, 0.03] {
            let mean_density: f64 = old_cells
                .iter()
                .zip(&mesh.areas)
                .map(|(cell, area)| cell.density_kg_m3 * area)
                .sum();
            let expected = plug_endpoint(mean_density, &loads, old_u, old_x, dt);
            let endpoint = |_: [f64; 2]| [expected[0], 0.0];
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&edge| mesh.is_boundary(edge))
                .map(|edge| {
                    if open.contains(&edge) {
                        Boundary::OpenReservoir {
                            edge,
                            traction: &traction,
                            velocity: &endpoint,
                        }
                    } else {
                        Boundary::NavierSlip {
                            edge,
                            value: &zero,
                            tangential_resistance_pa_s_per_m: 0.0,
                        }
                    }
                })
                .collect();
            // Zero spatial mean: compensate local old-density variation so
            // uniform acceleration sees the independently computed actual mass.
            let force = |cell: usize, _: [f64; 2], _: &Cell<2>| {
                [
                    (old_cells[cell].density_kg_m3 - mean_density) * (expected[0] - old_u) / dt,
                    0.0,
                ]
            };
            let control = controls(reference, depth);
            let result = system
                .density_momentum_coupled_bdf1_windkessel_step(
                    numerics(),
                    &|_, _, _| 0.2,
                    &force,
                    &boundaries,
                    &[&west, &east],
                    &loads,
                    History {
                        velocity: &old_velocity,
                        compliance_pressures_pa: &old_x,
                    },
                    &old_cells,
                    &inflow,
                    dt,
                    control,
                    nonlinear(),
                    2e-9,
                    2e-11,
                    || false,
                )
                .unwrap();
            let coupled = &result.flow;
            let step = &coupled.step;
            assert!(coupled.iterations >= 2);
            assert!(coupled.nonlinear_residual <= nonlinear().relative_tolerance);
            let (error, divergence) = system.velocity_error(&step.flow.field.x, &endpoint);
            assert!(
                error < 3e-10 && divergence < 3e-10,
                "{error:e}, {divergence:e}"
            );
            let q = [-expected[0], expected[0]];
            let old_q = [-old_u, old_u];
            let pressure: [f64; 2] = core::array::from_fn(|i| {
                loads[i].reservoir_pressure_pa
                    + loads[i].proximal_resistance_pa_s_per_m2 * q[i]
                    + loads[i].inertance_pa_s2_per_m2 * (q[i] - old_q[i]) / dt
                    + expected[i + 1]
            });
            for i in 0..2 {
                near(result.outward_volume_rates_m2_s[i], q[i], 2e-10);
                near(result.pressure_offsets_pa[i], pressure[i], 2e-7);
                let b = result.balances[i];
                let l = loads[i];
                let x = expected[i + 1];
                near(b.compliance_pressure_pa, x, 2e-9);
                near(
                    b.compliance_energy_j_per_m,
                    0.5 * l.compliance_m2_per_pa * x * x,
                    2e-9,
                );
                near(
                    b.inertance_energy_j_per_m,
                    0.5 * l.inertance_pa_s2_per_m2 * q[i] * q[i],
                    2e-10,
                );
                near(
                    b.proximal_dissipation_w_per_m,
                    l.proximal_resistance_pa_s_per_m2 * q[i] * q[i],
                    2e-9,
                );
                near(
                    b.distal_dissipation_w_per_m,
                    x * x / l.distal_resistance_pa_s_per_m2,
                    2e-9,
                );
                let old_energy = 0.5
                    * (l.compliance_m2_per_pa * old_x[i] * old_x[i]
                        + l.inertance_pa_s2_per_m2 * old_q[i] * old_q[i]);
                near(
                    b.stored_energy_change_rate_w_per_m,
                    (b.compliance_energy_j_per_m + b.inertance_energy_j_per_m - old_energy) / dt,
                    2e-9,
                );
                let temporal = (l.compliance_m2_per_pa * (x - old_x[i]).powi(2)
                    + l.inertance_pa_s2_per_m2 * (q[i] - old_q[i]).powi(2))
                    / (2.0 * dt);
                near(b.bdf_energy_defect_w_per_m, temporal, 2e-9);
                near(b.port_power_w_per_m, pressure[i] * q[i], 2e-8);
                near(b.power_balance_residual_w_per_m, 0.0, 2e-9);
                near(b.pressure_law_residual_pa, 0.0, 2e-8);
                near(b.compliance_flow_residual_m2_s, 0.0, 2e-11);
                let selected = if i == 0 { &west } else { &east };
                let reaction = step
                    .flow
                    .reactions
                    .resultant(selected, [0.0; 2], &[None])
                    .unwrap();
                let external_power = reaction
                    .imposed
                    .rigid_power(endpoint([0.0; 2]), 0.0)
                    .unwrap()
                    * reference
                    * depth;
                // Original common 7 Pa traction plus the added RCR pressure.
                near(
                    external_power,
                    (7.0 * q[i] + b.port_power_w_per_m) * depth,
                    5e-8,
                );
            }
            for (&p, center) in step.pressure_pa.iter().zip(&mesh.centroids) {
                near(
                    p,
                    7.0 + pressure[0] * (1.0 - center[0]) + pressure[1] * center[0],
                    2e-7,
                );
            }
            let plan = system
                .prepare_density_transport(
                    &step.flow.field.x[..system.n_u],
                    depth,
                    control.relative_volume_divergence_tolerance,
                    control.maximum_transport_storage_bytes,
                    || false,
                )
                .unwrap();
            let mut workspace = plan
                .workspace(control.maximum_transport_storage_bytes)
                .unwrap();
            let mut transported = vec![Cell::<2>::default(); old_cells.len()];
            let accounting = plan
                .advance_into(
                    &old_cells,
                    &inflow,
                    dt,
                    control.maximum_courant,
                    &mut transported,
                    &mut workspace,
                    || false,
                )
                .unwrap();
            assert_eq!(step.cells, transported);
            assert_eq!(step.transport, accounting);
            assert!(step.cells.iter().zip(&old_cells).any(|(a, b)| a != b));
            near(accounting.mass_balance_residual_kg, 0.0, 2e-11);
            for residual in accounting.carried_balance_residual {
                near(residual, 0.0, 4e-9);
            }
            let actual_mass: f64 = step
                .cells
                .iter()
                .zip(&mesh.areas)
                .map(|(cell, area)| cell.density_kg_m3 * area * depth)
                .sum();
            near(
                step.observations.kinetic_energy_after_j,
                0.5 * actual_mass * expected[0] * expected[0],
                2e-9,
            );
            assert!((actual_mass - reference * depth).abs() > 1.0);
            old_x = [
                result.balances[0].compliance_pressure_pa,
                result.balances[1].compliance_pressure_pa,
            ];
            old_u = expected[0];
            old_velocity = step.flow.field.x[..system.n_u].to_vec();
            old_cells = step.cells.clone();
        }
    }
    assert_eq!(initial[0].specific, [100.0, 0.2]);
}

#[test]
fn precharged_compliance_drives_reverse_flow_without_heating_carried_inventory() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let (west, east) = ports(&mesh);
    let loads = [Port {
        reservoir_pressure_pa: 0.0,
        proximal_resistance_pa_s_per_m2: 2.0,
        distal_resistance_pa_s_per_m2: 40.0,
        compliance_m2_per_pa: 0.03,
        inertance_pa_s2_per_m2: 0.3,
    }; 2];
    let old_x = [-8.0, 8.0];
    let old_velocity = vec![0.0; system.n_u];
    let cells = vec![
        Cell::<1> {
            density_kg_m3: 10.0,
            specific: [300.0]
        };
        mesh.tris.len()
    ];
    let dt = 0.04;
    let expected = plug_endpoint(10.0, &loads, 0.0, old_x, dt);
    assert!(expected[0] < 0.0);
    let endpoint = |_: [f64; 2]| [expected[0], 0.0];
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &endpoint,
                }
            } else {
                Boundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            }
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| mesh.is_boundary(edge).then_some(cells[0]))
        .collect();
    let run = |load: &[Port],
               history: &[f64],
               donor: &[Option<Cell<1>>],
               control,
               cancel: &mut dyn FnMut() -> bool| {
        system.density_momentum_coupled_bdf1_windkessel_step(
            numerics(),
            &|_, _, _| 0.2,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &[&west, &east],
            load,
            History {
                velocity: &old_velocity,
                compliance_pressures_pa: history,
            },
            &cells,
            donor,
            dt,
            control,
            nonlinear(),
            2e-9,
            2e-11,
            cancel,
        )
    };
    let result = run(&loads, &old_x, &inflow, controls(37.0, 2.5), &mut || false).unwrap();
    let (error, divergence) = system.velocity_error(&result.flow.step.flow.field.x, &endpoint);
    assert!(
        error < 3e-10 && divergence < 3e-10,
        "{error:e}, {divergence:e}"
    );
    for i in 0..2 {
        let b = result.balances[i];
        near(b.compliance_pressure_pa, expected[i + 1], 2e-9);
        assert!(
            b.port_power_w_per_m < 0.0,
            "precharged storage returns work"
        );
        assert!(b.stored_energy_change_rate_w_per_m < 0.0);
        assert!(b.proximal_dissipation_w_per_m > 0.0);
        assert!(b.distal_dissipation_w_per_m > 0.0);
        near(b.reservoir_power_w_per_m, 0.0, 0.0);
        near(
            b.port_power_w_per_m
                - b.stored_energy_change_rate_w_per_m
                - b.proximal_dissipation_w_per_m
                - b.distal_dissipation_w_per_m
                - b.bdf_energy_defect_w_per_m,
            0.0,
            2e-9,
        );
    }
    assert!(result.flow.step.observations.kinetic_energy_after_j > 0.0);
    for cell in &result.flow.step.cells {
        near(cell.specific[0], 300.0, 2e-10);
    }
    let only_west: Vec<_> = (0..mesh.edges.len())
        .map(|edge| west.contains(&edge).then_some(cells[0]))
        .collect();
    assert!(matches!(
        run(&loads, &old_x, &only_west, controls(37.0, 2.5), &mut || {
            false
        }),
        Err(RcrError::Momentum(Error::Transport(_)))
    ));
    let mut invalid = loads;
    invalid[0].compliance_m2_per_pa = 0.0;
    assert!(matches!(
        run(&invalid, &old_x, &inflow, controls(37.0, 2.5), &mut || {
            false
        }),
        Err(RcrError::Hydraulic(MaterialFluxError::Invalid(_)))
    ));
    assert!(run(&loads, &[0.0], &inflow, controls(37.0, 2.5), &mut || false).is_err());
    assert!(matches!(
        run(
            &loads,
            &old_x,
            &inflow,
            Controls {
                maximum_courant: 1e-6,
                ..controls(37.0, 2.5)
            },
            &mut || false
        ),
        Err(RcrError::Momentum(Error::Transport(
            DensityTransportError::Courant { .. }
        )))
    ));
    let mut polls = 0;
    run(&loads, &old_x, &inflow, controls(37.0, 2.5), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        run(&loads, &old_x, &inflow, controls(37.0, 2.5), &mut || {
            seen += 1;
            seen == polls
        }),
        Err(RcrError::Momentum(Error::Cancelled))
    ));
    assert!(old_velocity.iter().all(|u| *u == 0.0));
    assert_eq!(old_x, [-8.0, 8.0]);
    assert_eq!(cells[0].specific, [300.0]);
}

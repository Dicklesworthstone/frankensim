//! G1/G3/G4: native hydraulic continuation, accepted work and window rollback.

use super::*;
use crate::laminar_enthalpy::adaptive::{LaminarAdaptiveControls, LaminarAdaptiveFailure};
use crate::laminar_enthalpy::windkessel::adaptive::LaminarWindkesselAdaptiveControls;

fn adaptive_controls() -> LaminarWindkesselAdaptiveControls {
    LaminarWindkesselAdaptiveControls {
        coupled: LaminarAdaptiveControls {
            step: controls(),
            initial_step_s: 0.037,
            minimum_step_s: 1e-6,
            maximum_step_s: 0.037,
            maximum_native_steps: 300,
            maximum_accepted_steps: 100,
            maximum_endpoint_evaluations: 10_000,
            energy_absolute_tolerance_j: 1e-5,
            temperature_absolute_tolerance_k: 1e-5,
            velocity_absolute_tolerance_m_s: 1e-5,
            transfer_absolute_tolerance_j: 1e-5,
            relative_tolerance: 0.0,
            window_tolerance_j: 1e-7,
        },
        pressure_law_tolerance_pa: 1e-9,
        flow_law_tolerance_m2_s: 1e-10,
        pressure_absolute_tolerance_pa: 1e-4,
        flow_absolute_tolerance_m2_s: 1e-5,
        transport_absolute_tolerance_m3: 1e-5,
        storage_absolute_tolerance_j: 1e-4,
        transfer_absolute_tolerance_j: 1e-4,
        window_tolerance_j: 1e-7,
    }
}

#[test]
fn external_compliance_alone_requests_refinement_and_keeps_paired_be_history() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let original_cells = old.cells().clone();
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics::default(),
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let previous = moments(&mesh, &zero);
        let previous_x = [2.0; 2];
        let load = MaterialStokesWindkesselPort {
            reservoir_pressure_pa: 0.0,
            proximal_resistance_pa_s_per_m2: 0.3,
            distal_resistance_pa_s_per_m2: 2.1,
            compliance_m2_per_pa: 0.6,
            inertance_pa_s2_per_m2: 0.4,
        };
        let loads = [load; 2];
        let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
        let mut ports = [Vec::new(), Vec::new()];
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                let n = mesh.edges[edge].normal;
                if n[0].abs() > 0.5 {
                    ports[usize::from(n[0] > 0.0)].push(edge);
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &traction,
                        velocity: &zero,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity { edge, value: &zero }
                }
            })
            .collect();
        let port_edges = [ports[0].as_slice(), ports[1].as_slice()];
        let input = LaminarWindkesselInput {
            fluid: LaminarEnthalpyInput {
                history: LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                force_density: &zero,
                boundaries: &boundaries,
                slip_heat_routes: &[],
            },
            ports: &port_edges,
            loads: &loads,
            previous_compliance_pressures_pa: &previous_x,
        };
        let target: f64 = 0.111;
        let exact_x = previous_x[0]
            * (-target / (load.compliance_m2_per_pa * load.distal_resistance_pa_s_per_m2)).exp();
        let mut global_errors = Vec::new();
        for tight in [false, true] {
            let mut controls = adaptive_controls();
            if !tight {
                controls.pressure_absolute_tolerance_pa = 1.0;
                controls.storage_absolute_tolerance_j = 1.0;
                controls.transfer_absolute_tolerance_j = 1.0;
                controls.coupled.maximum_native_steps = 9;
                controls.coupled.maximum_accepted_steps = 6;
            }
            let solved = owner
                .advance_windkessel_adaptive(
                    cx,
                    &old,
                    target,
                    input,
                    controls,
                    |_, t, _| Ok(source(t, 0.4)),
                    |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
                )
                .unwrap();
            assert_eq!(solved.thermal().time_s(), target);
            if tight {
                assert!(solved.report().work.rejected_comparisons > 0);
            } else {
                assert_eq!(solved.report().work.native_steps, 9);
                assert_eq!(solved.report().steps.len(), 6);
            }
            let (mut x, mut loss, mut damping) = (previous_x[0], 0.0, 0.0);
            for step in &solved.report().steps {
                let dt = step.end_s - step.start_s;
                let next = x
                    / (1.0 + dt / (load.compliance_m2_per_pa * load.distal_resistance_pa_s_per_m2));
                let circuit = &step.report.exchange.circuit;
                for i in 0..2 {
                    near(circuit.balances[i].compliance_pressure_pa, next, 2e-8);
                    near(circuit.pressure_offsets_pa[i], next, 2e-8);
                    near(circuit.outward_volume_rates_m2_s[i], 0.0, 1e-9);
                }
                loss += next * next / load.distal_resistance_pa_s_per_m2 * dt * DEPTH;
                damping += 0.5 * load.compliance_m2_per_pa * (next - x).powi(2) * DEPTH;
                x = next;
            }
            let stored = 0.5 * load.compliance_m2_per_pa * (x * x - previous_x[0].powi(2)) * DEPTH;
            near(stored + loss + damping, 0.0, 1e-12);
            for balance in solved.report().hydraulic.as_ref().unwrap() {
                near(balance.compliance_energy_change_j, stored, 3e-8);
                near(balance.inertance_energy_change_j, 0.0, 1e-12);
                near(balance.distal_dissipation_j, loss, 3e-8);
                near(balance.proximal_dissipation_j, 0.0, 1e-12);
                near(balance.port_work_j, 0.0, 1e-9);
                near(balance.reservoir_work_j, 0.0, 1e-12);
                near(balance.bdf_energy_defect_j, damping, 3e-8);
                near(balance.balance_residual_j, 0.0, 1e-8);
                near(balance.transported_volume_m3, 0.0, 1e-9);
                near(balance.absolute_transported_volume_m3, 0.0, 1e-9);
            }
            for estimate in &solved.report().estimates {
                assert!(estimate.hydraulic.unwrap().maximum_error_ratio <= 1.0);
                assert!(estimate.velocity_error_ratio < 1e-3);
                assert!(estimate.energy_error_ratio < 1e-3);
                assert!(estimate.temperature_error_ratio < 1e-3);
            }
            near(solved.report().stored_energy_change_j, 0.0, 1e-10);
            near(solved.report().bulk_viscous_heat_j, 0.0, 1e-10);
            near(solved.report().slip_heat_j, 0.0, 1e-10);
            near(solved.report().external_stream_enthalpy_j, 0.0, 1e-10);
            for &id in &cells {
                near(solved.thermal().temperature_k(id).unwrap(), 300.0, 1e-10);
            }
            assert!(solved.report().suggested_step_s.is_finite());
            global_errors.push((x - exact_x).abs());
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(previous_x, [2.0; 2]);
        }
        assert!(global_errors[1] < global_errors[0]);
    });
}

#[test]
fn equilibrium_signed_throughput_work_and_late_window_refusal_remain_atomic() {
    for q in [-0.2, 0.2] {
        with_cx(|cx, gate| {
            let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
            let system = FluxSystem::new(&mesh);
            let law = Caloric { datum: 0.0 };
            let (old, cells) = inventory(&mesh, &law, &[300.0; 2], true);
            let original_cells = old.cells().clone();
            let material = reference(0.4);
            let owner = LaminarEnthalpy::new(
                cx,
                &system,
                &material,
                MaterialFluxNumerics::default(),
                &old,
                map(&cells),
                DEPTH,
                2 * mesh.edges.len(),
            )
            .unwrap();
            let velocity = |_| [q, 0.0];
            let previous = moments(&mesh, &velocity);
            let original_velocity = previous.clone();
            let previous_x = [-2.1 * q, 3.0 * q];
            let loads = [
                MaterialStokesWindkesselPort {
                    reservoir_pressure_pa: 6.1 * q,
                    proximal_resistance_pa_s_per_m2: 0.3,
                    distal_resistance_pa_s_per_m2: 2.1,
                    compliance_m2_per_pa: 0.6,
                    inertance_pa_s2_per_m2: 0.3,
                },
                MaterialStokesWindkesselPort {
                    reservoir_pressure_pa: 0.0,
                    proximal_resistance_pa_s_per_m2: 0.7,
                    distal_resistance_pa_s_per_m2: 3.0,
                    compliance_m2_per_pa: 0.8,
                    inertance_pa_s2_per_m2: 0.4,
                },
            ];
            let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
            let mut ports = [Vec::new(), Vec::new()];
            let mut routes = Vec::new();
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&e| mesh.is_boundary(e))
                .map(|edge| {
                    let n = mesh.edges[edge].normal;
                    if n[0].abs() > 0.5 {
                        ports[usize::from(n[0] > 0.0)].push(edge);
                        MaterialNavierStokesBoundary::OpenReservoir {
                            edge,
                            traction: &traction,
                            velocity: &velocity,
                        }
                    } else {
                        routes.push(StokesSlipHeatRoute {
                            edge,
                            wall_cell: WALL,
                            fluid_heat_fraction: 1.0,
                        });
                        MaterialNavierStokesBoundary::NavierSlip {
                            edge,
                            value: &zero,
                            tangential_resistance_pa_s_per_m: 0.0,
                        }
                    }
                })
                .collect();
            let port_edges = [ports[0].as_slice(), ports[1].as_slice()];
            let input = LaminarWindkesselInput {
                fluid: LaminarEnthalpyInput {
                    history: LaminarEnthalpyHistory {
                        velocity_m_s: &previous,
                        time_s: 0.0,
                    },
                    force_density: &zero,
                    boundaries: &boundaries,
                    slip_heat_routes: &routes,
                },
                ports: &port_edges,
                loads: &loads,
                previous_compliance_pressures_pa: &previous_x,
            };
            let target = 0.111;
            // Each load is exactly at x=Rd*Q. Both added port pressures are
            // 3.7*q, so the original fluid remains a uniform perfect-slip plug.
            let solved = owner
                .advance_windkessel_adaptive(
                    cx,
                    &old,
                    target,
                    input,
                    adaptive_controls(),
                    |_, t, _| Ok(source(t, 0.4)),
                    |_, _, _| Ok::<_, &'static str>(exchange(&mesh, law.h(300.0))),
                )
                .unwrap();
            assert_eq!(solved.thermal().time_s(), target);
            assert_eq!(solved.report().steps.len(), 6);
            assert_eq!(solved.report().work.native_steps, 9);
            let (error, divergence) =
                system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
            assert!(error < 2e-8 && divergence < 2e-8);
            for (i, balance) in solved
                .report()
                .hydraulic
                .as_ref()
                .unwrap()
                .iter()
                .enumerate()
            {
                let signed_q = if i == 0 { -q } else { q };
                near(
                    balance.transported_volume_m3,
                    signed_q * DEPTH * target,
                    2e-9,
                );
                near(
                    balance.absolute_transported_volume_m3,
                    q.abs() * DEPTH * target,
                    2e-9,
                );
                near(
                    balance.proximal_dissipation_j,
                    loads[i].proximal_resistance_pa_s_per_m2 * q * q * DEPTH * target,
                    2e-9,
                );
                near(
                    balance.distal_dissipation_j,
                    loads[i].distal_resistance_pa_s_per_m2 * q * q * DEPTH * target,
                    2e-9,
                );
                near(
                    balance.reservoir_work_j,
                    loads[i].reservoir_pressure_pa * signed_q * DEPTH * target,
                    2e-9,
                );
                near(
                    balance.port_work_j,
                    3.7 * q * signed_q * DEPTH * target,
                    2e-9,
                );
                near(balance.compliance_energy_change_j, 0.0, 2e-9);
                near(balance.inertance_energy_change_j, 0.0, 2e-9);
                near(balance.bdf_energy_defect_j, 0.0, 2e-9);
                near(balance.balance_residual_j, 0.0, 2e-9);
            }
            near(solved.report().stored_energy_change_j, 0.0, 1e-10);
            near(solved.report().bulk_viscous_heat_j, 0.0, 1e-10);
            near(solved.report().thermal_balance_residual_j, 0.0, 1e-9);
            for mode in 0..3 {
                let mut controls = adaptive_controls();
                if mode == 1 {
                    controls.coupled.maximum_native_steps = 3;
                }
                let attempted = owner.advance_windkessel_adaptive(
                    cx,
                    &old,
                    target,
                    input,
                    controls,
                    |_, t, _| Ok(source(t, 0.4)),
                    |_, candidate, _| {
                        if candidate.time_s() > 0.04 {
                            if mode == 0 {
                                return Err("late original heat source refused");
                            }
                            if mode == 2 {
                                gate.request();
                            }
                        }
                        Ok(exchange(&mesh, law.h(300.0)))
                    },
                );
                let error = match attempted {
                    Ok(_) => panic!("failed window returned a candidate"),
                    Err(error) => error,
                };
                match mode {
                    0 => {
                        assert!(matches!(
                            error.failure,
                            LaminarAdaptiveFailure::Native {
                                error: LaminarEnthalpyError::Exchange(_),
                                ..
                            }
                        ));
                        assert_eq!(error.work.native_steps, 4);
                    }
                    1 => {
                        assert!(matches!(error.failure, LaminarAdaptiveFailure::WorkBudget));
                        assert_eq!(error.work.native_steps, 3);
                    }
                    _ => {
                        assert!(matches!(
                            error.failure,
                            LaminarAdaptiveFailure::Native {
                                error: LaminarEnthalpyError::Cancelled,
                                ..
                            }
                        ));
                        assert_eq!(error.work.native_steps, 4);
                    }
                }
                assert_eq!(old.cells(), &original_cells);
                assert_eq!(old.time_s(), 0.0);
                assert_eq!(previous, original_velocity);
                assert_eq!(previous_x, [-2.1 * q, 3.0 * q]);
            }
        });
    }
}

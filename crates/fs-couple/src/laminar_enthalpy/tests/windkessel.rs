//! G1/G3/G4: original thermal rheology, native storage and atomic paired history.

use super::*;
use crate::laminar_enthalpy::windkessel::{LaminarWindkesselControls, LaminarWindkesselInput};
use fs_flux::MaterialStokesWindkesselPort;

fn policy() -> LaminarWindkesselControls {
    LaminarWindkesselControls {
        coupled: controls(),
        absolute_pressure_tolerance_pa: 1e-9,
        absolute_flow_tolerance_m2_s: 1e-10,
    }
}

#[test]
fn all_port_rcr_inertance_matches_independent_momentum_storage_and_reverse_heat() {
    with_cx(|cx, _| {
        for (drive, datum) in [(-4.0, 0.0), (4.0, 1e6)] {
            let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
            let system = FluxSystem::new(&mesh);
            let law = Caloric { datum };
            let (old, cells) = inventory(&mesh, &law, &[298.0, 308.0], true);
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
            let (dt, qp) = (0.2, 0.2);
            let previous_x = [-0.3, 0.4];
            let loads = [
                MaterialStokesWindkesselPort {
                    reservoir_pressure_pa: drive,
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
            // Independent elimination of the three physical BE equations:
            // If*dQ/dt = p_left-p_right; C_i*dx_i/dt = s_i*Q-x_i/Rd_i.
            let b =
                loads.map(|l| l.compliance_m2_per_pa / dt + 1.0 / l.distal_resistance_pa_s_per_m2);
            let a: [f64; 2] =
                std::array::from_fn(|i| loads[i].compliance_m2_per_pa * previous_x[i] / dt / b[i]);
            let total_inertia = RHO + loads.iter().map(|l| l.inertance_pa_s2_per_m2).sum::<f64>();
            let q = (total_inertia * qp / dt + drive + a[0] - a[1])
                / (total_inertia / dt + 1.0 + b[0].recip() + b[1].recip());
            let expected_x = [a[0] - q / b[0], a[1] + q / b[1]];
            let velocity = |_| [q, 0.0];
            let previous = moments(&mesh, &|_| [qp, 0.0]);
            let original_velocity = previous.clone();
            let original_cells = old.cells().clone();
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
            let observed_q = Cell::new(f64::NAN);
            let solved = owner
                .advance_windkessel(
                    cx,
                    &old,
                    dt,
                    LaminarWindkesselInput {
                        fluid: LaminarEnthalpyInput {
                            history: LaminarEnthalpyHistory {
                                velocity_m_s: &previous,
                                time_s: 0.0,
                            },
                            force_density: &zero,
                            boundaries: &boundaries,
                            slip_heat_routes: &routes,
                        },
                        ports: &[ports[0].as_slice(), ports[1].as_slice()],
                        loads: &loads,
                        previous_compliance_pressures_pa: &previous_x,
                    },
                    policy(),
                    &[0.0; 3],
                    |_, t, _| Ok(source(t, 0.4)),
                    |_, _, endpoint| {
                        observed_q.set(endpoint.outward_volume_rates_m2_s[1]);
                        for (balance, expected) in endpoint.balances.iter().zip(expected_x) {
                            near(balance.compliance_pressure_pa, expected, 2e-8);
                        }
                        Ok::<_, ()>(exchange(&mesh, law.h(330.0)))
                    },
                )
                .unwrap();
            let circuit = &solved.report().exchange.circuit;
            assert!(solved.report().endpoint_evaluations > 1);
            near(observed_q.get(), q, 2e-8);
            near(circuit.duration_s, dt, 0.0);
            near(circuit.depth_m, DEPTH, 0.0);
            near(circuit.outward_volume_rates_m2_s[0], -q, 2e-8);
            near(circuit.outward_volume_rates_m2_s[1], q, 2e-8);
            let next_x = circuit.compliance_pressures_pa();
            let mut physical_change = 0.5 * RHO * (q * q - qp * qp);
            let mut damping = 0.5 * RHO * (q - qp).powi(2);
            let mut resistor_heat = 0.0;
            for i in 0..2 {
                let load = loads[i];
                let sign = if i == 0 { -1.0 } else { 1.0 };
                let x = expected_x[i];
                let p = load.reservoir_pressure_pa
                    + x
                    + load.proximal_resistance_pa_s_per_m2 * sign * q
                    + load.inertance_pa_s2_per_m2 * sign * (q - qp) / dt;
                let balance = circuit.balances[i];
                near(next_x[i], x, 2e-8);
                near(circuit.pressure_offsets_pa[i], p, 3e-8);
                near(balance.port_power_w_per_m, p * sign * q, 3e-8);
                near(
                    balance.reservoir_power_w_per_m,
                    load.reservoir_pressure_pa * sign * q,
                    3e-8,
                );
                near(
                    balance.compliance_energy_j_per_m,
                    0.5 * load.compliance_m2_per_pa * x * x,
                    2e-8,
                );
                near(
                    balance.inertance_energy_j_per_m,
                    0.5 * load.inertance_pa_s2_per_m2 * q * q,
                    2e-8,
                );
                near(balance.power_balance_residual_w_per_m, 0.0, 3e-8);
                let stored = 0.5 * load.compliance_m2_per_pa * (x * x - previous_x[i].powi(2))
                    + 0.5 * load.inertance_pa_s2_per_m2 * (q * q - qp * qp);
                let defect = 0.5 * load.compliance_m2_per_pa * (x - previous_x[i]).powi(2)
                    + 0.5 * load.inertance_pa_s2_per_m2 * (q - qp).powi(2);
                let loss = dt
                    * (load.proximal_resistance_pa_s_per_m2 * q * q
                        + x * x / load.distal_resistance_pa_s_per_m2);
                near(balance.stored_energy_change_rate_w_per_m * dt, stored, 3e-8);
                near(balance.bdf_energy_defect_w_per_m * dt, defect, 3e-8);
                near(
                    (balance.proximal_dissipation_w_per_m + balance.distal_dissipation_w_per_m)
                        * dt,
                    loss,
                    3e-8,
                );
                physical_change += stored;
                damping += defect;
                resistor_heat += loss;
            }
            near(
                drive * q * dt,
                physical_change + damping + resistor_heat,
                3e-8,
            );
            assert!(resistor_heat > 0.0);
            let (error, divergence) =
                system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
            assert!(error < 2e-8 && divergence < 2e-8);
            let first = mesh
                .edges
                .iter()
                .find(|e| e.tris.1 == usize::MAX && q * e.normal[0] < 0.0)
                .unwrap()
                .tris
                .0;
            let second = 1 - first;
            let alpha = 2.0 * q.abs() * dt;
            let first_t =
                (old.temperature_k(cells[first]).unwrap() + alpha * 330.0) / (1.0 + alpha);
            let second_t =
                (old.temperature_k(cells[second]).unwrap() + alpha * first_t) / (1.0 + alpha);
            near(
                solved.thermal().temperature_k(cells[first]).unwrap(),
                first_t,
                3e-7,
            );
            near(
                solved.thermal().temperature_k(cells[second]).unwrap(),
                second_t,
                3e-7,
            );
            near(solved.thermal().temperature_k(WALL).unwrap(), 300.0, 1e-10);
            near(
                solved.report().external_stream_enthalpy_j,
                RHO * DEPTH * q.abs() * dt * CP * (330.0 - second_t),
                3e-7,
            );
            near(solved.report().bulk_viscous_heat_j, 0.0, 1e-12);
            near(solved.report().slip_heat_j, 0.0, 1e-12);
            near(solved.report().thermal_balance_residual_j, 0.0, 2e-8);
            for (i, sample) in solved.flow().samples().iter().enumerate() {
                near(
                    sample.material().source().temperature_k(),
                    solved.thermal().temperature_k(cells[i]).unwrap(),
                    1e-12,
                );
            }
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(previous, original_velocity);
            assert_eq!(previous_x, [-0.3, 0.4]);
        }
    });
}

#[test]
fn circuit_source_stream_budget_and_cancel_refusals_leave_all_histories_unchanged() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
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
        let (dt, q, xp) = (0.04, 0.2, 0.1);
        let velocity = |_| [q, 0.0];
        let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
        let previous = moments(&mesh, &velocity);
        let previous_x = [xp];
        let original_velocity = previous.clone();
        let original_cells = old.cells().clone();
        let x = (0.6 * xp / dt + q) / (0.6 / dt + 1.0 / 2.1);
        let loads = [MaterialStokesWindkesselPort {
            reservoir_pressure_pa: -0.3 * q - x,
            proximal_resistance_pa_s_per_m2: 0.3,
            distal_resistance_pa_s_per_m2: 2.1,
            compliance_m2_per_pa: 0.6,
            inertance_pa_s2_per_m2: 0.4,
        }];
        let outlet: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].normal[0] > 0.5)
            .collect();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                if mesh.edges[edge].normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &traction,
                        velocity: &velocity,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                }
            })
            .collect();
        let ports = [outlet.as_slice()];
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
            ports: &ports,
            loads: &loads,
            previous_compliance_pressures_pa: &previous_x,
        };
        for mode in 0..6 {
            let mut attempted_input = input;
            let mut controls = policy();
            if mode == 1 {
                attempted_input.previous_compliance_pressures_pa = &[];
            }
            if mode == 3 {
                controls.coupled.maximum_endpoint_evaluations = 1;
            }
            let exchanges = Cell::new(0usize);
            let attempted = owner.advance_windkessel(
                cx,
                &old,
                dt,
                attempted_input,
                controls,
                &[0.0; 2],
                |_, t, _| {
                    if mode == 0 {
                        return Err(material_error("original temperature source refused"));
                    }
                    Ok(source(t, 0.4))
                },
                |_, _, endpoint| {
                    near(endpoint.balances[0].compliance_pressure_pa, x, 2e-8);
                    exchanges.set(exchanges.get() + 1);
                    if mode == 2 || (mode == 5 && exchanges.get() > 1) {
                        return Err("original stream/heat source refused");
                    }
                    let mut heat = exchange(&mesh, law.h(320.0));
                    if mode == 4 {
                        heat.inlet_enthalpies.fill(None);
                    }
                    Ok(heat)
                },
            );
            match mode {
                0 | 1 => assert!(matches!(attempted, Err(LaminarEnthalpyError::Momentum(_)))),
                2 | 5 => assert!(matches!(attempted, Err(LaminarEnthalpyError::Exchange(_)))),
                3 => assert!(matches!(
                    attempted,
                    Err(LaminarEnthalpyError::EndpointBudget)
                )),
                _ => assert!(matches!(attempted, Err(LaminarEnthalpyError::Invalid(_)))),
            }
            assert_eq!(previous, original_velocity);
            assert_eq!(previous_x, [xp]);
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
        }
        let attempted = owner.advance_windkessel(
            cx,
            &old,
            dt,
            input,
            policy(),
            &[0.0; 2],
            |_, t, _| Ok(source(t, 0.4)),
            |_, _, endpoint| {
                near(endpoint.balances[0].compliance_pressure_pa, x, 2e-8);
                gate.request();
                Ok::<_, ()>(exchange(&mesh, law.h(320.0)))
            },
        );
        assert!(matches!(attempted, Err(LaminarEnthalpyError::Cancelled)));
        assert_eq!(previous, original_velocity);
        assert_eq!(previous_x, [xp]);
        assert_eq!(old.cells(), &original_cells);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn original_thermal_rheology_changes_circuit_endpoint_and_slip_heat_together() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (target_t, dt, wall_speed, beta): (f64, f64, f64, f64) = (310.0, 2.0, 5.0, 1.0);
        let curve = ThermalCarreauYasuda::new(
            CarreauYasuda::new(0.4, 0.2, 1.2, 1.0, 2.0, [0.0, 1e4]).unwrap(),
            TemperatureShift::new(
                TemperatureShiftModel::Arrhenius {
                    activation_energy_j_mol: 12_000.0,
                },
                300.0,
                [280.0, 360.0],
            )
            .unwrap(),
            ViscosityShiftMode::VerticalAndRate,
        );
        // Independent Arrhenius equation and affine Navier-slip wall balance.
        let mu = 0.4
            * (12_000.0 / fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * (target_t.recip() - 1.0 / 300.0))
                .exp();
        let gamma = wall_speed / (1.0 + mu / beta);
        let q = gamma / 2.0;
        let slip_w_per_m = beta * (wall_speed - gamma).powi(2);
        let top = mesh
            .edges
            .iter()
            .position(|e| e.tris.1 == usize::MAX && e.normal[1] > 0.5)
            .unwrap();
        let top_cell = mesh.edges[top].tris.0;
        // At the uniform target temperature advection cancels cell by cell.
        // Backward-Euler initial energies are manufactured from independent
        // bulk and actual relative-slip heat, not from a native flow result.
        let initial_t: Vec<_> = mesh
            .areas
            .iter()
            .enumerate()
            .map(|(i, &area)| {
                target_t
                    - dt * (mu * gamma * gamma * area
                        + if i == top_cell { slip_w_per_m } else { 0.0 })
                        / (RHO * area * CP)
            })
            .collect();
        let (old, cells) = inventory(&mesh, &law, &initial_t, true);
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
        let velocity = |p: [f64; 2]| [gamma * p[1], 0.0];
        let wall = |_| [wall_speed, 0.0];
        let traction = |_: [f64; 2], n: [f64; 2]| [mu * gamma * n[1], mu * gamma * n[0]];
        let previous = moments(&mesh, &velocity);
        let previous_x = [-0.3];
        let (rp, rd, compliance) = (0.3, 2.1, 0.6);
        let x = (compliance * previous_x[0] / dt + q) / (compliance / dt + 1.0 / rd);
        let load = MaterialStokesWindkesselPort {
            reservoir_pressure_pa: -rp * q - x,
            proximal_resistance_pa_s_per_m2: rp,
            distal_resistance_pa_s_per_m2: rd,
            compliance_m2_per_pa: compliance,
            inertance_pa_s2_per_m2: 0.4,
        };
        let outlet: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].normal[0] > 0.5)
            .collect();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                if edge == top {
                    MaterialNavierStokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_pa_s_per_m: beta,
                    }
                } else if mesh.edges[edge].normal[1] < -0.5 {
                    MaterialNavierStokesBoundary::Velocity { edge, value: &zero }
                } else {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &traction,
                        velocity: &velocity,
                    }
                }
            })
            .collect();
        let routes = [StokesSlipHeatRoute {
            edge: top,
            wall_cell: WALL,
            fluid_heat_fraction: 1.0,
        }];
        let off_endpoint_queries = Cell::new(0usize);
        let solved = owner
            .advance_windkessel(
                cx,
                &old,
                dt,
                LaminarWindkesselInput {
                    fluid: LaminarEnthalpyInput {
                        history: LaminarEnthalpyHistory {
                            velocity_m_s: &previous,
                            time_s: 0.0,
                        },
                        force_density: &zero,
                        boundaries: &boundaries,
                        slip_heat_routes: &routes,
                    },
                    ports: &[&outlet],
                    loads: &[load],
                    previous_compliance_pressures_pa: &previous_x,
                },
                policy(),
                &[0.0; 3],
                |point, t, _| {
                    if (t - target_t).abs() > 0.1 {
                        off_endpoint_queries.set(off_endpoint_queries.get() + 1);
                    }
                    curve
                        .evaluate(t, point.shear_rate_per_s)?
                        .bind_thermal(source(t, 0.4))
                },
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(target_t))),
            )
            .unwrap();
        assert!(off_endpoint_queries.get() > 0);
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
        assert!(error < 3e-7 && divergence < 3e-7, "{error}, {divergence}");
        let circuit = &solved.report().exchange.circuit;
        near(circuit.outward_volume_rates_m2_s[0], q, 3e-7);
        near(circuit.pressure_offsets_pa[0], 0.0, 3e-7);
        near(circuit.balances[0].compliance_pressure_pa, x, 3e-7);
        for (i, sample) in solved.flow().samples().iter().enumerate() {
            near(
                solved.thermal().temperature_k(cells[i]).unwrap(),
                target_t,
                3e-7,
            );
            near(sample.transport().dynamic_viscosity_pa_s(), mu, 3e-8);
            let original = sample.material().source();
            near(
                original.rheology().temperature().temperature_k(),
                target_t,
                3e-7,
            );
            near(original.rheology().shear_rate_per_s(), gamma, 3e-7);
            near(original.source().temperature_k(), target_t, 3e-7);
        }
        near(
            solved.report().bulk_viscous_heat_j,
            mu * gamma * gamma * DEPTH * dt,
            3e-7,
        );
        near(solved.report().slip_heat_j, slip_w_per_m * DEPTH * dt, 3e-7);
        near(
            solved.report().thermal.stored_energy_change_j,
            (mu * gamma * gamma + slip_w_per_m) * DEPTH * dt,
            3e-7,
        );
        near(solved.report().external_stream_enthalpy_j, 0.0, 3e-7);
        near(solved.report().other_external_heat_j, 0.0, 0.0);
        near(solved.report().thermal_balance_residual_j, 0.0, 2e-8);
        near(solved.thermal().temperature_k(WALL).unwrap(), 300.0, 1e-10);
        assert!(circuit.balances[0].distal_dissipation_w_per_m > 0.0);
        assert_eq!(old.cells(), &original_cells);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(previous_x, [-0.3]);
    });
}

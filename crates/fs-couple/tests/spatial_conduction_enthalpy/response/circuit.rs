//! Independent nonlinear RC/thermal and RL midpoint oracles.
use super::*;
use fs_couple::circuit_enthalpy::conduction::spatial::circuit::response::SpatialConductionResponseCircuit;
use fs_em::conduction::response::AdaptiveConductionRelaxation;
use fs_phs::circuit::{Branch, CircuitGraph};

fn matrix_controls() -> SpatialConductanceControls {
    SpatialConductanceControls {
        maximum_terminals: 3,
        maximum_field_solves: 2,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-10,
    }
}

fn terminals(problem: &SpatialConductionResponse) -> Vec<Vec<usize>> {
    [0.0, 2.0]
        .map(|x| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0])
                        == [x; 3]
                })
                .collect()
        })
        .to_vec()
}

fn adaptive(problem: SpatialConductionResponse) -> SpatialConductionResponse {
    problem
        .with_adaptive_relaxation(AdaptiveConductionRelaxation {
            initial_weight: 0.5,
            minimum_weight: 0.05,
            maximum_weight: 1.0,
        })
        .unwrap()
}

#[test]
fn nonlinear_3d_capacitor_current_source_and_enthalpy_share_one_accepted_history() {
    let r = [1.0, 2.0];
    let beta = [0.08, 0.06];
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        for sign in [1.0, -1.0] {
            let (a, ca) = nonlinear_material(r[0], CP[0], ALPHA[0], beta[0], coordinate);
            let (b, cb) = nonlinear_material(r[1], CP[1], ALPHA[1], beta[1], coordinate);
            let (cards, curves) = ([a, b], [ca, cb]);
            let context = query(300.0, coordinate);
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let cells = bindings(&cards, &curves, &context, &key);
            let regions = regions();
            let network = SpatialConductionThermalNetwork {
                cells: &cells,
                links: &LINKS,
                regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let mut thermal = initial(&curves, 999.0, None);
            let field = adaptive(problem(&cards, coordinate));
            let faces = terminals(&field);
            let graph = CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Capacitor { farads: 1.0 }),
                    (1, 0, Branch::CurrentSource { port: 0 }),
                ],
                transformers: vec![],
            };
            let (mut charge, ports) = (sign * 4.0, [sign * 0.25]);
            let mut circuit = SpatialConductionResponseCircuit::new(
                field,
                graph,
                faces,
                vec![0, 1],
                0.0,
                vec![],
                vec![charge],
                &ports,
                matrix_controls(),
                limits(),
                field_controls(),
                controls(),
                || false,
            )
            .unwrap();
            let dt = 0.25;
            let mut old_theta = [0.0; 2];
            let mut solve = config();
            solve.maximum_circuit_trials = 4096;
            for interval in 1..=2 {
                // For C=1 F, Vmid = q_old + h*u - h*J. Independently
                // eliminate the two original thermal balances at every J.
                let h = dt / 2.0;
                let (mut lo, mut hi) = (0.0, 4.0);
                for _ in 0..90 {
                    let j = (lo + hi) / 2.0;
                    let (_, e) = endpoint(
                        j,
                        dt,
                        old_theta,
                        r,
                        beta,
                        coordinate == ConductorMagnitude::ElectricField,
                    );
                    if e.iter().sum::<f64>() + h * j < charge.abs() + h * ports[0].abs() {
                        lo = j;
                    } else {
                        hi = j;
                    }
                }
                let j = (lo + hi) / 2.0;
                let (theta, e) = endpoint(
                    j,
                    dt,
                    old_theta,
                    r,
                    beta,
                    coordinate == ConductorMagnitude::ElectricField,
                );
                let voltage = sign * e.iter().sum::<f64>();
                let final_charge = 2.0 * voltage - charge;
                let expected_heat = dt * j * voltage.abs();
                let before = thermal.energy_j();
                let target = interval as f64 * dt;
                if interval == 1 && sign > 0.0 {
                    let mut starved = solve;
                    starved.maximum_circuit_trials = 1;
                    assert!(matches!(
                        network.step_response_circuit_candidate(
                            &thermal,
                            target,
                            &circuit,
                            &ports,
                            &[],
                            matrix_controls(),
                            limits(),
                            field_controls(),
                            controls(),
                            starved,
                            &[0.0; 2],
                            || false,
                        ),
                        Err(CircuitEnthalpyError::CircuitBudget)
                    ));
                    assert!(matches!(
                        network.step_response_circuit_candidate(
                            &thermal,
                            target,
                            &circuit,
                            &ports,
                            &[],
                            matrix_controls(),
                            limits(),
                            field_controls(),
                            controls(),
                            solve,
                            &[0.0; 2],
                            || true,
                        ),
                        Err(CircuitEnthalpyError::Cancelled)
                    ));
                    let mut one_field = controls();
                    one_field.maximum_field_solves = 1;
                    assert!(matches!(
                        circuit.step_candidate(
                            target,
                            &ports,
                            matrix_controls(),
                            limits(),
                            field_controls(),
                            one_field,
                            solve.circuit,
                            || false,
                        ),
                        Err(CircuitEnthalpyError::Source(
                            ConductionThermalError::SpatialResponse(
                                SpatialConductionResponseError::FieldSolveBudget { .. }
                            )
                        ))
                    ));
                    near(circuit.charges_c()[0], charge);
                    assert_eq!(circuit.time_s(), 0.0);
                    near(thermal.energy_j(), before);
                }
                let next = network
                    .step_response_circuit_candidate(
                        &thermal,
                        target,
                        &circuit,
                        &ports,
                        &[],
                        matrix_controls(),
                        limits(),
                        field_controls(),
                        controls(),
                        solve,
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                near(next.circuit().charges_c()[0], final_charge);
                let report = next.report();
                let step = &report.electrical.circuit;
                near(step.terminal_voltage_v[1], voltage);
                near(step.terminal_current_a[1], sign * j);
                near(step.terminal_current_a[0], -sign * j);
                near(step.heat_j, expected_heat);
                near(
                    step.stored_change_j,
                    0.5 * (final_charge * final_charge - charge * charge),
                );
                near(step.supplied_work_j, dt * ports[0] * voltage);
                near(step.energy_residual_j, 0.0);
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(
                        next.thermal().temperature_k(cell).unwrap(),
                        300.0 + theta[i],
                    );
                    for tet in 6 * i..6 * i + 6 {
                        near(step.spatial.field.current_density_a_m2[tet][0], -sign * j);
                        near(step.spatial.field.electric_field_v_m[tet][0], -sign * e[i]);
                        near(step.spatial.field.cell_power_w[tet], j * e[i] / 6.0);
                        let source = &report.electrical.evaluated[tet].conductor;
                        near(source.material().query().axes()["T"], 300.0 + theta[i]);
                        near(
                            source.material().query().axes()["operating"],
                            if coordinate == ConductorMagnitude::ElectricField {
                                e[i]
                            } else {
                                j
                            },
                        );
                        assert!(source.has_same_electrical_source(&circuit.conductors()[tet]));
                    }
                }
                near(next.thermal().energy_j() - before, expected_heat);
                near(report.deposited_joule_j, expected_heat);
                near(report.total_energy_residual_j, 0.0);
                assert_eq!(report.electrical.field_solves, report.electrical.iterations);
                assert!(report.electrical.iterations > 1);
                assert_eq!(report.circuit_trials % 2, 0);
                assert!(report.circuit_trials >= 2 * report.electrical.iterations);
                assert!(next.circuit().adaptive_relaxation().is_some());
                assert_eq!(context.axes()["operating"], 0.0);
                (thermal, circuit, _) = next.into_parts();
                charge = final_charge;
                old_theta = theta;
            }
        }
    }
}

#[test]
fn nonlinear_3d_inductor_discharge_closes_actual_midpoint_resistance_in_both_polarities() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let r = [1.0, 2.0];
    let beta = [0.2, 0.1];
    let (a, _) = nonlinear_material(r[0], CP[0], ALPHA[0], beta[0], coordinate);
    let (b, _) = nonlinear_material(r[1], CP[1], ALPHA[1], beta[1], coordinate);
    let cards = [a, b];
    for sign in [1.0, -1.0] {
        let problem = adaptive(problem(&cards, coordinate));
        let faces = terminals(&problem);
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Inductor { henries: 2.0 })],
            transformers: vec![],
        };
        let mut flux = sign * 4.0;
        let mut circuit = SpatialConductionResponseCircuit::new(
            problem,
            graph,
            faces,
            vec![0, 1],
            0.0,
            vec![flux],
            vec![],
            &[],
            matrix_controls(),
            limits(),
            field_controls(),
            controls(),
            || false,
        )
        .unwrap();
        let dt = 0.2;
        for interval in 1..=2 {
            // L*Jmid + h*Jmid*(R0 + beta*|Jmid|) = |flux_old|.
            let h = dt / 2.0;
            let a = h * beta.iter().sum::<f64>();
            let b = 2.0 + h * r.iter().sum::<f64>();
            let j = 2.0 * flux.abs() / (b + (b * b + 4.0 * a * flux.abs()).sqrt());
            let resistance = r.iter().sum::<f64>() + beta.iter().sum::<f64>() * j;
            let voltage = -sign * j * resistance;
            let final_flux = flux + dt * voltage;
            let next = circuit
                .step_candidate(
                    interval as f64 * dt,
                    &[],
                    matrix_controls(),
                    limits(),
                    field_controls(),
                    controls(),
                    config().circuit,
                    || false,
                )
                .unwrap();
            near(next.circuit().fluxes_wb()[0], final_flux);
            near(next.step().circuit.terminal_voltage_v[1], voltage);
            near(next.step().circuit.terminal_current_a[1], -sign * j);
            near(next.step().circuit.heat_j, dt * j * j * resistance);
            near(
                next.step().circuit.stored_change_j,
                (final_flux * final_flux - flux * flux) / 4.0,
            );
            near(next.step().circuit.energy_residual_j, 0.0);
            near(next.step().circuit.supplied_work_j, 0.0);
            assert!(next.step().iterations > 1);
            (circuit, _) = next.into_parts();
            flux = final_flux;
        }
    }
}

#[test]
fn nonlinear_three_electrode_rlc_uses_full_native_response_and_counts_every_basis() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let (a, _) = nonlinear_material(1.0, CP[0], ALPHA[0], 0.2, coordinate);
    let (b, _) = nonlinear_material(2.0, CP[1], ALPHA[1], 0.1, coordinate);
    let cards = [a, b];
    let original = problem(&cards, coordinate);
    // Two original volumes meet at an equipotential electrode; disconnected
    // material vertices retain each side's exterior face and current reaction.
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (other_positions, other_tets) = slab(&[1.0, 2.0]);
    let offset = positions.len() as u32;
    positions.extend(other_positions);
    tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    let field = adaptive(
        SpatialConductionResponse::new(
            positions,
            tets,
            original.sources().to_vec(),
            limits(),
            &mut || false,
        )
        .unwrap(),
    );
    let faces = [0.0, 1.0, 2.0]
        .map(|x| {
            field
                .exterior_faces()
                .filter(|&face| {
                    field.mesh().faces[face].map(|v| field.positions_m()[v as usize][0]) == [x; 3]
                })
                .collect::<Vec<_>>()
        })
        .to_vec();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Inductor { henries: 2.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let (dt, h, voltage, current, port) = (0.2, 0.1, 2.0, 0.5, 0.2);
    // Manufacture the accepted C/L history from independently chosen midpoint
    // voltage/current. J + 0.2*J²=2 on the left; right R=2+0.1*0.5.
    let left_j = 4.0_f64 / (1.0 + 2.6_f64.sqrt());
    let right_r = 2.0 + 0.1 * current;
    let right_v = voltage - right_r * current;
    let charge = voltage + h * (left_j + current - port);
    let flux = 2.0 * current - h * right_v;
    let circuit = SpatialConductionResponseCircuit::new(
        field,
        graph,
        faces,
        vec![0, 1, 2],
        0.0,
        vec![flux],
        vec![charge],
        &[port],
        matrix_controls(),
        limits(),
        field_controls(),
        controls(),
        || false,
    )
    .unwrap();
    let mut starved = controls();
    starved.maximum_field_solves = 1;
    assert!(matches!(
        circuit.step_candidate(
            dt,
            &[port],
            matrix_controls(),
            limits(),
            field_controls(),
            starved,
            config().circuit,
            || false
        ),
        Err(CircuitEnthalpyError::Source(
            ConductionThermalError::SpatialResponse(
                SpatialConductionResponseError::FieldSolveBudget {
                    completed: 0,
                    required: 2,
                    limit: 1
                }
            )
        ))
    ));
    let next = circuit
        .step_candidate(
            dt,
            &[port],
            matrix_controls(),
            limits(),
            field_controls(),
            controls(),
            config().circuit,
            || false,
        )
        .unwrap();
    let report = next.step();
    let step = &report.circuit;
    near(step.terminal_voltage_v[1], voltage);
    near(step.terminal_voltage_v[2], right_v);
    near(step.terminal_current_a[0], -left_j);
    near(step.terminal_current_a[1], left_j + current);
    near(step.terminal_current_a[2], -current);
    near(next.circuit().charges_c()[0], 2.0 * voltage - charge);
    near(next.circuit().fluxes_wb()[0], 4.0 * current - flux);
    near(
        step.spatial.field.cell_power_w[..6].iter().sum::<f64>(),
        voltage * left_j,
    );
    near(
        step.spatial.field.cell_power_w[6..].iter().sum::<f64>(),
        right_r * current * current,
    );
    near(
        step.heat_j,
        dt * (voltage * left_j + right_r * current * current),
    );
    near(step.supplied_work_j, dt * port * voltage);
    near(step.stored_change_j + step.heat_j, step.supplied_work_j);
    near(step.energy_residual_j, 0.0);
    assert_eq!(report.field_solves, 2 * report.iterations);
    assert!(report.iterations > 1);
    near(circuit.charges_c()[0], charge);
    near(circuit.fluxes_wb()[0], flux);
    assert_eq!(circuit.time_s(), 0.0);
}

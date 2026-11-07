//! G1/G3/G4: actual electrode circuits, magnetic entropy and finite heat.
use super::*;
use fs_couple::circuit_enthalpy::ResistorCell;
use fs_em::spatial::vector_potential::induction::terminals::{
    SpatialTerminalInduction,
    circuit::{
        SpatialTerminalCircuit, SpatialTerminalCircuitConfig, SpatialTerminalCircuitControls,
        SpatialTerminalCircuitInitial,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};

fn circuit_controls() -> SpatialTerminalCircuitControls {
    SpatialTerminalCircuitControls {
        field: field_controls(),
        voltage_probe_v: 0.1,
        current_probe_a: 0.1,
        maximum_field_trials: 4,
        relative_terminal_tolerance: 1e-9,
        voltage_tolerance_v: 2e-8,
        current_tolerance_a: 2e-8,
        energy_tolerance_j: 2e-8,
    }
}
fn field_limits() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 1000,
        maximum_conducting_cells: 100,
        maximum_assembly_terms: 10000,
    }
}
fn electrodes(problem: &SpatialVectorMagnetostatics) -> Vec<Vec<usize>> {
    [0.0, 1.0]
        .into_iter()
        .map(|x| {
            problem
                .exterior_faces()
                .into_iter()
                .filter(|&face| {
                    problem.mesh().faces[face]
                        .iter()
                        .all(|&v| problem.positions_m()[v as usize][0] == x)
                })
                .collect()
        })
        .collect()
}
fn connect<'m>(
    field: SpatialInduction<'m>,
    graph: CircuitGraph,
    charges_c: Vec<f64>,
    external_inputs: Vec<f64>,
    reference_ohm: f64,
) -> SpatialTerminalCircuit<'m> {
    let faces = electrodes(field.problem());
    let terminals =
        SpatialTerminalInduction::new(field, &faces, field_limits(), &mut || false).unwrap();
    SpatialTerminalCircuit::new(
        terminals,
        graph,
        vec![1, 0],
        vec![],
        SpatialTerminalCircuitInitial {
            fluxes_wb: vec![],
            charges_c,
            external_inputs,
        },
        SpatialTerminalCircuitConfig {
            reference_resistance_ohm: reference_ohm,
            maximum_native_states: 32,
            maximum_interface_scalars: 10000,
        },
        &mut || false,
    )
    .unwrap()
}
fn series_graph(resistance: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (2, 1, Branch::Resistor { ohms: resistance }),
        ],
        transformers: vec![],
    }
}
fn pair_thermal(curve: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    let events = [CellId(7), CellId(8)]
        .into_iter()
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: id.0,
                body_id: id.0,
                material_id: MaterialId(9),
                mass_kg: 1.0,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bindings<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a fs_matdb::QueryPoint,
    pair: bool,
) -> Vec<MaterialCell<'a>> {
    let ids = if pair {
        vec![CellId(7), CellId(8)]
    } else {
        vec![CellId(7)]
    };
    ids.into_iter()
        .map(|cell| MaterialCell {
            cell,
            material_id: MaterialId(9),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect()
}

#[test]
fn circuit_driven_sigma_and_mu_feedback_matches_independent_whitney_thermal_elimination() {
    // Unit-cube integration: int W_x=1/6, int W²=1/5, K(T)=4/mu(T).
    // The resistor is an independently heated second body. The magnetic
    // entropy is S=2*mu_T*a²/mu², so Qrev=-P-delta(T*S).
    let (rho, alpha, mu0, resistance, dt): (f64, f64, f64, f64, f64) = (0.2, 0.01, 1.0, 0.4, 0.1);
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for slope in [-0.0004, 0.0, 0.0004] {
            let card = source::absolute_affine_magnetic_card(rho, alpha, mu0, slope, basis);
            let curve = curve(&card);
            let problem = problem(&card, basis);
            let mut circuit = connect(
                initial(&problem, basis, 0.0, mu0),
                series_graph(resistance),
                vec![],
                vec![0.0],
                0.7,
            );
            let context = source::absolute_point(300.0);
            let cells = bindings(&card, &curve, &context, true);
            let regions = regions();
            let laws = laws(&card, &problem, 5.0);
            let network = SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let mut heat = pair_thermal(&curve);
            let (mut old_a, mut old_t, mut load_t) = (0.0, 300.0, 300.0);
            for (index, source_v) in [12.0, 0.0, -8.0].into_iter().enumerate() {
                let mu_old = mu0 + slope * (old_t - 300.0);
                let at = |t: f64| {
                    let mu = mu0 + slope * (t - 300.0);
                    let s = sigma(t, rho, alpha, basis);
                    let f = s / 6.0;
                    let denominator = s / 5.0 + 2.0 * dt / mu;
                    let magnetic = 2.0 * (mu_old.recip() + mu.recip()) * old_a;
                    let y = s - f * f / denominator;
                    let i0 = f * magnetic / denominator;
                    let voltage = (source_v - resistance * i0) / (1.0 + resistance * y);
                    let current = i0 + y * voltage;
                    let rate = (f * voltage - magnetic) / denominator;
                    let a = old_a + dt * rate;
                    let joule =
                        dt * s * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
                    let parameter = 2.0 * old_a * a * (mu.recip() - mu_old.recip());
                    let ts = 2.0
                        * slope
                        * (t * a * a / (mu * mu) - old_t * old_a * old_a / (mu_old * mu_old));
                    (a, voltage, current, joule, -parameter - ts, parameter)
                };
                let t = bracket(old_t, old_t + 20.0, |t| {
                    let (_, _, _, joule, reversible, _) = at(t);
                    4.0 * (t - old_t) - joule - reversible
                });
                let (a, voltage, current, joule, reversible, parameter) = at(t);
                let resistor_heat = dt * resistance * current * current;
                load_t += resistor_heat / 4.0;
                let proposal = network
                    .step_terminal_circuit_candidate(
                        &heat,
                        (index + 1) as f64 * dt,
                        &circuit,
                        &[source_v],
                        &[[0.0; 3]; 6],
                        &[[0.0; 3]; 8],
                        &[ResistorCell {
                            branch: 1,
                            cell: CellId(8),
                        }],
                        &laws,
                        thermal_controls(),
                        circuit_controls(),
                        assembly(),
                        &[0.0, 0.0],
                        || false,
                    )
                    .unwrap();
                let report = proposal.report();
                let step = &report.coupled.electrical;
                near(
                    proposal.thermal().temperature_k(CellId(7)).unwrap(),
                    t,
                    8e-8,
                );
                near(
                    proposal.thermal().temperature_k(CellId(8)).unwrap(),
                    load_t,
                    8e-8,
                );
                near(
                    proposal.circuit().field().potential_wb()[diagonal(&problem)],
                    a,
                    8e-8,
                );
                near(step.terminal_voltage_v[0], voltage, 8e-8);
                near(step.terminal_current_a[0], current, 2e-7);
                near(step.field.heat_j, joule, 2e-7);
                near(step.resistor_heat_j, resistor_heat, 2e-7);
                near(step.field.parameter_work_j, parameter, 2e-9);
                near(report.coupled.source.reversible_heat_j, reversible, 2e-9);
                near(proposal.deposited_heat_j(), joule + resistor_heat, 2e-7);
                near(step.external_work_j, dt * source_v * current, 2e-7);
                near(report.total_energy_residual_j, 0.0, 2e-8);
                assert!(step.cancelled_reference_heat_j > 0.0);
                near(heat.temperature_k(CellId(7)).unwrap(), old_t, 8e-8);
                let parts = proposal.into_parts();
                heat = parts.0;
                circuit = parts.1;
                old_a = a;
                old_t = t;
            }
        }
    }
}

#[test]
fn charged_capacitor_preserves_native_charge_while_heating_changes_field_current() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.03, 2.0, 0.0, basis);
    let curve = curve(&card);
    let problem = problem(&card, basis);
    let capacitance = 0.4;
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(
            1,
            0,
            Branch::Capacitor {
                farads: capacitance,
            },
        )],
        transformers: vec![],
    };
    let charge_index = assemble_circuit(&graph).unwrap().charge_index[0];
    let mut circuit = connect(
        initial(&problem, basis, 0.0, 2.0),
        graph,
        vec![8.0 * capacitance],
        vec![],
        1.3,
    );
    let context = source::absolute_point(300.0);
    let cells = bindings(&card, &curve, &context, false);
    let regions = regions();
    let laws = laws(&card, &problem, 5.0);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let mut heat = thermal(&curve, 0.0);
    let (mut old_a, mut old_t, mut charge) = (0.0, 300.0, 8.0 * capacitance);
    let dt = 0.05;
    for index in 0..8 {
        let at = |t: f64| {
            let s = sigma(t, 0.2, 0.03, basis);
            let f = s / 6.0;
            let d = s / 5.0 + dt;
            let y = s - f * f / d;
            let i0 = 2.0 * f * old_a / d;
            let voltage = (charge / capacitance - dt * i0 / (2.0 * capacitance))
                / (1.0 + dt * y / (2.0 * capacitance));
            let current = i0 + y * voltage;
            let rate = (f * voltage - 2.0 * old_a) / d;
            let a = old_a + dt * rate;
            let q = charge - dt * current;
            let joule = dt * s * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
            (a, q, current, joule)
        };
        let t = bracket(old_t, old_t + 10.0, |t| 4.0 * (t - old_t) - at(t).3);
        let (a, q, current, joule) = at(t);
        let candidate = network
            .step_terminal_circuit_candidate(
                &heat,
                (index + 1) as f64 * dt,
                &circuit,
                &[],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &[],
                &laws,
                thermal_controls(),
                circuit_controls(),
                assembly(),
                &[0.0],
                || false,
            )
            .unwrap();
        let step = &candidate.report().coupled.electrical;
        near(
            candidate.thermal().temperature_k(CellId(7)).unwrap(),
            t,
            6e-8,
        );
        near(candidate.circuit().native_state()[charge_index], q, 6e-8);
        near(
            candidate.circuit().field().potential_wb()[diagonal(&problem)],
            a,
            6e-8,
        );
        near(step.terminal_current_a[0], current, 2e-7);
        near(candidate.deposited_heat_j(), joule, 2e-7);
        near(
            step.circuit_stored_after_j,
            q * q / (2.0 * capacitance),
            8e-8,
        );
        near(step.external_work_j, 0.0, 1e-15);
        near(step.resistor_heat_j, 0.0, 1e-15);
        near(
            a * a + q * q / (2.0 * capacitance) + 4.0 * (t - 300.0),
            12.8,
            3e-7,
        );
        near(circuit.native_state()[charge_index], charge, 8e-8);
        let parts = candidate.into_parts();
        heat = parts.0;
        circuit = parts.1;
        old_a = a;
        old_t = t;
        charge = q;
    }
}

#[test]
fn missing_heat_receivers_source_changes_and_exhausted_trials_preserve_all_histories() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.01, 1.0, 0.0, basis);
    let curve = curve(&card);
    let problem = problem(&card, basis);
    let circuit = connect(
        initial(&problem, basis, 0.0, 1.0),
        series_graph(0.4),
        vec![],
        vec![0.0],
        0.7,
    );
    let heat = pair_thermal(&curve);
    let context = source::absolute_point(300.0);
    let cells = bindings(&card, &curve, &context, true);
    let regions = regions();
    let laws = laws(&card, &problem, 5.0);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let before = (
        heat.cells().clone(),
        circuit.native_state().to_vec(),
        circuit.field().potential_wb().to_vec(),
    );
    let destinations = [ResistorCell {
        branch: 1,
        cell: CellId(8),
    }];
    for bad in [
        vec![],
        vec![ResistorCell {
            branch: 0,
            cell: CellId(8),
        }],
        vec![ResistorCell {
            branch: 1,
            cell: CellId(99),
        }],
    ] {
        assert!(
            network
                .step_terminal_circuit_candidate(
                    &heat,
                    0.1,
                    &circuit,
                    &[12.0],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &bad,
                    &laws,
                    thermal_controls(),
                    circuit_controls(),
                    assembly(),
                    &[0.0, 0.0],
                    || false,
                )
                .is_err()
        );
    }
    let mut limited = thermal_controls();
    limited.maximum_circuit_trials = 2;
    assert!(matches!(
        network.step_terminal_circuit_candidate(
            &heat,
            0.1,
            &circuit,
            &[12.0],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &destinations,
            &laws,
            limited,
            circuit_controls(),
            assembly(),
            &[0.0, 0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(matches!(
        network.step_terminal_circuit_candidate(
            &heat,
            0.1,
            &circuit,
            &[12.0],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &destinations,
            &laws,
            thermal_controls(),
            circuit_controls(),
            assembly(),
            &[0.0, 0.0],
            || true,
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let other_card = source::absolute_affine_magnetic_card(0.3, 0.01, 1.0, 0.0, basis);
    let other_cells = bindings(&other_card, &curve, &context, true);
    let changed = SpatialMagneticThermalNetwork {
        cells: &other_cells,
        ..network
    };
    assert!(
        changed
            .step_terminal_circuit_candidate(
                &heat,
                0.1,
                &circuit,
                &[12.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &destinations,
                &laws,
                thermal_controls(),
                circuit_controls(),
                assembly(),
                &[0.0, 0.0],
                || false,
            )
            .is_err()
    );
    assert_eq!(heat.cells(), &before.0);
    assert_eq!(circuit.native_state(), before.1);
    assert_eq!(circuit.field().potential_wb(), before.2);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(circuit.time_s(), 0.0);
}

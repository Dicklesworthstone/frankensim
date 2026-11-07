//! Saturated physical electrodes and authored graph heat share finite receivers.
#[path = "../../../../../fs-em/tests/spatial_induction/terminals/circuit/common.rs"]
mod common;

use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell, spatial_induction::terminals::circuit::SpatialTerminalCircuitThermalNetwork,
};
use fs_em::spatial::vector_potential::induction::nonlinear::terminals::circuit::SpatialNonlinearTerminalCircuit;
use fs_phs::circuit::{Branch, CircuitGraph};

#[test]
fn saturated_electrode_circuit_heats_separate_bodies_and_releases_only_physical_stored_energy() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = card(0.2, 0.005, basis, true);
    let curve = curve(&card);
    let problem = problem(&card, basis);
    let faces: Vec<Vec<usize>> = [0.0, 1.0]
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
        .collect();
    let field = SpatialNonlinearTerminalInduction::new(
        initial(&problem, &card, basis, 0.0, true),
        &faces,
        SpatialInductionLimits {
            maximum_unknowns: 1000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 10000,
        },
        &mut || false,
    )
    .unwrap();
    let graph_r = 0.4;
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (2, 1, Branch::Resistor { ohms: graph_r }),
        ],
        transformers: vec![],
    };
    let mut circuit = SpatialNonlinearTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        common::initial(vec![], vec![], vec![0.0]),
        common::config(0.7),
        &mut || false,
    )
    .unwrap();

    // Reuse the original one-kilogram, four-J/K caloric fixture for both
    // bodies; the authored resistor's finite receiver is outside the FE body.
    let seed = thermal(&curve, 0.0);
    let core = seed.cells()[&CellId(7)].clone();
    let mut graph_cell = core.clone();
    graph_cell.id = CellId(8);
    graph_cell.body_id = 1;
    graph_cell.geometry.start_m = [2.0, 0.0, 0.0];
    graph_cell.geometry.end_m = [3.0, 0.0, 0.0];
    let mut heat = DepositionSimulation::new(
        [core, graph_cell]
            .into_iter()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect(),
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_events: 2,
            ..seed.budget()
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
    .0;
    let query = context(300.0);
    let [core_binding] = cells(&card, &curve, &query);
    let [mut graph_binding] = cells(&card, &curve, &query);
    graph_binding.cell = CellId(8);
    let bindings = [core_binding, graph_binding];
    let regions = regions();
    let graph_routes = [ResistorCell {
        branch: 1,
        cell: CellId(8),
    }];
    let network = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &bindings,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &graph_routes,
    };
    let (mut old_a, mut old_t, mut old_graph_t) = (0.0, 300.0, 300.0);
    let (mut physical, mut numerical, mut supplied) = (0.0, 0.0, 0.0);
    let dt = 0.1;
    for (index, drive) in [36.0, 0.0].into_iter().enumerate() {
        let endpoint = |t: f64| {
            let rho = 0.2 * (1.0 + 0.005 * (t - 300.0));
            let sigma = rho.recip();
            let mass = 31.0 * sigma / 180.0;
            // Exact cube integration gives 2H + M*a' = I/6 and
            // U = rho*I + a'/6. The native series graph imposes U+R*I=drive.
            let a = bracket(0.0, 1.999, |a| {
                let rate = (a - old_a) / dt;
                let current = (drive - rate / 6.0) / (graph_r + rho);
                2.0 * h(2.0 * a, true) + mass * rate - current / 6.0
            });
            let rate = (a - old_a) / dt;
            let current = (drive - rate / 6.0) / (graph_r + rho);
            let voltage = drive - graph_r * current;
            let core_heat = dt * (rho * current * current + mass * rate * rate);
            (a, current, voltage, core_heat, sigma)
        };
        let t = bracket(old_t, 500.0, |t| 4.0 * (t - old_t) - endpoint(t).3);
        let (a, current, voltage, core_heat, sigma) = endpoint(t);
        let graph_heat = dt * graph_r * current * current;
        let graph_t = old_graph_t + graph_heat / 4.0;
        let stored_change = energy(a, true) - energy(old_a, true);
        let numerical_loss = 2.0 * h(2.0 * a, true) * (a - old_a) - stored_change;
        let before = (
            heat.cells().clone(),
            circuit.field().potential_wb().to_vec(),
            circuit.time_s(),
        );
        let candidate = network
            .step_nonlinear_candidate(
                &heat,
                (index + 1) as f64 * dt,
                &circuit,
                &[drive],
                &[[0.0; 3]; 6],
                None,
                thermal_controls(),
                common::controls(field_controls()),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        let step = &report.coupled.electrical;
        near(
            candidate.thermal().temperature_k(CellId(7)).unwrap(),
            t,
            3e-8,
        );
        near(
            candidate.thermal().temperature_k(CellId(8)).unwrap(),
            graph_t,
            3e-8,
        );
        near(
            candidate.circuit().field().potential_wb()[diagonal(&problem)],
            a,
            3e-8,
        );
        near(step.terminal_voltage_v[0], voltage, 3e-7);
        near(step.terminal_current_a[0], current, 3e-7);
        near(step.terminal_current_a[1], -current, 3e-7);
        near(step.field.terminal_work_j, dt * voltage * current, 3e-7);
        near(step.field.heat_j, core_heat, 3e-7);
        assert_eq!(step.graph_resistors.len(), 1);
        assert_eq!(step.graph_resistors[0].branch, 1);
        near(step.graph_resistors[0].heat_j, graph_heat, 3e-7);
        near(step.resistor_heat_j, graph_heat, 3e-7);
        near(step.heat_j, core_heat + graph_heat, 3e-7);
        near(candidate.deposited_heat_j(), core_heat + graph_heat, 3e-7);
        near(
            report.coupled.thermal.stored_energy_change_j,
            core_heat + graph_heat,
            3e-7,
        );
        near(step.stored_change_j, stored_change, 3e-8);
        near(step.numerical_dissipation_j, numerical_loss, 3e-8);
        near(step.supplied_work_j, dt * drive * current, 3e-7);
        near(
            report.total_energy_residual_j,
            0.0,
            thermal_controls().total_energy_tolerance_j,
        );
        near(report.external_heat_j, 0.0, 1e-12);
        assert!(step.cancelled_reference_heat_j > 0.0);
        assert!(step.graph_materials.is_empty());
        for source in candidate.circuit().field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m().recip(), sigma, 3e-10);
            near(source.material().query().axes()["T"], t, 3e-8);
        }
        if index == 0 {
            assert!(2.0 * a > 1.0 && t > 310.0 && graph_t > t);
            assert!(current < endpoint(300.0).1);
            assert!(numerical_loss > 1.0);
        } else {
            assert!(stored_change < 0.0 && core_heat > 0.0 && graph_heat > 0.0);
            near(step.supplied_work_j, 0.0, 1e-12);
        }
        physical += candidate.deposited_heat_j();
        numerical += step.numerical_dissipation_j;
        supplied += step.supplied_work_j;
        near(step.stored_after_j + physical + numerical, supplied, 5e-7);
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(circuit.field().potential_wb(), before.1);
        assert_eq!(circuit.time_s(), before.2);
        old_a = a;
        old_t = t;
        old_graph_t = graph_t;
        let (next_heat, next_circuit, _) = candidate.into_parts();
        heat = next_heat;
        circuit = next_circuit;
    }
}

//! Solid-conductor saturation closes against authored native graph elements.
#[path = "../circuit/common.rs"]
mod common;
use super::*;
use common::{bracket, config, controls, initial};
use fs_em::spatial::vector_potential::induction::nonlinear::terminals::circuit::SpatialNonlinearTerminalCircuit;
use fs_phs::circuit::{Branch, CircuitGraph};

const SIGMA: f64 = 5.0;
const MASS: f64 = 31.0 * SIGMA / 180.0;

fn response(old: f64, next: f64, dt: f64) -> (f64, f64, f64) {
    let b = 2.0 * next.abs();
    let h = next.signum()
        * if b <= 1.0 {
            0.5 * b
        } else if b <= 2.0 {
            0.5 + 2.0 * (b - 1.0)
        } else {
            2.5 + 4.0 * (b - 2.0)
        };
    let rate = (next - old) / dt;
    // Independent Whitney integration: 2H + (31 sigma/180) adot = I/6,
    // with U = I/sigma + adot/6. These constants are not read from assembly.
    let current = 12.0 * h + 6.0 * MASS * rate;
    (h, current, current / SIGMA + rate / 6.0)
}

#[test]
fn saturated_series_resistor_matches_coupled_field_equations_and_reference_cancels() {
    let (material, source) = saturation(SIGMA);
    let problem = model(1, material);
    let field = SpatialNonlinearTerminalInduction::new(
        history(&mode(&problem, 0.0, 0.0), &source),
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (2, 1, Branch::Resistor { ohms: 0.4 }),
        ],
        transformers: vec![],
    };
    let dt = 0.2;
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    for reference in [0.001, 1.0, 1000.0] {
        let mut circuit = SpatialNonlinearTerminalCircuit::new(
            field.clone(),
            graph.clone(),
            vec![1, 0],
            vec![],
            initial(vec![], vec![], vec![0.0]),
            config(reference),
            &mut || false,
        )
        .unwrap();
        let mut old = 0.0;
        // Cross both source slope changes, then reverse the actual field.
        for next in [0.8_f64, 1.25, 0.4, -0.9, -1.3, -0.2] {
            let (h, current, voltage) = response(old, next, dt);
            let supply = voltage + 0.4 * current;
            let proposal = circuit
                .step_candidate(
                    &[supply],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    dt,
                    controls(nonlinear_controls()),
                    &mut || false,
                )
                .unwrap();
            let step = proposal.step();
            let change = energy(next.abs()) - energy(old.abs());
            let numerical = 2.0 * h * (next - old) - change;
            let joule = dt * (current * current / SIGMA + MASS * ((next - old) / dt).powi(2));
            near(proposal.circuit().field().potential_wb()[edge], next, 5e-8);
            near(step.terminal_voltage_v[0], voltage, 5e-8);
            near(step.terminal_current_a[0], current, 5e-8);
            near(step.terminal_current_a[1], -current, 5e-8);
            near(step.field.heat_j, joule, 5e-7);
            near(step.resistor_heat_j, dt * 0.4 * current * current, 5e-7);
            near(step.field.stored_after_j, energy(next.abs()), 5e-8);
            near(step.stored_change_j, change, 5e-8);
            near(step.numerical_dissipation_j, numerical, 5e-8);
            near(step.heat_j, joule + dt * 0.4 * current * current, 5e-7);
            near(step.supplied_work_j, dt * supply * current, 5e-7);
            near(step.energy_residual_j, 0.0, 2e-8);
            assert!(step.maximum_kirchhoff_current_defect_a <= 2e-8);
            assert!(step.report.as_ref().is_some_and(|r| r.converged));
            assert_eq!(step.graph_resistors.len(), 1);
            assert!(!step.reused_response);
            old = next;
            circuit = proposal.into_parts().0;
        }
    }
}

#[test]
fn original_inductor_flux_obeys_native_midpoint_law_beside_saturated_field_storage() {
    let (material, source) = saturation(SIGMA);
    let problem = model(1, material);
    let field = SpatialNonlinearTerminalInduction::new(
        history(&mode(&problem, 0.0, 0.0), &source),
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let inductance = 0.3;
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (
                2,
                1,
                Branch::Inductor {
                    henries: inductance,
                },
            ),
        ],
        transformers: vec![],
    };
    let mut flux = 0.15;
    let mut circuit = SpatialNonlinearTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![flux], vec![], vec![0.0]),
        config(0.7),
        &mut || false,
    )
    .unwrap();
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    let dt = 0.125;
    let mut old = 0.0;
    let initial_energy = flux * flux / (2.0 * inductance);
    let (mut physical_heat, mut numerical_loss, mut supplied) = (0.0, 0.0, 0.0);
    for drive in [12.0, 12.0, 24.0, 0.0, -8.0, -8.0, 0.0] {
        // I = (phi_old+phi_new)/(2L), phi_new-phi_old = dt*(E-U).
        let next = bracket(-1.95, 1.95, |a| {
            let (_, current, voltage) = response(old, a, dt);
            2.0 * (inductance * current - flux) - dt * (drive - voltage)
        });
        let (h, current, voltage) = response(old, next, dt);
        let next_flux = 2.0 * inductance * current - flux;
        let proposal = circuit
            .step_candidate(
                &[drive],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                controls(nonlinear_controls()),
                &mut || false,
            )
            .unwrap();
        let step = proposal.step();
        near(proposal.circuit().field().potential_wb()[edge], next, 5e-8);
        near(step.terminal_voltage_v[0], voltage, 5e-8);
        near(step.terminal_current_a[0], current, 5e-8);
        near(step.branches[1].current_a, current, 5e-8);
        near(
            proposal.circuit().native_state()[proposal.circuit().dae().flux_index[0]],
            next_flux,
            5e-8,
        );
        near(
            step.circuit_stored_change_j,
            current * (next_flux - flux),
            5e-8,
        );
        near(
            step.circuit_stored_after_j,
            next_flux * next_flux / (2.0 * inductance),
            5e-8,
        );
        near(step.field.stored_after_j, energy(next.abs()), 5e-8);
        near(
            step.numerical_dissipation_j,
            2.0 * h * (next - old) - energy(next.abs()) + energy(old.abs()),
            5e-8,
        );
        near(step.resistor_heat_j, 0.0, 1e-14);
        near(step.external_work_j, dt * drive * current, 5e-8);
        physical_heat += dt * (current * current / SIGMA + MASS * ((next - old) / dt).powi(2));
        numerical_loss += step.numerical_dissipation_j;
        supplied += dt * drive * current;
        near(
            step.circuit_stored_after_j + energy(next.abs()) + physical_heat + numerical_loss,
            initial_energy + supplied,
            3e-7,
        );
        old = next;
        flux = next_flux;
        circuit = proposal.into_parts().0;
    }
}

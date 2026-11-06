//! Massive conductor electrodes retain the actual reversible magnetic field.
use super::super::{near, preparation, source};
use crate::nonlinear::{energy, history, mode, model, nonlinear_controls, saturation};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        induction::nonlinear::terminals::SpatialNonlinearTerminalInduction,
        nonlinear::SpatialNonlinearVectorMagnetostatics,
    },
};

fn electrodes(problem: &SpatialNonlinearVectorMagnetostatics) -> Vec<Vec<usize>> {
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

#[test]
fn reversible_linear_limit_solves_total_current_and_inductive_voltage_with_exact_ohmic_heat() {
    let (material, source) = source(2.0, Some(5.0), None);
    let SpatialMagneticMaterial::Linear(material) = material else {
        unreachable!()
    };
    let problem = model(1, SpatialNonlinearMaterial::Linear(material));
    let source = source.unwrap();
    let mut state = SpatialNonlinearTerminalInduction::new(
        history(&mode(&problem, 0.0, 0.0), &source),
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap()
    .into_current_driven(preparation(), &mut || false)
    .unwrap();
    let (sigma, mass, stiffness, dt) = (5.0, 31.0 / 36.0, 2.0, 0.125);
    let mut old = 0.0_f64;
    for incoming in [1.7, 1.7, 0.0, -0.6, 0.0] {
        // int W_x=1/6 and int |W|²=1/5; eliminating the unknown
        // equipotential voltage gives mass=sigma*(1/5-1/36).
        let expected = (mass * old + dt * incoming / 6.0) / (mass + dt * stiffness);
        let rate = (expected - old) / dt;
        let voltage = incoming / sigma + rate / 6.0;
        let candidate = state
            .step_current_candidate(
                &[incoming, -incoming],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                nonlinear_controls(),
                &mut || false,
            )
            .unwrap();
        let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
        let step = candidate.step();
        near(candidate.state().potential_wb()[edge], expected, 3e-9);
        near(step.terminals[0].incoming_current_a, incoming, 3e-9);
        near(step.terminals[1].incoming_current_a, -incoming, 3e-9);
        near(step.terminals[0].voltage_to_reference_v, 0.0, 1e-14);
        near(step.terminals[1].voltage_to_reference_v, -voltage, 3e-9);
        near(step.terminal_work_j, dt * voltage * incoming, 3e-9);
        near(
            step.heat_j,
            dt * (incoming * incoming / sigma + mass * rate * rate),
            3e-9,
        );
        near(
            step.stored_after_j,
            0.5 * stiffness * expected * expected,
            3e-9,
        );
        near(
            step.numerical_dissipation_j,
            0.5 * stiffness * (expected - old).powi(2),
            3e-9,
        );
        near(step.impressed_work_j, 0.0, 1e-14);
        near(step.energy_residual_j, 0.0, 3e-9);
        old = expected;
        state = candidate.into_parts().0;
    }
}

#[test]
fn saturated_electrode_response_matches_piecewise_energy_and_voltage_current_duality() {
    let (material, source) = saturation(5.0);
    let problem = model(1, material);
    let mut state = SpatialNonlinearTerminalInduction::new(
        history(&mode(&problem, 0.0, 0.0), &source),
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let (sigma, dt) = (5.0, 0.2);
    let mut old = 0.0_f64;
    for (expected, h) in [(0.8_f64, 1.7_f64), (1.25, 4.5), (0.4, 0.4)] {
        let rate = (expected - old) / dt;
        let voltage = 12.0 * h / sigma + 6.0 * rate / 5.0;
        let incoming = sigma * (voltage - rate / 6.0);
        let candidate = state
            .step_voltage_candidate(
                &[voltage, 0.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                nonlinear_controls(),
                &mut || false,
            )
            .unwrap();
        let current = state
            .clone()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
        let dual = current
            .step_current_candidate(
                &[incoming, -incoming],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                nonlinear_controls(),
                &mut || false,
            )
            .unwrap();
        let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
        for value in [
            candidate.state().potential_wb()[edge],
            dual.state().potential_wb()[edge],
        ] {
            near(value, expected, 3e-8);
        }
        let step = candidate.step();
        let heat = dt * sigma * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
        near(step.terminals[0].incoming_current_a, incoming, 3e-8);
        near(
            dual.step().terminals[1].voltage_to_reference_v,
            -voltage,
            3e-8,
        );
        near(step.heat_j, heat, 3e-8);
        near(dual.step().heat_j, heat, 3e-8);
        near(step.stored_before_j, energy(old), 3e-8);
        near(step.stored_after_j, energy(expected), 3e-8);
        near(
            step.numerical_dissipation_j,
            2.0 * h * (expected - old) - energy(expected) + energy(old),
            3e-8,
        );
        near(step.terminal_work_j, dt * voltage * incoming, 3e-8);
        near(step.energy_residual_j, 0.0, 3e-8);
        assert!(
            step.maximum_charge_residual_a <= nonlinear_controls().induction.current_tolerance_a
        );
        // Publish the current-driven candidate through the voltage conversion:
        // the next step must retain exactly the same magnetic field and time.
        state = dual
            .into_parts()
            .0
            .into_voltage_driven(preparation(), &mut || false)
            .unwrap();
        old = expected;
    }
}

//! Voltage-driven solid electrodes close saturation and finite sigma(T) feedback.
#[path = "terminals/circuit.rs"]
mod circuit;

use super::*;
use fs_em::spatial::vector_potential::induction::nonlinear::terminals::SpatialNonlinearTerminalInduction;

fn bracket(mut lower: f64, mut upper: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lower) <= 0.0 && residual(upper) >= 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn finite_temperature_and_saturation_control_physical_electrode_current_and_release_stored_energy()
{
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
    let mut state = SpatialNonlinearTerminalInduction::new(
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
    let query = context(300.0);
    let cells = cells(&card, &curve, &query);
    let regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let mut heat = thermal(&curve, 0.0);
    let (mut old_a, mut old_t) = (0.0, 300.0);
    let dt = 0.1;
    for (index, voltage) in [15.0, 0.0].into_iter().enumerate() {
        let endpoint = |t: f64| {
            let sigma = (0.2 * (1.0 + 0.005 * (t - 300.0))).recip();
            // Six-tet cube: 2H(2a)+sigma*a'/5=sigma*U/6.
            let a = bracket(0.0, 1.999, |a| {
                2.0 * h(2.0 * a, true) + sigma * (a - old_a) / (5.0 * dt) - sigma * voltage / 6.0
            });
            let rate = (a - old_a) / dt;
            let incoming = sigma * (voltage - rate / 6.0);
            let joule = dt * sigma * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
            (a, incoming, joule, sigma)
        };
        let t = bracket(old_t, 500.0, |t| 4.0 * (t - old_t) - endpoint(t).2);
        let (a, incoming, joule, sigma) = endpoint(t);
        let before = (
            heat.cells().clone(),
            state.potential_wb().to_vec(),
            state.time_s(),
        );
        let mut polls = 0;
        let candidate = network
            .step_nonlinear_voltage_candidate(
                &heat,
                (index + 1) as f64 * dt,
                &state,
                &[voltage, 0.0],
                &[[0.0; 3]; 6],
                None,
                thermal_controls(),
                field_controls(),
                &[0.0],
                || {
                    polls += 1;
                    false
                },
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
            candidate.field().potential_wb()[diagonal(&problem)],
            a,
            3e-8,
        );
        near(step.terminals[0].incoming_current_a, incoming, 3e-7);
        near(step.terminal_work_j, dt * voltage * incoming, 3e-7);
        near(step.heat_j, joule, 3e-7);
        near(candidate.deposited_heat_j(), joule, 3e-7);
        near(report.coupled.thermal.stored_energy_change_j, joule, 3e-7);
        near(
            step.stored_change_j,
            energy(a, true) - energy(old_a, true),
            3e-8,
        );
        near(
            step.numerical_dissipation_j,
            2.0 * h(2.0 * a, true) * (a - old_a) - energy(a, true) + energy(old_a, true),
            3e-8,
        );
        near(
            report.total_energy_residual_j,
            0.0,
            thermal_controls().total_energy_tolerance_j,
        );
        for source in candidate.field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m().recip(), sigma, 3e-10);
            near(source.material().query().axes()["T"], t, 3e-8);
        }
        if index == 0 {
            assert!(2.0 * a > 1.0 && t > 310.0);
            assert!(incoming < endpoint(300.0).1);
            let mut reached = 0;
            assert!(matches!(
                network.step_nonlinear_voltage_candidate(
                    &heat,
                    dt,
                    &state,
                    &[voltage, 0.0],
                    &[[0.0; 3]; 6],
                    None,
                    thermal_controls(),
                    field_controls(),
                    &[0.0],
                    || {
                        reached += 1;
                        reached >= polls - 1
                    },
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
        } else {
            assert!(step.stored_change_j < 0.0 && joule > 0.0);
            near(step.terminal_work_j, 0.0, 1e-12);
        }
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb(), before.1);
        assert_eq!(state.time_s(), before.2);
        old_a = a;
        old_t = t;
        let (next_heat, next_state, _) = candidate.into_parts();
        heat = next_heat;
        state = next_state;
    }
}

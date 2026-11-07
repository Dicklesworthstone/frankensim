//! Authored terminal currents retain motional voltage and finite thermal feedback.
use super::*;

#[test]
fn terminal_controls_follow_motion_and_refreshed_resistance_on_one_thermal_clock() {
    for (current, speed) in [(2.0, 3.0), (2.0, -10.0), (-2.0, 3.0), (0.0, 3.0)] {
        for mode in 0..3 {
            let (a, ca) = material(1.0, 1.0, 0.002);
            let (b, cb) = material(3.0, 2.0, 0.0);
            let cards = [a, b];
            let curves = [ca, cb];
            let context = point(300.0);
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let cells = bindings(&cards, &curves, &context, &key);
            let mut thermal = initial(&curves, 999.0, None);
            let (mut problem, regions) = model(&cards, &context, 2);
            let network = SpatialConductionThermalNetwork {
                cells: &cells,
                links: &LINKS,
                regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let [positive, negative] = electrodes(&problem, 2.0);
            let terminals = [positive.as_slice(), negative.as_slice()];
            let currents = [current, -current];
            let controlled = [0];
            let voltage = [7.0];
            let current_only = [current];
            let response = SpatialConductanceControls {
                maximum_terminals: 2,
                maximum_field_solves: 1,
                maximum_response_values: 1000,
                maximum_dense_work: 1000,
                reference_voltage_v: 0.5,
                coefficient_tolerance_s: 1e-10,
            };
            let drive = match mode {
                0 => SpatialConductionDrive::total_current(&positive, &negative, current),
                1 => SpatialConductionDrive::TerminalCurrents {
                    terminals: &terminals,
                    currents_a: &currents,
                    response,
                },
                _ => SpatialConductionDrive::MixedTerminals {
                    terminals: &terminals,
                    current_terminals: &controlled,
                    voltages_v: &voltage,
                    currents_a: &current_only,
                    response,
                    mixed: SpatialMixedConductanceControls {
                        maximum_current_terminals: 1,
                        maximum_response_values: 100,
                        maximum_dense_work: 100,
                    },
                },
            };
            let velocities = vec![[0.0, speed, 0.0]; regions.len()];
            let fields = vec![[0.0, 0.0, 1.0]; regions.len()];
            let motion = SpatialConductionMotion {
                velocity_m_s: &velocities,
                flux_density_t: &fields,
            };
            let dt = 0.5;
            let mut theta = [0.0; 2];
            for step in 1..=2 {
                let before = thermal.energy_j();
                let next = network
                    .step_moving_candidate(
                        &thermal,
                        step as f64 * dt,
                        &problem,
                        drive,
                        motion,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                // Independent two-volume BE balance. Source-current control
                // fixes I²rho(T) regardless of the sign of mechanical work.
                let i2 = current * current;
                let aa = 1.5 - dt * i2 * 0.002;
                let rhs = [theta[0] + dt * i2, 2.0 * theta[1] + dt * i2 * 3.0];
                let det = 2.5 * aa - 0.25;
                theta = [
                    (2.5 * rhs[0] + 0.5 * rhs[1]) / det,
                    (0.5 * rhs[0] + aa * rhs[1]) / det,
                ];
                let rho = [1.0 + 0.002 * theta[0], 3.0];
                near(
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    300.0 + theta[0],
                );
                near(
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    300.0 + theta[1],
                );
                let report = next.report();
                let field = &report.electrical.field;
                let potential =
                    |faces: &[usize]| field.voltage_v()[problem.mesh().faces[faces[0]][0] as usize];
                let terminal_voltage = current * (rho[0] + rho[1]) + 2.0 * speed;
                near(
                    potential(&positive) - potential(&negative),
                    terminal_voltage,
                );
                if mode == 2 {
                    // The native response removes common bias and reports
                    // voltage relative to original terminal zero.
                    near(potential(&positive), 0.0);
                    near(potential(&negative), -terminal_voltage);
                }
                for (tet, density) in field.current_density_a_m2().iter().enumerate() {
                    near(density[0], -current);
                    near(density[1], 0.0);
                    near(density[2], 0.0);
                    near(field.cell_joule_power_w()[tet], i2 * rho[tet / 12] / 12.0);
                }
                let heat = dt * i2 * (rho[0] + rho[1]);
                near(
                    report.electrical.terminal_work_j,
                    dt * current * terminal_voltage,
                );
                near(
                    report.electrical.mechanical_work_j,
                    dt * 2.0 * current * speed,
                );
                near(report.electrical.heat_j, heat);
                near(report.electrical.net_supplied_work_j, heat);
                near(next.thermal().energy_j() - before, heat);
                near(report.total_energy_residual_j, 0.0);
                // Both the stationary response and actual motional solve count.
                assert!(report.circuit_trials >= 2);
                assert_eq!(report.circuit_trials % 2, 0);
                if current == 2.0 && speed < 0.0 {
                    assert!(report.electrical.terminal_work_j < 0.0);
                }
                let (next_thermal, next_problem, _) = next.into_parts();
                thermal = next_thermal;
                problem = next_problem;
            }
        }
    }
}

#[test]
fn moving_terminal_work_budget_and_late_cancellation_preserve_both_accepted_owners() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let (problem, regions) = model(&cards, &context, 1);
    let network = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let [positive, negative] = electrodes(&problem, 2.0);
    let drive = SpatialConductionDrive::total_current(&positive, &negative, 2.0);
    let velocities = vec![[0.0, 3.0, 0.0]; regions.len()];
    let fields = vec![[0.0, 0.0, 1.0]; regions.len()];
    let motion = SpatialConductionMotion {
        velocity_m_s: &velocities,
        flux_density_t: &fields,
    };
    let before_cells = thermal.cells().clone();
    let before_conductors = problem.conductors().to_vec();
    let mut bounded = config();
    bounded.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            drive,
            motion,
            limits(),
            field_controls(),
            bounded,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(matches!(
        network.step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            SpatialConductionDrive::Power {
                positive_faces: &positive,
                negative_faces: &negative,
                source: fs_em::conduction::ConductionPowerSource {
                    power_w: 1.0,
                    maximum_voltage_v: 100.0,
                    maximum_current_a: 100.0,
                },
                reference_voltage_v: 1.0,
            },
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::Source(
            fs_couple::circuit_enthalpy::conduction::ConductionThermalError::Invalid(
                "held terminal power is not admitted with prescribed motional EMF"
            )
        ))
    ));
    let calls = Cell::new(0usize);
    let accepted = network
        .step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            drive,
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    assert!(accepted.thermal().energy_j() > thermal.energy_j());
    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        network.step_moving_candidate(
            &thermal,
            0.5,
            &problem,
            drive,
            motion,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                calls.get() >= last
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before_cells);
    assert_eq!(problem.conductors(), &before_conductors);
}

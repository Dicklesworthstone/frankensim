//! G1/G3: nonlinear floating-electrode drive follows one original phase chart.
use super::*;
use fs_em::conduction::{
    axisymmetric::response::AxisymmetricConductionResponse,
    response::terminals::{
        ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits,
    },
};

fn split_phase_mesh() -> TriComplex2 {
    // Radius=1, total length=1, sweep=2: volume=mass=1 and terminal area=1.
    // Two disconnected half-length pieces share only their ideal middle disk.
    TriComplex2::from_indexed_triangles(
        "confined phase conductor with floating middle electrode",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.5],
            [0.0, 0.5],
            [0.0, 0.5],
            [1.0, 0.5],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(2.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn nonlinear_terminal_drive_retains_original_latent_and_sensible_history() {
    let terminals: [&[usize]; 3] = [&[0, 1], &[2, 3, 4, 5], &[6, 7]];
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let field_law = coordinate == ConductorMagnitude::ElectricField;
        for slope in [1.0, 3.0] {
            let card = source_card(slope, 1.0, coordinate);
            let caloric = curve(&card, 1.0);
            let mut thermal = initial(&caloric, 10.0);
            let context = field_point(305.0, 0.0);
            let sources = nonlinear(&card, &context).sources().to_vec();
            let mut problem = AxisymmetricConductionResponse::new(
                split_phase_mesh(),
                sources,
                limits(),
                &mut || false,
            )
            .unwrap();
            let cells = [binding(&card, &caloric, &context)];
            let layout = network(&cells);
            let mut old_h = 10.0;
            for interval in 1..=2 {
                // Eliminate rho and use the independent one-kilogram phase
                // chart: h-old_h=25*(1+s*f) for E, or 24+20*s*f for J.
                let (base_heat, slope_heat) = if field_law {
                    (25.0, 25.0 * slope)
                } else {
                    (24.0, 20.0 * slope)
                };
                let trial_fraction = (old_h + base_heat - 20.0) / (100.0 - slope_heat);
                let (h, t, fraction) = if trial_fraction <= 1.0 {
                    (20.0 + 100.0 * trial_fraction, 310.0, trial_fraction)
                } else {
                    let h = old_h + base_heat + slope_heat;
                    (h, 310.0 + (h - 120.0) / 2.0, 1.0)
                };
                let rho = if field_law {
                    (1.0 + slope * fraction) / 0.8
                } else {
                    1.2 + slope * fraction
                };
                let drive = ConductionTerminalDrive {
                    terminals: &terminals,
                    control: if interval == 1 {
                        ConductionTerminalControl::Currents(&[-2.0, 0.0, 2.0])
                    } else {
                        ConductionTerminalControl::Mixed {
                            current_terminals: &[2, 1],
                            currents_a: &[2.0, 0.0],
                            voltages_v: &[1e250],
                        }
                    },
                    limits: ConductionTerminalLimits {
                        maximum_terminals: 3,
                        maximum_workspace_values: 4096,
                    },
                };
                let next = layout
                    .step_axisymmetric_response_terminal_candidate(
                        &thermal,
                        5.0 * f64::from(interval),
                        &problem,
                        drive,
                        &[],
                        limits(),
                        controls(),
                        response_controls(),
                        config(),
                        &[h - old_h],
                        || false,
                    )
                    .unwrap();
                near(next.thermal().energy_j(), h);
                near(next.thermal().temperature_k(CellId(7)).unwrap(), t);
                near(
                    next.report().source.cells[0].phase.liquid_mass_fraction(),
                    fraction,
                );
                let step = &next.report().electrical;
                let terminal = step.terminals.as_ref().unwrap();
                near(terminal.voltage_v[0], 0.0);
                near(terminal.voltage_v[1], rho);
                near(terminal.voltage_v[2], 2.0 * rho);
                near(terminal.current_a[0], -2.0);
                near(terminal.current_a[1], 0.0);
                near(terminal.current_a[2], 2.0);
                near(step.heat_j, h - old_h);
                near(step.supplied_work_j, h - old_h);
                near(step.supplied_work_j, 5.0 * terminal.supplied_power_w);
                near(next.report().deposited_joule_j, h - old_h);
                for (triangle, conductor) in next.problem().conductors().iter().enumerate() {
                    near(conductor.material().query().axes()["T"], t);
                    near(conductor.material().query().axes()[FRACTION], fraction);
                    near(conductor.material().query().axes()["field"], 2.0 * rho);
                    near(conductor.material().query().axes()["current"], 2.0);
                    assert_eq!(conductor.liquid_mass_fraction_axis(), Some(FRACTION));
                    assert_eq!(
                        conductor.source_property().answer().receipt.selected,
                        problem.conductors()[triangle]
                            .source_property()
                            .answer()
                            .receipt
                            .selected
                    );
                    near(step.field.current_density_a_m2[triangle][1], -2.0);
                    near(
                        step.field.cell_power_w[triangle],
                        4.0 * rho * problem.mesh().face_measure(triangle).unwrap(),
                    );
                }
                assert!(
                    next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
                );
                // Mixing the new thermal history with the stale source owner
                // must refuse before publishing another interval.
                assert!(
                    layout
                        .step_axisymmetric_response_terminal_candidate(
                            next.thermal(),
                            5.0 * f64::from(interval + 1),
                            &problem,
                            drive,
                            &[],
                            limits(),
                            controls(),
                            response_controls(),
                            config(),
                            &[1.0],
                            || false,
                        )
                        .is_err()
                );
                old_h = h;
                (thermal, problem, _) = next.into_parts();
            }
            near(context.axes()[FRACTION], 0.0);
            near(context.axes()["field"], 0.0);
        }
    }
}

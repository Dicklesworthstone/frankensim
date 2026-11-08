//! Source-density gravity in the original pressure-selected thermal gas owner.
use super::*;
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::SeriesGravity;
use fs_airflow::graph::thermal::medium::steady_gas::pressure::{
    GasGravityBudget, PressureDrivenGasSolution, solve_pressure_gravity_driven_material_gas,
};
use fs_material::fluid_transport::FluidTransport;

fn gravity_budget() -> GasGravityBudget {
    GasGravityBudget {
        flow: pressure_budget(),
        maximum_neglected_potential_power_w: 0.1,
    }
}

fn solve(
    cx: &Cx<'_>,
    p: &SteadyAxialGas,
    transport: &[GasCellTransport],
    drive: f64,
    path: &SeriesGravity,
    budget: GasGravityBudget,
) -> Result<PressureDrivenGasSolution, GasPressureError> {
    solve_pressure_gravity_driven_material_gas(
        cx,
        p,
        transport,
        Pressure::new(drive),
        path,
        controls(),
        feedback(),
        budget,
    )
}

#[test]
fn zero_gravity_preserves_the_original_pressure_solution() {
    with_cx(|cx, _| {
        let mut p = problem(false);
        for cell in &mut p.cells {
            cell.wall_temperature_k = 400.0;
        }
        let gas = p
            .source
            .evaluate(&source::point(400.0, p.static_pressure_pa))
            .unwrap();
        let drive = 64.0 * LENGTH * gas.state().dynamic_viscosity * FLOW
            / (gas.state().density * area() * DIAMETER * DIAMETER);
        let original = solve_pressure_driven_material_gas(
            cx,
            &p,
            &laws(),
            Pressure::new(drive),
            controls(),
            feedback(),
            pressure_budget(),
        )
        .unwrap();
        let path = SeriesGravity::try_new(
            [0.0; 3],
            vec![[0.0; 3], [LENGTH, 0.0, 0.0], [2.0 * LENGTH, 0.0, 0.0]],
        )
        .unwrap();
        let with_gravity = solve(
            cx,
            &p,
            &laws(),
            drive,
            &path,
            GasGravityBudget {
                maximum_neglected_potential_power_w: 0.0,
                ..gravity_budget()
            },
        )
        .unwrap();
        assert_eq!(original.mass_flow_kg_per_s, with_gravity.mass_flow_kg_per_s);
        assert_eq!(
            original.pressure_residual_pa,
            with_gravity.pressure_residual_pa
        );
        assert_eq!(
            original.common_pressure_bound_pa,
            with_gravity.common_pressure_bound_pa
        );
        assert_eq!(original.flow_evaluations, with_gravity.flow_evaluations);
        assert_eq!(original.substeps, with_gravity.substeps);
        assert_eq!(original.gas.thermal, with_gravity.gas.thermal);
        assert!(original.gravity.is_none());
        let observed = with_gravity.gravity.unwrap();
        assert_eq!(observed.geometry, path);
        assert_eq!(observed.potential_power_change_w, 0.0);
        assert_eq!(observed.neglected_potential_power_bound_w, 0.0);
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(0.0),
                controls(),
                feedback(),
                pressure_budget()
            ),
            Err(GasPressureError::Invalid(_))
        ));
    });
}

#[test]
fn downhill_pressure_flow_scales_with_actual_density_and_signed_drive() {
    with_cx(|cx, _| {
        let mut p = problem(false);
        for cell in &mut p.cells {
            cell.wall_temperature_k = 400.0;
        }
        let g = 9.81;
        let path = SeriesGravity::try_new(
            [0.0, 0.0, -g],
            vec![[0.0; 3], [0.0, 0.0, -LENGTH], [0.0, 0.0, -2.0 * LENGTH]],
        )
        .unwrap();
        let mut reference_mass = 0.0;
        for pressure_scale in [1.0, 2.0] {
            p.static_pressure_pa = 100000.0 * pressure_scale;
            let gas = p
                .source
                .evaluate(&source::point(400.0, p.static_pressure_pa))
                .unwrap();
            let rho = gas.state().density;
            let head = rho * g * 2.0 * LENGTH;
            for drive_fraction in [0.0, -0.25] {
                let drive = drive_fraction * head;
                // Independent isothermal Poiseuille closure: hydrostatic head
                // scales with rho and pressure loss per unit mass with 1/rho.
                let expected = (head + drive) * rho * area() * DIAMETER * DIAMETER
                    / (64.0 * LENGTH * gas.state().dynamic_viscosity);
                let mut budget = gravity_budget();
                budget.flow.pressure.mass_flow_bracket_kg_per_s =
                    (0.75 * expected, 1.25 * expected);
                let solution = solve(cx, &p, &laws(), drive, &path, budget).unwrap();
                near(solution.mass_flow_kg_per_s, expected, 1e-6 * expected);
                let report = solution.gravity.as_ref().unwrap();
                near(report.hydrostatic.pressure_gain_pa(), head, 1e-10);
                near(
                    report.potential_power_change_w,
                    -solution.mass_flow_kg_per_s * g * 2.0 * LENGTH,
                    1e-13,
                );
                near(
                    report.neglected_potential_power_bound_w,
                    -report.potential_power_change_w,
                    1e-13,
                );
                assert_eq!(solution.driving_pressure_pa, drive);
                if drive_fraction == 0.0 {
                    if pressure_scale == 1.0 {
                        reference_mass = solution.mass_flow_kg_per_s;
                    } else {
                        near(
                            solution.mass_flow_kg_per_s,
                            4.0 * reference_mass,
                            1e-6 * expected,
                        );
                    }
                }
            }
        }
    });
}

#[test]
fn heated_stack_uses_an_explicit_source_column_and_has_no_neutral_flow() {
    with_cx(|cx, gate| {
        let mut p = problem(true);
        p.inlet = primitive(300.0, 0.75);
        for cell in &mut p.cells {
            cell.initial_temperature_k = 300.0;
            cell.wall_temperature_k = 480.0;
        }
        let g = 9.81;
        let path = SeriesGravity::try_new(
            [0.0, 0.0, -g],
            vec![[0.0; 3], [0.0, 0.0, LENGTH], [0.0, 0.0, 2.0 * LENGTH]],
        )
        .unwrap();
        let ambient = p
            .source
            .evaluate(&source::point(300.0, p.static_pressure_pa))
            .unwrap();
        let outside = path
            .evaluate(
                &[ambient.transport_state().unwrap(); 2],
                MassFlowRate::new(0.0),
                gate,
            )
            .unwrap();
        let drive = -outside.pressure_gain_pa();
        let solution = solve(cx, &p, &laws(), drive, &path, gravity_budget()).unwrap();
        let report = solution.gravity.as_ref().unwrap();
        let density_sum: f64 = solution
            .gas
            .cells
            .iter()
            .map(|cell| cell.gas.state().density)
            .sum();
        let buoyancy = g * LENGTH * (2.0 * ambient.state().density - density_sum);
        assert!(buoyancy > 0.0);
        assert!(
            solution
                .gas
                .cells
                .iter()
                .all(|cell| cell.gas.state().density < ambient.state().density)
        );
        near(
            solution.gas.total_pressure_drop_pa + solution.acceleration_pressure_pa,
            buoyancy,
            2e-7,
        );
        near(
            report.hydrostatic.pressure_gain_pa(),
            -g * LENGTH * density_sum,
            1e-11,
        );
        near(
            report.potential_power_change_w,
            solution.mass_flow_kg_per_s * g * 2.0 * LENGTH,
            1e-13,
        );
        assert!(solution.substeps <= gravity_budget().flow.maximum_substeps);
        // A source-identical outside column and unheated gas cancel hydrostatic
        // head; positive friction must not manufacture a natural-circulation root.
        for cell in &mut p.cells {
            cell.wall_temperature_k = 300.0;
        }
        assert!(matches!(
            solve(cx, &p, &laws(), drive, &path, gravity_budget()),
            Err(GasPressureError::NotBracketed { .. })
        ));
    });
}

#[test]
fn closed_path_buoyancy_has_zero_net_work_but_regional_energy_and_pressure_limits() {
    with_cx(|cx, gate| {
        let mut p = problem(true);
        let cell = p.cells[0].clone();
        p.cells = (0..4)
            .map(|i| AxialGasCell {
                end_m: (i + 1) as f64 * LENGTH,
                wall_temperature_k: if i < 2 { 480.0 } else { 300.0 },
                ..cell.clone()
            })
            .collect();
        let transport = [laws()[0]; 4];
        let g = 9.81;
        let vertices = vec![
            [0.0; 3],
            [0.0, 0.0, LENGTH],
            [LENGTH, 0.0, LENGTH],
            [LENGTH, 0.0, 0.0],
            [0.0; 3],
        ];
        let path = SeriesGravity::try_new([0.0, 0.0, -g], vertices.clone()).unwrap();
        let solution = solve(cx, &p, &transport, 0.0, &path, gravity_budget()).unwrap();
        let report = solution.gravity.as_ref().unwrap();
        let hot = solution.gas.cells[0].gas.state().density;
        let cold = solution.gas.cells[2].gas.state().density;
        let head = (cold - hot) * g * LENGTH;
        assert!(head > 0.0);
        near(report.hydrostatic.pressure_gain_pa(), head, 1e-11);
        near(
            solution.gas.total_pressure_drop_pa + solution.acceleration_pressure_pa,
            head,
            2e-7,
        );
        assert_eq!(report.hydrostatic.gravitational_power_w(), 0.0);
        assert_eq!(report.potential_power_change_w, 0.0);
        near(
            report.neglected_potential_power_bound_w,
            2.0 * solution.mass_flow_kg_per_s * g * LENGTH,
            1e-13,
        );
        near(
            report.hydrostatic_pressure_bound_pa,
            (cold + hot) * g * LENGTH,
            1e-11,
        );
        assert!(matches!(
            solve(
                cx,
                &p,
                &transport,
                0.0,
                &path,
                GasGravityBudget {
                    maximum_neglected_potential_power_w: 0.0,
                    ..gravity_budget()
                }
            ),
            Err(GasPressureError::NeglectedPotentialPower { .. })
        ));

        // Admit the first trial's friction while refusing its much larger
        // absolute hydrostatic excursions; a small NET head must not hide them.
        let mut lower = p.clone();
        lower.inlet_mass_flow_kg_s = pressure_budget().pressure.mass_flow_bracket_kg_per_s.0;
        let first =
            solve_material_gas_transport(cx, &lower, &transport, controls(), feedback()).unwrap();
        let inlet = first.thermal.field().state_at(0.0).unwrap();
        let mass_flux = lower.inlet_mass_flow_kg_s / area();
        let inlet_speed = mass_flux / inlet.state().density;
        let maximum_acceleration = first
            .cells
            .iter()
            .map(|cell| {
                (pressure_budget().momentum_flux_factor
                    * mass_flux
                    * (mass_flux / cell.gas.state().density - inlet_speed))
                    .abs()
            })
            .fold(0.0, f64::max);
        let absolute_head =
            g * LENGTH * (first.cells[0].gas.state().density + first.cells[2].gas.state().density);
        let allowance = first.total_pressure_drop_pa + maximum_acceleration + 0.5 * absolute_head;
        assert!(matches!(
            solve_pressure_gravity_driven_material_gas(
                cx,
                &p,
                &transport,
                Pressure::new(0.0),
                &path,
                controls(),
                GasTransportFeedback {
                    maximum_pressure_drop_fraction: allowance / p.static_pressure_pa,
                    ..feedback()
                },
                gravity_budget()
            ),
            Err(GasPressureError::PressureVariation { .. })
        ));
        let mut too_tall = vertices;
        too_tall[1][2] = 2.0 * LENGTH;
        let incompatible = SeriesGravity::try_new([0.0, 0.0, -g], too_tall).unwrap();
        assert!(matches!(
            solve(cx, &p, &transport, 0.0, &incompatible, gravity_budget()),
            Err(GasPressureError::Invalid(_))
        ));
        let saved = p.clone();
        gate.request();
        assert!(matches!(
            solve(cx, &p, &transport, 0.0, &path, gravity_budget()),
            Err(GasPressureError::Transport(
                MaterialSteadyGasError::Cancelled
            ))
        ));
        assert_eq!(p, saved);
    });
}

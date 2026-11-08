//! G1/G3: original-law pressure split, density and conservative manifold enthalpy.
use super::*;
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::SeriesGravity;
use fs_airflow::graph::thermal::medium::steady_gas::pressure::GasGravityBudget;
use fs_airflow::graph::thermal::medium::steady_gas::pressure::parallel::*;

fn budget() -> ParallelGasBudget {
    ParallelGasBudget {
        maximum_flow_evaluations: 128,
        maximum_substeps: 2_000_000,
        maximum_neglected_hydraulic_power_w: 2.0,
        maximum_neglected_kinetic_power_w: 2.0,
        maximum_neglected_potential_power_w: 0.2,
    }
}
fn branch<'a>(
    name: &'a str,
    p: &'a SteadyAxialGas,
    laws: &'a [GasCellTransport],
) -> ParallelGasBranch<'a> {
    ParallelGasBranch {
        name,
        problem: p,
        transport: laws,
        gravity: None,
        budget: GasGravityBudget {
            flow: pressure_budget(),
            maximum_neglected_potential_power_w: 0.0,
        },
    }
}
fn solve(
    cx: &Cx<'_>,
    branches: &[ParallelGasBranch<'_>],
    drive: f64,
    work: ParallelGasBudget,
) -> Result<ParallelGasSolution, ParallelGasError> {
    solve_parallel_pressure_driven_material_gas(
        cx,
        branches,
        Pressure::new(drive),
        controls(),
        feedback(),
        work,
    )
}
fn isothermal() -> SteadyAxialGas {
    let mut p = problem(false);
    for cell in &mut p.cells {
        cell.wall_temperature_k = 400.0;
    }
    p
}
fn isothermal_drive(p: &SteadyAxialGas) -> f64 {
    let gas = p
        .source
        .evaluate(&source::point(400.0, p.static_pressure_pa))
        .unwrap();
    64.0 * LENGTH * gas.state().dynamic_viscosity * FLOW
        / (gas.state().density * area() * DIAMETER * DIAMETER)
}

#[test]
fn parallel_mass_split_has_independent_poiseuille_and_density_oracles() {
    with_cx(|cx, _| {
        let short = isothermal();
        let mut long = short.clone();
        for cell in &mut long.cells {
            cell.end_m *= 2.0;
            cell.volume_m3 *= 2.0;
            cell.wall_area_m2 *= 2.0;
        }
        let short_laws = laws();
        let mut long_laws = laws();
        for law in &mut long_laws {
            law.settings.length_over_hydraulic_diameter = Some(2.0 * LENGTH / DIAMETER);
            law.resistance = MaterialFlowResistance::try_new(
                64.0 * LENGTH / (area() * DIAMETER * DIAMETER),
                0.0,
                Area::new(area()),
                Length::new(DIAMETER),
                (0.0, 2300.0),
            )
            .unwrap();
        }
        let drive = isothermal_drive(&short);
        let paths = [
            branch("short", &short, &short_laws),
            branch("long", &long, &long_laws),
        ];
        let value = solve(cx, &paths, drive, budget()).unwrap();
        near(
            value.branches[0].solution.mass_flow_kg_per_s,
            FLOW,
            1e-6 * FLOW,
        );
        near(
            value.branches[1].solution.mass_flow_kg_per_s,
            0.5 * FLOW,
            1e-6 * FLOW,
        );
        near(value.total_mass_flow_kg_per_s, 1.5 * FLOW, 2e-6 * FLOW);
        near(value.mixed_outlet.state().temperature, 400.0, 1e-9);
        near(
            value.inlet_volume_flow_m3_per_s,
            value.outlet_volume_flow_m3_per_s,
            1e-12,
        );
        assert_eq!(
            value.flow_evaluations,
            value
                .branches
                .iter()
                .map(|b| b.solution.flow_evaluations)
                .sum::<usize>()
        );
        assert_eq!(
            value.substeps,
            value
                .branches
                .iter()
                .map(|b| b.solution.substeps)
                .sum::<usize>()
        );
        // One branch retains the original pressure-selected thermal solution.
        let single = solve(cx, &paths[..1], drive, budget()).unwrap();
        let native = solve_pressure_driven_material_gas(
            cx,
            &short,
            &short_laws,
            Pressure::new(drive),
            controls(),
            feedback(),
            pressure_budget(),
        )
        .unwrap();
        assert_eq!(
            single.branches[0].solution.mass_flow_kg_per_s,
            native.mass_flow_kg_per_s
        );
        assert_eq!(single.branches[0].solution.gas.thermal, native.gas.thermal);
        assert_eq!(single.flow_evaluations, native.flow_evaluations);
        // At fixed temperature/viscosity and imposed dp, mass scales with rho
        // while inlet/outlet VOLUME flow remains unchanged.
        let mut dense_short = short.clone();
        let mut dense_long = long.clone();
        dense_short.static_pressure_pa *= 2.0;
        dense_long.static_pressure_pa *= 2.0;
        let dense = solve(
            cx,
            &[
                branch("short", &dense_short, &short_laws),
                branch("long", &dense_long, &long_laws),
            ],
            drive,
            budget(),
        )
        .unwrap();
        near(
            dense.total_mass_flow_kg_per_s,
            2.0 * value.total_mass_flow_kg_per_s,
            2e-6 * FLOW,
        );
        near(
            dense.inlet_volume_flow_m3_per_s,
            value.inlet_volume_flow_m3_per_s,
            2e-10,
        );
        // Isothermal branches have zero regional kinetic change, yet finite
        // inlet/outlet kinetic flux cannot disappear from the manifold allowance.
        assert_eq!(
            single.branches[0].solution.neglected_kinetic_power_bound_w,
            0.0
        );
        assert!(single.manifold_kinetic_power_bound_w > 0.0);
        assert!(matches!(
            solve(
                cx,
                &paths[..1],
                drive,
                ParallelGasBudget {
                    maximum_neglected_kinetic_power_w: 0.0,
                    ..budget()
                }
            ),
            Err(ParallelGasError::OmittedPower { .. })
        ));
    });
}

#[test]
fn thermal_density_redistributes_flow_and_outlets_mix_original_nonlinear_enthalpy() {
    with_cx(|cx, _| {
        let mut hot = problem(true);
        hot.source = source::model_transport_with_heat_capacity_slope(0.75, 0.01);
        let mut cold = hot.clone();
        for cell in &mut hot.cells {
            cell.wall_temperature_k = 480.0;
        }
        for cell in &mut cold.cells {
            cell.wall_temperature_k = 300.0;
        }
        let transport = laws();
        let value = solve(
            cx,
            &[
                branch("hot", &hot, &transport),
                branch("cold", &cold, &transport),
            ],
            isothermal_drive(&hot),
            budget(),
        )
        .unwrap();
        let hot_flow = &value.branches[0].solution;
        let cold_flow = &value.branches[1].solution;
        assert!(cold_flow.mass_flow_kg_per_s > hot_flow.mass_flow_kg_per_s);
        assert!(
            cold_flow.gas.cells.last().unwrap().gas.state().density
                > hot_flow.gas.cells.last().unwrap().gas.state().density
        );
        // Independent primitive NASA polynomial: Cp_molar/R = 3.25 + .01*T.
        // The mass-weighted source enthalpy has a quadratic T term, so mass-
        // weighted temperature is also physically wrong.
        let mut enthalpy_coordinate = 0.0;
        let mut temperature_average = 0.0;
        let mut gain = 0.0;
        for branch in &value.branches {
            let flow = &branch.solution;
            let m = flow.mass_flow_kg_per_s;
            let t = flow.gas.cells.last().unwrap().gas.state().temperature;
            enthalpy_coordinate += m * (3.25 * t + 0.005 * t * t);
            temperature_average += m * t;
            gain += m
                * (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.022)
                * (3.25 * (t - 400.0) + 0.005 * (t - 400.0) * (t + 400.0));
            near(flow.pressure_residual_pa, 0.0, 1.1e-7);
        }
        let h = enthalpy_coordinate / value.total_mass_flow_kg_per_s;
        let expected = (-3.25 + (3.25_f64 * 3.25 + 0.02 * h).sqrt()) / 0.01;
        near(value.mixed_outlet.state().temperature, expected, 2e-8);
        assert!((expected - temperature_average / value.total_mass_flow_kg_per_s).abs() > 0.5);
        near(value.total_enthalpy_gain_w, gain, 1e-9);
        near(value.total_wall_heat_w, gain, 5e-7);
        assert!(
            value.mixing_enthalpy_residual_w.abs()
                <= value.total_mass_flow_kg_per_s * controls().caloric.energy_tolerance_j_kg
        );
        near(
            value.outlet_volume_flow_m3_per_s,
            value.total_mass_flow_kg_per_s / value.mixed_outlet.state().density,
            1e-15,
        );
        assert_eq!(value.mixed_outlet.source_packs(), hot.source.source_packs());
        assert!(value.neglected_kinetic_power_bound_w >= value.manifold_kinetic_power_bound_w);
        assert!(value.substeps <= budget().maximum_substeps);
    });
}

#[test]
fn parallel_gravity_uses_shared_physical_ports_and_real_branch_density() {
    with_cx(|cx, _| {
        let mut hot = problem(true);
        hot.inlet = primitive(300.0, 0.75);
        for cell in &mut hot.cells {
            cell.initial_temperature_k = 300.0;
            cell.wall_temperature_k = 480.0;
        }
        let mut warm = hot.clone();
        for cell in &mut warm.cells {
            cell.wall_temperature_k = 400.0;
        }
        let path = SeriesGravity::try_new(
            [0.0, 0.0, -9.81],
            vec![[0.0; 3], [0.0, 0.0, LENGTH], [0.0, 0.0, 2.0 * LENGTH]],
        )
        .unwrap();
        let ambient = hot
            .source
            .evaluate(&source::point(300.0, hot.static_pressure_pa))
            .unwrap();
        let drive = 2.0 * LENGTH * 9.81 * ambient.state().density;
        let transport = laws();
        let paths = [("hot", &hot), ("warm", &warm)].map(|(name, p)| {
            let mut b = branch(name, p, &transport);
            b.gravity = Some(&path);
            b.budget.maximum_neglected_potential_power_w = 0.1;
            b
        });
        let value = solve(cx, &paths, drive, budget()).unwrap();
        for branch in &value.branches {
            let flow = &branch.solution;
            let gravity = flow.gravity.as_ref().unwrap();
            assert_eq!(gravity.geometry, path);
            let density: f64 = flow.gas.cells.iter().map(|c| c.gas.state().density).sum();
            let buoyancy = 9.81 * LENGTH * (2.0 * ambient.state().density - density);
            near(
                flow.gas.total_pressure_drop_pa + flow.acceleration_pressure_pa,
                buoyancy,
                2e-7,
            );
        }
        near(
            value.neglected_potential_power_bound_w,
            value.total_mass_flow_kg_per_s * 9.81 * 2.0 * LENGTH,
            1e-13,
        );
        let moved = SeriesGravity::try_new(
            [0.0, 0.0, -9.81],
            vec![
                [1.0, 0.0, 0.0],
                [1.0, 0.0, LENGTH],
                [1.0, 0.0, 2.0 * LENGTH],
            ],
        )
        .unwrap();
        let mismatch = [
            paths[0],
            ParallelGasBranch {
                gravity: Some(&moved),
                ..paths[1]
            },
        ];
        assert!(matches!(
            solve(cx, &mismatch, drive, budget()),
            Err(ParallelGasError::Invalid(_))
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                drive,
                ParallelGasBudget {
                    maximum_neglected_potential_power_w: 0.0,
                    ..budget()
                }
            ),
            Err(ParallelGasError::OmittedPower { .. })
        ));
    });
}

#[test]
fn incompatible_sources_and_shared_work_fail_without_publishing_partial_branches() {
    with_cx(|cx, gate| {
        let hot = problem(true);
        let cold = hot.clone();
        let transport = laws();
        let paths = [
            branch("first", &hot, &transport),
            branch("second", &cold, &transport),
        ];
        let drive = isothermal_drive(&hot);
        let value = solve(cx, &paths, drive, budget()).unwrap();
        let first = &value.branches[0].solution;
        assert!(first.flow_evaluations >= 2 && first.substeps > 0);
        assert!(matches!(
            solve(
                cx,
                &paths,
                drive,
                ParallelGasBudget {
                    maximum_flow_evaluations: first.flow_evaluations,
                    ..budget()
                }
            ),
            Err(ParallelGasError::WorkBudget { .. })
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                drive,
                ParallelGasBudget {
                    maximum_substeps: first.substeps,
                    ..budget()
                }
            ),
            Err(ParallelGasError::WorkBudget { .. })
        ));
        let mut different = cold.clone();
        different.source = source::model_transport(0.75, false);
        let mismatched = [
            paths[0],
            ParallelGasBranch {
                problem: &different,
                ..paths[1]
            },
        ];
        assert!(matches!(
            solve(cx, &mismatched, drive, budget()),
            Err(ParallelGasError::Invalid(_))
        ));
        assert!(matches!(
            solve(cx, &[], drive, budget()),
            Err(ParallelGasError::Invalid(_))
        ));
        let duplicate = [
            paths[0],
            ParallelGasBranch {
                name: "first",
                ..paths[1]
            },
        ];
        assert!(matches!(
            solve(cx, &duplicate, drive, budget()),
            Err(ParallelGasError::Invalid(_))
        ));
        let mut invalid_budget = paths[0];
        invalid_budget.budget.flow.pressure.max_flow_evaluations = 129;
        // The shared allowance is 128: clipping must not heal an invalid 129.
        assert!(matches!(
            solve(cx, &[invalid_budget], drive, budget()),
            Err(ParallelGasError::Invalid(_))
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                drive,
                ParallelGasBudget {
                    maximum_neglected_hydraulic_power_w: 0.0,
                    ..budget()
                }
            ),
            Err(ParallelGasError::OmittedPower { .. })
        ));
        let saved = hot.clone();
        gate.request();
        assert!(matches!(
            solve(cx, &paths, drive, budget()),
            Err(ParallelGasError::Cancelled)
        ));
        assert_eq!(hot, saved);
    });
}

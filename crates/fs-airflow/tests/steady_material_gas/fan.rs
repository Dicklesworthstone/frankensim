//! G1/G3: native fan affinity selects original thermal/density branch flow.
use super::*;
use fs_airflow::graph::thermal::medium::{
    FrozenFlowMedium,
    fan::{FanGasDomain, MaterialFanError},
    steady_gas::pressure::{
        GasGravityBudget,
        parallel::{fan::*, *},
    },
};
use fs_airflow::{
    AirflowError, FanArrangement, FanBank, FanCurve, FanPoint, SourceProvenance, ToleranceBasis,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::VolumetricFlowRate;

fn work() -> ParallelGasBudget {
    ParallelGasBudget {
        maximum_flow_evaluations: 4096,
        maximum_substeps: 16_000_000,
        maximum_neglected_hydraulic_power_w: 2.0,
        maximum_neglected_kinetic_power_w: 2.0,
        maximum_neglected_potential_power_w: 0.2,
    }
}
fn match_controls() -> GasFanControls {
    GasFanControls {
        pressure_bracket_pa: (2.0, 5.5),
        absolute_pressure_tolerance_pa: 1e-6,
        relative_pressure_tolerance: 1e-8,
        maximum_network_evaluations: 48,
        maximum_neglected_shaft_power_w: 0.1,
    }
}
fn branch<'a>(
    name: &'a str,
    p: &'a SteadyAxialGas,
    transport: &'a [GasCellTransport],
) -> ParallelGasBranch<'a> {
    ParallelGasBranch {
        name,
        problem: p,
        transport,
        gravity: None,
        budget: GasGravityBudget {
            flow: pressure_budget(),
            maximum_neglected_potential_power_w: 0.0,
        },
    }
}
fn bank(speed: f64, stall: f64) -> FanBank {
    let curve = FanCurve::new(
        "manufactured static fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(6.0)),
            FanPoint::new(VolumetricFlowRate::new(1e-4), Pressure::new(0.0)),
        ],
        SourceProvenance::new(
            "manufactured linear regression fixture; not measured",
            "gas-fan-linear-fixture",
        ),
        0.05,
        ToleranceBasis::EngineeringAllowance,
        VolumetricFlowRate::new(stall),
        (0.5, 2.0),
    )
    .unwrap();
    FanBank::new(curve, 1, FanArrangement::Series, speed).unwrap()
}
fn source_fan<'a>(
    reference_bank: &'a FanBank,
    reference: &'a ResolvedThermochemicalGasState,
) -> GasFanSource<'a> {
    GasFanSource {
        reference_bank,
        reference_medium: FrozenFlowMedium::ThermochemicalGas(reference),
        gas_domain: FanGasDomain {
            density_ratio: [0.5, 3.0],
            reynolds_ratio: [0.3, 4.0],
            maximum_pressure_fraction: 0.01,
        },
        minimum_static_efficiency: 0.3,
    }
}
fn isothermal() -> SteadyAxialGas {
    let mut p = problem(false);
    for cell in &mut p.cells {
        cell.wall_temperature_k = 400.0;
    }
    p
}
fn heated() -> (SteadyAxialGas, SteadyAxialGas) {
    let mut hot = problem(true);
    hot.source = source::model_transport_with_heat_capacity_slope(0.75, 0.01);
    let mut cold = hot.clone();
    for cell in &mut hot.cells {
        cell.wall_temperature_k = 480.0;
    }
    for cell in &mut cold.cells {
        cell.wall_temperature_k = 300.0;
    }
    (hot, cold)
}
fn solve<'a>(
    cx: &Cx<'_>,
    paths: &[ParallelGasBranch<'_>],
    source: GasFanSource<'a>,
    budget: ParallelGasBudget,
    limits: GasFanControls,
) -> Result<GasFanOperatingPoint<'a>, GasFanError> {
    solve_parallel_material_gas_with_fan(cx, paths, source, controls(), feedback(), budget, limits)
}

#[test]
fn fan_pipe_intersection_preserves_actual_density_speed_and_source() {
    with_cx(|cx, _| {
        let original = isothermal();
        let reference = original
            .source
            .evaluate(&source::point(400.0, 100000.0))
            .unwrap();
        let transport = laws();
        for (density_ratio, speed) in [(1.0, 1.0), (2.0, 1.0), (1.0, 0.75)] {
            let mut p = original.clone();
            p.static_pressure_pa *= density_ratio;
            let bank = bank(speed, 0.0);
            let retained = bank.clone();
            let paths = [
                branch("left", &p, &transport),
                branch("right", &p, &transport),
            ];
            let result = solve(
                cx,
                &paths,
                source_fan(&bank, &reference),
                work(),
                match_controls(),
            )
            .unwrap();
            let inlet = result.network.inlet.state();
            // Two identical parallel paths, each two Poiseuille cells:
            // dp = R_system * Q_total and dp_fan = lambda*(6*s^2 - 60000*s*Q).
            let resistance = 32.0 * LENGTH * reference.state().dynamic_viscosity
                / (area() * DIAMETER * DIAMETER);
            let expected_q = density_ratio * speed * speed * 6.0
                / (resistance + density_ratio * speed * 60000.0);
            near(result.network.inlet_volume_flow_m3_per_s, expected_q, 3e-11);
            near(
                result.network.total_mass_flow_kg_per_s,
                inlet.density * expected_q,
                1e-10,
            );
            for branch in &result.network.branches {
                near(
                    branch.solution.mass_flow_kg_per_s,
                    inlet.density * expected_q / 2.0,
                    5e-11,
                );
            }
            near(
                result.network.driving_pressure_pa,
                resistance * result.network.inlet_volume_flow_m3_per_s,
                2e-7,
            );
            near(
                result.fan_pressure.value(),
                density_ratio
                    * (6.0 * speed * speed
                        - 60000.0 * speed * result.network.inlet_volume_flow_m3_per_s),
                1e-10,
            );
            near(result.pressure_residual.value(), 0.0, 1.1e-6);
            let binding = result.fan.gas_density_binding().unwrap();
            assert_eq!(binding.reference(), reference.state());
            assert_eq!(binding.inlet(), inlet);
            assert_eq!(result.fan.inlet_gas(), Some(inlet));
            assert_eq!(result.source.reference_bank, &retained);
            assert_eq!(result.fan.curve(), bank.curve());
            assert_eq!(result.fan.speed_ratio(), speed);
            assert_eq!(result.fan.count(), bank.count());
            assert_eq!(result.fan.arrangement(), bank.arrangement());
            assert_eq!(result.network.inlet.source_packs(), p.source.source_packs());
            assert_eq!(bank, retained);
            assert!(!bank.has_material_density());
            assert!(result.total_flow_evaluations > result.network.flow_evaluations);
            assert!(result.total_substeps >= result.network.substeps);
            near(result.network.total_wall_heat_w, 0.0, 2e-7);
            near(result.network.total_enthalpy_gain_w, 0.0, 2e-7);
            assert!(result.static_pressure_power_w > 0.0);
            near(
                result.static_pressure_power_w,
                result.fan_pressure.value() * result.network.inlet_volume_flow_m3_per_s,
                1e-15,
            );
            near(
                result.declared_shaft_power_upper_w,
                result.static_pressure_power_w * 1.05 / 0.3,
                1e-15,
            );
        }
    });
}

#[test]
fn fan_uses_inlet_volume_while_heating_changes_flow_split_and_enthalpy() {
    with_cx(|cx, _| {
        let (hot, cold) = heated();
        let reference = hot
            .source
            .evaluate(&source::point(400.0, hot.static_pressure_pa))
            .unwrap();
        let transport = laws();
        let paths = [
            branch("hot", &hot, &transport),
            branch("cold", &cold, &transport),
        ];
        let bank = bank(1.0, 0.0);
        let result = solve(
            cx,
            &paths,
            source_fan(&bank, &reference),
            work(),
            match_controls(),
        )
        .unwrap();
        let network = &result.network;
        assert!(
            network.branches[1].solution.mass_flow_kg_per_s
                > network.branches[0].solution.mass_flow_kg_per_s
        );
        assert_eq!(result.fan.inlet_gas(), Some(network.inlet.state()));
        assert_eq!(network.inlet.source_packs(), hot.source.source_packs());
        near(
            result
                .fan
                .pressure_at(VolumetricFlowRate::new(network.inlet_volume_flow_m3_per_s))
                .unwrap()
                .value(),
            result.fan_pressure.value(),
            1e-12,
        );
        let wrong_volume_pressure = result
            .fan
            .pressure_at(VolumetricFlowRate::new(network.outlet_volume_flow_m3_per_s))
            .unwrap()
            .value();
        assert!((wrong_volume_pressure - result.fan_pressure.value()).abs() > 0.01);
        near(result.pressure_residual.value(), 0.0, 1.1e-6);
        let mut original_enthalpy_gain = 0.0;
        for branch in &network.branches {
            let flow = &branch.solution;
            near(flow.pressure_residual_pa, 0.0, 1.1e-7);
            let t = flow.gas.cells.last().unwrap().gas.state().temperature;
            original_enthalpy_gain += flow.mass_flow_kg_per_s
                * (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.022)
                * (3.25 * (t - 400.0) + 0.005 * (t - 400.0) * (t + 400.0));
        }
        near(network.total_enthalpy_gain_w, original_enthalpy_gain, 1e-9);
        near(network.total_wall_heat_w, original_enthalpy_gain, 5e-7);
        assert!(result.total_flow_evaluations <= work().maximum_flow_evaluations);
        assert!(result.total_substeps <= work().maximum_substeps);
        assert!(
            result.declared_shaft_power_upper_w <= match_controls().maximum_neglected_shaft_power_w
        );
        // Declared fan work is independently reported, never silently deposited.
        assert!(result.declared_shaft_power_upper_w > 0.0);
        assert!(network.neglected_hydraulic_power_w > 0.0);
    });
}

#[test]
fn outer_work_and_authored_fan_energy_bounds_are_not_reset_between_trials() {
    with_cx(|cx, gate| {
        let (hot, cold) = heated();
        let saved = hot.clone();
        let reference = hot
            .source
            .evaluate(&source::point(400.0, hot.static_pressure_pa))
            .unwrap();
        let transport = laws();
        let paths = [
            branch("hot", &hot, &transport),
            branch("cold", &cold, &transport),
        ];
        let bank = bank(1.0, 0.0);
        let source = source_fan(&bank, &reference);
        let first = solve_parallel_pressure_driven_material_gas(
            cx,
            &paths,
            Pressure::new(2.0),
            controls(),
            feedback(),
            work(),
        )
        .unwrap();
        assert!(first.substeps > 0);
        for limited in [
            ParallelGasBudget {
                maximum_flow_evaluations: first.flow_evaluations,
                ..work()
            },
            ParallelGasBudget {
                maximum_substeps: first.substeps,
                ..work()
            },
        ] {
            assert!(matches!(
                solve(cx, &paths, source, limited, match_controls()),
                Err(GasFanError::WorkBudget {
                    network_evaluations: 1,
                    ..
                })
            ));
        }
        assert!(matches!(
            solve(
                cx,
                &paths,
                source,
                work(),
                GasFanControls {
                    maximum_network_evaluations: 2,
                    ..match_controls()
                }
            ),
            Err(GasFanError::NotConverged { evaluations: 2, .. })
        ));
        // Reject the original aggregate allowance before any remaining-work clipping.
        assert!(matches!(
            solve(
                cx,
                &paths,
                source,
                ParallelGasBudget {
                    maximum_flow_evaluations: 128 * MAX_PARALLEL_GAS_BRANCHES + 1,
                    ..work()
                },
                match_controls()
            ),
            Err(GasFanError::Network(ParallelGasError::Invalid(_)))
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                GasFanSource {
                    minimum_static_efficiency: 0.0,
                    ..source
                },
                work(),
                match_controls()
            ),
            Err(GasFanError::Invalid(_))
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                source,
                work(),
                GasFanControls {
                    maximum_neglected_shaft_power_w: 0.0,
                    ..match_controls()
                }
            ),
            Err(GasFanError::NeglectedShaftPower { .. })
        ));
        assert_eq!(hot, saved);
        gate.request();
        assert!(matches!(
            solve(cx, &paths, source, work(), match_controls()),
            Err(GasFanError::Network(ParallelGasError::Cancelled))
        ));
        assert_eq!(hot, saved);
    });
}

#[test]
fn unsupported_fan_source_and_stall_refuse_without_extrapolating() {
    with_cx(|cx, _| {
        let p = isothermal();
        let reference = p
            .source
            .evaluate(&source::point(400.0, p.static_pressure_pa))
            .unwrap();
        let transport = laws();
        let paths = [
            branch("left", &p, &transport),
            branch("right", &p, &transport),
        ];
        let stalled = bank(1.0, 8e-5);
        assert!(matches!(
            solve(
                cx,
                &paths,
                source_fan(&stalled, &reference),
                work(),
                match_controls()
            ),
            Err(GasFanError::Fan(AirflowError::StallRegion { .. }))
        ));
        let bank = bank(1.0, 0.0);
        let mut source = source_fan(&bank, &reference);
        source.gas_domain.density_ratio = [0.5, 0.75];
        assert!(matches!(
            solve(cx, &paths, source, work(), match_controls()),
            Err(GasFanError::Binding(MaterialFanError::OutsideDomain { .. }))
        ));
        source = source_fan(&bank, &reference);
        let rebound = FrozenFlowMedium::ThermochemicalGas(&reference)
            .bind_fan_bank(&bank, source.reference_medium, source.gas_domain)
            .unwrap();
        assert!(matches!(
            solve(
                cx,
                &paths,
                GasFanSource {
                    reference_bank: rebound.bank(),
                    ..source
                },
                work(),
                match_controls()
            ),
            Err(GasFanError::Binding(MaterialFanError::Invalid(_)))
        ));
        assert!(matches!(
            solve(
                cx,
                &paths,
                source,
                work(),
                GasFanControls {
                    pressure_bracket_pa: (4.0, 5.5),
                    ..match_controls()
                }
            ),
            Err(GasFanError::NotBracketed { .. })
        ));
    });
}

//! G1/G3: common source enthalpy, convection and hydraulic feedback.
#[path = "finite_gas_convection/source.rs"]
#[allow(dead_code)]
mod source;

#[path = "steady_material_gas/network.rs"]
mod network;

#[path = "steady_material_gas/gravity.rs"]
mod gravity;

use fs_airflow::graph::thermal::medium::convection::ForcedConvectionSettings;
use fs_airflow::graph::thermal::medium::series::hydraulic::MaterialFlowResistance;
use fs_airflow::graph::thermal::medium::steady_gas::*;
use fs_airflow::graph::thermal::medium::steady_gas::pressure::{
    GasPressureBudget, GasPressureError, solve_pressure_driven_material_gas,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::PressureDrivenBudget;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::CorrelationId;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::QueryPoint;
use fs_material::gas_field::PrimitiveGasState;
use fs_material::gas_field::steady::{AxialGasCell, GasInjection, SteadyAxialGas, SteadyGasSolve};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_qty::{Area, Length, MassFlowRate, Pressure};

const LENGTH: f64 = 1.0;
const DIAMETER: f64 = 0.01;
const FLOW: f64 = 2e-5;
fn area() -> f64 {
    core::f64::consts::PI * DIAMETER * DIAMETER / 4.0
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 32,
                    kernel_id: 735,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn primitive(t: f64, fraction: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k: t,
        composition: Composition::new(
            CompositionBasis::MoleFraction,
            vec![fraction, 1.0 - fraction],
        )
        .unwrap(),
    }
}
fn problem(variable: bool) -> SteadyAxialGas {
    SteadyAxialGas {
        source: source::model_transport(0.75, variable),
        source_conditions: QueryPoint::new(),
        static_pressure_pa: 100000.0,
        phase: ThermochemicalPhase::FrozenGas,
        inlet: primitive(400.0, 0.75),
        inlet_mass_flow_kg_s: FLOW,
        injections: vec![],
        cells: (0..2)
            .map(|index| AxialGasCell {
                end_m: (index + 1) as f64 * LENGTH,
                volume_m3: area() * LENGTH,
                wall_area_m2: core::f64::consts::PI * DIAMETER * LENGTH,
                wall_conductance_w_k: 0.0,
                wall_temperature_k: 300.0,
                initial_temperature_k: 400.0,
            })
            .collect(),
    }
}
fn laws() -> [GasCellTransport; 2] {
    [GasCellTransport {
        characteristic_length: Length::new(DIAMETER),
        correlation: CorrelationId::CircularDuctLaminarCwt,
        settings: ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(LENGTH / DIAMETER),
            ..ForcedConvectionSettings::default()
        },
        resistance: MaterialFlowResistance::try_new(
            32.0 * LENGTH / (area() * DIAMETER * DIAMETER),
            0.0,
            Area::new(area()),
            Length::new(DIAMETER),
            (0.0, 2300.0),
        )
        .unwrap(),
    }; 2]
}
fn controls() -> SteadyGasSolve {
    SteadyGasSolve {
        caloric: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 500.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 128,
        },
        maximum_dt_s: 0.05,
        maximum_relaxations: 2000,
        maximum_substeps: 10000,
        steady_power_tolerance_w: 1e-7,
        energy_tolerance_j: 1e-9,
        maximum_mean_mach: 0.05,
    }
}
fn feedback() -> GasTransportFeedback {
    GasTransportFeedback {
        maximum_iterations: 32,
        maximum_substeps: 100000,
        conductance_relaxation: 1.0,
        maximum_pressure_drop_fraction: 0.01,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn original_enthalpy_balance_and_same_material_hydraulics_close() {
    with_cx(|cx, _| {
        let p = problem(false);
        let solved = solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()).unwrap();
        let inlet = p
            .source
            .evaluate(&source::point(400.0, p.static_pressure_pa))
            .unwrap();
        let capacity = FLOW * inlet.state().specific_heat_cp;
        let mut previous = 400.0;
        let mut heat = 0.0;
        let mut pressure = 0.0;
        let mut power = 0.0;
        for cell in &solved.cells {
            let gas = cell.gas.state();
            let expected = (capacity * previous + cell.conductance_w_k * 300.0)
                / (capacity + cell.conductance_w_k);
            near(gas.temperature, expected, 2e-5);
            near(
                cell.wall_heat_w,
                capacity * (gas.temperature - previous),
                2e-7,
            );
            near(cell.velocity_m_s, FLOW / (gas.density * area()), 1e-12);
            let dp =
                32.0 * gas.dynamic_viscosity * LENGTH * cell.velocity_m_s / (DIAMETER * DIAMETER);
            near(cell.hydraulic.pressure_drop_pa, dp, 1e-12);
            assert_eq!(cell.hydraulic.transport, cell.transport);
            assert_eq!(cell.gas.source_packs(), p.source.source_packs());
            assert_eq!(
                cell.gas,
                solved
                    .thermal
                    .field()
                    .state_at(p.cells[cell.cell].end_m)
                    .unwrap()
            );
            previous = gas.temperature;
            heat += cell.wall_heat_w;
            pressure += dp;
            power += dp * FLOW / gas.density;
        }
        near(heat, capacity * (previous - 400.0), 3e-7);
        near(solved.total_pressure_drop_pa, pressure, 1e-12);
        near(solved.hydraulic_dissipation_w, power, 1e-12);
        assert!(solved.iterations >= 2);
        assert!(solved.substeps <= feedback().maximum_substeps);
        assert!(solved.maximum_power_residual_bound_w <= controls().steady_power_tolerance_w);
        assert!(solved.cells[1].velocity_m_s < solved.cells[0].velocity_m_s);
        assert_eq!(p, problem(false));
    });
}

#[test]
fn temperature_curves_and_transported_species_drive_live_heat_and_pressure() {
    with_cx(|cx, _| {
        let p = problem(true);
        let base = solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()).unwrap();
        assert!(
            base.iterations > 2,
            "temperature-dependent h must feed back"
        );
        assert!(base.cells[1].conductance_w_k < base.cells[0].conductance_w_k);
        let mut injected = p.clone();
        injected.injections.push(GasInjection {
            cell: 1,
            mass_flow_kg_s: FLOW,
            state: primitive(400.0, 0.0),
        });
        let mixed =
            solve_material_gas_transport(cx, &injected, &laws(), controls(), feedback()).unwrap();
        near(
            mixed.cells[0].gas.state().temperature,
            base.cells[0].gas.state().temperature,
            2e-5,
        );
        assert_eq!(mixed.thermal.mass_flow_kg_s(), &[FLOW, 2.0 * FLOW]);
        assert_ne!(
            mixed.cells[1].gas.thermodynamics(),
            base.cells[1].gas.thermodynamics()
        );
        assert_ne!(
            mixed.cells[1].coefficient.value,
            base.cells[1].coefficient.value
        );
        assert_ne!(
            mixed.cells[1].hydraulic.pressure_drop_pa,
            base.cells[1].hydraulic.pressure_drop_pa
        );
        assert_eq!(mixed.cells[1].gas.source_packs(), p.source.source_packs());
        assert!(mixed.maximum_power_residual_bound_w <= controls().steady_power_tolerance_w);
        // Common pressure changes velocity and resistance through the same EOS,
        // while the constant-mass source/correlation laws remain authoritative.
        let mut pressurized = p.clone();
        pressurized.static_pressure_pa *= 2.0;
        let high = solve_material_gas_transport(cx, &pressurized, &laws(), controls(), feedback())
            .unwrap();
        near(
            high.cells[1].gas.state().temperature,
            base.cells[1].gas.state().temperature,
            2e-5,
        );
        near(
            high.total_pressure_drop_pa,
            0.5 * base.total_pressure_drop_pa,
            1e-5,
        );
        near(
            high.cells[1].velocity_m_s,
            0.5 * base.cells[1].velocity_m_s,
            1e-6,
        );
    });
}

#[test]
fn hydraulic_admission_uses_the_converged_state_not_a_thermal_guess() {
    with_cx(|cx, _| {
        let p = problem(true);
        let base = solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()).unwrap();
        let uncooled = p.solve(controls(), || false).unwrap();
        let initial = p
            .cells
            .iter()
            .zip(laws())
            .map(|(cell, law)| {
                law.resistance
                    .evaluate(
                        &uncooled.field().state_at(cell.end_m).unwrap(),
                        MassFlowRate::new(FLOW),
                    )
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let initial_pressure = initial
            .iter()
            .map(|cell| cell.pressure_drop_pa)
            .sum::<f64>();
        let limit = (initial_pressure + base.total_pressure_drop_pa) / (2.0 * p.static_pressure_pa);
        assert!(base.total_pressure_drop_pa / p.static_pressure_pa < limit);
        assert!(initial_pressure / p.static_pressure_pa > limit);
        // Cooling reduces mu(T), raising the final Reynolds number. The
        // hydraulic card can legitimately admit that final state and exclude
        // the initial no-wall pseudo-time guess, which has no physical status.
        let final_min_re = base
            .cells
            .iter()
            .map(|cell| cell.hydraulic.reynolds)
            .fold(f64::INFINITY, f64::min);
        let min_re = 0.5 * (initial[0].reynolds + final_min_re);
        assert!(initial[0].reynolds < min_re && min_re < final_min_re);
        let mut admitted = laws();
        for law in &mut admitted {
            law.resistance = MaterialFlowResistance::try_new(
                32.0 * LENGTH / (area() * DIAMETER * DIAMETER),
                0.0,
                Area::new(area()),
                Length::new(DIAMETER),
                (min_re, 2300.0),
            )
            .unwrap();
        }
        let solved = solve_material_gas_transport(
            cx,
            &p,
            &admitted,
            controls(),
            GasTransportFeedback {
                maximum_pressure_drop_fraction: limit,
                ..feedback()
            },
        )
        .unwrap();
        near(
            solved.total_pressure_drop_pa,
            base.total_pressure_drop_pa,
            1e-12,
        );
        assert!(solved.maximum_power_residual_bound_w <= controls().steady_power_tolerance_w);
    });
}

#[test]
fn duplicate_owners_geometry_work_pressure_and_cancellation_refuse_transactionally() {
    with_cx(|cx, gate| {
        let p = problem(true);
        let saved = p.clone();
        let mut duplicate = p.clone();
        duplicate.cells[0].wall_conductance_w_k = 1.0;
        assert!(matches!(
            solve_material_gas_transport(cx, &duplicate, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Invalid(_))
        ));
        let mut wrong_geometry = p.clone();
        wrong_geometry.cells[0].volume_m3 *= 2.0;
        assert!(matches!(
            solve_material_gas_transport(cx, &wrong_geometry, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Invalid(_))
        ));
        assert!(matches!(
            solve_material_gas_transport(
                cx,
                &p,
                &laws(),
                controls(),
                GasTransportFeedback {
                    maximum_iterations: 1,
                    ..feedback()
                }
            ),
            Err(MaterialSteadyGasError::NotConverged { .. })
        ));
        assert!(matches!(
            solve_material_gas_transport(
                cx,
                &p,
                &laws(),
                controls(),
                GasTransportFeedback {
                    maximum_pressure_drop_fraction: 1e-12,
                    ..feedback()
                }
            ),
            Err(MaterialSteadyGasError::PressureDrop { .. })
        ));
        assert_eq!(p, saved);
        gate.request();
        assert!(matches!(
            solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Cancelled)
        ));
        assert_eq!(p, saved);
    });
}

fn pressure_budget() -> GasPressureBudget {
    GasPressureBudget {
        pressure: PressureDrivenBudget {
            mass_flow_bracket_kg_per_s: (0.2 * FLOW, 2.75 * FLOW),
            pressure_absolute_tolerance_pa: 1e-7,
            pressure_relative_tolerance: 1e-8,
            max_flow_evaluations: 48,
            maximum_neglected_hydraulic_power_w: 1.0,
        },
        momentum_flux_factor: 4.0 / 3.0,
        kinetic_energy_flux_factor: 2.0,
        maximum_substeps: 1_000_000,
        maximum_neglected_kinetic_power_w: 1.0,
    }
}

#[test]
fn pressure_selected_isothermal_flow_uses_eos_density_and_the_supplied_drive() {
    with_cx(|cx, _| {
        let mut p = problem(false);
        for cell in &mut p.cells {
            cell.wall_temperature_k = 400.0;
        }
        let reference = p
            .source
            .evaluate(&source::point(400.0, p.static_pressure_pa))
            .unwrap();
        // Independent two-cell Poiseuille oracle; zero thermal acceleration.
        let drive = 64.0 * LENGTH * reference.state().dynamic_viscosity * FLOW
            / (reference.state().density * area() * DIAMETER * DIAMETER);
        p.inlet_mass_flow_kg_s = 0.3 * FLOW; // A template value, not the answer.
        for multiplier in [1.0, 2.0] {
            p.static_pressure_pa = 100000.0 * multiplier;
            let saved = p.clone();
            let solved = solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                GasPressureBudget {
                    momentum_flux_factor: 1.0,
                    kinetic_energy_flux_factor: 1.0,
                    ..pressure_budget()
                },
            )
            .unwrap();
            near(
                solved.mass_flow_kg_per_s,
                multiplier * FLOW,
                1e-6 * multiplier * FLOW,
            );
            near(solved.acceleration_pressure_pa, 0.0, 1e-10);
            assert!(
                solved.pressure_residual_pa.abs()
                    <= pressure_budget().pressure.pressure_absolute_tolerance_pa
            );
            for cell in &solved.gas.cells {
                near(
                    cell.gas.state().density,
                    multiplier * reference.state().density,
                    1e-11,
                );
                near(
                    cell.velocity_m_s,
                    FLOW / (reference.state().density * area()),
                    1e-6,
                );
                assert_eq!(cell.gas.source_packs(), p.source.source_packs());
            }
            assert!(solved.flow_evaluations <= pressure_budget().pressure.max_flow_evaluations);
            assert!(solved.substeps <= pressure_budget().maximum_substeps);
            assert_eq!(p, saved);
        }
    });
}

#[test]
fn pressure_selected_heating_and_cooling_close_density_momentum_and_enthalpy() {
    with_cx(|cx, _| {
        let p = problem(true);
        let given = solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()).unwrap();
        let inlet = given.thermal.field().state_at(0.0).unwrap();
        let inlet_speed = FLOW / (inlet.state().density * area());
        let outlet_speed = FLOW / (given.cells.last().unwrap().gas.state().density * area());
        let drive = given.total_pressure_drop_pa
            + pressure_budget().momentum_flux_factor
                * (FLOW / area())
                * (outlet_speed - inlet_speed);
        let mut previous_mass = 0.0;
        for wall in [300.0, 480.0] {
            let mut trial = p.clone();
            for cell in &mut trial.cells {
                cell.wall_temperature_k = wall;
            }
            let solved = solve_pressure_driven_material_gas(
                cx,
                &trial,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                pressure_budget(),
            )
            .unwrap();
            let mass = solved.mass_flow_kg_per_s;
            let inlet = solved.gas.thermal.field().state_at(0.0).unwrap();
            let cp = inlet.state().specific_heat_cp;
            let mut previous_temperature = inlet.state().temperature;
            let mut previous_speed = mass / (inlet.state().density * area());
            let initial_speed = previous_speed;
            let mut pressure_loss = 0.0;
            let mut kinetic = 0.0;
            let mut kinetic_bound = 0.0;
            for cell in &solved.gas.cells {
                let gas = cell.gas.state();
                let speed = mass / (gas.density * area());
                near(gas.density * area() * cell.velocity_m_s, mass, 1e-14);
                near(
                    cell.wall_heat_w,
                    mass * cp * (gas.temperature - previous_temperature),
                    3e-7,
                );
                pressure_loss +=
                    32.0 * gas.dynamic_viscosity * LENGTH * speed / (DIAMETER * DIAMETER);
                let change = pressure_budget().kinetic_energy_flux_factor
                    * 0.5
                    * mass
                    * (speed * speed - previous_speed * previous_speed);
                kinetic += change;
                kinetic_bound += change.abs();
                previous_speed = speed;
                previous_temperature = gas.temperature;
                assert_eq!(
                    cell.gas,
                    solved
                        .gas
                        .thermal
                        .field()
                        .state_at(trial.cells[cell.cell].end_m)
                        .unwrap()
                );
            }
            let acceleration = pressure_budget().momentum_flux_factor
                * (mass / area())
                * (previous_speed - initial_speed);
            assert_eq!(
                solved.momentum_flux_factor,
                pressure_budget().momentum_flux_factor
            );
            assert_eq!(
                solved.kinetic_energy_flux_factor,
                pressure_budget().kinetic_energy_flux_factor
            );
            near(solved.acceleration_pressure_pa, acceleration, 1e-12);
            near(pressure_loss + acceleration, drive, 2e-7);
            near(solved.kinetic_power_change_w, kinetic, 1e-15);
            near(solved.neglected_kinetic_power_bound_w, kinetic_bound, 1e-15);
            assert!(
                solved.substeps > solved.gas.substeps,
                "all pressure trials share one work account"
            );
            if wall == 300.0 {
                near(mass, FLOW, 1e-6 * FLOW);
                assert!(acceleration < 0.0);
                previous_mass = mass;
            } else {
                assert!(acceleration > 0.0);
                assert!(
                    mass < previous_mass,
                    "heating reduces density and increases required pressure"
                );
            }
        }
    });
}

#[test]
fn pressure_gas_shared_work_energy_and_pressure_bounds_refuse_without_mutation() {
    with_cx(|cx, gate| {
        let p = problem(true);
        let saved = p.clone();
        let base = solve_material_gas_transport(cx, &p, &laws(), controls(), feedback()).unwrap();
        let inlet = base.thermal.field().state_at(0.0).unwrap();
        let acceleration = pressure_budget().momentum_flux_factor
            * (FLOW / area())
            * (FLOW / (base.cells[1].gas.state().density * area())
                - FLOW / (inlet.state().density * area()));
        let drive = base.total_pressure_drop_pa + acceleration;
        let mut budget = pressure_budget();
        budget.pressure.mass_flow_bracket_kg_per_s = (FLOW, 2.0 * FLOW);
        // The known solution is the first endpoint, making each refusal local.
        let mut no_hydraulic_heat = budget;
        no_hydraulic_heat
            .pressure
            .maximum_neglected_hydraulic_power_w = 0.0;
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                no_hydraulic_heat
            ),
            Err(GasPressureError::NeglectedHydraulicPower { .. })
        ));
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                GasPressureBudget {
                    maximum_neglected_kinetic_power_w: 0.0,
                    ..budget
                }
            ),
            Err(GasPressureError::NeglectedKineticPower { .. })
        ));
        // Friction alone fits, but the same source density's momentum term does not.
        let fraction =
            (base.total_pressure_drop_pa + 0.5 * acceleration.abs()) / p.static_pressure_pa;
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                GasTransportFeedback {
                    maximum_pressure_drop_fraction: fraction,
                    ..feedback()
                },
                budget
            ),
            Err(GasPressureError::PressureVariation { .. })
        ));
        assert!(base.substeps > 0);
        // The first trial has its original work plus one spare substep;
        // the second full thermal solve cannot reset the shared allowance.
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(1.1 * drive),
                controls(),
                feedback(),
                GasPressureBudget {
                    maximum_substeps: base.substeps + 1,
                    ..budget
                }
            ),
            Err(GasPressureError::SubstepBudget { .. })
                | Err(GasPressureError::Transport(
                    MaterialSteadyGasError::NotConverged { .. }
                ))
                | Err(GasPressureError::Transport(
                    MaterialSteadyGasError::Thermal(_)
                ))
        ));
        assert_eq!(p, saved);
        gate.request();
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                budget
            ),
            Err(GasPressureError::Transport(
                MaterialSteadyGasError::Cancelled
            ))
        ));
        assert_eq!(p, saved);
    });
}

#[test]
fn pressure_gas_requires_an_admitted_constant_area_sign_bracket() {
    with_cx(|cx, _| {
        let mut p = problem(false);
        for cell in &mut p.cells {
            cell.wall_temperature_k = 400.0;
        }
        let material = p
            .source
            .evaluate(&source::point(400.0, p.static_pressure_pa))
            .unwrap();
        let drive = 64.0 * LENGTH * material.state().dynamic_viscosity * FLOW
            / (material.state().density * area() * DIAMETER * DIAMETER);
        let mut budget = pressure_budget();
        budget.pressure.mass_flow_bracket_kg_per_s = (2.0 * FLOW, 2.75 * FLOW);
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                budget
            ),
            Err(GasPressureError::NotBracketed { .. })
        ));
        budget = pressure_budget();
        budget.pressure.max_flow_evaluations = 2;
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                budget
            ),
            Err(GasPressureError::NotConverged { evaluations: 2, .. })
        ));
        for (beta, alpha) in [
            (0.5, 1.0),
            (2.0, 2.0),
            (f64::NAN, 2.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(matches!(
                solve_pressure_driven_material_gas(
                    cx,
                    &p,
                    &laws(),
                    Pressure::new(drive),
                    controls(),
                    feedback(),
                    GasPressureBudget {
                        momentum_flux_factor: beta,
                        kinetic_energy_flux_factor: alpha,
                        ..pressure_budget()
                    }
                ),
                Err(GasPressureError::Invalid(_))
            ));
        }
        let mut varying_area = laws();
        varying_area[1].resistance = MaterialFlowResistance::try_new(
            1.0,
            0.0,
            Area::new(2.0 * area()),
            Length::new(DIAMETER),
            (0.0, 2300.0),
        )
        .unwrap();
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &varying_area,
                Pressure::new(drive),
                controls(),
                feedback(),
                pressure_budget()
            ),
            Err(GasPressureError::Invalid(_))
        ));
        p.injections.push(GasInjection {
            cell: 1,
            mass_flow_kg_s: FLOW,
            state: primitive(400.0, 0.75),
        });
        assert!(matches!(
            solve_pressure_driven_material_gas(
                cx,
                &p,
                &laws(),
                Pressure::new(drive),
                controls(),
                feedback(),
                pressure_budget()
            ),
            Err(GasPressureError::Invalid(_))
        ));
    });
}

//! The actual thermal-only nonlinear owner, starting without a solved endpoint.
use super::*;

#[path = "implicit/adaptive.rs"]
mod adaptive;
use fs_couple::circuit_enthalpy::contacts::heat::{
    PreparedHeatNetwork,
    thermal::{HeatAssemblySolve, step_heat_assembly_candidate},
};
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_solver::NewtonKrylovConfig;
use std::cell::Cell;

fn controls() -> HeatAssemblySolve {
    HeatAssemblySolve {
        maximum_cells: 2,
        maximum_evaluations: 600,
        maximum_newton_iterations: 40,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7,
        closed_energy_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn with_fixture(
    scale: f64,
    reversed: bool,
    run: impl FnOnce(
        &fs_exec::Cx<'_>,
        &fs_exec::CancelGate,
        &PreparedHeatNetwork<'_, '_>,
        &PreparedSourceConjugateLinks<'_, '_>,
        &LiquidTemperatureLaw<'_>,
        &DepositionSimulation<'_>,
    ),
) {
    with_scheduled_fixture(scale, reversed, &[], run);
}

fn with_scheduled_fixture(
    scale: f64,
    reversed: bool,
    scheduled: &[TimedEvent],
    run: impl FnOnce(
        &fs_exec::Cx<'_>,
        &fs_exec::CancelGate,
        &PreparedHeatNetwork<'_, '_>,
        &PreparedSourceConjugateLinks<'_, '_>,
        &LiquidTemperatureLaw<'_>,
        &DepositionSimulation<'_>,
    ),
) {
    let (solid_card, solid_caloric) = wall_source(scale);
    let solid_point = wall_point(310.0);
    let layer = wall_layer(
        &solid_card,
        solid_caloric.curve(),
        &solid_point,
        [280.0, 320.0],
    );
    let (fluid_card, fluid_caloric) = fluid_source(1.0);
    let liquid = liquid_law(&fluid_card);
    let fluid_point = context(290.0);
    let thermal = pair(solid_caloric.curve(), &fluid_caloric, reversed, 0.0);
    let thermal = if scheduled.is_empty() {
        thermal
    } else {
        let mut events = thermal.cells().values().cloned().map(|cell| TimedEvent {
            time_s: 0.0, event: Event::Birth(cell),
        }).collect::<Vec<_>>();
        events.extend(thermal.links().values().cloned().map(|link| TimedEvent {
            time_s: 0.0, event: Event::Contact(link),
        }));
        events.extend_from_slice(scheduled);
        DepositionSimulation::new(events, &[
            (MaterialId(0), solid_caloric.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid_caloric.curve() as &dyn ThermalLaw),
        ], DepositionBudget {
            max_cells: 2, max_links: 1, max_events: 3 + scheduled.len(), max_substeps: 200,
        }).unwrap().advance_to(0.0, StepConfig {
            max_dt_s: 1.0, energy_tolerance_j: 1e-10,
        }, || Ok(())).unwrap().0
    };
    let owner =
        PreparedSourceConjugateLinks::new(&thermal, &[binding(&layer)], 1, || Ok(())).unwrap();
    let initial = owner
        .resolve(
            &thermal,
            config(),
            |_, _, t| liquid.at_temperature(t),
            || Ok(()),
        )
        .unwrap();
    let cells = [
        MaterialCell {
            cell: IDS[0],
            material_id: MaterialId(0),
            curve: solid_caloric.curve(),
            point: ThermalMaterialPoint {
                material: &solid_card,
                context: &solid_point,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        },
        MaterialCell {
            cell: IDS[1],
            material_id: MaterialId(1),
            curve: fluid_caloric.curve(),
            point: ThermalMaterialPoint {
                material: &fluid_card,
                context: &fluid_point,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        },
    ];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 173,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
            .unwrap();
        run(&cx, &gate, &prepared, &owner, &liquid, &thermal);
    });
}

#[test]
fn cold_start_solves_both_finite_stores_and_the_actual_surface() {
    let dt = 0.2;
    let expected = endpoint(1.0, dt);
    for reversed in [false, true] {
        with_fixture(1.0, reversed, |cx, _, prepared, owner, liquid, thermal| {
            let saved = thermal.cells().clone();
            let assembly = SourceConjugateHeatAssembly {
                heat: prepared,
                links: owner,
                layer_config: config(),
                resolve_fluid: |_: &SourceConvectiveLink,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    liquid.at_temperature(t)
                },
            };
            let (next, report) =
                step_heat_assembly_candidate(&assembly, cx, thermal, dt, controls(), &[0.0; 2])
                    .unwrap();
            assert!(report.nonlinear.converged && report.evaluations > 1);
            for (cell, temperature) in IDS.into_iter().zip(expected) {
                near(next.temperature_k(cell).unwrap(), temperature, 2e-7);
            }
            let sample = &report.source.conjugate.samples()[0];
            let (surface, power) = surface_oracle(expected[0], expected[1], 1.0);
            near(sample.layer.surface_temperature_k, surface, 2e-7);
            near(sample.applied_heat_into_fluid_w, power, 2e-7);
            near(
                sample
                    .layer
                    .film
                    .surface
                    .as_ref()
                    .unwrap()
                    .thermal
                    .temperature_k(),
                surface,
                2e-7,
            );
            assert_eq!(
                sample.wall_temperature_k,
                next.temperature_k(IDS[0]).unwrap()
            );
            assert_eq!(
                sample.fluid_temperature_k,
                next.temperature_k(IDS[1]).unwrap()
            );
            near(
                MASS * CP * (next.temperature_k(IDS[0]).unwrap() - 310.0)
                    + dt * sample.applied_heat_into_fluid_w,
                0.0,
                1e-7,
            );
            near(
                MASS * FLUID_CP * (next.temperature_k(IDS[1]).unwrap() - 290.0)
                    - dt * sample.applied_heat_into_fluid_w,
                0.0,
                1e-7,
            );
            near(next.energy_j(), thermal.energy_j(), 1e-7);
            near(report.external_heat_j, 0.0, 1e-12);
            near(report.energy_residual_j, 0.0, 1e-7);
            assert!(report.maximum_cell_residual_j <= controls().cell_energy_tolerance_j);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.time_s(), 0.0);
        });
    }
}

#[test]
fn stiff_endpoint_responds_to_source_conductivity_and_cannot_use_frozen_initial_heat() {
    let dt = 2.0;
    let mut accepted_walls = Vec::new();
    for scale in [1.0, 4.0] {
        with_fixture(scale, false, |cx, _, prepared, owner, liquid, thermal| {
            let assembly = SourceConjugateHeatAssembly {
                heat: prepared,
                links: owner,
                layer_config: config(),
                resolve_fluid: |_: &SourceConvectiveLink,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    liquid.at_temperature(t)
                },
            };
            let (next, report) =
                step_heat_assembly_candidate(&assembly, cx, thermal, dt, controls(), &[0.0; 2])
                    .unwrap();
            let expected = endpoint(scale, dt);
            let initial_g = surface_oracle(310.0, 290.0, scale).1 / 20.0;
            let stiffness = dt * initial_g * (1.0 / (MASS * CP) + 1.0 / (MASS * FLUID_CP));
            assert!(
                stiffness > 8.0,
                "this interval exceeds the explicit pair rate scale"
            );
            let frozen_wall = 294.0 + 16.0 / (1.0 + stiffness);
            let wall = next.temperature_k(IDS[0]).unwrap();
            near(wall, expected[0], 2e-7);
            near(next.temperature_k(IDS[1]).unwrap(), expected[1], 2e-7);
            assert!(
                (wall - frozen_wall).abs() > 0.05,
                "initial source/surface coefficients must not be frozen across the solve"
            );
            near(report.energy_residual_j, 0.0, 1e-7);
            accepted_walls.push(wall);
        });
    }
    assert!(
        accepted_walls[1] < accepted_walls[0] - 0.5,
        "changing the original conductivity source must change actual stiff cooling"
    );
}

#[test]
fn nonlinear_work_source_and_cancellation_refusals_publish_no_state() {
    with_fixture(1.0, false, |cx, gate, prepared, owner, liquid, thermal| {
        let saved = thermal.cells().clone();
        let assembly = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |_: &SourceConvectiveLink,
                            _: &DepositionSimulation<'_>,
                            t: Temperature| liquid.at_temperature(t),
        };
        let mut exhausted = controls();
        exhausted.maximum_evaluations = 1;
        assert!(matches!(
            step_heat_assembly_candidate(&assembly, cx, thermal, 0.2, exhausted, &[0.0; 2]),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        assert!(
            step_heat_assembly_candidate(&assembly, cx, thermal, 0.0, controls(), &[0.0; 2],)
                .is_err()
        );
        let later_trial_seen = Cell::new(false);
        let limited = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |binding: &SourceConvectiveLink,
                            state: &DepositionSimulation<'_>,
                            t: Temperature| {
                if state.temperature_k(binding.fluid_cell).unwrap() > 290.5 {
                    later_trial_seen.set(true);
                    return Err("authored finite-fluid source support ended");
                }
                liquid
                    .at_temperature(t)
                    .map_err(|_| "original liquid source refused")
            },
        };
        assert!(matches!(
            step_heat_assembly_candidate(&limited, cx, thermal, 0.2, controls(), &[0.0; 2]),
            Err(CircuitEnthalpyError::Source(
                SourceConjugateHeatError::Conjugate(SourceConjugateError::Convection(
                    SourceCellConvectionError::Source {
                        source: "authored finite-fluid source support ended",
                        ..
                    }
                ))
            ))
        ));
        assert!(later_trial_seen.get());
        // A refusal must not poison the same original material/heat owners.
        let (accepted, _) =
            step_heat_assembly_candidate(&assembly, cx, thermal, 0.2, controls(), &[0.0; 2])
                .unwrap();
        near(
            accepted.temperature_k(IDS[0]).unwrap(),
            endpoint(1.0, 0.2)[0],
            2e-7,
        );
        let queries = Cell::new(0usize);
        let interrupted = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |_: &SourceConvectiveLink,
                            _: &DepositionSimulation<'_>,
                            t: Temperature| {
                queries.set(queries.get() + 1);
                if queries.get() == 3 {
                    gate.request();
                }
                liquid.at_temperature(t)
            },
        };
        assert!(matches!(
            step_heat_assembly_candidate(&interrupted, cx, thermal, 0.2, controls(), &[0.0; 2]),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert!(
            queries.get() >= 3,
            "cancel while resolving an actual layer/film trial"
        );
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

//! Manufactured heterogeneous segments on the real material/thermal owners.

use super::{ElectricalPropertyBasis, OhmicDrive, source_fixture};
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::constitutive_graph::MaterialGraphSolve;
use fs_couple::constitutive_graph::enthalpy::GraphThermalSolve;
use fs_couple::incandescent::filament::{FilamentError, FilamentSegment, PreparedSeriesFilament};
use fs_exec::{CancelGate, Cx};
use fs_matdb::PropertyKey;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const CELLS: [CellId; 2] = [CellId(17), CellId(4)];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 81,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn timeline<'a>(
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    h: [f64; 2],
    length: f64,
) -> DepositionSimulation<'a> {
    let mut events = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELLS[i],
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(i as u64),
                mass_kg: 0.1,
                specific_h_j_kg: h[i],
                geometry: CellGeometry {
                    start_m: [i as f64 * length, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * length, 0.0, 0.0],
                    volume_m3: 0.0001,
                    surface_area_m2: 0.03,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 8,
            a: CELLS[0],
            b: CELLS[1],
            area_m2: 0.01,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn controls() -> GraphThermalSolve {
    GraphThermalSolve {
        maximum_cells: 2,
        maximum_thermal_trials: 400,
        maximum_graph_evaluations: 4000,
        maximum_iterations: 24,
        difference_step_j: 1e-4,
        cell_tolerance_j: 1e-7,
        total_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..Default::default()
        },
        graph: MaterialGraphSolve {
            maximum_unknowns: 2,
            maximum_terms: 8,
            maximum_evaluations: 400,
            maximum_iterations: 24,
            equation_tolerance: 1e-10,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: 4,
                ..Default::default()
            },
        },
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-7 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn exercise(drive: OhmicDrive, phase: bool) {
    with_cx(|cx, _| {
        let (a, ca) = source_fixture(ElectricalPropertyBasis::Resistivity, phase);
        let (b, cb) = source_fixture(ElectricalPropertyBasis::Conductivity, phase);
        let nodes = [a, b];
        let curves = [ca, cb];
        let initial_h = if phase {
            [3500.0, 4500.0]
        } else {
            [500.0, 3500.0]
        };
        let initial = timeline(&curves, initial_h, 0.01);
        let prepared = PreparedSeriesFilament::try_new(
            cx,
            &initial,
            (0..2)
                .map(|i| FilamentSegment {
                    cell: CELLS[i],
                    curve: &curves[i],
                    source: nodes[i].clone(),
                })
                .collect(),
            drive,
            2,
        )
        .unwrap();
        let contexts = nodes
            .each_ref()
            .map(|node| node.source().original().material().query());
        let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = [0, 1].map(|i| MaterialCell {
            cell: CELLS[i],
            material_id: MaterialId(i as u64),
            curve: &curves[i],
            point: ThermalMaterialPoint {
                material: nodes[i].source().card(),
                context: &contexts[i],
                temperature_axis: "T",
                liquid_fraction_axis: if phase { Some("f") } else { None },
            },
            conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 8,
            distance_a_m: 0.005,
            distance_b_m: 0.005,
        }];
        let bulk = MaterialNetwork {
            cells: &bindings,
            resistors: &[],
            links: &links,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let network = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let heat = network.prepare(cx, &initial).unwrap();
        let dt = if phase { 0.01 } else { 1e-7 };
        let result = heat
            .graph_step_candidate(
                cx,
                &initial,
                dt,
                prepared.thermal_input(cx, &initial).unwrap(),
                controls(),
                &[0.0; 2],
            )
            .unwrap();
        let temperatures = CELLS.map(|id| result.thermal().temperature_k(id).unwrap());
        let resistances = [0, 1].map(|i| {
            if phase {
                let f = curves[i]
                    .values_at_specific_enthalpy(
                        result.thermal().cells()[&CELLS[i]].specific_h_j_kg,
                    )
                    .unwrap()
                    .liquid_mass_fraction;
                1.0 + 2.0 * f
            } else {
                temperatures[i] / 300.0
            }
        });
        // Independent authored fixture values: rho_a=q, sigma_b=q; L/A=1.
        let resistances = [resistances[0], 1.0 / resistances[1]];
        let current = match drive {
            OhmicDrive::Current(i) => i,
            OhmicDrive::Voltage(v) => v / resistances.iter().sum::<f64>(),
        };
        let total_power = current * current * resistances.iter().sum::<f64>();
        close(result.report().graph_heat_j, dt * total_power);
        close(
            result.thermal().energy_j() - initial.energy_j(),
            dt * total_power,
        );
        close(result.report().other_external_heat_j, 0.0);
        let conductance = result.report().source.material.links[0].conductance_w_k;
        close(conductance, 1e6);
        for i in 0..2 {
            let node = prepared.node_id(CELLS[i]).unwrap().to_owned();
            close(
                result.graph().external[&(node.clone(), "current".into())],
                current,
            );
            close(
                result.graph().graph.outputs[&(node.clone(), "voltage".into())],
                current * resistances[i],
            );
            let power = current * current * resistances[i];
            close(
                result.graph().graph.outputs[&(node.clone(), "heat".into())],
                power,
            );
            let stored = 0.1 * (result.thermal().cells()[&CELLS[i]].specific_h_j_kg - initial_h[i]);
            close(
                stored,
                dt * (power + conductance * (temperatures[1 - i] - temperatures[i])),
            );
            assert_eq!(
                result.thermal().cells()[&CELLS[i]].geometry,
                initial.cells()[&CELLS[i]].geometry
            );
            assert_eq!(result.thermal().cells()[&CELLS[i]].mass_kg, 0.1);
            if phase {
                close(temperatures[i], 600.0);
                let fraction = curves[i]
                    .values_at_specific_enthalpy(
                        result.thermal().cells()[&CELLS[i]].specific_h_j_kg,
                    )
                    .unwrap()
                    .liquid_mass_fraction;
                assert!(fraction > if i == 0 { 0.25 } else { 0.75 });
                close(
                    result.graph().external[&(node, "liquid_mass_fraction_change".into())],
                    fraction,
                );
            }
        }
        if !phase {
            assert!(temperatures[0] > 350.0 && temperatures[1] < 650.0);
        }
        assert!(prepared.thermal_input(cx, result.thermal()).is_ok());
        assert_eq!(initial.cells()[&CELLS[0]].specific_h_j_kg, initial_h[0]);
    });
}

#[test]
fn signed_drives_split_native_heat_and_preserve_axial_exchange() {
    for drive in [
        OhmicDrive::Current(200.0),
        OhmicDrive::Current(-200.0),
        OhmicDrive::Voltage(400.0),
        OhmicDrive::Voltage(-400.0),
    ] {
        exercise(drive, false);
    }
}

#[test]
fn different_phase_cells_update_native_resistance_at_constant_temperature() {
    for drive in [OhmicDrive::Current(-2.0), OhmicDrive::Voltage(6.0)] {
        exercise(drive, true);
    }
}

#[test]
fn original_cell_chart_geometry_work_cap_and_cancellation_are_enforced() {
    with_cx(|cx, gate| {
        let (node, curve) = source_fixture(ElectricalPropertyBasis::Resistivity, false);
        let (_, other) = source_fixture(ElectricalPropertyBasis::Conductivity, false);
        let curves = [curve, other];
        let thermal = timeline(&curves, [500.0, 500.0], 0.01);
        let segment = || FilamentSegment {
            cell: CELLS[0],
            curve: &curves[0],
            source: node.clone(),
        };
        let prepare = |segments, cap| {
            PreparedSeriesFilament::try_new(cx, &thermal, segments, OhmicDrive::Current(1.0), cap)
        };
        assert!(matches!(
            prepare(vec![segment(), segment()], 2),
            Err(FilamentError::Invalid(_))
        ));
        assert!(matches!(
            prepare(vec![segment()], 0),
            Err(FilamentError::Invalid(_))
        ));
        let copy = curves[0].clone();
        assert!(matches!(
            prepare(
                vec![FilamentSegment {
                    curve: &copy,
                    ..segment()
                }],
                1
            ),
            Err(FilamentError::Invalid(_))
        ));
        assert!(matches!(
            prepare(
                vec![FilamentSegment {
                    cell: CELLS[1],
                    curve: &curves[1],
                    source: node.clone()
                }],
                1
            ),
            Err(FilamentError::Electrical(_))
        ));
        let changed_geometry = timeline(&curves, [500.0, 500.0], 0.02);
        assert!(matches!(
            PreparedSeriesFilament::try_new(
                cx,
                &changed_geometry,
                vec![segment()],
                OhmicDrive::Current(1.0),
                1
            ),
            Err(FilamentError::Invalid(_))
        ));
        let prepared = prepare(vec![segment()], 1).unwrap();
        assert!(matches!(
            prepared.thermal_input(cx, &changed_geometry),
            Err(FilamentError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            prepare(vec![segment()], 1),
            Err(FilamentError::Cancelled)
        ));
        assert!(matches!(
            prepared.thermal_input(cx, &thermal),
            Err(FilamentError::Cancelled)
        ));
    });
}

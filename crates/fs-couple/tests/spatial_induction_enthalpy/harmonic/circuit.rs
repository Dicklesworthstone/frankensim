//! Loaded induction and an original material resistor heat two finite bodies.
use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell,
    material::CellResistor,
    spatial_induction::harmonic::circuit::{
        SpatialHarmonicCircuitThermalCandidate, SpatialHarmonicCircuitThermalError,
        heat_assembly_candidate as circuit_heat_candidate,
    },
};
use fs_em::spatial::vector_potential::induction::{
    harmonic::terminals::{HarmonicCircuitControls, SpatialHarmonicTerminalControls},
    terminals::{SpatialTerminalInduction, current::SpatialCurrentInduction},
};
use fs_phs::{
    circuit::{Branch, CircuitGraph},
    harmonic::HarmonicBudget,
};

const LOAD_L: f64 = 0.004;
const LOAD_C: f64 = 0.025;

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    near(actual.re, expected.re, tolerance);
    near(actual.im, expected.im, tolerance);
}

fn excitation() -> [[C64; 3]; 6] {
    [[C64::new(6.4, -4.8); 3]; 6]
}

fn template(problem: &SpatialVectorMagnetostatics) -> SpatialCurrentInduction<'_> {
    let faces: Vec<Vec<_>> = [0.0, 1.0]
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
    let limits = SpatialInductionLimits {
        maximum_unknowns: 1000,
        maximum_conducting_cells: 100,
        maximum_assembly_terms: 10000,
    };
    SpatialTerminalInduction::new(
        initial(problem, BASIS, 0.7, 1.0),
        &faces,
        limits,
        &mut || false,
    )
    .unwrap()
    .into_current_driven(limits, &mut || false)
    .unwrap()
}

fn graph(driven: bool) -> CircuitGraph {
    // Source-owned R deliberately differs from either material temperature;
    // it must be lowered again at the accepted state and every proposal.
    let mut branches = if driven {
        vec![
            (4, 3, Branch::Resistor { ohms: 999.0 }),
            (3, 2, Branch::Inductor { henries: LOAD_L }),
            (2, 1, Branch::Capacitor { farads: LOAD_C }),
        ]
    } else {
        vec![
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 3, Branch::Inductor { henries: LOAD_L }),
            (3, 0, Branch::Capacitor { farads: LOAD_C }),
        ]
    };
    branches.push((1, 0, Branch::Inductor { henries: f64::NAN }));
    if driven {
        branches.push((4, 0, Branch::VoltageSource { port: 0 }));
    }
    CircuitGraph {
        node_count: if driven { 5 } else { 4 },
        branches,
        transformers: vec![],
    }
}

fn with_resistor<R>(
    heat: &MaterialHeatNetwork<'_>,
    destination: CellId,
    run: impl FnOnce(&MaterialHeatNetwork<'_>) -> R,
) -> R {
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 0,
            cell: destination,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let bulk = MaterialNetwork {
        cells: heat.electrical.bulk.cells,
        resistors: &resistors,
        links: heat.electrical.bulk.links,
    };
    let bases = [(0, BASIS)];
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: heat.electrical.contacts,
        electrical_bases: &bases,
    };
    let loaded = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: heat.thermal_contacts,
        radiation: heat.radiation,
    };
    run(&loaded)
}

fn terminal_controls() -> SpatialHarmonicTerminalControls {
    SpatialHarmonicTerminalControls {
        field: ac_controls(),
        maximum_ports: 4,
        maximum_response_scalars: 100000,
        current_probe_a: 0.25,
    }
}

fn circuit_controls() -> HarmonicCircuitControls {
    HarmonicCircuitControls {
        native: HarmonicBudget {
            max_dimension: 128,
            max_dense_bytes: 4_000_000,
            max_dense_work: 10_000_000,
            relative_residual_tolerance: 1e-12,
        },
        max_interface_scalars: 100000,
        reference_resistance_scale: 1.0,
        relative_terminal_residual_tolerance: 1e-12,
        voltage_tolerance_v: 1e-9,
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_loaded<'t, A: HeatAssembly>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    template: &SpatialCurrentInduction<'_>,
    graph: &CircuitGraph,
    inputs: &[C64],
    impressed: &[[C64; 3]],
    controls: CircuitEnthalpySolve,
) -> Result<
    SpatialHarmonicCircuitThermalCandidate<'t, A::Sample>,
    CircuitEnthalpyError<SpatialHarmonicCircuitThermalError<A::Error>>,
> {
    circuit_heat_candidate(
        heat,
        cx,
        thermal,
        target,
        template,
        graph,
        &[3],
        &[],
        inputs,
        impressed,
        OMEGA,
        &regions(),
        1e-10,
        controls,
        terminal_controls(),
        circuit_controls(),
        assembly(),
        &[0.0; 2],
    )
}

struct Oracle {
    current: C64,
    voltage: C64,
    resistance: f64,
    field_power: f64,
    resistor_power: f64,
    circuit_supply: C64,
    impressed_supply: C64,
    magnetic_energy: f64,
    network_energy: f64,
}

fn oracle(ta: f64, tb: f64, drive: C64, frozen: [Option<f64>; 2]) -> Oracle {
    // Independent original cube integrals: K=4, int W_x=1/6, int |W|²=1/5.
    // Positive physical port current enters x=1; F=J/2. The external resistor
    // occupies the OTHER complete caloric cell, so its R uses Tb, not Ta.
    let sigma = 1.0 / (RHO * (1.0 + ALPHA * (frozen[0].unwrap_or(ta) - 300.0)));
    let resistance = RHO * (1.0 + ALPHA * (frozen[1].unwrap_or(tb) - 300.0));
    let f = excitation()[0][0].scale(0.5);
    let d = C64::new(4.0, -OMEGA * sigma * (1.0 / 5.0 - 1.0 / 36.0));
    let open = (C64::new(0.0, OMEGA) * f / d).scale(1.0 / 6.0);
    let z = C64::from_re(1.0 / sigma) + C64::new(0.0, -OMEGA / 36.0) / d;
    let load = C64::new(resistance, 1.0 / (OMEGA * LOAD_C) - OMEGA * LOAD_L);
    let current = (drive - open) / (z + load);
    let a = (f - current.scale(1.0 / 6.0)) / d;
    let induced = C64::new(0.0, OMEGA) * a;
    let voltage = current.scale(1.0 / sigma) + induced.scale(1.0 / 6.0);
    Oracle {
        current,
        voltage,
        resistance,
        field_power: 0.5
            * sigma
            * (voltage.norm_sq() - (voltage * induced.conj()).re / 3.0 + induced.norm_sq() / 5.0),
        resistor_power: 0.5 * resistance * current.norm_sq(),
        circuit_supply: (drive * current.conj()).scale(0.5),
        impressed_supply: (-induced * f.conj()).scale(0.5),
        magnetic_energy: a.norm_sq(),
        network_energy: 0.25 * current.norm_sq() * (LOAD_L + 1.0 / (OMEGA * OMEGA * LOAD_C)),
    }
}

fn reference(old: [f64; 2], drive: C64, frozen: [Option<f64>; 2]) -> [f64; 2] {
    let sink = |ta| {
        root(200.0, 1000.0, |tb| {
            4.0 * (tb - old[1])
                - DT * (oracle(ta, tb, drive, frozen).resistor_power - transfer(ta, tb, 0.5))
        })
    };
    let ta = root(200.0, 1000.0, |ta| {
        let tb = sink(ta);
        4.0 * (ta - old[0])
            - DT * (oracle(ta, tb, drive, frozen).field_power + transfer(ta, tb, 0.5))
    });
    [ta, sink(ta)]
}

#[test]
fn passive_and_driven_rlc_feedback_match_independent_two_temperature_energy_balances() {
    for driven in [false, true] {
        with_case(false, 0.5, |cx, _, heat, mut thermal, card| {
            with_resistor(heat, SINK, |loaded| {
                let prepared = loaded.prepare(cx, &thermal).unwrap();
                let original_resistor = prepared
                    .resolve_heat(cx, &thermal)
                    .unwrap()
                    .bulk
                    .resistors
                    .into_iter()
                    .next()
                    .unwrap();
                let problem = problem(card, BASIS);
                let template = template(&problem);
                let old_a = template.potential_wb().to_vec();
                let old_sources = template.conductors().to_vec();
                let graph = graph(driven);
                let drive = if driven {
                    C64::new(4.0, -3.0)
                } else {
                    C64::ZERO
                };
                let inputs = if driven { vec![drive] } else { vec![] };
                let mut temperatures = [300.0, 340.0];
                for interval in 1..=2 {
                    let expected = reference(temperatures, drive, [None; 2]);
                    let frozen_field =
                        reference(temperatures, drive, [Some(temperatures[0]), None]);
                    let frozen_load = reference(temperatures, drive, [None, Some(temperatures[1])]);
                    assert!((expected[0] - frozen_field[0]).abs() > 1e-4);
                    assert!((expected[0] - frozen_load[0]).abs() > 1e-4);
                    let exact = oracle(expected[0], expected[1], drive, [None; 2]);
                    let old_cells = thermal.cells().clone();
                    let next = advance_loaded(
                        &prepared,
                        cx,
                        &thermal,
                        f64::from(interval) * DT,
                        &template,
                        &graph,
                        &inputs,
                        &excitation(),
                        thermal_controls(),
                    )
                    .unwrap();
                    for i in 0..2 {
                        near(
                            next.thermal().temperature_k([BODY, SINK][i]).unwrap(),
                            expected[i],
                            3e-8,
                        );
                    }
                    let report = next.report();
                    let step = &report.electrical;
                    let circuit = &step.circuit;
                    let field = &circuit.field;
                    complex_near(field.terminals[1].incoming_current_a, exact.current, 3e-9);
                    complex_near(
                        field.terminals[1].voltage_to_reference_v,
                        exact.voltage,
                        3e-9,
                    );
                    complex_near(
                        circuit.circuit_supplied_power_va,
                        exact.circuit_supply,
                        5e-9,
                    );
                    complex_near(
                        circuit.impressed_supplied_power_va,
                        exact.impressed_supply,
                        5e-9,
                    );
                    near(circuit.mean_field_power_w, exact.field_power, 5e-9);
                    near(circuit.mean_resistor_power_w, exact.resistor_power, 5e-9);
                    near(
                        field.balance.mean_magnetic_energy_j,
                        exact.magnetic_energy,
                        1e-10,
                    );
                    near(
                        circuit.mean_network_stored_energy_j,
                        exact.network_energy,
                        1e-10,
                    );
                    if !driven {
                        assert!(field.balance.mean_terminal_power_w < 0.0);
                        assert!(circuit.sources.is_empty());
                    }
                    assert_eq!(field.cells.len(), 6);
                    assert!(field.cells.iter().all(|cell| cell.cell == BODY));
                    near(field.heat_j, DT * exact.field_power, 1e-8);
                    assert_eq!(step.resistors.len(), 1);
                    let resistor = &step.resistors[0];
                    assert_eq!(resistor.material.branch, 0);
                    let resistor_current = if driven {
                        exact.current
                    } else {
                        -exact.current
                    };
                    complex_near(resistor.current_a, resistor_current, 3e-9);
                    complex_near(
                        resistor.voltage_v,
                        resistor_current.scale(exact.resistance),
                        3e-9,
                    );
                    near(resistor.mean_power_w, exact.resistor_power, 5e-9);
                    near(resistor.heat_j, DT * exact.resistor_power, 1e-8);
                    let accepted_resistor_temperature = next.thermal().temperature_k(SINK).unwrap();
                    near(
                        resistor.material.conductor.resistance_ohm(),
                        RHO * (1.0 + ALPHA * (accepted_resistor_temperature - 300.0)),
                        1e-12,
                    );
                    near(
                        *resistor
                            .material
                            .conductor
                            .material()
                            .query_point()
                            .get("T")
                            .unwrap(),
                        accepted_resistor_temperature,
                        1e-12,
                    );
                    assert_eq!(
                        resistor.material.conductor.source_key(),
                        original_resistor.conductor.source_key()
                    );
                    assert_eq!(
                        resistor.material.conductor.material().card_identity(),
                        original_resistor.conductor.material().card_identity()
                    );
                    let total_heat = DT * (exact.field_power + exact.resistor_power);
                    near(step.heat_j, total_heat, 2e-8);
                    near(report.deposited_joule_j, total_heat, 2e-8);
                    near(
                        step.supplied_work_j,
                        DT * (exact.circuit_supply.re + exact.impressed_supply.re),
                        2e-8,
                    );
                    near(
                        next.thermal().energy_j() - thermal.energy_j(),
                        total_heat,
                        3e-8,
                    );
                    near(report.non_electrical_heat_j, 0.0, 2e-8);
                    near(report.total_energy_residual_j, 0.0, 2e-8);
                    assert_eq!(step.field_trials, 3);
                    assert_eq!(report.circuit_trials % 3, 0);
                    assert_eq!(thermal.cells(), &old_cells);
                    assert_eq!(template.potential_wb(), old_a);
                    assert_eq!(template.conductors(), old_sources);
                    assert_eq!(template.time_s(), 0.0);
                    temperatures = expected;
                    thermal = next.into_parts().0;
                }
            });
        });
    }
}

#[test]
fn missing_resistor_sources_and_double_electrical_heat_ownership_refuse_atomically() {
    with_case(false, 0.5, |cx, _, heat, thermal, card| {
        let problem = problem(card, BASIS);
        let template = template(&problem);
        let graph = graph(false);
        let old_cells = thermal.cells().clone();
        let old_a = template.potential_wb().to_vec();
        let missing = heat.prepare(cx, &thermal).unwrap();
        assert!(
            advance_loaded(
                &missing,
                cx,
                &thermal,
                DT,
                &template,
                &graph,
                &[],
                &excitation(),
                thermal_controls(),
            )
            .is_err()
        );
        // A full resolved conductor already owns BODY. A whole-cell graph
        // resistor cannot claim that same material volume and deposit it twice.
        with_resistor(heat, BODY, |duplicate| {
            let duplicate = duplicate.prepare(cx, &thermal).unwrap();
            assert!(
                advance_loaded(
                    &duplicate,
                    cx,
                    &thermal,
                    DT,
                    &template,
                    &graph,
                    &[],
                    &excitation(),
                    thermal_controls(),
                )
                .is_err()
            );
        });
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(template.potential_wb(), old_a);
        assert_eq!(template.time_s(), 0.0);
    });
}

#[test]
fn actual_field_work_budget_and_typed_heat_refusal_cancel_without_consuming_histories() {
    with_case(false, 0.5, |cx, gate, heat, thermal, card| {
        with_resistor(heat, SINK, |loaded| {
            let prepared = loaded.prepare(cx, &thermal).unwrap();
            let problem = problem(card, BASIS);
            let template = template(&problem);
            let graph = graph(true);
            let inputs = [C64::new(4.0, -3.0)];
            let old_cells = thermal.cells().clone();
            let old_a = template.potential_wb().to_vec();
            let tagged = Tagged {
                inner: &prepared,
                gate,
                calls: Cell::new(0),
                fail_at: Cell::new(None),
                cancel_at: Cell::new(None),
            };
            let mut limited = thermal_controls();
            limited.maximum_circuit_trials = 2;
            assert!(matches!(
                advance_loaded(
                    &tagged,
                    cx,
                    &thermal,
                    DT,
                    &template,
                    &graph,
                    &inputs,
                    &excitation(),
                    limited,
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            assert_eq!(
                tagged.calls.get(),
                1,
                "budget must refuse before the first proposal"
            );
            tagged.calls.set(0);
            tagged.fail_at.set(Some(3));
            assert!(matches!(
                advance_loaded(
                    &tagged,
                    cx,
                    &thermal,
                    DT,
                    &template,
                    &graph,
                    &inputs,
                    &excitation(),
                    thermal_controls(),
                ),
                Err(CircuitEnthalpyError::Source(
                    SpatialHarmonicCircuitThermalError::Heat(TaggedError::Refused)
                ))
            ));
            tagged.calls.set(0);
            tagged.fail_at.set(None);
            let retry = advance_loaded(
                &tagged,
                cx,
                &thermal,
                DT,
                &template,
                &graph,
                &inputs,
                &excitation(),
                thermal_controls(),
            )
            .unwrap();
            assert_eq!(retry.report().electrical.field_trials, 3);
            assert_eq!(retry.report().circuit_trials, (tagged.calls.get() - 1) * 3);
            assert_eq!(retry.report().source.network.resolution, tagged.calls.get());
            let expected = reference([300.0, 340.0], inputs[0], [None; 2]);
            near(
                retry.thermal().temperature_k(BODY).unwrap(),
                expected[0],
                3e-8,
            );
            tagged.calls.set(0);
            let zero = advance_loaded(
                &tagged,
                cx,
                &thermal,
                DT,
                &template,
                &graph,
                &inputs,
                &[[C64::ZERO; 3]; 6],
                thermal_controls(),
            )
            .unwrap();
            assert_eq!(zero.report().electrical.field_trials, 2);
            assert_eq!(zero.report().circuit_trials, (tagged.calls.get() - 1) * 2);
            tagged.calls.set(0);
            tagged.cancel_at.set(Some(4));
            assert!(matches!(
                advance_loaded(
                    &tagged,
                    cx,
                    &thermal,
                    DT,
                    &template,
                    &graph,
                    &inputs,
                    &excitation(),
                    thermal_controls(),
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(template.potential_wb(), old_a);
            assert_eq!(template.time_s(), 0.0);
        });
    });
}

//! Bulk induction and a separately heated source resistor close one real circuit.
#[path = "../../../../../fs-em/tests/spatial_induction/terminals/circuit/common.rs"]
mod common;
use super::*;
use common::{config, controls, initial};
use fs_couple::circuit_enthalpy::{
    ResistorCell, material::CellResistor,
    spatial_induction::terminals::circuit::SpatialTerminalCircuitThermalNetwork,
};
use fs_em::{
    spatial::vector_potential::induction::hysteretic::terminals::circuit::SpatialHystereticTerminalCircuit,
    winding::nonlinear::resistors::GraphResistorMaterial,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph};

const CAPACITY: f64 = 10.0;
const BULK: CellId = CellId(7);
const RESISTOR: CellId = CellId(8);

fn two_bodies(curve: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        [BULK, RESISTOR]
            .into_iter()
            .map(|id| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: 0,
                    body_id: id.0,
                    material_id: MaterialId(9),
                    mass_kg: 1.0,
                    specific_h_j_kg: CAPACITY * (300.0 - 200.0),
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [1.0, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 1.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect(),
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn endpoint_at_bulk_temperature(
    old: Endpoint,
    old_resistor_k: f64,
    bulk_k: f64,
    supplied_v: f64,
) -> (Endpoint, f64) {
    let sigma = resistivity(bulk_k).recip();
    // Exact graph-body caloric equation at fixed I for this manufactured
    // affine rho(T): C delta(T_R) = dt rho(T_R) I², with L=A=1.
    let resistor_temperature = |current: f64| {
        let capacity = CAPACITY - DT * current * current * 0.005 / 18.0;
        assert!(capacity > 0.0);
        old_resistor_k + DT * current * current * resistivity(old_resistor_k) / capacity
    };
    let h = common::bracket(-29.0, 29.0, |h| {
        let next = reference::at_field(old, sigma, DT, h);
        next.voltage + resistivity(resistor_temperature(next.incoming)) * next.incoming - supplied_v
    });
    let next = reference::at_field(old, sigma, DT, h);
    (next, resistor_temperature(next.incoming))
}

#[test]
fn bulk_and_graph_material_feedback_heat_distinct_finite_receivers_without_numerical_loss() {
    let card = data::card_with_heat_capacity(500.0, CAPACITY);
    let curve = caloric(&card);
    let problem = cube(&card);
    let query = point(300.0);
    let [bulk] = cells(&card, &curve, &query);
    let [mut graph_body] = cells(&card, &curve, &query);
    graph_body.cell = RESISTOR;
    let material_cells = [bulk, graph_body];
    let regions = regions();
    let routes = [ResistorCell {
        branch: 1,
        cell: RESISTOR,
    }];
    let resistors = [CellResistor {
        destination: routes[0],
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let network = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &material_cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &routes,
    };
    let field = SpatialHystereticTerminalInduction::new(
        data::initial(&problem, Some(conductor(&card, 300.0))),
        &electrodes(&problem),
        SpatialInductionLimits {
            maximum_unknowns: 2000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 100_000,
        },
        &mut || false,
    )
    .unwrap();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (
                2,
                1,
                Branch::Resistor {
                    ohms: resistivity(300.0),
                },
            ),
        ],
        transformers: vec![],
    };
    let mut circuit = SpatialHystereticTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![], vec![], vec![0.0]),
        config(0.7),
        &mut || false,
    )
    .unwrap()
    .with_graph_resistors(
        vec![GraphResistorMaterial {
            branch: 1,
            conductor: conductor(&card, 300.0)
                .into_uniform_conductor(1.0, 1.0)
                .unwrap(),
        }],
        &mut || false,
    )
    .unwrap();
    let mut thermal = two_bodies(&curve);
    let mut old = Endpoint::default();
    let (mut old_bulk_k, mut old_resistor_k) = (300.0, 300.0);
    for (index, drive) in [28.0, 0.0, -15.0].into_iter().enumerate() {
        let bulk_k = common::bracket(old_bulk_k, 480.0, |t| {
            let (next, _) = endpoint_at_bulk_temperature(old, old_resistor_k, t, drive);
            CAPACITY * (t - old_bulk_k)
                - joule(old, next, resistivity(t).recip(), DT)
                - PIN * (next.z - old.z).abs()
        });
        let (next, resistor_k) = endpoint_at_bulk_temperature(old, old_resistor_k, bulk_k, drive);
        let bulk_joule = joule(old, next, resistivity(bulk_k).recip(), DT);
        let pinning = PIN * (next.z - old.z).abs();
        let graph_heat = DT * resistivity(resistor_k) * next.incoming.powi(2);
        let change = stored(next.h, next.z) - stored(old.h, old.z);
        let numerical = next.h * (next.b - old.b) - change - pinning;
        let target = (index + 1) as f64 * DT;
        let before = (
            thermal.cells().clone(),
            circuit.native_state().to_vec(),
            circuit.field().potential_wb().to_vec(),
            circuit.field().cell_states().to_vec(),
            circuit.time_s(),
        );
        let mut polls = 0;
        let proposal = network
            .step_hysteretic_material_candidate(
                &thermal,
                target,
                &circuit,
                &[drive],
                &[[0.0; 3]; 6],
                &resistors,
                &bases,
                thermal_controls(),
                controls(field_controls()),
                &[0.0; 2],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let step = &proposal.report().coupled.electrical;
        near(
            proposal.thermal().temperature_k(BULK).unwrap(),
            bulk_k,
            5e-7,
        );
        near(
            proposal.thermal().temperature_k(RESISTOR).unwrap(),
            resistor_k,
            5e-7,
        );
        assert_memory(proposal.circuit().field().cell_states(), next);
        near(step.terminal_current_a[0], next.incoming, 5e-7);
        near(step.terminal_voltage_v[0], next.voltage, 5e-7);
        near(step.field.joule_heat_j, bulk_joule, 6e-7);
        near(step.field.hysteresis_heat_j, pinning, 5e-8);
        near(step.field.stored_change_j, change, 5e-8);
        near(step.numerical_dissipation_j, numerical, 5e-8);
        near(step.graph_resistors[0].heat_j, graph_heat, 6e-7);
        near(step.resistor_heat_j, graph_heat, 6e-7);
        near(step.heat_j, bulk_joule + pinning + graph_heat, 8e-7);
        near(CAPACITY * (bulk_k - old_bulk_k), bulk_joule + pinning, 5e-8);
        near(CAPACITY * (resistor_k - old_resistor_k), graph_heat, 5e-8);
        near(
            proposal.deposited_heat_j(),
            bulk_joule + pinning + graph_heat,
            8e-7,
        );
        near(
            proposal.report().coupled.thermal.stored_energy_change_j,
            bulk_joule + pinning + graph_heat,
            8e-7,
        );
        near(proposal.report().total_energy_residual_j, 0.0, 2e-8);
        near(step.external_work_j, DT * drive * next.incoming, 8e-7);
        for source in proposal.circuit().field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m(), resistivity(bulk_k), 1e-11);
            near(source.material().query().axes()["T"], bulk_k, 5e-7);
        }
        assert_eq!(step.graph_materials.len(), 1);
        near(
            step.graph_materials[0].conductor.resistance_ohm(),
            resistivity(resistor_k),
            1e-11,
        );
        near(
            step.graph_materials[0].conductor.material().query().axes()["T"],
            resistor_k,
            5e-7,
        );
        if index == 0 {
            assert!(bulk_k > 320.0 && resistor_k > 320.0);
            assert!((bulk_k - resistor_k).abs() > 1.0);
            assert!(pinning > 1.0 && numerical > 1.0 && graph_heat > 200.0);
            assert!(step.cancelled_reference_heat_j > 1.0);
            // Neither numerical endpoint loss nor the canceled reference
            // resistor may be paid into either finite caloric receiver.
            assert!((proposal.deposited_heat_j() - step.heat_j - numerical).abs() > 1.0);
            assert!(polls > 20);
            let mut seen = 0;
            assert!(matches!(
                network.step_hysteretic_material_candidate(
                    &thermal,
                    target,
                    &circuit,
                    &[drive],
                    &[[0.0; 3]; 6],
                    &resistors,
                    &bases,
                    thermal_controls(),
                    controls(field_controls()),
                    &[0.0; 2],
                    || {
                        seen += 1;
                        seen >= polls - 1
                    },
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let wrong = circuit
                .with_graph_resistors(
                    vec![GraphResistorMaterial {
                        branch: 1,
                        conductor: conductor(&card, 301.0)
                            .into_uniform_conductor(1.0, 1.0)
                            .unwrap(),
                    }],
                    &mut || false,
                )
                .unwrap();
            assert!(
                network
                    .step_hysteretic_material_candidate(
                        &thermal,
                        target,
                        &wrong,
                        &[drive],
                        &[[0.0; 3]; 6],
                        &resistors,
                        &bases,
                        thermal_controls(),
                        controls(field_controls()),
                        &[0.0; 2],
                        || false,
                    )
                    .is_err()
            );
        }
        assert_eq!(thermal.cells(), &before.0);
        assert_eq!(circuit.native_state(), before.1);
        assert_eq!(circuit.field().potential_wb(), before.2);
        assert_eq!(circuit.field().cell_states(), before.3);
        assert_eq!(circuit.time_s(), before.4);
        old = next;
        old_bulk_k = bulk_k;
        old_resistor_k = resistor_k;
        let (next_thermal, next_circuit, _) = proposal.into_parts();
        thermal = next_thermal;
        circuit = next_circuit;
    }
    assert!(old.z < 0.0 && old_bulk_k > old_resistor_k);
}

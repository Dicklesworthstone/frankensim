//! Manufactured source-backed common-core heat, not measured material data.
use super::*;
use crate::{
    circuit_enthalpy::{
        CircuitEnthalpyError,
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        material::ConductionLink,
    },
    material_circuit::hysteretic::shared::{
        SharedHystereticCircuit, SharedHystereticCircuitInitial,
        thermal::ThermalSharedHystereticCircuit,
    },
};
use fs_exec::{CancelGate, Cx};
use fs_material::{
    field_storage::hysteresis::shared::{SharedHystereticCore, SharedHystereticWinding},
    state_point::THERMAL_CONDUCTIVITY_DIMS,
};

const CORE: CellId = CellId(10);
const COILS: [CellId; 2] = [CellId(11), CellId(12)];
const PORTS: [f64; 3] = [5.0, 17.0, -3.0];
const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity),
    (4, ElectricalPropertyBasis::Resistivity),
];

fn conductive(card: &MaterialCard, value: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = card.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                200.0,
                2000.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured shared-core finite exchange".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn with_system<R>(
    pinning: f64,
    alpha: f64,
    core_upper: f64,
    conductivity: Option<f64>,
    initial_temperatures: [f64; 3],
    run: impl FnOnce(
        &DepositionSimulation<'_>,
        &ThermalSharedHystereticCircuit<'_, '_>,
        &MaterialNetwork<'_>,
    ) -> R,
) -> R {
    let (a, ca) = source::core_material(pinning, core_upper);
    let (b, cb) = material(0.7, alpha);
    let (c, cc) = material(1.1, alpha);
    let mut cards = [a, b, c];
    let mut curves = [ca, cb, cc];
    if let Some(k) = conductivity {
        for i in 0..3 {
            let (card, curve) = conductive(&cards[i], k);
            cards[i] = card;
            curves[i] = curve;
        }
    }
    let mut events = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(10 + i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: initial_temperatures[i] - 200.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
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
    let links = if conductivity.is_some() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 0,
                a: CORE,
                b: COILS[0],
                area_m2: 1.0,
                conductance_w_k: 999.0,
            }),
        });
        vec![ConductionLink {
            link: 0,
            distance_a_m: 0.5,
            distance_b_m: 0.5,
        }]
    } else {
        vec![]
    };
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
            (MaterialId(2), curves[2].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 1,
            max_events: 4,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1, 2].map(|i| MaterialCell {
        cell: CellId(10 + i as u64),
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: conductivity.map(|_| (&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [1, 4].map(|branch| CellResistor {
        destination: ResistorCell {
            branch,
            cell: if branch == 1 { COILS[0] } else { COILS[1] },
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    });
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let sources = network.resolve(&thermal, || false).unwrap().resistors;
    let core = SharedHystereticCore::new(
        source::sourced_core(&cards[0]),
        vec![
            SharedHystereticWinding {
                turns_ratio: 2.0,
                leakage_h: 0.3,
            },
            SharedHystereticWinding {
                turns_ratio: -1.0,
                leakage_h: 0.4,
            },
        ],
    )
    .unwrap();
    let initial = SharedHystereticCircuitInitial {
        magnetic: core.initial_state(),
        remaining_fluxes_wb: vec![],
        charges_c: vec![],
        external_inputs: vec![0.0; 3],
        time_s: 0.0,
    };
    let graph = CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: f64::NAN }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 2 }),
            (3, 4, Branch::Resistor { ohms: f64::NAN }),
            (4, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    };
    let circuit = SharedHystereticCircuit::new(graph, vec![2, 5], core, sources, initial).unwrap();
    let driver = ThermalSharedHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &BASES,
        CORE,
        [290.0, core_upper],
        || false,
    )
    .unwrap();
    run(&thermal, &driver, &network)
}

fn thermal_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 3,
        ..heat_controls()
    }
}

#[test]
fn one_common_core_heat_and_each_original_winding_heat_advance_actual_temperatures() {
    let run = |pinning, alpha| {
        with_system(
            pinning,
            alpha,
            310.0,
            None,
            [300.0; 3],
            |thermal, driver, _| {
                let before = driver.circuit().clone();
                let cells_before = thermal.cells().clone();
                let next = driver
                    .step_candidate(
                        thermal,
                        0.1,
                        &PORTS,
                        thermal_controls(),
                        controls(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                let report = next.report();
                let step = &report.electrical;
                let core_heat = step.magnetic.core.core_heat_j;
                assert!(core_heat > 0.0);
                near(
                    next.thermal().temperature_k(CORE).unwrap() - 300.0,
                    core_heat,
                    2e-8,
                );
                assert_eq!(
                    step.circuit
                        .resistors
                        .iter()
                        .map(|r| r.material.branch)
                        .collect::<Vec<_>>(),
                    vec![1, 4]
                );
                for (loss, cell) in step.circuit.resistors.iter().zip(COILS) {
                    near(
                        next.thermal().temperature_k(cell).unwrap() - 300.0,
                        loss.heat_j,
                        2e-8,
                    );
                }
                near(
                    report.deposited_joule_j,
                    core_heat + step.circuit.heat_j,
                    1e-10,
                );
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    report.deposited_joule_j,
                    3e-8,
                );
                near(report.total_energy_residual_j, 0.0, 2e-8);
                assert_eq!(step.numerical_dissipation_j, 0.0);
                assert_eq!(next.circuit().circuit().branch_indices(), &[2, 5]);
                let source = next.circuit().circuit().core().core().source().unwrap();
                assert_eq!(
                    source.material().query(),
                    point(next.thermal().temperature_k(CORE).unwrap())
                );
                assert_eq!(
                    source.material().card_identity(),
                    before
                        .core()
                        .core()
                        .source()
                        .unwrap()
                        .material()
                        .card_identity()
                );
                assert_eq!(driver.circuit(), &before);
                assert_eq!(thermal.cells(), &cells_before);

                // A second interval keeps the original common polarization. Its
                // core receiver gets only the new interval's physical pinning heat.
                let memory = next.circuit().circuit().magnetic_state().clone();
                let continued = next
                    .circuit()
                    .step_candidate(
                        next.thermal(),
                        0.2,
                        &PORTS,
                        thermal_controls(),
                        controls(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                assert_eq!(next.circuit().circuit().magnetic_state(), &memory);
                assert!(
                    continued
                        .circuit()
                        .circuit()
                        .magnetic_state()
                        .core_state()
                        .polarizations_t()[0]
                        > memory.core_state().polarizations_t()[0]
                );
                near(
                    continued.thermal().temperature_k(CORE).unwrap() - 300.0,
                    core_heat + continued.report().electrical.magnetic.core.core_heat_j,
                    4e-8,
                );
                (
                    continued.circuit().circuit().magnetic_state().currents_a()[0],
                    core_heat,
                    continued.thermal().temperature_k(COILS[0]).unwrap(),
                )
            },
        )
    };
    let first = run(1.0, 0.005);
    let harder = run(2.0, 0.005);
    let fixed_resistance = run(1.0, 0.0);
    assert!((harder.0 - first.0).abs() > 1e-4);
    assert!((harder.1 - first.1).abs() > 1e-4);
    assert!((harder.2 - first.2).abs() > 1e-4);
    assert!(fixed_resistance.0 > first.0);
}

#[test]
fn native_heat_assembly_cools_common_core_and_changes_winding_feedback() {
    let run = |conductivity| {
        with_system(
            1.0,
            0.005,
            310.0,
            Some(conductivity),
            [300.0, 290.0, 300.0],
            |thermal, driver, network| {
                let electrical = ElectricalContactNetwork {
                    bulk: network,
                    contacts: &[],
                    electrical_bases: &BASES,
                };
                let assembly = MaterialHeatNetwork {
                    electrical: &electrical,
                    thermal_contacts: &[],
                    radiation: None,
                };
                let gate = CancelGate::new();
                let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
                pool.scope(|arena| {
                    let cx = Cx::new(
                        &gate,
                        arena,
                        fs_exec::StreamKey {
                            seed: 43,
                            kernel_id: 1,
                            tile: 0,
                            iteration: 0,
                        },
                        fs_exec::Budget::INFINITE,
                        fs_exec::ExecMode::Deterministic,
                    );
                    let prepared = assembly.prepare(&cx, thermal).unwrap();
                    let next = driver
                        .step_candidate_with_heat_assembly(
                            &cx,
                            &prepared,
                            thermal,
                            0.05,
                            &PORTS,
                            thermal_controls(),
                            controls(),
                            &[0.0; 3],
                        )
                        .unwrap();
                    let [tc, t1, t2] = [CORE, COILS[0], COILS[1]]
                        .map(|cell| next.thermal().temperature_k(cell).unwrap());
                    let step = &next.report().electrical;
                    let transfer = 0.05 * conductivity * (tc - t1);
                    near(tc - 300.0, step.magnetic.core.core_heat_j - transfer, 3e-8);
                    near(
                        t1 - 290.0,
                        step.circuit.resistors[0].heat_j + transfer,
                        3e-8,
                    );
                    near(t2 - 300.0, step.circuit.resistors[1].heat_j, 3e-8);
                    near(next.report().non_electrical_heat_j, 0.0, 2e-8);
                    near(next.report().total_energy_residual_j, 0.0, 2e-8);
                    assert_eq!(next.report().source.material.links.len(), 1);
                    near(
                        next.report().source.material.links[0].conductance_w_k,
                        conductivity,
                        1e-12,
                    );
                    let source = next.circuit().circuit().core().core().source().unwrap();
                    assert_eq!(source.material().query(), point(tc));
                    (tc, t1, next.circuit().circuit().reference_flux_linkage_wb())
                })
            },
        )
    };
    let weak = run(0.25);
    let strong = run(4.0);
    assert!(strong.0 < weak.0 && strong.1 > weak.1);
    assert!((strong.2 - weak.2).abs() > 1e-8);
}

#[test]
fn late_core_temperature_refusal_and_cancellation_preserve_all_accepted_histories() {
    with_system(
        1.0,
        0.005,
        300.01,
        None,
        [300.0; 3],
        |thermal, driver, _| {
            let before = driver.circuit().clone();
            let cells = thermal.cells().clone();
            assert!(
                driver
                    .step_candidate(
                        thermal,
                        0.05,
                        &PORTS,
                        thermal_controls(),
                        controls(),
                        &[0.0; 3],
                        || false
                    )
                    .is_err()
            );
            assert_eq!(driver.circuit(), &before);
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(thermal.time_s(), 0.0);
            let checks = std::cell::Cell::new(0usize);
            let retry = driver
                .step_candidate(
                    thermal,
                    1e-5,
                    &PORTS,
                    thermal_controls(),
                    controls(),
                    &[0.0; 3],
                    || {
                        checks.set(checks.get() + 1);
                        false
                    },
                )
                .unwrap();
            let finish = checks.get();
            assert!(finish > 10);
            checks.set(0);
            assert!(matches!(
                driver.step_candidate(
                    thermal,
                    1e-5,
                    &PORTS,
                    thermal_controls(),
                    controls(),
                    &[0.0; 3],
                    || {
                        checks.set(checks.get() + 1);
                        checks.get() >= finish
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(driver.circuit(), &before);
            assert_eq!(thermal.cells(), &cells);
            let repeated = driver
                .step_candidate(
                    thermal,
                    1e-5,
                    &PORTS,
                    thermal_controls(),
                    controls(),
                    &[0.0; 3],
                    || false,
                )
                .unwrap();
            assert_eq!(retry.thermal().cells(), repeated.thermal().cells());
            assert_eq!(retry.circuit().circuit(), repeated.circuit().circuit());
        },
    );
}

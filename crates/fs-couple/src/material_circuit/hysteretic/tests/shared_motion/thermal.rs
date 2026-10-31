//! Finite source-resolved heat, common magnetic memory and native gap motion.
//! Authored values test integration and physical accounting, not alloy fidelity.
use super::*;
use crate::{
    circuit_enthalpy::{
        CircuitEnthalpyError,
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        material::ConductionLink,
    },
    material_circuit::hysteretic::shared::motion::thermal::{
        SharedHystereticGapHeatInput, ThermalSharedHystereticGapMotion,
    },
};
use fs_exec::{CancelGate, Cx};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;

const CORE: CellId = CellId(10);
const COILS: [CellId; 2] = [CellId(11), CellId(12)];
const DRIVE: [f64; 3] = [5.0, 17.0, -3.0];
const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity),
    (4, ElectricalPropertyBasis::Resistivity),
];

#[derive(Clone, Copy)]
struct Case {
    pinning: f64,
    alpha: f64,
    core_upper: f64,
    conductivity: Option<f64>,
    temperature: [f64; 3],
    supported: bool,
}
impl Default for Case {
    fn default() -> Self {
        Self {
            pinning: 1.0,
            alpha: 0.005,
            core_upper: 310.0,
            conductivity: None,
            temperature: [300.0; 3],
            supported: false,
        }
    }
}

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
                source: "authored moving common-core heat transport".into(),
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
    case: Case,
    run: impl FnOnce(
        &DepositionSimulation<'_>,
        &ThermalSharedHystereticGapMotion<'_, '_>,
        &MaterialNetwork<'_>,
    ) -> R,
) -> R {
    let (a, ca) = source::core_material(case.pinning, case.core_upper);
    let (b, cb) = material(RESISTANCE[0], case.alpha);
    let (c, cc) = material(RESISTANCE[1], case.alpha);
    let (mut cards, mut curves) = ([a, b, c], [ca, cb, cc]);
    if let Some(k) = case.conductivity {
        for i in 0..3 {
            let (card, curve) = conductive(&cards[i], k);
            cards[i] = card;
            curves[i] = curve;
        }
    }
    let mut events: Vec<_> = (0..3)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(10 + i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: case.temperature[i] - 200.0,
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
        .collect();
    let links = if case.conductivity.is_some() {
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
            max_dt_s: 0.1,
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
        conductivity: case
            .conductivity
            .map(|_| (&key, MaterialPropertySelection::SingleClaimOnly)),
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
    let core = with_gap(shared(source::sourced_core(&cards[0])), [0.2, 2.0]);
    let initial = SharedHystereticCircuitInitial {
        magnetic: core.initial_state(),
        remaining_fluxes_wb: vec![],
        charges_c: vec![],
        external_inputs: vec![0.0; 3],
        time_s: 0.0,
    };
    let (graph, branches) = graph(true);
    let circuit = SharedHystereticCircuit::new(graph, branches, core, sources, initial).unwrap();
    let motion =
        SharedHystereticGapMotion::new(circuit, body(MASS), case.supported.then(support), 0.2)
            .unwrap();
    let driver = ThermalSharedHystereticGapMotion::bind(
        motion,
        &thermal,
        &network,
        &BASES,
        CORE,
        [290.0, case.core_upper],
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
fn input(force: f64) -> SharedHystereticGapHeatInput {
    SharedHystereticGapHeatInput {
        external_force_n: force,
        controls: policy(),
    }
}
fn temperatures(thermal: &DepositionSimulation<'_>) -> [f64; 3] {
    [CORE, COILS[0], COILS[1]].map(|cell| thermal.temperature_k(cell).unwrap())
}

// Independent source formulas: one scalar play branch, series gap, physical
// winding leakage and a sourced guided body. No production energy evaluator.
fn independent_storage(motion: &SharedHystereticGapMotion) -> f64 {
    let magnetic = motion.circuit().magnetic_state();
    let h = magnetic.core_state().field_a_m();
    let z = magnetic.core_state().polarizations_t()[0];
    let phi = 0.01 * h + z;
    let core = 0.005 * h * h + z * z / 0.4;
    let leakage: f64 = magnetic
        .currents_a()
        .iter()
        .zip(LEAKAGE)
        .map(|(i, l)| 0.5 * l * i * i)
        .sum();
    let gap = 1.25 * motion.gap_m() * phi * phi;
    let (support_mass, elastic) = if motion.support().is_some() {
        (
            3.0 * 0.02 * 2.0 / 3.0,
            0.5 * (motion.gap_m() - 0.68).powi(2),
        )
    } else {
        (0.0, 0.0)
    };
    near(motion.generalized_mass_kg(), MASS + support_mass, 1e-12);
    core + leakage + gap + 0.5 * (MASS + support_mass) * motion.velocity_m_s().powi(2) + elastic
}

#[test]
fn finite_winding_heat_changes_retained_gap_motion_with_one_common_core_loss() {
    let run = |case: Case| {
        with_system(case, |thermal, driver, _| {
            let initial = driver.circuit().clone();
            let first = driver
                .step_candidate(
                    thermal,
                    0.1,
                    &DRIVE,
                    thermal_controls(),
                    input(0.4),
                    &[0.0; 3],
                    || false,
                )
                .unwrap();
            let second = first
                .circuit()
                .step_candidate(
                    first.thermal(),
                    0.2,
                    &DRIVE,
                    thermal_controls(),
                    input(0.4),
                    &[0.0; 3],
                    || false,
                )
                .unwrap();
            for (before, after, thermal_before, thermal_after, report) in [
                (
                    &initial,
                    first.circuit().circuit(),
                    thermal,
                    first.thermal(),
                    first.report(),
                ),
                (
                    first.circuit().circuit(),
                    second.circuit().circuit(),
                    first.thermal(),
                    second.thermal(),
                    second.report(),
                ),
            ] {
                let t0 = temperatures(thermal_before);
                let t1 = temperatures(thermal_after);
                let step = &report.electrical;
                let z0 = before
                    .circuit()
                    .magnetic_state()
                    .core_state()
                    .polarizations_t()[0];
                let z1 = after
                    .circuit()
                    .magnetic_state()
                    .core_state()
                    .polarizations_t()[0];
                let pinning_heat = case.pinning * (z1 - z0).abs();
                assert!(pinning_heat > 0.0);
                near(step.magnetic.core.core_heat_j, pinning_heat, 2e-10);
                near(t1[0] - t0[0], pinning_heat, 2e-8);
                for (i, loss) in step.circuit.resistors.iter().enumerate() {
                    let resistance = RESISTANCE[i] * (1.0 + case.alpha * (t1[i + 1] - 300.0));
                    near(loss.material.conductor.resistance_ohm(), resistance, 2e-12);
                    assert_eq!(loss.material.conductor.material().query(), point(t1[i + 1]));
                    near(
                        loss.heat_j,
                        0.1 * resistance * loss.current_a.powi(2),
                        2e-10,
                    );
                    near(t1[i + 1] - t0[i + 1], loss.heat_j, 2e-8);
                }
                let dg = after.gap_m() - before.gap_m();
                assert!(dg > 0.0 && step.field_work_on_body_j < 0.0);
                near(step.mechanical_supply_j, 0.4 * dg, 2e-12);
                near(
                    step.kinetic_energy_j[1] - step.kinetic_energy_j[0] + step.elastic_energy_j[1]
                        - step.elastic_energy_j[0],
                    0.4 * dg + step.field_work_on_body_j,
                    2e-8,
                );
                near(step.stored_before_j, independent_storage(before), 2e-10);
                near(step.stored_after_j, independent_storage(after), 2e-10);
                near(
                    report.deposited_joule_j,
                    pinning_heat + step.circuit.heat_j,
                    2e-10,
                );
                near(report.non_electrical_heat_j, 0.0, 2e-12);
                near(
                    independent_storage(after) - independent_storage(before)
                        + thermal_after.energy_j()
                        - thermal_before.energy_j(),
                    step.electrical_supply_j + 0.4 * dg,
                    2e-8,
                );
                near(report.total_energy_residual_j, 0.0, 2e-8);
                assert!(report.circuit_trials > step.circuit_trials);
                assert_eq!(step.circuit_trials, step.mechanical_trials + 3);
                let source = after.circuit().core().core().source().unwrap();
                assert_eq!(source.material().query(), point(t1[0]));
                assert_eq!(
                    source.material().card_identity(),
                    before
                        .circuit()
                        .core()
                        .core()
                        .source()
                        .unwrap()
                        .material()
                        .card_identity()
                );
                assert_eq!(after.body(), before.body());
                assert_eq!(after.support(), before.support());
            }
            assert_eq!(driver.circuit(), &initial);
            assert_eq!(temperatures(thermal), case.temperature);
            assert_eq!(thermal.time_s(), 0.0);
            (
                second.circuit().circuit().gap_m(),
                second.circuit().circuit().velocity_m_s(),
                second
                    .circuit()
                    .circuit()
                    .circuit()
                    .magnetic_state()
                    .currents_a()[0],
            )
        })
    };
    let heated = run(Case::default());
    let fixed = run(Case {
        alpha: 0.0,
        ..Case::default()
    });
    let harder = run(Case {
        pinning: 2.0,
        ..Case::default()
    });
    let supported = run(Case {
        supported: true,
        ..Case::default()
    });
    assert!((heated.2 - fixed.2).abs() > 1e-6);
    assert!((heated.1 - fixed.1).abs() > 1e-9);
    assert!((heated.1 - harder.1).abs() > 1e-7);
    assert!((heated.0 - supported.0).abs() > 1e-6);
}

#[test]
fn complete_heat_assembly_cools_the_original_core_and_changes_actual_motion() {
    let run = |k: f64| {
        with_system(
            Case {
                conductivity: Some(k),
                temperature: [300.0, 290.0, 300.0],
                ..Case::default()
            },
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
                            seed: 47,
                            kernel_id: 3,
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
                            0.1,
                            &DRIVE,
                            thermal_controls(),
                            input(0.4),
                            &[0.0; 3],
                        )
                        .unwrap();
                    let t = temperatures(next.thermal());
                    let step = &next.report().electrical;
                    let transfer = 0.1 * k * (t[0] - t[1]);
                    near(
                        t[0] - 300.0,
                        step.magnetic.core.core_heat_j - transfer,
                        2e-8,
                    );
                    near(
                        t[1] - 290.0,
                        step.circuit.resistors[0].heat_j + transfer,
                        2e-8,
                    );
                    near(t[2] - 300.0, step.circuit.resistors[1].heat_j, 2e-8);
                    near(
                        next.report().source.material.links[0].conductance_w_k,
                        k,
                        1e-12,
                    );
                    near(next.report().deposited_joule_j, step.heat_j, 2e-10);
                    near(next.report().non_electrical_heat_j, 0.0, 2e-12);
                    near(next.report().total_energy_residual_j, 0.0, 2e-8);
                    let motion = next.circuit().circuit();
                    assert_eq!(
                        motion
                            .circuit()
                            .core()
                            .core()
                            .source()
                            .unwrap()
                            .material()
                            .query(),
                        point(t[0])
                    );
                    near(
                        independent_storage(motion) - independent_storage(driver.circuit())
                            + next.thermal().energy_j()
                            - thermal.energy_j(),
                        step.electrical_supply_j + step.mechanical_supply_j,
                        2e-8,
                    );
                    (t, motion.velocity_m_s())
                })
            },
        )
    };
    let weak = run(0.25);
    let strong = run(4.0);
    assert!(strong.0[0] < weak.0[0] && strong.0[1] > weak.0[1]);
    assert!((strong.1 - weak.1).abs() > 1e-9);
}

#[test]
fn source_refusal_budget_and_late_cancel_preserve_momentum_memory_and_heat() {
    with_system(Case::default(), |thermal, driver, network| {
        let original = driver.circuit().clone();
        let cells = thermal.cells().clone();
        let stale = original.with_core_source_query(&point(301.0)).unwrap();
        assert!(
            ThermalSharedHystereticGapMotion::bind(
                stale,
                thermal,
                network,
                &BASES,
                CORE,
                [290.0, 310.0],
                || false,
            )
            .is_err()
        );
        let mut polls = 0usize;
        let complete = driver
            .step_candidate(
                thermal,
                0.1,
                &DRIVE,
                thermal_controls(),
                input(0.4),
                &[0.0; 3],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let mut limited = thermal_controls();
        limited.maximum_circuit_trials = 4;
        assert!(matches!(
            driver.step_candidate(thermal, 0.1, &DRIVE, limited, input(0.4), &[0.0; 3], || {
                false
            },),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        let mut visited = 0usize;
        assert!(polls > 10);
        assert!(matches!(
            driver.step_candidate(
                thermal,
                0.1,
                &DRIVE,
                thermal_controls(),
                input(0.4),
                &[0.0; 3],
                || {
                    visited += 1;
                    visited >= polls
                },
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        let exact = CircuitEnthalpySolve {
            maximum_circuit_trials: complete.report().circuit_trials,
            ..thermal_controls()
        };
        let retry = driver
            .step_candidate(thermal, 0.1, &DRIVE, exact, input(0.4), &[0.0; 3], || false)
            .unwrap();
        assert_eq!(retry.circuit().circuit(), complete.circuit().circuit());
        assert_eq!(retry.thermal().cells(), complete.thermal().cells());
        assert_eq!(driver.circuit(), &original);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
    with_system(
        Case {
            core_upper: 300.01,
            ..Case::default()
        },
        |thermal, driver, _| {
            let original = driver.circuit().clone();
            let cells = thermal.cells().clone();
            assert!(
                driver
                    .step_candidate(
                        thermal,
                        0.1,
                        &DRIVE,
                        thermal_controls(),
                        input(0.4),
                        &[0.0; 3],
                        || false,
                    )
                    .is_err()
            );
            assert_eq!(driver.circuit(), &original);
            assert_eq!(thermal.cells(), &cells);
        },
    );
}

//! Real source cards, finite electrical storage and both independent gas owners.
use super::material::{MaterialContact, step_candidate as material_step};
use super::*;
use crate::cavity_phs::piston::tests::{controls as piston_controls, drive, gas};
use crate::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::material::*;
use crate::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use crate::material_circuit::MaterialCircuitSolve;
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}
fn material(k: f64, r: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, a, b) in [
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            r * 0.0025 * 0.6,
            r * 0.0025 * 1.8,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            k,
            2.0 * k,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            500.0,
            500.0,
        ),
        ("density", Density::DIMS, 4000.0, 4000.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, a), (800.0, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 200.0, 800.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic multi-chamber material".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    (card, curve)
}
fn cell(i: usize) -> CellId {
    CellId(if i == 0 { 30 } else { 10 })
}
fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    p: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: cell(i),
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: p,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn resistors() -> [CellResistor; 2] {
    [0, 1].map(|i| CellResistor {
        destination: ResistorCell {
            branch: i + 1,
            cell: cell(i),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 0.01,
        area_m2: 0.000025,
    })
}
fn initial<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    temperatures: [f64; 2],
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell(i),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 0.001,
                specific_h_j_kg: 500.0 * (temperatures[i] - 200.0),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 0.01, 0.0, 0.0],
                    volume_m3: 0.001 / 4000.0,
                    surface_area_m2: 0.04,
                },
                boundary: Boundary {
                    ambient_k: 400.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: cell(0),
            b: cell(1),
            area_m2: 0.002,
            conductance_w_k: 999.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
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
            max_dt_s: 0.001,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn circuit(
    binding: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    charge: f64,
) -> MaterialCircuit {
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 0.01 }),
                (1, 2, Branch::Resistor { ohms: 0.0 }),
                (2, 0, Branch::Resistor { ohms: 0.0 }),
            ],
            transformers: vec![],
        },
        binding.resolve(thermal, || false).unwrap().resistors,
        vec![],
        vec![charge],
        &[],
    )
    .unwrap()
}
fn solve() -> NetworkSolve {
    NetworkSolve {
        network: CircuitEnthalpySolve {
            maximum_cells: 2,
            maximum_circuit_trials: 300,
            maximum_newton_iterations: 32,
            difference_step_j: 1e-4,
            cell_energy_tolerance_j: 1e-7,
            thermal_energy_tolerance_j: 1e-7,
            total_energy_tolerance_j: 1e-7,
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: 1e-7,
                power_tolerance_w: 1e-7,
            },
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-8,
                relative_tolerance: 0.0,
                forcing_minimum: 1e-4,
                linear_restart: 4,
                ..NewtonKrylovConfig::default()
            },
        },
        maximum_chambers: 2,
        maximum_contacts: 4,
        maximum_boundary_solves: 600,
        total_energy_tolerance_j: 1e-6,
    }
}
fn chamber(id: u64, piston: &GasPiston, dt: f64) -> Chamber<'_> {
    Chamber {
        id,
        piston,
        input: PistonInput {
            duration_s: dt,
            ..drive(0.0)
        },
        piston_solve: piston_controls(),
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-1.0, 1.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
    }
}
fn patches() -> [MaterialContact; 3] {
    [(7, 10, cell(0)), (3, 10, cell(1)), (11, 20, cell(0))].map(|(id, chamber, cell)| {
        MaterialContact {
            contact: Contact { id, chamber, cell },
            geometry: GasConductionGeometry {
                area_m2: 0.008,
                gas_distance_m: 0.001,
                cell_distance_m: 0.02,
            },
        }
    })
}
fn pistons() -> [GasPiston; 2] {
    [0.75, 0.25].map(|x| GasPiston::new(gas(x, 0.001), 1.0, 0.01, 0.0).unwrap())
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-6, "{a} != {b}");
}

#[test]
fn no_chambers_or_isolated_chambers_reuse_exact_original_owners() {
    let (a, ca) = material(0.2, 2.0);
    let (b, cb) = material(0.4, 4.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }];
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [400.0; 2]);
    let circuit = circuit(&binding, &thermal, 1.0);
    let original = binding
        .step_candidate(
            &thermal,
            0.00005,
            &circuit,
            &[],
            solve().network,
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let empty = material_step(
        &binding,
        &thermal,
        0.00005,
        &circuit,
        &[],
        &[],
        &[],
        solve(),
        &[0.0; 2],
        || false,
    )
    .unwrap();
    assert_eq!(empty.thermal().cells(), original.thermal().cells());
    assert_eq!(empty.circuit(), original.circuit());
    let gases = pistons();
    let mut chambers = [
        chamber(10, &gases[0], 0.00005),
        chamber(20, &gases[1], 0.00005),
    ];
    chambers[1].input.heat_j = 0.02;
    let isolated = material_step(
        &binding,
        &thermal,
        0.00005,
        &circuit,
        &[],
        &chambers,
        &[],
        solve(),
        &[0.0; 2],
        || false,
    )
    .unwrap();
    assert_eq!(isolated.thermal().cells(), original.thermal().cells());
    assert_eq!(isolated.circuit(), original.circuit());
    for (i, c) in chambers.iter().enumerate() {
        let mut expected = gases[i].clone();
        expected.step(c.input, c.piston_solve, || false).unwrap();
        assert_eq!(isolated.chambers()[i].boundary.piston(), &expected);
    }
    near(isolated.step().external_network_heat_j, 0.0);
}

#[test]
fn source_only_change_reaches_both_gases_bodies_and_finite_electrical_history() {
    let run = |k| {
        let (a, ca) = material(k, 2.0);
        let (b, cb) = material(0.4, 4.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let p = point(400.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &p, &key);
        let resistors = resistors();
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.001,
            distance_b_m: 0.001,
        }];
        let binding = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let mut thermal = initial(&curves, [400.0; 2]);
        let mut circuit = circuit(&binding, &thermal, 1.0);
        let mut gases = pistons();
        let old = gases.clone();
        let energy0 = thermal.energy_j() + 50.0;
        let mut external_work = 0.0;
        for n in 1..=20 {
            let target = f64::from(n) * 0.00005;
            let dt = target - thermal.time_s();
            let chambers = [chamber(10, &gases[0], dt), chamber(20, &gases[1], dt)];
            let candidate = material_step(
                &binding,
                &thermal,
                target,
                &circuit,
                &[],
                &chambers,
                &patches(),
                solve(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
            assert!(
                candidate.step().total_energy_residual_j.abs() <= solve().total_energy_tolerance_j
            );
            near(candidate.step().external_network_heat_j, 0.0);
            for chamber in candidate.chambers() {
                for contact in chamber.boundary.source() {
                    let material = candidate
                        .material()
                        .cells
                        .iter()
                        .find(|s| s.cell == contact.source.0.cell)
                        .unwrap();
                    assert_eq!(material.conductivity.as_ref(), Some(&contact.source.1));
                    near(
                        contact.port.temperature_k,
                        candidate
                            .thermal()
                            .temperature_k(contact.source.0.cell)
                            .unwrap(),
                    );
                }
                external_work += chamber.boundary.step().piston.exterior_pressure_work_j;
            }
            let (t, c, report, _) = candidate.into_parts();
            for (i, chamber) in report.source.chambers.into_iter().enumerate() {
                gases[i] = chamber.boundary.into_parts().0;
            }
            let gas_energy: f64 = gases
                .iter()
                .zip(&old)
                .map(|(g, o)| {
                    g.gas().mass_kg()
                        * (g.gas().state().specific_internal_energy_j_kg()
                            - o.gas().state().specific_internal_energy_j_kg())
                        + 0.5 * g.mass_kg() * g.velocity_m_s().powi(2)
                })
                .sum();
            near(
                t.energy_j() + report.electrical.stored_after_j + gas_energy + external_work,
                energy0,
            );
            thermal = t;
            circuit = c;
        }
        (
            thermal.temperature_k(cell(0)).unwrap(),
            thermal.temperature_k(cell(1)).unwrap(),
            circuit.charges_c()[0],
            gases,
        )
    };
    let a = run(0.2);
    let b = run(2.0);
    assert!((a.0 - b.0).abs() > 1e-5 && (a.1 - b.1).abs() > 1e-5 && (a.2 - b.2).abs() > 1e-10);
    for i in 0..2 {
        assert!((a.3[i].position_m() - b.3[i].position_m()).abs() > 1e-10);
        assert!(
            (a.3[i].gas().state().state().sound_speed - b.3[i].gas().state().state().sound_speed)
                .abs()
                > 1e-6
        );
    }
}

#[test]
fn finite_hot_and_cold_cells_exchange_through_one_gas_at_zero_net_heat() {
    let (a, ca) = material(0.2, 2.0);
    let (b, cb) = material(0.4, 4.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }];
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [500.0, 300.0]);
    let circuit = circuit(&binding, &thermal, 0.0);
    let gases = pistons();
    let destinations: Vec<_> = resistors.iter().map(|r| r.destination).collect();
    let ports = [
        Contact {
            id: 1,
            chamber: 10,
            cell: cell(0),
        },
        Contact {
            id: 2,
            chamber: 10,
            cell: cell(1),
        },
    ];
    let result = super::step_candidate(
        &thermal,
        0.00005,
        CircuitEnthalpyInput {
            circuit: &circuit,
            ports: &[],
            destinations: &destinations,
        },
        &[chamber(10, &gases[0], 0.00005)],
        &ports,
        solve(),
        &[0.0; 2],
        |s| {
            let mut c = binding.resolve(s, || false)?;
            c.conductances_w_k.insert(7, 0.0); // explicit isolation of direct solid/solid path in this reference
            Ok::<_, MaterialNetworkError>(c)
        },
        |_, _, _, _| {
            Ok::<_, MaterialNetworkError>(ConductanceSample {
                conductance_w_k: 1.0,
                source: (),
            })
        },
        || false,
    )
    .unwrap();
    let expected = 100.0 / (1.0 + 0.00005 / 0.5);
    near(
        result.thermal().temperature_k(cell(0)).unwrap(),
        400.0 + expected,
    );
    near(
        result.thermal().temperature_k(cell(1)).unwrap(),
        400.0 - expected,
    );
    let gas = &result.chambers()[0].boundary;
    near(gas.piston().position_m(), gases[0].position_m());
    near(gas.step().exchanged_heat_j, 0.0);
    assert!(gas.source()[0].heat_j > 0.04 && gas.source()[1].heat_j < -0.04);
    near(gas.source()[0].heat_j + gas.source()[1].heat_j, 0.0);
    near(result.thermal().energy_j(), thermal.energy_j());
}

#[test]
fn ordering_replay_and_aggregate_surface_admission_preserve_real_owners() {
    let (a, ca) = material(0.2, 2.0);
    let (b, cb) = material(0.4, 4.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }];
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [400.0; 2]);
    let circuit = circuit(&binding, &thermal, 1.0);
    let gases = pistons();
    let chambers = [
        chamber(10, &gases[0], 0.00005),
        chamber(20, &gases[1], 0.00005),
    ];
    let patches = patches();
    let run = |c: &[Chamber<'_>], p: &[MaterialContact]| {
        material_step(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            c,
            p,
            solve(),
            &[0.0; 2],
            || false,
        )
    };
    let a = run(&chambers, &patches).unwrap();
    let b = run(
        &[chambers[1], chambers[0]],
        &[patches[2], patches[0], patches[1]],
    )
    .unwrap();
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(a.circuit(), b.circuit());
    for (a, b) in a.chambers().iter().zip(b.chambers()) {
        assert_eq!(a.id, b.id);
        assert_eq!(a.boundary, b.boundary);
    }
    let mut over = patches;
    over[0].geometry.area_m2 = 0.024;
    over[2].geometry.area_m2 = 0.024;
    assert!(run(&chambers, &over).is_err()); // each alone fits, their sum does not
    assert!(
        run(
            &[
                chambers[0],
                Chamber {
                    id: 20,
                    ..chambers[0]
                }
            ],
            &patches
        )
        .is_err()
    );
    assert!(run(&chambers, &[patches[0], patches[0]]).is_err());
    let bad = Chamber {
        input: PistonInput {
            duration_s: 0.001,
            ..chambers[1].input
        },
        ..chambers[1]
    };
    assert!(run(&[chambers[0], bad], &patches).is_err());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), &[1.0]);
}

#[test]
fn last_chamber_failure_or_final_work_and_cancellation_publish_nothing() {
    let (a, ca) = material(0.2, 2.0);
    let (b, cb) = material(0.4, 4.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }];
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, [400.0; 2]);
    let circuit = circuit(&binding, &thermal, 1.0);
    let gases = pistons();
    let saved = gases.clone();
    let chambers = [
        chamber(10, &gases[0], 0.00005),
        chamber(20, &gases[1], 0.00005),
    ];
    let mut count = 0;
    let accepted = material_step(
        &binding,
        &thermal,
        0.00005,
        &circuit,
        &[],
        &chambers,
        &patches(),
        solve(),
        &[0.0; 2],
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let calls = accepted.step().boundary_solves;
    let mut actual = 0;
    assert!(matches!(
        material_step(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            &chambers,
            &patches(),
            solve(),
            &[0.0; 2],
            || {
                actual += 1;
                actual == count
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        material_step(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            &chambers,
            &patches(),
            NetworkSolve {
                maximum_boundary_solves: calls - 1,
                ..solve()
            },
            &[0.0; 2],
            || false
        )
        .is_err()
    );
    let destinations: Vec<_> = resistors.iter().map(|r| r.destination).collect();
    let ports: Vec<_> = patches().iter().map(|p| p.contact).collect();
    let rejected = super::step_candidate(
        &thermal,
        0.00005,
        CircuitEnthalpyInput {
            circuit: &circuit,
            ports: &[],
            destinations: &destinations,
        },
        &chambers,
        &ports,
        solve(),
        &[0.0; 2],
        |s| binding.resolve(s, || false),
        |contact, _, _, _| {
            if contact.chamber == 20 {
                Err("second chamber interface unavailable")
            } else {
                Ok(ConductanceSample {
                    conductance_w_k: 0.2,
                    source: (),
                })
            }
        },
        || false,
    );
    assert!(matches!(
        rejected,
        Err(CircuitEnthalpyError::Source(NetworkError::Boundary {
            chamber: 20,
            ..
        }))
    ));
    assert_eq!(gases, saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), &[1.0]);
}

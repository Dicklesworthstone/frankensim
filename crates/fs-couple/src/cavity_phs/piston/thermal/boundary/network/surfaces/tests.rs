//! Real source cards, finite circuit, two walls and ONE shared gas/motion owner.
use super::material::*;
use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
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
fn card(k: f64, r: f64, slope: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, a, b) in [
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            r * 0.0025 * (1.0 - 200.0 * slope),
            r * 0.0025 * (1.0 + 400.0 * slope),
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k, k),
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
                    source: "authored multi-wall fixture".into(),
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
            process: "authored".into(),
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
    linked: bool,
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
    if linked {
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
    }
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
fn solve() -> NetworkGasSolve {
    NetworkGasSolve {
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
        maximum_boundary_solves: 300,
        total_energy_tolerance_j: 1e-6,
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn surfaces() -> [MaterialGasSurface; 2] {
    [0, 1].map(|i| MaterialGasSurface {
        port: NetworkThermalPort {
            id: i as u64 + 100,
            cell: cell(i),
        },
        geometry: GasConductionGeometry {
            area_m2: 0.01,
            gas_distance_m: 0.001,
            cell_distance_m: 0.02,
        },
    })
}
fn gas_drive(surfaces: &[MaterialGasSurface], dt: f64) -> MaterialSurfacesDrive<'_> {
    MaterialSurfacesDrive {
        surfaces,
        input: PistonInput {
            duration_s: dt,
            ..drive(0.0)
        },
        piston: controls(),
        exchange: exchange(),
    }
}
fn piston() -> GasPiston {
    GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-6, "{a} != {b}");
}

#[test]
fn balanced_finite_walls_exchange_energy_without_heating_or_moving_the_gas() {
    let (a, ca) = card(0.2, 2.0, 0.0);
    let (b, cb) = card(0.2, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curves, [500.0, 300.0], false);
    let circuit = circuit(&binding, &thermal, 0.0);
    let old = piston();
    let surfaces = surfaces();
    let candidate = old
        .material_network_surfaces_candidate(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(&surfaces, 0.00005),
            solve(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let (next, electrical, gas, report) = candidate.into_parts();
    let g = report.ports[0].conductance_w_k;
    // Closed-form backward Euler on two equal 0.5 J/K stores with a 400 K gas.
    let heat = 0.00005 * g * 100.0 / (1.0 + 0.00005 * g / 0.5);
    assert!(heat > 0.0);
    near(report.ports[0].heat_j, heat);
    near(report.ports[1].heat_j, -heat);
    near(next.temperature_k(cell(0)).unwrap(), 500.0 - heat / 0.5);
    near(next.temperature_k(cell(1)).unwrap(), 300.0 + heat / 0.5);
    near(gas.gas().state().state().temperature, 400.0);
    assert!((gas.position_m() - old.position_m()).abs() < 1e-10);
    near(next.energy_j(), thermal.energy_j());
    assert_eq!(electrical.charges_c(), &[0.0]);
    near(report.coupled.external_network_heat_j, 0.0);
}

#[test]
fn conductivity_only_swap_reaches_both_walls_circuit_and_shared_gas_motion() {
    let run = |k| {
        let (a, ca) = card(k, 2.0, 0.002);
        let (b, cb) = card(0.4, 4.0, 0.001);
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
        let mut thermal = initial(&curves, [400.0; 2], true);
        let mut circuit = circuit(&binding, &thermal, 1.0);
        let mut gas = piston();
        let old = gas.clone();
        let initial_heat = thermal.energy_j();
        let mut exterior = 0.0;
        let surfaces = surfaces();
        for i in 1..=40 {
            let target = f64::from(i) * 0.00005;
            let dt = target - thermal.time_s();
            let candidate = gas
                .material_network_surfaces_candidate(
                    &binding,
                    &thermal,
                    target,
                    &circuit,
                    &[],
                    gas_drive(&surfaces, dt),
                    solve(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let (t, c, g, report) = candidate.into_parts();
            assert_eq!(report.ports.len(), 2);
            for (port, destination) in report.ports.iter().zip(&report.destinations) {
                let sample = report
                    .network
                    .source
                    .cells
                    .iter()
                    .find(|s| s.cell == destination.cell)
                    .unwrap();
                assert_eq!(sample.conductivity.as_ref(), Some(&port.source));
                near(
                    port.port.temperature_k,
                    t.temperature_k(destination.cell).unwrap(),
                );
            }
            near(report.coupled.external_network_heat_j, 0.0);
            exterior += report.boundary.piston.exterior_pressure_work_j;
            near(
                t.energy_j() - initial_heat + report.network.electrical.stored_after_j - 50.0
                    + g.gas().mass_kg()
                        * (g.gas().state().specific_internal_energy_j_kg()
                            - old.gas().state().specific_internal_energy_j_kg())
                    + 0.5 * g.mass_kg() * g.velocity_m_s().powi(2)
                    + exterior,
                0.0,
            );
            thermal = t;
            circuit = c;
            gas = g;
        }
        (
            thermal.temperature_k(cell(0)).unwrap(),
            thermal.temperature_k(cell(1)).unwrap(),
            circuit.charges_c()[0],
            gas,
        )
    };
    let a = run(0.2);
    let b = run(2.0);
    assert!((a.0 - b.0).abs() > 1e-4 && (a.1 - b.1).abs() > 1e-4);
    assert!((a.2 - b.2).abs() > 1e-9);
    assert!((a.3.position_m() - b.3.position_m()).abs() > 1e-9);
    assert!(
        (a.3.gas().state().state().sound_speed - b.3.gas().state().state().sound_speed).abs()
            > 1e-6
    );
}

#[test]
fn splitting_one_cell_surface_retains_motion_and_checks_total_coverage() {
    let (a, ca) = card(0.2, 2.0, 0.002);
    let (b, cb) = card(0.4, 4.0, 0.001);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = initial(&curves, [500.0, 400.0], false);
    let circuit = circuit(&binding, &thermal, 0.0);
    let gas = piston();
    let original = surfaces()[0];
    let run = |surfaces: &[MaterialGasSurface]| {
        gas.material_network_surfaces_candidate(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(surfaces, 0.00005),
            solve(),
            &[0.0; 2],
            || false,
        )
    };
    let one = run(&[original]).unwrap();
    let mut split = [original, original];
    split[1].port.id += 10;
    for s in &mut split {
        s.geometry.area_m2 *= 0.5;
    }
    let two = run(&split).unwrap();
    near(
        one.boundary().piston().position_m(),
        two.boundary().piston().position_m(),
    );
    near(
        one.network().thermal().energy_j(),
        two.network().thermal().energy_j(),
    );
    near(
        two.boundary().source()[0].heat_j,
        two.boundary().source()[1].heat_j,
    );
    let reverse = run(&[split[1], split[0]]).unwrap();
    assert_eq!(two.boundary(), reverse.boundary());
    assert_eq!(
        two.network().thermal().cells(),
        reverse.network().thermal().cells()
    );
    for s in &mut split {
        s.geometry.area_m2 = 0.03;
    }
    assert!(
        run(&split).is_err(),
        "two individually admissible patches exceed total exposed area"
    );
}

#[test]
fn late_failure_and_missing_cold_surface_data_return_no_partial_owners() {
    let (a, ca) = card(0.2, 2.0, 0.002);
    let (b, cb) = card(0.4, 4.0, 0.001);
    let cards = [a, b];
    let curves = [ca, cb];
    let p = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut cells = cells(&cards, &curves, &p, &key);
    let resistors = resistors();
    let thermal = initial(&curves, [500.0, 300.0], false);
    let gas = piston();
    let saved = gas.clone();
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let circuit = circuit(&binding, &thermal, 0.0);
    let old_cells = thermal.cells().clone();
    let old_circuit = circuit.clone();
    let surfaces = surfaces();
    let mut polls = 0;
    let accepted = gas
        .material_network_surfaces_candidate(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(&surfaces, 0.00005),
            solve(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut actual = 0;
    assert!(matches!(
        gas.material_network_surfaces_candidate(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(&surfaces, 0.00005),
            solve(),
            &[0.0; 2],
            || {
                actual += 1;
                actual == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        gas.material_network_surfaces_candidate(
            &binding,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(&surfaces, 0.00005),
            NetworkGasSolve {
                maximum_boundary_solves: accepted.step().boundary_solves - 1,
                ..solve()
            },
            &[0.0; 2],
            || false
        )
        .is_err()
    );
    cells[1].conductivity = None;
    let incomplete = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let mut isolated = surfaces;
    for s in &mut isolated {
        s.geometry.area_m2 = 0.0;
    }
    assert!(
        gas.material_network_surfaces_candidate(
            &incomplete,
            &thermal,
            0.00005,
            &circuit,
            &[],
            gas_drive(&isolated, 0.00005),
            solve(),
            &[0.0; 2],
            || false
        )
        .is_err()
    );
    assert_eq!(gas, saved);
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(circuit, old_circuit);
}

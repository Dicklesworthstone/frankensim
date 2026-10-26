//! Independent two-body finite-capacitor/current/thermal closure.
use super::*;
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialNetwork};
use fs_em::winding::nonlinear::{
    NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindings,
    resistors::GraphResistorMaterial,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{QuantitySpec, Temperature};

const RESISTOR_CELL: CellId = CellId(8);
const RESISTOR_MATERIAL: MaterialId = MaterialId(10);
const CAPACITANCE: f64 = 1.0;
const CHARGE: f64 = 8.0;
const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 17,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];
const DESTINATIONS: [ResistorCell; 1] = [ResistorCell {
    branch: 2,
    cell: RESISTOR_CELL,
}];

fn graph_card(alpha: f64, upper_t: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let basis = ElectricalPropertyBasis::Resistivity;
    for (name, dims, value) in [
        (
            basis.property_name(),
            basis.dimensions(),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, 0.5),
                    (300.0, 0.5),
                    (upper_t, 0.5 * (1.0 + alpha * (upper_t - 300.0))),
                ],
                dims: basis.dimensions(),
            },
        ),
        (
            "density",
            Density::DIMS,
            PropertyValue::Scalar {
                value: 2.0,
                dims: Density::DIMS,
            },
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            PropertyValue::Scalar {
                value: 3.0,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            PropertyValue::Scalar {
                value: 1.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
        ),
    ] {
        let interpolation = if matches!(value, PropertyValue::Curve { .. }) {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    if name == basis.property_name() {
                        upper_t
                    } else {
                        500.0
                    },
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured ordinary resistor body".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic resistor".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn two_bodies<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'a> {
    let mut events = [CELL, RESISTOR_CELL]
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: if i == 0 { MATERIAL } else { RESISTOR_MATERIAL },
                mass_kg: 2.0 * VOLUME,
                specific_h_j_kg: 300.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: VOLUME,
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
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 17,
            a: CELL,
            b: RESISTOR_CELL,
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MATERIAL, curves[0].curve() as &dyn ThermalLaw),
            (RESISTOR_MATERIAL, curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
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
fn bindings<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: if i == 0 { CELL } else { RESISTOR_CELL },
        material_id: if i == 0 { MATERIAL } else { RESISTOR_MATERIAL },
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn resistors() -> [CellResistor; 1] {
    [CellResistor {
        destination: DESTINATIONS[0],
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: VOLUME,
    }]
}
fn native<'a>(
    problem: &'a fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
) -> NonlinearWindingCircuit<'a> {
    let profiles = field::profiles(1);
    NonlinearWindingCircuit::new(
        field::initial_field(problem, &profiles, &[0.0]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: CAPACITANCE,
                    },
                ),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: 1.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        field::conductors(card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![CHARGE],
            external_inputs: vec![],
        },
        field::config(),
        &mut || false,
    )
    .unwrap()
}
fn solve_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        ..controls()
    }
}

// The P1 fixture has lambda(I)=2*depth*B(I), with the independently authored
// piecewise H(B) inverse below. For a trial current the TWO material heat
// equations are a linear 2x2 system, including the exact shared G=1 W/K.
// Only a scalar bisection remains; no production constitutive/solver is called.
fn oracle(q0: f64, i0: f64, t0: [f64; 2], alpha: f64) -> (f64, [f64; 2]) {
    let temperature = |i: f64| {
        let heat = DT * i * i;
        let a = CAPACITY + DT - heat * 0.005;
        let d = CAPACITY + DT - heat * alpha;
        let b = -DT;
        let rhs = [
            CAPACITY * (t0[0] - 300.0) + heat,
            CAPACITY * (t0[1] - 300.0) + heat,
        ];
        let determinant = a * d - b * b;
        assert!(determinant > 0.0);
        [
            300.0 + (d * rhs[0] - b * rhs[1]) / determinant,
            300.0 + (a * rhs[1] - b * rhs[0]) / determinant,
        ]
    };
    let residual = |i: f64| {
        let t = temperature(i);
        2.0 * DEPTH * (flux_density(i, true) - flux_density(i0, true)) / DT
            + (winding_resistance(t[0]) + 1.0 + alpha * (t[1] - 300.0)) * i
            + DT * i / (2.0 * CAPACITANCE)
            - q0 / CAPACITANCE
    };
    let (mut lo, mut hi) = (0.0, 8.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let current = f64::midpoint(lo, hi);
    (current, temperature(current))
}

#[test]
fn finite_capacitor_and_nonlinear_field_use_ordinary_resistor_temperature_and_its_own_heat_cell() {
    let mut outcomes = Vec::new();
    for alpha in [0.0, 0.08] {
        let cards = [field::source(true, 18.0), graph_card(alpha, 500.0)];
        let curves = cards.each_ref().map(caloric);
        let mut thermal = two_bodies(&curves);
        let problem = field::problem(&cards[0], DEPTH);
        let mut circuit = native(&problem, &cards[0]);
        let context = field::material::point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = routes();
        let resistors = resistors();
        let layout = WindingThermalNetwork {
            cells: &cells,
            links: &LINKS,
            magnetic_regions: &regions,
            graph_resistors: &DESTINATIONS,
            volume_tolerance_m3: 1e-12,
        };
        let bases = [(2, ElectricalPropertyBasis::Resistivity)];
        let original = thermal.energy_j() + CHARGE * CHARGE / (2.0 * CAPACITANCE);
        let mut numerical = 0.0;
        for step in 1..=2 {
            let q0 = circuit.state()[circuit.dae().charge_index[0]];
            let i0 = circuit.currents_a()[0];
            let t0 = [
                thermal.temperature_k(CELL).unwrap(),
                thermal.temperature_k(RESISTOR_CELL).unwrap(),
            ];
            let (expected, t) = oracle(q0, i0, t0, alpha);
            let candidate = layout
                .step_nonlinear_material_candidate(
                    &thermal,
                    step as f64 * DT,
                    &circuit,
                    &[],
                    &resistors,
                    &bases,
                    solve_controls(),
                    field::controls(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(candidate.circuit().currents_a()[0], expected);
            near(
                candidate.circuit().state()[candidate.circuit().dae().charge_index[0]],
                q0 - DT * expected,
            );
            near(candidate.thermal().temperature_k(CELL).unwrap(), t[0]);
            near(
                candidate.thermal().temperature_k(RESISTOR_CELL).unwrap(),
                t[1],
            );
            let report = candidate.report();
            let electrical = &report.electrical;
            assert_eq!(electrical.supplied_work_j, 0.0);
            assert_eq!(electrical.graph_materials.len(), 1);
            let source = &electrical.graph_materials[0];
            assert_eq!(source.branch, 2);
            assert_eq!(
                source.conductor.material().card_identity(),
                cards[1].content_hash()
            );
            near(
                source.conductor.resistance_ohm(),
                1.0 + alpha * (t[1] - 300.0),
            );
            near(
                electrical.graph_resistors[0].heat_j,
                DT * expected * expected * source.conductor.resistance_ohm(),
            );
            near(
                CAPACITY * (t[0] - t0[0]),
                DT * winding_resistance(t[0]) * expected * expected - DT * (t[0] - t[1]),
            );
            near(
                CAPACITY * (t[1] - t0[1]),
                electrical.graph_resistors[0].heat_j + DT * (t[0] - t[1]),
            );
            near(report.source.network.links[0].conductance_w_k, 1.0);
            near(report.deposited_joule_j, electrical.heat_j);
            near(report.thermal.stored_energy_change_j, electrical.heat_j);
            assert!(electrical.numerical_dissipation_j > 0.0);
            numerical += electrical.numerical_dissipation_j;
            near(
                candidate.thermal().energy_j() + electrical.stored_after_j + numerical,
                original,
            );
            assert_eq!(circuit.time_s(), (step - 1) as f64 * DT);
            let (next_thermal, next_circuit, _) = candidate.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
        }
        // The legacy thermal entry cannot silently freeze an attached source.
        assert!(
            layout
                .step_nonlinear_candidate(
                    &thermal,
                    3.0 * DT,
                    &circuit,
                    &[],
                    solve_controls(),
                    field::controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        outcomes.push((
            circuit.currents_a()[0],
            thermal.temperature_k(RESISTOR_CELL).unwrap(),
        ));
    }
    assert!(
        outcomes[1].0 < outcomes[0].0 - 1e-3,
        "ordinary R(T) must react back on the finite circuit and field"
    );
    assert!((outcomes[1].1 - outcomes[0].1).abs() > 1e-3);
}

#[test]
fn resistor_refresh_and_thermal_refusals_preserve_every_accepted_history() {
    for upper in [500.0, 300.01] {
        let cards = [field::source(true, 18.0), graph_card(0.08, upper)];
        let curves = cards.each_ref().map(caloric);
        let thermal = two_bodies(&curves);
        let energy = thermal.energy_j();
        let problem = field::problem(&cards[0], DEPTH);
        let original = native(&problem, &cards[0]);
        let context = field::material::point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = routes();
        let resistors = resistors();
        let material_network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &LINKS,
        };
        let initial = material_network.resolve(&thermal, || false).unwrap();
        let sources = initial
            .resistors
            .iter()
            .map(|s| GraphResistorMaterial {
                branch: s.branch,
                conductor: s.conductor.clone(),
            })
            .collect::<Vec<_>>();
        let circuit = original
            .with_graph_resistors(sources.clone(), &mut || false)
            .unwrap();
        assert_eq!(circuit.state(), original.state());
        assert_eq!(
            circuit.field().potential_wb_m(),
            original.field().potential_wb_m()
        );
        assert_eq!(circuit.flux_linkages_wb(), original.flux_linkages_wb());
        assert_eq!(circuit.currents_a(), original.currents_a());
        assert!(circuit.with_graph_resistors(vec![], &mut || false).is_err());
        let mut wrong = sources.clone();
        wrong[0].branch = 1;
        assert!(circuit.with_graph_resistors(wrong, &mut || false).is_err());
        assert!(circuit.with_graph_resistors(sources, &mut || true).is_err());
        let layout = WindingThermalNetwork {
            cells: &cells,
            links: &LINKS,
            magnetic_regions: &regions,
            graph_resistors: &DESTINATIONS,
            volume_tolerance_m3: 1e-12,
        };
        let bases = [(2, ElectricalPropertyBasis::Resistivity)];
        if upper < 400.0 {
            assert!(
                layout
                    .step_nonlinear_material_candidate(
                        &thermal,
                        DT,
                        &circuit,
                        &[],
                        &resistors,
                        &bases,
                        solve_controls(),
                        field::controls(),
                        &[0.0; 2],
                        || false
                    )
                    .is_err()
            );
        } else {
            let mut calls = 0;
            layout
                .step_nonlinear_material_candidate(
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    &resistors,
                    &bases,
                    solve_controls(),
                    field::controls(),
                    &[0.0; 2],
                    || {
                        calls += 1;
                        false
                    },
                )
                .unwrap();
            let mut actual = 0;
            assert!(matches!(
                layout.step_nonlinear_material_candidate(
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    &resistors,
                    &bases,
                    solve_controls(),
                    field::controls(),
                    &[0.0; 2],
                    || {
                        actual += 1;
                        actual == calls
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert!(
                layout
                    .step_nonlinear_material_candidate(
                        &thermal,
                        DT,
                        &circuit,
                        &[],
                        &[],
                        &[],
                        solve_controls(),
                        field::controls(),
                        &[0.0; 2],
                        || false
                    )
                    .is_err()
            );
            let mut limited = solve_controls();
            limited.maximum_circuit_trials = 1;
            assert!(
                layout
                    .step_nonlinear_material_candidate(
                        &thermal,
                        DT,
                        &circuit,
                        &[],
                        &resistors,
                        &bases,
                        limited,
                        field::controls(),
                        &[0.0; 2],
                        || false
                    )
                    .is_err()
            );
        }
        assert_eq!(circuit.state(), original.state());
        assert_eq!(
            circuit.field().potential_wb_m(),
            original.field().potential_wb_m()
        );
        assert_eq!(circuit.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), energy);
        assert_eq!(thermal.temperature_k(RESISTOR_CELL).unwrap(), 300.0);
    }
}

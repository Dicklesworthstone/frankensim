//! G1/G3: source-backed nonlinear conduction and finite heat on the real owners.
//! The two authored one-axis laws below are synthetic numerical fixtures.
#[path = "circuit_enthalpy/common.rs"]
mod common;

use common::{DESTINATIONS, config, initial, near, point};
use fs_couple::circuit_enthalpy::conductivity::step_candidate;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitThermalCoefficients};
use fs_couple::material_circuit::conductivity::{
    ConductionAxis, ConductionCoordinate, ConductivitySolve, ConductorResponse,
};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitError, MaterialResistor};
use fs_matdb::*;
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_PROPERTY, ResolvedConductorState, resolve_conductor_state,
};
use fs_material::phase::{
    EquilibriumPhaseState, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::collections::BTreeMap;

const FIELD: Dims = Dims([1, 1, -3, 0, -1, 0]);
const AXES: [ConductionAxis<'static>; 1] = [ConductionAxis {
    name: "E",
    coordinate: ConductionCoordinate::ElectricFieldMagnitude,
}];

fn query(t: f64) -> QueryPoint {
    point(t)
        .with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0)
        .unwrap()
}

fn conductivity_controls() -> ConductivitySolve {
    ConductivitySolve {
        maximum_iterations: 64,
        absolute_resistance_tolerance_ohm: 1e-12,
        relative_resistance_tolerance: 1e-12,
    }
}

fn field_material(
    slope: f64,
    maximum_temperature: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = common::material(1.0, 1.0, 0.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == ELECTRICAL_RESISTIVITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "E".into(),
                abscissa_dims: FIELD,
                knots: vec![(0.0, 1.0), (20.0, 1.0 + 20.0 * slope)],
                dims: claim.key.quantity().dims(),
            };
            claim.validity = fs_evidence::ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    maximum_temperature,
                )
                .with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0, 20.0)
                .with("lot", 7.0, 7.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &query(200.0),
            upper: &query(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn circuit(cards: &[MaterialCard; 2]) -> MaterialCircuit {
    let resistors = cards
        .iter()
        .enumerate()
        .map(|(i, card)| MaterialResistor {
            branch: DESTINATIONS[i].branch,
            conductor: resolve_conductor_state(
                card,
                &if i == 0 { query(300.0) } else { point(300.0) },
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
            .into_uniform_conductor(1.0, 1.0)
            .unwrap(),
        })
        .collect();
    MaterialCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Resistor { ohms: 999.0 }),
            ],
            transformers: vec![],
        },
        resistors,
        vec![],
        vec![10.0],
        &[],
    )
    .unwrap()
}

fn responses(cards: &[MaterialCard; 2]) -> Vec<ConductorResponse<'_>> {
    cards
        .iter()
        .enumerate()
        .map(|(i, card)| ConductorResponse {
            branch: DESTINATIONS[i].branch,
            card,
            axes: if i == 0 { &AXES } else { &[] },
        })
        .collect()
}

fn coefficients(
    state: &DepositionSimulation<'_>,
    original: &MaterialCircuit,
    cards: &[MaterialCard; 2],
    caloric: &[ResolvedHeatCapacityCurve; 2],
) -> Result<CircuitThermalCoefficients<Vec<EquilibriumPhaseState>>, String> {
    let mut resistors = Vec::new();
    let mut phases = Vec::new();
    for (i, old) in original.resistors().iter().enumerate() {
        let phase = caloric[i]
            .curve()
            .state_at_specific_enthalpy(state.cells()[&DESTINATIONS[i].cell].specific_h_j_kg)
            .map_err(|e| e.to_string())?;
        let context = old.conductor.material().query();
        let physical = ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        }
        .at(phase)
        .map_err(|e| e.to_string())?;
        let bundle = old
            .conductor
            .material()
            .at_query_point(&cards[i], &physical)
            .map_err(|e| e.to_string())?;
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            bundle,
            old.conductor.source_basis(),
        )
        .map_err(|e| e.to_string())?
        .into_uniform_conductor(old.conductor.length_m(), old.conductor.area_m2())
        .map_err(|e| e.to_string())?;
        resistors.push(MaterialResistor {
            branch: old.branch,
            conductor,
        });
        phases.push(phase);
    }
    Ok(CircuitThermalCoefficients {
        resistors,
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        source: phases,
    })
}

// Independent series-RC and two-cell BE equations, analytically reduced to I.
// R_a=1+beta*E_a, E_a=I*R_a, R_b=3*(1+alpha*(T_b-300)).
fn reference(beta: f64, alpha: f64, dt: f64) -> (f64, f64, f64) {
    let values = |current: f64| {
        let ra = 1.0 / (1.0 - beta * current);
        let determinant = (1.0 + dt) * (2.0 + dt) - dt * dt;
        let tb = dt * current * current * (dt * ra + (1.0 + dt) * 3.0)
            / (determinant - dt * current * current * (1.0 + dt) * 3.0 * alpha);
        let rb = 3.0 * (1.0 + alpha * tb);
        let ta = dt * current * current * ((2.0 + dt) * ra + dt * rb) / determinant;
        (ra, rb, ta, tb)
    };
    let (mut low, mut high) = (0.0, 3.0);
    for _ in 0..80 {
        let current = 0.5 * (low + high);
        let (ra, rb, _, _) = values(current);
        if current * (ra + rb + 0.5 * dt) > 10.0 {
            high = current;
        } else {
            low = current;
        }
    }
    let current = 0.5 * (low + high);
    let (_, _, ta, tb) = values(current);
    (10.0 - dt * current, 300.0 + ta, 300.0 + tb)
}

#[test]
fn g1_field_and_temperature_feedback_close_the_same_finite_circuit_and_heat_stores() {
    let (a, ca) = field_material(0.05, 2000.0);
    let (b, cb) = common::material(3.0, 2.0, 0.005);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 1.0, None);
    let electrical = circuit(&cards);
    let expected = reference(0.05, 0.005, 0.2);
    let run = || {
        step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            config(),
            conductivity_controls(),
            &[0.0; 2],
            |s| coefficients(s, &electrical, &cards, &curves),
            || false,
        )
        .unwrap()
    };
    let result = run();
    let replay = run();
    assert_eq!(result.circuit(), replay.circuit());
    assert_eq!(result.thermal().cells(), replay.thermal().cells());
    near(result.circuit().charges_c()[0], expected.0);
    near(
        result.thermal().temperature_k(CellId(30)).unwrap(),
        expected.1,
    );
    near(
        result.thermal().temperature_k(CellId(10)).unwrap(),
        expected.2,
    );
    let report = result.report();
    near(
        report.electrical.electrical.stored_after_j + result.thermal().energy_j(),
        50.0 + thermal.energy_j(),
    );
    near(report.non_electrical_heat_j, 0.0);
    near(
        report.deposited_joule_j,
        report.electrical.electrical.heat_j,
    );
    assert!(report.circuit_trials >= report.electrical.iterations);
    assert!(report.electrical.iterations > 1);
    for (i, actual) in report.electrical.evaluated.iter().enumerate() {
        near(
            actual.conductor.material().state_coordinate("T").unwrap(),
            result
                .thermal()
                .temperature_k(DESTINATIONS[i].cell)
                .unwrap(),
        );
        if i == 0 {
            near(
                actual.conductor.material().state_coordinate("E").unwrap(),
                actual.point.electric_field_v_m.abs(),
            );
        } else {
            // A thermal-only member needs no fictitious electrical source axis.
            assert!(actual.conductor.material().state_coordinate("E").is_none());
        }
        near(
            report.source[i].temperature_k(),
            result
                .thermal()
                .temperature_k(DESTINATIONS[i].cell)
                .unwrap(),
        );
        assert!(actual.power_defect_w.abs() <= config().circuit.power_tolerance_w);
    }
    assert_eq!(electrical.charges_c(), &[10.0]);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g3_material_curve_swaps_change_both_electrical_discharge_and_temperature() {
    let run = |beta: f64, alpha: f64| {
        let (a, ca) = field_material(beta, 2000.0);
        let (b, cb) = common::material(3.0, 2.0, alpha);
        let cards = [a, b];
        let curves = [ca, cb];
        let thermal = initial(&curves, 1.0, None);
        let electrical = circuit(&cards);
        let result = step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            config(),
            conductivity_controls(),
            &[0.0; 2],
            |s| coefficients(s, &electrical, &cards, &curves),
            || false,
        )
        .unwrap();
        let expected = reference(beta, alpha, 0.2);
        let actual = (
            result.circuit().charges_c()[0],
            result.thermal().temperature_k(CellId(30)).unwrap(),
            result.thermal().temperature_k(CellId(10)).unwrap(),
        );
        near(actual.0, expected.0);
        near(actual.1, expected.1);
        near(actual.2, expected.2);
        actual
    };
    let baseline = run(0.0, 0.0);
    let field = run(0.1, 0.0);
    let both = run(0.1, 0.008);
    assert!(field.0 > baseline.0);
    assert!(both.0 > field.0);
    assert_ne!(field.1, baseline.1);
    assert_ne!(both.2, field.2);
}

#[test]
fn g4_coupled_source_budget_and_cancellation_refusals_preserve_both_histories() {
    let (a, ca) = field_material(0.05, 300.1);
    let (b, cb) = common::material(3.0, 2.0, 0.005);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 1.0, None);
    let electrical = circuit(&cards);
    let original = electrical.clone();
    let cells = thermal.cells().clone();
    let proposed = step_candidate(
        &thermal,
        0.2,
        common::input(&electrical),
        &responses(&cards),
        config(),
        conductivity_controls(),
        &[0.0; 2],
        |s| coefficients(s, &electrical, &cards, &curves),
        || false,
    );
    assert!(
        proposed.is_err(),
        "heated endpoint exceeds the original field-law temperature support"
    );
    let mut limited = config();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            limited,
            conductivity_controls(),
            &[0.0; 2],
            |s| coefficients(s, &electrical, &cards, &curves),
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    // The outer budget must not silently make an invalid authored policy valid.
    let mut invalid = conductivity_controls();
    invalid.maximum_iterations = 129;
    assert!(matches!(
        step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            limited,
            invalid,
            &[0.0; 2],
            |s| coefficients(s, &electrical, &cards, &curves),
            || false
        ),
        Err(CircuitEnthalpyError::Circuit(
            MaterialCircuitError::Invalid(_)
        ))
    ));
    let mut polls = 0;
    assert!(matches!(
        step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            config(),
            conductivity_controls(),
            &[0.0; 2],
            |s| coefficients(s, &electrical, &cards, &curves),
            || {
                polls += 1;
                polls >= 12
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(matches!(
        step_candidate(
            &thermal,
            0.2,
            common::input(&electrical),
            &responses(&cards),
            config(),
            conductivity_controls(),
            &[0.0; 2],
            |_| Err::<CircuitThermalCoefficients<()>, _>("source unavailable"),
            || false
        ),
        Err(CircuitEnthalpyError::Source("source unavailable"))
    ));
    assert_eq!(electrical, original);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
}

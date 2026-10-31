//! One native circuit/enthalpy trajectory with E-dependent bulk and T-dependent
//! interface sources. Independent equations check the distinct physical losses.
use super::*;
use fs_couple::circuit_enthalpy::contacts::{ContactHeatCells, ContactThermalCoefficients};
use fs_couple::circuit_enthalpy::conductivity::step_with_contacts_candidate;
use fs_couple::material_circuit::contacts::SeriesElectricalContact;

#[path = "../electrical_contacts/source.rs"]
mod source;

const PARTITION: f64 = 0.25;
const ROUTES: [ContactHeatCells; 1] = [ContactHeatCells {
    id: 7,
    cells: [CellId(30), CellId(10)],
    heat_fraction_to_first: PARTITION,
}];

fn contact_coefficients(
    state: &DepositionSimulation<'_>,
    original: &MaterialCircuit,
    cards: &[MaterialCard; 2],
    caloric: &[ResolvedHeatCapacityCurve; 2],
    interface: &InterfaceSystemCard,
) -> Result<ContactThermalCoefficients<Vec<EquilibriumPhaseState>>, String> {
    let bulk = coefficients(state, original, cards, caloric)?;
    let old = &original.contacts()[0];
    let point = old.source().resolved().query()
        .with("Ta", bulk.source[0].temperature_k()).map_err(|e| e.to_string())?
        .with("Tb", bulk.source[1].temperature_k()).map_err(|e| e.to_string())?;
    let contact = SeriesElectricalContact::new(old.id(), old.branch(),
        old.source().at_query_point(interface, &point).map_err(|e| e.to_string())?,
        old.area_m2()).map_err(|e| e.to_string())?;
    Ok(ContactThermalCoefficients { bulk, contacts: vec![contact] })
}

// Independent reduction of two backward-Euler heat equations and midpoint RC.
// R_a=1+beta E_a, R_c=scale(1+gamma(T_a-300)), C=1, Cp=(1,2), R_b=3.
fn reference(beta: f64, scale: f64, gamma: f64, conductance: f64, dt: f64)
    -> (f64, f64, f64)
{
    let q = dt * conductance;
    let determinant = (1.0 + q) * (2.0 + q) - q * q;
    let values = |current: f64| {
        let ra = 1.0 / (1.0 - beta * current);
        let heat = dt * current * current;
        let base_ta = heat * ((2.0 + q) * ra + q * 3.0) / determinant;
        let contact_ta = heat * ((2.0 + q) * PARTITION + q * (1.0 - PARTITION)) / determinant;
        let rc = scale * (1.0 + gamma * base_ta) / (1.0 - scale * gamma * contact_ta);
        let ta = base_ta + contact_ta * rc;
        let tb = (heat * (3.0 + (1.0 - PARTITION) * rc) + q * ta) / (2.0 + q);
        (ra, rc, ta, tb)
    };
    let (mut low, mut high) = (0.0, 2.5);
    for _ in 0..80 {
        let current = 0.5 * (low + high);
        let (ra, rc, _, _) = values(current);
        if current * (ra + 3.0 + rc + 0.5 * dt) > 10.0 { high = current; }
        else { low = current; }
    }
    let current = 0.5 * (low + high);
    let (_, _, ta, tb) = values(current);
    (10.0 - dt * current, 300.0 + ta, 300.0 + tb)
}

#[test]
fn g1_nonlinear_bulk_and_contact_heat_share_the_actual_thermal_endpoint() {
    let beta = 0.05;
    let gamma = 0.005;
    let dt = 0.2;
    let mut charges = Vec::new();
    for scale in [0.5, 2.0] {
        let (a, ca) = field_material(beta, 2000.0);
        let (b, cb) = common::material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let mut thermal = initial(&curves, 1.0, None);
        let interface = source::card(&cards[0], &cards[1], scale, gamma);
        let mut electrical = circuit(&cards).rebind_contacts_candidate(
            vec![source::contact(7, 1, &interface, 1.0)], &[],
        ).unwrap();
        let total0 = 50.0 + thermal.energy_j();
        for n in 1..=3 {
            let result = step_with_contacts_candidate(
                &thermal, f64::from(n) * dt, common::input(&electrical),
                &responses(&cards), &ROUTES, config(), conductivity_controls(), &[0.0; 2],
                |state| contact_coefficients(state, &electrical, &cards, &curves, &interface),
                || false,
            ).unwrap();
            let report = result.report();
            let native = &report.electrical.electrical;
            let ta = result.thermal().temperature_k(CellId(30)).unwrap();
            let tb = result.thermal().temperature_k(CellId(10)).unwrap();
            let surface = &native.contacts[0];
            let bulk = &native.resistors[0];
            near(surface.contact.resistance_ohm(), scale * (1.0 + gamma * (ta - 300.0)));
            near(surface.contact.source().resolved().state_coordinate("Ta").unwrap(), ta);
            near(surface.contact.source().resolved().state_coordinate("Tb").unwrap(), tb);
            near(report.electrical.evaluated[0].point.electric_field_v_m, bulk.voltage_v);
            near(report.electrical.evaluated[0].conductor.resistance_ohm(), 1.0 + beta * bulk.voltage_v);
            let ha = result.thermal().cells()[&CellId(30)].specific_h_j_kg
                - thermal.cells()[&CellId(30)].specific_h_j_kg;
            let hb = result.thermal().cells()[&CellId(10)].specific_h_j_kg
                - thermal.cells()[&CellId(10)].specific_h_j_kg;
            let exchange = dt * (tb - ta);
            near(ha, bulk.heat_j + PARTITION * surface.heat_j + exchange);
            near(hb, native.resistors[1].heat_j + (1.0 - PARTITION) * surface.heat_j - exchange);
            near(native.stored_after_j + result.thermal().energy_j(), total0);
            near(report.deposited_joule_j, native.heat_j);
            near(report.non_electrical_heat_j, 0.0);
            assert_eq!(report.source.destinations, ROUTES);
            assert!(report.electrical.iterations > 1);
            if n == 1 {
                let expected = reference(beta, scale, gamma, 1.0, dt);
                near(result.circuit().charges_c()[0], expected.0);
                near(ta, expected.1);
                near(tb, expected.2);
                charges.push(result.circuit().charges_c()[0]);
            }
            let (next_thermal, next_electrical, _) = result.into_parts();
            thermal = next_thermal;
            electrical = next_electrical;
        }
    }
    assert!(charges[1] > charges[0], "source-only contact change alters actual discharge");
}

#[test]
fn g4_contact_source_budget_route_and_cancellation_failures_preserve_both_histories() {
    let (a, ca) = field_material(0.05, 2000.0);
    let (b, cb) = common::material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 1.0, None);
    let interface = source::card(&cards[0], &cards[1], 2.0, 0.005);
    let electrical = circuit(&cards).rebind_contacts_candidate(
        vec![source::contact(7, 1, &interface, 1.0)], &[],
    ).unwrap();
    let saved_cells = thermal.cells().clone();
    let saved_circuit = electrical.clone();
    let run = |card: &InterfaceSystemCard, routes: &[ContactHeatCells], controls,
               cancel: &mut dyn FnMut() -> bool| {
        step_with_contacts_candidate(&thermal, 0.2, common::input(&electrical),
            &responses(&cards), routes, controls, conductivity_controls(), &[0.0; 2],
            |state| contact_coefficients(state, &electrical, &cards, &curves, card), cancel)
    };
    assert!(run(&interface, &[], config(), &mut || false).is_err());
    let other = source::card(&cards[0], &cards[1], 3.0, 0.005);
    assert!(run(&other, &ROUTES, config(), &mut || false).is_err());
    assert!(matches!(run(&interface, &ROUTES, fs_couple::circuit_enthalpy::CircuitEnthalpySolve {
        maximum_circuit_trials: 1, ..config()
    }, &mut || false), Err(CircuitEnthalpyError::CircuitBudget)));
    let mut polls = 0;
    run(&interface, &ROUTES, config(), &mut || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(run(&interface, &ROUTES, config(), &mut || { count += 1; count == polls }),
        Err(CircuitEnthalpyError::Cancelled)));
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(electrical, saved_circuit);
}

use super::*;
use crate::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use crate::section_temperature::tests::{card, section};
use crate::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::Density;

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let lo = QueryPoint::new().with("T", 300.0).unwrap().with("pressure", 100000.0).unwrap();
    let hi = lo.clone().with("T", 400.0).unwrap();
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T", lower: &lo, upper: &hi,
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 0.01).unwrap()
}

fn budget() -> DepositionBudget {
    DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 10000 }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w }
}
fn interface() -> SectionThermalInterface {
    SectionThermalInterface { area_m2: 0.001, resistance_geometry_per_m: [1.0, 1.0],
        contact_resistance_k_w: 0.0 }
}
fn config() -> StepConfig {
    StepConfig { max_dt_s: 0.002, energy_tolerance_j: 2e-8 }
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-8 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn electrical_heat_advances_enthalpy_and_refreshes_the_entire_material_section() {
    let (card, req) = card(1.0);
    let chart = caloric(&card);
    let reference = section(&card, &req);
    let law = ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[&card, &card], &["T", "T"]).unwrap();
    let body = HeatedConcentricSection::try_new(law, &[&chart, &chart], 1.0,
        &[boundary(0.0), boundary(0.0)], vec![interface()], budget()).unwrap();
    let start = body.observer.electrical().unwrap();
    let (heated, report) = body.advance_to(1.0, OhmicDrive::Voltage(0.2), config(), || Ok(())).unwrap();
    assert!(heated.thermal().temperature_k(CellId(0)).unwrap() > 300.0);
    assert!(heated.section().axial_stiffness_n() < reference.axial_stiffness_n());
    assert!(heated.section().bending_stiffness_n_m2() < reference.bending_stiffness_n_m2());
    assert_eq!(heated.section().line_density_kg_m(), reference.line_density_kg_m());
    assert!(heated.observer.link_conductances().unwrap()[0] > body.observer.link_conductances().unwrap()[0]);
    assert!(report.electrical_energy_j > 0.0);
    near(report.thermal.stored_energy_change_j, report.electrical_energy_j);
    near(report.thermal.external_heat_j, report.electrical_energy_j);
    assert!(report.thermal.substeps > 1);
    assert_eq!(start.joule_power_w, 0.0);
    assert_eq!(body.thermal().time_s(), 0.0);
    assert_eq!(body.section(), &reference);
    for cell in heated.thermal().cells().values() { assert_eq!(cell.boundary.source_w, 0.0); }
    let (coasted, idle) = heated.advance_to(1.2, OhmicDrive::Current(0.0), config(), || Ok(())).unwrap();
    assert_eq!(idle.electrical_energy_j, 0.0);
    near(coasted.thermal().energy_j(), heated.thermal().energy_j());
}

#[test]
fn a_data_only_resistivity_swap_changes_the_actual_thermal_and_mechanical_history() {
    let run = |scale| {
        let (card, req) = card(scale);
        let chart = caloric(&card);
        let law = ConcentricSectionTemperatureLaw::try_new(section(&card, &req), &[&card, &card], &["T", "T"]).unwrap();
        let body = HeatedConcentricSection::try_new(law, &[&chart, &chart], 1.0,
            &[boundary(0.0), boundary(0.0)], vec![interface()], budget()).unwrap();
        let (next, report) = body.advance_to(1.0, OhmicDrive::Voltage(0.2), config(), || Ok(())).unwrap();
        (next.thermal().temperature_k(CellId(0)).unwrap(), next.section().axial_stiffness_n(), report.electrical_energy_j)
    };
    let a = run(1.0);
    let b = run(2.0);
    assert!(a.0 > b.0);
    assert!(a.1 < b.1);
    assert!(a.2 > b.2);
}

#[test]
fn shared_links_redistribute_layer_heat_without_duplicating_it() {
    let (card, req) = card(1.0);
    let chart = caloric(&card);
    let law = ConcentricSectionTemperatureLaw::try_new(section(&card, &req), &[&card, &card], &["T", "T"]).unwrap();
    let body = HeatedConcentricSection::try_new(law, &[&chart, &chart], 1.0,
        &[boundary(0.0), boundary(1.0)], vec![interface()], budget()).unwrap();
    let (next, report) = body.advance_to(1.0, OhmicDrive::Current(0.0), config(), || Ok(())).unwrap();
    let inner = next.thermal().temperature_k(CellId(0)).unwrap();
    let outer = next.thermal().temperature_k(CellId(1)).unwrap();
    assert!(inner > 300.0 && outer > inner);
    near(report.thermal.stored_energy_change_j, 1.0);
    near(report.thermal.external_heat_j, 1.0);
    assert_eq!(report.electrical_energy_j, 0.0);
    assert_eq!(next.thermal().cells()[&CellId(1)].boundary.source_w, 1.0);
}

#[test]
fn unsupported_endpoint_and_cancellation_preserve_thermal_and_material_state() {
    let (card, req) = card(1.0);
    let chart = caloric(&card);
    let law = ConcentricSectionTemperatureLaw::try_new(section(&card, &req), &[&card, &card], &["T", "T"]).unwrap();
    let body = HeatedConcentricSection::try_new(law, &[&chart, &chart], 1.0,
        &[boundary(0.0), boundary(0.0)], vec![interface()], budget()).unwrap();
    let old_section = body.section().clone();
    let old_cells = body.thermal().cells().clone();
    assert!(body.advance_to(1.0, OhmicDrive::Voltage(1e6), config(), || Ok(())).is_err());
    let mut calls = 0;
    let cancelled = body.advance_to(1.0, OhmicDrive::Voltage(0.2), config(), || {
        calls += 1;
        if calls >= 15 { Err(DepositionError::Cancelled) } else { Ok(()) }
    });
    assert!(matches!(cancelled, Err(SectionHeatingError::Thermal(DepositionError::Cancelled))));
    assert_eq!(body.section(), &old_section);
    assert_eq!(body.thermal().cells(), &old_cells);
    assert_eq!(body.thermal().time_s(), 0.0);
    assert!(body.advance_to(0.02, OhmicDrive::Voltage(0.2), config(), || Ok(())).is_ok());
}

#[test]
fn foreign_caloric_sources_and_invalid_interfaces_are_not_silently_admitted() {
    let (card_a, req) = card(1.0);
    let (card_b, _) = card(2.0);
    let chart_a = caloric(&card_a);
    let chart_b = caloric(&card_b);
    let law = ConcentricSectionTemperatureLaw::try_new(section(&card_a, &req), &[&card_a, &card_a], &["T", "T"]).unwrap();
    let make = |charts: &[&ResolvedHeatCapacityCurve], path| HeatedConcentricSection::try_new(
        law.clone(), charts, 1.0, &[boundary(0.0), boundary(0.0)], path, budget()).is_err();
    assert!(make(&[&chart_a, &chart_b], vec![interface()]));
    assert!(make(&[&chart_a, &chart_a], vec![]));
    let mut invalid = interface();
    invalid.resistance_geometry_per_m[0] = 0.0;
    assert!(make(&[&chart_a, &chart_a], vec![invalid]));
    invalid = interface();
    invalid.area_m2 = 100.0;
    assert!(make(&[&chart_a, &chart_a], vec![invalid]));
}

mod axial;

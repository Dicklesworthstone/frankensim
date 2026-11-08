//! Analytic feedback oracles with synthetic sources, not measured qualifications.
use super::*;
use crate::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, PhaseStateError, SolidLiquidPhase};
use crate::section::{LayerIsotropicParameters, thermal::AxialConstraint};
use crate::state_point::{
    IsotropicElasticParameters, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point,
};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::{DepositionError, Reservoir, StepConfig};

const LENGTH: f64 = 0.5;
const CP: f64 = 100.0;
const PRESSURE: f64 = 100_000.0;
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), PRESSURE).unwrap()
}
fn fixture(qualified: bool) -> (MaterialCard, Vec<ScalarPropertyRequirement>, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, values) in [
        ("mass-density", Density::DIMS, [1000.0, 512.0, 480.0]),
        ("axial-E", Pressure::DIMS, [200e9; 3]),
        ("shear-G", Pressure::DIMS, [80e9; 3]),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, [CP; 3]),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, [20.0; 3]),
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, [1e-6; 3]),
    ] {
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0, 400.0);
        if qualified {
            validity = validity.with_quantity("p", quantity(QuantityKind::Pressure), PRESSURE, PRESSURE);
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, values[0]), (350.0, values[1]), (400.0, values[2])], dims },
            validity, interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "synthetic density feedback oracle".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "free isotropic specimen".into(), revision: 0 },
        claims, vec![]).unwrap();
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new("mass-density", Density::DIMS);
    let (lower, upper) = (point(300.0), point(400.0));
    let chart = ResolvedHeatCapacityCurve::try_new_exact(&card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T", lower: &lower, upper: &upper,
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0).unwrap();
    (card, requirements, chart)
}
fn law<'a>(card: &'a MaterialCard, req: &[ScalarPropertyRequirement], count: usize)
    -> ConcentricSectionTemperatureLaw<'a>
{
    let material = resolve_material_state_point(card, &point(300.0), req,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let parameters = LayerIsotropicParameters { density: req[0].clone(),
        parameters: IsotropicElasticParameters::YoungShear {
            young_modulus: req[1].clone(), shear_modulus: req[2].clone() } };
    let layers = (0..count).map(|i| ConcentricMaterialLayer {
        outer_radius_m: 0.01 * (i + 1) as f64, material: material.clone(),
    }).collect();
    let section = ConcentricElasticSection::try_new_with_isotropic_parameters(
        0.005, layers, &vec![parameters; count],
    ).unwrap();
    ConcentricSectionTemperatureLaw::try_new(section, &vec![card; count], &vec!["T"; count]).unwrap()
}
fn budget() -> DepositionBudget {
    DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 10_000 }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w }
}
fn config(max_dt_s: f64) -> StepConfig {
    StepConfig { max_dt_s, energy_tolerance_j: 1e-8 }
}
fn body<'a>(card: &'a MaterialCard, req: &[ScalarPropertyRequirement],
    chart: &'a ResolvedHeatCapacityCurve, boundary: Boundary) -> HeatedConcentricSection<'a>
{
    HeatedConcentricSection::try_new_with_free_density_expansion(
        law(card, req, 1), chart, LENGTH, boundary, "p", budget(),
    ).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn density_changes_dimensions_resistance_and_work_without_changing_mass_or_source_roles() {
    let (card, req, chart) = fixture(true);
    let reference = law(&card, &req, 1).reference().clone();
    let mass = reference.line_density_kg_m() * LENGTH;
    let initial = body(&card, &req, &chart, boundary(mass * CP * 50.0));
    let (hot, report) = initial.advance_to(1.0, OhmicDrive::Current(0.0), config(0.25), || Ok(())).unwrap();
    let g = hot.density_geometry().unwrap();
    near(hot.thermal().temperature_k(CellId(0)).unwrap(), 350.0);
    near(g.density_kg_m3, 512.0);
    near(g.isotropic_scale, 1.25);
    near(g.length_m, LENGTH * 1.25);
    near(hot.section().inner_radius_m(), 0.005 * 1.25);
    near(hot.section().layers()[0].outer_radius_m, 0.01 * 1.25);
    near(g.volume_m3, mass / 512.0);
    near(g.surface_area_m2, initial.density_geometry().unwrap().surface_area_m2 * 1.25_f64.powi(2));
    near(hot.section().line_density_kg_m() * hot.length_m(), mass);
    assert_eq!(g.mass_kg, mass);
    assert_eq!(hot.thermal().cells()[&CellId(0)].geometry, initial.thermal().cells()[&CellId(0)].geometry);
    assert_eq!(hot.section().isotropic_parameters(), reference.isotropic_parameters());
    assert_eq!(hot.section().density_source_property(0).unwrap().answer().receipt.selected,
        reference.density_source_property(0).unwrap().answer().receipt.selected);
    let electrical = |b: &HeatedConcentricSection<'_>| b.section().thermal_transport().unwrap()
        .with_ohmic_transport(b.length_m()).unwrap().evaluate(OhmicDrive::Current(1.0), 0.0).unwrap().voltage_drop_v;
    near(electrical(&hot), electrical(&initial) / 1.25);
    let account = report.expansion.unwrap();
    near(account.boundary_work_j, PRESSURE * (g.volume_m3 - mass / 1000.0));
    near(account.internal_energy_change_j + account.boundary_work_j, mass * CP * 50.0);
    near(account.supplied_heat_j, report.thermal.external_heat_j);
    near(account.heat_residual_j, 0.0);
    assert_eq!(report.thermal.substeps, 4);
}

#[test]
fn voltage_and_current_feedback_follow_independent_closed_form_heating_oracles() {
    let (card, req, chart) = fixture(true);
    let initial = body(&card, &req, &chart, boundary(0.0));
    let resistance = 1e-6 * LENGTH / initial.section().area_m2();
    let rate = 100.0_f64.powi(2) * resistance / (initial.density_geometry().unwrap().mass_kg * CP);
    // rho/rho0 = 1-a*(T-300). Integrating dT/dt = rate*s or rate/s
    // gives powers 4/3 (voltage) and 2/3 (current) of this density ratio.
    let a = 0.00976;
    let duration = 2.0;
    let exact_v = 300.0 + (1.0 - (1.0 - 4.0 / 3.0 * a * rate * duration).powf(0.75)) / a;
    let exact_i = 300.0 + (1.0 - (1.0 - 2.0 / 3.0 * a * rate * duration).powf(1.5)) / a;
    let (voltage, rv) = initial.advance_to(duration, OhmicDrive::Voltage(100.0 * resistance),
        config(0.01), || Ok(())).unwrap();
    let (current, ri) = initial.advance_to(duration, OhmicDrive::Current(100.0),
        config(0.01), || Ok(())).unwrap();
    assert!((voltage.thermal().temperature_k(CellId(0)).unwrap() - exact_v).abs() < 0.002);
    assert!((current.thermal().temperature_k(CellId(0)).unwrap() - exact_i).abs() < 0.002);
    assert!(rv.electrical_energy_j > ri.electrical_energy_j);
    assert!(rv.thermal.substeps > 100 && ri.thermal.substeps > 100);
    for report in [rv, ri] {
        near(report.electrical_energy_j, report.thermal.external_heat_j);
        let account = report.expansion.unwrap();
        near(account.internal_energy_change_j + account.boundary_work_j, report.electrical_energy_j);
    }
}

#[test]
fn current_surface_area_drives_both_native_convection_and_the_cfl_limit() {
    let (card, req, chart) = fixture(true);
    let mut bath = boundary(0.0);
    bath.ambient_k = 360.0;
    bath.htc_w_m2_k = 1000.0;
    let initial = body(&card, &req, &chart, bath);
    let (first, _) = initial.advance_to(0.01, OhmicDrive::Current(0.0), config(0.01), || Ok(())).unwrap();
    let temperature = first.thermal().temperature_k(CellId(0)).unwrap();
    let expected = 1000.0 * first.density_geometry().unwrap().surface_area_m2 * (360.0 - temperature) * 0.01;
    let (_, second) = first.advance_to(0.02, OhmicDrive::Current(0.0), config(0.01), || Ok(())).unwrap();
    near(second.thermal.external_heat_j, expected);
    let frozen = 1000.0 * initial.density_geometry().unwrap().surface_area_m2 * (360.0 - temperature) * 0.01;
    assert!(expected - frozen > 0.1);
    let (_, cfl) = initial.advance_to(0.35, OhmicDrive::Current(0.0), config(1.0), || Ok(())).unwrap();
    // Fixed area would admit two steps; expanded area requires a third.
    assert_eq!(cfl.thermal.substeps, 3);
    assert_eq!(cfl.electrical_energy_j, 0.0);
    near(cfl.expansion.unwrap().supplied_heat_j, cfl.thermal.external_heat_j);
}

#[test]
fn missing_source_pressure_and_incompatible_constraints_are_refused() {
    let (card, req, chart) = fixture(true);
    let make = |count, b, axis| HeatedConcentricSection::try_new_with_free_density_expansion(
        law(&card, &req, count), &chart, LENGTH, b, axis, budget());
    assert!(make(2, boundary(0.0), "p").is_err());
    assert!(make(1, boundary(0.0), "T").is_err());
    let mut contact = boundary(0.0);
    contact.reservoir = Some(Reservoir { area_m2: 0.001, conductance_w_k: 1.0, temperature_k: 300.0 });
    assert!(make(1, contact, "p").is_err());
    assert!(matches!(make(1, boundary(0.0), "p").unwrap()
        .with_axial_constraint(AxialConstraint::FixedStrain(0.0), 0.01),
        Err(SectionHeatingError::Invalid(_))));
    let (unqualified, other_req, other_chart) = fixture(false);
    assert!(HeatedConcentricSection::try_new_with_free_density_expansion(
        law(&unqualified, &other_req, 1), &other_chart, LENGTH, boundary(0.0), "p", budget()).is_err());
    assert!(HeatedConcentricSection::try_new_with_free_density_expansion(
        law(&card, &req, 1), &other_chart, LENGTH, boundary(0.0), "p", budget()).is_err());
    assert!(HeatedConcentricSection::try_new(law(&card, &req, 1), &[&chart], LENGTH,
        &[boundary(0.0)], vec![], budget()).is_err());
}

#[test]
fn late_source_failure_and_cancellation_roll_back_thermal_material_and_geometry() {
    let (card, req, chart) = fixture(true);
    let mass = law(&card, &req, 1).reference().line_density_kg_m() * LENGTH;
    let initial = body(&card, &req, &chart, boundary(mass * CP * 50.0));
    let (cells, section, geometry) = (initial.thermal().cells().clone(),
        initial.section().clone(), *initial.density_geometry().unwrap());
    let mut evaluated = 0;
    assert!(initial.advance_to(3.0, OhmicDrive::Current(0.0), config(0.1),
        || { evaluated += 1; Ok(()) }).is_err());
    assert!(evaluated > 50);
    let mut checks = 0;
    assert!(matches!(initial.advance_to(1.0, OhmicDrive::Current(0.0), config(0.1), || {
        checks += 1;
        if checks == 40 { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(SectionHeatingError::Thermal(DepositionError::Cancelled))));
    assert_eq!(initial.thermal().cells(), &cells);
    assert_eq!(initial.section(), &section);
    assert_eq!(initial.density_geometry(), Some(&geometry));
    assert_eq!(initial.thermal().time_s(), 0.0);
    assert!(initial.advance_to(0.1, OhmicDrive::Current(0.0), config(0.1), || Ok(())).is_ok());
}

#[test]
fn isobaric_accounting_keeps_accepted_endpoints_and_exposes_unmatched_heat() {
    let (_, _, source) = fixture(true);
    let chart = source.curve();
    let initial = chart.state_at_specific_enthalpy(123.456).unwrap();
    let final_state = chart.state_at_specific_enthalpy(876.543).unwrap();
    let account = chart.account_isobaric_transition(initial, final_state, 2.0, PRESSURE, -17.0).unwrap();
    assert_eq!(account.initial, initial);
    assert_eq!(account.state, final_state);
    near(account.enthalpy_change_j, 2.0 * (876.543 - 123.456));
    near(account.heat_residual_j, -17.0 - account.enthalpy_change_j);
    near(account.boundary_work_j, PRESSURE * 2.0
        * (final_state.bulk_density_kg_m3().recip() - initial.bulk_density_kg_m3().recip()));
    let reversed = chart.account_isobaric_transition(final_state, initial, 2.0, PRESSURE, 17.0).unwrap();
    near(reversed.boundary_work_j, -account.boundary_work_j);
    assert_eq!(chart.account_isobaric_transition(initial, initial, 2.0, PRESSURE, 1.0)
        .unwrap().heat_residual_j, 1.0);
    let foreign_chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x4e; 32]), SolidLiquidPhase::Solid, chart.knots().to_vec()).unwrap();
    let foreign = foreign_chart.state_at_specific_enthalpy(123.456).unwrap();
    for (before, after) in [(foreign, final_state), (initial, foreign)] {
        assert_eq!(chart.account_isobaric_transition(before, after, 2.0, PRESSURE, 0.0),
            Err(PhaseStateError::StateCurveMismatch));
    }
}

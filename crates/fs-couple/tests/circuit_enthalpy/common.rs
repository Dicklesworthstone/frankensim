use std::collections::BTreeMap;
use fs_couple::circuit_enthalpy::*;
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_couple::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use fs_evidence::ValidityDomain;
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::phase::{EquilibriumPhaseState, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

pub fn point(t: f64) -> QueryPoint { QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),t).unwrap().with("lot",7.0).unwrap() }
/// Authored synthetic source curves, not qualified physical material cards.
pub fn material(resistance: f64, cp: f64, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    material_with_conductivity(resistance,cp,alpha,1.0)
}
pub fn material_with_conductivity(resistance: f64, cp: f64, alpha: f64, k: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density",Density::DIMS);
    let resistivity = PropertyKey::new("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS);
    let mut claims = ClaimSet::new();
    for (key,a,b) in [(heat.clone(),cp,cp),(density.clone(),1.0,1.0),
        (resistivity,resistance*(1.0-100.0*alpha),resistance*(1.0+1700.0*alpha)),
        (PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),k,k)] {
        let dims=key.quantity().dims();
        claims.insert_claim(PropertyClaim {key,value:PropertyValue::Curve {
            abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots:vec![(200.0,a),(2000.0,b)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),200.0,2000.0).with("lot",7.0,7.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic circuit thermal fixture".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"authored test".into(),revision:0},claims,vec![]).unwrap();
    let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&heat,density:&density,temperature_axis:"T",lower:&point(200.0),upper:&point(2000.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,caloric)
}
pub fn config() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {maximum_cells:8,maximum_circuit_trials:1000,maximum_newton_iterations:24,
        difference_step_j:1e-3,cell_energy_tolerance_j:1e-6,thermal_energy_tolerance_j:1e-6,
        total_energy_tolerance_j:2e-6,circuit:MaterialCircuitSolve {energy_tolerance_j:1e-8,power_tolerance_w:1e-8},
        newton:NewtonKrylovConfig {absolute_tolerance:1e-8,relative_tolerance:0.0,
            forcing_minimum:1e-4,linear_restart:8,..NewtonKrylovConfig::default()} }
}
pub const DESTINATIONS: [ResistorCell;2] = [
    ResistorCell {branch:1,cell:CellId(30)},ResistorCell {branch:2,cell:CellId(10)}];
pub fn boundary() -> Boundary { Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0} }
pub fn initial<'a>(curves:&'a [ResolvedHeatCapacityCurve;2],conductance:f64,endpoint_event:Option<f64>)
    -> DepositionSimulation<'a> {
    let mut events=Vec::new();
    for (index,id) in [CellId(30),CellId(10)].into_iter().enumerate() {
        // Constant Cp inferred only for this exact authored test fixture.
        let cp=(curves[index].curve().knots().last().unwrap().specific_enthalpy_j_kg
            -curves[index].curve().knots()[0].specific_enthalpy_j_kg)/1800.0;
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id,segment_id:index as u64,body_id:index as u64,material_id:MaterialId(index as u64),
            mass_kg:1.0,specific_h_j_kg:100.0*cp,geometry:CellGeometry {
                start_m:[index as f64,0.0,0.0],end_m:[index as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:6.0},
            boundary:boundary()})});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:CellId(30),b:CellId(10),area_m2:1.0,conductance_w_k:conductance})});
    if let Some(time_s)=endpoint_event { events.push(TimedEvent {time_s,event:Event::Boundary(CellId(30),Boundary {
        source_w:4.0,..boundary()})}); }
    DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],DepositionBudget {
            max_cells:8,max_links:8,max_events:16,max_substeps:1000}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-7},||Ok(())).unwrap().0
}
pub fn coefficients(state:&DepositionSimulation<'_>,cards:&[MaterialCard;2],curves:&[ResolvedHeatCapacityCurve;2])
    -> Result<CircuitThermalCoefficients<Vec<EquilibriumPhaseState>>,JouleHeatingError> {
    let context=point(300.0); let mut resistors=Vec::new(); let mut phases=Vec::new();
    for (index,binding) in DESTINATIONS.iter().enumerate() {
        let cell=&state.cells()[&binding.cell];
        let phase=curves[index].curve().state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|_|JouleHeatingError::Invalid("fixture phase chart refused"))?;
        let source=MaterialOhmicBody {material:ThermalMaterialPoint {material:&cards[index],context:&context,
            temperature_axis:"T",liquid_fraction_axis:None},selection:MaterialPropertySelection::SingleClaimOnly,
            length_m:1.0,area_m2:1.0};
        resistors.push(MaterialResistor {branch:binding.branch,conductor:source.resolve_conductor(phase,cell.geometry.volume_m3)?});
        phases.push(phase);
    }
    Ok(CircuitThermalCoefficients {resistors,sources_w:BTreeMap::new(),conductances_w_k:BTreeMap::new(),source:phases})
}
pub fn circuit(state:&DepositionSimulation<'_>,cards:&[MaterialCard;2],curves:&[ResolvedHeatCapacityCurve;2]) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:1.0}),(1,2,Branch::Resistor {ohms:999.0}),
        (2,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        coefficients(state,cards,curves).unwrap().resistors,vec![],vec![10.0],&[]).unwrap()
}
pub fn input<'a>(circuit:&'a MaterialCircuit) -> CircuitEnthalpyInput<'a> {
    CircuitEnthalpyInput {circuit,ports:&[],destinations:&DESTINATIONS}
}
pub fn near(a:f64,b:f64) { assert!((a-b).abs()<3e-6,"{a} != {b}"); }

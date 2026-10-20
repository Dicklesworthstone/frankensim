use std::collections::BTreeMap;
use core::convert::Infallible;
use fs_couple::strip_enthalpy::{LayerThermalBinding, StripEnthalpyConfig, ThermalCoefficients};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermovisco::{Limits, source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery}};
use fs_qty::{Density, Dims, Pressure, Temperature, Time};
use fs_solid::deposition::{self, LayerGeometry, LayeredStrip};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

pub const LENGTH: f64 = 0.1;
pub const AREA: f64 = 1e-5;
pub const VOLUME: f64 = LENGTH * AREA;
pub const DENSITY: f64 = 1000.0;
pub const MASS: f64 = DENSITY * VOLUME;
pub const LAYERS: [u64; 2] = [7, 11];
pub const CELLS: [CellId; 2] = [CellId(30), CellId(10)];

pub struct Material {
    pub mechanical: ResolvedThermoMaxwell,
    pub caloric: ResolvedHeatCapacityCurve,
}
pub fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap() }
pub fn material(cp: f64, branch_pa: f64, hot_shift: f64, free_at_320: f64) -> Material {
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        ("equilibrium", Pressure::DIMS, vec![(280.0,1e8),(320.0,1e8)]),
        ("branch", Pressure::DIMS, vec![(280.0,branch_pa),(320.0,branch_pa)]),
        ("time", Time::DIMS, vec![(280.0,1.0),(320.0,1.0)]),
        ("modulus_scale", Dims([0;6]), vec![(280.0,1.0),(320.0,1.0)]),
        ("shift", Dims([0;6]), vec![(280.0,1.0),(300.0,1.0),(320.0,hot_shift)]),
        ("free", Dims([0;6]), vec![(280.0,0.0),(320.0,free_at_320)]),
        ("density", Density::DIMS, vec![(280.0,DENSITY),(320.0,DENSITY)]),
        ("cp", SPECIFIC_HEAT_CAPACITY_DIMS, vec![(280.0,cp),(320.0,cp)]),
    ] {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic layered heat-feedback test".into(),
                license:"CC0-1.0".into(),artifact:None} }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"declared test data, not measured material".into(),revision:0},claims,vec![]).unwrap();
    let mechanical = ResolvedThermoMaxwell::try_new(&card,ThermoMaxwellSourceQuery {
        equilibrium_modulus:&PropertyKey::new("equilibrium",Pressure::DIMS),
        branches:&[(PropertyKey::new("branch",Pressure::DIMS),PropertyKey::new("time",Time::DIMS))],
        modulus_scale:&PropertyKey::new("modulus_scale",Dims([0;6])),
        relaxation_scale:&PropertyKey::new("shift",Dims([0;6])),
        free_strain:&PropertyKey::new("free",Dims([0;6])),temperature_axis:"T",
        lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly,
    },Limits {max_strain:0.04,max_strain_rate_s:10.0},||false).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("cp",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,-20.0*cp,1e-8).unwrap();
    Material {mechanical,caloric}
}
pub fn geometry(index: usize) -> LayerGeometry {
    LayerGeometry {area_m2:AREA,centroid_z_m:(index as f64-0.5)*0.001,thickness_m:0.001}
}
pub fn boundary() -> Boundary { Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0} }
pub fn cell(id: CellId, index: usize) -> DepositedCell {
    let z = geometry(index).centroid_z_m;
    DepositedCell {id,segment_id:id.0,body_id:1,material_id:MaterialId(index as u64),
        mass_kg:MASS,specific_h_j_kg:0.0,geometry:CellGeometry {start_m:[0.0,0.0,z],
            end_m:[LENGTH,0.0,z],volume_m3:VOLUME,surface_area_m2:0.00222},boundary:boundary()}
}
pub fn initial(materials: &[Material;2], conductance: f64, later: Vec<TimedEvent>)
    -> (DepositionSimulation<'_>, LayeredStrip<'_>) {
    let mut events:Vec<_> = CELLS.into_iter().enumerate().map(|(i,id)|
        TimedEvent {time_s:0.0,event:Event::Birth(cell(id,i))}).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:5,a:CELLS[0],b:CELLS[1],area_m2:0.001,conductance_w_k:conductance})});
    events.extend(later);
    let heat = DepositionSimulation::new(events,&[
        (MaterialId(0),materials[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),materials[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:1000}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let mut strip = LayeredStrip::new(LENGTH,deposition::Budget {max_layers:8,max_branch_states:32}).unwrap();
    for i in 0..2 {
        strip = strip.deposit(LAYERS[i],geometry(i),materials[i].mechanical.law(),
            heat.temperature_k(CELLS[i]).unwrap()).unwrap();
    }
    (heat,strip)
}
pub fn bindings(m: &[Material;2]) -> [LayerThermalBinding<'_>;2] {
    [0,1].map(|i|LayerThermalBinding {layer:LAYERS[i],cell:CELLS[i],
        material:&m[i].mechanical,caloric:&m[i].caloric})
}
pub fn config() -> StripEnthalpyConfig {
    StripEnthalpyConfig {maximum_cells:8,maximum_structural_trials:1000,maximum_newton_iterations:32,
        difference_step_j:1e-7,cell_energy_tolerance_j:1e-9,thermal_energy_tolerance_j:1e-9,
        total_energy_tolerance_j:2e-9,mechanics:deposition::Config {force_tolerance_n:1e-8,
            moment_tolerance_nm:1e-10,energy_tolerance_j:1e-10,max_strain:0.04,max_rotation_rad:0.1},
        newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-14,
            forcing_minimum:1e-4,linear_restart:8,..NewtonKrylovConfig::default()}}
}
pub fn coefficients(_: &DepositionSimulation<'_>) -> Result<ThermalCoefficients<()>,Infallible> {
    Ok(ThermalCoefficients {sources_w:BTreeMap::new(),conductances_w_k:BTreeMap::new(),source:()})
}
pub fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual-expected).abs() <= tolerance,"actual={actual}, expected={expected}");
}
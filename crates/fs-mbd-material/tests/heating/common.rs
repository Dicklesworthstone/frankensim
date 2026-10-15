use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    InterpolationPolicy, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, DENSITY_PROPERTY,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_mbd_material::thermal::heating::{HeatedBody, ThermalRegion};
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::*;

pub fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("pressure", pressure).unwrap()
}
pub fn card(cp: f64, alpha: f64, alpha_upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name,dims,value,upper) in [
        (DENSITY_PROPERTY,Density::DIMS,2.0,500.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS,cp,500.0),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS,alpha,alpha_upper),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name,dims),
            value: PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.0,value),(upper,value)],dims },
            validity:ValidityDomain::unconstrained().with("T",300.0,upper).with("pressure",1e5,2e5),
            interpolation:InterpolationPolicy::LinearInside,uncertainty:UncertaintyModel::Unstated,
            observations:vec![],provenance:Provenance { source:"synthetic heating fixture, not measurements".into(),
                license:"CC0-1.0".into(),artifact:None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(),phase:"solid".into(),
        process:"fixture".into(),revision:0 },claims,vec![]).unwrap()
}
pub fn caloric(card: &MaterialCard, pressure: f64, phase: SolidLiquidPhase) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
        lower:&point(300.0,pressure),upper:&point(500.0,pressure),selection:SelectionPolicy::SingleClaimOnly,
    },phase,0.0,1e-8).unwrap()
}
pub fn reference(card: &MaterialCard, count: usize) -> ResolvedBody {
    let p = point(300.0,1e5);
    let keys = ["a","b"];
    let regions: Vec<_> = keys[..count].iter().enumerate().map(|(i,&key)| {
        let cube = VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap();
        RegionBinding { key,material:card,point:&p,selection:MaterialPropertySelection::SingleClaimOnly,
            geometry:VolumeMoments::new(1.0,Vec3::new(i as f64,0.0,0.0),cube.inertia_m5()).unwrap() }
    }).collect();
    resolve_body(&regions,count,||false).unwrap()
}
pub fn boundary(source_w: f64) -> Boundary {
    Boundary { ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w }
}
pub fn region<'a>(card:&'a MaterialCard, caloric:&'a ResolvedHeatCapacityCurve,
    key:&'a str, id:u64, watts:f64) -> ThermalRegion<'a> {
    ThermalRegion { key,card,caloric,cell_id:CellId(id),material_id:MaterialId(id+10),
        segment_id:id+20,body_id:1,geometry:CellGeometry { start_m:[id as f64-0.5,0.0,0.0],
            end_m:[id as f64+0.5,0.0,0.0],volume_m3:1.0,surface_area_m2:6.0 },
        initial_specific_enthalpy_j_kg:0.0,boundary:boundary(watts),
        expansion_selection:MaterialPropertySelection::SingleClaimOnly,fixed_point_m:Vec3::default() }
}
pub fn budget() -> DepositionBudget {
    DepositionBudget { max_cells:8,max_links:8,max_events:32,max_substeps:1000 }
}
pub fn config() -> StepConfig { StepConfig { max_dt_s:0.25,energy_tolerance_j:1e-9 } }
pub fn heated<'a>(reference:&'a ResolvedBody, card:&'a MaterialCard, caloric:&'a ResolvedHeatCapacityCurve,
    watts:f64,events:Vec<TimedEvent>) -> HeatedBody<'a> {
    HeatedBody::new(reference,vec![region(card,caloric,"a",0,watts)],events,budget(),config(),||Ok(())).unwrap()
}
pub fn close(a:f64,b:f64) { assert!((a-b).abs()<2e-10*a.abs().max(b.abs()).max(1.0),"{a} != {b}"); }

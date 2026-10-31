//! One real source card supplies both thermal storage and boundary conductivity.
use super::*;
use crate::cavity_phs::piston::{GasPiston, thermal::HeatExchangeSolve};
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim,
    Provenance, SelectionPolicy, UncertaintyModel};
use fs_qty::Density;

pub(in crate::cavity_phs::piston::thermal) fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T",semantic(QuantityKind::AbsoluteTemperature),t).unwrap()
        .with("lot",7.0).unwrap().with("melt",0.0).unwrap()
}
pub(in crate::cavity_phs::piston::thermal) fn material(
    k: Option<f64>,
    typed: bool,
    phase_curve: bool,
) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let quantity = if typed { semantic(QuantityKind::ThermalConductivity) }
        else { QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS) };
    let conductivity = PropertyKey::with_quantity("thermal_conductivity",quantity);
    let cp = PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density",Density::DIMS);
    let mut claims = ClaimSet::new();
    let mut insert = |key:PropertyKey, a:f64, b:f64, phase:bool| {
        let dims = key.quantity().dims();
        let (axis, axis_dims, knots) = if phase {
            ("melt",Dims([0;6]),vec![(0.0,a),(1.0,b)])
        } else { ("T",Temperature::DIMS,vec![(200.0,a),(800.0,b)]) };
        claims.insert_claim(PropertyClaim { key,
            value:PropertyValue::Curve { abscissa:axis.into(),abscissa_dims:axis_dims,knots,dims },
            validity:ValidityDomain::unconstrained()
                .with_quantity("T",semantic(QuantityKind::AbsoluteTemperature),200.0,800.0)
                .with("lot",7.0,7.0).with("melt",0.0,1.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance { source:"synthetic body source".into(),
                license:"CC0-1.0".into(),artifact:None } }).unwrap();
    };
    insert(cp.clone(),500.0,500.0,false);
    insert(density.clone(),4000.0,4000.0,false);
    if let Some(k) = k { insert(conductivity.clone(),k,2.0*k,phase_curve); }
    let card = MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(),
        phase:"declared in thermal chart".into(),process:"authored test".into(),revision:0 },claims,vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&cp,density:&density,temperature_axis:"T",lower:&point(200.0),upper:&point(800.0),
        selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0,1e-5).unwrap();
    (card,conductivity,caloric)
}
fn layers<'a>(card:&'a MaterialCard,key:&'a PropertyKey,point:&'a QueryPoint) -> MaterialConductionLayers<'a> {
    MaterialConductionLayers { material:card,context:point,temperature_axis:"T",liquid_fraction_axis:None,
        conductivity:key,selection:MaterialPropertySelection::SingleClaimOnly,
        area_m2:0.02,gas_distance_m:0.001,body_distance_m:0.02 }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j:[-1.0,1.0],source_tolerance_j:1e-9,
        energy_tolerance_j:1e-7,maximum_iterations:128 }
}

#[test]
fn exact_series_law_keeps_typed_source_receipts_and_other_conditions() {
    for typed in [false,true] {
        let (card,key,caloric)=material(Some(0.2),typed,false);
        let context=point(300.0);
        let body=caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
        let gas=gas(0.75,0.001);
        let sample=layers(&card,&key,&context).resolve(gas.state(),body).unwrap();
        let expected=0.02/(0.001/gas.state().state().thermal_conductivity+0.02/0.3);
        assert!((sample.conductance_w_k/expected-1.0).abs()<1e-12);
        let query=context.clone().with("T",500.0).unwrap();
        let requirement=ScalarPropertyRequirement::try_with_key(&key,ScalarAdmissibility::NonNegative).unwrap();
        let original=resolve_material_state_point(&card,&query,&[requirement],MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert_eq!(sample.source,original);
        assert!(!caloric.receipts().is_empty());
        assert_eq!(body.material_card_identity(),card.content_hash());
    }
}

#[test]
fn conductivity_only_card_swap_changes_the_coupled_trajectory() {
    let run=|k| {
        let (card,key,caloric)=material(Some(k),true,false);
        let context=point(300.0);
        let source=layers(&card,&key,&context);
        let mut body=UniformEnthalpyStepInput { curve:caloric.curve(),
            initial:caloric.curve().state_at_specific_enthalpy(150000.0).unwrap(),
            mass_kg:0.001,volume_m3:0.001/4000.0,surface_area_m2:0.02,
            internal_heat_j:0.0,duration_s:drive(0.0).duration_s };
        let mut piston=GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap();
        let mut heat=0.0;
        for _ in 0..80 {
            let candidate=piston.thermal_candidate(drive(0.0),body,controls(),exchange(),
                |g,b| source.resolve(g,b),||false).unwrap();
            let (next,thermal,report)=candidate.into_parts();
            assert_eq!(thermal.state.material_card_identity(),card.content_hash());
            assert!(report.source_residual_j.abs()<=exchange().source_tolerance_j);
            assert!(report.total_energy_residual_j.abs()<=exchange().energy_tolerance_j);
            assert_eq!(thermal.report,source.resolve(&report.piston.midpoint,report.thermal_midpoint).unwrap().source);
            heat+=report.exchanged_heat_j; body.initial=thermal.state; piston=next;
        }
        (piston,heat)
    };
    let (a,qa)=run(0.2); let (b,qb)=run(2.0);
    assert!(qb>qa && qa>0.0);
    assert!((a.position_m()-b.position_m()).abs()>1e-8);
    assert!((a.gas().state().state().sound_speed-b.gas().state().state().sound_speed).abs()>1e-3);
}

#[test]
fn latent_fraction_coordinate_changes_source_without_changing_temperature() {
    let (card,key,_)=material(Some(0.2),true,true);
    let knot=|h,t,x| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,temperature_k:t,
        liquid_mass_fraction:x,bulk_density_kg_m3:4000.0 };
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![knot(0.0,200.0,0.0),
        knot(100000.0,350.0,0.0),knot(160000.0,350.0,1.0),knot(400000.0,800.0,1.0)]).unwrap();
    let context=point(300.0);
    let mut source=layers(&card,&key,&context);
    source.liquid_fraction_axis=Some("melt");
    let a=curve.state_at_specific_enthalpy(115000.0).unwrap();
    let b=curve.state_at_specific_enthalpy(145000.0).unwrap();
    let gas=gas(0.75,0.001);
    let first=source.resolve(gas.state(),a).unwrap();
    let second=source.resolve(gas.state(),b).unwrap();
    assert_eq!(a.temperature_k(),b.temperature_k());
    assert!(second.conductance_w_k>first.conductance_w_k);
    assert!((first.source.property(key.name()).unwrap().value_si()-0.25).abs()<1e-12);
    assert!((second.source.property(key.name()).unwrap().value_si()-0.35).abs()<1e-12);
    source.liquid_fraction_axis=None;
    assert!(source.resolve(gas.state(),b).is_err());
}

#[test]
fn identity_domains_quantity_and_missing_data_do_not_become_insulation() {
    let (card,key,caloric)=material(Some(0.2),false,false);
    let context=point(300.0); let gas=gas(0.75,0.001);
    let state=caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
    let other=material(Some(2.0),false,false).0;
    assert!(layers(&other,&key,&context).resolve(gas.state(),state).is_err());
    let wrong=PropertyKey::new(key.name(),Temperature::DIMS);
    assert!(layers(&card,&wrong,&context).resolve(gas.state(),state).is_err());
    let outside=context.clone().with("lot",9.0).unwrap();
    assert!(layers(&card,&key,&outside).resolve(gas.state(),state).is_err());
    let mut isolated=layers(&card,&key,&context); isolated.area_m2=0.0;
    assert_eq!(isolated.resolve(gas.state(),state).unwrap().conductance_w_k,0.0);
    let (missing,key,caloric)=material(None,false,false);
    let state=caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
    let mut isolated=layers(&missing,&key,&context); isolated.area_m2=0.0;
    assert!(isolated.resolve(gas.state(),state).is_err());
    let (zero,key,caloric)=material(Some(0.0),false,false);
    let state=caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
    assert_eq!(layers(&zero,&key,&context).resolve(gas.state(),state).unwrap().conductance_w_k,0.0);
}

mod snapshot;

//! Data-only Ohmic source swaps through actual thermal/phase/gas/momentum owners.
use super::*;
use crate::cavity_phs::piston::{GasPiston, thermal::{ConductanceSample, HeatExchangeSolve}};
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use fs_evidence::ValidityDomain;
use fs_material::conductor::resolve_uniform_conductor;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with("lot", 7.0).unwrap().with("melt", 0.0).unwrap().with("omega", 1.0).unwrap()
}
fn fixture(scale: Option<f64>, axis: &str) -> (MaterialCard, ResolvedHeatCapacityCurve, PropertyKey) {
    fixture_with_basis(scale,axis,ElectricalPropertyBasis::Resistivity)
}
fn fixture_with_basis(scale: Option<f64>, axis: &str, basis:ElectricalPropertyBasis)
    -> (MaterialCard, ResolvedHeatCapacityCurve, PropertyKey) {
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new("density", Density::DIMS);
    let k = PropertyKey::with_quantity("thermal_conductivity", quantity(QuantityKind::ThermalConductivity));
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value: PropertyValue| {
        // A scalar claim is a plateau; only a curve has knots to interpolate.
        let interpolation = match &value {
            PropertyValue::Scalar { .. } => InterpolationPolicy::ConstantWithinValidity,
            PropertyValue::Curve { .. } => InterpolationPolicy::LinearInside,
        };
        claims.insert_claim(PropertyClaim { key, value,
            validity: ValidityDomain::unconstrained().with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 200.0, 800.0)
                .with("lot",7.0,7.0).with("melt",0.0,1.0).with("omega",0.0,100.0),
            interpolation, uncertainty: UncertaintyModel::Unstated,
            observations:vec![], provenance: Provenance { source:"synthetic electrothermal data".into(),
                license:"CC0-1.0".into(), artifact:None } }).unwrap();
    };
    insert(cp.clone(), PropertyValue::Scalar { value:500.0, dims:SPECIFIC_HEAT_CAPACITY_DIMS });
    insert(rho.clone(), PropertyValue::Scalar { value:4000.0, dims:Density::DIMS });
    insert(k.clone(), PropertyValue::Curve { abscissa:"T".into(), abscissa_dims:Temperature::DIMS,
        knots:vec![(200.0,0.2),(800.0,0.4)], dims:THERMAL_CONDUCTIVITY_DIMS });
    if let Some(scale) = scale {
        let (dims, knots) = match axis {
            "T" => (Temperature::DIMS, vec![(200.0,0.0125*scale),(800.0,0.05*scale)]),
            "melt" => (Dims([0;6]), vec![(0.0,0.025*scale),(1.0,0.05*scale)]),
            _ => (Dims([0,0,-1,0,0,0]), vec![(0.0,0.025*scale),(100.0,0.05*scale)]),
        };
        insert(PropertyKey::new(basis.property_name(),basis.dimensions()),
            PropertyValue::Curve { abscissa:axis.into(),abscissa_dims:dims,knots,dims:basis.dimensions() });
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry:"authored resistor".into(),
        phase:"declared by chart".into(),process:"synthetic".into(),revision:0 },claims,vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity:&cp,density:&rho,temperature_axis:"T",lower:&point(200.0),upper:&point(800.0),
        selection:SelectionPolicy::SingleClaimOnly },SolidLiquidPhase::Solid,0.0,1e-5).unwrap();
    (card,caloric,k)
}
fn heater<'a>(card: &'a MaterialCard, context: &'a QueryPoint, drive: OhmicDrive) -> MaterialJouleHeating<'a> {
    MaterialJouleHeating { material:ThermalMaterialPoint {material:card,context,temperature_axis:"T",liquid_fraction_axis:None},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:0.01,area_m2:0.000025,drive }
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve,initial:curve.state_at_specific_enthalpy(h).unwrap(),mass_kg:0.001,
        volume_m3:0.00000025,surface_area_m2:0.02,internal_heat_j:0.0,duration_s:drive(0.0).duration_s }
}
fn bounds() -> HeatExchangeSolve {
    HeatExchangeSolve {heat_bracket_j:[-1.0,1.0],source_tolerance_j:1e-9,energy_tolerance_j:1e-7,maximum_iterations:128}
}

#[test]
fn original_resistivity_receipts_and_ideal_drive_signs_are_retained() {
    let (card, caloric, _) = fixture(Some(1.0),"T"); let context=point(300.0);
    let state=caloric.curve().state_at_specific_enthalpy(100000.0).unwrap();
    for drive in [OhmicDrive::Current(5.0),OhmicDrive::Current(-5.0),OhmicDrive::Voltage(50.0),OhmicDrive::Voltage(-50.0)] {
        let heat=heater(&card,&context,drive).resolve(state,0.00000025).unwrap();
        let actual_point=context.clone().with("T",400.0).unwrap();
        let expected=resolve_uniform_conductor(&card,&actual_point,MaterialPropertySelection::SingleClaimOnly,0.01,0.000025).unwrap();
        assert_eq!(heat.source.conductor,expected);
        assert!((heat.power_w-250.0).abs()<1e-10);
        assert_eq!(heat.source.body,state);
        assert!(heat.source.electrical_power_residual_w.abs()<1e-10);
        assert!(heat.source.current_a*heat.source.voltage_v>0.0);
    }
}

fn run(scale: f64, drive_mode: OhmicDrive) -> (GasPiston, f64, f64) {
    let (card,caloric,key)=fixture(Some(scale),"T"); let context=point(300.0);
    let source=heater(&card,&context,drive_mode);
    let layers=source.conduction_layers(&key,MaterialPropertySelection::SingleClaimOnly,0.02,0.001,0.02);
    let mut thermal=body(caloric.curve(),100000.0);
    let mut piston=GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap();
    let initial=piston.clone(); let h0=thermal.initial.specific_enthalpy_j_kg(); let mut energy=0.0;
    for i in 0..80 {
        let heater=MaterialJouleHeating {drive:if i<60 {drive_mode}else{OhmicDrive::Current(0.0)},..source.clone()};
        let proposal=piston.powered_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
            |g,b| layers.resolve(g,b).map_err(JouleHeatingError::Thermal),
            |_,b| heater.resolve(b,thermal.volume_m3),||false).unwrap();
        assert_eq!(proposal.power().source.body,proposal.coupled().step().thermal_midpoint);
        assert_eq!(proposal.power().source.conductor.material(),
            heater.resolve(proposal.power().source.body,thermal.volume_m3).unwrap().source.conductor.material());
        energy+=proposal.step().evaluated_energy_j;
        assert!(proposal.step().total_energy_residual_j.abs()<=bounds().energy_tolerance_j);
        let (coupled,_,_)=proposal.into_parts(); let (next,accepted,_)=coupled.into_parts();
        piston=next; thermal.initial=accepted.state;
        assert_eq!(thermal.initial.material_card_identity(),card.content_hash());
    }
    let total=piston.gas().mass_kg()*(piston.gas().state().specific_internal_energy_j_kg()-initial.gas().state().specific_internal_energy_j_kg())
        +0.5*piston.mass_kg()*piston.velocity_m_s().powi(2)+thermal.mass_kg*(thermal.initial.specific_enthalpy_j_kg()-h0)
        +100000.0*(piston.gas().volume_m3()-initial.gas().volume_m3())-energy;
    assert!(total.abs()<80.0*bounds().energy_tolerance_j);
    assert!(thermal.initial.temperature_k()>400.0);
    (piston,energy,thermal.initial.temperature_k())
}
#[test]
fn resistivity_only_swaps_change_heat_and_motion_oppositely_under_current_and_voltage_control() {
    let (a,qa,ta)=run(1.0,OhmicDrive::Current(5.0)); let (b,qb,tb)=run(2.0,OhmicDrive::Current(5.0));
    assert!(qb>qa && tb>ta);
    assert!(b.position_m()>a.position_m());
    assert!(b.gas().state().state().sound_speed>a.gas().state().state().sound_speed);
    let (a,qa,ta)=run(1.0,OhmicDrive::Voltage(50.0)); let (b,qb,tb)=run(2.0,OhmicDrive::Voltage(50.0));
    assert!(qb<qa && tb<ta);
    assert!(b.position_m()<a.position_m());
    assert!(b.gas().state().state().sound_speed<a.gas().state().state().sound_speed);
}

#[test]
fn latent_phase_resistivity_feedback_matches_independent_implicit_formula_at_fixed_temperature() {
    let (card,_,_)=fixture(Some(1.0),"melt"); let context=point(300.0);
    let knot=|h,t,x| EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,liquid_mass_fraction:x,bulk_density_kg_m3:4000.0};
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![knot(0.0,200.0,0.0),
        knot(100000.0,400.0,0.0),knot(160000.0,400.0,1.0),knot(400000.0,800.0,1.0)]).unwrap();
    let thermal=body(&curve,130000.0); let mut heater=heater(&card,&context,OhmicDrive::Current(5.0));
    heater.material.liquid_fraction_axis=Some("melt");
    let old=GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap();
    let result=old.powered_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        |_,b| Ok(ConductanceSample {conductance_w_k:0.0,source:b}),
        |_,b| heater.resolve(b,thermal.volume_m3),||false).unwrap();
    // R=10*(1+x), W=dt*25*10*(1+x0+W/(2*m*latent)).
    let coefficient=thermal.duration_s*25.0*10.0;
    let expected=coefficient*1.5/(1.0-coefficient/(2.0*thermal.mass_kg*60000.0));
    assert!((result.step().applied_energy_j-expected).abs()<1e-9);
    assert_eq!(result.coupled().thermal().state.temperature_k(),400.0);
    assert!(result.coupled().thermal().state.liquid_mass_fraction()>0.5);
    assert!(result.power().source.conductor.resistance_ohm()>15.0);
    assert_eq!(result.coupled().piston(),&old);
    heater.material.liquid_fraction_axis=None;
    assert!(heater.resolve(result.power().source.body,thermal.volume_m3).is_err());
}

#[test]
fn zero_drive_cannot_hide_missing_sources_wrong_identity_geometry_or_dynamic_coefficients() {
    let (card,caloric,_)=fixture(Some(1.0),"T"); let context=point(300.0);
    let body=body(caloric.curve(),100000.0); let source=heater(&card,&context,OhmicDrive::Current(0.0));
    assert_eq!(source.resolve(body.initial,body.volume_m3).unwrap().power_w,0.0);
    assert!(source.resolve(body.initial,2.0*body.volume_m3).is_err());
    let (other,_,_)=fixture(Some(2.0),"T");
    assert!(heater(&other,&context,OhmicDrive::Current(0.0)).resolve(body.initial,body.volume_m3).is_err());
    let (missing,caloric,_)=fixture(None,"T");
    assert!(heater(&missing,&context,OhmicDrive::Current(0.0)).resolve(caloric.curve().state_at_specific_enthalpy(100000.0).unwrap(),body.volume_m3).is_err());
    let outside=context.clone().with("lot",9.0).unwrap();
    assert!(heater(&card,&outside,OhmicDrive::Current(0.0)).resolve(body.initial,body.volume_m3).is_err());
    let wrong=QueryPoint::new().with_quantity("T",quantity(QuantityKind::Pressure),300.0).unwrap();
    assert!(heater(&card,&wrong,OhmicDrive::Current(0.0)).resolve(body.initial,body.volume_m3).is_err());
    let (dynamic,caloric,_)=fixture(Some(1.0),"omega");
    assert!(heater(&dynamic,&context,OhmicDrive::Current(0.0)).resolve(caloric.curve().state_at_specific_enthalpy(100000.0).unwrap(),body.volume_m3).is_err());
}

mod circuit;

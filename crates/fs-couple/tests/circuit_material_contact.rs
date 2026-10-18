//! Contact-source-only changes propagate through material heat and finite circuit history.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::{self,CircuitThermalCoefficients};
use fs_couple::circuit_enthalpy::material::*;
use fs_conduction::transient::deposition::contact::*;
use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::{interface::ThermalInterfacePoint,source::ThermalMaterialPoint};
use fs_matdb::*;
use fs_qty::{Temperature,QuantitySpec};
use fs_thermal::enthalpy::*;

#[test]
fn changing_only_contact_resistance_changes_both_temperatures_and_capacitor_discharge() {
    let run=|resistance:f64| {
        let (a,ca)=material_with_conductivity(1.0,1.0,0.005,1.0);
        let (b,cb)=material_with_conductivity(3.0,2.0,0.0,2.0);
        let cards=[a,b]; let curves=[ca,cb]; let context=point(300.0);
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=[0,1].map(|i|MaterialCell {cell:DESTINATIONS[i].cell,material_id:MaterialId(i as u64),curve:curves[i].curve(),
            point:ThermalMaterialPoint {material:&cards[i],context:&context,temperature_axis:"T",liquid_fraction_axis:None},
            conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
        let resistors=[0,1].map(|i|CellResistor {destination:DESTINATIONS[i],selection:MaterialPropertySelection::SingleClaimOnly,
            length_m:1.0,area_m2:1.0});
        let links=[ConductionLink {link:7,distance_a_m:0.5,distance_b_m:0.5}];
        let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
        let mut claims=ClaimSet::new(); let pin=claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(RP,RD),value:PropertyValue::Curve {abscissa:"aT".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(200.0,0.8*resistance),(2000.0,4.4*resistance)],dims:RD},
            validity:fs_evidence::ValidityDomain::unconstrained().with("aT",200.0,2000.0).with("bT",200.0,2000.0).with("pressure",200000.0,200000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
            provenance:Provenance {source:"synthetic contact substitution".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
        let source=InterfaceSystemCard::assemble(
            SurfaceSpec {material:cards[0].id().clone(),texture_frame:"side-a".into()},
            SurfaceSpec {material:cards[1].id().clone(),texture_frame:"side-b".into()},
            SystemContext {medium:"fixed contact".into(),third_body:None,environment:"fixed".into(),history:"initial".into()},claims,vec![]).unwrap();
        let cp=QueryPoint::new().with_quantity("aT",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap()
            .with_quantity("bT",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap().with("pressure",200000.0).unwrap();
        let contacts=[MaterialContactLink {link:7,sides:[0,1].map(|i|ContactCell {cell:DESTINATIONS[i].cell,
            material_id:MaterialId(i as u64),curve:curves[i].curve()}),
            point:ThermalInterfacePoint {interface:&source,materials:[&cards[0],&cards[1]],context:&cp,
                temperature_axes:["aT","bT"],liquid_fraction_axes:[None,None],mean_temperature_axis:None},selection:ClaimSelection::Pinned(pin)}];
        let mut thermal=initial(&curves,999.0,None); let mut electrical=circuit(&thermal,&cards,&curves);
        let total=thermal.energy_j()+50.0;
        for i in 1..=40 {
            let next=circuit_enthalpy::step_candidate(&thermal,f64::from(i)*0.05,input(&electrical),config(),&[0.0;2],|state| {
                let mut coefficients=network.resolve(state,||false).map_err(|e|e.to_string())?;
                let contact=apply_material_contacts(state,&contacts,&mut coefficients.conductances_w_k,||false).map_err(|e|e.to_string())?;
                Ok::<_,String>(CircuitThermalCoefficients {resistors:coefficients.resistors,sources_w:coefficients.sources_w,
                    conductances_w_k:coefficients.conductances_w_k,source:(coefficients.source,contact)})
            },||false).unwrap();
            let r=next.report(); let contact=&r.source.1[0];
            near(contact.bulk_conductance_w_k,1.0/0.75);
            let t=next.thermal().temperature_k(CellId(30)).unwrap();
            near(contact.interface.value_m2_k_per_w(),resistance*(1.0+0.002*(t-300.0)));
            near(contact.conductance_w_k,1.0/(0.75+contact.interface.value_m2_k_per_w()));
            assert_eq!(contact.interface.receipt().selected,pin);
            near(next.thermal().energy_j()+r.electrical.stored_after_j,total);
            let (t,c,_)=next.into_parts(); thermal=t; electrical=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),electrical.charges_c()[0])
    };
    let low=run(0.1); let high=run(10.0);
    assert!(low.0>high.0 && low.1<high.1);
    assert!((low.2-high.2).abs()>1e-5,"contact data must alter the temperature-dependent circuit history");
}

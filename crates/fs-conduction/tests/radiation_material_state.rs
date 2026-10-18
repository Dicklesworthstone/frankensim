//! Full physical material conditions reach the existing radiosity owner.
mod support;
#[path="radiation_material/common.rs"]
mod common;
use common::*;
use fs_conduction::{SurfaceEmissivity, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::SelectionPolicy;

#[test]
fn phase_emissivity_retains_pressure_lot_and_original_claim_receipt() {
    let card=card(0.2,0.8); let curve=curve(&card,100.0); let context=context();
    let binding=binding(&card,&context);
    assert!(SurfaceEmissivity::from_card("wall",&card,400.0,SelectionPolicy::SingleClaimOnly).is_err());
    for (h,expected) in [(25.0,0.35),(75.0,0.65)] {
        let phase=curve.state_at_specific_enthalpy(h).unwrap();
        let (emissivity,source)=SurfaceEmissivity::from_thermal_state("wall",binding,phase,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        near(emissivity.value(),expected); assert_eq!(emissivity.temperature_k(),400.0);
        let point=binding.at(phase).unwrap();
        assert_eq!(point.axes()["lot"],7.0); assert_eq!(point.axes()["pressure"],100000.0);
        let original=card.claims().query(SURFACE_EMISSIVITY_PROPERTY,&point,SelectionPolicy::SingleClaimOnly).unwrap();
        assert_eq!(emissivity.receipt(),&original.receipt);
        assert_eq!(source.property(SURFACE_EMISSIVITY_PROPERTY).unwrap().answer().receipt,original.receipt);
        assert_eq!(emissivity.card_identity(),card.content_hash());
    }
}

#[test]
fn refreshed_phase_changes_radiosity_without_reusing_stale_factors_or_geometry() {
    let a=card(0.2,0.8); let ca=curve(&a,100.0); let b=card(0.6,0.6); let cb=curve(&b,100.0);
    let old=enclosure([emissivity(&a,&ca,25.0),emissivity(&b,&cb,-50.0)]);
    support::with_cx(|cx| {
        let prepared=old.prepare_radiosity(cx).unwrap();
        let first=prepared.solve(cx,&[400.0,350.0]).unwrap();
        let new=old.with_emissivities(vec![emissivity(&a,&ca,75.0),emissivity(&b,&cb,-50.0)]).unwrap();
        assert_eq!(old.view_factors(),new.view_factors());
        for (a,b) in old.surfaces().iter().zip(new.surfaces()) {
            assert_eq!(a.face_slots(),b.face_slots()); assert_eq!(a.name(),b.name()); assert_eq!(a.area_m2(),b.area_m2());
        }
        let final_state=new.solve(cx,&[400.0,350.0]).unwrap();
        let expected=STEFAN_BOLTZMANN_W_M2_K4*(400.0_f64.powi(4)-350.0_f64.powi(4))/(1.0/0.65+1.0/0.6-1.0);
        near(final_state.net_outward_heat_w[0],expected);
        near(final_state.net_outward_heat_w[1],-expected);
        assert!(final_state.net_outward_heat_w[0]>first.net_outward_heat_w[0]);
        assert_eq!(first,prepared.solve(cx,&[400.0,350.0]).unwrap());
        assert!(old.with_emissivities(vec![]).is_err());
    });
}

#[test]
fn missing_conditions_foreign_caloric_card_and_nonphysical_emissivity_refuse() {
    let a=card(0.2,0.8); let curve=curve(&a,100.0); let phase=curve.state_at_specific_enthalpy(50.0).unwrap();
    let missing=fs_matdb::QueryPoint::new().with("T",400.0).unwrap().with("melt",0.5).unwrap();
    assert!(SurfaceEmissivity::from_thermal_state("wall",binding(&a,&missing),phase,
        MaterialPropertySelection::SingleClaimOnly).is_err());
    let b=card(0.3,0.8);
    assert!(SurfaceEmissivity::from_thermal_state("wall",binding(&b,&context()),phase,
        MaterialPropertySelection::SingleClaimOnly).is_err());
    for value in [0.0,1.2] {
        let c=card(value,value); let law=common::curve(&c,100.0);
        assert!(SurfaceEmissivity::from_thermal_state("wall",binding(&c,&context()),
            law.state_at_specific_enthalpy(50.0).unwrap(),MaterialPropertySelection::SingleClaimOnly).is_err());
    }
}

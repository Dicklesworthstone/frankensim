//! Same-specimen caloric/transport selection, not composition-derived Cp.
use super::*;
use fs_material::sheet::thermal::*;
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_PROPERTY as CP, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY as K, THERMAL_CONDUCTIVITY_DIMS};

fn thermal_card(interpolation:InterpolationPolicy) -> MaterialCard {
    let mut claims=card(&entries()[..10]).claims().clone();
    for (name,dims,values) in [(CP,SPECIFIC_HEAT_CAPACITY_DIMS,vec![(290.0,1000.0),(300.0,1200.0),(310.0,1300.0)]),
        (K,THERMAL_CONDUCTIVITY_DIMS,vec![(290.0,0.04),(310.0,0.08)])] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Dims([0,0,0,1,0,0]),knots:values,dims},
            interpolation,validity:ValidityDomain::unconstrained().with("T",290.0,310.0).with("specimen",7.0,7.0),
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"synthetic thermal specimen".into(),license:"test-only".into(),artifact:None}
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),phase:"solid".into(),
        process:"same specimen".into(),revision:1},claims,vec![]).unwrap()
}
fn bundle(c:&MaterialCard,specimen:f64) -> ResolvedSheet {
    let mut requirements=Vec::new();
    for (key,_,dims) in entries().into_iter().take(10) {
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let q=QueryPoint::new().with("T",300.0).unwrap().with("specimen",specimen).unwrap();
    ResolvedSheet::try_from_material(resolve_material_state_point(c,&q,&requirements,Pick::SingleClaimOnly).unwrap()).unwrap()
}
fn claim(c:&MaterialCard,key:&str,dims:Dims)->fs_matdb::ClaimId {
    let requirement=ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::StrictlyPositive).unwrap();
    resolve_material_state_point(c,&QueryPoint::new().with("T",300.0).unwrap().with("specimen",7.0).unwrap(),
        &[requirement],Pick::SingleClaimOnly).unwrap().property(key).unwrap().answer().receipt.selected
}
#[test]
fn integrates_selected_cp_and_retains_complete_sheet_and_sample_context() {
    let c=thermal_card(InterpolationPolicy::LinearInside);let original=bundle(&c,7.0);
    let caloric=SheetCaloric::try_new(original.clone(),&c,claim(&c,CP,SPECIFIC_HEAT_CAPACITY_DIMS),
        [290.0,310.0],-30.0,0.01).unwrap();
    assert_eq!(caloric.sheet(),&original);
    assert_eq!(caloric.heat_capacity_samples().len(),3);
    for s in caloric.heat_capacity_samples() {
        assert_eq!(s.card_identity(),c.content_hash());assert_eq!(s.query_point()["specimen"],7.0);
    }
    // Exact source integral: 10*(1000+1200)/2 + 10*(1200+1300)/2.
    let a=caloric.state_at_temperature(290.0).unwrap();let b=caloric.state_at_temperature(310.0).unwrap();
    close(b.specific_enthalpy_j_kg()-a.specific_enthalpy_j_kg(),23500.0);
    assert_eq!(b.liquid_mass_fraction(),0.0);
    assert!(caloric.state_at_temperature(310.1).is_err());
    let k=conductivity_at_temperature(&original,&c,claim(&c,K,THERMAL_CONDUCTIVITY_DIMS),305.0).unwrap();
    close(k.property(K).unwrap().value_si(),0.07);
    assert_eq!(k.query_point()["specimen"],7.0);
}
#[test]
fn exact_knots_do_not_authorize_transient_interpolation_or_other_specimens() {
    let exact=thermal_card(InterpolationPolicy::TabulatedOnly);let cp=claim(&exact,CP,SPECIFIC_HEAT_CAPACITY_DIMS);
    assert!(SheetCaloric::try_new(bundle(&exact,7.0),&exact,cp,[290.0,310.0],0.0,0.01).is_err());
    let c=thermal_card(InterpolationPolicy::LinearInside);let cp=claim(&c,CP,SPECIFIC_HEAT_CAPACITY_DIMS);
    assert!(SheetCaloric::try_new(bundle(&c,8.0),&c,cp,[290.0,310.0],0.0,0.01).is_err());
    assert!(SheetCaloric::try_new(bundle(&c,7.0),&c,cp,[289.0,310.0],0.0,0.01).is_err());
    assert!(SheetCaloric::try_new(bundle(&c,7.0),&exact,cp,[290.0,310.0],0.0,0.01).is_err());
}
#[test]
fn conductivity_is_not_a_substitute_for_caloric_storage() {
    let c=thermal_card(InterpolationPolicy::LinearInside);let k=claim(&c,K,THERMAL_CONDUCTIVITY_DIMS);
    let sheet=bundle(&c,7.0);
    assert!(conductivity_at_temperature(&sheet,&c,k,300.0).is_ok());
    assert!(SheetCaloric::try_new(sheet.clone(),&c,k,[290.0,310.0],0.0,0.01).is_err());
    assert!(conductivity_at_temperature(&sheet,&c,k,320.0).is_err());
}

//! Synthetic interface fixtures, never material-name defaults.
use fs_evidence::ValidityDomain;
use fs_material::electrical_interface::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::*;
use fs_qty::{Dims,Temperature};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;

pub fn point() -> QueryPoint {
    QueryPoint::new().with("Ta",300.0).unwrap().with("Tb",300.0).unwrap()
        .with("pressure",100000.0).unwrap().with("melt_a",0.0).unwrap().with("melt_b",0.0).unwrap()
}
pub fn card(a:&MaterialCard,b:&MaterialCard,r:f64,slope:f64) -> InterfaceSystemCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve {abscissa:"Ta".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,r*(1.0-100.0*slope)),(2000.0,r*(1.0+1700.0*slope))],dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained().with("Ta",200.0,2000.0).with("Tb",200.0,2000.0)
            .with("pressure",100000.0,100000.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"synthetic electrical contact".into(),license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    InterfaceSystemCard::assemble(SurfaceSpec {material:a.id().clone(),texture_frame:"first finish".into()},
        SurfaceSpec {material:b.id().clone(),texture_frame:"second finish".into()},
        SystemContext {medium:"authored contact".into(),third_body:None,environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap()
}
pub fn contact(id:u64,branch:usize,card:&InterfaceSystemCard,area:f64) -> SeriesElectricalContact {
    SeriesElectricalContact::new(id,branch,resolve_electrical_interface_state(card,&point(),
        MaterialPropertySelection::SingleClaimOnly).unwrap(),area).unwrap()
}
pub fn scalar_card(property:&str,dims:Dims,value:f64) -> MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(property,dims),value:PropertyValue::Scalar {value,dims},
        validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance {source:"synthetic field storage".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"test".into(),phase:"solid".into(),process:"authored".into(),revision:0},claims,vec![]).unwrap()
}

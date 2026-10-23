#![allow(dead_code)]
use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{Dims, QuantitySpec, Temperature};

pub fn bulk(name: &str) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "supplied chart".into(),
        process: "synthetic".into(), revision: 0 }, ClaimSet::new(), vec![]).unwrap()
}
pub fn curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),
        [(-100.0,300.0,0.0),(0.0,400.0,0.0),(100.0,400.0,1.0),(500.0,800.0,1.0)]
        .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
            liquid_mass_fraction:f, bulk_density_kg_m3:1.0 }).to_vec()).unwrap()
}
pub fn point() -> QueryPoint {
    QueryPoint::new().with_quantity("aT",QuantitySpec::dimensional(Temperature::DIMS),400.0).unwrap()
        .with_quantity("bT",QuantitySpec::dimensional(Temperature::DIMS),350.0).unwrap()
        .with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),375.0).unwrap()
        .with_quantity("fa",QuantitySpec::dimensional(Dims::NONE),0.0).unwrap()
        .with_quantity("fb",QuantitySpec::dimensional(Dims::NONE),0.0).unwrap()
        .with("pressure",200000.0).unwrap().with("lot",7.0).unwrap()
}
pub fn claim(axis: &str, dims: Dims, scale: f64) -> PropertyClaim {
    PropertyClaim { key:PropertyKey::new(RP,RD), value:PropertyValue::Curve {
        abscissa:axis.into(), abscissa_dims:dims, knots:vec![(0.0,scale),(1.0,2.0*scale)], dims:RD },
        validity:ValidityDomain::unconstrained()
            .with_quantity("aT",QuantitySpec::dimensional(Temperature::DIMS),300.0,800.0)
            .with_quantity("bT",QuantitySpec::dimensional(Temperature::DIMS),300.0,800.0)
            .with_quantity("fa",QuantitySpec::dimensional(Dims::NONE),0.0,1.0)
            .with_quantity("fb",QuantitySpec::dimensional(Dims::NONE),0.0,1.0)
            .with("pressure",200000.0,200000.0).with("lot",7.0,7.0),
        uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::LinearInside,
        observations:vec![], provenance:Provenance { source:"synthetic two-sided contact".into(),
            license:"CC0-1.0".into(), artifact:None } }
}
pub fn interface(cards: &[MaterialCard;2], claims: ClaimSet) -> InterfaceSystemCard {
    InterfaceSystemCard::assemble(
        SurfaceSpec { material:cards[0].id().clone(), texture_frame:"first surface".into() },
        SurfaceSpec { material:cards[1].id().clone(), texture_frame:"second surface".into() },
        SystemContext { medium:"declared contact".into(), third_body:None,
            environment:"declared environment".into(), history:"initial".into() },claims,vec![]).unwrap()
}
pub fn source(cards: &[MaterialCard;2], scale: f64) -> (InterfaceSystemCard, ClaimId) {
    let mut claims=ClaimSet::new(); let id=claims.insert_claim(claim("fa",Dims::NONE,scale)).unwrap();
    (interface(cards,claims),id)
}
pub fn binding<'a>(cards:&'a [MaterialCard;2], source:&'a InterfaceSystemCard, point:&'a QueryPoint)
    -> ThermalInterfacePoint<'a> {
    ThermalInterfacePoint { interface:source, materials:[&cards[0],&cards[1]], context:point,
        temperature_axes:["aT","bT"], liquid_fraction_axes:[Some("fa"),Some("fb")],
        mean_temperature_axis:Some("T") }
}

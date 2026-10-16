//! Synthetic cards exercised through the real material and interface resolvers.
#![allow(dead_code)]
use std::collections::BTreeMap;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::source::*;
use fs_contact::normal_patch::InputUncertainty;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterfaceSystemCard,
    InterpolationPolicy, LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point, resolve_isotropic_solid_state_point};
use fs_qty::{Density, Dims, Pressure};
use fs_tribo::InputAuthority;

fn provenance(source: &str) -> Provenance {
    Provenance { source: format!("synthetic point-contact fixture: {source}"),
        license: "CC0-1.0".into(), artifact: None }
}
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim { key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T",280.0,320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: provenance(name) }
}
pub fn solid(name: &str, density: f64, young: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key,value,dims) in [("density",density,Density::DIMS),
        ("young_modulus",young,Pressure::DIMS),("poisson_ratio",0.25,Dims::NONE),
        ("yield_stress",1e10,Pressure::DIMS)] {
        claims.insert_claim(claim(key,value,dims)).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry:name.into(),phase:"solid".into(),
        process:"synthetic test".into(),revision:0 },claims,vec![]).unwrap()
}
pub fn query() -> QueryPoint { QueryPoint::new().with("T",293.15).unwrap() }
pub fn bind(a: &MaterialCard, b: &MaterialCard, chi: Option<f64>, adhesion: f64) -> BoundNormalContactModel {
    let mut parameters = BTreeMap::new();
    for (name,value,dims) in [("characteristic-rate",1.0,Dims([1,0,-1,0,0,0])),
        ("max-patch-to-radius",0.1,Dims::NONE),("max-strain",0.01,Dims::NONE),
        ("max-patch-to-depth",0.1,Dims::NONE),("max-patch-to-layer",0.1,Dims::NONE),
        ("max-pressure-to-yield",0.2,Dims::NONE),("max-rate-ratio",1.0,Dims::NONE)] {
        parameters.insert(name.into(),LawParameter{value,dims});
    }
    if let Some(value)=chi { parameters.insert("dissipation".into(),
        LawParameter{value,dims:Dims([-1,0,1,0,0,0])}); }
    let model = ConstitutiveModelCard {
        law:LawId(if chi.is_some(){NORMAL_HUNT_CROSSLEY_POINT_LAW_ID}else{NORMAL_HERTZ_LAW_ID}.into()),
        law_version:NORMAL_CONTACT_LAW_VERSION, parameters,
        state_schema_version:NORMAL_CONTACT_STATE_SCHEMA_VERSION,
        initial_state:InitialStatePolicy::ZeroInternalState,
        validity:ValidityDomain::unconstrained().with("T",280.0,320.0),
        sources:vec![],provenance:provenance("normal model"),
    };
    let mut claims=ClaimSet::new();
    claims.insert_claim(claim(ADHESION_ENERGY_PROPERTY,adhesion,ADHESION_ENERGY_DIMS)).unwrap();
    let card=InterfaceSystemCard::assemble(
        SurfaceSpec{material:a.id().clone(),texture_frame:"surface-a".into()},
        SurfaceSpec{material:b.id().clone(),texture_frame:"surface-b".into()},
        SystemContext{medium:"dry".into(),third_body:None,environment:"test".into(),history:"initial".into()},
        claims,vec![model]).unwrap();
    let p=query();
    let sa=resolve_isotropic_solid_state_point(a,&p,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let sb=resolve_isotropic_solid_state_point(b,&p,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let dry=bind_dry_interface_system_card(&card,InputAuthority::SyntheticFixture).unwrap();
    let elastic=bind_isotropic_elastic_interface(&dry,&sa,&sb).unwrap();
    let interface=resolve_interface_state_point(&card,&p,&[
        ScalarPropertyRequirement::try_new(ADHESION_ENERGY_PROPERTY,ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative).unwrap()],MaterialPropertySelection::SingleClaimOnly).unwrap();
    let normal=bind_normal_interface_state(&elastic,&interface).unwrap();
    bind_normal_contact_model(&normal,NormalContactModelSelection::SingleSupported).unwrap()
}
pub fn envelope() -> PointContactEnvelope {
    PointContactEnvelope { maximum_indentation_m:1e-4,maximum_absolute_rate_m_per_s:0.2,
        maximum_step_s:0.001,half_space_depth_m:1.0,layer_thickness_m:0.1,
        uncertainty:InputUncertainty{radius_relative:0.01,modulus_relative:0.02,load_relative:0.03} }
}
pub fn source(young: f64, chi: Option<f64>) -> SourcePointContact {
    let a=solid("A",1000.0,young); let b=solid("B",2000.0,2e7);
    SourcePointContact::new(bind(&a,&b,chi,0.0),
        PointContactGeometry::SpherePlane{effective_radius_m:0.1},envelope()).unwrap()
}
pub fn close(a:f64,b:f64) {
    assert!((a-b).abs() <= 2e-11*a.abs().max(b.abs()).max(1e-20),"{a} != {b}");
}

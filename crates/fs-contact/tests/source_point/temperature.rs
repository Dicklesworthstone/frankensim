//! Separate material temperatures flow through the original contact evaluator.
use super::*;
use fs_contact::interface_binding::{InterfaceBindingError, NormalContactModelSelection,
    bind_dry_interface_system_card, bind_isotropic_elastic_interface};
use fs_contact::normal_patch::source::state::{PointContactMaterialCards,ContactMaterialResolveError};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,MaterialCard,MaterialStateId,InterfaceSystemCard,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,UncertaintyModel,InterpolationPolicy,QueryPoint};
use fs_material::state_point::{MaterialPropertySelection,resolve_isotropic_solid_state_point};
use fs_qty::{Density,Pressure,Temperature,Dims};

fn thermal_card(name:&str,e280:f64,e340:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),
        ("young_modulus",e280,Pressure::DIMS),("poisson_ratio",0.25,Dims::NONE),
        ("yield_stress",1e10,Pressure::DIMS)] {
        let curve=name=="young_modulus";
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if curve {PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,e280),(340.0,e340)],dims}} else {PropertyValue::Scalar {value,dims}},
            validity:ValidityDomain::unconstrained().with("T",280.0,340.0),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:if curve {InterpolationPolicy::LinearInside} else {InterpolationPolicy::ConstantWithinValidity},
            observations:vec![],provenance:Provenance {source:"synthetic separate-temperature fixture".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:name.into(),phase:"solid".into(),
        process:"synthetic test".into(),revision:0},claims,vec![]).unwrap()
}
fn cards<'a>(a:&'a MaterialCard,b:&'a MaterialCard,i:&'a InterfaceSystemCard)->PointContactMaterialCards<'a> {
    PointContactMaterialCards {surface_a:a,surface_b:b,interface:i,
        bulk_selection:[MaterialPropertySelection::SingleClaimOnly,MaterialPropertySelection::SingleClaimOnly],
        interface_selection:MaterialPropertySelection::SingleClaimOnly,
        model_selection:NormalContactModelSelection::SingleSupported,authority:InputAuthority::SyntheticFixture}
}
fn point(t:f64)->QueryPoint {QueryPoint::new().with("T",t).unwrap()}
fn geometry()->PointContactGeometry {PointContactGeometry::SpherePlane {effective_radius_m:0.1}}
fn interface(a:&MaterialCard,b:&MaterialCard)->InterfaceSystemCard {
    bind(a,b,Some(0.4),0.0).normal_state().elastic().interface().card().clone()
}

#[test]
fn common_temperature_limit_preserves_original_source_and_strict_entry_point() {
    let a=thermal_card("A",4e7,1e7);let b=thermal_card("B",2e7,2.6e7);
    let model=bind(&a,&b,Some(0.4),0.0);let i=interface(&a,&b);
    let expected=SourcePointContact::new(model,geometry(),envelope()).unwrap();
    let c=cards(&a,&b,&i);let p=query();
    assert_eq!(c.resolve_at_temperatures([&p,&p],&p,geometry(),envelope()).unwrap(),expected);
    assert_eq!(c.resolve(&p,geometry(),envelope()).unwrap(),expected);
    let sa=resolve_isotropic_solid_state_point(&a,&point(290.0),MaterialPropertySelection::SingleClaimOnly).unwrap();
    let sb=resolve_isotropic_solid_state_point(&b,&point(310.0),MaterialPropertySelection::SingleClaimOnly).unwrap();
    let dry=bind_dry_interface_system_card(&i,InputAuthority::SyntheticFixture).unwrap();
    assert!(matches!(bind_isotropic_elastic_interface(&dry,&sa,&sb),Err(InterfaceBindingError::StatePointMismatch)));
}

#[test]
fn each_bulk_retains_its_temperature_and_original_force_pressure_and_loss() {
    let a=thermal_card("A",4e7,1e7);let b=thermal_card("B",2e7,2.6e7);let i=interface(&a,&b);
    let c=cards(&a,&b,&i);let pa=point(290.0);let pb=point(310.0);let pi=point(300.0);
    let law=c.resolve_at_temperatures([&pa,&pb],&pi,geometry(),envelope()).unwrap();
    let e=law.model().normal_state().elastic();
    close(e.surface_a_state().young_modulus_pa(),3.5e7);
    close(e.surface_b_state().young_modulus_pa(),2.3e7);
    assert_eq!(e.surface_a_state().resolved().query_point(),&[("T".into(),290.0)]);
    assert_eq!(e.surface_b_state().resolved().query_point(),&[("T".into(),310.0)]);
    assert_eq!(e.state_coordinate("T"),Some(300.0));
    let reduced=1.0/(0.9375/3.5e7+0.9375/2.3e7);
    close(e.reduced_modulus_pa(),reduced);
    let k=(4.0/3.0)*0.1_f64.sqrt()*reduced;let d=2e-5_f64;let rate=0.05;
    let response=law.evaluate("hot-cold",d,rate,1e-4).unwrap();
    close(response.normal_force_n,k*d.powf(1.5)*(1.0+0.4*rate));
    close(response.reversible_energy_j,0.4*k*d.powf(2.5));
    close(response.dissipated_power_w,0.4*k*d.powf(1.5)*rate*rate);
    close(response.patch_radius_m,(0.1*d).sqrt());
    assert_eq!(response.authority,InputAuthority::SyntheticFixture);
}

#[test]
fn same_mean_temperature_is_not_an_effective_material_or_interface_temperature() {
    let a=thermal_card("A",4e7,1e7);let b=thermal_card("B",2e7,2.6e7);let i=interface(&a,&b);
    let c=cards(&a,&b,&i);let cold=point(290.0);let hot=point(310.0);let mean=point(300.0);
    let x=c.resolve_at_temperatures([&cold,&hot],&mean,geometry(),envelope()).unwrap();
    let y=c.resolve_at_temperatures([&hot,&cold],&mean,geometry(),envelope()).unwrap();
    let averaged=c.resolve(&mean,geometry(),envelope()).unwrap();
    assert_ne!(x.stiffness_n_per_m_3_2(),y.stiffness_n_per_m_3_2());
    assert_ne!(x.stiffness_n_per_m_3_2(),averaged.stiffness_n_per_m_3_2());
    let interface_changed=c.resolve_at_temperatures([&cold,&hot],&point(305.0),geometry(),envelope()).unwrap();
    assert_eq!(x.stiffness_n_per_m_3_2(),interface_changed.stiffness_n_per_m_3_2());
    assert_ne!(x.model().identity(),interface_changed.model().identity());
    // True elliptic geometry still goes through the existing finite-patch owner.
    let ellipse=PointContactGeometry::EllipticParaboloid {maximum_curvature_m_inverse:20.0,minimum_curvature_m_inverse:5.0};
    let ex=c.resolve_at_temperatures([&cold,&hot],&mean,ellipse,envelope()).unwrap();
    let ey=c.resolve_at_temperatures([&hot,&cold],&mean,ellipse,envelope()).unwrap();
    let rx=ex.evaluate("ellipse",1e-5,0.0,1e-4).unwrap();let ry=ey.evaluate("ellipse",1e-5,0.0,1e-4).unwrap();
    close(rx.normal_force_n/ry.normal_force_n,x.stiffness_n_per_m_3_2()/y.stiffness_n_per_m_3_2());
    let ax=rx.elliptic_patch_axes.unwrap();let ay=ry.elliptic_patch_axes.unwrap();
    close(ax.semi_major_axis_m,ay.semi_major_axis_m);
    close(ax.semi_minor_axis_m,ay.semi_minor_axis_m);
}

#[test]
fn valid_interface_temperature_cannot_hide_a_bulk_outside_the_model_domain() {
    // Bulk data supports 340 K, but the selected normal model supports only 320 K.
    let a=thermal_card("A",4e7,1e7);let b=thermal_card("B",2e7,2.6e7);let i=interface(&a,&b);
    let c=cards(&a,&b,&i);
    for (ta,tb) in [(330.0,300.0),(300.0,330.0)] {
        assert!(matches!(c.resolve_at_temperatures([&point(ta),&point(tb)],&point(300.0),geometry(),envelope()),
            Err(ContactMaterialResolveError::Interface(InterfaceBindingError::NormalContactModelOutsideValidity {..}))));
    }
    assert!(c.resolve_at_temperatures([&point(350.0),&point(300.0)],&point(300.0),geometry(),envelope()).is_err());
    assert!(c.resolve_at_temperatures([&point(290.0),&point(310.0)],&point(330.0),geometry(),envelope()).is_err());
}

#[test]
fn role_and_nonthermal_context_mismatch_remain_errors() {
    let a=thermal_card("A",4e7,1e7);let b=thermal_card("B",2e7,2.6e7);let i=interface(&a,&b);
    let c=cards(&a,&b,&i);let p=point(300.0);
    let extra=p.clone().with("pressure",100_000.0).unwrap();
    assert!(c.resolve_at_temperatures([&p,&p],&extra,geometry(),envelope()).is_err());
    assert!(cards(&b,&a,&i).resolve_at_temperatures([&p,&p],&p,geometry(),envelope()).is_err());
    assert!(c.resolve_at_temperatures([&p,&p],&QueryPoint::new(),geometry(),envelope()).is_err());
    assert!(c.resolve_at_temperatures([&p,&p],&point(0.0),geometry(),envelope()).is_err());
}

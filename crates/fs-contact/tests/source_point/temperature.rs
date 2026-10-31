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
    assert_eq!(c.resolve_at_fixed_environments([&p,&p],&p,geometry(),envelope()).unwrap(),expected);
    let sa=resolve_isotropic_solid_state_point(&a,&point(290.0),MaterialPropertySelection::SingleClaimOnly).unwrap();
    let sb=resolve_isotropic_solid_state_point(&b,&point(310.0),MaterialPropertySelection::SingleClaimOnly).unwrap();
    let dry=bind_dry_interface_system_card(&i,InputAuthority::SyntheticFixture).unwrap();
    assert!(matches!(bind_isotropic_elastic_interface(&dry,&sa,&sb),Err(InterfaceBindingError::StatePointMismatch)));
}

#[test]
fn selected_bulk_pairs_retain_complete_sources_in_native_point_contact() {
    use fs_material::state_point::{IsotropicElasticParameters,IsotropicSolidStatePoint,
        ScalarAdmissibility as Admit,ScalarPropertyRequirement as Requirement,resolve_material_state_point};
    let a=thermal_card("A",4e7,4e7);
    let base=thermal_card("B",2e7,2e7);
    let mut claims=base.claims().clone();
    for (name,value,dims) in [("shear_modulus",8e6,Pressure::DIMS),
        ("bulk_modulus",2e7/1.5,Pressure::DIMS),("retained_viscosity",10.0,fs_qty::DynViscosity::DIMS)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,340.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic contact source pair".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    let b=MaterialCard::assemble(base.id().clone(),claims,vec![]).unwrap();
    let i=interface(&a,&b);let c=cards(&a,&b,&i);let p=query();
    let expected=c.resolve(&p,geometry(),envelope()).unwrap();
    let sa=resolve_isotropic_solid_state_point(&a,&p,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let sb=resolve_isotropic_solid_state_point(&b,&p,MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert_eq!(sa.identity(),sa.resolved().identity());
    assert_eq!(c.resolve_from_bulk_states_at_temperatures([&sa,&sb],&p,geometry(),envelope()).unwrap(),expected);
    let requirement=|name,dims|Requirement::try_new(name,dims,Admit::StrictlyPositive).unwrap();
    let density=requirement("density",Density::DIMS);
    let yield_requirement=requirement("yield_stress",Pressure::DIMS);
    for parameters in [IsotropicElasticParameters::YoungShear {
        young_modulus:requirement("young_modulus",Pressure::DIMS),
        shear_modulus:requirement("shear_modulus",Pressure::DIMS),
    },IsotropicElasticParameters::BulkShear {
        bulk_modulus:requirement("bulk_modulus",Pressure::DIMS),
        shear_modulus:requirement("shear_modulus",Pressure::DIMS),
    }] {
        let mut requirements=vec![density.clone(),yield_requirement.clone(),
            requirement("retained_viscosity",fs_qty::DynViscosity::DIMS)];
        requirements.extend(parameters.requirements().into_iter().cloned());
        let resolved=resolve_material_state_point(&b,&p,&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap();
        let selected=IsotropicSolidStatePoint::try_from_resolved_with_parameters(
            resolved.clone(),&density,&parameters,&yield_requirement).unwrap();
        let law=c.resolve_from_bulk_states_at_temperatures([&sa,&selected],&p,geometry(),envelope()).unwrap();
        close(law.stiffness_n_per_m_3_2(),expected.stiffness_n_per_m_3_2());
        let response=law.evaluate("selected-pair",2e-5,0.05,1e-4).unwrap();
        let reference=expected.evaluate("canonical",2e-5,0.05,1e-4).unwrap();
        close(response.normal_force_n,reference.normal_force_n);
        close(response.pressure.peak_pressure_pa,reference.pressure.peak_pressure_pa);
        let elastic=law.model().normal_state().elastic();
        assert_eq!(elastic.surface_b_state().resolved(),&resolved);
        assert_eq!(elastic.surface_b_state_identity(),selected.identity());
        assert!(resolved.property("poisson_ratio").is_none());
        assert_eq!(c.resolve_from_bulk_states_at_fixed_environments([&sa,&selected],&p,
            geometry(),envelope()).unwrap(),law);
    }
    let foreign=resolve_isotropic_solid_state_point(&base,&p,MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(c.resolve_from_bulk_states_at_temperatures([&sa,&foreign],&p,geometry(),envelope()).is_err(),
        "same material ID must not substitute a different immutable card");
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
        assert!(matches!(c.resolve_at_fixed_environments([&point(ta),&point(tb)],&point(300.0),geometry(),envelope()),
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


mod environment {
    use super::*;
    use fs_qty::{QuantitySpec, semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm}};

    fn moisture(basis: MoistureBasis) -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Moisture(basis), ValueForm::Static))
    }
    fn material(name: &str, axis: &str, quantity: QuantitySpec, e0: f64, e1: f64) -> MaterialCard {
        let base = thermal_card(name, e0, e0);
        let mut claims = ClaimSet::new();
        for id in base.claims().observation_ids() {
            claims.register_observation(base.claims().observation(id).unwrap().clone()).unwrap();
        }
        for (_, original) in base.claims().claims_ordered() {
            let mut claim = original.clone();
            claim.validity = claim.validity.with_quantity(axis, quantity, 0.0, 0.2);
            if claim.key.name() == "young_modulus" {
                claim.value = PropertyValue::Curve { abscissa: axis.into(), abscissa_dims: quantity.dims(),
                    knots: vec![(0.0, e0), (0.2, e1)], dims: Pressure::DIMS };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(claim).unwrap();
        }
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
    }
    fn ordered() -> InterfaceSystemCard {
        let base = interface(&thermal_card("A",4e7,4e7), &thermal_card("B",2e7,2e7));
        let mut claims = ClaimSet::new();
        for id in base.claims().observation_ids() {
            claims.register_observation(base.claims().observation(id).unwrap().clone()).unwrap();
        }
        for (_, original) in base.claims().claims_ordered() {
            let mut claim = original.clone();
            claim.validity = claim.validity.with_quantity("humidity",
                moisture(MoistureBasis::RelativeHumidity), 0.0, 1.0);
            claims.insert_claim(claim).unwrap();
        }
        InterfaceSystemCard::assemble(base.surface_a().clone(), base.surface_b().clone(),
            base.context().clone(), claims, base.models().to_vec()).unwrap()
    }

    #[test]
    fn independent_material_and_interface_states_reach_original_force_pressure_and_loss() {
        let qa = moisture(MoistureBasis::DryMass);
        let qb = moisture(MoistureBasis::WetMass);
        let a = material("A", "water_a", qa, 4e7, 2e7);
        let b = material("B", "water_b", qb, 2e7, 1e7);
        let i = ordered();
        let c = cards(&a,&b,&i);
        let pa = point(290.0).with_quantity("water_a",qa,0.2).unwrap();
        let pb = point(310.0).with_quantity("water_b",qb,0.1).unwrap();
        let pi = point(300.0).with_quantity("humidity",moisture(MoistureBasis::RelativeHumidity),0.4).unwrap();
        assert!(c.resolve_at_temperatures([&pa,&pb],&pi,geometry(),envelope()).is_err());
        let law = c.resolve_at_fixed_environments([&pa,&pb],&pi,geometry(),envelope()).unwrap();
        let elastic = law.model().normal_state().elastic();
        close(elastic.surface_a_state().young_modulus_pa(),2e7);
        close(elastic.surface_b_state().young_modulus_pa(),1.5e7);
        assert_eq!(elastic.surface_a_state().resolved().query(),pa);
        assert_eq!(elastic.surface_b_state().resolved().query(),pb);
        assert_eq!(law.model().normal_state().interface_state().query(),pi);
        let reduced = 1.0/(0.9375/2e7+0.9375/1.5e7);
        close(elastic.reduced_modulus_pa(),reduced);
        let k = (4.0/3.0)*0.1_f64.sqrt()*reduced;
        let d = 2e-5_f64;
        let response = law.evaluate("separate-environments",d,0.05,1e-4).unwrap();
        close(response.normal_force_n,k*d.powf(1.5)*(1.0+0.4*0.05));
        close(response.dissipated_power_w,0.4*k*d.powf(1.5)*0.05*0.05);
        let elastic_response = law.evaluate("elastic-patch",d,0.0,1e-4).unwrap();
        close(elastic_response.pressure.peak_pressure_pa,
            1.5*k*d.powf(1.5)/(core::f64::consts::PI*elastic_response.patch_radius_m.powi(2)));
        let changed = c.resolve_at_fixed_environments([&pa.clone().with("water_a",0.0).unwrap(),&pb],
            &pi,geometry(),envelope()).unwrap();
        assert_ne!(law.stiffness_n_per_m_3_2(),changed.stiffness_n_per_m_3_2());
        let interface_changed = c.resolve_at_fixed_environments([&pa,&pb],
            &pi.clone().with("humidity",0.6).unwrap(),geometry(),envelope()).unwrap();
        assert_eq!(law.stiffness_n_per_m_3_2(),interface_changed.stiffness_n_per_m_3_2());
        assert_ne!(law.model().identity(),interface_changed.model().identity());
        assert!(c.resolve_at_fixed_environments([&pa.clone().with("water_a",0.3).unwrap(),&pb],
            &pi,geometry(),envelope()).is_err());
        assert!(c.resolve_at_fixed_environments([&pa.clone().with("slip",0.1).unwrap(),&pb],
            &pi,geometry(),envelope()).is_err());
        assert!(cards(&b,&a,&i).resolve_at_fixed_environments([&pb,&pa],&pi,geometry(),envelope()).is_err());
    }

    #[test]
    fn a_frequency_sample_cannot_become_frozen_contact_elasticity() {
        let q = QuantitySpec::dimensional(fs_qty::Frequency::DIMS);
        let a = material("A","omega",q,4e7,2e7);
        let b = thermal_card("B",2e7,2e7);
        let i = interface(&thermal_card("A",4e7,4e7),&b);
        let pa = point(300.0).with_quantity("omega",q,0.1).unwrap();
        let p = point(300.0).with("omega",0.1).unwrap();
        assert!(matches!(cards(&a,&b,&i).resolve_at_fixed_environments([&pa,&p],&p,geometry(),envelope()),
            Err(ContactMaterialResolveError::Interface(InterfaceBindingError::InvalidDerived {..}))));
    }
}

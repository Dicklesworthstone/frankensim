//! G0/G1/G3: one material card drives the existing Maxwell and heat-feedback laws.
//! All source data are synthetic; these tests do not qualify a real specimen.
use core::convert::Infallible;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, SourceError, ThermoMaxwellSourceQuery};
use fs_material::thermovisco::thermal::{FeedbackConfig, MechanicalControl};
use fs_qty::{Density, Dims, Pressure, Temperature, Time};

const ONE: Dims = Dims([0;6]);
fn limits() -> Limits { Limits { max_strain: 0.04, max_strain_rate_s: 10.0 } }
fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap() }
fn add(claims: &mut ClaimSet, name: &str, dims: Dims, knots: Vec<(f64,f64)>) {
    let lo=knots[0].0; let hi=knots.last().unwrap().0;
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
        value: PropertyValue::Curve { abscissa:"T".into(), abscissa_dims:Temperature::DIMS, knots, dims },
        validity:ValidityDomain::unconstrained().with("T",lo,hi),
        uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::LinearInside,
        observations:vec![], provenance:Provenance { source:"synthetic Maxwell source".into(),
            license:"CC0-1.0".into(), artifact:None } }).unwrap();
}
fn claims(cp:f64) -> ClaimSet {
    let mut claims=ClaimSet::new();
    for (name,dims,value) in [("E_inf",Pressure::DIMS,1e8),("E0",Pressure::DIMS,3e8),
        ("tau0",Time::DIMS,1.0),("scale",ONE,1.0),("free",ONE,0.0),
        ("density",Density::DIMS,1000.0),("cp",Dims([2,0,-2,-1,0,0]),cp)] {
        add(&mut claims,name,dims,vec![(280.0,value),(320.0,value)]);
    }
    add(&mut claims,"shift",ONE,vec![(280.0,1.0),(300.0,1.0),(310.0,0.4),(320.0,0.05)]);
    claims
}
fn card(claims:ClaimSet) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(),phase:"solid".into(),
        process:"source coupling test".into(),revision:0 },claims,vec![]).unwrap()
}
fn resolve(card:&MaterialCard, lower:&QueryPoint, upper:&QueryPoint) -> Result<ResolvedThermoMaxwell,SourceError> {
    ResolvedThermoMaxwell::try_new(card,ThermoMaxwellSourceQuery {
        equilibrium_modulus:&PropertyKey::new("E_inf",Pressure::DIMS),
        branches:&[(PropertyKey::new("E0",Pressure::DIMS),PropertyKey::new("tau0",Time::DIMS))],
        modulus_scale:&PropertyKey::new("scale",ONE),relaxation_scale:&PropertyKey::new("shift",ONE),
        free_strain:&PropertyKey::new("free",ONE),temperature_axis:"T",lower,upper,
        selection:SelectionPolicy::SingleClaimOnly },limits(),|| false)
}
fn caloric(card:&MaterialCard, cp:f64, lower:&QueryPoint, upper:&QueryPoint) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("cp",Dims([2,0,-2,-1,0,0])),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower,upper,
        selection:SelectionPolicy::SingleClaimOnly },SolidLiquidPhase::Solid,-20.0*cp,1e-7).unwrap()
}

#[test]
fn all_source_breakpoints_drive_the_original_history_update() {
    let card=card(claims(5.0));
    let source=resolve(&card,&point(280.0),&point(320.0)).unwrap();
    assert_eq!(source.samples().len(),4);
    let direct=ThermoMaxwell::new(1e8,vec![(3e8,1.0)],vec![
        ThermalKnot {temperature_k:280.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0},
        ThermalKnot {temperature_k:300.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0},
        ThermalKnot {temperature_k:310.0,modulus_scale:1.0,relaxation_scale:0.4,free_strain:0.0},
        ThermalKnot {temperature_k:320.0,modulus_scale:1.0,relaxation_scale:0.05,free_strain:0.0},
    ],limits()).unwrap();
    assert_eq!(source.law(),&direct);
    let mut a=source.law().birth(300.0,0.0).unwrap();
    let mut b=direct.birth(300.0,0.0).unwrap();
    for (t,strain) in [(307.0,0.01),(313.0,0.02),(303.0,0.01),(300.0,0.0)] {
        let sa=a.advance(t,strain,0.2,1e-6).unwrap();
        let sb=b.advance(t,strain,0.2,1e-6).unwrap();
        assert_eq!(sa.point,sb.point); assert_eq!(sa.stress_pa,sb.stress_pa);
        a=sa.point; b=sb.point;
    }
    for sample in source.samples() { assert_eq!(sample.card_identity(),card.content_hash()); }
}

#[test]
fn replacing_only_caloric_source_data_changes_actual_stress_and_self_heating() {
    let mut output=Vec::new();
    for cp in [5.0,10.0] {
        let card=card(claims(cp));
        let source=resolve(&card,&point(280.0),&point(320.0)).unwrap();
        let thermal=caloric(&card,cp,&point(280.0),&point(320.0));
        let old=source.caloric_point(&thermal,300.0,0.0,1e-6,0.001).unwrap();
        let step=old.advance(0.02,0.2,FeedbackConfig { max_iterations:100,
            temperature_tolerance_k:1e-10,energy_tolerance_j:1e-9,relaxation:1.0 },|input| {
            let h=input.initial.specific_enthalpy_j_kg()+input.internal_heat_j/input.mass_kg;
            Ok::<_,Infallible>(UniformEnthalpyStep { state:input.curve.state_at_specific_enthalpy(h).unwrap(),
                external_heat_j:0.0,energy_residual_tolerance_j:1e-10,report:() })
        },|| false).unwrap();
        output.push((step.next.caloric_state().temperature_k(),step.next.mechanical().stress_pa()));
        assert!(step.combined_residual_j.abs()<1e-9);
        assert_eq!(old.time_s(),0.0);
    }
    assert!(output[0].0>output[1].0); assert!(output[0].1<output[1].1);
}

#[test]
fn interior_conflicts_and_hidden_reference_coefficient_changes_are_not_averaged() {
    let mut conflict=claims(5.0);
    add(&mut conflict,"E0",Pressure::DIMS,vec![(305.0,2e8),(315.0,2e8)]);
    assert!(resolve(&card(conflict),&point(280.0),&point(320.0)).is_err());
    // The envelope endpoints alone would falsely appear constant.
    let mut hidden=ClaimSet::new();
    for (name,dims,value) in [("E_inf",Pressure::DIMS,1e8),("tau0",Time::DIMS,1.0),
        ("scale",ONE,1.0),("shift",ONE,1.0),("free",ONE,0.0)] {
        add(&mut hidden,name,dims,vec![(280.0,value),(320.0,value)]);
    }
    add(&mut hidden,"E0",Pressure::DIMS,vec![(280.0,3e8),(309.0,2e8),(320.0,3e8)]);
    assert!(matches!(resolve(&card(hidden),&point(280.0),&point(320.0)),Err(SourceError::Invalid(_))));
}

#[test]
fn foreign_caloric_cards_or_fixed_contexts_cannot_be_bound() {
    let a=card(claims(5.0)); let b=card(claims(10.0));
    let lower=point(280.0).with("p",100000.0).unwrap();
    let upper=point(320.0).with("p",100000.0).unwrap();
    let source=resolve(&a,&lower,&upper).unwrap();
    let same=caloric(&a,5.0,&lower,&upper);
    assert!(source.caloric_point(&same,300.0,0.0,1e-6,0.001).is_ok());
    let other=caloric(&b,10.0,&lower,&upper);
    assert!(source.caloric_point(&other,300.0,0.0,1e-6,0.001).is_err());
    let wrong=caloric(&a,5.0,&point(280.0).with("p",200000.0).unwrap(),
        &point(320.0).with("p",200000.0).unwrap());
    assert!(source.caloric_point(&wrong,300.0,0.0,1e-6,0.001).is_err());
    assert!(resolve(&a,&lower,&point(320.0).with("p",200000.0).unwrap()).is_err());
}

#[test]
fn source_binding_does_not_relax_caloric_physics_or_source_units() {
    let mut c=ClaimSet::new();
    for (name,dims,value) in [("E_inf",Pressure::DIMS,1e8),("E0",Pressure::DIMS,3e8),
        ("tau0",Time::DIMS,1.0),("scale",ONE,1.0),("shift",ONE,1.0),
        ("density",Density::DIMS,1000.0),("cp",Dims([2,0,-2,-1,0,0]),5.0)] {
        add(&mut c,name,dims,vec![(280.0,value),(320.0,value)]);
    }
    add(&mut c,"free",ONE,vec![(280.0,0.0),(320.0,0.001)]);
    let c=card(c); let source=resolve(&c,&point(280.0),&point(320.0)).unwrap();
    let thermal=caloric(&c,5.0,&point(280.0),&point(320.0));
    assert!(source.caloric_point(&thermal,300.0,0.0,1e-6,0.001).is_err());
    assert!(source.law().birth(280.0,0.0).unwrap().advance(320.0,0.0,1.0,1e-6).unwrap().stress_pa<0.0);
    let wrong=PropertyKey::new("E0",Time::DIMS);
    assert!(ResolvedThermoMaxwell::try_new(&c,ThermoMaxwellSourceQuery {
        equilibrium_modulus:&wrong,branches:&[],modulus_scale:&PropertyKey::new("scale",ONE),
        relaxation_scale:&PropertyKey::new("shift",ONE),free_strain:&PropertyKey::new("free",ONE),
        temperature_axis:"T",lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly,
    },limits(),|| false).is_err());
}

#[test]
fn cancellation_refuses_before_publication_and_zero_branches_still_need_declared_laws() {
    let c=card(claims(5.0)); let mut calls=0;
    let query=ThermoMaxwellSourceQuery { equilibrium_modulus:&PropertyKey::new("E_inf",Pressure::DIMS),
        branches:&[],modulus_scale:&PropertyKey::new("scale",ONE),relaxation_scale:&PropertyKey::new("shift",ONE),
        free_strain:&PropertyKey::new("free",ONE),temperature_axis:"T",
        lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly };
    assert!(matches!(ResolvedThermoMaxwell::try_new(&c,query,limits(),|| {calls+=1;calls==5}),Err(SourceError::Cancelled)));
    let source=ResolvedThermoMaxwell::try_new(&c,query,limits(),|| false).unwrap();
    assert_eq!(source.law().branch_count(),0);
    let thermal=caloric(&c,5.0,&point(280.0),&point(320.0));
    let state=source.caloric_point(&thermal,300.0,0.0,1e-6,0.001).unwrap();
    let control=MechanicalControl::Stress { stress_pa:2e6,tolerance_pa:1e-6 };
    let step=state.advance_controlled(control,0.2,FeedbackConfig {max_iterations:5,
        temperature_tolerance_k:1e-8,energy_tolerance_j:1e-9,relaxation:1.0},|input| {
        assert_eq!(input.internal_heat_j,0.0);
        Ok::<_,Infallible>(UniformEnthalpyStep {state:input.initial,external_heat_j:0.0,
            energy_residual_tolerance_j:1e-10,report:()})
    },|| false).unwrap();
    assert_eq!(step.next.caloric_state().temperature_k(),300.0);
}

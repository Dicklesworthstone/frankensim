use super::*;
use crate::section::thermal::AxialConstraint;
use crate::state_point::{INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, ScalarAdmissibility, ScalarPropertyRequirement};
use fs_matdb::{MaterialStateId, PropertyClaim, PropertyValue, InterpolationPolicy, Provenance, UncertaintyModel};
use fs_evidence::ValidityDomain;

fn expanding_card(alpha: f64) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
    let (base,mut req)=card(1.0);let mut claims=base.claims().clone();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ALPHA,INVERSE_TEMPERATURE_DIMS),
        value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:fs_qty::Temperature::DIMS,
            knots:vec![(300.0,alpha),(400.0,2.0*alpha)],dims:INVERSE_TEMPERATURE_DIMS},
        validity:ValidityDomain::unconstrained().with("T",300.0,400.0).with("pressure",100000.0,100000.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"synthetic heated constraint".into(),license:"test-only".into(),artifact:None},
    }).unwrap();
    req.push(ScalarPropertyRequirement::try_new(ALPHA,INVERSE_TEMPERATURE_DIMS,ScalarAdmissibility::Finite).unwrap());
    (MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),process:"test".into(),revision:0},claims,vec![]).unwrap(),req)
}
fn body<'a>(card:&'a MaterialCard,req:&[ScalarPropertyRequirement],chart:&'a ResolvedHeatCapacityCurve)
    ->HeatedConcentricSection<'a> {
    let law=ConcentricSectionTemperatureLaw::try_new(section(card,req),&[card,card],&["T","T"]).unwrap();
    HeatedConcentricSection::try_new(law,&[chart,chart],1.0,
        &[boundary(0.0),boundary(1.0)],vec![interface()],budget()).unwrap()
}

#[test]
fn fixed_and_free_supports_follow_the_same_actual_heat_history() {
    let (card,req)=expanding_card(1e-5);let chart=caloric(&card);
    let original=body(&card,&req,&chart);
    let fixed=original.clone().with_axial_constraint(AxialConstraint::FixedStrain(0.001),0.01).unwrap();
    let free=original.clone().with_axial_constraint(AxialConstraint::FixedForce(0.0),0.01).unwrap();
    let (a,ra)=fixed.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
    let (b,rb)=free.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
    let (plain,rp)=original.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
    assert_eq!(a.thermal().cells(),b.thermal().cells());
    assert_eq!(a.thermal().cells(),plain.thermal().cells());
    assert_eq!(ra.thermal,rp.thermal);assert_eq!(ra.thermal,rb.thermal);
    assert_eq!(ra.electrical_energy_j,rp.electrical_energy_j);
    let held=a.axial_state().unwrap();let released=b.axial_state().unwrap();
    assert_eq!(held.section().section(),a.section());
    assert!(held.response().force_n()<fixed.axial_state().unwrap().response().force_n());
    near(released.response().force_n(),0.0);
    assert!(released.response().total_strain()>0.0);
    assert!(released.response().layers().iter().any(|l|l.stress_pa()!=0.0));
    assert_eq!(original.thermal().time_s(),0.0);
    assert!(original.axial_state().is_none());
}

#[test]
fn changing_only_expansion_data_changes_force_not_heat_or_electrical_work() {
    let run=|alpha| {
        let (card,req)=expanding_card(alpha);let chart=caloric(&card);
        let initial=body(&card,&req,&chart).with_axial_constraint(AxialConstraint::FixedStrain(0.001),0.01).unwrap();
        let (next,report)=initial.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
        (next.thermal().temperature_k(CellId(0)).unwrap(),report.electrical_energy_j,
            next.axial_state().unwrap().response().force_n())
    };
    let a=run(1e-5);let b=run(2e-5);
    assert_eq!(a.0,b.0);assert_eq!(a.1,b.1);assert!(b.2<a.2);
}

#[test]
fn late_strain_refusal_and_cancellation_do_not_publish_thermal_or_axial_history() {
    let (card,req)=expanding_card(0.01);let chart=caloric(&card);
    let initial=body(&card,&req,&chart).with_axial_constraint(AxialConstraint::FixedStrain(0.0),1e-4).unwrap();
    let old_cells=initial.thermal().cells().clone();let old_axial=initial.axial_state().unwrap().clone();
    assert!(initial.advance_to(0.01,OhmicDrive::Voltage(0.2),config(),||Ok(())).is_ok());
    let mut checks=0;
    assert!(matches!(initial.advance_to(1.0,OhmicDrive::Voltage(0.2),config(),||{checks+=1;Ok(())}),
        Err(SectionHeatingError::Material(_))));
    assert!(checks>20);
    assert_eq!(initial.thermal().cells(),&old_cells);assert_eq!(initial.axial_state(),Some(&old_axial));
    let mut checks=0;
    assert!(matches!(initial.advance_to(0.02,OhmicDrive::Voltage(0.2),config(),|| {
        checks+=1;if checks>20 {Err(DepositionError::Cancelled)} else {Ok(())}
    }),Err(SectionHeatingError::Thermal(DepositionError::Cancelled))));
    assert_eq!(initial.thermal().time_s(),0.0);assert_eq!(initial.axial_state(),Some(&old_axial));
}

#[test]
fn attaching_after_heating_preserves_original_stress_free_reference_and_restart() {
    let (card,req)=expanding_card(1e-5);let chart=caloric(&card);
    let initial=body(&card,&req,&chart);
    let (hot,_)=initial.advance_to(0.1,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap();
    let attached=hot.clone().with_axial_constraint(AxialConstraint::FixedStrain(0.0),0.01).unwrap();
    assert!(attached.axial_state().unwrap().response().force_n()<0.0);
    assert_eq!(attached.thermal().cells(),hot.thermal().cells());
    let a=attached.advance_to(0.2,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap().0;
    let b=attached.clone().advance_to(0.2,OhmicDrive::Voltage(0.2),config(),||Ok(())).unwrap().0;
    assert_eq!(a.axial_state(),b.axial_state());assert_eq!(a.thermal().cells(),b.thermal().cells());
    let (plain,req)=card_without_expansion();let chart=caloric(&plain);
    assert!(body(&plain,&req,&chart).with_axial_constraint(AxialConstraint::FixedStrain(0.0),0.01).is_err());
}
fn card_without_expansion()->(MaterialCard,Vec<ScalarPropertyRequirement>) {card(1.0)}

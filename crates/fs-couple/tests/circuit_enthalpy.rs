//! Actual shared circuit and thermal owners, not copied RC or heat steppers.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::*;
use fs_thermal::enthalpy::{CellId,StepConfig};

#[test]
fn finite_capacitor_heat_is_partitioned_then_conducted_between_material_bodies() {
    let (a,ca)=material(1.0,1.0,0.0); let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let thermal=initial(&curves,2.0,None);
    let old=circuit(&thermal,&cards,&curves);
    let proposed=step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||false).unwrap();
    // Independently solved series RC midpoint and the two-cell BE linear system.
    let q1=10.0*(1.0-0.125)/(1.0+0.125);
    let current=0.5*(10.0+q1)/4.0;
    let qa=current*current; let qb=3.0*qa;
    near(proposed.circuit().charges_c()[0],q1);
    near(proposed.thermal().temperature_k(CellId(30)).unwrap(),300.0+(4.0*qa+2.0*qb)/8.0);
    near(proposed.thermal().temperature_k(CellId(10)).unwrap(),300.0+(2.0*qa+3.0*qb)/8.0);
    let r=proposed.report();near(r.electrical.resistors[0].heat_j,qa);near(r.electrical.resistors[1].heat_j,qb);
    near(r.deposited_joule_j,qa+qb);near(r.non_electrical_heat_j,0.0);
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    for (i,id) in [CellId(30),CellId(10)].into_iter().enumerate() {
        near(r.source[i].temperature_k(),proposed.thermal().temperature_k(id).unwrap());
        assert_eq!(r.electrical.resistors[i].material.conductor.material().property("electrical_resistivity").unwrap().value_si(),[1.0,3.0][i]);
    }
    assert_eq!(old.charges_c(),&[10.0]);assert_eq!(thermal.time_s(),0.0);
}

#[test]
fn changing_one_resistor_material_changes_both_thermal_bodies_and_finite_supply() {
    let run=|second_r| {
        let (a,ca)=material(1.0,1.0,0.0);let (b,cb)=material(second_r,2.0,0.0);
        let cards=[a,b];let curves=[ca,cb];let thermal=initial(&curves,0.0,None);
        let old=circuit(&thermal,&cards,&curves);
        let next=step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
            |s|coefficients(s,&cards,&curves),||false).unwrap();
        (next.thermal().temperature_k(CellId(30)).unwrap(),next.thermal().temperature_k(CellId(10)).unwrap(),
            next.circuit().charges_c()[0],next.report().deposited_joule_j)
    };
    let a=run(3.0);let b=run(6.0);
    assert!(b.0<a.0 && b.1<a.1 && b.2>a.2 && b.3<a.3);
}

#[test]
fn actual_endpoint_resistance_feedback_and_restart_keep_all_storage() {
    let (a,ca)=material(1.0,1.0,0.005);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let mut thermal=initial(&curves,1.0,None);
    let mut electrical=circuit(&thermal,&cards,&curves);let energy0=thermal.energy_j()+50.0;
    for k in 1..=8 {
        let solve=||step_candidate(&thermal,k as f64*0.25,input(&electrical),config(),&[0.0;2],
            |s|coefficients(s,&cards,&curves),||false).unwrap();
        let next=solve();let replay=solve();
        assert_eq!(next.circuit(),replay.circuit());assert_eq!(next.thermal().cells(),replay.thermal().cells());
        let r=next.report();let t=next.thermal().temperature_k(CellId(30)).unwrap();
        near(r.electrical.resistors[0].material.conductor.resistance_ohm(),1.0+0.005*(t-300.0));
        near(next.thermal().energy_j()+r.electrical.stored_after_j,energy0);
        let (t,c,_)=next.into_parts();thermal=t;electrical=c;
    }
}

#[test]
fn independent_sources_remain_distinct_from_internal_electrical_heat() {
    let (a,ca)=material(1.0,1.0,0.0);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let thermal=initial(&curves,2.0,None);let old=circuit(&thermal,&cards,&curves);
    let candidate=step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],|s| {
        let mut c=coefficients(s,&cards,&curves)?;
        c.sources_w.insert(CellId(30),2.0);c.sources_w.insert(CellId(10),-1.0);
        Ok::<_,fs_couple::cavity_phs::piston::thermal::power::joule::JouleHeatingError>(c)
    },||false).unwrap();
    near(candidate.report().non_electrical_heat_j,1.0);
    near(candidate.thermal().energy_j()-thermal.energy_j()+candidate.report().electrical.stored_after_j-50.0,1.0);
    assert_eq!(candidate.thermal().cells()[&CellId(30)].boundary,thermal.cells()[&CellId(30)].boundary);
}

#[test]
fn complete_mapping_final_budget_and_late_cancellation_are_joint_refusals() {
    let (a,ca)=material(1.0,1.0,0.0);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let thermal=initial(&curves,2.0,None);let old=circuit(&thermal,&cards,&curves);
    let saved=thermal.cells().clone();let electrical=old.clone();
    let duplicate=[DESTINATIONS[0];2];let inactive=[DESTINATIONS[0],ResistorCell {branch:2,cell:CellId(99)}];
    for map in [&DESTINATIONS[..1],&duplicate[..],&inactive[..]] {
        assert!(step_candidate(&thermal,1.0,CircuitEnthalpyInput {destinations:map,..input(&old)},config(),&[0.0;2],
            |s|coefficients(s,&cards,&curves),||false).is_err());
    }
    let mut polls=0;
    let normal=step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||{polls+=1;false}).unwrap();
    let limited=CircuitEnthalpySolve {maximum_circuit_trials:normal.report().circuit_trials-1,..config()};
    assert!(matches!(step_candidate(&thermal,1.0,input(&old),limited,&[0.0;2],
        |s|coefficients(s,&cards,&curves),||false),Err(CircuitEnthalpyError::CircuitBudget)));
    let mut actual=0;
    assert!(matches!(step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||{actual+=1;actual==polls}),Err(CircuitEnthalpyError::Cancelled)));
    assert!(matches!(step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
        |_|Err::<CircuitThermalCoefficients<()>,_>("missing resistor source"),||false),
        Err(CircuitEnthalpyError::Source("missing resistor source"))));
    assert_eq!(thermal.cells(),&saved);assert_eq!(old,electrical);
}

#[test]
fn thermal_events_are_not_implicitly_consumed_by_circuit_trials() {
    let (a,ca)=material(1.0,1.0,0.0);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let thermal=initial(&curves,1.0,Some(1.0));let old=circuit(&thermal,&cards,&curves);
    assert!(step_candidate(&thermal,2.0,input(&old),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||false).is_err());
    let left=step_candidate(&thermal,1.0,input(&old),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||false).unwrap();
    assert_eq!(left.thermal().cells()[&CellId(30)].boundary.source_w,0.0);
    assert!(step_candidate(left.thermal(),2.0,input(left.circuit()),config(),&[0.0;2],
        |s|coefficients(s,&cards,&curves),||false).is_err());
    let right=left.thermal().advance_to(1.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
    assert_eq!(right.cells()[&CellId(30)].boundary.source_w,4.0);
    assert_eq!(right.energy_j(),left.thermal().energy_j());
}

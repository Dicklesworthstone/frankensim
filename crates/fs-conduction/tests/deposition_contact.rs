//! Real shared thermal residual, independently derived finite-contact balance.
#[path = "contact_material_state/common.rs"]
mod common;
use common::*;
use std::collections::BTreeMap;
use fs_conduction::transient::deposition::*;
use fs_conduction::transient::deposition::contact::*;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_matdb::{ClaimSelection, ClaimSet};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};

#[path = "deposition_contact/observer.rs"]
mod observer;

fn initial<'a>(curves: &'a [EquilibriumEnthalpyPhaseCurve;2], area:f64) -> DepositionSimulation<'a> {
    let mut events=Vec::new();
    for (i,(id,h)) in [(CellId(30),25.0),(CellId(10),-50.0)].into_iter().enumerate() {
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id,segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
            specific_h_j_kg:h,geometry:CellGeometry {start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],
                volume_m3:1.0,surface_area_m2:6.0},boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:CellId(30),b:CellId(10),area_m2:area,conductance_w_k:999.0})});
    DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
        (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:10})
        .unwrap().advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0
}
fn sides(curves:&[EquilibriumEnthalpyPhaseCurve;2]) -> [ContactCell<'_>;2] {
    [0,1].map(|i|ContactCell {cell:CellId(if i==0 {30}else{10}),material_id:MaterialId(i as u64),curve:&curves[i]})
}

#[test]
fn phase_contact_is_in_series_and_feedback_closes_the_actual_implicit_cell_balance() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let old=initial(&curves,1.0); let (source,id)=source(&cards,1.0); let context=point();
    let contacts=[MaterialContactLink {link:7,sides:sides(&curves),point:binding(&cards,&source,&context),selection:ClaimSelection::Pinned(id)}];
    // Independent balance: q*(1/2 + 1 + (25-q)/100) = 0.1*(50-q).
    let q=10.0/(1.85+(1.85_f64*1.85-0.2).sqrt());
    let step=BackwardEulerStep::new(&old,0.1,||Ok(())).unwrap();
    assert_eq!(step.cell_ids(),&[CellId(10),CellId(30)]);
    let trial=step.evaluate(&[q,-q],|state| {
        let mut g=BTreeMap::from([(7,2.0)]);
        let samples=apply_material_contacts(state,&contacts,&mut g,||false).map_err(|_|DepositionError::MaterialDomain)?;
        Ok::<_,DepositionError>(StepCoefficients {sources_w:BTreeMap::new(),conductances_w_k:g,payload:samples})
    },||Ok(())).unwrap();
    assert!(trial.maximum_cell_residual_j()<1e-10);
    let (next,samples,audit)=trial.accept(1e-9,1e-9,||Ok(())).unwrap();
    assert!((next.energy_j()-old.energy_j()).abs()<1e-9);
    assert_eq!(next.temperature_k(CellId(30)).unwrap(),400.0);
    assert!((next.temperature_k(CellId(10)).unwrap()-(350.0+q)).abs()<1e-10);
    assert!((samples[0].interface.value_m2_k_per_w()-(1.25-q/100.0)).abs()<1e-12);
    assert!((samples[0].conductance_w_k-1.0/(1.75-q/100.0)).abs()<1e-12);
    assert_eq!(samples[0].series.as_ref().unwrap().terms().len(),2);
    assert!(audit.external_heat_j.abs()<1e-10);
    assert!((q-5.0/1.85).abs()>0.01,"lagging the contact's phase must not satisfy this case");
    assert_eq!(old.cells()[&CellId(30)].specific_h_j_kg,25.0);
}

#[test]
fn isolation_still_checks_contact_data_and_late_cancellation_preserves_all_coefficients() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let old=initial(&curves,1.0); let (source,id)=source(&cards,1.0); let context=point();
    let mut contact=MaterialContactLink {link:7,sides:sides(&curves),point:binding(&cards,&source,&context),selection:ClaimSelection::Pinned(id)};
    let saved=BTreeMap::from([(7,2.0),(99,123.0)]); let mut g=saved.clone(); let mut polls=0;
    assert!(apply_material_contacts(&old,&[contact],&mut g,||{polls+=1;polls==3}).is_err());
    assert_eq!(g,saved);
    let zero=initial(&curves,0.0);
    let samples=apply_material_contacts(&zero,&[contact],&mut g,||false).unwrap();
    assert_eq!(g[&7],0.0); assert_eq!(g[&99],123.0); assert!(samples[0].series.is_none());
    assert_eq!(samples[0].interface.receipt().selected,id);
    let absent=interface(&cards,ClaimSet::new()); contact.point.interface=&absent;
    let snapshot=g.clone(); assert!(apply_material_contacts(&zero,&[contact],&mut g,||false).is_err());
    assert_eq!(g,snapshot);
}

#[test]
fn mismatched_caloric_instances_order_and_missing_bulk_inputs_refuse() {
    let cards=[bulk("A"),bulk("B")]; let curves=[curve(&cards[0]),curve(&cards[1])];
    let old=initial(&curves,1.0); let (source,id)=source(&cards,1.0); let context=point();
    let mut contact=MaterialContactLink {link:7,sides:sides(&curves),point:binding(&cards,&source,&context),selection:ClaimSelection::Pinned(id)};
    let copied=curves[0].clone(); contact.sides[0].curve=&copied;
    let saved=BTreeMap::from([(7,2.0)]); let mut g=saved.clone();
    assert!(apply_material_contacts(&old,&[contact],&mut g,||false).is_err()); assert_eq!(g,saved);
    contact.sides=sides(&curves); contact.sides.swap(0,1);
    assert!(apply_material_contacts(&old,&[contact],&mut g,||false).is_err()); assert_eq!(g,saved);
    contact.sides=sides(&curves);
    assert!(apply_material_contacts(&old,&[contact,contact],&mut g,||false).is_err()); assert_eq!(g,saved);
    assert!(apply_material_contacts(&old,&[contact],&mut BTreeMap::new(),||false).is_err());
    g.insert(7,0.0); let report=apply_material_contacts(&old,&[contact],&mut g,||false).unwrap();
    assert_eq!(report[0].conductance_w_k,0.0); assert_eq!(g[&7],0.0);
}

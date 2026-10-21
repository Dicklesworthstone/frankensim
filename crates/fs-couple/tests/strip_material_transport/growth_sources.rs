//! Source continuity through real thermal births and subsequent mechanics.
use super::*;
use fs_couple::circuit_enthalpy::material::MaterialNetworkError;
use fs_couple::strip_enthalpy::{LayerThermalBinding, growth::{GrowingStrip, LayerBirth}};

#[test]
fn growth_retains_old_claims_when_hot_layers_gain_a_competing_source() {
    let m = [5.0,50.0].map(|cp| material_with_conductivity_alternative(cp,3e8,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],Some(vec![(300.01,0.02),(320.0,0.02)])));
    let mut incoming = cell(CellId(99),2); incoming.material_id=MaterialId(0);
    let (heat,strip) = initial(&m,9.0,vec![
        TimedEvent {time_s:0.4,event:Event::Birth(incoming)},
        TimedEvent {time_s:0.4,event:Event::Contact(ThermalLink {
            id:6,a:CELLS[1],b:CellId(99),area_m2:0.001,conductance_w_k:9.0})}]);
    let binding = LayerThermalBinding {layer:99,cell:CellId(99),
        material:&m[0].mechanical,caloric:&m[0].caloric};
    let (old,_) = GrowingStrip::try_new(heat,strip,1,bindings(&m).to_vec(),
        vec![LayerBirth {binding,geometry:geometry(2)}],3,1e-9,||false).unwrap();
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    let before=prepared.resolve(old.thermal(),||false).unwrap().source;
    let (born,_) = old.step_transport_candidate(&prepared,0.4,control(),config(),&[0.0;2],||false).unwrap();
    assert!(born.thermal().temperature_k(CELLS[0]).unwrap()>300.01);
    let mut expanded=cells(&m,&context,&key).into_iter().collect::<Vec<_>>();
    expanded.push(MaterialCell {cell:CellId(99),material_id:MaterialId(0),curve:m[0].caloric.curve(),
        point:ThermalMaterialPoint {material:&m[0].card,context:&context,
            temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let links=[LINKS[0],ConductionLink {link:6,distance_a_m:0.0005,distance_b_m:0.0005}];
    let successor=MaterialNetwork {cells:&expanded,resistors:&[],links:&links};
    // Re-selection is ambiguous, but the previously accepted source is still valid.
    assert!(successor.prepare(born.thermal(),||false).is_err());
    let continued=prepared.extend_thermal(&successor,born.thermal(),||false).unwrap();
    let sample=continued.resolve(born.thermal(),||false).unwrap().source;
    for old_sample in &before.cells {
        let now=sample.cells.iter().find(|s|s.cell==old_sample.cell).unwrap();
        let claim=|s:&fs_couple::circuit_enthalpy::material::CellMaterialSample|
            s.conductivity.as_ref().unwrap().property("conductivity").unwrap().answer().receipt.selected;
        assert_eq!(claim(now),claim(old_sample));
    }
    let stretch=Control {axial:AxialControl::Strain(0.015),bending:BendingControl::CurvaturePerM(0.0)};
    let (next,r)=born.step_transport_candidate(&continued,0.6,stretch,config(),&[0.0;3],||false).unwrap();
    assert!(r.interval.layer_energy_j[&99].dissipated>0.0);
    for link in &r.interval.source.links {
        for k in link.conductivity_w_m_k { near(k,0.002,1e-14); }
        near(link.conductance_w_k,0.002,1e-14);
    }
    assert_eq!(next.strip().layers()[&99].birth_time_s,0.4);
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    assert_eq!(old.thermal().time_s(),0.0);
}

#[test]
fn an_isolated_old_cell_can_admit_conductivity_when_a_real_contact_appears() {
    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,_)=initial(&m,9.0,vec![
        TimedEvent {time_s:0.0,event:Event::Separation(5)},
        TimedEvent {time_s:0.4,event:Event::Contact(ThermalLink {
            id:5,a:CELLS[0],b:CELLS[1],area_m2:0.001,conductance_w_k:9.0})}]);
    let context=point(300.0); let key=key(); let mut isolated=cells(&m,&context,&key);
    for c in &mut isolated {c.conductivity=None;}
    let network=MaterialNetwork {cells:&isolated,resistors:&[],links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();
    let (right,_)=heat.advance_to(0.4,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap();
    let bound=cells(&m,&context,&key);
    let successor=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let next=prepared.extend_thermal(&successor,&right,||false).unwrap();
    let sample=next.resolve(&right,||false).unwrap().source;
    assert_eq!(sample.links.len(),1);
    near(sample.links[0].conductance_w_k,0.002,1e-14);
    assert!(sample.cells.iter().all(|c|c.conductivity.is_some()));
    assert!(heat.links().is_empty());
}

#[test]
fn extension_refuses_source_retargeting_and_circuit_state_loss() {
    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,_)=initial(&m,0.002,vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let prepared=network.prepare(&heat,||false).unwrap();
    let other_context=point(300.0).with("lot",8.0).unwrap();
    let changed=cells(&m,&other_context,&key);
    let retargeted=MaterialNetwork {cells:&changed,resistors:&[],links:&LINKS};
    assert!(matches!(prepared.extend_thermal(&retargeted,&heat,||false),Err(MaterialNetworkError::Invalid(_))));
    let mut missing=cells(&m,&context,&key); missing[0].conductivity=None;
    let removed=MaterialNetwork {cells:&missing,resistors:&[],links:&LINKS};
    assert!(matches!(prepared.extend_thermal(&removed,&heat,||false),Err(MaterialNetworkError::Invalid(_))));
    let resistors=[CellResistor {destination:ResistorCell {branch:0,cell:CELLS[0]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}];
    let electrical=MaterialNetwork {cells:&bound,resistors:&resistors,links:&LINKS};
    assert!(prepared.extend_thermal(&electrical,&heat,||false).is_err());
    let electrical_prepared=electrical.prepare(&heat,||false).unwrap();
    assert!(electrical_prepared.extend_thermal(&network,&heat,||false).is_err());
    assert_eq!(heat.time_s(),0.0);
}

#[test]
fn extension_never_replaces_an_expired_pin_with_an_eligible_alternative() {
    let m=[5.0,50.0].map(|cp|material_with_conductivity_alternative(cp,3e8,0.05,0.0,
        vec![(280.0,0.002),(300.05,0.002)],Some(vec![(300.01,0.02),(320.0,0.02)])));
    let (heat,_)=initial(&m,9.0,vec![
        TimedEvent {time_s:0.0,event:Event::Separation(5)},
        TimedEvent {time_s:0.0,event:Event::Boundary(CELLS[0],Boundary {source_w:0.01,..boundary()})}]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();
    // Deliberately warm with the thermal owner alone, outside the pinned model.
    let (warm,_)=heat.advance_to(0.1,StepConfig {max_dt_s:0.05,energy_tolerance_j:1e-9},||Ok(())).unwrap();
    assert!(warm.temperature_k(CELLS[0]).unwrap()>300.05);
    assert!(network.prepare(&warm,||false).is_ok(),"a different claim is now eligible");
    assert!(prepared.extend_thermal(&network,&warm,||false).is_err());
    assert!(prepared.resolve(&heat,||false).is_ok(),"old accepted source choices remain unchanged");
    assert_eq!(heat.time_s(),0.0);
}

#[test]
fn late_cancelled_extension_preserves_old_choices_and_replays() {
    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,_)=initial(&m,0.002,vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let prepared=network.prepare(&heat,||false).unwrap();
    let before=prepared.resolve(&heat,||false).unwrap().source;
    let mut reordered=cells(&m,&context,&key); reordered.swap(0,1);
    let successor=MaterialNetwork {cells:&reordered,resistors:&[],links:&LINKS};
    let mut polls=0;
    let next=prepared.extend_thermal(&successor,&heat,||{polls+=1;false}).unwrap();
    assert_eq!(next.resolve(&heat,||false).unwrap().source,before);
    let mut current=0;
    assert!(matches!(prepared.extend_thermal(&successor,&heat,||{current+=1;current==polls}),
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))));
    assert_eq!(prepared.resolve(&heat,||false).unwrap().source,before);
    assert_eq!(prepared.extend_thermal(&successor,&heat,||false).unwrap().resolve(&heat,||false).unwrap().source,before);
}

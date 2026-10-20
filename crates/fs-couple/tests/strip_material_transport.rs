//! G0/G3 composition: original material, conduction, section and nonlinear owners.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::ResistorCell;
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::strip_enthalpy::{self, StripEnthalpyError};
use fs_couple::strip_enthalpy::transport::{self, TransportError};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{PropertyKey, QueryPoint};
use fs_solid::deposition::{AxialControl, BendingControl, Control};
use fs_thermal::enthalpy::*;

const LINKS: [ConductionLink;1] = [ConductionLink {
    link:5,distance_a_m:0.0005,distance_b_m:0.0005,
}];
fn control() -> Control {
    Control {axial:AxialControl::Strain(0.01),bending:BendingControl::MomentNm(0.0)}
}
fn cells<'a>(m: &'a [Material;2], context: &'a QueryPoint, key: &'a PropertyKey) -> [MaterialCell<'a>;2] {
    [0,1].map(|i|MaterialCell {cell:CELLS[i],material_id:MaterialId(i as u64),
        curve:m[i].caloric.curve(),point:ThermalMaterialPoint {
            material:&m[i].card,context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))})
}
fn key() -> PropertyKey { PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS) }
fn materials(k: Vec<(f64,f64)>) -> [Material;2] {
    [material_with_conductivity(5.0,3e8,0.05,0.0,k.clone()),
     material_with_conductivity(50.0,3e8,0.05,0.0,k)]
}

#[test]
fn endpoint_conductivity_drives_the_original_residual_and_retains_its_source() {
    let m=materials(vec![(280.0,0.0002),(300.0,0.002),(320.0,0.04)]);
    // This deliberately wrong authored coefficient must be replaced, not added.
    let (heat,strip)=initial(&m,9.0,vec![]);
    let context=point(300.0); let key=key(); let cells=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let pinned=network.prepare(&heat,||false).unwrap(); let map=bindings(&m);
    let run=||transport::step_candidate(&heat,&strip,&map,&pinned,0.2,control(),config(),&[0.0;2],||false).unwrap();
    let next=run(); let replay=run();
    let manual=strip_enthalpy::step_candidate(&heat,&strip,&map,0.2,control(),config(),&[0.0;2],
        |t|transport::coefficients(&pinned,t,&map,||false),||false).unwrap();
    assert_eq!(next.strip(),manual.strip()); assert_eq!(next.strip(),replay.strip());
    assert_eq!(next.thermal().cells(),manual.thermal().cells());
    let r=next.report(); let link=&r.source.links[0];
    for i in 0..2 {
        let t=next.thermal().temperature_k(CELLS[i]).unwrap();
        near(link.conductivity_w_m_k[i],0.002+0.0019*(t-300.0),1e-12);
        let sample=r.source.cells.iter().find(|s|s.cell==CELLS[i]).unwrap();
        near(sample.phase.temperature_k(),t,1e-12);
        let state=sample.conductivity.as_ref().unwrap();
        assert_eq!(state.card_identity(),m[i].card.content_hash());
        m[i].mechanical.admit_state_point(state).unwrap();
    }
    near(link.conductance_w_k,0.001/(0.0005/link.conductivity_w_m_k[0]
        +0.0005/link.conductivity_w_m_k[1]),1e-12);
    near(r.non_mechanical_heat_j,0.0,1e-12);
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    assert_eq!(next.thermal().links()[&5].conductance_w_k,9.0,"trial coefficient is not a persistent input");
    assert_eq!(heat.time_s(),0.0);
}

#[test]
fn changing_only_material_conductivity_changes_temperatures_and_compatible_curvature() {
    let run=|k| {
        let m=materials(vec![(280.0,k),(320.0,k)]);
        let (heat,strip)=initial(&m,9.0,vec![]);
        let context=point(300.0); let key=key(); let cells=cells(&m,&context,&key);
        let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
        let pinned=network.prepare(&heat,||false).unwrap();
        let next=transport::step_candidate(&heat,&strip,&bindings(&m),&pinned,0.2,
            control(),config(),&[0.0;2],||false).unwrap();
        (next.thermal().temperature_k(CELLS[0]).unwrap()-next.thermal().temperature_k(CELLS[1]).unwrap(),
            next.strip().curvature_per_m())
    };
    let low=run(0.0002); let high=run(0.02);
    assert!(low.0>0.1 && high.0<low.0);
    assert!(low.1.abs()>1e-4 && high.1.abs()<low.1.abs());
}

#[test]
fn same_card_does_not_allow_different_fixed_source_conditions() {
    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,strip)=initial(&m,0.002,vec![]);
    let context=point(300.0).with("lot",7.0).unwrap(); let key=key(); let cells=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let pinned=network.prepare(&heat,||false).unwrap();
    assert!(matches!(transport::step_candidate(&heat,&strip,&bindings(&m),&pinned,0.2,
        control(),config(),&[0.0;2],||false),
        Err(StripEnthalpyError::Source(TransportError::Context {..}))));
    assert_eq!(heat.time_s(),0.0); assert_eq!(strip.axial_strain(),0.0);
}

#[test]
fn pinned_source_domain_and_late_cancellation_preserve_both_histories() {
    let narrow=materials(vec![(280.0,0.002),(300.05,0.002)]);
    let (heat,strip)=initial(&narrow,0.002,vec![]);
    let before=heat.cells().clone(); let old=strip.clone();
    let context=point(300.0); let key=key(); let bound=cells(&narrow,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let pinned=network.prepare(&heat,||false).unwrap();
    assert!(transport::step_candidate(&heat,&strip,&bindings(&narrow),&pinned,0.2,
        control(),config(),&[0.0;2],||false).is_err());
    assert_eq!(heat.cells(),&before); assert_eq!(strip,old);

    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,strip)=initial(&m,0.002,vec![]); let before=heat.cells().clone(); let old=strip.clone();
    let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let pinned=network.prepare(&heat,||false).unwrap(); let map=bindings(&m);
    let mut polls=0;
    transport::step_candidate(&heat,&strip,&map,&pinned,0.2,control(),config(),&[0.0;2],
        ||{polls+=1;false}).unwrap();
    let mut current=0;
    assert!(matches!(transport::step_candidate(&heat,&strip,&map,&pinned,0.2,control(),config(),&[0.0;2],
        ||{current+=1;current==polls}),Err(StripEnthalpyError::Cancelled)));
    assert_eq!(heat.cells(),&before); assert_eq!(strip,old);
}

#[test]
fn an_electrical_branch_cannot_be_silently_discarded_by_the_thermal_adapter() {
    let m=materials(vec![(280.0,0.002),(320.0,0.002)]);
    let (heat,strip)=initial(&m,0.002,vec![]);
    let context=point(300.0); let key=key(); let cells=cells(&m,&context,&key);
    let resistors=[CellResistor {destination:ResistorCell {branch:0,cell:CELLS[0]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}];
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&LINKS};
    let pinned=network.prepare(&heat,||false).unwrap();
    assert!(matches!(transport::step_candidate(&heat,&strip,&bindings(&m),&pinned,0.2,
        control(),config(),&[0.0;2],||false),
        Err(StripEnthalpyError::Source(TransportError::Invalid(_)))));
}

#[test]
fn material_transport_continues_after_atomic_birth_and_an_explicit_topology_update() {
    use fs_couple::strip_enthalpy::{LayerThermalBinding, growth::{GrowingStrip, LayerBirth}};
    let m=materials(vec![(280.0,0.0002),(300.0,0.002),(320.0,0.04)]);
    let mut newborn=cell(CellId(99),2); newborn.material_id=MaterialId(0);
    let later=vec![TimedEvent {time_s:0.4,event:Event::Birth(newborn)},
        TimedEvent {time_s:0.4,event:Event::Contact(ThermalLink {
            id:6,a:CELLS[1],b:CellId(99),area_m2:0.001,conductance_w_k:9.0})}];
    let (heat,strip)=initial(&m,9.0,later);
    let newborn_binding=LayerThermalBinding {layer:99,cell:CellId(99),
        material:&m[0].mechanical,caloric:&m[0].caloric};
    let (old,_)=GrowingStrip::try_new(heat,strip,1,bindings(&m).to_vec(),
        vec![LayerBirth {binding:newborn_binding,geometry:geometry(2)}],3,1e-9,||false).unwrap();
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&bound,resistors:&[],links:&LINKS};
    let pinned=network.prepare(old.thermal(),||false).unwrap();
    let (born,r)=old.step_transport_candidate(&pinned,0.4,control(),config(),&[0.0;2],||false).unwrap();
    assert_eq!(r.born_layers,1); assert_eq!(r.events.events,2);
    assert_eq!(r.interval.source.cells.len(),2);
    assert_eq!(r.interval.source.links.len(),1);
    assert_eq!(born.thermal().cells().len(),3);
    assert_eq!(born.thermal().links().len(),2);
    assert!(born.step_transport_candidate(&pinned,0.6,control(),config(),&[0.0;3],||false).is_err(),
        "an obsolete transport topology must refuse, not silently ignore the new cell");
    let mut next_cells=cells(&m,&context,&key).into_iter().collect::<Vec<_>>();
    next_cells.push(MaterialCell {cell:CellId(99),material_id:MaterialId(0),
        curve:m[0].caloric.curve(),point:ThermalMaterialPoint {
            material:&m[0].card,context:&context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let links=[LINKS[0],ConductionLink {link:6,distance_a_m:0.0005,distance_b_m:0.0005}];
    let next_network=MaterialNetwork {cells:&next_cells,resistors:&[],links:&links};
    let next_pinned=next_network.prepare(born.thermal(),||false).unwrap();
    let stretch=Control {axial:AxialControl::Strain(0.015),bending:BendingControl::CurvaturePerM(0.0)};
    let (next,r)=born.step_transport_candidate(&next_pinned,0.6,stretch,config(),&[0.0;3],||false).unwrap();
    assert!(r.interval.layer_energy_j[&99].dissipated>0.0);
    assert_eq!(r.interval.source.links.len(),2);
    assert!(r.interval.source.links.iter().all(|s|s.conductance_w_k<0.1));
    assert_eq!(next.strip().layers()[&99].birth_time_s,0.4);
    assert!(!next.birth_equilibrium_pending());
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    assert_eq!(old.thermal().time_s(),0.0);
}

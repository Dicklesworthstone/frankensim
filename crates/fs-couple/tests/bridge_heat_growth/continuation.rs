//! G0/G3: growth and its actual post-event radiative network are one candidate.
use super::*;
use fs_conduction::transient::deposition::radiation::{MaterialRadiation, RadiatingCell};
use fs_couple::circuit_enthalpy::contacts::heat::{HeatNetworkError, MaterialHeatNetwork, PreparedHeatNetwork};
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::material_circuit::MaterialCircuit;
use fs_exec::{CancelGate, Cx};
use fs_thermal::enthalpy::{ThermalLink, TimedEvent};

fn with_successor<R>(small_surface: bool,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &GrowingBridge<'_,'_>, &MaterialCircuit,
        &PreparedHeatNetwork<'_,'_>, &MaterialHeatNetwork<'_>) -> R) -> R {
    let mut incoming=arrival(77,99,LENGTH,0.02,1);
    if small_surface {
        if let Event::Birth(cell)=&mut incoming.event { cell.geometry.surface_area_m2=2e-5; }
    }
    let contact=TimedEvent {time_s:DT,event:Event::Contact(ThermalLink {
        id:6,a:CELLS[0],b:CellId(77),area_m2:AREA,conductance_w_k:99.0})};
    with_case(0.7,0.5,true,true,vec![incoming,contact],|m,cx,gate,n,t,b,c| {
        let prepared=n.prepare(cx,t).unwrap();
        let (old,_)=GrowingBridge::try_new(t.clone(),b.clone(),1,map(m,true),
            vec![future(&m[1],77,99,0.02)],2,mechanics(),1e-9,||false).unwrap();
        let mut cells=n.electrical.bulk.cells.iter().map(|s|MaterialCell {
            cell:s.cell,material_id:s.material_id,curve:s.curve,point:s.point,
            conductivity:s.conductivity.clone()}).collect::<Vec<_>>();
        let source=&n.electrical.bulk.cells[1];
        cells.push(MaterialCell {cell:CellId(77),material_id:source.material_id,
            curve:source.curve,point:source.point,conductivity:source.conductivity.clone()});
        let mut links=n.electrical.bulk.links.to_vec();
        links.push(ConductionLink {link:6,distance_a_m:0.05,distance_b_m:0.01});
        let bulk=MaterialNetwork {cells:&cells,resistors:n.electrical.bulk.resistors,links:&links};
        let electrical=ElectricalContactNetwork {bulk:&bulk,contacts:n.electrical.contacts,
            electrical_bases:n.electrical.electrical_bases};
        let old_rad=n.radiation.unwrap();
        // Explicitly retire the passive receiver's old face and introduce a face
        // of the newly deposited material. The limiting enclosure is authored,
        // not inferred from the bent beam, and the passive body still holds heat.
        let surfaces=old_rad.cells.iter().map(|s|RadiatingCell {
            surface:s.surface,cell:if s.surface==1 {CellId(77)}else{s.cell},
            material_id:s.material_id,curve:s.curve,material:s.material,
            selection:s.selection.clone()}).collect::<Vec<_>>();
        let radiation=MaterialRadiation {enclosure:old_rad.enclosure,cells:&surfaces,
            power_tolerance_w:old_rad.power_tolerance_w,radiosity_tolerance_w_m2:old_rad.radiosity_tolerance_w_m2};
        let successor=MaterialHeatNetwork {electrical:&electrical,
            thermal_contacts:n.thermal_contacts,radiation:Some(&radiation)};
        run(cx,gate,&old,c,&prepared,&successor)
    })
}

#[test]
fn returned_successor_heats_new_material_then_advances_the_original_beam_and_supply() {
    with_successor(false,|cx,_,old,c,p,successor| {
        let (manual,manual_c,_) = old.step_heat_network_candidate(cx,c,p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let manual_p=p.extend_fixed_circuit(cx,successor,manual.thermal(),&[Some(0),None]).unwrap();
        let (born,next_c,ready,r) = old.step_heat_network_with_successor_candidate(
            cx,c,p,successor,&[Some(0),None],&[],DT,1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        assert_eq!(born.thermal().cells(),manual.thermal().cells());
        assert_eq!(born.thermal().links(),manual.thermal().links());
        assert_eq!(born.bridge(),manual.bridge()); assert_eq!(next_c,manual_c);
        let expected=manual_p.resolve(cx,manual.thermal()).unwrap();
        let actual=ready.resolve(cx,born.thermal()).unwrap();
        assert_eq!(actual.bulk.source.radiation,expected.bulk.source.radiation);
        assert_eq!(actual.bulk.source.thermal_contacts,expected.bulk.source.thermal_contacts);
        assert_eq!(actual.bulk.source.radiation.as_ref().unwrap().surfaces()[1].cell,CellId(77));
        assert_eq!(r.births.thermal.events,2); assert_eq!(r.births.thermal.substeps,0);
        assert!(r.birth_equilibrium_pending); assert_eq!(next_c.resistors().len(),1);
        assert!(p.resolve(cx,born.thermal()).is_err());
        let (after,final_c,second)=born.step_heat_network_candidate(cx,&next_c,&ready,&[],2.0*DT,
            1.0,&[0.02;2],controls(),&[0.0;3]).unwrap();
        assert!(second.interval.material.radiation.as_ref().unwrap().sources_w()[&CellId(77)]>0.0);
        assert!(after.thermal().temperature_k(CellId(77)).unwrap()>300.0);
        let tip=LENGTH+0.02;
        assert!((after.bridge().deformed_position_m(tip,0.0).unwrap()[1]
            -born.bridge().deformed_position_m(tip,0.0).unwrap()[1]).abs()>1e-10);
        assert!(final_c.charges_c()[0]<next_c.charges_c()[0]);
        assert!(!second.birth_equilibrium_pending);
        assert!(second.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
        assert_eq!(old.thermal().time_s(),0.0); assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn impossible_new_surface_rolls_back_a_successful_birth_and_finite_circuit_debit() {
    with_successor(true,|cx,gate,old,c,p,successor| {
        let (born,debit,_) = old.step_heat_network_candidate(cx,c,p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        assert_eq!(born.remaining_births(),0); assert_eq!(born.bridge().cells().len(),2);
        assert!(debit.charges_c()[0]<CHARGE,"the physical interval really used stored energy");
        let saved=old.clone(); let original=c.clone();
        assert!(matches!(old.step_heat_network_with_successor_candidate(
            cx,c,p,successor,&[Some(0),None],&[],DT,1.0,&[0.02],controls(),&[0.0;2]),
            Err(ElectricalGrowthError::Successor(HeatNetworkError::Radiation(_)))));
        assert_eq!(old.thermal().cells(),saved.thermal().cells());
        assert_eq!(old.thermal().links(),saved.thermal().links());
        assert_eq!(old.bridge(),saved.bridge()); assert_eq!(c,&original);
        assert_eq!(old.remaining_births(),1); assert!(p.resolve(cx,old.thermal()).is_ok());
        gate.request();
        assert!(matches!(old.step_heat_network_with_successor_candidate(
            cx,c,p,successor,&[Some(0),None],&[],DT,1.0,&[0.02],controls(),&[0.0;2]),
            Err(ElectricalGrowthError::Interval(electrical::ElectricalError::Cancelled))));
        assert_eq!(c,&original); assert_eq!(old.thermal().time_s(),0.0);
    });
}

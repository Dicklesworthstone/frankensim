//! Same original fixture: real chemistry, contacts, radiation, charge and motion.
use super::*;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::Contact;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::NetworkError;
use fs_couple::cavity_phs::piston::thermal::ThermalPistonError;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::{ChamberDrive, MaterialExchangeState};
use std::collections::BTreeMap;

fn surfaces(patches: &[MaterialContact]) -> Vec<GasSurface> {
    patches.iter().map(|p| GasSurface { contact: p.contact, area_m2: p.geometry.area_m2 }).collect()
}

#[test]
fn generic_heat_law_reproduces_original_conduction_and_radiation_exactly() {
    fixture([0.2,0.8],150.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let gases=gas_fixture::pistons();
        let chambers=[gas_fixture::chamber(10,&gases[0],0.00005),gas_fixture::chamber(20,&gases[1],0.00005)];
        let patches=gas_fixture::patches();let geometry=surfaces(&patches);
        let original=prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&chambers,
            &patches,controls(),&[0.0;2]).unwrap();
        let generic=prepared.step_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],
            &chambers,&geometry,controls(),&[0.0;2],|contact,state,source,gas| {
                let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
                prepared.resolve_conductive_gas_contact(patch,state,source,gas)
            }).unwrap();
        assert_eq!(generic.thermal().cells(),original.thermal().cells());
        assert_eq!(generic.circuit(),original.circuit());
        assert_eq!(generic.step(),original.step());
        for (a,b) in generic.chambers().iter().zip(original.chambers()) {
            assert_eq!(a.boundary,b.boundary);
        }
        assert_eq!(generic.material().material.radiation.as_ref().unwrap().net_heat_w(),
            original.material().material.radiation.as_ref().unwrap().net_heat_w());
    }));
}

#[test]
fn typed_sources_follow_one_chemical_stage_through_the_complete_feedback() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let law=reaction(50.0);let gases=[piston(&law),gas_fixture::pistons()[1].clone()];
        let chambers=[gas_fixture::chamber(10,&gases[0],0.00005),gas_fixture::chamber(20,&gases[1],0.00005)];
        let patches=gas_fixture::patches();let geometry=surfaces(&patches);
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        let expected=prepared.reactive_step_candidate(cx,&thermal,0.00005,&circuit,&[],
            &chambers,&reactions,&patches,controls(),&[0.0;2]).unwrap();
        let calls=std::cell::Cell::new(0usize);
        let result=prepared.reactive_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],
            &chambers,&reactions,&geometry,controls(),&[0.0;2],|contact,state,source,gas| {
                calls.set(calls.get()+1);
                if contact.chamber==10 { assert_eq!(gas.components().len(),2); }
                let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
                let sample=prepared.resolve_conductive_gas_contact(patch,state,source,gas)?;
                let phase=source.bulk.cells.iter().find(|c|c.cell==contact.cell).unwrap().phase;
                Ok::<_,MaterialNetworkError>(ConductanceSample {
                    conductance_w_k:sample.conductance_w_k,source:(sample.source,gas.clone(),phase),
                })
            }).unwrap();
        assert!(calls.get()>0);
        assert_eq!(result.thermal().cells(),expected.thermal().cells());
        assert_eq!(result.circuit(),expected.circuit());
        assert_eq!(result.reactions(),expected.reactions());
        assert_eq!(result.step(),expected.step());
        let chemistry=law.parcel_candidate(gases[0].gas(),0.00005,chemical(),||false).unwrap();
        assert_eq!(&result.reactions()[0].reaction,chemistry.step());
        assert!(result.report().electrical.contacts[0].heat_j>0.0);
        assert!(result.report().source.material.material.radiation.is_some());
        for (a,b) in result.chambers().iter().zip(expected.chambers()) {
            assert_eq!(a.piston,b.piston);
            assert_eq!(a.boundary,b.boundary);
            for surface in &a.surfaces {
                let (contact,(material,gas,phase))=&surface.source;
                assert_eq!(gas,&a.boundary.piston.midpoint);
                let cell=result.report().source.material.material.bulk.cells.iter()
                    .find(|c|c.cell==contact.cell).unwrap();
                assert_eq!(*phase,cell.phase);
                assert_eq!(Some(material),cell.conductivity.as_ref());
            }
        }
    }));
}

#[test]
fn all_laws_share_aggregate_exposure_and_zero_area_cannot_heat() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let gases=gas_fixture::pistons();
        let chambers=[gas_fixture::chamber(10,&gases[0],0.00005),gas_fixture::chamber(20,&gases[1],0.00005)];
        let surface=GasSurface {contact:Contact {id:1,chamber:10,cell:CellId(30)},area_m2:0.8};
        let other=GasSurface {contact:Contact {id:2,chamber:20,..surface.contact},..surface};
        let callback=|_:Contact,_:&DepositionSimulation<'_>,_:&MaterialExchangeSample,_:&fs_material::thermochemical_gas::ResolvedThermochemicalGasState|
            Ok::<_,&'static str>(ConductanceSample {conductance_w_k:0.3,source:()});
        assert!(prepared.step_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],&chambers,
            &[surface,other],controls(),&[0.0;2],callback).is_err());
        assert!(prepared.step_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],&chambers,
            &[GasSurface {area_m2:0.0,..surface}],controls(),&[0.0;2],callback).is_err());
        assert!(prepared.step_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],&chambers,
            &[surface,surface],controls(),&[0.0;2],callback).is_err());
        assert_eq!(thermal.time_s(),0.0);
    }));
}

#[test]
fn source_error_and_late_cancellation_do_not_publish_chemical_or_electrical_state() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit| with_cx(|cx,gate| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(50.0);let gas=piston(&law);
        let old=gas.clone();let cells=thermal.cells().clone();let electrical=circuit.clone();
        let chambers=[gas_fixture::chamber(10,&gas,0.00005)];
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        let geometry=surfaces(&gas_fixture::patches()[..1]);
        let failure=prepared.reactive_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],
            &chambers,&reactions,&geometry,controls(),&[0.0;2],|_,_,_,_| {
                Err::<ConductanceSample<()>,_>("deliberate original contact refusal")
            });
        assert!(matches!(failure,Err(ReactiveExchangeError::Exchange(CircuitEnthalpyError::Source(
            NetworkError::Boundary {error:ThermalPistonError::Interface("deliberate original contact refusal"),..})))));
        let cancelled=prepared.reactive_with_gas_law_candidate(cx,&thermal,0.00005,&circuit,&[],
            &chambers,&reactions,&geometry,controls(),&[0.0;2],|_,_,_,_| {
                gate.request();Ok::<_,&'static str>(ConductanceSample {conductance_w_k:0.0,source:()})
            });
        assert!(matches!(cancelled,Err(ReactiveExchangeError::Cancelled)));
        assert_eq!(gas,old);assert_eq!(thermal.cells(),&cells);assert_eq!(circuit,electrical);
        assert_eq!(thermal.time_s(),0.0);
    }));
}


fn live_drive(id: u64, gas: &GasPiston) -> ChamberDrive {
    let chamber = gas_fixture::chamber(id, gas, 0.00005);
    ChamberDrive {
        chamber: id,
        heat_j: chamber.input.heat_j,
        exterior_pressure_pa: chamber.input.exterior_pressure_pa,
        external_force_n: chamber.input.external_force_n,
        piston_solve: chamber.piston_solve,
        exchange: chamber.exchange,
    }
}

#[test]
fn live_state_preview_and_advance_publish_the_typed_heat_law_once() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let law=reaction(50.0);
        let mut state=MaterialExchangeState::from_parts(
            thermal,circuit,BTreeMap::from([
                (10,piston(&law)),
                (20,gas_fixture::pistons()[1].clone()),
            ]));
        let drives=[
            live_drive(10,&state.chambers()[&10]),
            live_drive(20,&state.chambers()[&20]),
        ];
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        let patches=gas_fixture::patches();
        let geometry=surfaces(&patches);
        let chemistry=law.parcel_candidate(
            state.chambers()[&10].gas(),0.00005,chemical(),||false).unwrap();
        let resolve=|contact:Contact,state:&DepositionSimulation<'_>,
            source:&MaterialExchangeSample,
            gas:&fs_material::thermochemical_gas::ResolvedThermochemicalGasState| {
            let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
            let sample=prepared.resolve_conductive_gas_contact(patch,state,source,gas)?;
            Ok::<_,MaterialNetworkError>(ConductanceSample {
                conductance_w_k:sample.conductance_w_k,
                source:(sample.source,gas.clone()),
            })
        };
        let (preview,preview_report)=state.step_with_gas_law_candidate(
            &prepared,cx,0.00005,&[],&drives,&reactions,&geometry,
            controls(),&[0.0;2],resolve).unwrap();
        assert_eq!(state.time_s(),0.0);
        assert_eq!(preview_report.report.source.reactions[0].reaction,*chemistry.step());

        let report=state.advance_with_gas_law(
            &prepared,cx,0.00005,&[],&drives,&reactions,&geometry,
            controls(),&[0.0;2],|contact,state,source,gas| {
                let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
                let sample=prepared.resolve_conductive_gas_contact(patch,state,source,gas)?;
                Ok::<_,MaterialNetworkError>(ConductanceSample {
                    conductance_w_k:sample.conductance_w_k,
                    source:(sample.source,gas.clone()),
                })
            }).unwrap();
        assert_eq!(state.time_s(),0.00005);
        assert_eq!(state.thermal().cells(),preview.thermal().cells());
        assert_eq!(state.circuit(),preview.circuit());
        assert_eq!(state.chambers(),preview.chambers());
        assert_eq!(report.step,preview_report.step);
        assert_eq!(report.report.source.reactions,preview_report.report.source.reactions);
        for chamber in &report.report.source.chambers {
            let accepted=&state.chambers()[&chamber.id];
            assert_eq!(accepted.last_window_duration_s(),Some(0.00005));
            assert_eq!(accepted.last_window_transfers(),preview.chambers()[&chamber.id].last_window_transfers());
            for surface in &chamber.surfaces {
                assert_eq!(&surface.source.1,&chamber.boundary.piston.midpoint);
            }
        }
    }));
}

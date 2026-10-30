//! Additional physical links retain their original two-sided caloric ownership.
use super::*;
use crate::circuit_enthalpy::contacts::{ContactThermalCoefficients, heat::{MaterialHeatNetwork, HeatNetworkError},
    material::ElectricalContactNetwork};
use crate::circuit_enthalpy::CircuitThermalCoefficients;
use fs_exec::Cx;

const LINK:u64=79;
fn linked<'a>(source:&MaterialNetwork<'a>, old:&DepositionSimulation<'_>)->DepositionSimulation<'a> {
    let mut events:Vec<_>=old.cells().values().map(|cell|TimedEvent {time_s:0.0,event:Event::Birth(cell.clone())}).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:LINK,a:CellId(0),b:CellId(1),area_m2:0.001,conductance_w_k:0.0})});
    let laws:Vec<_>=source.cells.iter().map(|s|(s.material_id,s.curve as &dyn ThermalLaw)).collect();
    DepositionSimulation::new(events,&laws,DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:128})
        .unwrap().advance_to(0.0,StepConfig {max_dt_s:1.0/256.0,energy_tolerance_j:1e-11},||Ok(())).unwrap().0
}

#[test]
fn fresh_reciprocal_link_heat_reaches_both_cells_without_replacing_electrical_or_material_history() {
    fixture_powered(1.0,false,|m,plain,powered,heat,bases| {
        let mut heat=linked(plain,&heat);let gate=CancelGate::new_clock_free();
        let links=BTreeMap::from([(LINK,1e-6)]);
        let (_,material,mechanical)=native(m);
        assert!(FrameRelaxationThermalSystem::new(material,mechanical,&heat,plain,&DEST,budget(),&gate).is_err());
        let (placement,material,mechanical)=native(m);
        let physical=placement.physical_state(mechanical.components(),&gate).unwrap();
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(
            vec![(material,DEST.to_vec())],mechanical,&heat,plain,&links,budget(),&gate).unwrap();
        let initial=powered.prepare_with_link_conductances(&heat,bases,&links,||false).unwrap()
            .resolve_with_link_conductances(&heat,&links,||false).unwrap();
        let circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
            (1,0,Branch::Capacitor {farads:CAPACITANCE}), (1,2,Branch::Resistor {ohms:999.0}),
            (2,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},initial.resistors,vec![],vec![INITIAL_CHARGE],&[]).unwrap();
        // Ordinary binding must not silently discard an independent link.
        assert!(owner.bind_material_circuit(&heat,powered,circuit.clone(),bases,&gate).is_err());
        assert!(owner.electrical().is_none());
        owner.bind_material_circuit_with_link_conductances(&heat,powered,circuit.clone(),bases,
            vec![],&links,&gate).unwrap();
        let external=vec![0.0;owner.mechanical().mode_count()];let saved=heat.cells().clone();
        let input=RelaxationElectricalInput {ports:&[],solve:electrical_controls()};
        assert!(owner.step_electrical(&mut heat,&external,input,controls(),&[0.0;2],&gate).is_err());
        assert_eq!(heat.cells(),&saved);assert_eq!(owner.electrical().unwrap(),&circuit);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        assert_eq!(owner.mechanical().samples_rendered(),0);
        let electrical=ElectricalContactNetwork {bulk:powered,contacts:&[],electrical_bases:bases};
        let network=MaterialHeatNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:13,kernel_id:52,tile:0,iteration:0},
                fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
            let prepared=network.prepare_with_link_conductances(&cx,&heat,&links).unwrap();
            for _ in 0..4 {
                let old=heat.cells().clone();
                let report=owner.step_heat_network_with(&mut heat,&external,RelaxationHeatInput {network:&prepared,cx:&cx},
                    Some(input),controls(),&[0.0;2],&gate,|trial| {
                        let t=trial.temperature_k(CellId(0)).unwrap();
                        let g=1e-6*(1.0+0.01*(t-300.0));
                        let ContactThermalCoefficients {bulk,contacts}=prepared.resolve_with_link_conductances(
                            &cx,trial,&BTreeMap::from([(LINK,g)]))?;
                        Ok::<_,HeatNetworkError>(ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
                            resistors:bulk.resistors,sources_w:bulk.sources_w,conductances_w_k:bulk.conductances_w_k,
                            source:(bulk.source,g)}})
                    }).unwrap();
                let g=report.transport_source.unwrap();
                let exchange=g*(heat.temperature_k(CellId(0)).unwrap()-heat.temperature_k(CellId(1)).unwrap())/256.0;
                for i in 0..2 {
                    let id=CellId(i as u64);let cell=&heat.cells()[&id];
                    let material=report.layers.iter().find(|l|l.destination.cell==id).unwrap().heat_j;
                    let joule=report.electrical.as_ref().unwrap().resistors.iter().find(|r|r.material.branch==i+1).unwrap().heat_j;
                    near(cell.mass_kg*(cell.specific_h_j_kg-old[&id].specific_h_j_kg),
                        material+joule+if i==0 {-exchange}else{exchange});
                }
                near(report.external_thermal_heat_j,0.0);near(report.total_energy_residual_j,0.0);
            }
        });
    });
}

#[test]
fn invalid_external_link_coverage_refuses_and_empty_maps_keep_the_original_path() {
    fixture_powered(1.0,false,|m,plain,powered,heat,bases| {
        let gate=CancelGate::new_clock_free();let linked=linked(plain,&heat);
        for bad in [BTreeMap::new(),BTreeMap::from([(LINK,-1.0)]),BTreeMap::from([(LINK,f64::NAN)]),
            BTreeMap::from([(LINK,1e-6),(999,0.0)])] {
            let (_,material,mechanical)=native(m);
            assert!(FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                mechanical,&linked,plain,&bad,budget(),&gate).is_err());
        }
        let mut accepted=Vec::new();
        for explicit in [false,true] {
            let mut heat=heat.clone();let (_,material,mechanical)=native(m);
            let mut owner=if explicit {
                FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                    mechanical,&heat,plain,&BTreeMap::new(),budget(),&gate)
            } else {FrameRelaxationThermalSystem::new(material,mechanical,&heat,plain,&DEST,budget(),&gate)}.unwrap();
            let c=circuit(powered,&heat,bases,None);
            if explicit {owner.bind_material_circuit_with_link_conductances(&heat,powered,c,bases,vec![],&BTreeMap::new(),&gate).unwrap();}
            else {owner.bind_material_circuit(&heat,powered,c,bases,&gate).unwrap();}
            let force=vec![0.0;owner.mechanical().mode_count()];
            owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                controls(),&[0.0;2],&gate).unwrap();
            accepted.push((heat.cells().clone(),owner.electrical().unwrap().clone(),owner.mechanical().material_memory().to_vec()));
        }
        assert_eq!(accepted[0],accepted[1]);
    });
}

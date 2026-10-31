//! Endpoint flow and heat are one evaluation, including discarded Newton trials.
use super::*;
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::transport::{
    RelaxationFlowCoefficients, RelaxationHeatInput,
};
use fs_couple::circuit_enthalpy::{CircuitThermalCoefficients, contacts::ContactThermalCoefficients};

#[test]
fn resolved_constant_flow_is_identical_to_borrowed_prescribed_flow() {
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let mut outcomes=Vec::new();
        for dynamic in [false,true] {
            let mut heat=heat.clone();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
            let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let graph=flow(0.2,true,false,4000.0);let force=vec![0.0;owner.mechanical().mode_count()];
            for _ in 0..3 {
                let resolve=|trial:&DepositionSimulation<'_>| {
                    let ContactThermalCoefficients {bulk,contacts}=network.resolve_with_links(cx,trial,&links(),2,
                        |_,_,t|law.at_temperature(t)).unwrap();
                    Ok::<_, &'static str>(ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
                        resistors:bulk.resistors,sources_w:bulk.sources_w,conductances_w_k:bulk.conductances_w_k,
                        source:(bulk.source.heat,0.2),
                    }})
                };
                let heat_input=RelaxationHeatInput {network:prepared,cx};
                let report=if dynamic {
                    owner.step_heat_network_with_resolved_advection(&mut heat,&force,heat_input,Some(input()),settings(4),
                        &[0.0;4],&gate,|trial|resolve(trial).map(|heat|RelaxationFlowCoefficients {heat,flow:graph.clone()}))
                } else {
                    owner.step_heat_network_with_advection(&mut heat,&force,heat_input,Some(input()),settings(4),
                        &[0.0;4],&gate,&graph,resolve)
                }.unwrap();
                assert_eq!(report.transport_source,Some(0.2));
            }
            assert_eq!(graph.audit(),Default::default());
            outcomes.push((heat.cells().clone(),owner.electrical().unwrap().clone(),
                owner.mechanical().material_memory().to_vec(),placement.physical_state(owner.mechanical().components(),&gate).unwrap()));
        }
        assert_eq!(outcomes[0],outcomes[1]);
    });
}

#[test]
fn final_flow_and_source_receipt_use_the_actual_accepted_caloric_endpoint() {
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (_,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit,
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
        let force=vec![0.0;owner.mechanical().mode_count()];let calls=std::cell::Cell::new(0);
        let final_query=std::cell::Cell::new(0.0);
        let resolve=|trial:&DepositionSimulation<'_>| {
            let temperature=trial.temperature_k(COOLANT).unwrap();
            let mass=0.2*(1.0+0.02*(temperature-300.0));
            calls.set(calls.get()+1);final_query.set(temperature);
            let ContactThermalCoefficients {bulk,contacts}=network.resolve_with_links(cx,trial,&links(),2,
                |_,_,t|law.at_temperature(t)).unwrap();
            Ok::<_, &'static str>(RelaxationFlowCoefficients {
                flow:flow(mass,true,false,4000.0),
                heat:ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
                    resistors:bulk.resistors,sources_w:bulk.sources_w,conductances_w_k:bulk.conductances_w_k,
                    source:(bulk.source.heat,(temperature,mass)),
                }},
            })
        };
        for _ in 0..3 {
            let before=heat.cells().clone();
            let report=owner.step_heat_network_with_resolved_advection(&mut heat,&force,
                RelaxationHeatInput {network:prepared,cx},Some(input()),settings(4),&[0.0;4],&gate,&resolve).unwrap();
            let (temperature,mass)=report.transport_source.unwrap();
            assert_eq!(temperature.to_bits(),heat.temperature_k(COOLANT).unwrap().to_bits());
            assert_eq!(temperature.to_bits(),final_query.get().to_bits());
            let other=heat.temperature_k(CellId(8)).unwrap();
            let transfer=mass*4000.0*(temperature-other)*DT;
            let cell=&heat.cells()[&CellId(8)];
            near(cell.mass_kg*(cell.specific_h_j_kg-before[&CellId(8)].specific_h_j_kg),transfer,3e-8);
            assert_eq!(report.enthalpy_flow.as_ref().unwrap().external_enthalpy_j,0.0);
        }
        assert!(calls.get()>3);
    });
}

#[test]
fn invalid_late_graph_or_typed_source_failure_rolls_back_all_histories_and_retries() {
    #[derive(Debug,PartialEq)] enum Failure { Late }
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
        let force=vec![0.0;owner.mechanical().mode_count()];
        let cells=heat.cells().clone();let memory=owner.mechanical().material_memory().to_vec();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let source=owner.material().members()[0].source.clone();
        for kind in 0..3 {
            let calls=std::cell::Cell::new(0);
            let result=owner.step_heat_network_with_resolved_advection(&mut heat,&force,
                RelaxationHeatInput {network:prepared,cx},Some(input()),settings(4),&[0.0;4],&gate,|trial| {
                    let call=calls.get();calls.set(call+1);
                    if call>0 && kind==0 {return Err(Failure::Late);}
                    let graph=if call==0 {flow(0.2,true,false,4000.0)} else if kind==1 {
                        // A graph may be admissible at the start but become
                        // forbidden only inside a rejected thermal stencil.
                        Advected::new((),vec![MassFlow {id:1,a:CellId(0),b:COOLANT,kg_s:0.2},
                            MassFlow {id:2,a:COOLANT,b:CellId(0),kg_s:0.2}],2).unwrap()
                    } else {Advected::new((),vec![MassFlow {id:1,a:COOLANT,b:CellId(8),kg_s:0.2}],1).unwrap()};
                    let ContactThermalCoefficients {bulk,contacts}=network.resolve_with_links(cx,trial,&links(),2,
                        |_,_,t|law.at_temperature(t)).unwrap();
                    Ok(RelaxationFlowCoefficients {flow:graph,heat:ContactThermalCoefficients {contacts,
                        bulk:CircuitThermalCoefficients {resistors:bulk.resistors,sources_w:bulk.sources_w,
                            conductances_w_k:bulk.conductances_w_k,source:(bulk.source.heat,())}}})
                });
            if kind==0 {assert!(matches!(result,Err(RelaxationTransportError::Source(Failure::Late))));}
            else {assert!(result.is_err());}
            assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);
            assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(owner.mechanical().samples_rendered(),0);
            assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,source);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        }
        network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),settings(4),
            &[0.0;4],&gate,&links(),2,|_,_,t|law.at_temperature(t),&flow(0.2,true,false,4000.0)).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

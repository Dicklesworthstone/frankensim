//! Original circulation changes the complete thermal/circuit/material endpoint.
use super::*;
use fs_thermal::enthalpy::advection::{Advected,MassFlow,MassBoundary};

fn settings(cells:usize)->FrameThermalSolve {FrameThermalSolve {maximum_cells:cells,..solve()}}
fn flow(rate:f64,closed:bool,reverse:bool,cp:f64)->Advected<()> {
    if closed {
        let sign=if reverse {-1.0}else{1.0};
        let pair=if reverse {[CellId(8),COOLANT]}else{[COOLANT,CellId(8)]};
        Advected::new((),vec![MassFlow {id:1,a:pair[0],b:pair[1],kg_s:sign*rate},
            MassFlow {id:2,a:pair[1],b:pair[0],kg_s:sign*rate}],2).unwrap()
    } else {
        Advected::new((),vec![],2).unwrap().with_boundaries(vec![
            MassBoundary::Inlet {id:1,cell:COOLANT,material_id:MaterialId(7),kg_s:rate,specific_h_j_kg:5.0*cp},
            MassBoundary::Outlet {id:2,cell:COOLANT,kg_s:rate}],2).unwrap()
    }
}

#[test]
fn closed_and_open_flow_share_local_balances_and_change_continuing_electrical_and_mechanical_response() {
    let mut endpoints=Vec::new();
    for (closed,rate,reverse,powered) in [(true,0.0,false,true),(true,0.2,false,true),
        (true,0.8,false,true),(true,0.2,true,true),(false,0.2,false,true),(true,0.2,false,false)] {
        let s=Spec {circulating:closed,powered,..spec()};
        fixture(s,|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            if powered {owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();}
            let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let graph=flow(rate,closed,reverse,s.cp);let guesses=vec![0.0;heat.cells().len()];
            let controls=settings(guesses.len());
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            let energy0=owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                +if powered {CHARGE*CHARGE/(2.0*CAP)}else{0.0};
            let (mut supply,mut outside)=(0.0,0.0);
            for step in 1..=12 {
                let before=heat.cells().clone();
                let report=network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,
                    powered.then_some(input()),controls,&guesses,&gate,&links(),2,|_,_,t|law.at_temperature(t),&graph).unwrap();
                let t=heat.temperature_k(COOLANT).unwrap();
                let other=if closed {heat.temperature_k(CellId(8)).unwrap()}else{295.0};
                let carried=rate*s.cp*(other-t);
                let advected=report.enthalpy_flow.as_ref().unwrap();
                near(advected.external_enthalpy_j,if closed {0.0}else{DT*carried},1e-9);
                near(report.external_thermal_heat_j,advected.external_enthalpy_j,1e-8);
                let mut received=0.0;
                for sample in report.transport_source.as_ref().unwrap().finite.samples() {
                    let i=sample.wall_cell.0 as usize;let wall=heat.temperature_k(sample.wall_cell).unwrap();
                    let power=PATCHES[i]*h(s,wall,t)*(wall-t);received+=DT*power;
                    near(sample.heat_into_fluid_w,power,1e-8);
                    let joule=report.electrical.as_ref().map_or(0.0,|e|e.resistors[i].heat_j);
                    let loss=report.layers.iter().find(|l|l.destination.cell==sample.wall_cell).unwrap().heat_j;
                    let cell=&heat.cells()[&sample.wall_cell];
                    near(cell.mass_kg*(cell.specific_h_j_kg-before[&sample.wall_cell].specific_h_j_kg),loss+joule-DT*power,2e-8);
                }
                let cell=&heat.cells()[&COOLANT];
                near(cell.mass_kg*(cell.specific_h_j_kg-before[&COOLANT].specific_h_j_kg),received+DT*carried,3e-8);
                if closed {
                    let cell=&heat.cells()[&CellId(8)];
                    near(cell.mass_kg*(cell.specific_h_j_kg-before[&CellId(8)].specific_h_j_kg),-DT*carried,3e-8);
                    near(advected.internal_at_cell_entropy_j_k,DT*carried*(1.0/t-1.0/other),1e-10);
                }
                let electrical=report.electrical.as_ref();
                supply+=report.mechanical.external_work_j+report.external_thermal_heat_j
                    +electrical.map_or(0.0,|e|e.native.supplied);
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                    +electrical.map_or(0.0,|e|e.stored_after_j)-energy0+outside-supply,0.0,4e-7);
                assert!(report.entropy.accounting_residual_j_k.abs()<1e-8);
                assert_eq!(heat.time_s(),step as f64*DT);
                for (id,cell) in heat.cells() {assert_eq!(cell.mass_kg,before[id].mass_kg);}
            }
            near(owner.mechanical().material_memory()[0],0.007*(-12.0*DT/0.7).exp(),1e-10);
            assert_eq!(graph.audit(),Default::default());
            endpoints.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(CellId(1)).unwrap(),heat.temperature_k(COOLANT).unwrap(),
                owner.electrical().map_or(0.0,|e|e.charges_c()[0])]);
        });
    }
    for i in [1,2] {for j in 0..4 {assert!((endpoints[i][j]-endpoints[0][j]).abs()>1e-12);}}
    for j in 0..4 {near(endpoints[1][j],endpoints[3][j],1e-10);}
}

#[test]
fn empty_graph_has_identical_finite_coolant_trajectory() {
    fixture(spec(),|m,plain,electric,circuit,heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let mut outcomes=Vec::new();
        for advected in [false,true] {
            let mut heat=heat.clone();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
            let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let graph=Advected::new((),vec![],0).unwrap();let force=vec![0.0;owner.mechanical().mode_count()];
            for _ in 0..3 {
                if advected {network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),
                    solve(),&[0.0;3],&gate,&links(),2,|_,_,t|law.at_temperature(t),&graph).unwrap();}
                else {network.step_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),solve(),
                    &[0.0;3],&gate,&links(),2,|_,_,t|law.at_temperature(t)).unwrap();}
            }
            outcomes.push((heat.cells().clone(),owner.electrical().unwrap().clone(),
                owner.mechanical().material_memory().to_vec(),placement.physical_state(owner.mechanical().components(),&gate).unwrap()));
        }
        assert_eq!(outcomes[0],outcomes[1]);
    });
}

#[test]
fn invalid_mass_path_and_late_failure_preserve_the_complete_circulating_history() {
    #[derive(Debug)] enum Failure { Late, Material(MaterialStatePointError) }
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t).map_err(Failure::Material)};
        let graph=flow(0.2,true,false,4000.0);let force=vec![0.0;owner.mechanical().mode_count()];
        let cells=heat.cells().clone();let memories=owner.mechanical().material_memory().to_vec();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let material=owner.material().members()[0].source.clone();
        let bad=Advected::new((),vec![],2).unwrap().with_boundaries(vec![
            MassBoundary::Inlet {id:1,cell:CellId(0),material_id:MaterialId(0),kg_s:1.0,specific_h_j_kg:10.0*CP},
            MassBoundary::Outlet {id:2,cell:CellId(0),kg_s:1.0}],2).unwrap();
        let unbalanced=Advected::new((),vec![MassFlow {id:1,a:COOLANT,b:CellId(8),kg_s:0.2}],1).unwrap();
        for bad in [&bad,&unbalanced] {
            assert!(network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),settings(4),
                &[0.0;4],&gate,&links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material),bad).is_err());
        }
        let calls=std::cell::Cell::new(0);
        let result=network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),settings(4),
            &[0.0;4],&gate,&links(),2,|_,_,t| {
                calls.set(calls.get()+1);if calls.get()>2 {return Err(Failure::Late);}
                law.at_temperature(t).map_err(Failure::Material)
            },&graph);
        assert!(matches!(result,Err(RelaxationTransportError::Source(ConvectionHeatError::Convection(
            SourceCellConvectionError::Source {source:Failure::Late,..})))));
        let mut exhausted=settings(4);exhausted.maximum_material_trials=1;
        assert!(network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),exhausted,
            &[0.0;4],&gate,&links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material),&graph).is_err());
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(matches!(network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),settings(4),
            &[0.0;4],&cancelled,&links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material),&graph),
            Err(RelaxationTransportError::Relaxation(CircuitEnthalpyError::Cancelled))));
        assert_eq!(heat.cells(),&cells);assert_eq!(owner.electrical().unwrap(),&circuit);
        assert_eq!(owner.mechanical().material_memory(),memories);assert_eq!(owner.material().members()[0].source,material);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        assert_eq!(heat.time_s(),0.0);assert_eq!(owner.mechanical().samples_rendered(),0);
        network.step_advected_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),settings(4),
            &[0.0;4],&gate,&links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material),&graph).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

mod resolved;

//! Independent pressure/RC/cell balances on a continuing original coolant loop.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::heat::relaxation::pressure::*;
use fs_airflow::graph::thermal::medium::series::hydraulic::{MaterialFlowResistance,
    pressure::SeriesPressureBudget};
use fs_qty::Area;
const GEOMETRY_A: f64=1e10;
const GEOMETRY_B: f64=2e6;
const FLOW_AREA: f64=0.0002;
fn sections()->[PressureCoolantSection;2] {
    [COOLANT,CellId(8)].map(|cell|PressureCoolantSection {cell,resistance:
        MaterialFlowResistance::try_new(GEOMETRY_A,GEOMETRY_B,Area::new(FLOW_AREA),
            Length::new(0.01),(0.0,1e6)).unwrap()})
}
fn patches()->[PressureConvectiveLink;2] {
    links().map(|link|PressureConvectiveLink {link:link.link,section:0,
        characteristic_length:link.characteristic_length,property_temperature:link.property_temperature,
        correlation:CorrelationId::FlatPlateLaminarAverage,settings:ForcedConvectionSettings::default()})
}
fn path<'a>(sections:&'a[PressureCoolantSection],links:&'a[PressureConvectiveLink],drive:f64)->PressureCoolantLoop<'a> {
    PressureCoolantLoop {sections,links,driving_pressure_pa:drive,maximum_links:2,
        pressure_budget:SeriesPressureBudget {maximum_sections:2,absolute_pressure_tolerance_pa:1e-8,relative_pressure_tolerance:1e-12}}
}
fn hydraulic(s:Spec,t:[f64;2],drive:f64)->(f64,[f64;2]) {
    let mu=t.map(|t|s.viscosity*0.001*(1.0-0.01*(t-300.0)));
    let a=GEOMETRY_A*(mu[0]+mu[1])/1000.0;
    let b=2.0*GEOMETRY_B/1000.0;
    // Independent ordinary quadratic inverse, not the scaled native algorithm.
    let mass=(2.0*drive.abs()/(a+(a*a+4.0*b*drive.abs()).sqrt())).copysign(drive);
    let q=mass/1000.0;
    let heat=mu.map(|mu|(GEOMETRY_A*mu*q+GEOMETRY_B*1000.0*q*q.abs())*q);
    (mass,heat)
}
fn controls4()->FrameThermalSolve {FrameThermalSolve {maximum_cells:4,..solve()}}

#[test]
fn temperature_resolved_pressure_speed_enthalpy_and_hydraulic_heat_close_original_balances() {
    let mut outcomes=Vec::new();
    for (s,drive) in [(Spec {circulating:true,..spec()},4000.0),
        (Spec {circulating:true,viscosity:2.0,..spec()},4000.0),
        (Spec {circulating:true,..spec()},8000.0),
        (Spec {circulating:true,..spec()},-4000.0),
        (Spec {circulating:true,powered:false,..spec()},4000.0)] {
        fixture(s,|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            if s.powered {owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();}
            let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let sections=sections();let patches=patches();
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            let e0=owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                +if s.powered {CHARGE*CHARGE/(2.0*CAP)}else{0.0};
            let (mut supplied,mut outside)=(0.0,0.0);let mut first_rate=0.0;let mut last_rate=0.0;
            for step in 1..=12 {
                let before=heat.cells().clone();let q0=owner.electrical().map_or(0.0,|c|c.charges_c()[0]);
                let report=network.step_pressure_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,
                    s.powered.then_some(input()),controls4(),&[0.0;4],&gate,path(&sections,&patches,drive),
                    |_,_,t|law.at_temperature(t)).unwrap();
                let payload=report.transport_source.as_ref().unwrap();
                let t=[heat.temperature_k(COOLANT).unwrap(),heat.temperature_k(CellId(8)).unwrap()];
                let (mass,power)=hydraulic(s,t,drive);let speed=mass.abs()/1000.0/FLOW_AREA;
                near(payload.hydraulics.mass_flow_kg_s(),mass,1e-12);
                if step==1 {first_rate=mass;}last_rate=mass;
                for i in 0..2 {
                    let actual=&payload.hydraulics.evaluations()[i];
                    let original=law.at_temperature(Temperature::new(t[i])).unwrap();
                    assert_eq!(payload.hydraulics.sections()[i].medium,original);
                    near(actual.dissipated_power_w,power[i],1e-12);
                    near(payload.hydraulic_heat_w[&sections[i].cell],power[i],1e-12);
                }
                near(payload.hydraulics.dissipated_power_w(),drive*mass/1000.0,1e-11);
                near(report.external_thermal_heat_j,DT*(power[0]+power[1]),1e-8);
                assert_eq!(report.enthalpy_flow.as_ref().unwrap().external_enthalpy_j,0.0);
                let mut received=0.0;
                for sample in payload.convection.finite.samples() {
                    let i=sample.wall_cell.0 as usize;let wall=heat.temperature_k(sample.wall_cell).unwrap();
                    let expected=PATCHES[i]*h(s,wall,t[0])*speed.sqrt()*(wall-t[0]);
                    near(sample.heat_into_fluid_w,expected,1e-8);received+=expected;
                    let cell=&heat.cells()[&sample.wall_cell];
                    let loss=report.layers.iter().find(|r|r.destination.cell==sample.wall_cell).unwrap().heat_j;
                    let joule=report.electrical.as_ref().map_or(0.0,|e|e.resistors[i].heat_j);
                    near(cell.mass_kg*(cell.specific_h_j_kg-before[&sample.wall_cell].specific_h_j_kg),
                        loss+joule-DT*expected,3e-8);
                }
                let carried=mass.abs()*s.cp*(t[1]-t[0]);
                for i in 0..2 {
                    let id=sections[i].cell;let cell=&heat.cells()[&id];
                    let exchange=if i==0 {received+carried}else{-carried};
                    near(cell.mass_kg*(cell.specific_h_j_kg-before[&id].specific_h_j_kg),DT*(exchange+power[i]),4e-8);
                }
                if let Some(e)=&report.electrical {
                    let r=e.resistors.iter().map(|r|r.material.conductor.resistance_ohm()).sum::<f64>();
                    let a=DT/(2.0*CAP*r);near(owner.electrical().unwrap().charges_c()[0],q0*(1.0-a)/(1.0+a),1e-12);
                }
                supplied+=report.mechanical.external_work_j+DT*(power[0]+power[1])
                    +report.electrical.as_ref().map_or(0.0,|e|e.native.supplied);
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                    +report.electrical.as_ref().map_or(0.0,|e|e.stored_after_j)-e0+outside-supplied,0.0,5e-7);
                assert_eq!(owner.mechanical().samples_rendered(),step);assert_eq!(heat.time_s(),step as f64*DT);
            }
            assert!((first_rate-last_rate).abs()>1e-6); // flow was not frozen at initialization
            outcomes.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(COOLANT).unwrap(),owner.electrical().map_or(0.0,|c|c.charges_c()[0]),last_rate]);
        });
    }
    for i in [1,2] {for k in 0..4 {assert!((outcomes[i][k]-outcomes[0][k]).abs()>1e-12);}}
    for k in 0..3 {near(outcomes[3][k],outcomes[0][k],1e-10);}
    near(outcomes[3][3],-outcomes[0][3],1e-12);
}

#[test]
fn unsupported_bulk_state_duplicate_inventory_and_late_failure_leave_all_histories_for_retry() {
    #[derive(Debug)] enum Failure { Late, Material(MaterialStatePointError) }
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t).map_err(Failure::Material)};
        let sections=sections();let patches=patches();let force=vec![0.0;owner.mechanical().mode_count()];
        let cells=heat.cells().clone();let memory=owner.mechanical().material_memory().to_vec();
        let material=owner.material().members()[0].source.clone();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        for kind in 0..4 {
            let calls=std::cell::Cell::new(0);let duplicate=[sections[0],sections[0]];
            let sections=if kind==0 {&duplicate}else{&sections};
            let mut controls=controls4();if kind==3 {controls.maximum_material_trials=1;}
            let result=network.step_pressure_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),
                controls,&[0.0;4],&gate,path(sections,&patches,4000.0),|_,_,temperature| {
                    let n=calls.get();calls.set(n+1);
                    if kind==2 && n>=4 {return Err(Failure::Late);}
                    let t=if kind==1 {Temperature::new(temperature.value()+1.0)}else{temperature};
                    law.at_temperature(t).map_err(Failure::Material)
                });
            assert!(result.is_err());
            if kind==2 {assert!(matches!(result,Err(RelaxationTransportError::Source(
                PressureCoolantError::Source {source:Failure::Late,..}))));}
            assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);
            assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(owner.mechanical().samples_rendered(),0);
            assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,material);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        }
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(network.step_pressure_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),controls4(),
            &[0.0;4],&cancelled,path(&sections,&patches,4000.0),|_,_,t|law.at_temperature(t).map_err(Failure::Material)).is_err());
        assert_eq!(owner.mechanical().samples_rendered(),0);
        network.step_pressure_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,Some(input()),controls4(),
            &[0.0;4],&gate,path(&sections,&patches,4000.0),|_,_,t|law.at_temperature(t).map_err(Failure::Material)).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

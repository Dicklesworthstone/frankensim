//! Existing pump matcher, original source queries and independent whole-cell work.
use super::*;
use fs_airflow::SourceProvenance;
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::{PressureSection, SeriesPressureFlow};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpHeadPoint, PumpOperatingControls, PumpOperatingError,
};
use fs_qty::VolumetricFlowRate;

fn curve(head: f64, viscosity: (f64, f64)) -> PumpHeadCurve {
    PumpHeadCurve::try_new(vec![
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.0), head: Length::new(head) },
        PumpHeadPoint { flow: VolumetricFlowRate::new(0.0005), head: Length::new(head-0.5) },
    ], 9.81, viscosity, SourceProvenance::new("manufactured fixed-speed head curve", "pump-series-test")).unwrap()
}
fn operating() -> PumpOperatingControls {
    PumpOperatingControls { pressure_bracket_pa: (500.0, 6000.0), absolute_pressure_tolerance_pa: 1e-8,
        relative_pressure_tolerance: 1e-12, max_network_evaluations: 64 }
}
fn pumped<'a>(sections: &'a [PressureCoolantSection], links: &'a [PressureConvectiveLink],
    pump: &'a PumpHeadCurve) -> PumpCoolantLoop<'a> {
    PumpCoolantLoop { sections, links, pump, pressure_budget: path(sections,links,0.0).pressure_budget,
        operating_point: operating(), maximum_links: 2 }
}
fn independent(s: Spec, t: [f64; 2], head: f64) -> (f64, f64) {
    let mu = t.map(|t| s.viscosity*0.001*(1.0-0.01*(t-300.0)));
    // Native matcher bisects pressure; this solves its scalar intersection
    // analytically in MASS flow, without reusing native reductions/matching.
    let a = GEOMETRY_A*(mu[0]+mu[1])/1000.0 + 9.81*1000.0;
    let b = 2.0*GEOMETRY_B/1000.0;
    let p0 = 1000.0*9.81*head;
    let mass = 2.0*p0/(a+(a*a+4.0*b*p0).sqrt());
    (mass, 1000.0*9.81*(head-1000.0*mass/1000.0))
}

#[test]
fn original_series_and_pump_sources_match_the_analytic_intersection_and_fixed_drive() {
    let s=Spec {circulating:true,..spec()};
    fixture(s,|_,_,_,_,_,_,law,_,cx| {
        let gate=CancelGate::new_clock_free(); let pump=curve(0.6,(1e-7,5e-6));
        let t=[310.0,295.0]; let geometry=sections(); let links=patches();
        let samples:Vec<_>=t.into_iter().map(|t|law.at_temperature(Temperature::new(t)).unwrap()).collect();
        let prepared=||geometry.iter().zip(&samples).map(|(g,s)|PressureSection {resistance:g.resistance,medium:s.clone()}).collect();
        let p=path(&geometry,&links,0.0).pressure_budget;
        let result=SeriesPressureFlow::solve_with_pump(prepared(),&pump,p,operating(),cx,&gate).unwrap();
        let (mass,pressure)=independent(s,t,0.6);
        near(result.state.mass_flow_kg_s(),mass,1e-11);
        near(result.network_pressure_pa,pressure,2e-8);
        assert!(core::ptr::eq(result.curve,&pump));
        assert_eq!(result.network_evaluations,result.total_flow_evaluations);
        assert_eq!(result.total_solid_solves,0);
        assert!(result.network_evaluations<=operating().max_network_evaluations);
        let fixed=SeriesPressureFlow::solve(prepared(),result.network_pressure_pa,p,&gate).unwrap();
        assert_eq!(result.state.evaluations(),fixed.evaluations());
        for (actual,original) in result.state.sections().iter().zip(&samples) {assert_eq!(&actual.medium,original);}
        near(result.supply.hydraulic_power_w,pressure*mass/1000.0,2e-11);
        // Another resistance section must not become the inlet-density authority.
        near(result.supply.inlet_kinematic_viscosity_m2_s,0.0009/1000.0,1e-18);
    });
}

#[test]
fn pump_selected_flow_and_hydraulic_work_share_the_thermal_circuit_material_endpoint() {
    let mut endpoints=Vec::new();
    for (s,head) in [(Spec {circulating:true,..spec()},0.6),
        (Spec {circulating:true,viscosity:2.0,..spec()},0.6),
        (Spec {circulating:true,..spec()},0.8),
        (Spec {circulating:true,powered:false,..spec()},0.6)] {
        fixture(s,|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            if s.powered {owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();}
            let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let geometry=sections();let links=patches();let pump=curve(head,(1e-7,5e-6));
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            let energy0=owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                +if s.powered {CHARGE*CHARGE/(2.0*CAP)}else{0.0};
            let (mut supplied,mut outside,mut first,mut last)=(0.0,0.0,0.0,0.0);
            for n in 1..=12 {
                let cells=heat.cells().clone();let q0=owner.electrical().map_or(0.0,|e|e.charges_c()[0]);
                let report=network.step_pump_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,
                    s.powered.then_some(input()),controls4(),&[0.0;4],&gate,pumped(&geometry,&links,&pump),
                    |_,_,t|law.at_temperature(t)).unwrap();
                let accepted=report.transport_source.as_ref().unwrap();
                let coolant=&accepted.coolant;let receipt=&accepted.pump;
                let t=[heat.temperature_k(COOLANT).unwrap(),heat.temperature_k(CellId(8)).unwrap()];
                let (mass,pressure)=independent(s,t,head);let (_,power)=hydraulic(s,t,pressure);
                near(coolant.hydraulics.mass_flow_kg_s(),mass,1e-11);
                near(receipt.supply.available_pressure_pa,pressure,2e-8);
                near(receipt.supply.hydraulic_power_w,pressure*mass/1000.0,2e-11);
                near(receipt.hydraulic_power_residual_w,receipt.supply.hydraulic_power_w
                    -coolant.hydraulics.dissipated_power_w(),1e-15);
                assert!(core::ptr::eq(receipt.curve,&pump));
                assert!(receipt.pressure_residual_pa.abs()<2e-8);
                near(report.external_thermal_heat_j,DT*receipt.supply.hydraulic_power_w,2e-8);
                assert_eq!(report.enthalpy_flow.as_ref().unwrap().external_enthalpy_j,0.0);
                if n==1 {first=mass;}last=mass;
                for i in 0..2 {
                    assert_eq!(coolant.hydraulics.sections()[i].medium,law.at_temperature(Temperature::new(t[i])).unwrap());
                    near(coolant.hydraulic_heat_w[&geometry[i].cell],power[i],2e-11);
                }
                let speed=mass/1000.0/FLOW_AREA;let mut received=0.0;
                for sample in coolant.convection.finite.samples() {
                    let i=sample.wall_cell.0 as usize;let wall=heat.temperature_k(sample.wall_cell).unwrap();
                    let flux=PATCHES[i]*h(s,wall,t[0])*speed.sqrt()*(wall-t[0]);
                    near(sample.heat_into_fluid_w,flux,2e-8);received+=flux;
                    let cell=&heat.cells()[&sample.wall_cell];
                    let relaxation=report.layers.iter().find(|r|r.destination.cell==sample.wall_cell).unwrap().heat_j;
                    let joule=report.electrical.as_ref().map_or(0.0,|e|e.resistors[i].heat_j);
                    near(cell.mass_kg*(cell.specific_h_j_kg-cells[&sample.wall_cell].specific_h_j_kg),
                        relaxation+joule-DT*flux,4e-8);
                }
                let carried=mass*s.cp*(t[1]-t[0]);
                for i in 0..2 {
                    let id=geometry[i].cell;let c=&heat.cells()[&id];
                    near(c.mass_kg*(c.specific_h_j_kg-cells[&id].specific_h_j_kg),
                        DT*(power[i]+if i==0 {received+carried}else{-carried}),5e-8);
                }
                if let Some(e)=&report.electrical {
                    let r=e.resistors.iter().map(|r|r.material.conductor.resistance_ohm()).sum::<f64>();
                    let a=DT/(2.0*CAP*r);near(owner.electrical().unwrap().charges_c()[0],q0*(1.0-a)/(1.0+a),1e-12);
                }
                supplied+=report.mechanical.external_work_j+DT*receipt.supply.hydraulic_power_w
                    +report.electrical.as_ref().map_or(0.0,|e|e.native.supplied);
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                    +report.electrical.as_ref().map_or(0.0,|e|e.stored_after_j)-energy0+outside-supplied,0.0,6e-7);
                assert_eq!(owner.mechanical().samples_rendered(),n);assert_eq!(heat.time_s(),n as f64*DT);
            }
            assert!((last-first).abs()>1e-6);
            endpoints.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(COOLANT).unwrap(),owner.electrical().map_or(0.0,|e|e.charges_c()[0]),last]);
        });
    }
    for i in [1,2] {for k in 0..4 {assert!((endpoints[i][k]-endpoints[0][k]).abs()>1e-12);}}
}

#[test]
fn unsupported_pump_or_late_trial_preserves_charge_cells_material_and_fluid_histories() {
    fixture(Spec {circulating:true,..spec()},|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
        let geometry=sections();let links=patches();let pump=curve(0.6,(1e-7,5e-6));
        let unsupported=curve(0.6,(1e-10,2e-10));
        let cells=heat.cells().clone();let memory=owner.mechanical().material_memory().to_vec();
        let material=owner.material().members()[0].source.clone();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let force=vec![0.0;owner.mechanical().mode_count()];
        for fault in 0..5 {
            let mut p=pumped(&geometry,&links,if fault==0 {&unsupported}else{&pump});
            let mut controls=controls4();
            match fault {
                1=>p.operating_point.pressure_bracket_pa=(200.0,500.0),
                2=>p.operating_point.max_network_evaluations=2,
                3=>controls.maximum_material_trials=1,
                // A loose pressure tolerance cannot hide missing hydraulic work.
                4=>p.operating_point.absolute_pressure_tolerance_pa=1e9,
                _=>{},
            }
            let result=network.step_pump_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,
                Some(input()),controls,&[0.0;4],&gate,p,|_,_,t|law.at_temperature(t));
            assert!(result.is_err());
            if fault==0 {assert!(matches!(result,Err(RelaxationTransportError::Source(
                PressureCoolantError::Pump(PumpOperatingError::OutsideCurve {..})))));}
            if fault==1 {assert!(matches!(result,Err(RelaxationTransportError::Source(
                PressureCoolantError::Pump(PumpOperatingError::NotBracketed {..})))));}
            if fault==2 {assert!(matches!(result,Err(RelaxationTransportError::Source(
                PressureCoolantError::Pump(PumpOperatingError::NotConverged {..})))));}
            assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);
            assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(owner.mechanical().samples_rendered(),0);
            assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,material);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        }
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(network.step_pump_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,
            Some(input()),controls4(),&[0.0;4],&cancelled,pumped(&geometry,&links,&pump),|_,_,t|law.at_temperature(t)).is_err());
        assert_eq!(owner.mechanical().samples_rendered(),0);
        network.step_pump_driven_relaxation_with_links(cx,&mut owner,&mut heat,&force,
            Some(input()),controls4(),&[0.0;4],&gate,pumped(&geometry,&links,&pump),|_,_,t|law.at_temperature(t)).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

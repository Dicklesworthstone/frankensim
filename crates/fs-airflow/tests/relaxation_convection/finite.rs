//! One finite coolant receives reciprocal heat from both original solid layers.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::links::{
    SourceConvectiveLink,evaluate_source_link_convection,
};
use std::collections::BTreeMap;
const COOLANT:CellId=CellId(7);
const PATCHES:[f64;2]=[0.01,0.04];
#[derive(Clone,Copy)] struct Spec {powered:bool,k:f64,cp:f64,mass:f64,reversed:bool,reservoir:bool,circulating:bool,
    viscosity:f64,density_slope:f64,coolant_initial_k:[f64;2],adiabatic:bool}
fn spec()->Spec {Spec {powered:true,k:1.0,cp:100.0,mass:0.5,reversed:false,reservoir:false,circulating:false,
    viscosity:1.0,density_slope:0.0,coolant_initial_k:[310.0,295.0],adiabatic:false}}
fn finite_source(s:Spec)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let mut claims=ClaimSet::new();
    for (name,dims,value) in [("density",Density::DIMS,1000.0),("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,s.cp)] {
        let varying=name=="density" && s.density_slope!=0.0;
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:if varying {
            PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*s.density_slope),(340.0,value+40.0*s.density_slope)],dims}
            }else{PropertyValue::Scalar {value,dims}},
            validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),280.0,340.0).with("p",100000.0,100000.0),
            interpolation:if varying {InterpolationPolicy::LinearInside}else{InterpolationPolicy::ConstantWithinValidity},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured finite coolant".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    for (name,dims,value,slope) in [("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,s.k*0.5,s.k*0.002),
        ("dynamic_viscosity",DYNAMIC_VISCOSITY_DIMS,s.viscosity*0.001,-s.viscosity*0.00001)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Curve {
            abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots:vec![(280.0,value-20.0*slope),(340.0,value+40.0*slope)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),280.0,340.0).with("p",100000.0,100000.0),
            interpolation:InterpolationPolicy::LinearInside,uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured finite coolant transport".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"finite coolant".into(),phase:"liquid".into(),
        process:"authored test law".into(),revision:0},claims,vec![]).unwrap();
    let query=|t|point(t).with("p",100000.0).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),density:&PropertyKey::new("density",Density::DIMS),
        temperature_axis:"T",lower:&query(280.0),upper:&query(340.0),selection:SelectionPolicy::SingleClaimOnly},
        SolidLiquidPhase::Liquid,0.0,1e-6).unwrap();
    (card,curve)
}
fn links()->[SourceConvectiveLink;2] {
    [0,1].map(|i|SourceConvectiveLink {link:20+i,fluid_cell:COOLANT,characteristic_length:Length::new(0.1),
        property_temperature:ConvectionPropertyTemperature::FilmMean,law:SourceCellConvectionLaw::Forced {
            correlation:CorrelationId::FlatPlateLaminarAverage,speed_m_s:1.0,settings:ForcedConvectionSettings::default()}})
}
fn h(s:Spec,wall:f64,fluid:f64)->f64 {
    let d=0.5*(wall+fluid)-300.0;let k=s.k*(0.5+0.002*d);let mu=s.viscosity*0.001*(1.0-0.01*d);
    0.664*(1000.0*0.1/mu).sqrt()*(mu*s.cp/k).cbrt()*k/0.1
}
fn solve()->FrameThermalSolve {FrameThermalSolve {maximum_cells:3,maximum_material_trials:256,..controls()}}
fn fixture<R>(s:Spec,run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,&MaterialNetwork<'_>,&MaterialNetwork<'_>,
    MaterialCircuit,DepositionSimulation<'_>,&fs_couple::circuit_enthalpy::contacts::heat::PreparedHeatNetwork<'_,'_>,
    &LiquidTemperatureLaw<'_>,&BTreeMap<u64,f64>,&Cx<'_>)->R)->R {
    let (a,ac)=source(0,true);let (b,bc)=source(1,true);let (fluid,fc)=finite_source(s);
    let cards=[a,b];let curves=[ac,bc];let f=frame(&cards);let m=model(&f);
    let solid_query=point(300.0);let fluid_query=point(310.0).with("p",100000.0).unwrap();
    let law=LiquidTemperatureLaw::try_new(&fluid,&fluid_query,"T",LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let mut cells:Vec<_>=(0..2).map(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context:&solid_query,temperature_axis:"T",liquid_fraction_axis:None},conductivity:None}).collect();
    cells.push(MaterialCell {cell:COOLANT,material_id:MaterialId(7),curve:fc.curve(),point:ThermalMaterialPoint {
        material:&fluid,context:&fluid_query,temperature_axis:"T",liquid_fraction_axis:None},conductivity:None});
    if s.circulating {cells.push(MaterialCell {cell:CellId(8),material_id:MaterialId(7),curve:fc.curve(),point:ThermalMaterialPoint {
        material:&fluid,context:&fluid_query,temperature_axis:"T",liquid_fraction_axis:None},conductivity:None});}
    let mut events:Vec<_>=(0..2).map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(i as u64),segment_id:0,body_id:0,
        material_id:MaterialId(i as u64),mass_kg:[3.0,5.0][i]*area(i),specific_h_j_kg:10.0*CP,
        geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area(i),surface_area_m2:exposed_area(i)},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})}).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:COOLANT,segment_id:0,body_id:1,material_id:MaterialId(7),
        mass_kg:s.mass,specific_h_j_kg:(s.coolant_initial_k[0]-280.0)*s.cp,geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],
            volume_m3:s.mass/1000.0,surface_area_m2:PATCHES.iter().sum()},boundary:Boundary {ambient_k:310.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});
    if s.circulating {events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(8),segment_id:0,body_id:1,material_id:MaterialId(7),
        mass_kg:s.mass,specific_h_j_kg:(s.coolant_initial_k[1]-280.0)*s.cp,geometry:CellGeometry {start_m:[1.0,0.0,0.0],end_m:[2.0,0.0,0.0],
            volume_m3:s.mass/1000.0,surface_area_m2:0.05},boundary:Boundary {ambient_k:295.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});}
    for i in 0..if s.adiabatic {0}else{2} {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:20+i as u64,
        a:if s.reversed {COOLANT}else{CellId(i as u64)},b:if s.reversed {CellId(i as u64)}else{COOLANT},
        area_m2:PATCHES[i],conductance_w_k:0.0})});}
    let heat=DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw),(MaterialId(7),fc.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:4,max_links:2,max_events:6,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:DT,energy_tolerance_j:1e-10},||Ok(())).unwrap().0;
    let bindings=if s.adiabatic {Vec::new()}else{links().to_vec()};
    let initial=evaluate_source_link_convection(&heat,&bindings,2,|_,_,t|law.at_temperature(t),||false).unwrap();
    let resistors:[_;2]=core::array::from_fn(|i|CellResistor {destination:ResistorCell {branch:i+1,cell:CellId(i as u64)},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area(i)});
    let plain=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let electric=MaterialNetwork {cells:&cells,resistors:&resistors,links:plain.links};
    let bases=[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let resolved=electric.prepare_with_link_conductances(&heat,&bases,initial.conductances_w_k(),||false).unwrap()
        .resolve_with_link_conductances(&heat,initial.conductances_w_k(),||false).unwrap();
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![(1,0,Branch::Capacitor {farads:CAP}),
        (1,2,Branch::Resistor {ohms:999.0}),(2,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},resolved.resistors,vec![],vec![CHARGE],&[]).unwrap();
    let electrical=ElectricalContactNetwork {bulk:if s.powered {&electric}else{&plain},contacts:&[],electrical_bases:if s.powered {&bases}else{&[]}};
    let assembly=MaterialHeatNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:14,kernel_id:42,tile:0,iteration:0},fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        let prepared=assembly.prepare_with_link_conductances(&cx,&heat,initial.conductances_w_k()).unwrap();
        run(&m,&plain,&electric,circuit,heat,&prepared,&law,initial.conductances_w_k(),&cx)
    })
}

#[test]
fn both_walls_share_original_finite_coolant_with_independent_cell_and_total_energy() {
    let mut outcomes=Vec::new();
    for s in [spec(),Spec {k:2.0,..spec()},Spec {cp:200.0,..spec()},Spec {mass:1.0,..spec()},
        Spec {powered:false,..spec()},Spec {reversed:true,..spec()},Spec {reservoir:true,..spec()}] {
        fixture(s,|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let initial_energy=native.total_energy_j().unwrap()+caloric(&heat)+if s.powered {CHARGE*CHARGE/(2.0*CAP)}else{0.0};
            let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
                native,&heat,plain,initial,budget(),&gate).unwrap();
            if s.powered {owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit,
                &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();}
            let reservoir=boundaries();let network=SourceConvectionHeatNetwork {heat:prepared,cells:if s.reservoir {&reservoir}else{&[]},maximum_cells:2,
                resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t)};
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();force.resize(owner.mechanical().mode_count(),0.0);
            let (mut supply,mut outside)=(0.0,0.0);
            for step in 1..=12 {
                let before=heat.cells().clone();let q0=owner.electrical().map_or(0.0,|c|c.charges_c()[0]);
                let report=network.step_relaxation_with_links(cx,&mut owner,&mut heat,&force,s.powered.then_some(input()),solve(),&[0.0;3],&gate,
                    &links(),2,|_,_,t|law.at_temperature(t)).unwrap();
                let payload=report.transport_source.as_ref().unwrap();let fluid=heat.temperature_k(COOLANT).unwrap();
                let mut joule=[0.0;2];let mut stored=0.0;let mut received=0.0;let mut entropy=0.0;
                if let Some(e)=&report.electrical {
                    let r:[f64;2]=core::array::from_fn(|i|[0.2,0.6][i]*(1.0+0.02*(heat.temperature_k(CellId(i as u64)).unwrap()-300.0))/area(i));
                    let a=DT/(2.0*CAP*(r[0]+r[1]));let q=q0*(1.0-a)/(1.0+a);
                    near(owner.electrical().unwrap().charges_c()[0],q,1e-12);
                    for i in 0..2 {joule[i]=DT*((q0+q)/(2.0*CAP*(r[0]+r[1]))).powi(2)*r[i];}
                    stored=e.stored_after_j;supply+=e.native.supplied;
                }
                for sample in payload.finite.samples() {
                    let i=sample.wall_cell.0 as usize;let wall=heat.temperature_k(sample.wall_cell).unwrap();
                    let power=PATCHES[i]*h(s,wall,fluid)*(wall-fluid);received+=DT*power;
                    entropy+=DT*power*(1.0/fluid-1.0/wall);
                    near(sample.heat_into_fluid_w,power,1e-8);
                    near(sample.property_temperature_k,0.5*(wall+fluid),1e-12);
                    assert_eq!(sample.material.material().card_identity(),plain.cells[2].point.material.content_hash());
                    assert_eq!(sample.material.material().query().axes()["p"],100000.0);
                    let loss=report.layers.iter().find(|l|l.destination.cell==sample.wall_cell).unwrap().heat_j;
                    let ambient=payload.reservoirs.convection().sources_w().get(&sample.wall_cell).copied().unwrap_or(0.0);
                    let cell=&heat.cells()[&sample.wall_cell];
                    near(cell.mass_kg*(cell.specific_h_j_kg-before[&sample.wall_cell].specific_h_j_kg),
                        loss+joule[i]-DT*power+DT*ambient,2e-8);
                }
                let cell=&heat.cells()[&COOLANT];
                near(cell.mass_kg*(cell.specific_h_j_kg-before[&COOLANT].specific_h_j_kg),received,2e-8);
                near(report.external_thermal_heat_j,DT*payload.reservoirs.convection().total_heat_into_cells_w(),1e-8);
                near(report.entropy.conduction_generation_j_k,entropy,1e-9);assert!(entropy>=0.0);
                assert_eq!(report.entropy.has_unclassified_sources,s.reservoir);
                supply+=report.mechanical.external_work_j+report.external_thermal_heat_j;
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+stored-initial_energy+outside-supply,0.0,3e-7);
                assert_eq!(heat.time_s(),step as f64*DT);assert_eq!(owner.mechanical().samples_rendered(),step);
            }
            near(owner.mechanical().material_memory()[0],0.007*(-12.0*DT/0.7).exp(),1e-10);
            assert!(heat.temperature_k(COOLANT).unwrap()<310.0 && heat.temperature_k(COOLANT).unwrap()>300.0);
            outcomes.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(CellId(1)).unwrap(),heat.temperature_k(COOLANT).unwrap(),owner.electrical().map_or(0.0,|c|c.charges_c()[0])]);
        });
    }
    for i in [1,2,3] {for j in 0..4 {assert!((outcomes[i][j]-outcomes[0][j]).abs()>1e-12);}}
    for j in 0..4 {near(outcomes[5][j],outcomes[0][j],1e-10);}
}

#[test]
fn incomplete_links_and_late_coolant_refusal_preserve_the_whole_finite_inventory_for_retry() {
    #[derive(Debug)] enum Failure { Late, Material(MaterialStatePointError) }
    fixture(spec(),|m,plain,electric,circuit,mut heat,prepared,law,initial,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new_shared_with_link_conductances(vec![(material,DEST.to_vec())],
            native,&heat,plain,initial,budget(),&gate).unwrap();
        owner.bind_material_circuit_with_link_conductances(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],vec![],initial,&gate).unwrap();
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
            resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t|law.at_temperature(t).map_err(Failure::Material)};
        let external=vec![0.0;owner.mechanical().mode_count()];let before=heat.cells().clone();
        let memory=owner.mechanical().material_memory().to_vec();let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let source=owner.material().members()[0].source.clone();
        assert!(network.step_relaxation_with_links(cx,&mut owner,&mut heat,&external,Some(input()),solve(),&[0.0;3],&gate,
            &links()[..1],2,|_,_,t|law.at_temperature(t).map_err(Failure::Material)).is_err());
        let calls=std::cell::Cell::new(0);
        let error=network.step_relaxation_with_links(cx,&mut owner,&mut heat,&external,Some(input()),solve(),&[0.0;3],&gate,
            &links(),2,|_,_,t| {
                calls.set(calls.get()+1);if calls.get()>2 {return Err(Failure::Late);}
                law.at_temperature(t).map_err(Failure::Material)
            });
        assert!(matches!(error,Err(RelaxationTransportError::Source(ConvectionHeatError::Convection(
            SourceCellConvectionError::Source {source:Failure::Late,..})))));
        let mut exhausted=solve();exhausted.maximum_material_trials=1;
        assert!(network.step_relaxation_with_links(cx,&mut owner,&mut heat,&external,Some(input()),exhausted,&[0.0;3],&gate,
            &links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material)).is_err());
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(matches!(network.step_relaxation_with_links(cx,&mut owner,&mut heat,&external,Some(input()),solve(),&[0.0;3],&cancelled,
            &links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material)),
            Err(RelaxationTransportError::Relaxation(CircuitEnthalpyError::Cancelled))));
        assert_eq!(heat.cells(),&before);assert_eq!(owner.electrical().unwrap(),&circuit);
        assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,source);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        assert_eq!(heat.time_s(),0.0);assert_eq!(owner.mechanical().samples_rendered(),0);
        network.step_relaxation_with_links(cx,&mut owner,&mut heat,&external,Some(input()),solve(),&[0.0;3],&gate,
            &links(),2,|_,_,t|law.at_temperature(t).map_err(Failure::Material)).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

mod advection;

mod pressure;

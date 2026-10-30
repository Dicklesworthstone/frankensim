//! Source coolant data changes the continuing native Maxwell/fluid/circuit trajectory.
//! All cards are manufactured numerical fixtures, not measured materials.
use fs_couple::air_path::cylinder::frame::modal::relaxation::{FrameRelaxationBinding};
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::{
    FrameRelaxationThermalSystem, transport::{RelaxationHeatInput, RelaxationTransportError},
    electrical::RelaxationElectricalInput,
};
use fs_couple::frame_vibration::thermal::{FrameLayerCell,FrameThermalSolve};
use fs_couple::render::schedule::force::coupled::CoupledModalSystem;
use fs_couple::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError,ResistorCell,
    material::{CellResistor,MaterialNetwork},contacts::{heat::MaterialHeatNetwork,
    material::ElectricalContactNetwork}};
use fs_couple::material_circuit::{MaterialCircuit,MaterialCircuitSolve};
use fs_airflow::graph::thermal::medium::convection::deposition::source::{
    SourceConvectiveCell,SourceCellConvectionLaw,SourceCellConvectionError,ConvectionPropertyTemperature,
    heat::{SourceConvectionHeatNetwork,ConvectionHeatError},
};
use fs_airflow::graph::thermal::medium::convection::ForcedConvectionSettings;
use fs_convection::CorrelationId;
use fs_material::liquid_temperature::{LiquidTemperatureLaw,LiquidViscosityBasis};
use fs_material::liquid::{ResolvedLiquidState,DYNAMIC_VISCOSITY_DIMS};
use fs_material::conductor::{ElectricalPropertyBasis,ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::state_point::{THERMAL_CONDUCTIVITY_DIMS,MaterialStatePointError};
use fs_exec::{CancelGate,Cx};
use fs_qty::Length;
use fs_phs::circuit::{Branch,CircuitGraph};
use fs_couple::air_path::cylinder::{CylinderFlowResolution,frame::{MaterialFluidFrame,FrameFluidBinding,FrameFluidBudget}};
use fs_couple::air_path::cylinder::frame::modal::FrameFluidPlacement;
use fs_couple::air_path::cylinder::frame::modal::relaxation::{FrameLongitudinalRelaxation,FrameRelaxationBudget};
use fs_couple::frame_vibration::FrameVibrationBudget;
use fs_couple::modal_acoustic_time::{ModalAcousticState,ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::circuit_enthalpy::material::MaterialCell;
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_evidence::ValidityDomain;
use fs_matdb::{MaterialCard,ClaimSet,MaterialStateId,PropertyClaim,PropertyKey,PropertyValue,
    Provenance,QueryPoint,InterpolationPolicy,UncertaintyModel,SelectionPolicy};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_material::state_point::{ScalarPropertyRequirement,ScalarAdmissibility,MaterialPropertySelection,
    resolve_material_state_point,SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_material::phase::{ResolvedHeatCapacityCurve,HeatCapacitySourceQuery,SolidLiquidPhase};
use fs_phs::RelaxationBranch;
use fs_qty::{Density,Pressure,Time,Temperature,QuantitySpec};
use fs_couple::frame_vibration::{Frame2d,FrameMassFormulation,Geometry};
use fs_thermal::enthalpy::*;
use fs_solver::NewtonKrylovConfig;

const CP:f64=2000.0; // Manufactured heat capacity, NOT measured material data.
const RADII:[f64;2]=[0.06,0.1];
fn point(t:f64)->QueryPoint {QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),t).unwrap().with("omega",1.0).unwrap()}
fn source(layer:usize,varying:bool)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let mut claims=ClaimSet::new();let tau=[0.025,0.12][layer];
    for (name,dims,value) in [("density",Density::DIMS,[3.0,5.0][layer]),
        ("young_modulus",Pressure::DIMS,[1000.0,1500.0][layer]),
        ("equilibrium",Pressure::DIMS,[1000.0,1500.0][layer]),
        ("increment",Pressure::DIMS,[600.0,900.0][layer]),("time",Time::DIMS,tau),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,CP),
        ("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS,[0.2,0.6][layer])] {
        let curve=(name=="time"&&varying)||name=="electrical_resistivity";
        let mut validity=ValidityDomain::unconstrained().with("T",290.0,340.0);
        if name=="time" {validity=validity.with("omega",0.0,10000.0);}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if curve {PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:if name=="time" {vec![(290.0,1.2*tau),(300.0,tau),(340.0,0.2*tau)]}
                    else {vec![(290.0,0.8*value),(300.0,value),(340.0,1.8*value)]},dims}}
                else {PropertyValue::Scalar {value,dims}},validity,
            interpolation:if curve {InterpolationPolicy::LinearInside}else{InterpolationPolicy::ConstantWithinValidity},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured self-heating section".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:format!("layer {layer}"),phase:"solid".into(),process:"synthetic".into(),revision:0},claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower:&point(290.0),upper:&point(340.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,curve)
}
fn frame(cards:&[MaterialCard;2])->Frame2d {
    let requirements=[("density",Density::DIMS),("young_modulus",Pressure::DIMS),("equilibrium",Pressure::DIMS),
        ("increment",Pressure::DIMS),("time",Time::DIMS)].map(|(n,d)|ScalarPropertyRequirement::try_new(n,d,ScalarAdmissibility::Finite).unwrap());
    let section=ConcentricElasticSection::try_new(0.0,cards.iter().enumerate().map(|(layer,card)|ConcentricMaterialLayer {
        outer_radius_m:RADII[layer],material:resolve_material_state_point(card,&point(300.0),&requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()}).collect()).unwrap();
    let mut f=Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    f.add_node(0.0,0.0);f.add_node(1.0,0.0);f.fix(0,[true;3]);f.fix(1,[true,false,true]);
    f.add_material_element(0,1,3,section).unwrap();f
}
fn model(frame:&Frame2d)->MaterialFluidFrame<'_,FluidTransportState> {
    MaterialFluidFrame::new(frame,1,vec![FrameFluidBinding {element:0,fluid:FluidTransportState::try_new(2.0,4000.0,0.5,0.02).unwrap()}],
        CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:2.0,intervals:1},
        FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},FrameFluidBudget {max_states:8,max_dense_entries:64,
            max_projection_terms:4096,mass_inverse_tolerance:1e-10},&CancelGate::new_clock_free()).unwrap()
}
fn native<'m,'f>(model:&'m MaterialFluidFrame<'f,FluidTransportState>)
    ->(FrameFluidPlacement<'m,'f,FluidTransportState>,FrameRelaxationBinding<'m,'f,FluidTransportState>,CoupledModalSystem) {
    let gate=CancelGate::new_clock_free();
    let x=model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:0.0003,velocity_m_sqrt_kg_per_s:-0.001}],&[0.0001;4]).unwrap();
    let cfg=ModalCouplingConfig {max_modes:16,max_connections:32,max_setup_terms:1_000_000,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9};
    let (placement,native)=model.coupled_component(&x,256,ModalAcousticTimeBudget::audible_reference(),cfg,&gate).unwrap().into_system(None,&gate).unwrap();
    let size=2*native.mode_count();
    let native=native.with_material_relaxation(vec![RelaxationBranch {projection:vec![0.0;size],stiffness:2.0,relaxation_time_s:0.7}],Some(vec![0.007]),&gate).unwrap();
    let selected=[FrameLongitudinalRelaxation {element:0,layers:vec![LayerRelaxationProperties {equilibrium_modulus:"equilibrium".into(),
        branches:vec![("increment".into(),"time".into())]};2]}];
    let (material,native)=placement.with_longitudinal_relaxation(native,&selected,Some(vec![0.0002,0.003,-0.001,0.0005,-0.002,0.0007]),
        FrameRelaxationBudget {max_members:1,max_source_arms:2,max_memory_coordinates:6,max_setup_terms:100000},&gate).unwrap();
    (placement,material,native)
}
fn thermal<'a>(curves:&'a [ResolvedHeatCapacityCurve;2],mass_scale:f64)->DepositionSimulation<'a> {
    let mut events=Vec::new();
    for layer in 0..2 {
        let inner=if layer==0 {0.0_f64}else{RADII[layer-1]};
        let area=core::f64::consts::PI*(RADII[layer].powi(2)-inner.powi(2));
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(layer as u64),segment_id:0,body_id:0,
            material_id:MaterialId(layer as u64),mass_kg:mass_scale*[3.0,5.0][layer]*area,specific_h_j_kg:10.0*CP,
            geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area,surface_area_m2:exposed_area(layer)},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});
    }
    DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),(MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:2,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:1.0/256.0,energy_tolerance_j:1e-11},||Ok(())).unwrap().0
}
fn budget()->MaterialHeatBudget {MaterialHeatBudget {max_groups:4,max_dimension:256,max_dense_bytes:64*1024*1024,max_multiply_terms:100_000_000_000}}
fn controls()->FrameThermalSolve {FrameThermalSolve {maximum_cells:2,maximum_material_trials:128,maximum_newton_iterations:12,
    difference_step_j:1e-4,cell_energy_tolerance_j:1e-8,thermal_energy_tolerance_j:1e-8,total_energy_tolerance_j:1e-8,
    newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,..NewtonKrylovConfig::default()}}}
const DEST:[FrameLayerCell;2]=[FrameLayerCell {element:0,layer:0,cell:CellId(0)},FrameLayerCell {element:0,layer:1,cell:CellId(1)}];

fn area(i:usize)->f64 {core::f64::consts::PI*(RADII[i].powi(2)-if i==0 {0.0}else{RADII[i-1].powi(2)})}
fn exposed_area(i:usize)->f64 {2.0*area(i)+if i==1 {2.0*core::f64::consts::PI*RADII[i]}else{0.0}}
fn caloric(t:&DepositionSimulation<'_>)->f64 {t.cells().values().map(|c|c.mass_kg*c.specific_h_j_kg).sum()}
fn near(a:f64,b:f64,tol:f64){assert!((a-b).abs()<=tol,"{a:e} != {b:e}");}
const CAP:f64=0.001;
const CHARGE:f64=0.0002;
const DT:f64=1.0/256.0;
fn input()->RelaxationElectricalInput<'static> {RelaxationElectricalInput {ports:&[],
    solve:MaterialCircuitSolve {energy_tolerance_j:1e-10,power_tolerance_w:1e-9}}}

fn fixture<R>(powered:bool,run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,
    &MaterialNetwork<'_>,&MaterialNetwork<'_>,MaterialCircuit,DepositionSimulation<'_>,
    &fs_couple::circuit_enthalpy::contacts::heat::PreparedHeatNetwork<'_,'_>,&Cx<'_>)->R)->R {
    let (a,ac)=source(0,true);let (b,bc)=source(1,true);let cards=[a,b];let curves=[ac,bc];
    let f=frame(&cards);let m=model(&f);let query=point(300.0);let heat=thermal(&curves,1.0);
    let cells:[_;2]=core::array::from_fn(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),
        curve:curves[i].curve(),point:ThermalMaterialPoint {material:&cards[i],context:&query,
        temperature_axis:"T",liquid_fraction_axis:None},conductivity:None});
    let resistors:[_;2]=core::array::from_fn(|i|CellResistor {destination:ResistorCell {branch:i+1,cell:CellId(i as u64)},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area(i)});
    let plain=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let electric=MaterialNetwork {cells:&cells,resistors:&resistors,links:plain.links};
    let bases=[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let resolved=electric.prepare_with_bases(&heat,&bases,||false).unwrap().resolve(&heat,||false).unwrap();
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:CAP}),(1,2,Branch::Resistor {ohms:999.0}),
        (2,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},resolved.resistors,vec![],vec![CHARGE],&[]).unwrap();
    let electrical=ElectricalContactNetwork {bulk:if powered {&electric}else{&plain},contacts:&[],
        electrical_bases:if powered {&bases}else{&[]}};
    let assembly=MaterialHeatNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:71,kernel_id:29,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        let prepared=assembly.prepare(&cx,&heat).unwrap();run(&m,&plain,&electric,circuit,heat,&prepared,&cx)
    })
}
fn coolant(scale:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,dims,value,slope) in [("density",Density::DIMS,1000.0,-0.5),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,4000.0,2.0),
        ("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,scale*0.5,scale*0.002),
        ("dynamic_viscosity",DYNAMIC_VISCOSITY_DIMS,0.001,-0.00001)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*slope),(340.0,value+40.0*slope)],dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,340.0).with("p",100000.0,100000.0),
            interpolation:InterpolationPolicy::LinearInside,uncertainty:UncertaintyModel::Unstated,
            observations:vec![],provenance:Provenance {source:"synthetic coolant substitution".into(),license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic coolant".into(),phase:"liquid".into(),
        process:"authored transport".into(),revision:0},claims,vec![]).unwrap()
}
fn coolant_law(card:&MaterialCard)->LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(card,&QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),300.0)
        .unwrap().with("p",100000.0).unwrap(),"T",LiquidViscosityBasis::Dynamic,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn boundaries()->[SourceConvectiveCell;2] {
    [0,1].map(|i|SourceConvectiveCell {cell:CellId(i),reference_temperature:Temperature::new(310.0),
        characteristic_length:Length::new(0.1),property_temperature:ConvectionPropertyTemperature::FilmMean,
        law:SourceCellConvectionLaw::Forced {correlation:CorrelationId::FlatPlateLaminarAverage,
            speed_m_s:1.0,settings:ForcedConvectionSettings::default()}})
}
fn h_reference(scale:f64,wall:f64)->f64 {
    let d=0.5*(wall+310.0)-300.0;let rho=1000.0-0.5*d;let cp=4000.0+2.0*d;
    let k=scale*(0.5+0.002*d);let mu=0.001*(1.0-0.01*d);
    0.664*(rho*0.1/mu).sqrt()*(mu*cp/k).cbrt()*k/0.1
}

#[test]
fn original_coolant_data_changes_heat_charge_rates_and_the_continuing_wet_frame() {
    for powered in [false,true] {
        fixture(powered,|m,plain,electric,circuit,heat,prepared,cx| {
            let gate=CancelGate::new_clock_free();let mut endpoints=Vec::new();
            for scale in [0.5,2.0] {
                let mut heat=heat.clone();let (placement,material,native)=native(m);
                let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
                if powered {owner.bind_material_circuit(&heat,electric,circuit.clone(),
                    &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],&gate).unwrap();}
                let card=coolant(scale);let law=coolant_law(&card);let boundaries=boundaries();
                let network=SourceConvectionHeatNetwork {heat:prepared,cells:&boundaries,maximum_cells:2,
                    resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t:Temperature|law.at_temperature(t)};
                let initial=owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                    +if powered {CHARGE*CHARGE/(2.0*CAP)}else{0.0};
                let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
                force.resize(owner.mechanical().mode_count(),0.0);
                let (mut work,mut outside,mut supplied_heat)=(0.0,0.0,0.0);
                for n in 1..=12 {
                    let cells=heat.cells().clone();let q0=owner.electrical().map_or(0.0,|c|c.charges_c()[0]);
                    let report=network.step_relaxation(cx,&mut owner,&mut heat,&force,
                        powered.then_some(input()),controls(),&[0.0;2],&gate).unwrap();
                    let sources=report.transport_source.as_ref().unwrap();
                    let mut external=0.0;let mut joule=[0.0;2];let mut stored=0.0;
                    if let Some(e)=&report.electrical {
                        let t=[heat.temperature_k(CellId(0)).unwrap(),heat.temperature_k(CellId(1)).unwrap()];
                        let r=[0.2*(1.0+0.02*(t[0]-300.0))/area(0),0.6*(1.0+0.02*(t[1]-300.0))/area(1)];
                        let a=DT/(2.0*CAP*(r[0]+r[1]));let q=q0*(1.0-a)/(1.0+a);
                        near(owner.electrical().unwrap().charges_c()[0],q,1e-12);
                        let current=(q0+q)/(2.0*CAP*(r[0]+r[1]));
                        for i in 0..2 {joule[i]=DT*current*current*r[i];near(e.resistors[i].heat_j,joule[i],1e-12);}
                        stored=e.stored_after_j;work+=e.native.supplied;
                    }
                    for (i,(fluid,sample)) in sources.materials().iter().zip(sources.convection().samples()).enumerate() {
                        let id=CellId(i as u64);let t=heat.temperature_k(id).unwrap();
                        near(fluid.property_temperature_k,0.5*(t+310.0),1e-12);
                        assert_eq!(fluid.material.material().card_identity(),card.content_hash());
                        assert_eq!(fluid.material.material().query().axes()["p"],100000.0);
                        near(sample.wall_temperature_k,t,0.0);
                        near(sample.reference_temperature_k,310.0,0.0);
                        near(sample.exposed_area_m2,exposed_area(i),1e-12);
                        let power=exposed_area(i)*h_reference(scale,t)*(310.0-t);
                        near(sample.heat_into_cell_w,power,1e-8);
                        let loss=report.layers.iter().find(|l|l.destination.cell==id).unwrap().heat_j;
                        let cell=&heat.cells()[&id];
                        near(cell.mass_kg*(cell.specific_h_j_kg-cells[&id].specific_h_j_kg),loss+joule[i]+DT*power,2e-8);
                        external+=DT*power;
                    }
                    near(report.external_thermal_heat_j,external,1e-8);
                    assert!(report.entropy.has_unclassified_sources); // Actual reservoir powers remain explicit.
                    work+=report.mechanical.external_work_j;
                    supplied_heat+=external;
                    outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                    near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+stored-initial+outside-work-supplied_heat,0.0,2e-7);
                    assert!(!report.rates.changes.is_empty());
                    assert_eq!(heat.time_s(),n as f64*DT);assert_eq!(owner.mechanical().samples_rendered(),n);
                }
                // The unselected constrained material history is never heated or reset by this binding.
                near(owner.mechanical().material_memory()[0],0.007*(-12.0*DT/0.7).exp(),1e-10);
                endpoints.push([placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                    heat.temperature_k(CellId(1)).unwrap(),owner.electrical().map_or(0.0,|c|c.charges_c()[0])]);
            }
            assert!((endpoints[1][0]-endpoints[0][0]).abs()>1e-12);
            assert!(endpoints[1][1]>endpoints[0][1]);
            if powered {assert!((endpoints[1][2]-endpoints[0][2]).abs()>1e-10);}
        });
    }
}

#[test]
fn empty_convection_is_the_same_powered_or_unpowered_heat_trajectory() {
    for powered in [false,true] {
        fixture(powered,|m,plain,electric,circuit,heat,prepared,cx| {
            let gate=CancelGate::new_clock_free();let mut endpoints=Vec::new();
            for composed in [false,true] {
                let mut heat=heat.clone();let (placement,material,native)=native(m);
                let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
                if powered {owner.bind_material_circuit(&heat,electric,circuit.clone(),
                    &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],&gate).unwrap();}
                let network=SourceConvectionHeatNetwork {heat:prepared,cells:&[],maximum_cells:0,
                    resolve_fluid:|_:&SourceConvectiveCell,_:&DepositionSimulation<'_>,_:Temperature|
                        ->Result<ResolvedLiquidState,MaterialStatePointError>{panic!("empty boundary invoked a fluid source")}};
                let force=vec![0.0;owner.mechanical().mode_count()];
                for _ in 0..3 {
                    if composed {let report=network.step_relaxation(cx,&mut owner,&mut heat,&force,
                        powered.then_some(input()),controls(),&[0.0;2],&gate).unwrap();
                        assert!(report.transport_source.unwrap().materials().is_empty());}
                    else {owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network:prepared,cx},
                        powered.then_some(input()),controls(),&[0.0;2],&gate).unwrap();}
                }
                endpoints.push((heat.cells().clone(),owner.electrical().cloned(),
                    owner.mechanical().material_memory().to_vec(),
                    placement.physical_state(owner.mechanical().components(),&gate).unwrap()));
            }
            assert_eq!(endpoints[0],endpoints[1]);
        });
    }
}

#[test]
fn late_coolant_failure_keeps_its_typed_cause_and_rolls_back_the_complete_assembly() {
    #[derive(Debug)]
    enum Failure { Domain(CellId), Original(MaterialStatePointError) }
    fixture(true,|m,plain,electric,circuit,mut heat,prepared,cx| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        owner.bind_material_circuit(&heat,electric,circuit.clone(),
            &[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)],&gate).unwrap();
        let card=coolant(1.0);let law=coolant_law(&card);let boundaries=boundaries();
        let calls=std::cell::Cell::new(0usize);let refuse=std::cell::Cell::new(true);
        let network=SourceConvectionHeatNetwork {heat:prepared,cells:&boundaries,maximum_cells:2,
            resolve_fluid:|binding:&SourceConvectiveCell,_:&DepositionSimulation<'_>,t:Temperature| {
                calls.set(calls.get()+1);
                if refuse.get() && calls.get()>2 {return Err(Failure::Domain(binding.cell));}
                law.at_temperature(t).map_err(Failure::Original)
            }};
        let cells=heat.cells().clone();let memories=owner.mechanical().material_memory().to_vec();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let material=owner.material().members()[0].source.clone();
        let force=vec![0.0;owner.mechanical().mode_count()];
        let error=network.step_relaxation(cx,&mut owner,&mut heat,&force,Some(input()),controls(),&[0.0;2],&gate);
        assert!(matches!(error,Err(RelaxationTransportError::Source(ConvectionHeatError::Convection(
            SourceCellConvectionError::Source {cell:CellId(0),source:Failure::Domain(CellId(0)),..})))));
        assert!(calls.get()>2);assert_eq!(heat.cells(),&cells);assert_eq!(owner.electrical().unwrap(),&circuit);
        assert_eq!(owner.mechanical().material_memory(),memories);
        assert_eq!(owner.material().members()[0].source,material);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        assert_eq!(owner.mechanical().samples_rendered(),0);assert_eq!(heat.time_s(),0.0);
        refuse.set(false);
        let mut exhausted=controls();exhausted.maximum_material_trials=1;
        assert!(network.step_relaxation(cx,&mut owner,&mut heat,&force,Some(input()),exhausted,&[0.0;2],&gate).is_err());
        assert_eq!(heat.cells(),&cells);assert_eq!(owner.electrical().unwrap(),&circuit);
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(matches!(network.step_relaxation(cx,&mut owner,&mut heat,&force,Some(input()),controls(),&[0.0;2],&cancelled),
            Err(RelaxationTransportError::Relaxation(CircuitEnthalpyError::Cancelled))));
        network.step_relaxation(cx,&mut owner,&mut heat,&force,Some(input()),controls(),&[0.0;2],&gate).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

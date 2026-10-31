//! Shared native motion, separately sourced layer heat, and cross-body conduction.
use super::*;
use crate::circuit_enthalpy::material::ConductionLink;
use crate::render::schedule::force::coupled::{ModalAttachment,ModalConnection};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;

fn conducting_source(layer:usize,varying:bool,k:f64)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let (original,_)=source(layer,varying);let mut claims=ClaimSet::new();
    for (_,claim) in original.claims().claims_ordered() {claims.insert_claim(claim.clone()).unwrap();}
    claims.insert_claim(PropertyClaim {key:PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),
        value:PropertyValue::Scalar {value:k,dims:THERMAL_CONDUCTIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),290.0,340.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"synthetic cross-body heat path".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(original.id().clone(),claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower:&point(290.0),upper:&point(340.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,curve)
}
fn destinations(body:usize)->Vec<FrameLayerCell> {
    (0..2).map(|layer|FrameLayerCell {element:0,layer,cell:CellId((2*body+layer) as u64)}).collect()
}
fn shared_budget()->MaterialHeatBudget {MaterialHeatBudget {max_groups:8,..budget()}}
fn shared_controls()->FrameThermalSolve {FrameThermalSolve {maximum_cells:4,maximum_material_trials:256,..controls()}}

fn coupled<'m,'f>(models:[&'m MaterialFluidFrame<'f,FluidTransportState>;2],connected:bool)
    ->(Vec<FrameFluidPlacement<'m,'f,FluidTransportState>>,Vec<FrameRelaxationBinding<'m,'f,FluidTransportState>>,CoupledModalSystem) {
    let gate=CancelGate::new_clock_free();let mut components=Vec::new();let mut connections=Vec::new();let mut placements=Vec::new();
    let cfg=ModalCouplingConfig {max_modes:16,max_connections:32,max_setup_terms:1_000_000,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9};
    for (body,model) in models.into_iter().enumerate() {
        let state=model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:if body==0 {0.0003}else{0.0},
            velocity_m_sqrt_kg_per_s:if body==0 {-0.001}else{0.0}}],&[0.0;4]).unwrap();
        let (placement,c,l)=model.coupled_component(&state,256,ModalAcousticTimeBudget::audible_reference(),cfg,&gate).unwrap()
            .into_parts(components.len(),connections.len(),&gate).unwrap();
        components.extend(c);connections.extend(l);placements.push(placement);
    }
    if connected {
        let point=[0.0,0.0,0.0,0.0,1.0,0.0];
        connections.push(ModalConnection {left:ModalAttachment {component:placements[0].structural_component(),
            shapes:placements[0].project_nodal_force(&point,&gate).unwrap()},
            right:ModalAttachment {component:placements[1].structural_component(),
                shapes:placements[1].project_nodal_force(&point,&gate).unwrap()},
            stiffness_n_m:4.0,damping_n_s_m:0.02,rest_extension_m:0.0});
    }
    let mut mechanical=CoupledModalSystem::new(components,connections,cfg,&gate).unwrap();
    let extent=2*mechanical.mode_count();
    let foreign=||vec![RelaxationBranch {projection:vec![0.0;extent],stiffness:2.0,relaxation_time_s:0.7}];
    let selected=[FrameLongitudinalRelaxation {element:0,layers:vec![LayerRelaxationProperties {equilibrium_modulus:"equilibrium".into(),
        branches:vec![("increment".into(),"time".into())]};2]}];
    let mut bindings=Vec::new();
    for (body,placement) in placements.iter().enumerate() {
        mechanical=mechanical.with_additional_material_relaxation(foreign(),Some(vec![0.007]),&gate).unwrap();
        let z=if body==0 {vec![0.0002,0.003,-0.001,0.0005,-0.002,0.0007]}else{vec![0.0;6]};
        let (binding,next)=placement.with_longitudinal_relaxation(mechanical,&selected,Some(z),
            FrameRelaxationBudget {max_members:1,max_source_arms:2,max_memory_coordinates:6,max_setup_terms:1_000_000},&gate).unwrap();
        bindings.push(binding);mechanical=next;
    }
    mechanical=mechanical.with_additional_material_relaxation(foreign(),Some(vec![0.007]),&gate).unwrap();
    (placements,bindings,mechanical)
}

fn shared_fixture<R>(varying:bool,receiver_conductivity:Option<f64>,
    run:impl FnOnce([&MaterialFluidFrame<'_,FluidTransportState>;2],&MaterialNetwork<'_>,DepositionSimulation<'_>)->R)->R {
    let (a,ac)=conducting_source(0,varying,1e-4);let (b,bc)=conducting_source(1,varying,1e-4);
    let (c,cc)=conducting_source(0,varying,receiver_conductivity.unwrap_or(1e-4));let (d,dc)=conducting_source(1,varying,1e-4);
    let cards=[[a,b],[c,d]];let curves=[ac,bc,cc,dc];let frames=[frame(&cards[0]),frame(&cards[1])];
    let models=[model(&frames[0]),model(&frames[1])];let query=point(300.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=core::array::from_fn::<_,4,_>(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i/2][i%2],context:&query,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let links=if receiver_conductivity.is_some() {vec![ConductionLink {link:17,distance_a_m:0.005,distance_b_m:0.005}]}else{vec![]};
    let source=MaterialNetwork {cells:&cells,resistors:&[],links:&links};let mut events=Vec::new();
    for i in 0..4 {
        let layer=i%2;let inner=if layer==0 {0.0_f64}else{RADII[layer-1]};
        let area=core::f64::consts::PI*(RADII[layer].powi(2)-inner.powi(2));
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {id:CellId(i as u64),segment_id:0,body_id:(i/2) as u64,
            material_id:MaterialId(i as u64),mass_kg:[3.0,5.0][layer]*area,specific_h_j_kg:10.0*CP,
            geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area,surface_area_m2:1.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});
    }
    if !links.is_empty() {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:17,a:CellId(0),b:CellId(2),area_m2:1e-3,conductance_w_k:999.0})});}
    let laws:Vec<_>=curves.iter().enumerate().map(|(i,c)|(MaterialId(i as u64),c.curve() as &dyn ThermalLaw)).collect();
    let thermal=DepositionSimulation::new(events,&laws,DepositionBudget {max_cells:4,max_links:1,max_events:5,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:1.0/256.0,energy_tolerance_j:1e-11},||Ok(())).unwrap().0;
    run([&models[0],&models[1]],&source,thermal)
}
fn bind<'m,'f>(materials:Vec<FrameRelaxationBinding<'m,'f,FluidTransportState>>,reverse:bool)
    ->Vec<(FrameRelaxationBinding<'m,'f,FluidTransportState>,Vec<FrameLayerCell>)> {
    let mut inputs:Vec<_>=materials.into_iter().enumerate().map(|(body,m)|(m,destinations(body))).collect();
    if reverse {inputs.reverse();}inputs
}

#[test]
fn connected_frames_heat_their_own_layers_and_preserve_foreign_prefix_gap_and_suffix() {
    let mut endpoints=Vec::new();
    for reverse in [false,true] {
        shared_fixture(true,None,|models,source,mut heat| {
            let gate=CancelGate::new_clock_free();let (placements,materials,mechanical)=coupled(models,true);
            let initial=mechanical.total_energy_j().unwrap()+caloric(&heat);
            let mut owner=FrameRelaxationThermalSystem::new_shared(bind(materials,reverse),mechanical,&heat,source,shared_budget(),&gate).unwrap();
            assert_eq!(owner.materials().len(),2);assert_eq!(owner.material_layer_ranges(),&[0..2,2..4]);
            let foreign:Vec<_>=owner.outside_groups.iter().map(|&i|owner.groups[i].clone()).collect();
            assert_eq!(foreign,vec![0..1,7..8,14..15]);
            let mut force=placements[0].project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();force.resize(owner.mechanical().mode_count(),0.0);
            let (mut work,mut outside)=(0.0,0.0);let mut received=[0.0;4];
            for step in 1..=16 {
                let report=owner.step(&mut heat,&force,shared_controls(),&[0.0;4],&gate).unwrap();
                assert_eq!(report.material_layer_ranges,vec![0..2,2..4]);
                assert_eq!(report.layers[0].destination.cell,CellId(if reverse {2}else{0}));
                for layer in report.layers {assert!(layer.heat_j>=0.0);received[layer.destination.cell.0 as usize]+=layer.heat_j;}
                work+=report.mechanical.external_work_j;outside+=report.external_mechanical_loss_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)-initial+outside-work,0.0);
                for i in [0,7,14] {near(owner.mechanical().material_memory()[i],0.007*(-f64::from(step)/256.0/0.7).exp());}
                for material in owner.materials() {material.states(owner.mechanical(),&gate).unwrap();}
            }
            assert!(received[0]+received[1]>0.0&&received[2]+received[3]>0.0);
            for i in 0..4 {assert!(heat.temperature_k(CellId(i)).unwrap()>300.0);}
            let receiver=placements[1].physical_state(owner.mechanical().components(),&gate).unwrap();assert!(receiver[0].abs()>1e-9);
            endpoints.push((receiver,heat.cells().values().map(|c|c.specific_h_j_kg).collect::<Vec<_>>()));
        });
    }
    for (a,b) in endpoints[0].0.iter().zip(&endpoints[1].0) {near(*a,*b);}
    for (a,b) in endpoints[0].1.iter().zip(&endpoints[1].1) {near(*a,*b);}
}

#[test]
fn original_cross_body_conductivity_changes_both_material_rates_and_driving_motion() {
    let mut ends=Vec::new();
    for k in [0.0,1e-4] {
        shared_fixture(true,Some(k),|models,source,mut heat| {
            let gate=CancelGate::new_clock_free();let (placements,materials,mechanical)=coupled(models,true);
            let mut owner=FrameRelaxationThermalSystem::new_shared(bind(materials,false),mechanical,&heat,source,shared_budget(),&gate).unwrap();
            let mut force=placements[0].project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();force.resize(owner.mechanical().mode_count(),0.0);
            for _ in 0..24 {
                let report=owner.step(&mut heat,&force,shared_controls(),&[0.0;4],&gate).unwrap();
                near(report.source.links[0].conductance_w_k,if k==0.0 {0.0}else{1e-5});
                near(report.external_thermal_heat_j,0.0);
                near(report.total_energy_residual_j,0.0);
            }
            ends.push((heat.temperature_k(CellId(0)).unwrap(),heat.temperature_k(CellId(2)).unwrap(),
                placements[0].physical_state(owner.mechanical().components(),&gate).unwrap()[0]));
        });
    }
    assert!(ends[1].0<ends[0].0&&ends[1].1>ends[0].1);
    assert!((ends[1].2-ends[0].2).abs()>1e-13);
}

#[test]
fn disconnected_unheated_receiver_stays_quiescent_without_resetting_its_foreign_history() {
    shared_fixture(true,None,|models,source,mut heat| {
        let gate=CancelGate::new_clock_free();let (placements,materials,mechanical)=coupled(models,false);
        let mut owner=FrameRelaxationThermalSystem::new_shared(bind(materials,false),mechanical,&heat,source,shared_budget(),&gate).unwrap();
        let force=vec![0.0;owner.mechanical().mode_count()];
        for _ in 0..8 {owner.step(&mut heat,&force,shared_controls(),&[0.0;4],&gate).unwrap();}
        assert!(placements[1].physical_state(owner.mechanical().components(),&gate).unwrap().iter().all(|x|*x==0.0));
        near(heat.temperature_k(CellId(2)).unwrap(),300.0);near(heat.temperature_k(CellId(3)).unwrap(),300.0);
        assert!(owner.mechanical().material_memory()[7]>0.0);
    });
}

#[test]
fn shared_late_refusal_keeps_both_bodies_rates_and_all_four_cells_for_retry() {
    shared_fixture(true,Some(1e-4),|models,source,mut heat| {
        let gate=CancelGate::new_clock_free();let (placements,materials,mechanical)=coupled(models,true);
        let mut owner=FrameRelaxationThermalSystem::new_shared(bind(materials,true),mechanical,&heat,source,shared_budget(),&gate).unwrap();
        let states:Vec<_>=placements.iter().map(|p|p.physical_state(owner.mechanical().components(),&gate).unwrap()).collect();
        let memories=owner.mechanical().material_memory().to_vec();let cells=heat.cells().clone();
        let sources:Vec<_>=owner.materials().iter().map(|b|b.members()[0].source.clone()).collect();
        let force=vec![0.0;owner.mechanical().mode_count()];let mut exhausted=shared_controls();exhausted.maximum_material_trials=1;
        assert!(owner.step(&mut heat,&force,exhausted,&[0.0;4],&gate).is_err());
        assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);assert_eq!(owner.mechanical().samples_rendered(),0);
        assert_eq!(owner.mechanical().material_memory(),memories);
        for (i,binding) in owner.materials().iter().enumerate() {assert_eq!(binding.members()[0].source,sources[i]);}
        for (i,p) in placements.iter().enumerate() {assert_eq!(p.physical_state(owner.mechanical().components(),&gate).unwrap(),states[i]);}
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(owner.step(&mut heat,&force,shared_controls(),&[0.0;4],&cancelled).is_err());
        owner.step(&mut heat,&force,shared_controls(),&[0.0;4],&gate).unwrap();assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

#[test]
fn duplicate_cross_body_inventory_and_incomplete_aggregate_heat_budget_refuse() {
    shared_fixture(true,None,|models,source,heat| {
        let gate=CancelGate::new_clock_free();
        for duplicate in [false,true] {
            let (_,materials,mechanical)=coupled(models,true);let mut inputs=bind(materials,false);
            if duplicate {inputs[1].1=destinations(0);}
            let limits=if duplicate {shared_budget()}else{MaterialHeatBudget {max_groups:6,..shared_budget()}};
            assert!(FrameRelaxationThermalSystem::new_shared(inputs,mechanical,&heat,source,limits,&gate).is_err());
        }
        let (_,materials,mechanical)=coupled(models,true);let mut inputs=bind(materials,false);
        let overlapping=inputs[0].0.memory_range();inputs[1].0.memory=overlapping;
        assert!(FrameRelaxationThermalSystem::new_shared(inputs,mechanical,&heat,source,shared_budget(),&gate).is_err());
        assert_eq!(heat.time_s(),0.0);
    });
}

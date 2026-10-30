use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution,frame::{MaterialFluidFrame,FrameFluidBinding,FrameFluidBudget}};
use crate::air_path::cylinder::frame::modal::FrameFluidPlacement;
use crate::air_path::cylinder::frame::modal::relaxation::{FrameLongitudinalRelaxation,FrameRelaxationBudget};
use crate::frame_vibration::FrameVibrationBudget;
use crate::modal_acoustic_time::{ModalAcousticState,ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use crate::circuit_enthalpy::material::MaterialCell;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
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
use fs_solid::frame2d::{Frame2d,FrameMassFormulation,Geometry};
use fs_thermal::enthalpy::*;
use fs_solver::NewtonKrylovConfig;

const CP:f64=0.0001; // Manufactured heat capacity, NOT measured material data.
const RADII:[f64;2]=[0.06,0.1];
fn point(t:f64)->QueryPoint {QueryPoint::new().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),t).unwrap().with("omega",1.0).unwrap()}
fn source(layer:usize,varying:bool)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let mut claims=ClaimSet::new();let tau=[0.025,0.12][layer];
    for (name,dims,value) in [("density",Density::DIMS,[3.0,5.0][layer]),
        ("young_modulus",Pressure::DIMS,[1000.0,1500.0][layer]),
        ("equilibrium",Pressure::DIMS,[1000.0,1500.0][layer]),
        ("increment",Pressure::DIMS,[600.0,900.0][layer]),("time",Time::DIMS,tau),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,CP)] {
        let curve=name=="time"&&varying;
        let mut validity=ValidityDomain::unconstrained().with("T",290.0,340.0);
        if name=="time" {validity=validity.with("omega",0.0,10000.0);}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if curve {PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(290.0,1.2*tau),(300.0,tau),(340.0,0.2*tau)],dims}}
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
            geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area,surface_area_m2:1.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})});
    }
    DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),(MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:2,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:1.0/256.0,energy_tolerance_j:1e-11},||Ok(())).unwrap().0
}
fn budget()->MaterialHeatBudget {MaterialHeatBudget {max_groups:4,max_dimension:256,max_dense_bytes:64*1024*1024,max_multiply_terms:100_000_000_000}}
fn controls()->FrameThermalSolve {FrameThermalSolve {maximum_cells:2,maximum_material_trials:128,maximum_newton_iterations:12,
    difference_step_j:1e-8,cell_energy_tolerance_j:1e-10,thermal_energy_tolerance_j:1e-10,total_energy_tolerance_j:1e-9,
    newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,..NewtonKrylovConfig::default()}}}
const DEST:[FrameLayerCell;2]=[FrameLayerCell {element:0,layer:0,cell:CellId(0)},FrameLayerCell {element:0,layer:1,cell:CellId(1)}];
fn fixture<R>(varying:bool,run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,&MaterialNetwork<'_>,DepositionSimulation<'_>)->R)->R {
    let (a,ac)=source(0,varying);let (b,bc)=source(1,varying);let cards=[a,b];let curves=[ac,bc];
    let frame=frame(&cards);let model=model(&frame);let query=point(300.0);
    let cells=core::array::from_fn::<_,2,_>(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context:&query,temperature_axis:"T",liquid_fraction_axis:None},conductivity:None});
    let source=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    run(&model,&source,thermal(&curves,1.0))
}
fn caloric(t:&DepositionSimulation<'_>)->f64 {t.cells().values().map(|c|c.mass_kg*c.specific_h_j_kg).sum()}
fn near(a:f64,b:f64){assert!((a-b).abs()<1e-9,"{a:e} != {b:e}");}

#[test]
fn own_layer_heat_changes_rates_and_motion_while_all_energy_and_foreign_history_survive() {
    let mut outcomes=Vec::new();
    for varying in [false,true] {
        fixture(varying,|model,source,mut heat| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(model);
            let initial=native.total_energy_j().unwrap()+caloric(&heat);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,source,&DEST,budget(),&gate).unwrap();
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();force.resize(owner.mechanical().mode_count(),0.0);
            let (mut work,mut outside,mut deposited)=(0.0,0.0,0.0);
            let mut changes=0;
            for _ in 0..12 {
                let report=owner.step(&mut heat,&force,controls(),&[0.0;2],&gate).unwrap();
                changes+=report.rates.changes.len();work+=report.mechanical.external_work_j;
                outside+=report.external_mechanical_loss_j;deposited+=report.deposited_material_heat_j;
                assert!(report.layers.iter().all(|l|l.heat_j>=0.0));
                near(report.layers.iter().map(|l|l.heat_j).sum(),report.deposited_material_heat_j);
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)-initial+outside-work,0.0);
                near(report.total_energy_residual_j,0.0);
            }
            assert!(deposited>0.0&&outside>0.0);
            assert_eq!(changes>0,varying);
            for i in 0..2 {assert!(heat.temperature_k(CellId(i)).unwrap()>300.0);}
            near(owner.mechanical().material_memory()[0],0.007*(-12.0_f64/256.0/0.7).exp());
            assert_eq!(owner.mechanical().samples_rendered(),12);
            outcomes.push(placement.physical_state(owner.mechanical().components(),&gate).unwrap());
        });
    }
    assert!((outcomes[0][0]-outcomes[1][0]).abs()>1e-12);
}

#[test]
fn late_caloric_refusal_and_cancellation_preserve_rates_cells_and_every_momentum_for_retry() {
    fixture(true,|model,source,mut heat| {
        let gate=CancelGate::new_clock_free();let (placement,material,native)=native(model);
        let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,source,&DEST,budget(),&gate).unwrap();
        let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let memories=owner.mechanical().material_memory().to_vec();let cells=heat.cells().clone();
        let original=owner.material().members()[0].source.clone();
        let force=vec![0.0;owner.mechanical().mode_count()];let mut exhausted=controls();exhausted.maximum_material_trials=1;
        assert!(owner.step(&mut heat,&force,exhausted,&[0.0;2],&gate).is_err());
        assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);
        assert_eq!(owner.mechanical().material_memory(),memories);
        assert_eq!(owner.material().members()[0].source,original);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(owner.step(&mut heat,&force,controls(),&[0.0;2],&cancelled).is_err());
        assert_eq!(owner.mechanical().samples_rendered(),0);
        owner.step(&mut heat,&force,controls(),&[0.0;2],&gate).unwrap();assert_eq!(owner.mechanical().samples_rendered(),1);
    });
}

#[test]
fn missing_or_duplicate_layer_inventory_is_not_replaced_by_arbitrary_heat_fractions() {
    fixture(true,|model,source,heat| {
        let gate=CancelGate::new_clock_free();
        for destinations in [vec![DEST[0]],vec![DEST[0],DEST[0]],vec![DEST[0],FrameLayerCell {cell:CellId(0),..DEST[1]}]] {
            let (_,material,native)=native(model);
            assert!(FrameRelaxationThermalSystem::new(material,native,&heat,source,&destinations,budget(),&gate).is_err());
        }
        assert_eq!(heat.time_s(),0.0);
    });
}

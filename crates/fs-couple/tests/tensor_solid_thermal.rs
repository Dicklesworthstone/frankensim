//! Synthetic Kelvin/Maxwell self-heating; no empirical material qualification.
use std::collections::BTreeMap;
use fs_blake3::ContentHash;
use fs_couple::tensor_material::{InitialMemory, TensorMaterialBudget, solid::TensorSolidBody};
use fs_couple::tensor_material::solid::thermal::{
    ElementThermalCell, SolidThermalConfig, SolidThermalError,
};
use fs_exec::{CancelGate, Cx};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::{TetAssemblyBudget, tet10::{TetQuadrature, promote_p1_to_p2},
    visco::{SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial,
        ViscoElements, ViscoStateBudget}};
use fs_thermal::enthalpy::*;

const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const FIXED: [usize;11] = [0,1,2,3,4,5,6,7,8,9,10];
const ZERO: [[f64;3];4] = [[0.;3];4];
const MASS: f64 = 1.0/3.0;
fn scope<R>(f: impl FnOnce(&Cx<'_>, &CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena|
        f(&Cx::new(&gate,arena,fs_exec::StreamKey {seed:87,kernel_id:87,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn diagonal(value:f64)->Mandel6 {
    core::array::from_fn(|i|core::array::from_fn(|j|if i==j {value} else {0.}))
}
fn material(kelvin_s:f64, maxwell:f64)->TetViscoMaterial {
    let branches=if maxwell==0. {vec![]} else {
        vec![PronyBranch::new(diagonal(120.*maxwell),0.08).unwrap()]
    };
    TetViscoMaterial {law:TensorProny::new(diagonal(120.),Some(diagonal(120.*kelvin_s)),
        branches,(0.,1e5)).unwrap(),density_kg_m3:2.,material_state_identity:[87;32]}
}
fn source(mesh:SolidMeshRef<'_>,field:TetViscoField<'_>,fixed:&[usize],cx:&Cx<'_>)->ViscoElements {
    ViscoElements::new(mesh,field,fixed,TetAssemblyBudget::standard(),
        ViscoStateBudget {max_internal_variables:10000},cx).unwrap()
}
fn body<'s>(source:&'s ViscoElements,basis:&[Vec<f64>],cx:&Cx<'_>,gate:&CancelGate)->TensorSolidBody<'s> {
    let projection=source.project_basis(basis,SolidProjectionBudget {max_coordinates:8,
        max_basis_scalars:4096,max_material_points:64,max_scalar_terms:1_000_000},cx).unwrap();
    TensorSolidBody::new(projection,TensorMaterialBudget {max_state_coordinates:112,
        max_material_points:64},gate).unwrap()
}
fn curve(cp:f64)->EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([87;32]),vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.,temperature_k:300.,liquid_mass_fraction:0.,bulk_density_kg_m3:2.},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:100.*cp,temperature_k:400.,liquid_mass_fraction:0.,bulk_density_kg_m3:2.},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:100.*cp+1.,temperature_k:400.,liquid_mass_fraction:1.,bulk_density_kg_m3:2.},
    ]).unwrap()
}
fn boundary(source_w:f64)->Boundary {
    Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w}
}
fn thermal<'a>(curve:&'a EquilibriumEnthalpyPhaseCurve,ids:&[u64],extra:Vec<TimedEvent>,steps:usize)
    ->DepositionSimulation<'a> {
    let mut events:Vec<_>=ids.iter().map(|&id|TimedEvent {time_s:0.,event:Event::Birth(DepositedCell {
        id:CellId(id),segment_id:id,body_id:1,material_id:MaterialId(0),mass_kg:MASS,
        specific_h_j_kg:0.,geometry:CellGeometry {start_m:[0.;3],end_m:[1.,0.,0.],
            volume_m3:1./6.,surface_area_m2:3.},boundary:boundary(0.),
    })}).collect();
    events.extend(extra);
    DepositionSimulation::new(events,&[(MaterialId(0),curve as &dyn ThermalLaw)],
        DepositionBudget {max_cells:4,max_links:4,max_events:16,max_substeps:steps}).unwrap()
        .advance_to(0.,StepConfig {max_dt_s:0.03125,energy_tolerance_j:1e-11},||Ok(())).unwrap().0
}
fn cells<'a>(curve:&'a EquilibriumEnthalpyPhaseCurve,ids:&[u64])->Vec<ElementThermalCell<'a>> {
    ids.iter().map(|&id|ElementThermalCell {cell:CellId(id),curve,
        temperature_interval_k:[300.,350.]}).collect()
}
fn controls()->SolidThermalConfig {
    SolidThermalConfig {thermal:StepConfig {max_dt_s:0.03125,energy_tolerance_j:1e-11},energy_tolerance_j:1e-10}
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a:e} != {b:e}");}

#[test]
fn kelvin_motion_heats_the_actual_caloric_cell_and_continues_without_reset() {
    scope(|cx,gate| {
        let mat=material(0.01,0.);
        let source=source(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},TetViscoField::Uniform(&mat),&FIXED,cx);
        let body=body(&source,&[vec![1.]],cx,gate);
        let initial=body.history(&[0.001],&[0.1],InitialMemory::Relaxed,gate).unwrap();
        let curve=curve(1.);
        let heat=thermal(&curve,&[30],vec![],16);
        let bindings=cells(&curve,&[30]);
        let mut native=initial.clone(); let mut caloric=heat.clone();
        let initial_energy=0.5*(1./30.)*0.1_f64.powi(2)+0.5*20.*0.001_f64.powi(2);
        for target in [0.125,0.25] {
            let candidate=native.advance_thermal_candidate(&caloric,&(),target,&ZERO,&bindings,controls(),gate).unwrap();
            let w=591.0_f64.sqrt(); let decay=(-3.*target).exp();
            let shape=0.001*(w*target).cos()+0.103/w*(w*target).sin();
            let q=decay*shape;
            let v=decay*(-3.*shape-0.001*w*(w*target).sin()+0.103*(w*target).cos());
            let energy=0.5*(1./30.)*v*v+0.5*20.*q*q;
            near(candidate.solid().material_history().displacement_m()[0],q,1e-11);
            near(candidate.solid().material_history().velocity_m_s()[0],v,1e-10);
            near(candidate.thermal().temperature_k(CellId(30)).unwrap(),300.+(initial_energy-energy)/MASS,1e-8);
            near(candidate.thermal().energy_j()+candidate.solid().material_history().stored_energy_j(),initial_energy,1e-10);
            assert_eq!(candidate.report().thermal.substeps,4);
            near(candidate.report().external_thermal_heat_j,0.,1e-12);
            near(candidate.report().combined_residual_j,0.,1e-10);
            let (next,thermal,_,_)=candidate.into_parts(); native=next; caloric=thermal;
        }
        assert_eq!(native.material_history().steps(),8);
        assert_eq!(initial.material_history().time_s(),0.);
        assert_eq!(heat.cells()[&CellId(30)].specific_h_j_kg,0.);
    });
}

#[test]
fn viscosity_and_caloric_data_change_heat_and_temperature_without_substituted_physics() {
    scope(|cx,gate| {
        let mut results=Vec::new();
        for (viscosity,cp) in [(0.01,1.),(0.03,1.),(0.01,2.),(0.,1.)] {
            let mat=material(viscosity,0.);
            let source=source(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},TetViscoField::Uniform(&mat),&FIXED,cx);
            let body=body(&source,&[vec![1.]],cx,gate);
            let history=body.history(&[0.001],&[0.1],InitialMemory::Relaxed,gate).unwrap();
            let chart=curve(cp); let heat=thermal(&chart,&[30],vec![],16);
            let next=history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&chart,&[30]),controls(),gate).unwrap();
            results.push([next.solid().material_history().displacement_m()[0],
                next.thermal().temperature_k(CellId(30)).unwrap()-300.,next.report().material_heat_j]);
        }
        assert!((results[0][0]-results[1][0]).abs()>1e-6);
        assert!((results[0][1]-results[1][1]).abs()>1e-6);
        assert_eq!(results[0][0],results[2][0]);
        assert_eq!(results[0][2],results[2][2]);
        near(results[0][1],2.*results[2][1],1e-10);
        assert_eq!(results[3][1],0.); near(results[3][2],0.,1e-14);
    });
}

#[derive(Clone,Default)]
struct OtherPhysics {
    sources:usize, conductances:usize, transfers:usize, limits:usize, events:usize,
    intervals:Vec<(f64,f64)>, steps:usize, fail_after:Option<f64>,
}
impl DepositionObserver for OtherPhysics {
    type Error=DepositionError;
    fn sources(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {
        self.sources+=1; Ok(BTreeMap::from([(CellId(30),0.001)]))
    }
    fn conductances(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<u64,f64>,Self::Error> {
        self.conductances+=1; Ok(BTreeMap::from([(55,0.5)]))
    }
    fn heat_transfers(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<Vec<HeatTransfer>,Self::Error> {
        self.transfers+=1; Ok(vec![HeatTransfer {a:CellId(10),b:CellId(30),power_w:0.0002}])
    }
    fn limit_timestep_s(&mut self,_:&DepositionSimulation<'_>,dt:f64,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<f64,Self::Error> {self.limits+=1; Ok(dt.min(0.015625))}
    fn requires_step_interval(&self)->bool {true}
    fn step_sources(&mut self,_:&DepositionSimulation<'_>,start:f64,dt:f64,_:f64,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error> {
        self.intervals.push((start,dt)); Ok(BTreeMap::from([(CellId(10),0.002)]))
    }
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error> {self.events+=1;Ok(())}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,thermal:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        self.steps+=1;
        if self.fail_after.is_some_and(|t|thermal.time_s()>t) {Err(DepositionError::InvalidInput("synthetic late observer refusal"))}
        else {Ok(())}
    }
}

#[test]
fn regional_maxwell_heat_shares_conduction_events_and_every_existing_observer_hook() {
    scope(|cx,gate| {
        let mut nodes=NODES.to_vec();nodes.extend(NODES.map(|mut p|{p[0]+=2.;p}));
        let tets=[[0,1,2,3],[4,5,6,7]];
        let fixed:Vec<_>=(0..24).filter(|&d|d!=11 && d!=23).collect();
        let materials=[material(0.01,0.),material(0.,0.5)];
        let source=source(SolidMeshRef::Tet4 {nodes_m:&nodes,elements:&tets},TetViscoField::PerElement(&materials),&fixed,cx);
        let body=body(&source,&[vec![1.,0.],vec![0.,1.]],cx,gate);
        let initial=body.history(&[0.001,0.001],&[0.1,0.1],InitialMemory::Instantaneous,gate).unwrap();
        let chart=curve(1.); let switch=0.0234375;
        let heat=thermal(&chart,&[30,10],vec![
            TimedEvent {time_s:0.,event:Event::Contact(ThermalLink {id:55,a:CellId(30),b:CellId(10),
                area_m2:0.1,conductance_w_k:0.})},
            TimedEvent {time_s:switch,event:Event::Boundary(CellId(30),boundary(0.004))},
        ],16);
        let mut force=vec![[0.;3];8];force[3][2]=0.01;force[7][2]=-0.02;
        let observer=OtherPhysics::default();
        let next=initial.advance_thermal_candidate(&heat,&observer,0.0625,&force,&cells(&chart,&[30,10]),controls(),gate).unwrap();
        let obs=next.observer(); let count=next.report().thermal.substeps;
        assert_eq!(count,5); assert_eq!(obs.events,1);
        assert_eq!((obs.sources,obs.conductances,obs.transfers,obs.limits,obs.steps),(count,count,count,count,count));
        assert_eq!(obs.intervals.len(),count);assert!(obs.intervals.iter().any(|&(t,_)|t==switch));
        let mut direct=initial.clone(); let mut temperatures=[300.,300.]; let mut heats=[0.,0.];
        for &(start,dt) in &obs.intervals {
            let step=direct.advance(&force,dt,controls().energy_tolerance_j,gate).unwrap();
            let transfer=0.5*(temperatures[1]-temperatures[0])+0.0002;
            let external_a=0.001+if start>=switch {0.004} else {0.};
            temperatures[0]+=(step.element_heat_j[0]+dt*(external_a+transfer))/MASS;
            temperatures[1]+=(step.element_heat_j[1]+dt*(0.002-transfer))/MASS;
            for i in 0..2 {heats[i]+=step.element_heat_j[i];}
        }
        for (index,id) in [30,10].into_iter().enumerate() {
            near(next.thermal().temperature_k(CellId(id)).unwrap(),temperatures[index],1e-9);
            near(next.report().element_heat_j[index],heats[index],1e-11);
        }
        assert!((heats[0]-heats[1]).abs()>1e-8);
        assert_eq!(next.solid().material_history().state(),direct.material_history().state());
        near(next.report().external_thermal_heat_j,0.003*0.0625+0.004*(0.0625-switch),1e-12);
        near(next.report().combined_residual_j,0.,1e-10);
        assert!(observer.intervals.is_empty());assert_eq!(observer.steps,0);
        // All original authored boundaries/conductances survive the overrides.
        assert_eq!(next.thermal().cells()[&CellId(30)].boundary.source_w,0.004);
        assert_eq!(next.thermal().links()[&55].conductance_w_k,0.);
        let failed=OtherPhysics {fail_after:Some(0.02),..Default::default()};
        assert!(matches!(initial.advance_thermal_candidate(&heat,&failed,0.0625,&force,
            &cells(&chart,&[30,10]),controls(),gate),Err(SolidThermalError::Observer(_))));
        assert_eq!(initial.material_history().time_s(),0.);assert_eq!(failed.steps,0);
        assert_eq!(heat.cells()[&CellId(30)].specific_h_j_kg,0.);
    });
}

#[test]
fn late_temperature_phase_and_budget_refusals_do_not_consume_mechanical_or_thermal_history() {
    scope(|cx,gate| {
        let mat=material(0.01,0.5);
        let source=source(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},TetViscoField::Uniform(&mat),&FIXED,cx);
        let body=body(&source,&[vec![1.]],cx,gate);
        let history=body.history(&[0.001],&[0.1],InitialMemory::Instantaneous,gate).unwrap();
        let saved=history.material_history().state().to_vec();
        let chart=curve(1.);let heat=thermal(&chart,&[30],vec![],16);
        let mut limited=cells(&chart,&[30]);limited[0].temperature_interval_k=[300.,300.00000001];
        assert!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&limited,controls(),gate).is_err());
        let phase_chart=curve(1e-7); let phase_heat=thermal(&phase_chart,&[30],vec![],16);
        let mut phase_binding=cells(&phase_chart,&[30]);phase_binding[0].temperature_interval_k=[300.,450.];
        assert!(history.advance_thermal_candidate(&phase_heat,&(),0.125,&ZERO,&phase_binding,controls(),gate).is_err());
        let budget_heat=thermal(&chart,&[30],vec![],1);
        assert!(matches!(history.advance_thermal_candidate(&budget_heat,&(),0.125,&ZERO,&cells(&chart,&[30]),
            controls(),gate),Err(SolidThermalError::Thermal(DepositionError::Budget(_)))));
        assert_eq!(history.material_history().state(),saved.as_slice());
        assert_eq!(history.material_history().time_s(),0.);
        assert_eq!(history.material_history().dissipated_heat_j(),0.);
        assert_eq!(heat.cells()[&CellId(30)].specific_h_j_kg,0.);
        let retry=history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&chart,&[30]),controls(),gate).unwrap();
        assert!(retry.report().deposited_heat_j>0.);
        gate.request();
        assert!(matches!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&chart,&[30]),
            controls(),gate),Err(SolidThermalError::Cancelled)));
        assert_eq!(history.material_history().state(),saved.as_slice());
    });
}

#[test]
fn exact_caloric_binding_clock_and_complete_routes_are_required() {
    scope(|cx,gate| {
        let mat=material(0.01,0.);
        let source=source(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},TetViscoField::Uniform(&mat),&FIXED,cx);
        let body=body(&source,&[vec![1.]],cx,gate);
        let mut history=body.history(&[0.001],&[0.1],InitialMemory::Relaxed,gate).unwrap();
        let chart=curve(1.);let heat=thermal(&chart,&[30],vec![],16);let copy=chart.clone();
        assert!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&[],controls(),gate).is_err());
        assert!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&copy,&[30]),controls(),gate).is_err());
        assert!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&chart,&[31]),controls(),gate).is_err());
        history.advance(&ZERO,0.03125,1e-10,gate).unwrap();
        assert!(history.advance_thermal_candidate(&heat,&(),0.125,&ZERO,&cells(&chart,&[30]),controls(),gate).is_err());
        assert_eq!(heat.time_s(),0.);
    });
}

#[test]
fn tet10_quadrature_heat_reaches_one_original_element_without_duplication() {
    scope(|cx,gate| {
        let mat=material(0.01,0.);
        let (nodes,tets)=promote_p1_to_p2(&NODES,&TETS,None,cx).unwrap();
        let source=source(SolidMeshRef::Tet10 {nodes_m:&nodes,elements:&tets,
            quadrature:TetQuadrature {stiffness_points_per_dir:2,mass_points_per_dir:4}},
            TetViscoField::Uniform(&mat),&[],cx);
        let basis=vec![nodes.iter().flat_map(|p|[0.,0.,p[2]]).collect()];
        let body=body(&source,&basis,cx,gate);
        assert!(body.projection().material_points().len()>1);
        let history=body.history(&[0.001],&[0.1],InitialMemory::Relaxed,gate).unwrap();
        let chart=curve(1.);let heat=thermal(&chart,&[30],vec![],16);
        let next=history.advance_thermal_candidate(&heat,&(),0.125,&vec![[0.;3];nodes.len()],
            &cells(&chart,&[30]),controls(),gate).unwrap();
        assert_eq!(next.report().element_heat_j.len(),1);
        near(next.report().element_heat_j[0],next.report().thermal.stored_energy_change_j,1e-10);
        near(next.report().mechanical_energy_change_j+next.report().deposited_heat_j,0.,1e-10);
    });
}

//! Independent annular momentum and regional enthalpy balances.

use super::*;
use crate::graph::rheology::{CarreauPipeBranch, CarreauPipeGraph};
use crate::graph::rheology::transient::LinearHydraulicStorage;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::rheology::pipe::LaminarCarreauPipe;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::{Density, VolumetricFlowRate};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::collections::BTreeMap;

const R: f64 = 0.01;
const L: f64 = 0.5;
const RHO: f64 = 1000.0;
const CP: f64 = 500.0;

struct Caloric;
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8;32]> { Some([67;32]) }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t=300.0+h/CP;
        if t.is_finite() && (280.0..=360.0).contains(&t) { Ok(t) }
        else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0/CP }
}
fn inventory<'a>(law: &'a Caloric, temperatures: &[f64]) -> DepositionSimulation<'a> {
    // Two explicit junction inventories, followed by independent real parallel
    // pipe volumes. Native stable CellId order is not a transport surrogate.
    let events=temperatures.iter().enumerate().map(|(i,&t)| {
        let (id,mass,start,end)=match i {
            0 => (CellId(0),0.04,[-0.1,0.0,0.0],[0.0;3]),
            1 => (CellId(1),0.03,[L,0.0,0.0],[L+0.1,0.0,0.0]),
            _ => (CellId(18+i as u64),RHO*core::f64::consts::PI*R*R*L,[0.0;3],[L,0.0,0.0]),
        };
        TimedEvent { time_s:0.0, event:Event::Birth(DepositedCell {
            id,segment_id:i as u64,body_id:1,material_id:MaterialId(1),mass_kg:mass,
            specific_h_j_kg:CP*(t-300.0),
            geometry:CellGeometry { start_m:start,end_m:end,volume_m3:mass/RHO,surface_area_m2:1.0 },
            boundary:Boundary { ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0 },
        }) }
    }).collect();
    DepositionSimulation::new(events,&[(MaterialId(1),law as &dyn ThermalLaw)],
        DepositionBudget { max_cells:8,max_links:8,max_events:10,max_substeps:100 }).unwrap()
        .advance_to(0.0,StepConfig { max_dt_s:1.0,energy_tolerance_j:1e-9 },|| Ok(())).unwrap().0
}
fn parts(cx:&Cx<'_>,thermal:&DepositionSimulation<'_>,parallel:usize,shear:bool)
    -> (RegionalCarreauTransport,TransientCarreauGraph) {
    let law=ThermalCarreauYasuda::new(
        CarreauYasuda::new(0.2,if shear {0.04}else{0.2},0.02,0.7,1.5,[0.0,1e7]).unwrap(),
        TemperatureShift::new(TemperatureShiftModel::Arrhenius { activation_energy_j_mol:9000.0 },
            300.0,[280.0,360.0]).unwrap(),ViscosityShiftMode::VerticalAndRate);
    let pipe=LaminarCarreauPipe::new(law,300.0,R,L,1e5).unwrap();
    let graph=CarreauPipeGraph::new(2,(0..parallel).map(|i| CarreauPipeBranch {
        name:format!("pipe-{i}"),from:0,to:1,pipe,
    }).collect()).unwrap();
    let cells:Vec<_>=(0..parallel).map(|i|CellId(20+i as u64)).collect();
    let regional=RegionalCarreauTransport::new(cx,graph,Density::new(RHO),
        &[CellId(0),CellId(1)],&cells,thermal,16).unwrap();
    let momentum=TransientCarreauGraph::new(regional.graph().clone(),RHO,&vec![2;parallel],vec![None;2]).unwrap();
    (regional,momentum)
}
fn model(cx:&Cx<'_>,thermal:&DepositionSimulation<'_>,parallel:usize,shear:bool)->TransientRegionalCarreauTransport {
    let (regional,momentum)=parts(cx,thermal,parallel,shear);
    TransientRegionalCarreauTransport::new(cx,regional,momentum).unwrap()
}
fn controls()->TransientRegionalCarreauControls {
    let mut hydraulics=TransientCarreauControls {
        absolute_flow_tolerance:VolumetricFlowRate::new(1e-20),relative_flow_tolerance:2e-15,
        ..TransientCarreauControls::default()
    };
    hydraulics.radial.absolute_velocity_tolerance_m_s=1e-13;
    hydraulics.radial.relative_velocity_tolerance=1e-13;
    TransientRegionalCarreauControls {
        hydraulics,
        thermal:ImplicitTransportControls {
            solver:NewtonKrylovConfig { absolute_tolerance:1e-10,relative_tolerance:0.0,
                ..NewtonKrylovConfig::default() },maximum_cells:8,maximum_newton_iterations:30,
            maximum_endpoint_evaluations:1024,energy_scale_j:100.0,finite_difference_step:1e-5,
            cell_tolerance_j:1e-7,total_tolerance_j:3e-7,
        },maximum_total_branch_solves:10000,mechanical_energy_tolerance_j:1e-10,
    }
}
fn pressures(dp:f64)->[FixedPressure;2] {
    [FixedPressure { node:0,pressure:Pressure::new(dp) },FixedPressure { node:1,pressure:Pressure::new(0.0) }]
}
fn exchange(flow:&TransientCarreauStep,inlet:f64)->CarreauTransportExchange<()> {
    let mass_boundaries=flow.node_outflows.iter().enumerate().filter_map(|(node,flow)| {
        let kg_s=RHO*flow.value();
        if kg_s>0.0 { Some(MassBoundary::Inlet { id:100+node as u64,cell:CellId(node as u64),kg_s,
            material_id:MaterialId(1),specific_h_j_kg:CP*(inlet-300.0) }) }
        else if kg_s<0.0 { Some(MassBoundary::Outlet { id:100+node as u64,cell:CellId(node as u64),kg_s:-kg_s }) }
        else { None }
    }).collect();
    CarreauTransportExchange { mass_boundaries,
        heat:StepCoefficients { sources_w:BTreeMap::new(),conductances_w_k:BTreeMap::new(),payload:() }.into() }
}
fn near(a:f64,b:f64,absolute:f64,relative:f64) {
    assert!(a.is_finite() && b.is_finite() && (a-b).abs()<=absolute+relative*a.abs().max(b.abs()),"{a:.16e} != {b:.16e}");
}
fn snapshot(state:&TransientCarreauState)->(f64,Vec<Vec<f64>>) {
    (state.accepted_time_s(),state.branches().iter().map(|branch|branch.velocities_m_s().to_vec()).collect())
}

struct Reference {
    velocity:[f64;2],mean_flow:f64,heat:f64,kinetic:f64,temperatures:[f64;3],
}
fn mu(t:f64)->f64 { 0.2*(9000.0/8.31446261815324*(1.0/t-1.0/300.0)).exp() }
// Independent two-annulus FV solve: masses are exact annular volumes. The
// internal face is at R/2; no-slip wall distance is half a radial cell.
fn momentum_reference(old:[f64;2],dt:f64,dp:f64,t:f64)->([f64;2],f64,f64,f64) {
    let area=core::f64::consts::PI*R*R;
    let a0=area/4.0; let a1=3.0*area/4.0;
    let m0=RHO*L*a0; let m1=RHO*L*a1;
    let k=2.0*core::f64::consts::PI*mu(t)*L; let wall=4.0*k;
    let aa=m0+0.5*dt*k; let bb=-0.5*dt*k; let dd=m1+0.5*dt*(k+wall);
    let f0=(m0-0.5*dt*k)*old[0]+0.5*dt*k*old[1]+dt*dp*a0;
    let f1=0.5*dt*k*old[0]+(m1-0.5*dt*(k+wall))*old[1]+dt*dp*a1;
    let determinant=aa*dd-bb*bb;
    let next=[(dd*f0-bb*f1)/determinant,(aa*f1-bb*f0)/determinant];
    let mid=[0.5*(old[0]+next[0]),0.5*(old[1]+next[1])];
    let q=a0*mid[0]+a1*mid[1];
    let heat=dt*(k*(mid[1]-mid[0]).powi(2)+wall*mid[1]*mid[1]);
    let kinetic=m0*mid[0]*(next[0]-old[0])+m1*mid[1]*(next[1]-old[1]);
    (next,q,heat,kinetic)
}
fn reference(old_velocity:[f64;2],old:&DepositionSimulation<'_>,dt:f64,dp:f64,inlet:f64)->Reference {
    let ids=[CellId(0),CellId(20),CellId(1)];
    let old_t:Vec<_>=ids.iter().map(|&id|old.temperature_k(id).unwrap()).collect();
    let masses:Vec<_>=ids.iter().map(|id|old.cells()[id].mass_kg).collect();
    let residual=|pipe_t:f64| {
        let (_,q,heat,_)=momentum_reference(old_velocity,dt,dp,pipe_t);
        let upstream=if q>=0.0 {0}else{2}; let moved=dt*RHO*q.abs();
        let t=(masses[upstream]*old_t[upstream]+moved*inlet)/(masses[upstream]+moved);
        masses[1]*CP*(pipe_t-old_t[1])-moved*CP*(t-pipe_t)-heat
    };
    let (mut lo,mut hi)=(280.0,360.0);
    assert!(residual(lo)<0.0 && residual(hi)>0.0);
    for _ in 0..90 { let mid=0.5*(lo+hi); if residual(mid)>0.0 {hi=mid}else{lo=mid} }
    let pipe_t=0.5*(lo+hi);
    let (velocity,q,heat,kinetic)=momentum_reference(old_velocity,dt,dp,pipe_t);
    let (up,down)=if q>=0.0 {(0,2)}else{(2,0)}; let moved=dt*RHO*q.abs();
    let mut temperatures=[0.0;3]; temperatures[1]=pipe_t;
    temperatures[up]=(masses[up]*old_t[up]+moved*inlet)/(masses[up]+moved);
    temperatures[down]=(masses[down]*old_t[down]+moved*pipe_t)/(masses[down]+moved);
    Reference {velocity,mean_flow:q,heat,kinetic,temperatures}
}
fn check(step:&TransientRegionalCarreauStep<'_,()>,expected:&Reference) {
    for (&actual,&expected) in step.momentum.branches()[0].velocities_m_s().iter().zip(&expected.velocity) {
        near(actual,expected,1e-11,2e-8);
    }
    for (&id,&t) in [CellId(0),CellId(20),CellId(1)].iter().zip(&expected.temperatures) {
        near(step.thermal.temperature_k(id).unwrap(),t,2e-7,0.0);
    }
    let flow=&step.report.hydraulics.branches[0];
    near(flow.midpoint_volume_flow_m3_s,expected.mean_flow,1e-14,2e-8);
    near(flow.viscous_heat_j,expected.heat,1e-11,2e-8);
    near(flow.kinetic_energy_change_j,expected.kinetic,1e-11,2e-8);
    near(flow.pressure_work_j,expected.heat+expected.kinetic,1e-11,2e-8);
    near(step.next.momentum().pipes()[0].pipe().temperature_k(),step.thermal.temperature_k(CellId(20)).unwrap(),0.0,0.0);
    assert!(step.report.combined_energy_residual_j.abs()<=step.report.combined_energy_tolerance_j);
    near(step.report.thermal.thermal.stored_energy_change_j,
        step.report.external_stream_enthalpy_j+step.report.other_heat_j+expected.heat,3e-7,1e-9);
    assert_eq!(step.momentum.accepted_time_s(),step.thermal.time_s());
}

#[test]
fn g1_startup_coastdown_and_recovering_head_match_annular_and_advected_energy_oracles() {
    with_cx(|cx| {
        let caloric=Caloric;
        let initial=inventory(&caloric,&[300.0,295.0,310.0]);
        let owner=model(cx,&initial,1,false);
        let state=owner.momentum().rest_state(&[None;2],cx).unwrap();
        let before=snapshot(&state); let cells=initial.cells().clone();
        let mut work=owner.momentum().workspace().unwrap();
        let expected=reference([0.0;2],&initial,0.04,5000.0,340.0);
        let startup=owner.advance(cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,340.0))).unwrap();
        check(&startup,&expected);
        assert!(startup.report.hydraulics.kinetic_energy_change_j>startup.report.hydraulics.viscous_heat_j);
        assert_eq!(snapshot(&state),before); assert_eq!(initial.cells(),&cells);
        let expected=reference(expected.velocity,&startup.thermal,0.03,0.0,300.0);
        let coast=startup.next.advance(cx,&startup.thermal,&startup.momentum,&mut work,0.07,&pressures(0.0),&[],controls(),&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,300.0))).unwrap();
        check(&coast,&expected);
        assert_eq!(coast.report.hydraulics.reservoir_work_j,0.0);
        assert!(coast.report.hydraulics.viscous_heat_j>0.0 && coast.report.hydraulics.kinetic_energy_change_j<0.0);
        let expected=reference(expected.velocity,&coast.thermal,0.01,-100.0,300.0);
        let recovered=coast.next.advance_with_heads(cx,&coast.thermal,&coast.momentum,&mut work,0.08,&pressures(0.0),&[],
            &[Pressure::new(-100.0)],controls(),&[0.0;3],|_,_,flow|Ok::<_,()>(exchange(flow,300.0))).unwrap();
        check(&recovered,&expected);
        assert!(recovered.report.hydraulics.branch_head_work_j<0.0);
        assert!(recovered.report.hydraulics.viscous_heat_j>0.0);
    });
}

#[test]
fn g3_original_nonlinear_branches_keep_distinct_transported_temperatures_and_radial_histories() {
    with_cx(|cx| {
        let caloric=Caloric; let initial=inventory(&caloric,&[300.0,300.0,300.0,320.0]);
        let owner=model(cx,&initial,2,true);
        let state=owner.momentum().rest_state(&[None;2],cx).unwrap();
        let mut work=owner.momentum().workspace().unwrap();
        let actual=owner.advance(cx,&initial,&state,&mut work,0.1,&pressures(5000.0),&[],controls(),&[0.0;4],
            |_,_,flow|Ok::<_,()>(exchange(flow,345.0))).unwrap();
        let temperatures:Vec<_>=[CellId(20),CellId(21)].into_iter().map(|id|actual.thermal.temperature_k(id).unwrap()).collect();
        assert!((temperatures[0]-temperatures[1]).abs()>10.0);
        let original=owner.momentum().at_temperatures(&temperatures,cx).unwrap();
        let mut reference_state=state.clone(); let mut reference_work=owner.momentum().workspace().unwrap();
        let direct=original.step(&mut reference_state,&mut reference_work,&pressures(5000.0),&[],0.1,controls().hydraulics,cx).unwrap();
        assert_eq!(snapshot(&reference_state),snapshot(&actual.momentum));
        assert_eq!(direct.branches,actual.report.hydraulics.branches);
        let mut frozen_state=state.clone();
        let frozen=owner.momentum().step(&mut frozen_state,&mut reference_work,&pressures(5000.0),&[],0.1,controls().hydraulics,cx).unwrap();
        assert_ne!(frozen.branches,actual.report.hydraulics.branches);
        for (i,&t) in temperatures.iter().enumerate() {
            assert_eq!(actual.next.momentum().pipes()[i].pipe().temperature_k(),t);
            assert_eq!(actual.next.momentum().pipes()[i].pipe().law(),owner.momentum().pipes()[i].pipe().law());
        }
        assert!(actual.report.total_branch_solves>actual.report.hydraulics.branch_solves);
        assert!(actual.report.total_branch_solves<=controls().maximum_total_branch_solves);
        assert!(actual.report.combined_energy_residual_j.abs()<=actual.report.combined_energy_tolerance_j);
    });
}

#[test]
fn g1_zero_and_balanced_steady_limits_retain_original_finite_storage() {
    with_cx(|cx| {
        let caloric=Caloric; let initial=inventory(&caloric,&[300.0;3]);
        let owner=model(cx,&initial,1,false); let mut work=owner.momentum().workspace().unwrap();
        let rest=owner.momentum().rest_state(&[None;2],cx).unwrap();
        let heated=owner.advance(cx,&initial,&rest,&mut work,0.5,&pressures(0.0),&[],controls(),&[0.0;3],
            |_,_,flow| {let mut response=exchange(flow,300.0);response.heat.bulk.sources_w.insert(CellId(20),20.0);Ok::<_,()>(response)}).unwrap();
        assert_eq!(heated.report.hydraulics.viscous_heat_j,0.0);
        assert_eq!(heated.report.hydraulics.kinetic_energy_change_j,0.0);
        assert_eq!(heated.report.external_stream_enthalpy_j,0.0);
        near(heated.report.other_heat_j,10.0,1e-12,0.0);
        near(heated.report.thermal.thermal.stored_energy_change_j,10.0,1e-7,0.0);
        // Exact discrete steady vector K u = pressure * annular area.
        let area=core::f64::consts::PI*R*R;let k=2.0*core::f64::consts::PI*mu(300.0)*L;
        let outer=5000.0*area/(4.0*k);let inner=outer+5000.0*area/(4.0*k);
        let steady=owner.momentum().state_from_velocities(&[vec![inner,outer]],&[None;2],0.0,cx).unwrap();
        let next=owner.advance(cx,&initial,&steady,&mut work,0.2,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow| {let mut response=exchange(flow,300.0);
                response.heat.bulk.sources_w.insert(CellId(20),-flow.viscous_heat_j/flow.duration_s);Ok::<_,()>(response)}).unwrap();
        for (&a,&b) in next.momentum.branches()[0].velocities_m_s().iter().zip(steady.branches()[0].velocities_m_s()) {near(a,b,1e-12,1e-12);}
        for id in [CellId(0),CellId(1),CellId(20)] {near(next.thermal.temperature_k(id).unwrap(),300.0,1e-9,0.0);}
        near(next.report.other_heat_j,-next.report.hydraulics.viscous_heat_j,1e-12,1e-12);
        near(next.report.hydraulics.kinetic_energy_change_j,0.0,1e-12,0.0);
    });
}

#[test]
fn g3_equal_initial_volume_flows_retain_distinct_radial_heat_and_thermal_response() {
    with_cx(|cx| {
        let caloric=Caloric;let initial=inventory(&caloric,&[300.0;3]);let owner=model(cx,&initial,1,false);
        let mut outputs=Vec::new();
        for velocity in [[0.8,0.2],[0.2,0.4]] {
            let state=owner.momentum().state_from_velocities(&[velocity.to_vec()],&[None;2],0.0,cx).unwrap();
            let mut work=owner.momentum().workspace().unwrap();
            let expected=reference(velocity,&initial,0.01,0.0,300.0);
            let result=owner.advance(cx,&initial,&state,&mut work,0.01,&pressures(0.0),&[],controls(),&[0.0;3],
                |_,_,flow|Ok::<_,()>(exchange(flow,300.0))).unwrap();check(&result,&expected);outputs.push(result);
        }
        near(outputs[0].report.hydraulics.branches[0].initial_volume_flow_m3_s,
            outputs[1].report.hydraulics.branches[0].initial_volume_flow_m3_s,1e-18,0.0);
        assert!((outputs[0].report.hydraulics.viscous_heat_j-outputs[1].report.hydraulics.viscous_heat_j).abs()>1e-5);
        assert!((outputs[0].thermal.temperature_k(CellId(20)).unwrap()-outputs[1].thermal.temperature_k(CellId(20)).unwrap()).abs()>1e-6);
    });
}

#[test]
fn g0_compliance_clock_work_and_stream_refusals_preserve_both_states_and_allow_retry() {
    with_cx(|cx| {
        let caloric=Caloric;let initial=inventory(&caloric,&[300.0;3]);let (regional,momentum)=parts(cx,&initial,1,false);
        let compliant=TransientCarreauGraph::new(regional.graph().clone(),RHO,&[2],vec![None,Some(LinearHydraulicStorage {
            compliance_m3_pa:1e-10,minimum_pressure:Pressure::new(-1e6),maximum_pressure:Pressure::new(1e6),
        })]).unwrap();
        assert!(matches!(TransientRegionalCarreauTransport::new(cx,regional.clone(),compliant),Err(ImplicitTransportError::Invalid(_))));
        let wrong_density=TransientCarreauGraph::new(regional.graph().clone(),900.0,&[2],vec![None;2]).unwrap();
        assert!(TransientRegionalCarreauTransport::new(cx,regional.clone(),wrong_density).is_err());
        let owner=TransientRegionalCarreauTransport::new(cx,regional,momentum).unwrap();
        let state=owner.momentum().rest_state(&[None;2],cx).unwrap();let before=snapshot(&state);let cells=initial.cells().clone();
        let mut work=owner.momentum().workspace().unwrap();
        let wrong_time=owner.momentum().state_from_velocities(&[vec![0.0;2]],&[None;2],0.25,cx).unwrap();
        assert!(owner.advance(cx,&initial,&wrong_time,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,340.0))).is_err());
        let short=TransientRegionalCarreauControls {maximum_total_branch_solves:1,..controls()};
        assert!(matches!(owner.advance(cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],short,&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,340.0))),Err(ImplicitTransportError::Resolve(TransientCarreauTransportError::Momentum(TransientCarreauError::BranchSolveBudget)))));
        assert!(owner.advance(cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|{let mut response=exchange(flow,340.0);response.mass_boundaries.clear();Ok::<_,()>(response)}).is_err());
        assert!(matches!(owner.advance(cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,_|Err::<CarreauTransportExchange<()>,_>("original source refused")),
            Err(ImplicitTransportError::Resolve(TransientCarreauTransportError::Exchange("original source refused")))));
        assert_eq!(snapshot(&state),before);assert_eq!(initial.cells(),&cells);assert_eq!(initial.time_s(),0.0);
        let retry=owner.advance(cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,340.0))).unwrap();
        let mut fresh=owner.momentum().workspace().unwrap();
        let reference=owner.advance(cx,&initial,&state,&mut fresh,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|Ok::<_,()>(exchange(flow,340.0))).unwrap();
        assert_eq!(snapshot(&retry.momentum),snapshot(&reference.momentum));assert_eq!(retry.thermal.cells(),reference.thermal.cells());
    });
}

#[test]
fn g4_cancellation_after_native_momentum_keeps_radial_and_thermal_histories() {
    use fs_alloc::{ArenaConfig,ArenaPool};use fs_exec::{Budget,CancelGate,ExecMode,StreamKey};
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:17,kernel_id:711,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let caloric=Caloric;let initial=inventory(&caloric,&[300.0;3]);let owner=model(&cx,&initial,1,true);
        let state=owner.momentum().rest_state(&[None;2],&cx).unwrap();let before=snapshot(&state);let cells=initial.cells().clone();
        let mut work=owner.momentum().workspace().unwrap();
        assert!(owner.advance(&cx,&initial,&state,&mut work,0.04,&pressures(5000.0),&[],controls(),&[0.0;3],
            |_,_,flow|{gate.request();Ok::<_,()>(exchange(flow,340.0))}).is_err());
        assert_eq!(snapshot(&state),before);assert_eq!(initial.cells(),&cells);assert_eq!(initial.time_s(),0.0);
    });
}

//! Real P1 FEM coupled to a hydraulic bypass/mixing graph. At high NTU the
//! stationary tangent and transpose iterations stall within the same budget.
//! Closed-form slab/NTU balances are independent of either interface solver.
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_airflow::conjugate::{AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::graph::thermal::coupled_transport::solve_coupled_transport_iqn;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::{
    CoupledLinearization, CoupledSensitivityError, InterfaceSolveConfig,
};
use fs_airflow::graph::thermal::transport::{
    BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    Nonlinearity, ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::adjoint::robin::RobinLinearization;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_couple::iqn_ils::IqnIlsConfig;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate,arena,
        StreamKey {seed:29,kernel_id:732,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn budget() -> InterfaceSolveConfig {
    InterfaceSolveConfig {max_iterations:12,absolute_tolerance:1e-9,relative_tolerance:1e-10,relaxation:1.0}
}
const H: [f64;2]=[250.0,35.0];
fn fixture<R>(cx: &Cx<'_>, f: impl FnOnce(&TransportNetwork<'_>,&RobinLinearization,&ConjugateConfig)->R)->R {
    let edge=|name,from,to,q:f64| GraphBranch {from,to,
        loss:LossElement::new(name,LossResistance::new(1.0/(q*q)),0.0,
            SourceProvenance::new("analytic fixture","iqn-adjoint-v1"),ToleranceBasis::Analytic).unwrap()};
    let flow=LossGraph::new(4,vec![edge("first",0,2,0.005),edge("bypass",1,2,0.002),edge("last",2,3,0.007)])
        .unwrap().solve(&[
            FixedPressure {node:0,pressure:Pressure::new(2.0)},
            FixedPressure {node:1,pressure:Pressure::new(2.0)},
            FixedPressure {node:3,pressure:Pressure::new(0.0)},
        ],GraphSolveConfig {max_sweeps:4096,max_node_iterations:80,
            absolute_flow_tolerance:VolumetricFlowRate::new(1e-13),relative_flow_tolerance:1e-12},cx).unwrap();
    let network=TransportNetwork::new(cx,&flow,
        TransportAir {density:Density::new(1.0),specific_heat_j_kg_k:1000.0},
        vec![BranchThermalModel::Exchange(vec![AirSegment::new("left",1.0,H[0]).unwrap()]),
            BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(vec![AirSegment::new("right",1.0,H[1]).unwrap()])],
        &[TransportInlet {node:0,temperature:Temperature::new(330.0)},
            TransportInlet {node:1,temperature:Temperature::new(290.0)}],
        TransportConfig {absolute_flow_tolerance:VolumetricFlowRate::new(1e-12),relative_flow_tolerance:1e-10,
            absolute_heat_tolerance_w:1e-7,relative_heat_tolerance:1e-7}).unwrap();
    let (complex,positions)=box_grid([2,1,1],[0.1,1.0,1.0]);
    let mesh=ConductionMesh::new(complex,positions).unwrap();
    let material=ConductivityModel::isotropic_declared(10.0).unwrap();
    let source=ScalarField::Uniform(0.0);
    let gate=ConjugateConfig {temperature_tolerance_k:1e-9,max_iterations:12,..ConjugateConfig::default()};
    let mut retained=None;
    solve_coupled_transport_iqn(cx,&network,&gate,IqnIlsConfig::default(),|cx,refs| {
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("left",|face| on_box_face(face.centroid[0],0.0),ThermalBc::robin(H[0],refs[0]).unwrap()).unwrap()
            .region("right",|face| on_box_face(face.centroid[0],0.1),ThermalBc::robin(H[1],refs[1]).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let mut config=SolveConfig::default();
        config.nonlinearity=Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:8};
        config.initial=InitialGuess::Uniform(310.0);
        config.linear.tolerance=1e-12; config.stop.residual_rtol=1e-12; config.stop.step_atol=0.0;
        let linear=RobinLinearization::new(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
            element_materials:None,source:&source},config,&["left","right"]).unwrap();
        let states=linear.ports().iter().map(|port| SolidRegionState::from_robin_flux(
            linear.primal().report.robin_fluxes.iter().find(|flux| flux.region==port.name).unwrap())).collect();
        retained=Some(linear); Ok(states)
    }).unwrap();
    f(&network,&retained.unwrap(),&gate)
}
fn exact(h:[f64;2],inlet:[f64;2])->[f64;2] {
    let r0=1.0/(5.0*(-(-h[0]/5.0).exp_m1()));
    let r1=1.0/(7.0*(-(-h[1]/7.0).exp_m1()));
    let mixed=(5.0*inlet[0]+2.0*inlet[1])/7.0;
    let q=(inlet[0]-mixed)/(r0+0.01+r1-1.0/7.0);
    [inlet[0]-q*r0,mixed+q*(r1-1.0/7.0)]
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<3e-6*b.abs().max(1.0),"{a:.15e} != {b:.15e}");}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}

#[test]
fn high_ntu_primal_and_gradients_close_under_the_unchanged_sweep_budget() {
    with_gate(&CancelGate::new_clock_free(),|cx| fixture(cx,|network,solid,gate| {
        let c=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let walls=solid.wall_means(cx,&solid.primal().temperature).unwrap();
        for (actual,expected) in walls.iter().zip(exact(H,[330.0,290.0])) {close(*actual,expected);}
        let mut d=c.zero_direction(); d.inlets_k[0]=1.0; d.inlets_k[1]=-1.0;
        let mut o=c.zero_objective(); o.wall_temperatures[0]=1.0;
        assert!(matches!(c.apply(cx,&d,budget()),Err(CoupledSensitivityError::DidNotConverge {..})));
        assert!(matches!(c.pullback(cx,&o,budget()),Err(CoupledSensitivityError::DidNotConverge {..})));
        let tangent=c.apply_iqn(cx,&d,budget(),IqnIlsConfig::default()).unwrap();
        let gradient=c.pullback_iqn(cx,&o,budget(),IqnIlsConfig::default()).unwrap();
        assert!(tangent.iterations<=12 && gradient.iterations<=12);
        let base=exact(H,[330.0,290.0]); let shifted=exact(H,[331.0,289.0]);
        for i in 0..2 {close(tangent.solid.mean_wall_temperatures_k[i],shifted[i]-base[i]);}
        for i in 0..2 {
            let mut supplies=[330.0,290.0]; supplies[i]+=1.0;
            close(gradient.inlets[i],exact(H,supplies)[0]-base[0]);
            let delta=1e-4_f64; let mut plus=H; let mut minus=H;
            plus[i]*=delta.exp(); minus[i]*=(-delta).exp();
            close(gradient.log_htc[i],(exact(plus,[330.0,290.0])[0]-exact(minus,[330.0,290.0])[0])/(2.0*delta));
        }
        close(gradient.inlets[0]+gradient.inlets[1],1.0);
    }));
}

#[test]
fn accelerated_transpose_identity_keeps_load_and_heat_objective_terms() {
    with_gate(&CancelGate::new_clock_free(),|cx| fixture(cx,|network,solid,gate| {
        let c=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let mut d=c.zero_direction(); d.inlets_k[0]=0.4; d.inlets_k[1]=-0.1; d.log_htc=vec![0.3,-0.2];
        for (i,load) in d.nodal_load_w.iter_mut().enumerate() {*load=(i%3) as f64*0.03;}
        let mut o=c.zero_objective(); o.nodal_temperatures.fill(0.02); o.wall_temperatures=vec![0.7,-0.3];
        o.solid_heat_rates=vec![0.03,-0.02]; o.air.branch_outlets[0]=0.8; o.air.branch_outlets[2]=-0.2;
        o.air.references=vec![-0.1,0.4]; o.air.wall_heat_rate=0.02; o.air.external_heat_gain=-0.01;
        let jvp=c.apply_iqn(cx,&d,budget(),IqnIlsConfig::default()).unwrap();
        let vjp=c.pullback_iqn(cx,&o,budget(),IqnIlsConfig::default()).unwrap();
        let mut lhs=dot(&o.nodal_temperatures,&jvp.solid.temperature_k)
            +dot(&o.wall_temperatures,&jvp.solid.mean_wall_temperatures_k)
            +dot(&o.solid_heat_rates,&jvp.solid.heat_rates_w)+dot(&o.air.references,&jvp.air.reference_temperatures_k)
            +o.air.wall_heat_rate*jvp.air.wall_heat_rate_w+o.air.external_heat_gain*jvp.air.external_heat_gain_w;
        for (weight,value) in o.air.branch_outlets.iter().zip(&jvp.air.branch_outlets_k) {lhs+=weight*value.unwrap_or(0.0);}
        close(lhs,dot(&d.inlets_k,&vjp.inlets)+dot(&d.log_htc,&vjp.log_htc)+dot(&d.nodal_load_w,&vjp.nodal_load));
    }));
}

#[test]
fn zero_invalid_budget_and_cancellation_keep_the_original_derivative_contract() {
    with_gate(&CancelGate::new_clock_free(),|cx| fixture(cx,|network,solid,gate| {
        let c=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let zero=c.pullback_iqn(cx,&c.zero_objective(),budget(),IqnIlsConfig::default()).unwrap();
        assert_eq!(zero.iterations,1); assert!(zero.nodal_load.iter().all(|v| *v==0.0));
        let invalid=IqnIlsConfig {max_history:0,..IqnIlsConfig::default()};
        assert!(matches!(c.apply_iqn(cx,&c.zero_direction(),budget(),invalid),Err(CoupledSensitivityError::Acceleration(_))));
        let mut d=c.zero_direction(); d.inlets_k[0]=1.0;
        let short=InterfaceSolveConfig {max_iterations:1,..budget()};
        assert!(matches!(c.apply_iqn(cx,&d,short,IqnIlsConfig::default()),Err(CoupledSensitivityError::DidNotConverge {..})));
        let blocked=CancelGate::new_clock_free(); blocked.request();
        with_gate(&blocked,|cx| assert!(matches!(
            c.pullback_iqn(cx,&c.zero_objective(),budget(),IqnIlsConfig::default()),Err(CoupledSensitivityError::Interrupted))));
    }));
}

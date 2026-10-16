use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    Nonlinearity, ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};
use crate::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use crate::conjugate::AirSegment;
use crate::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use crate::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet};
use super::super::solve_coupled_transport;

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 39, kernel_id: 815, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R { with_gate(&CancelGate::new(), f) }
fn budget() -> InterfaceSolveConfig { InterfaceSolveConfig { max_iterations: 300,
    absolute_tolerance: 1e-10, relative_tolerance: 1e-10, relaxation: 1.0 } }
fn edge(name: &str, from: usize, to: usize, q: f64) -> GraphBranch {
    GraphBranch { from, to, loss: LossElement::new(name, LossResistance::new(1.0/(q*q)),
        0.0, SourceProvenance::new("analytic fixture", "coupled-gradient-v1"), ToleranceBasis::Analytic).unwrap() }
}
fn fixture<R>(cx: &Cx<'_>, h: [f64;2], supplies: [f64;2], mixed: bool, wrong_point: bool,
    f: impl FnOnce(&TransportNetwork<'_>, &RobinLinearization, &ConjugateConfig) -> R) -> R {
    let edges = if mixed { vec![edge("first",0,2,0.005),edge("bypass",1,2,0.002),edge("last",2,3,0.007)] }
        else { vec![edge("first",0,2,0.005),edge("last",1,3,0.007)] };
    let graph = LossGraph::new(4,edges).unwrap();
    let mut fixed = vec![FixedPressure { node:0,pressure:Pressure::new(if mixed {2.0} else {1.0}) },
        FixedPressure { node:1,pressure:Pressure::new(if mixed {2.0} else {1.0}) },
        FixedPressure { node:3,pressure:Pressure::new(0.0) }];
    if !mixed { fixed.push(FixedPressure { node:2,pressure:Pressure::new(0.0) }); }
    let flow = graph.solve(&fixed,GraphSolveConfig { max_sweeps:4096,max_node_iterations:80,
        absolute_flow_tolerance:VolumetricFlowRate::new(1e-13),relative_flow_tolerance:1e-12 },cx).unwrap();
    let mut models = vec![BranchThermalModel::Exchange(vec![AirSegment::new("left",1.0,h[0]).unwrap()])];
    if mixed { models.push(BranchThermalModel::Adiabatic); }
    models.push(BranchThermalModel::Exchange(vec![AirSegment::new("right",1.0,h[1]).unwrap()]));
    let network = TransportNetwork::new(cx,&flow,TransportAir { density:Density::new(1.0),specific_heat_j_kg_k:1000.0 },models,
        &[TransportInlet {node:0,temperature:Temperature::new(supplies[0])},TransportInlet {node:1,temperature:Temperature::new(supplies[1])}],
        TransportConfig {absolute_flow_tolerance:VolumetricFlowRate::new(1e-12),relative_flow_tolerance:1e-10,
            absolute_heat_tolerance_w:1e-7,relative_heat_tolerance:1e-7}).unwrap();
    let (complex, positions) = box_grid([2,1,1],[0.1,1.0,1.0]);
    let mesh = ConductionMesh::new(complex,positions).unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let build = |refs: &[f64]| {
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("left",|face| on_box_face(face.centroid[0],0.0),ThermalBc::robin(h[0],refs[0]).unwrap()).unwrap()
            .region("right",|face| on_box_face(face.centroid[0],0.1),ThermalBc::robin(h[1],refs[1]).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let mut config = SolveConfig::default();
        config.nonlinearity = Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:8};
        config.initial = InitialGuess::Uniform(310.0);
        config.linear.tolerance=1e-10; config.stop.residual_rtol=1e-11; config.stop.step_atol=0.0;
        RobinLinearization::new(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
            element_materials:None,source:&source},config,&["left","right"]).unwrap()
    };
    let gate = ConjugateConfig { temperature_tolerance_k:1e-8, max_iterations:300, ..ConjugateConfig::default() };
    let mut retained = None;
    solve_coupled_transport(cx,&network,&gate,|_,refs| {
        let linear=build(refs);
        let states=linear.ports().iter().map(|port| SolidRegionState::from_robin_flux(
            linear.primal().report.robin_fluxes.iter().find(|flux| flux.region==port.name).unwrap())).collect();
        retained=Some(linear);
        Ok(states)
    }).unwrap();
    let solid=if wrong_point { build(&[310.0,310.0]) } else { retained.unwrap() };
    f(&network,&solid,&gate)
}
fn exact(h:[f64;2],inlet:[f64;2],mixed:bool)->[f64;2] {
    let r0=1.0/(5.0*(-(-h[0]/5.0).exp_m1()));
    let r1=1.0/(7.0*(-(-h[1]/7.0).exp_m1()));
    let downstream=if mixed {(5.0*inlet[0]+2.0*inlet[1])/7.0} else {inlet[1]};
    let q=(inlet[0]-downstream)/(r0+0.01+r1-if mixed {1.0/7.0} else {0.0});
    [inlet[0]-q*r0,downstream+q*(r1-if mixed {1.0/7.0} else {0.0})]
}
fn close(a:f64,b:f64) { assert!((a-b).abs()<3e-6*b.abs().max(1.0),"{a:.14e} != {b:.14e}"); }
fn dot(a:&[f64],b:&[f64])->f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }

#[test]
fn fully_coupled_wall_gradients_match_independent_closed_forms_with_and_without_mixing() {
    with_cx(|cx| for mixed in [false,true] {
        fixture(cx,[2.0,3.0],[330.0,290.0],mixed,false,|network,solid,gate| {
            let coupled=CoupledLinearization::new(cx,network,solid,gate).unwrap();
            for wall in 0..2 {
                let mut objective=coupled.zero_objective(); objective.wall_temperatures[wall]=1.0;
                let gradient=coupled.pullback(cx,&objective,budget()).unwrap();
                let eps=1e-4_f64;
                for control in 0..2 {
                    let mut plus=[2.0,3.0]; let mut minus=plus;
                    plus[control]*=eps.exp(); minus[control]*=(-eps).exp();
                    close(gradient.log_htc[control],(exact(plus,[330.0,290.0],mixed)[wall]-exact(minus,[330.0,290.0],mixed)[wall])/(2.0*eps));
                    let mut plus=[330.0,290.0]; let mut minus=plus;
                    plus[control]+=eps; minus[control]-=eps;
                    close(gradient.inlets[control],(exact([2.0,3.0],plus,mixed)[wall]-exact([2.0,3.0],minus,mixed)[wall])/(2.0*eps));
                }
                close(gradient.inlets[0]+gradient.inlets[1],1.0);
            }
        });
    });
}

#[test]
fn coupled_tangent_matches_perturbed_shared_solid_fem_not_a_thermostat() {
    with_cx(|cx| {
        let (actual, frozen)=fixture(cx,[2.0,3.0],[330.0,290.0],true,false,|network,solid,gate| {
            let coupled=CoupledLinearization::new(cx,network,solid,gate).unwrap();
            let mut direction=coupled.zero_direction(); direction.log_htc=vec![0.3,-0.2];
            direction.inlets_k[0]=0.4; direction.inlets_k[1]=-0.1;
            let result=coupled.apply(cx,&direction,budget()).unwrap();
            let mut objective=coupled.zero_objective(); objective.air.branch_outlets[0]=1.0;
            let total=coupled.pullback(cx,&objective,budget()).unwrap();
            let frozen=coupled.air.pullback(cx,&objective.air).unwrap();
            assert!((total.log_htc[0]-frozen.log_conductances[0]).abs()>0.01);
            (result.solid.temperature_k,frozen.log_conductances[0])
        });
        assert!(frozen.is_finite());
        let eps=1e-4_f64;
        let plus=fixture(cx,[2.0*(0.3*eps).exp(),3.0*(-0.2*eps).exp()],[330.0+0.4*eps,290.0-0.1*eps],true,false,
            |_,solid,_| solid.primal().temperature.clone());
        let minus=fixture(cx,[2.0*(-0.3*eps).exp(),3.0*(0.2*eps).exp()],[330.0-0.4*eps,290.0+0.1*eps],true,false,
            |_,solid,_| solid.primal().temperature.clone());
        for ((a,b),expected) in plus.iter().zip(&minus).zip(actual) { close((a-b)/(2.0*eps),expected); }
    });
}

#[test]
fn full_coupled_transpose_identity_includes_source_loads_and_heat_objectives() {
    with_cx(|cx| fixture(cx,[2.0,3.0],[330.0,290.0],true,false,|network,solid,gate| {
        let coupled=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let mut d=coupled.zero_direction(); d.log_htc=vec![0.3,-0.2]; d.inlets_k[0]=0.4; d.inlets_k[1]=-0.1;
        for (i,load) in d.nodal_load_w.iter_mut().enumerate() { *load=(i%3) as f64*0.03; }
        let mut o=coupled.zero_objective(); o.nodal_temperatures.fill(0.02); o.wall_temperatures=vec![0.7,-0.3];
        o.solid_heat_rates=vec![0.03,-0.02]; o.air.branch_outlets[0]=0.8; o.air.branch_outlets[2]=-0.2;
        o.air.references=vec![-0.1,0.4]; o.air.wall_heat_rate=0.02; o.air.external_heat_gain=-0.01;
        let jvp=coupled.apply(cx,&d,budget()).unwrap(); let vjp=coupled.pullback(cx,&o,budget()).unwrap();
        let mut lhs=dot(&o.nodal_temperatures,&jvp.solid.temperature_k)+dot(&o.wall_temperatures,&jvp.solid.mean_wall_temperatures_k)
            +dot(&o.solid_heat_rates,&jvp.solid.heat_rates_w)+dot(&o.air.references,&jvp.air.reference_temperatures_k)
            +o.air.wall_heat_rate*jvp.air.wall_heat_rate_w+o.air.external_heat_gain*jvp.air.external_heat_gain_w;
        for (weight,value) in o.air.branch_outlets.iter().zip(&jvp.air.branch_outlets_k) { lhs+=weight*value.unwrap_or(0.0); }
        close(lhs,dot(&d.inlets_k,&vjp.inlets)+dot(&d.log_htc,&vjp.log_htc)+dot(&d.nodal_load_w,&vjp.nodal_load));
    }));
}

#[test]
fn shared_solid_cannot_change_total_outlet_enthalpy_without_a_heat_source() {
    with_cx(|cx| fixture(cx,[2.0,3.0],[330.0,290.0],true,false,|network,solid,gate| {
        let coupled=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let mut o=coupled.zero_objective(); o.air.branch_outlets[2]=1.0;
        let g=coupled.pullback(cx,&o,budget()).unwrap();
        close(g.log_htc[0],0.0); close(g.log_htc[1],0.0);
        close(g.inlets[0],5.0/7.0); close(g.inlets[1],2.0/7.0);
    }));
}

#[test]
fn nonconverged_primal_invalid_budget_and_budget_exhaustion_refuse() {
    with_cx(|cx| {
        fixture(cx,[2.0,3.0],[330.0,290.0],true,true,|n,s,g| assert!(CoupledLinearization::new(cx,n,s,g).is_err()));
        fixture(cx,[2.0,3.0],[330.0,290.0],true,false,|n,s,g| {
            let c=CoupledLinearization::new(cx,n,s,g).unwrap();
            let mut d=c.zero_direction(); d.inlets_k[0]=1.0;
            let mut b=budget(); b.max_iterations=1;
            assert!(matches!(c.apply(cx,&d,b),Err(CoupledSensitivityError::DidNotConverge{..})));
            b=budget(); b.relaxation=0.0; assert!(c.apply(cx,&d,b).is_err());
            d.log_htc.pop(); assert!(c.apply(cx,&d,budget()).is_err());
        });
    });
}

#[test]
fn cancellation_and_zero_objective_do_not_publish_partial_gradients() {
    with_cx(|cx| fixture(cx,[2.0,3.0],[330.0,290.0],true,false,|n,s,g| {
        let c=CoupledLinearization::new(cx,n,s,g).unwrap();
        let zero=c.pullback(cx,&c.zero_objective(),budget()).unwrap();
        assert_eq!(zero.iterations,1); assert!(zero.nodal_load.iter().all(|v| *v==0.0));
        let gate=CancelGate::new(); gate.request();
        with_gate(&gate,|blocked| assert!(matches!(c.pullback(blocked,&c.zero_objective(),budget()),Err(CoupledSensitivityError::Interrupted))));
    }));
}

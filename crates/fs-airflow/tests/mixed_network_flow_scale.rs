use fs_airflow::conjugate::{AirSegment, ConjugateConfig, IqnIlsConfig, SolidRegionState};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::graph::thermal::coupled_transport::solve_coupled_transport_iqn;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::{
    CoupledLinearization, CoupledSensitivityError, InterfaceSolveConfig,
};
use fs_airflow::graph::thermal::transport::{
    BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork,
};
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::adjoint::robin::RobinLinearization;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 43, kernel_id: 817, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn budget() -> InterfaceSolveConfig {
    InterfaceSolveConfig { max_iterations: 40, absolute_tolerance: 1e-11,
        relative_tolerance: 1e-11, relaxation: 1.0 }
}
fn edge(name: &str, from: usize, to: usize, q: f64) -> GraphBranch {
    GraphBranch { from, to, loss: LossElement::new(name, LossResistance::new(1.0/(q*q)),
        0.0, SourceProvenance::new("analytic fixture", "flow-scale-v1"), ToleranceBasis::Analytic).unwrap() }
}
fn fixture<R>(cx: &Cx<'_>, speed: f64, source: f64,
    f: impl FnOnce(&TransportNetwork<'_>, &RobinLinearization, &ConjugateConfig) -> R) -> R {
    let graph = LossGraph::new(4, vec![edge("first",0,2,0.005),
        edge("bypass",1,2,0.002), edge("last",2,3,0.007)]).unwrap();
    let flow = graph.solve(&[
        FixedPressure { node: 0, pressure: Pressure::new(2.0*speed*speed) },
        FixedPressure { node: 1, pressure: Pressure::new(2.0*speed*speed) },
        FixedPressure { node: 3, pressure: Pressure::new(0.0) },
    ], GraphSolveConfig { max_sweeps:4096, max_node_iterations:80,
        absolute_flow_tolerance:VolumetricFlowRate::new(1e-14),relative_flow_tolerance:1e-13 },cx).unwrap();
    let network = TransportNetwork::new(cx, &flow,
        TransportAir { density:Density::new(1.0),specific_heat_j_kg_k:1000.0 },
        vec![BranchThermalModel::Exchange(vec![AirSegment::new("left",1.0,2.0).unwrap()]),
            BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(vec![AirSegment::new("right",1.0,3.0).unwrap()])],
        &[TransportInlet {node:0,temperature:Temperature::new(330.0)},
            TransportInlet {node:1,temperature:Temperature::new(290.0)}],
        TransportConfig {absolute_flow_tolerance:VolumetricFlowRate::new(1e-12),relative_flow_tolerance:1e-10,
            absolute_heat_tolerance_w:1e-7,relative_heat_tolerance:1e-7}).unwrap();
    let (complex,positions) = box_grid([2,1,1],[0.1,1.0,1.0]);
    let mesh = ConductionMesh::new(complex,positions).unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(source);
    let gate = ConjugateConfig { max_iterations:40,temperature_tolerance_k:1e-9,
        ..ConjugateConfig::default() };
    let mut retained = None;
    solve_coupled_transport_iqn(cx,&network,&gate,IqnIlsConfig::default(),|_,refs| {
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("left",|face|on_box_face(face.centroid[0],0.0),ThermalBc::robin(2.0,refs[0]).unwrap()).unwrap()
            .region("right",|face|on_box_face(face.centroid[0],0.1),ThermalBc::robin(3.0,refs[1]).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let mut config=SolveConfig::default();
        config.initial=InitialGuess::Uniform(310.0);
        config.linear.tolerance=1e-12; config.stop.residual_rtol=1e-12; config.stop.step_atol=0.0;
        let linear=RobinLinearization::new(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,
            material:&material,element_materials:None,source:&source},config,&["left","right"]).unwrap();
        let response=linear.ports().iter().map(|port|SolidRegionState::from_robin_flux(
            linear.primal().report.robin_fluxes.iter().find(|flux|flux.region==port.name).unwrap())).collect();
        retained=Some(linear);
        Ok(response)
    }).unwrap();
    f(&network,&retained.unwrap(),&gate)
}
fn exact_left(speed:f64)->f64 {
    let c0=5.0*speed; let c1=7.0*speed;
    let r0=1.0/(c0*(-(-2.0/c0).exp_m1()));
    let r1=1.0/(c1*(-(-3.0/c1).exp_m1()));
    let mixed=(5.0*330.0+2.0*290.0)/7.0;
    let q=(330.0-mixed)/(r0+0.01+r1-1.0/c1);
    330.0-q*r0
}

#[test]
fn flow_scale_gradient_matches_the_independent_mixed_slab_solution() {
    with_cx(&CancelGate::new_clock_free(),|cx| {
        for speed in [0.5,1.0,1.6] {
            fixture(cx,speed,0.0,|network,solid,gate| {
                let binding=CoupledLinearization::new(cx,network,solid,gate).unwrap();
                let mut objective=binding.zero_objective(); objective.wall_temperatures[0]=1.0;
                let result=binding.pullback_flow_scale_iqn(cx,&objective,budget(),IqnIlsConfig::default()).unwrap();
                let eps=1e-5_f64;
                let expected=(exact_left(speed*eps.exp())-exact_left(speed*(-eps).exp()))/(2.0*eps);
                assert!((result.log_flow_scale-expected).abs()<2e-5,
                    "speed={speed}: {} != {expected}",result.log_flow_scale);
                let old=binding.pullback_iqn(cx,&objective,budget(),IqnIlsConfig::default()).unwrap();
                assert_eq!(old.log_htc,result.thermal.log_htc);
                assert_eq!(old.nodal_load,result.thermal.nodal_load);
            });
        }
    });
}

#[test]
fn heat_functionals_keep_the_explicit_capacity_term() {
    with_cx(&CancelGate::new_clock_free(),|cx| fixture(cx,1.0,100.0,|network,solid,gate| {
        let binding=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        // The 0.1 m^3 solid supplies 10 W at every flow. Omitting the
        // degree-one heat term from the chain rule would give roughly -10.
        let mut objective=binding.zero_objective(); objective.air.external_heat_gain=1.0;
        let heat=binding.pullback_flow_scale_iqn(cx,&objective,budget(),IqnIlsConfig::default()).unwrap();
        assert!(heat.log_flow_scale.abs()<1e-7);
        let mut objective=binding.zero_objective(); objective.air.branch_outlets[2]=1.0;
        let outlet=binding.pullback_flow_scale_iqn(cx,&objective,budget(),IqnIlsConfig::default()).unwrap();
        assert!((outlet.log_flow_scale+10.0/7.0).abs()<1e-7);
    }));
}

#[test]
fn cancelled_or_unconverged_adjoint_cannot_publish_flow_gradients() {
    with_cx(&CancelGate::new_clock_free(),|cx| fixture(cx,1.0,0.0,|network,solid,gate| {
        let binding=CoupledLinearization::new(cx,network,solid,gate).unwrap();
        let mut objective=binding.zero_objective(); objective.wall_temperatures[0]=1.0;
        let mut config=budget(); config.max_iterations=1;
        assert!(matches!(binding.pullback_flow_scale_iqn(cx,&objective,config,IqnIlsConfig::default()),
            Err(CoupledSensitivityError::DidNotConverge { .. })));
        let cancel=CancelGate::new_clock_free(); cancel.request();
        with_cx(&cancel,|blocked| assert!(matches!(
            binding.pullback_flow_scale_iqn(blocked,&objective,budget(),IqnIlsConfig::default()),
            Err(CoupledSensitivityError::Interrupted))));
    }));
}

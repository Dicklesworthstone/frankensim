use super::*;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    Nonlinearity, ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_qty::{Density, Temperature, VolumetricFlowRate};
use crate::conjugate::{AirSegment, SolidRegionState};
use crate::graph::{GraphSolveConfig, LossGraph};
use crate::graph::tests::{boundary, edge, with_cx};
use crate::graph::thermal::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet};
use crate::graph::thermal::coupled_transport::solve_coupled_transport;

fn budget() -> InterfaceSolveConfig { InterfaceSolveConfig { max_iterations:300,
    absolute_tolerance:1e-10,relative_tolerance:1e-10,relaxation:1.0 } }
fn fixture<R>(cx: &Cx<'_>, q: [f64;2], f: impl FnOnce(&CoupledLinearization<'_, '_>) -> R) -> R {
    let graph=LossGraph::new(4,vec![edge("first",0,2,1.0/(q[0]*q[0])),
        edge("bypass",1,2,1.0/(q[1]*q[1])),edge("last",2,3,1.0/((q[0]+q[1])*(q[0]+q[1]))) ]).unwrap();
    let flow=graph.solve(&[boundary(0,2.0),boundary(1,2.0),boundary(3,0.0)],GraphSolveConfig {
        max_sweeps:4096,max_node_iterations:80,absolute_flow_tolerance:VolumetricFlowRate::new(1e-13),relative_flow_tolerance:1e-12 },cx).unwrap();
    let models=vec![BranchThermalModel::Exchange(vec![AirSegment::new("left",1.0,2.0).unwrap()]),
        BranchThermalModel::Adiabatic,BranchThermalModel::Exchange(vec![AirSegment::new("right",1.0,3.0).unwrap()])];
    let network=TransportNetwork::new(cx,&flow,TransportAir {density:Density::new(1.0),specific_heat_j_kg_k:1000.0},models,
        &[TransportInlet {node:0,temperature:Temperature::new(330.0)},TransportInlet {node:1,temperature:Temperature::new(290.0)}],
        TransportConfig {absolute_flow_tolerance:VolumetricFlowRate::new(1e-12),relative_flow_tolerance:1e-10,
            absolute_heat_tolerance_w:1e-7,relative_heat_tolerance:1e-7}).unwrap();
    let (complex,positions)=box_grid([2,1,1],[0.1,1.0,1.0]);
    let mesh=ConductionMesh::new(complex,positions).unwrap();
    let material=ConductivityModel::isotropic_declared(10.0).unwrap();
    let source=ScalarField::Uniform(0.0);
    let gate=ConjugateConfig {temperature_tolerance_k:1e-8,max_iterations:300,..ConjugateConfig::default()};
    let mut retained=None;
    solve_coupled_transport(cx,&network,&gate,|_,refs| {
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("left",|face|on_box_face(face.centroid[0],0.0),ThermalBc::robin(2.0,refs[0]).unwrap()).unwrap()
            .region("right",|face|on_box_face(face.centroid[0],0.1),ThermalBc::robin(3.0,refs[1]).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let mut config=SolveConfig::default();
        config.nonlinearity=Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:8};
        config.initial=InitialGuess::Uniform(310.0);config.linear.tolerance=1e-10;
        config.stop.residual_rtol=1e-11;config.stop.step_atol=0.0;
        let linear=RobinLinearization::new(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
            element_materials:None,source:&source},config,&["left","right"]).unwrap();
        let states=linear.ports().iter().map(|port|SolidRegionState::from_robin_flux(
            linear.primal().report.robin_fluxes.iter().find(|flux|flux.region==port.name).unwrap())).collect();
        retained=Some(linear);Ok(states)
    }).unwrap();
    let solid=retained.unwrap();
    let binding=CoupledLinearization::new(cx,&network,&solid,&gate).unwrap();
    f(&binding)
}
fn exact(q:[f64;2])->[f64;2] {
    let c0=1000.0*q[0];let cb=1000.0*q[1];let c=c0+cb;
    let r0=1.0/(c0*(-(-2.0/c0).exp_m1()));let r1=1.0/(c*(-(-3.0/c).exp_m1()));
    let mixed=(c0*330.0+cb*290.0)/c;
    let heat=(330.0-mixed)/(r0+0.01+r1-1.0/c);
    [330.0-heat*r0,mixed+heat*(r1-1.0/c)]
}
fn close(a:f64,b:f64) { assert!((a-b).abs()<3e-5*b.abs().max(1.0),"{a} != {b}"); }

#[test]
fn total_flow_controls_match_independent_shared_solid_and_mixing_solution() {
    with_cx(|cx|fixture(cx,[0.005,0.002],|binding| {
        for wall in 0..2 {
            let mut objective=binding.zero_objective();objective.wall_temperatures[wall]=1.0;
            let (_,flows)=binding.pullback_branch_flows(cx,&objective,budget()).unwrap();
            for source in 0..2 {
                let h=1e-7;let mut plus=[0.005,0.002];let mut minus=plus;
                plus[source]+=h;minus[source]-=h;
                // The downstream branch changes too: this is a conserved direction.
                close(flows[source]+flows[2],(exact(plus)[wall]-exact(minus)[wall])/(2.0*h));
            }
            let common=binding.pullback_flow_scale(cx,&objective,budget()).unwrap();
            close(common.log_flow_scale,flows.iter().zip(&binding.hydraulics().branches)
                .map(|(g,b)|g*b.flow.value()).sum());
        }
    }));
}

#[test]
fn total_flow_control_matches_perturbed_production_fem_and_refuses_zero_budget() {
    with_cx(|cx| {
        let gradient=fixture(cx,[0.005,0.002],|binding| {
            let mut objective=binding.zero_objective();objective.wall_temperatures[0]=1.0;
            let mut invalid=budget();invalid.max_iterations=0;
            assert!(binding.pullback_branch_flows(cx,&objective,invalid).is_err());
            let (_,g)=binding.pullback_branch_flows(cx,&objective,budget()).unwrap();g[0]+g[2]
        });
        let h=1e-6;
        let run=|q|fixture(cx,q,|binding|binding.solid.wall_means(cx,binding.solid.temperature()).unwrap()[0]);
        close(gradient,(run([0.005+h,0.002])-run([0.005-h,0.002]))/(2.0*h));
    });
}

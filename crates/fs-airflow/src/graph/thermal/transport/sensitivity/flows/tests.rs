use super::*;
use crate::conjugate::AirSegment;
use crate::graph::{GraphSolution, LossGraph};
use crate::graph::tests::{boundary, bridge, config as hydraulic_config, edge, with_cx};
use crate::graph::thermal::transport::{TransportAir, TransportConfig, TransportInlet};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_qty::{Density, Temperature, VolumetricFlowRate};

fn exchange(name: &str, u: f64) -> BranchThermalModel {
    BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, u).unwrap()])
}
fn inlet(node: usize, t: f64) -> TransportInlet {
    TransportInlet { node, temperature: Temperature::new(t) }
}
fn network<'a>(cx: &Cx<'_>, flow: &'a GraphSolution, models: &[BranchThermalModel], inlets: &[TransportInlet]) -> TransportNetwork<'a> {
    // Deliberately admit small off-conservation coordinate perturbations:
    // each raw energy/defect derivative is checked, not repaired to zero.
    TransportNetwork::new(cx, flow,
        TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 2.0 },
        models.to_vec(), inlets, TransportConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(0.01), relative_flow_tolerance: 0.0,
            absolute_heat_tolerance_w: 1.0, relative_heat_tolerance: 0.0,
        }).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a:.16e} != {b:.16e}");
}
fn value(m: &TransportMarch, w: &TransportObjective) -> f64 {
    m.node_temperatures_k.iter().zip(&w.node_temperatures).map(|(t,w)| t.unwrap_or(0.0)*w).sum::<f64>()
        + m.branches.iter().zip(&w.branch_outlets).map(|(b,w)| b.outlet_temperature_k.unwrap_or(0.0)*w).sum::<f64>()
        + m.reference_temperatures_k.iter().zip(&w.references).map(|(t,w)| t*w).sum::<f64>()
        + w.wall_heat_rate*m.wall_heat_rate_w + w.external_heat_gain*m.external_heat_gain_w
        + w.heat_imbalance*m.heat_imbalance_w + w.hydraulic_energy_defect*m.hydraulic_energy_defect_w
}
fn check_columns(cx: &Cx<'_>, flow: &GraphSolution, models: &[BranchThermalModel], supplies: &[TransportInlet], walls: &[f64]) {
    let net = network(cx, flow, models, supplies);
    let lin = net.linearize(cx, walls).unwrap();
    let mut objectives = Vec::new();
    for kind in 0..7 {
        let mut w = lin.zero_objective();
        match kind {
            0 => w.node_temperatures.fill(0.4),
            1 => w.branch_outlets.fill(-0.7),
            2 => w.references.fill(0.3),
            3 => w.wall_heat_rate = 1.0,
            4 => w.external_heat_gain = 1.0,
            5 => w.heat_imbalance = 1.0,
            _ => w.hydraulic_energy_defect = 1.0,
        }
        objectives.push(w);
    }
    for objective in objectives {
        let gradient = lin.pullback_branch_flows(cx, &objective).unwrap();
        for (branch, &actual) in gradient.iter().enumerate() {
            let h = 1e-5 * flow.branches[branch].flow.value().abs();
            let run = |step: f64| {
                let mut changed = flow.clone();
                changed.branches[branch].flow = VolumetricFlowRate::new(flow.branches[branch].flow.value() + step);
                value(&network(cx, &changed, models, supplies).march(cx, walls).unwrap(), &objective)
            };
            close(actual, (run(h)-run(-h))/(2.0*h), 2e-6);
        }
    }
}

#[test]
fn single_channel_flow_columns_match_closed_form() {
    with_cx(|cx| {
        let flow = LossGraph::new(2, vec![edge("duct",0,1,1.0)]).unwrap()
            .solve(&[boundary(0,4.0),boundary(1,0.0)],hydraulic_config(),cx).unwrap();
        let net = network(cx, &flow, &[exchange("wall",4.0)], &[inlet(0,300.0)]);
        let lin = net.linearize(cx, &[340.0]).unwrap();
        let a = (-1.0_f64).exp();
        let mut w = lin.zero_objective(); w.branch_outlets[0]=1.0;
        close(lin.pullback_branch_flows(cx,&w).unwrap()[0],-20.0*a,1e-12);
        w.branch_outlets[0]=0.0; w.references[0]=1.0;
        close(lin.pullback_branch_flows(cx,&w).unwrap()[0],-20.0*(1.0-2.0*a),1e-12);
        w.references[0]=0.0; w.wall_heat_rate=1.0;
        close(lin.pullback_branch_flows(cx,&w).unwrap()[0],80.0*(1.0-2.0*a),1e-12);
    });
}

#[test]
fn reverse_branch_mixing_and_every_energy_output_match_perturbed_marches() {
    with_cx(|cx| {
        let flow = bridge().solve(&[boundary(0,25.0),boundary(3,0.0)],hydraulic_config(),cx).unwrap();
        let models = vec![exchange("feed",1.0),BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(vec![AirSegment::new("first",1.0,1.2).unwrap(),
                AirSegment::new("last",1.0,0.7).unwrap()]),
            BranchThermalModel::Adiabatic,exchange("out",2.0)];
        assert!(flow.branches[2].flow.value()<0.0);
        check_columns(cx,&flow,&models,&[inlet(0,290.0)],&[345.0,320.0,335.0,330.0]);
    });
}

#[test]
fn intermediate_reservoir_injection_and_extraction_are_not_frozen() {
    with_cx(|cx| {
        for middle in [4.0,8.0] {
            let flow = LossGraph::new(3,vec![edge("a",0,1,1.0),edge("b",1,2,1.0)]).unwrap()
                .solve(&[boundary(0,9.0),boundary(1,middle),boundary(2,0.0)],hydraulic_config(),cx).unwrap();
            let mut supplies=vec![inlet(0,310.0)];
            if middle==8.0 { supplies.push(inlet(1,285.0)); }
            check_columns(cx,&flow,&[exchange("a",1.0),exchange("b",2.0)],&supplies,&[345.0,330.0]);
        }
    });
}

#[test]
fn common_scale_matches_euler_identity_with_the_raw_residual_retained() {
    with_cx(|cx| {
        let mut flow = bridge().solve(&[boundary(0,25.0),boundary(3,0.0)],hydraulic_config(),cx).unwrap();
        flow.branches[1].flow = VolumetricFlowRate::new(flow.branches[1].flow.value()*1.00001);
        let net = network(cx,&flow,&[exchange("a",1.0),exchange("b",2.0),exchange("c",3.0),
            exchange("d",4.0),exchange("e",5.0)],&[inlet(0,290.0)]);
        let lin=net.linearize(cx,&[345.0,320.0,335.0,325.0,330.0]).unwrap();
        let mut w=lin.zero_objective(); w.references.fill(0.4); w.wall_heat_rate=0.3;
        w.external_heat_gain=-0.2; w.heat_imbalance=0.7; w.hydraulic_energy_defect=-0.5;
        let flows=lin.pullback_branch_flows(cx,&w).unwrap();
        let conditional=lin.pullback(cx,&w).unwrap();
        let degree_one=w.wall_heat_rate*lin.primal.wall_heat_rate_w+w.external_heat_gain*lin.primal.external_heat_gain_w
            +w.heat_imbalance*lin.primal.heat_imbalance_w+w.hydraulic_energy_defect*lin.primal.hydraulic_energy_defect_w;
        close(flows.iter().zip(&flow.branches).map(|(g,b)|g*b.flow.value()).sum::<f64>()
            +conditional.log_conductances.iter().sum::<f64>(),degree_one,1e-10);
    });
}

#[test]
fn zero_flow_zero_reservoir_exchange_and_interruption_refuse() {
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:17,kernel_id:713,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let graph=LossGraph::new(3,vec![edge("a",0,1,1.0),edge("b",1,2,1.0)]).unwrap();
        let flow=graph.solve(&[boundary(0,9.0),boundary(1,4.5),boundary(2,0.0)],hydraulic_config(),&cx).unwrap();
        let net=network(&cx,&flow,&[exchange("a",1.0),exchange("b",2.0)],&[inlet(0,300.0)]);
        let lin=net.linearize(&cx,&[340.0,330.0]).unwrap();
        assert!(matches!(lin.pullback_branch_flows(&cx,&lin.zero_objective()),Err(TransportError::Node {node:1,..})));
        let idle=LossGraph::new(2,vec![edge("idle",0,1,1.0)]).unwrap()
            .solve(&[boundary(0,0.0),boundary(1,0.0)],hydraulic_config(),&cx).unwrap();
        let idle_net=network(&cx,&idle,&[BranchThermalModel::Adiabatic],&[]);
        let idle_lin=idle_net.linearize(&cx,&[]).unwrap();
        assert!(idle_lin.pullback_branch_flows(&cx,&idle_lin.zero_objective()).is_err());
        gate.request();
        assert!(matches!(lin.pullback_branch_flows(&cx,&lin.zero_objective()),Err(TransportError::Interrupted)));
    });
}

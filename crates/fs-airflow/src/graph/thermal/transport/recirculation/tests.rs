use super::*;
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use crate::graph::LossGraph;
use super::super::sensitivity::{TransportDifferential, TransportObjective};

fn config() -> TransportConfig {
    TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(1e-10),
        relative_flow_tolerance: 0.0, absolute_heat_tolerance_w: 1e-7, relative_heat_tolerance: 0.0 }
}
fn air() -> TransportAir { TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 1.0 } }
fn inlet(node: usize, value: f64) -> TransportInlet { TransportInlet { node, temperature: Temperature::new(value) } }
fn heated(name: &str, value: f64) -> BranchThermalModel {
    BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, value).unwrap()])
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a:.16e} vs {b:.16e}"); }
fn link(supply_node: usize, return_node: usize, fraction: f64) -> RecirculationLink {
    RecirculationLink { supply_node, return_node, fraction }
}

#[test]
fn returned_single_exchanger_matches_independent_geometric_feedback_root() {
    with_cx(|cx| {
        let graph = LossGraph::new(2, vec![edge("channel",0,1,1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0,1.0),boundary(1,0.0)],hydraulic_config(),cx).unwrap();
        for fraction in [0.0,0.25,0.8,0.9999] {
            let network = TransportNetwork::new(cx,&flow,air(),vec![heated("wall",0.8)],
                &[inlet(0,300.0)],config()).unwrap()
                .with_recirculation(cx,vec![link(0,1,fraction)],1e-9).unwrap();
            let result = network.march(cx,&[350.0]).unwrap();
            let a = (-0.8_f64).exp();
            let intake = 300.0 + fraction * (1.0-a) * 50.0 / (1.0-fraction*a);
            let exhaust = a*intake+(1.0-a)*350.0;
            near(result.node_temperatures_k[0].unwrap(),intake,2e-8);
            near(result.node_temperatures_k[1].unwrap(),exhaust,2e-8);
            near(result.wall_heat_rate_w,(1.0-fraction)*(exhaust-300.0),2e-8);
            near(network.initial_references(cx).unwrap()[0],300.0,2e-8);
            if fraction > 0.0 {
                let report = result.recirculation.unwrap();
                near(report.external_heat_gain_w,result.wall_heat_rate_w,1e-7);
                near(report.supplies[0].mixed_temperature_k,intake,2e-8);
                assert!(report.max_mixing_residual_k <= 1e-9);
            } else { assert!(result.recirculation.is_none()); }
        }
    });
}

#[test]
fn multi_supply_return_mixes_by_makeup_capacity_and_has_zero_unheated_gain() {
    with_cx(|cx| {
        let graph = LossGraph::new(3,vec![edge("a",0,2,1.0),edge("b",1,2,1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0,9.0),boundary(1,4.0),boundary(2,0.0)],hydraulic_config(),cx).unwrap();
        let build = |links| TransportNetwork::new(cx,&flow,air(),vec![BranchThermalModel::Adiabatic;2],
            &[inlet(0,290.0),inlet(1,340.0)],config()).unwrap()
            .with_recirculation(cx,links,1e-9).unwrap();
        let a = build(vec![link(0,2,0.5),link(1,2,0.25)]).march(cx,&[]).unwrap();
        let b = build(vec![link(1,2,0.25),link(0,2,0.5)]).march(cx,&[]).unwrap();
        assert_eq!(a,b);
        // Both FRESH streams are 1.5 W/K; return temperature is 315 K.
        near(a.node_temperatures_k[2].unwrap(),315.0,1e-8);
        near(a.node_temperatures_k[0].unwrap(),302.5,1e-8);
        near(a.node_temperatures_k[1].unwrap(),333.75,1e-8);
        near(a.recirculation.unwrap().external_heat_gain_w,0.0,1e-7);
    });
}

fn value(march: &TransportMarch, w: &TransportObjective) -> f64 {
    march.node_temperatures_k.iter().zip(&w.node_temperatures).map(|(v,w)|v.unwrap_or(0.0)*w).sum::<f64>()
        + march.branches.iter().zip(&w.branch_outlets).map(|(v,w)|v.outlet_temperature_k.unwrap_or(0.0)*w).sum::<f64>()
        + march.reference_temperatures_k.iter().zip(&w.references).map(|(v,w)|v*w).sum::<f64>()
        + w.wall_heat_rate*march.wall_heat_rate_w + w.external_heat_gain*march.external_heat_gain_w
        + w.heat_imbalance*march.heat_imbalance_w + w.hydraulic_energy_defect*march.hydraulic_energy_defect_w
}
fn tangent_value(d: &TransportDifferential, w: &TransportObjective) -> f64 {
    d.node_temperatures_k.iter().zip(&w.node_temperatures).map(|(v,w)|v.unwrap_or(0.0)*w).sum::<f64>()
        + d.branch_outlets_k.iter().zip(&w.branch_outlets).map(|(v,w)|v.unwrap_or(0.0)*w).sum::<f64>()
        + d.reference_temperatures_k.iter().zip(&w.references).map(|(v,w)|v*w).sum::<f64>()
        + w.wall_heat_rate*d.wall_heat_rate_w + w.external_heat_gain*d.external_heat_gain_w
        + w.heat_imbalance*d.heat_imbalance_w + w.hydraulic_energy_defect*d.hydraulic_energy_defect_w
}

#[test]
fn feedback_tangent_and_transpose_match_full_perturbations_of_all_controls() {
    with_cx(|cx| {
        let graph = LossGraph::new(3,vec![edge("a",0,2,1.0),edge("b",1,2,1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0,9.0),boundary(1,4.0),boundary(2,0.0)],hydraulic_config(),cx).unwrap();
        let build = |step: f64| TransportNetwork::new(cx,&flow,air(),
            vec![heated("a-wall",1.2*(0.2*step).exp()),heated("b-wall",0.7*(-0.4*step).exp())],
            &[inlet(0,290.0+step*0.9),inlet(1,310.0-step*0.6)],config()).unwrap()
            .with_recirculation(cx,vec![link(0,2,0.6),link(1,2,0.3)],1e-9).unwrap();
        let network = build(0.0);
        let lin = network.linearize(cx,&[340.0,320.0]).unwrap();
        let mut direction = lin.zero_direction();
        direction.walls_k = vec![0.3,-0.7]; direction.inlets_k = vec![0.9,-0.6,0.0];
        direction.log_conductances = vec![0.2,-0.4];
        let mut objective = lin.zero_objective();
        objective.node_temperatures = vec![0.2,0.1,0.7];
        objective.branch_outlets = vec![-0.2,0.3]; objective.references = vec![0.4,-0.1];
        objective.wall_heat_rate = 0.12; objective.external_heat_gain = -0.07;
        objective.heat_imbalance = 0.03; objective.hydraulic_energy_defect = 0.09;
        let tangent = lin.apply(cx,&direction).unwrap();
        let gradient = lin.pullback(cx,&objective).unwrap();
        let dual = gradient.walls.iter().zip(&direction.walls_k).map(|(a,b)|a*b).sum::<f64>()
            + gradient.inlets.iter().zip(&direction.inlets_k).map(|(a,b)|a*b).sum::<f64>()
            + gradient.log_conductances.iter().zip(&direction.log_conductances).map(|(a,b)|a*b).sum::<f64>();
        near(tangent_value(&tangent,&objective),dual,1e-10);
        let eps = 1e-4;
        let plus = build(eps).march(cx,&[340.0+0.3*eps,320.0-0.7*eps]).unwrap();
        let minus = build(-eps).march(cx,&[340.0-0.3*eps,320.0+0.7*eps]).unwrap();
        near(dual,(value(&plus,&objective)-value(&minus,&objective))/(2.0*eps),2e-6);
        let mut uniform = lin.zero_direction(); uniform.walls_k.fill(1.0);
        uniform.inlets_k[0]=1.0; uniform.inlets_k[1]=1.0;
        for t in lin.apply(cx,&uniform).unwrap().node_temperatures_k.into_iter().flatten() { near(t,1.0,1e-12); }
    });
}

#[test]
fn return_admission_refuses_bad_ownership_overdraw_and_fully_closed_loops() {
    with_cx(|cx| {
        let graph = LossGraph::new(3,vec![edge("small",0,1,1.0),edge("large",0,2,1.0/9.0)]).unwrap();
        let flow = graph.solve(&[boundary(0,1.0),boundary(1,0.0),boundary(2,0.0)],hydraulic_config(),cx).unwrap();
        let build = || TransportNetwork::new(cx,&flow,air(),vec![BranchThermalModel::Adiabatic;2],
            &[inlet(0,300.0)],config()).unwrap();
        for links in [vec![link(0,1,0.5)], vec![link(1,0,0.1)], vec![link(0,5,0.1)],
            vec![link(0,1,-0.1)], vec![link(0,1,f64::NAN)], vec![link(0,2,1.0)],
            vec![link(0,1,0.25),link(0,2,0.75)],vec![link(0,2,0.1);2],
            vec![link(0,2,0.0);MAX_RECIRCULATION_LINKS+1]] {
            assert!(build().with_recirculation(cx,links,1e-9).is_err());
        }
        let normal = build().march(cx,&[]).unwrap();
        let zero = build().with_recirculation(cx,vec![link(0,2,0.0)],1e-9).unwrap().march(cx,&[]).unwrap();
        assert_eq!(normal,zero);
    });
}

#[test]
fn conductance_rebinding_and_reversed_branches_retain_the_return_model() {
    with_cx(|cx| {
        let graph = LossGraph::new(2,vec![edge("reverse",1,0,1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0,1.0),boundary(1,0.0)],hydraulic_config(),cx).unwrap();
        let network = TransportNetwork::new(cx,&flow,air(),vec![BranchThermalModel::Exchange(vec![
            AirSegment::new("last",1.0,0.5).unwrap(),AirSegment::new("first",1.0,0.3).unwrap()])],
            &[inlet(0,300.0)],config()).unwrap().with_recirculation(cx,vec![link(0,1,0.4)],1e-9).unwrap();
        let scaled = network.scaled_conductance(cx,2.0).unwrap();
        assert_eq!(scaled.recirculation_links(),network.recirculation_links());
        assert_eq!(scaled.regions(),vec!["first","last"]);
        let result = scaled.march(cx,&[350.0,350.0]).unwrap();
        let a = (-1.6_f64).exp();
        let intake = 300.0+0.4*(1.0-a)*50.0/(1.0-0.4*a);
        near(result.node_temperatures_k[1].unwrap(),a*intake+(1.0-a)*350.0,1e-8);
    });
}

#[test]
fn reduced_solve_pivots_and_transpose_preserve_duality_and_refuse_singular_data() {
    with_cx(|cx| {
        let a = vec![vec![0.0,2.0],vec![3.0,4.0]];
        let x = solve_feedback(cx,&a,&[6.0,18.0],false).unwrap();
        near(x[0],2.0,1e-12); near(x[1],3.0,1e-12);
        let y = solve_feedback(cx,&a,&[0.7,-0.2],true).unwrap();
        near(0.7*x[0]-0.2*x[1],6.0*y[0]+18.0*y[1],1e-12);
        assert!(solve_feedback(cx,&[vec![0.0]],&[1.0],false).is_err());
        assert!(solve_feedback(cx,&[vec![1.0]],&[f64::NAN],false).is_err());
    });
}

#[test]
fn cancellation_is_checked_before_reduced_solve_work() {
    let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate,arena,fs_exec::StreamKey {seed:1,kernel_id:711,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        assert!(matches!(solve_feedback(&cx,&[vec![1.0]],&[1.0],false),Err(TransportError::Interrupted)));
    });
}

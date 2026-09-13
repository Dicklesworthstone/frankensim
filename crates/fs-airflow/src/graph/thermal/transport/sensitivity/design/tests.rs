use super::*;
use fs_qty::{Density, VolumetricFlowRate};
use crate::graph::{GraphSolution, LossGraph};
use crate::graph::tests::{boundary, config as flow_config, edge, with_cx};
use crate::graph::thermal::transport::{TransportAir, TransportConfig};

fn request() -> UniformCoolingRequest {
    UniformCoolingRequest { node: 1, wall_temperature: Temperature::new(290.0), outlet_limit: Temperature::new(310.0),
        minimum_scale: 0.1, maximum_scale: 20.0, temperature_tolerance_k: 1.0e-7,
        log_scale_tolerance: 1.0e-8, max_evaluations: 100 }
}
fn flow(reverse: bool) -> GraphSolution {
    let (a,b) = if reverse { (1,0) } else { (0,1) };
    with_cx(|cx| LossGraph::new(2,vec![edge("duct",a,b,1.0)]).unwrap()
        .solve(&[boundary(0,4.0),boundary(1,0.0)],flow_config(),cx)).unwrap()
}
fn network<'a>(cx: &Cx<'_>, flow: &'a GraphSolution, source: f64, model: Vec<BranchThermalModel>) -> TransportNetwork<'a> {
    TransportNetwork::new(cx,flow,TransportAir {density:Density::new(1.0),specific_heat_j_kg_k:2.0},model,
        &[TransportInlet {node:0,temperature:Temperature::new(source)}],
        TransportConfig {absolute_flow_tolerance:VolumetricFlowRate::new(1e-8),relative_flow_tolerance:1e-8,
            absolute_heat_tolerance_w:1e-5,relative_heat_tolerance:1e-8}).unwrap()
}
fn models() -> Vec<BranchThermalModel> {
    vec![BranchThermalModel::Exchange(vec![AirSegment::new("cooler",1.0,0.5).unwrap()])]
}

#[test]
fn sizes_single_exchanger_to_closed_form_and_returns_the_passing_side() {
    let flow=flow(false);
    with_cx(|cx| {
        let base=network(cx,&flow,340.0,models()); let r=request();
        let sized=base.size_uniform_cooling(cx,r).unwrap();
        let exact=8.0*(50.0_f64/20.0).ln();
        assert!((sized.scale-exact).abs()<2e-7);
        assert!(sized.outlet_temperature.value()<=r.outlet_limit.value());
        assert!(sized.lower_outlet_temperature.value()>r.outlet_limit.value());
        assert!((sized.scale/sized.lower_scale).ln()<=1.01*r.log_scale_tolerance);
        assert!(310.0-sized.outlet_temperature.value()<=r.temperature_tolerance_k);
        assert!(sized.evaluations<=r.max_evaluations && !sized.at_lower_bound);
        let expected_slope=-(sized.outlet_temperature.value()-290.0)*sized.scale/8.0;
        assert!((sized.slope_per_log_scale_k-expected_slope).abs()<1e-9);
        assert_eq!(base.scaled_conductance(cx,sized.scale).unwrap().march(cx,&[290.0]).unwrap(),sized.march);
    });
}

#[test]
fn lower_bound_success_unattainable_limit_and_work_budget_are_distinct() {
    let flow=flow(false);
    with_cx(|cx| {
        let base=network(cx,&flow,340.0,models());
        let easy=base.size_uniform_cooling(cx,UniformCoolingRequest {outlet_limit:Temperature::new(340.0),..request()}).unwrap();
        assert!(easy.at_lower_bound); assert_eq!(easy.scale,request().minimum_scale); assert_eq!(easy.evaluations,1);
        assert!(matches!(base.size_uniform_cooling(cx,UniformCoolingRequest {outlet_limit:Temperature::new(291.0),..request()}),
            Err(CoolingSizingError::Unattainable {..})));
        assert!(matches!(base.size_uniform_cooling(cx,UniformCoolingRequest {max_evaluations:2,..request()}),
            Err(CoolingSizingError::BudgetExhausted {evaluations:2,..})));
    });
}

#[test]
fn bypass_floor_is_not_mistaken_for_a_reachable_temperature() {
    let graph=LossGraph::new(3,vec![edge("cooler",0,1,1.0),edge("bypass",0,1,4.0),edge("exhaust",1,2,5.0/9.0)]).unwrap();
    let flow=with_cx(|cx| graph.solve(&[boundary(0,9.0),boundary(2,0.0)],flow_config(),cx)).unwrap();
    with_cx(|cx| {
        let base=network(cx,&flow,340.0,vec![models().remove(0),BranchThermalModel::Adiabatic,BranchThermalModel::Adiabatic]);
        let failed=base.size_uniform_cooling(cx,UniformCoolingRequest {node:2,maximum_scale:1.0e4,
            outlet_limit:Temperature::new(300.0),..request()});
        match failed { Err(CoolingSizingError::Unattainable {outlet_temperature,..}) => {
            assert!((outlet_temperature.value()-(290.0+50.0/3.0)).abs()<1e-6);
        },other=>panic!("expected an unattainable bypass floor, got {other:?}") }
        let target=base.size_uniform_cooling(cx,UniformCoolingRequest {node:2,maximum_scale:100.0,..request()}).unwrap();
        let exact=8.0*(10.0_f64).ln();
        assert!((target.scale-exact).abs()<1e-6);
    });
}

#[test]
fn scaling_keeps_reverse_flow_region_order_and_changes_actual_exchanger_values() {
    let flow=flow(true);
    with_cx(|cx| {
        let base=network(cx,&flow,340.0,vec![BranchThermalModel::Exchange(vec![
            AirSegment::new("declared-first",1.0,0.2).unwrap(),AirSegment::new("declared-last",1.0,0.3).unwrap()])]);
        assert_eq!(base.regions(),vec!["declared-last","declared-first"]);
        let scaled=base.scaled_conductance(cx,2.0).unwrap(); assert_eq!(scaled.regions(),base.regions());
        let m=scaled.march(cx,&[300.0,290.0]).unwrap();
        let mid=300.0+40.0*(-0.6_f64/4.0).exp(); let out=290.0+(mid-290.0)*(-0.4_f64/4.0).exp();
        assert!((m.node_temperatures_k[1].unwrap()-out).abs()<1e-10);
        assert_eq!(base.scaled_conductance(cx,1.0).unwrap().march(cx,&[300.0,290.0]).unwrap(),base.march(cx,&[300.0,290.0]).unwrap());
    });
}

#[test]
fn cooling_assumptions_and_invalid_requests_refuse_before_search() {
    let flow=flow(false);
    with_cx(|cx| {
        let warm=network(cx,&flow,340.0,models()); let cold=network(cx,&flow,280.0,models());
        assert!(matches!(cold.size_uniform_cooling(cx,request()),Err(CoolingSizingError::SupplyBelowWall {node:0})));
        assert!(matches!(warm.size_uniform_cooling(cx,UniformCoolingRequest {node:9,..request()}),Err(CoolingSizingError::UnknownTarget {node:9})));
        for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
            assert!(warm.scaled_conductance(cx,bad).is_err());
            assert!(warm.size_uniform_cooling(cx,UniformCoolingRequest {log_scale_tolerance:bad,..request()}).is_err());
        }
        assert!(warm.size_uniform_cooling(cx,UniformCoolingRequest {minimum_scale:25.0,..request()}).is_err());
        assert!(warm.size_uniform_cooling(cx,UniformCoolingRequest {max_evaluations:1,..request()}).is_err());
    });
}

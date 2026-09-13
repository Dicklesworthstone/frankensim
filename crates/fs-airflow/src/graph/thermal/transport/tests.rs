//! G1/G3 transport checks through the real hydraulic solver, with independent
//! exponential-law and enthalpy-mixing oracles. No mocked transport solver.

use super::*;
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use crate::graph::{GraphSolution, LossGraph};

fn air() -> TransportAir {
    TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 1.0 }
}
fn budgets() -> TransportConfig {
    TransportConfig {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-9), relative_flow_tolerance: 1e-9,
        absolute_heat_tolerance_w: 1e-7, relative_heat_tolerance: 1e-8,
    }
}
fn inlet(node: usize, kelvin: f64) -> TransportInlet {
    TransportInlet { node, temperature: Temperature::new(kelvin) }
}
fn heated(name: &str, conductance: f64) -> BranchThermalModel {
    BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, conductance).unwrap()])
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-7, "actual {actual:.16e}, expected {expected:.16e}");
}
fn merge_flow(cx: &Cx<'_>) -> GraphSolution {
    // p=(2,1,0), parallel capacities 3 and 1, downstream capacity 4.
    LossGraph::new(3, vec![edge("left", 0, 1, 1.0 / 9.0),
        edge("right", 0, 1, 1.0), edge("downstream", 1, 2, 1.0 / 16.0)])
        .unwrap().solve(&[boundary(0, 2.0), boundary(2, 0.0)], hydraulic_config(), cx).unwrap()
}

#[test]
fn split_merge_and_downstream_heating_match_closed_form() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let net = TransportNetwork::new(cx, &flow, air(),
            vec![heated("left-wall", 3.0), heated("right-wall", 2.0), heated("tail-wall", 4.0)],
            &[inlet(0, 300.0)], budgets()).unwrap();
        let result = net.march(cx, &[340.0, 360.0, 320.0]).unwrap();
        let left = 340.0 - 40.0 * (-1.0_f64).exp();
        let right = 360.0 - 60.0 * (-2.0_f64).exp();
        let mixed = (3.0 * left + right) / 4.0;
        let outlet = 320.0 + (mixed - 320.0) * (-1.0_f64).exp();
        close(result.node_temperatures_k[1].unwrap(), mixed);
        close(result.node_temperatures_k[2].unwrap(), outlet);
        close(result.branches[2].inlet_temperature_k.unwrap(), mixed);
        close(result.wall_heat_rate_w, 4.0 * (outlet - 300.0));
        close(result.external_heat_gain_w, result.wall_heat_rate_w);
        assert_eq!(net.regions(), vec!["left-wall", "right-wall", "tail-wall"]);
        assert_eq!(net.initial_references(cx).unwrap(), vec![300.0; 3]);
    });
}

#[test]
fn separate_sources_mix_by_flow_not_arithmetic_mean() {
    with_cx(|cx| {
        let graph = LossGraph::new(3, vec![edge("cold", 0, 2, 1.0), edge("hot", 1, 2, 1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0, 9.0), boundary(1, 4.0), boundary(2, 0.0)], hydraulic_config(), cx).unwrap();
        let net = TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Adiabatic; 2],
            &[inlet(0, 290.0), inlet(1, 340.0)], budgets()).unwrap();
        let result = net.march(cx, &[]).unwrap();
        close(result.node_temperatures_k[2].unwrap(), 310.0);
        close(result.external_heat_gain_w, 0.0);
        assert!((result.node_temperatures_k[2].unwrap() - 315.0).abs() > 4.0);
    });
}

#[test]
fn reservoir_injection_mixes_with_arriving_air_but_extraction_does_not_reset_it() {
    with_cx(|cx| {
        for (resistances, inlets, expected) in [
            ([1.0, 1.0 / 9.0], vec![inlet(0, 330.0), inlet(1, 300.0)], 310.0),
            ([1.0 / 9.0, 1.0], vec![inlet(0, 330.0)], 330.0),
        ] {
            let graph = LossGraph::new(3, vec![edge("in", 0, 1, resistances[0]), edge("out", 1, 2, resistances[1])]).unwrap();
            let flow = graph.solve(&[boundary(0, 2.0), boundary(1, 1.0), boundary(2, 0.0)], hydraulic_config(), cx).unwrap();
            let net = TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Adiabatic; 2], &inlets, budgets()).unwrap();
            let result = net.march(cx, &[]).unwrap();
            close(result.node_temperatures_k[1].unwrap(), expected);
            close(result.node_temperatures_k[2].unwrap(), expected);
            close(result.external_heat_gain_w, 0.0);
        }
    });
}

#[test]
fn reverse_flow_reverses_the_physical_segment_order() {
    with_cx(|cx| {
        let graph = LossGraph::new(2, vec![edge("reversed", 1, 0, 1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0, 1.0), boundary(1, 0.0)], hydraulic_config(), cx).unwrap();
        let models = vec![BranchThermalModel::Exchange(vec![
            AirSegment::new("near-sink", 1.0, 0.7).unwrap(),
            AirSegment::new("near-source", 1.0, 0.3).unwrap(),
        ])];
        let net = TransportNetwork::new(cx, &flow, air(), models, &[inlet(0, 300.0)], budgets()).unwrap();
        assert_eq!(net.regions(), vec!["near-source", "near-sink"]);
        let result = net.march(cx, &[360.0, 320.0]).unwrap();
        let after_first = 360.0 - 60.0 * (-0.3_f64).exp();
        let expected = 320.0 + (after_first - 320.0) * (-0.7_f64).exp();
        close(result.node_temperatures_k[1].unwrap(), expected);
        assert_eq!(result.branches[0].march.as_ref().unwrap().segments[0].region, "near-source");
    });
}

#[test]
fn leakage_bypass_dilutes_the_heated_stream_before_downstream_cooling() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let net = TransportNetwork::new(cx, &flow, air(),
            vec![heated("heated-stream", 3.0), BranchThermalModel::Adiabatic, heated("tail", 4.0)],
            &[inlet(0, 300.0)], budgets()).unwrap();
        let result = net.march(cx, &[340.0, 310.0]).unwrap();
        let hot = 340.0 - 40.0 * (-1.0_f64).exp();
        let mixed = (3.0 * hot + 300.0) / 4.0;
        close(result.node_temperatures_k[1].unwrap(), mixed);
        close(result.node_temperatures_k[2].unwrap(), 310.0 + (mixed - 310.0) * (-1.0_f64).exp());
        assert!(result.branches[1].march.is_none());
        assert_eq!(result.branches[1].outlet_temperature_k, Some(300.0));
    });
}

#[test]
fn one_branch_reuses_the_existing_channel_law_exactly() {
    with_cx(|cx| {
        let graph = LossGraph::new(2, vec![edge("single", 0, 1, 1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0, 1.0), boundary(1, 0.0)], hydraulic_config(), cx).unwrap();
        let segment = AirSegment::new("wall", 1.0, 0.8).unwrap();
        let net = TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Exchange(vec![segment.clone()])],
            &[inlet(0, 300.0)], budgets()).unwrap();
        let expected = AirPath::new(300.0, flow.branches[0].flow.value(), 1.0, vec![segment]).unwrap().march(&[350.0]).unwrap();
        assert_eq!(net.march(cx, &[350.0]).unwrap().branches[0].march.as_ref(), Some(&expected));
    });
}

#[test]
fn stagnant_adiabatic_air_remains_unknown_and_stagnant_heat_exchange_refuses() {
    with_cx(|cx| {
        let graph = LossGraph::new(2, vec![edge("stagnant", 0, 1, 1.0)]).unwrap();
        let flow = graph.solve(&[boundary(0, 2.0), boundary(1, 2.0)], hydraulic_config(), cx).unwrap();
        let net = TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Adiabatic], &[], budgets()).unwrap();
        let result = net.march(cx, &[]).unwrap();
        assert_eq!(result.node_temperatures_k, vec![None, None]);
        assert_eq!(result.branches[0].outlet_temperature_k, None);
        assert!(matches!(TransportNetwork::new(cx, &flow, air(), vec![heated("wall", 1.0)], &[], budgets()),
            Err(TransportError::Branch { branch: 0, .. })));
    });
}

#[test]
fn incomplete_boundaries_duplicate_ownership_and_invalid_wall_values_refuse() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        for inlets in [vec![], vec![inlet(0, 300.0), inlet(0, 300.0)],
            vec![inlet(0, 300.0), inlet(2, 300.0)], vec![inlet(0, f64::NAN)]] {
            assert!(TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Adiabatic; 3], &inlets, budgets()).is_err());
        }
        assert!(TransportNetwork::new(cx, &flow, air(),
            vec![heated("same", 1.0), heated("same", 1.0), BranchThermalModel::Adiabatic],
            &[inlet(0, 300.0)], budgets()).is_err());
        let net = TransportNetwork::new(cx, &flow, air(),
            vec![heated("left", 1.0), BranchThermalModel::Adiabatic, BranchThermalModel::Adiabatic],
            &[inlet(0, 300.0)], budgets()).unwrap();
        for walls in [vec![], vec![f64::NAN], vec![0.0], vec![-1.0], vec![f64::INFINITY]] {
            assert!(net.march(cx, &walls).is_err());
        }
    });
}

#[test]
fn branch_flows_not_public_summary_fields_determine_admission() {
    with_cx(|cx| {
        let mut flow = merge_flow(cx);
        flow.branches[2].flow = VolumetricFlowRate::new(flow.branches[2].flow.value() * 1.01);
        let result = TransportNetwork::new(cx, &flow, air(), vec![BranchThermalModel::Adiabatic; 3], &[inlet(0, 300.0)], budgets());
        assert!(matches!(result, Err(TransportError::FlowImbalance { node: 1, .. })));
    });
}

#[test]
fn hydraulic_energy_defect_is_disclosed_but_never_corrected_away() {
    with_cx(|cx| {
        let mut flow = merge_flow(cx);
        flow.branches[2].flow = VolumetricFlowRate::new(flow.branches[2].flow.value() * 1.01);
        let mut config = budgets();
        config.absolute_flow_tolerance = VolumetricFlowRate::new(0.1);
        let net = TransportNetwork::new(cx, &flow, air(),
            vec![heated("left", 3.0), heated("right", 1.0), BranchThermalModel::Adiabatic],
            &[inlet(0, 300.0)], config).unwrap();
        match net.march(cx, &[350.0, 350.0]) {
            Err(TransportError::HeatImbalance { residual_w, hydraulic_defect_w, .. }) => {
                assert!(residual_w.abs() > 1.0);
                close(residual_w, hydraulic_defect_w);
            }
            other => panic!("unbalanced heat transport must refuse, got {other:?}"),
        }
    });
}

#[test]
fn declaration_permutation_and_repeated_marches_preserve_physical_results() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let net = TransportNetwork::new(cx, &flow, air(), vec![heated("a", 3.0), heated("b", 2.0), heated("c", 4.0)],
            &[inlet(0, 300.0)], budgets()).unwrap();
        let first = net.march(cx, &[340.0, 360.0, 320.0]).unwrap();
        assert_eq!(first, net.march(cx, &[340.0, 360.0, 320.0]).unwrap());
        let mut permuted = flow.clone();
        permuted.branches.swap(0, 1);
        let net = TransportNetwork::new(cx, &permuted, air(), vec![heated("b", 2.0), heated("a", 3.0), heated("c", 4.0)],
            &[inlet(0, 300.0)], budgets()).unwrap();
        let second = net.march(cx, &[360.0, 340.0, 320.0]).unwrap();
        for (a, b) in first.node_temperatures_k.iter().zip(&second.node_temperatures_k) {
            close(a.unwrap(), b.unwrap());
        }
    });
}

#[test]
fn cancellation_is_checked_even_for_an_empty_transport_pass() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 1, kernel_id: 712, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        let flow = merge_flow(&cx);
        let net = TransportNetwork::new(&cx, &flow, air(), vec![BranchThermalModel::Adiabatic; 3], &[inlet(0, 300.0)], budgets()).unwrap();
        gate.request();
        assert!(matches!(net.march(&cx, &[]), Err(TransportError::Interrupted)));
    });
}

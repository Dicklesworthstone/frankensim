//! Focused G0/G1/G3/G5 graph checks; numerical oracles are analytic.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

use super::*;
use crate::{LossResistance, SourceProvenance, ToleranceBasis};

pub(super) fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new(), f)
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate, arena,
            StreamKey { seed: 17, kernel_id: 711, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ))
    })
}

pub(super) fn config() -> GraphSolveConfig {
    GraphSolveConfig {
        max_sweeps: 4096,
        max_node_iterations: 80,
        absolute_flow_tolerance: VolumetricFlowRate::new(1.0e-10),
        relative_flow_tolerance: 1.0e-10,
    }
}

pub(super) fn edge(name: &str, from: usize, to: usize, resistance: f64) -> GraphBranch {
    GraphBranch {
        from, to,
        loss: LossElement::new(
            name, LossResistance::new(resistance), 0.15,
            SourceProvenance::new("analytic graph fixture", "graph-fixture-v1"),
            ToleranceBasis::EngineeringAllowance,
        ).unwrap(),
    }
}

pub(super) fn boundary(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure { node, pressure: Pressure::new(pressure) }
}

pub(super) fn bridge() -> LossGraph {
    LossGraph::new(4, vec![
        edge("feed", 0, 1, 1.0), edge("left", 1, 3, 4.0),
        edge("cross", 2, 1, 7.0), edge("right-feed", 0, 2, 16.0),
        edge("right", 2, 3, 2.25),
    ]).unwrap()
}

pub(super) fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2.0e-7 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}");
}

#[test]
fn g1_non_series_parallel_bridge_has_manufactured_solution() {
    let graph = bridge();
    let result = with_cx(|cx| graph.solve(&[boundary(0, 25.0), boundary(3, 0.0)], config(), cx)).unwrap();
    for (actual, expected) in result.pressures.iter().zip([25.0, 16.0, 9.0, 0.0]) {
        close(actual.value(), expected);
    }
    for (actual, expected) in result.branches.iter().zip([3.0, 2.0, -1.0, 1.0, 2.0]) {
        close(actual.flow.value(), expected);
        close(actual.pressure_drop.value(), actual.loss.resistance.value()
            * actual.flow.value() * actual.flow.value().abs());
    }
    close(result.node_outflows[0].value(), 4.0);
    close(result.node_outflows[3].value(), -4.0);
    assert!(result.max_node_imbalance.value() < 1.0e-8);
    assert!(result.boundary_imbalance.value().abs() < 1.0e-8);
    assert_eq!(result.branches[2].loss, graph.branches()[2].loss);
}

#[test]
fn g1_series_and_parallel_agree_with_closed_forms() {
    let fixtures = [
        (LossGraph::new(3, vec![edge("a", 0, 1, 1.0), edge("b", 1, 2, 4.0)]).unwrap(), 2, 5.0, 1.0),
        (LossGraph::new(2, vec![edge("a", 0, 1, 4.0), edge("b", 0, 1, 9.0)]).unwrap(), 1, 36.0, 5.0),
        (LossGraph::new(2, vec![edge("a", 0, 1, 4.0)]).unwrap(), 1, 9.0, 1.5),
    ];
    for (graph, outlet, pressure, flow) in fixtures {
        let result = with_cx(|cx| graph.solve(&[boundary(0, pressure), boundary(outlet, 0.0)], config(), cx)).unwrap();
        close(result.node_outflows[0].value(), flow);
    }
}

#[test]
fn g1_balanced_bridge_and_dead_leg_allow_exact_zero_flow() {
    let graph = LossGraph::new(5, vec![
        edge("a", 0, 1, 1.0), edge("b", 1, 3, 1.0),
        edge("c", 0, 2, 1.0), edge("d", 2, 3, 1.0),
        edge("cross", 1, 2, 1.0), edge("dead-leg", 1, 4, 1.0),
    ]).unwrap();
    let result = with_cx(|cx| graph.solve(&[boundary(0, 4.0), boundary(3, 0.0)], config(), cx)).unwrap();
    assert_eq!(result.branches[4].flow.value(), 0.0);
    assert_eq!(result.branches[5].flow.value(), 0.0);
    close(result.node_outflows[0].value(), 2.0 * 2.0_f64.sqrt());
}

#[test]
fn g3_pressure_scale_and_reference_shift_preserve_the_expected_flows() {
    let graph = bridge();
    for (pressure, shift, factor) in [(25.0, 5.0, 1.0), (100.0, 0.0, 2.0), (25.0, -30.0, 1.0)] {
        let result = with_cx(|cx| graph.solve(&[boundary(0, pressure + shift), boundary(3, shift)], config(), cx)).unwrap();
        for (actual, expected) in result.branches.iter().zip([3.0, 2.0, -1.0, 1.0, 2.0]) {
            close(actual.flow.value(), factor * expected);
        }
    }
}

#[test]
fn g3_reversing_edge_orientation_only_changes_the_flow_sign() {
    let original = bridge();
    let mut edges = original.branches().to_vec();
    for edge in &mut edges { std::mem::swap(&mut edge.from, &mut edge.to); }
    let reversed = LossGraph::new(4, edges).unwrap();
    let boundaries = [boundary(0, 25.0), boundary(3, 0.0)];
    let a = with_cx(|cx| original.solve(&boundaries, config(), cx)).unwrap();
    let b = with_cx(|cx| reversed.solve(&boundaries, config(), cx)).unwrap();
    for (a, b) in a.branches.iter().zip(&b.branches) { close(a.flow.value(), -b.flow.value()); }
}

#[test]
fn g0_every_component_needs_a_pressure_reference() {
    let graph = LossGraph::new(4, vec![edge("a", 0, 1, 1.0), edge("b", 2, 3, 4.0)]).unwrap();
    let refused = with_cx(|cx| graph.solve(&[boundary(0, 1.0)], config(), cx));
    assert!(matches!(refused, Err(GraphError::UnreferencedNode(2))));
    let result = with_cx(|cx| graph.solve(&[boundary(0, 1.0), boundary(2, 8.0)], config(), cx)).unwrap();
    assert!(result.branches.iter().all(|edge| edge.flow.value() == 0.0));
}

#[test]
fn g1_three_pressure_reservoirs_balance_at_the_junction() {
    let graph = LossGraph::new(4, vec![
        edge("a", 0, 3, 1.0), edge("b", 1, 3, 1.0), edge("c", 3, 2, 1.0),
    ]).unwrap();
    let result = with_cx(|cx| graph.solve(&[boundary(0, 9.0), boundary(1, 4.0), boundary(2, 0.0)], config(), cx)).unwrap();
    let p = result.pressures[3].value();
    close((9.0 - p).sqrt() - (p - 4.0).sqrt() - p.sqrt(), 0.0);
    close(result.boundary_imbalance.value(), 0.0);
}

#[test]
fn g0_invalid_topology_and_mutated_coefficients_refuse() {
    for branches in [vec![edge("a", 0, 0, 1.0)], vec![edge("a", 0, 2, 1.0)],
        vec![edge("a", 0, 1, 1.0), edge("a", 1, 0, 2.0)]] {
        assert!(LossGraph::new(2, branches).is_err());
    }
    for resistance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut branch = edge("mutated", 0, 1, 1.0);
        branch.loss.resistance = LossResistance::new(resistance);
        assert!(LossGraph::new(2, vec![branch]).is_err());
    }
}

#[test]
fn g0_invalid_boundaries_and_budgets_refuse() {
    let graph = bridge();
    for boundaries in [vec![], vec![boundary(4, 0.0)], vec![boundary(0, f64::NAN)],
        vec![boundary(0, 1.0), boundary(0, 1.0)]] {
        assert!(with_cx(|cx| graph.solve(&boundaries, config(), cx)).is_err());
    }
    for bad in [0.0, f64::NAN, f64::INFINITY] {
        let mut budget = config();
        budget.absolute_flow_tolerance = VolumetricFlowRate::new(bad);
        assert!(with_cx(|cx| graph.solve(&[boundary(0, 25.0), boundary(3, 0.0)], budget, cx)).is_err());
    }
}

#[test]
fn g0_exhausted_budget_does_not_publish_a_partial_solution() {
    let mut budget = config();
    budget.max_sweeps = 1;
    let result = with_cx(|cx| bridge().solve(&[boundary(0, 25.0), boundary(3, 0.0)], budget, cx));
    assert!(matches!(result, Err(GraphError::DidNotConverge { sweeps: 1, .. })));
}

#[test]
fn g4_pre_requested_cancellation_refuses() {
    let gate = CancelGate::new();
    gate.request();
    let result = with_gate(&gate, |cx| bridge().solve(&[boundary(0, 25.0), boundary(3, 0.0)], config(), cx));
    assert!(matches!(result, Err(GraphError::Interrupted)));
}

#[test]
fn g5_same_input_repeats_bit_for_bit() {
    let graph = bridge();
    let boundaries = [boundary(0, 25.0), boundary(3, 0.0)];
    let a = with_cx(|cx| graph.solve(&boundaries, config(), cx)).unwrap();
    let b = with_cx(|cx| graph.solve(&boundaries, config(), cx)).unwrap();
    assert_eq!(a.sweeps, b.sweeps);
    for (a, b) in a.pressures.iter().zip(&b.pressures) { assert_eq!(a.value().to_bits(), b.value().to_bits()); }
    for (a, b) in a.branches.iter().zip(&b.branches) { assert_eq!(a.flow.value().to_bits(), b.flow.value().to_bits()); }
}

#[test]
fn g0_separate_square_roots_avoid_avoidable_ratio_overflow() {
    let graph = LossGraph::new(2, vec![edge("extreme", 0, 1, 1.0e-200)]).unwrap();
    let result = with_cx(|cx| graph.solve(&[boundary(0, 1.0e200), boundary(1, 0.0)], config(), cx)).unwrap();
    close(result.branches[0].flow.value() / 1.0e200, 1.0);
}

#[test]
fn g1_asymmetric_dead_legs_do_not_stall_point_relaxation() {
    for (graph, inlet, outlet, branch_index) in [
        (LossGraph::new(4, vec![edge("a", 0, 1, 1.0), edge("b", 1, 2, 4.0), edge("dead", 1, 3, 1.0)]).unwrap(), 0, 2, 2),
        (LossGraph::new(4, vec![edge("a", 1, 2, 1.0), edge("b", 2, 3, 4.0), edge("dead", 2, 0, 1.0)]).unwrap(), 1, 3, 2),
    ] {
        let result = with_cx(|cx| graph.solve(&[boundary(inlet, 5.0), boundary(outlet, 0.0)], config(), cx)).unwrap();
        close(result.node_outflows[inlet].value(), 1.0);
        assert_eq!(result.branches[branch_index].flow.value(), 0.0);
        assert!(result.sweeps < 100);
    }
}

use super::*;
use super::super::tests::{boundary, bridge, config, edge, with_cx};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use crate::LossResistance;

fn sensitivity() -> GraphSensitivityConfig {
    GraphSensitivityConfig {
        minimum_absolute_flow: VolumetricFlowRate::new(1e-12),
        linear: CheckedCgConfig { tolerance: 1e-12, max_iterations: 128, max_corrections: 2 },
    }
}

fn primal() -> GraphSolveConfig {
    let mut value = config();
    value.absolute_flow_tolerance = VolumetricFlowRate::new(1e-12);
    value.relative_flow_tolerance = 1e-12;
    value
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}");
}

fn objective(solution: &GraphSolution, wp: &[PressureWeight], wq: &[FlowWeight]) -> f64 {
    solution.pressures.iter().zip(wp).map(|(p, w)| p.value() * w.value()).sum::<f64>()
        + solution.branches.iter().zip(wq).map(|(q, w)| q.flow.value() * w.value()).sum::<f64>()
}

fn perturbed(
    graph: &LossGraph, boundaries: &[FixedPressure], z: &[f64], dp: &[Pressure],
    step: f64, cx: &Cx<'_>,
) -> GraphSolution {
    let mut branches = graph.branches().to_vec();
    for (branch, &direction) in branches.iter_mut().zip(z) {
        branch.loss.resistance = LossResistance::new(branch.loss.resistance.value() * (1.0 + step * direction));
    }
    let changed = LossGraph::new(graph.node_count(), branches).unwrap();
    let boundaries: Vec<_> = boundaries.iter().zip(dp).map(|(b, d)|
        boundary(b.node, b.pressure.value() + step * d.value())).collect();
    changed.solve(&boundaries, primal(), cx).unwrap()
}

#[test]
fn all_fixed_single_branch_has_analytic_tangent_and_gradient() {
    let graph = LossGraph::new(2, vec![edge("duct", 0, 1, 4.0)]).unwrap();
    with_cx(|cx| {
        let linear = graph.solve_linearized(&[boundary(0, 9.0), boundary(1, 0.0)], primal(), sensitivity(), cx).unwrap();
        let tangent = linear.tangent(&[0.2], &[Pressure::new(2.0), Pressure::new(-1.0)], cx).unwrap();
        close(tangent.flows[0].value(), 0.1, 1e-12);
        close(tangent.node_outflows[0].value(), 0.1, 1e-12);
        close(tangent.node_outflows[1].value(), -0.1, 1e-12);
        assert_eq!(tangent.max_node_imbalance.value(), 0.0);
        let gradient = linear.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx).unwrap();
        close(gradient.log_resistances[0], -0.75, 1e-12);
        close(gradient.boundary_pressures[0].value(), 1.0 / 12.0, 1e-12);
        close(gradient.boundary_pressures[1].value(), -1.0 / 12.0, 1e-12);
        assert!(gradient.linear_report.converged_euclidean());
        assert_eq!(gradient.linear_report.iters, 0);
        assert_eq!(linear.solution().branches[0].loss, graph.branches()[0].loss);
    });
}

#[test]
fn series_pressure_objective_includes_indirect_and_boundary_derivatives() {
    let graph = LossGraph::new(3, vec![edge("a", 0, 1, 1.0), edge("b", 1, 2, 4.0)]).unwrap();
    with_cx(|cx| {
        let linear = graph.solve_linearized(&[boundary(0, 5.0), boundary(2, 0.0)], primal(), sensitivity(), cx).unwrap();
        let tangent = linear.tangent(&[1.0, 0.0], &[Pressure::new(0.0); 2], cx).unwrap();
        close(tangent.pressures[1].value(), -0.8, 1e-10);
        for flow in &tangent.flows { close(flow.value(), -0.1, 1e-10); }
        let gradient = linear.adjoint(&[PressureWeight::new(0.0), PressureWeight::new(1.0), PressureWeight::new(0.0)],
            &[FlowWeight::new(0.0); 2], cx).unwrap();
        close(gradient.log_resistances[0], -0.8, 1e-10);
        close(gradient.log_resistances[1], 0.8, 1e-10);
        close(gradient.boundary_pressures[0].value(), 0.8, 1e-10);
        close(gradient.boundary_pressures[1].value(), 0.2, 1e-10);
    });
}

#[test]
fn reverse_flow_bridge_passes_jvp_vjp_and_full_finite_difference_checks() {
    let graph = bridge();
    // Deliberately not in node order: boundary gradient order is part of the API.
    let boundaries = [boundary(3, 0.0), boundary(0, 25.0)];
    let z = [0.2, -0.4, 0.1, 0.5, -0.3];
    let dp = [Pressure::new(0.6), Pressure::new(-0.2)];
    let wp = [0.2, -0.3, 0.1, 0.4].map(PressureWeight::new);
    let wq = [1.0, -2.0, 0.7, 0.25, 3.0].map(FlowWeight::new);
    with_cx(|cx| {
        let linear = graph.solve_linearized(&boundaries, primal(), sensitivity(), cx).unwrap();
        assert!(linear.solution().branches[2].flow.value() < 0.0);
        let tangent = linear.tangent(&z, &dp, cx).unwrap();
        let gradient = linear.adjoint(&wp, &wq, cx).unwrap();
        let forward = tangent.pressures.iter().zip(&wp).map(|(p, w)| p.value()*w.value()).sum::<f64>()
            + tangent.flows.iter().zip(&wq).map(|(q, w)| q.value()*w.value()).sum::<f64>();
        let reverse = gradient.log_resistances.iter().zip(z).map(|(g, d)| g*d).sum::<f64>()
            + gradient.boundary_pressures.iter().zip(dp).map(|(g, d)| g.value()*d.value()).sum::<f64>();
        close(forward, reverse, 1e-10);
        assert!(tangent.max_node_imbalance.value() < 1e-10);
        let h = 1e-4;
        let plus = perturbed(&graph, &boundaries, &z, &dp, h, cx);
        let minus = perturbed(&graph, &boundaries, &z, &dp, -h, cx);
        for node in 0..graph.node_count() {
            close(tangent.pressures[node].value(),
                (plus.pressures[node].value()-minus.pressures[node].value())/(2.0*h), 2e-7);
        }
        for edge in 0..graph.branches().len() {
            close(tangent.flows[edge].value(),
                (plus.branches[edge].flow.value()-minus.branches[edge].flow.value())/(2.0*h), 2e-7);
            let mut direction = vec![0.0; graph.branches().len()];
            direction[edge] = 1.0;
            let a = perturbed(&graph, &boundaries, &direction, &[Pressure::new(0.0); 2], h, cx);
            let b = perturbed(&graph, &boundaries, &direction, &[Pressure::new(0.0); 2], -h, cx);
            close(gradient.log_resistances[edge], (objective(&a, &wp, &wq)-objective(&b, &wp, &wq))/(2.0*h), 2e-7);
        }
        for index in 0..boundaries.len() {
            let mut direction = vec![Pressure::new(0.0); boundaries.len()];
            direction[index] = Pressure::new(1.0);
            let a = perturbed(&graph, &boundaries, &[0.0; 5], &direction, h, cx);
            let b = perturbed(&graph, &boundaries, &[0.0; 5], &direction, -h, cx);
            close(gradient.boundary_pressures[index].value(), (objective(&a, &wp, &wq)-objective(&b, &wp, &wq))/(2.0*h), 2e-7);
        }
    });
}

#[test]
fn reference_shift_and_uniform_loss_scaling_have_exact_structural_derivatives() {
    let graph = bridge();
    with_cx(|cx| {
        let linear = graph.solve_linearized(&[boundary(0, 25.0), boundary(3, 0.0)], primal(), sensitivity(), cx).unwrap();
        let shift = linear.tangent(&[0.0; 5], &[Pressure::new(1.0); 2], cx).unwrap();
        for p in shift.pressures { close(p.value(), 1.0, 1e-10); }
        for q in shift.flows { close(q.value(), 0.0, 1e-10); }
        let scale = linear.tangent(&[1.0; 5], &[Pressure::new(0.0); 2], cx).unwrap();
        for p in scale.pressures { close(p.value(), 0.0, 1e-10); }
        for (q, primal) in scale.flows.iter().zip(&linear.solution().branches) {
            close(q.value(), -0.5*primal.flow.value(), 1e-10);
        }
        let gradient = linear.adjoint(&[PressureWeight::new(0.0); 4], &[FlowWeight::new(1.0); 5], cx).unwrap();
        close(gradient.boundary_pressures.iter().map(|g| g.value()).sum(), 0.0, 1e-10);
    });
}

#[test]
fn orientation_change_preserves_physical_objective_gradient() {
    let graph = bridge();
    let mut branches = graph.branches().to_vec();
    for edge in &mut branches { std::mem::swap(&mut edge.from, &mut edge.to); }
    let reversed = LossGraph::new(4, branches).unwrap();
    with_cx(|cx| {
        let boundaries = [boundary(0, 25.0), boundary(3, 0.0)];
        let a = graph.solve_linearized(&boundaries, primal(), sensitivity(), cx).unwrap();
        let b = reversed.solve_linearized(&boundaries, primal(), sensitivity(), cx).unwrap();
        let ga = a.adjoint(&[PressureWeight::new(0.2); 4], &[FlowWeight::new(1.0); 5], cx).unwrap();
        let gb = b.adjoint(&[PressureWeight::new(0.2); 4], &[FlowWeight::new(-1.0); 5], cx).unwrap();
        for (a, b) in ga.log_resistances.iter().zip(gb.log_resistances) { close(*a, b, 1e-10); }
        for (a, b) in ga.boundary_pressures.iter().zip(gb.boundary_pressures) { close(a.value(), b.value(), 1e-10); }
    });
}

#[test]
fn disconnected_grounded_components_do_not_cross_couple_gradients() {
    let graph = LossGraph::new(5, vec![edge("a", 0, 1, 1.0), edge("b", 2, 3, 1.0), edge("c", 3, 4, 3.0)]).unwrap();
    with_cx(|cx| {
        let boundaries = [boundary(0, 9.0), boundary(1, 0.0), boundary(2, 16.0), boundary(4, 0.0)];
        let linear = graph.solve_linearized(&boundaries, primal(), sensitivity(), cx).unwrap();
        let gradient = linear.adjoint(&[PressureWeight::new(0.0); 5],
            &[FlowWeight::new(1.0), FlowWeight::new(0.0), FlowWeight::new(0.0)], cx).unwrap();
        close(gradient.log_resistances[0], -1.5, 1e-12);
        assert_eq!(&gradient.log_resistances[1..], &[0.0, 0.0]);
        assert_eq!(gradient.boundary_pressures[2].value(), 0.0);
        assert_eq!(gradient.boundary_pressures[3].value(), 0.0);
    });
}

#[test]
fn stagnation_is_a_derivative_refusal_not_a_primal_failure() {
    let graph = LossGraph::new(2, vec![edge("idle", 0, 1, 4.0)]).unwrap();
    with_cx(|cx| {
        let boundaries = [boundary(0, 0.0), boundary(1, 0.0)];
        assert_eq!(graph.solve(&boundaries, primal(), cx).unwrap().branches[0].flow.value(), 0.0);
        let mut cfg = sensitivity();
        cfg.minimum_absolute_flow = VolumetricFlowRate::new(0.0);
        assert!(matches!(graph.solve_linearized(&boundaries, primal(), cfg, cx),
            Err(GraphSensitivityError::NonDifferentiableBranch { branch: 0, .. })));
        cfg.minimum_absolute_flow = VolumetricFlowRate::new(2.0);
        assert!(matches!(graph.solve_linearized(&[boundary(0, 9.0), boundary(1, 0.0)], primal(), cfg, cx),
            Err(GraphSensitivityError::NonDifferentiableBranch { branch: 0, .. })));
    });
}

#[test]
fn malformed_directions_weights_and_budgets_refuse() {
    let graph = bridge();
    with_cx(|cx| {
        let boundaries = [boundary(0, 25.0), boundary(3, 0.0)];
        for value in [-1.0, f64::NAN, f64::INFINITY] {
            let mut cfg = sensitivity(); cfg.minimum_absolute_flow = VolumetricFlowRate::new(value);
            assert!(graph.solve_linearized(&boundaries, primal(), cfg, cx).is_err());
        }
        for value in [0.0, 1.0, f64::NAN, f64::INFINITY] {
            let mut cfg = sensitivity(); cfg.linear.tolerance = value;
            assert!(graph.solve_linearized(&boundaries, primal(), cfg, cx).is_err());
        }
        let mut cfg = sensitivity(); cfg.linear.max_iterations = 0;
        assert!(graph.solve_linearized(&boundaries, primal(), cfg, cx).is_err());
        let linear = graph.solve_linearized(&boundaries, primal(), sensitivity(), cx).unwrap();
        assert!(linear.tangent(&[], &[Pressure::new(0.0); 2], cx).is_err());
        assert!(linear.tangent(&[0.0; 5], &[], cx).is_err());
        assert!(linear.tangent(&[f64::NAN; 5], &[Pressure::new(0.0); 2], cx).is_err());
        assert!(linear.tangent(&[0.0; 5], &[Pressure::new(f64::INFINITY); 2], cx).is_err());
        assert!(linear.adjoint(&[], &[FlowWeight::new(0.0); 5], cx).is_err());
        assert!(linear.adjoint(&[PressureWeight::new(0.0); 4], &[], cx).is_err());
        assert!(linear.adjoint(&[PressureWeight::new(f64::NAN); 4], &[FlowWeight::new(0.0); 5], cx).is_err());
        assert!(linear.adjoint(&[PressureWeight::new(0.0); 4], &[FlowWeight::new(f64::INFINITY); 5], cx).is_err());
    });
}

#[test]
fn linear_budget_exhaustion_cannot_publish_a_gradient() {
    let graph = bridge();
    with_cx(|cx| {
        let mut cfg = sensitivity(); cfg.linear.max_iterations = 1; cfg.linear.max_corrections = 0;
        let linear = graph.solve_linearized(&[boundary(0, 25.0), boundary(3, 0.0)], primal(), cfg, cx).unwrap();
        let result = linear.adjoint(&[PressureWeight::new(0.0), PressureWeight::new(1.0),
            PressureWeight::new(0.0), PressureWeight::new(0.0)], &[FlowWeight::new(0.0); 5], cx);
        assert!(matches!(result, Err(GraphSensitivityError::LinearDidNotConverge(_))));
    });
}

#[test]
fn cancellation_after_primal_acceptance_refuses_both_derivative_directions() {
    let graph = bridge();
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 17, kernel_id: 711, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let linear = graph.solve_linearized(&[boundary(0, 25.0), boundary(3, 0.0)], primal(), sensitivity(), &cx).unwrap();
        gate.request();
        assert!(matches!(linear.tangent(&[0.0; 5], &[Pressure::new(0.0); 2], &cx),
            Err(GraphSensitivityError::Graph(GraphError::Interrupted))));
        assert!(matches!(linear.adjoint(&[PressureWeight::new(0.0); 4], &[FlowWeight::new(1.0); 5], &cx),
            Err(GraphSensitivityError::Graph(GraphError::Interrupted))));
    });
}

#[test]
fn repeating_a_gradient_is_bit_identical() {
    let graph = bridge();
    with_cx(|cx| {
        let linear = graph.solve_linearized(&[boundary(0, 25.0), boundary(3, 0.0)], primal(), sensitivity(), cx).unwrap();
        let wp = [PressureWeight::new(0.2); 4];
        let wq = [FlowWeight::new(0.7); 5];
        let a = linear.adjoint(&wp, &wq, cx).unwrap();
        let b = linear.adjoint(&wp, &wq, cx).unwrap();
        assert_eq!(a.linear_report.iters, b.linear_report.iters);
        for (a, b) in a.log_resistances.iter().zip(b.log_resistances) { assert_eq!(a.to_bits(), b.to_bits()); }
        for (a, b) in a.boundary_pressures.iter().zip(b.boundary_pressures) { assert_eq!(a.value().to_bits(), b.value().to_bits()); }
    });
}

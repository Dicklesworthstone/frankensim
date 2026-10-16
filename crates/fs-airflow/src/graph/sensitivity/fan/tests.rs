use super::*;
use super::super::super::tests::{bridge, config, edge, with_cx};
use crate::{FanArrangement, FanCurve, FanPoint, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solver::CheckedCgConfig;

fn primal() -> GraphSolveConfig {
    let mut result = config();
    result.absolute_flow_tolerance = VolumetricFlowRate::new(1e-12);
    result.relative_flow_tolerance = 1e-12;
    result
}
fn scalar() -> GraphFanConfig {
    GraphFanConfig { max_iterations: 128, absolute_flow_tolerance: VolumetricFlowRate::new(1e-10),
        absolute_pressure_tolerance: Pressure::new(1e-10) }
}
fn sensitivity() -> GraphSensitivityConfig {
    GraphSensitivityConfig { minimum_absolute_flow: VolumetricFlowRate::new(1e-12),
        linear: CheckedCgConfig { tolerance: 1e-12, max_iterations: 128, max_corrections: 2 } }
}
fn bank(count: usize, arrangement: FanArrangement, speed: f64, samples: &[(f64, f64)]) -> FanBank {
    FanBank::new(FanCurve::new("sensitivity fixture", samples.iter().map(|&(q, p)|
        FanPoint::new(VolumetricFlowRate::new(q), Pressure::new(p))).collect(),
        SourceProvenance::new("analytic fixture", "fan-sensitivity-v1"), 0.01,
        ToleranceBasis::EngineeringAllowance, VolumetricFlowRate::new(0.0), (0.25, 2.0),
    ).unwrap(), count, arrangement, speed).unwrap()
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
    graph: &LossGraph, fan: &FanBank, z: &[f64], ds: f64, dc: &[f64], h: f64, cx: &Cx<'_>,
) -> GraphSolution {
    let mut branches = graph.branches().to_vec();
    for (branch, &direction) in branches.iter_mut().zip(z) {
        branch.loss.resistance = LossResistance::new(branch.loss.resistance.value() * (1.0 + h * direction));
    }
    let changed = LossGraph::new(graph.node_count(), branches).unwrap();
    let samples: Vec<_> = fan.curve().points().iter().zip(dc).map(|(p, &d)|
        (p.flow.value(), p.pressure.value() + h * d)).collect();
    let fan = bank(fan.count(), fan.arrangement(), fan.speed_ratio() * (1.0 + h * ds), &samples);
    changed.solve_with_fan(0, graph.node_count() - 1, &fan, primal(), scalar(), cx).unwrap().network
}

#[test]
fn analytic_single_loss_includes_feedback_and_fan_law_scaling() {
    let graph = LossGraph::new(2, vec![edge("duct", 0, 1, 4.0)]).unwrap();
    for (count, arrangement, speed) in [
        (1, FanArrangement::Series, 1.0), (2, FanArrangement::Series, 1.0),
        (2, FanArrangement::Parallel, 1.0), (1, FanArrangement::Series, 0.7),
    ] {
        let fan = bank(count, arrangement, speed, &[(0.0, 100.0), (10.0, 10.0)]);
        with_cx(|cx| {
            let linear = graph.solve_with_fan_linearized(0, 1, &fan, primal(), scalar(), sensitivity(), cx).unwrap();
            let q = linear.solution().node_outflows[0].value();
            let slope = 9.0 * fan.pressure_factor() / fan.flow_factor();
            let gradient = linear.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx).unwrap();
            close(gradient.log_resistances[0], -4.0 * q * q / (8.0 * q + slope), 1e-10);
            close(gradient.log_speed_ratio, q, 1e-10);
            close(gradient.curve_pressures.iter().map(|w| w.value()).sum::<f64>(),
                fan.pressure_factor() / (8.0 * q + slope), 1e-10);
            let tangent = linear.tangent(&[0.0], 1.0, &[Pressure::new(0.0); 2], cx).unwrap();
            close(tangent.fan_flow.value(), q, 1e-10);
            close(tangent.fan_pressure.value(), 8.0 * q * q, 1e-10);
            assert!(gradient.objective_linear_report.converged_euclidean());
            assert!(gradient.inlet_linear_report.converged_euclidean());
            assert_eq!(linear.fan(), &fan);
        });
    }
}

#[test]
fn bridge_total_gradients_match_every_control_and_forward_reverse_products() {
    let graph = bridge();
    let fan = bank(1, FanArrangement::Series, 1.0, &[(0.0, 100.0), (5.0, 60.0), (10.0, 10.0)]);
    let wp = [0.2, -0.3, 0.1, 0.4].map(PressureWeight::new);
    let wq = [1.0, -2.0, 0.7, 0.25, 3.0].map(FlowWeight::new);
    with_cx(|cx| {
        let linear = graph.solve_with_fan_linearized(0, 3, &fan, primal(), scalar(), sensitivity(), cx).unwrap();
        assert!(linear.solution().branches[2].flow.value() < 0.0);
        let gradient = linear.adjoint(&wp, &wq, cx).unwrap();
        let z = [0.2, -0.4, 0.1, 0.5, -0.3];
        let dc = [0.6, -0.2, 0.1];
        let ds = 0.15;
        let tangent = linear.tangent(&z, ds, &dc.map(Pressure::new), cx).unwrap();
        let forward = tangent.network.pressures.iter().zip(wp).map(|(p, w)| p.value() * w.value()).sum::<f64>()
            + tangent.network.flows.iter().zip(wq).map(|(q, w)| q.value() * w.value()).sum::<f64>();
        let reverse = gradient.log_resistances.iter().zip(z).map(|(g, d)| g * d).sum::<f64>()
            + gradient.log_speed_ratio * ds
            + gradient.curve_pressures.iter().zip(dc).map(|(g, d)| g.value() * d).sum::<f64>();
        close(forward, reverse, 1e-10);
        assert!(tangent.network.max_node_imbalance.value() < 1e-10);
        let h = 1e-4;
        let plus = perturbed(&graph, &fan, &z, ds, &dc, h, cx);
        let minus = perturbed(&graph, &fan, &z, ds, &dc, -h, cx);
        close(forward, (objective(&plus, &wp, &wq) - objective(&minus, &wp, &wq)) / (2.0 * h), 2e-7);
        for control in 0..9 {
            let mut z = [0.0; 5];
            let mut dc = [0.0; 3];
            let mut ds = 0.0;
            let expected = if control < 5 {
                z[control] = 1.0; gradient.log_resistances[control]
            } else if control == 5 {
                ds = 1.0; gradient.log_speed_ratio
            } else {
                dc[control - 6] = 1.0; gradient.curve_pressures[control - 6].value()
            };
            let plus = perturbed(&graph, &fan, &z, ds, &dc, h, cx);
            let minus = perturbed(&graph, &fan, &z, ds, &dc, -h, cx);
            close(expected, (objective(&plus, &wp, &wq) - objective(&minus, &wp, &wq)) / (2.0 * h), 2e-7);
        }
        assert_eq!(gradient.curve_pressures.iter().filter(|w| w.value() == 0.0).count(), 1);
    });
}

#[test]
fn knots_and_retained_domain_endpoints_refuse_two_sided_derivatives() {
    with_cx(|cx| {
        for (resistance, speed, samples) in [
            (1.0, 1.0, vec![(0.0, 100.0), (5.0, 25.0), (10.0, 10.0)]),
            (1.0, 1.0, vec![(0.0, 100.0), (10.0, 100.0)]),
            (4.0, 2.0, vec![(0.0, 100.0), (10.0, 10.0)]),
        ] {
            let graph = LossGraph::new(2, vec![edge("duct", 0, 1, resistance)]).unwrap();
            let fan = bank(1, FanArrangement::Series, speed, &samples);
            assert!(graph.solve_with_fan(0, 1, &fan, primal(), scalar(), cx).is_ok());
            assert!(matches!(graph.solve_with_fan_linearized(0, 1, &fan, primal(), scalar(), sensitivity(), cx),
                Err(GraphFanSensitivityError::NonDifferentiableFan(_))));
        }
    });
}

#[test]
fn flat_fan_curve_has_no_resistance_pressure_feedback() {
    let graph = LossGraph::new(2, vec![edge("duct", 0, 1, 4.0)]).unwrap();
    let fan = bank(1, FanArrangement::Series, 1.0, &[(0.0, 100.0), (10.0, 100.0)]);
    with_cx(|cx| {
        let linear = graph.solve_with_fan_linearized(0, 1, &fan, primal(), scalar(), sensitivity(), cx).unwrap();
        let gradient = linear.adjoint(&[PressureWeight::new(1.0), PressureWeight::new(0.0)],
            &[FlowWeight::new(0.0)], cx).unwrap();
        assert_eq!(gradient.log_resistances, vec![0.0]);
        close(gradient.log_speed_ratio, 200.0, 1e-10);
        let tangent = linear.tangent(&[1.0], 0.0, &[Pressure::new(0.0); 2], cx).unwrap();
        assert_eq!(tangent.fan_pressure.value(), 0.0);
        close(tangent.fan_flow.value(), -2.5, 1e-10);
    });
}

#[test]
fn malformed_fan_controls_and_cancelled_derivatives_refuse() {
    let graph = LossGraph::new(2, vec![edge("duct", 0, 1, 4.0)]).unwrap();
    let fan = bank(1, FanArrangement::Series, 1.0, &[(0.0, 100.0), (10.0, 10.0)]);
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 17, kernel_id: 712, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let linear = graph.solve_with_fan_linearized(0, 1, &fan, primal(), scalar(), sensitivity(), &cx).unwrap();
        assert!(linear.tangent(&[0.0], 0.0, &[], &cx).is_err());
        assert!(linear.tangent(&[0.0], f64::NAN, &[Pressure::new(0.0); 2], &cx).is_err());
        assert!(linear.tangent(&[0.0], 0.0, &[Pressure::new(f64::INFINITY); 2], &cx).is_err());
        gate.request();
        assert!(matches!(linear.tangent(&[0.0], 1.0, &[Pressure::new(0.0); 2], &cx),
            Err(GraphFanSensitivityError::Hydraulic(GraphSensitivityError::Graph(GraphError::Interrupted)))));
        assert!(matches!(linear.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], &cx),
            Err(GraphFanSensitivityError::Hydraulic(GraphSensitivityError::Graph(GraphError::Interrupted)))));
    });
}

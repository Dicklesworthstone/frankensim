//! Analytic G0/G1/G3 checks of Darcy drag in the production graph owner.

mod forchheimer;

use super::*;
use crate::graph::{GraphFanConfig, LossGraph};
use crate::graph::sensitivity::{FlowWeight, GraphSensitivityConfig, PressureWeight};
use crate::graph::tests::{boundary, close, config, edge, with_cx};
use crate::{FanArrangement, FanBank, FanCurve, FanPoint};
use fs_material::gas::GasSpec;
use fs_material::porous::permeability::{
    IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use fs_qty::{Pressure, VolumetricFlowRate};
use fs_solver::CheckedCgConfig;

fn source() -> SourceProvenance {
    SourceProvenance::new("analytic Darcy fixture", "darcy-fixture-v1")
}

fn linear(name: &str, resistance: f64) -> LinearLossElement {
    LinearLossElement::new(name, LinearLossResistance::new(resistance), 0.0,
        source(), ToleranceBasis::Analytic).unwrap()
}

fn passive(name: &str, from: usize, to: usize, a: f64, b: f64) -> PassiveBranch {
    if a == 0.0 { return edge(name, from, to, b).into(); }
    let loss = linear(name, a);
    if b == 0.0 { PassiveBranch::Linear { from, to, loss } }
    else { PassiveBranch::LinearQuadratic { branch: edge(name, from, to, b), linear_loss: loss } }
}

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 100_000.0).unwrap()
}

fn barrier(permeability: f64, slip: f64, max_speed: f64) -> DarcyGasBarrier {
    DarcyGasBarrier::try_new(
        vec![PorousGasLayer::try_new(IntrinsicPermeability::isotropic(permeability).unwrap(),
            0.01, slip, "explicit Darcy fixture").unwrap()],
        [1.0, 0.0, 0.0],
        PorousGasSupport { temperature_k: (250.0, 800.0), pressure_pa: (90_000.0, 110_000.0),
            maximum_superficial_speed_m_s: max_speed },
    ).unwrap()
}

fn material_loss(permeability: f64, temperature: f64, max_speed: f64) -> LinearLossElement {
    LinearLossElement::from_darcy("filter", barrier(permeability, 0.0, max_speed),
        gas(temperature), 0.02, 0.01, 0.0, source(), ToleranceBasis::Analytic).unwrap()
}

fn fan() -> FanBank {
    FanBank::new(FanCurve::new(
        "linear analytic fan",
        vec![FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(100.0)),
            FanPoint::new(VolumetricFlowRate::new(10.0), Pressure::new(0.0))],
        source(), 0.0, ToleranceBasis::Analytic, VolumetricFlowRate::new(0.0), (0.5, 2.0),
    ).unwrap(), 1, FanArrangement::Series, 1.0).unwrap()
}

fn fan_config() -> GraphFanConfig {
    GraphFanConfig { max_iterations: 128,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-9),
        absolute_pressure_tolerance: Pressure::new(1e-8) }
}

fn sensitivity() -> GraphSensitivityConfig {
    GraphSensitivityConfig { minimum_absolute_flow: VolumetricFlowRate::new(0.0),
        linear: CheckedCgConfig { tolerance: 1e-12, max_iterations: 128, max_corrections: 2 } }
}

#[test]
fn linear_flow_reverses_exactly_and_mixed_loss_matches_analytic_root() {
    for (a, b, pressure, expected) in [(10.0, 0.0, 30.0, 3.0), (2.0, 3.0, 16.0, 2.0)] {
        let graph = LossGraph::new_mixed(2, vec![passive("filter", 0, 1, a, b)]).unwrap();
        for sign in [-1.0, 0.0, 1.0] {
            let result = with_cx(|cx| graph.solve(
                &[boundary(0, sign * pressure), boundary(1, 0.0)], config(), cx)).unwrap();
            close(result.branches[0].flow.value(), sign * expected);
            assert_eq!(result.branches[0].linear_loss, graph.linear_losses()[0]);
            if sign == 0.0 { assert_eq!(result.branches[0].flow.value(), 0.0); }
        }
    }
}

#[test]
fn mixed_inversion_preserves_linear_and_quadratic_extreme_limits() {
    for (a, b, pressure, expected) in [
        (1e200, 1e-200, 1e200, 1.0),
        (1e-200, 1e200, 1e200, 1.0),
        (1e300, 1e-300, 1e-5, 1e-305),
    ] {
        let graph = LossGraph::new_mixed(2, vec![passive("scale", 0, 1, a, b)]).unwrap();
        let result = with_cx(|cx| graph.solve(
            &[boundary(0, pressure), boundary(1, 0.0)], config(), cx)).unwrap();
        close(result.branches[0].flow.value() / expected, 1.0);
    }
}

#[test]
fn mixed_series_and_parallel_networks_conserve_volume() {
    let series = LossGraph::new_mixed(3, vec![
        passive("darcy", 0, 1, 2.0, 0.0), passive("orifice", 1, 2, 0.0, 3.0),
    ]).unwrap();
    let result = with_cx(|cx| series.solve(
        &[boundary(0, 16.0), boundary(2, 0.0)], config(), cx)).unwrap();
    close(result.pressures[1].value(), 12.0);
    for branch in &result.branches { close(branch.flow.value(), 2.0); }
    assert!(result.max_node_imbalance.value() < 1e-8);
    let parallel = LossGraph::new_mixed(2, vec![
        passive("darcy", 0, 1, 2.0, 0.0), passive("reverse-orifice", 1, 0, 0.0, 4.0),
    ]).unwrap();
    let result = with_cx(|cx| parallel.solve(
        &[boundary(0, 16.0), boundary(1, 0.0)], config(), cx)).unwrap();
    close(result.branches[0].flow.value(), 8.0);
    close(result.branches[1].flow.value(), -2.0);
    close(result.node_outflows[0].value(), 10.0);
}

#[test]
fn mixed_loop_has_manufactured_pressure_and_reverse_flow_solution() {
    let graph = LossGraph::new_mixed(4, vec![
        passive("feed", 0, 1, 1.0, 2.0 / 3.0),
        passive("left", 1, 3, 2.0, 3.0),
        passive("cross", 2, 1, 4.0, 3.0),
        passive("right-feed", 0, 2, 16.0, 0.0),
        passive("right", 2, 3, 0.0, 2.25),
    ]).unwrap();
    let result = with_cx(|cx| graph.solve(
        &[boundary(0, 25.0), boundary(3, 0.0)], config(), cx)).unwrap();
    for (actual, expected) in result.pressures.iter().zip([25.0, 16.0, 9.0, 0.0]) {
        close(actual.value(), expected);
    }
    for (actual, expected) in result.branches.iter().zip([3.0, 2.0, -1.0, 1.0, 2.0]) {
        close(actual.flow.value(), expected);
        let a = actual.linear_loss.as_ref().map_or(0.0, |loss| loss.resistance.value());
        let q = actual.flow.value();
        close(actual.pressure_drop.value(), a * q + actual.loss.resistance.value() * q * q.abs());
    }
    close(result.boundary_imbalance.value(), 0.0);
}

#[test]
fn material_permeability_and_actual_viscosity_change_solved_flow() {
    for (permeability, temperature) in [(1e-9, 300.0), (2e-9, 300.0), (1e-9, 600.0)] {
        let loss = material_loss(permeability, temperature, 5.0);
        let expected = permeability * 0.02 * 100.0 / (gas(temperature).dynamic_viscosity * 0.01);
        let graph = LossGraph::new_mixed(2, vec![PassiveBranch::Linear {
            from: 0, to: 1, loss: loss.clone(),
        }]).unwrap();
        let result = with_cx(|cx| graph.solve(
            &[boundary(0, 100.0), boundary(1, 0.0)], config(), cx)).unwrap();
        close(result.branches[0].flow.value(), expected);
        assert_eq!(result.branches[0].linear_loss.as_ref().unwrap().material(), loss.material());
    }
}

#[test]
fn material_velocity_pressure_and_sealed_paths_refuse_explicitly() {
    for (loss, drop) in [(material_loss(1e-9, 300.0, 0.1), 100.0),
        (material_loss(1e-9, 300.0, 100.0), 2000.0)]
    {
        let graph = LossGraph::new_mixed(2, vec![PassiveBranch::Linear { from: 0, to: 1, loss }]).unwrap();
        assert!(matches!(with_cx(|cx| graph.solve(
            &[boundary(0, drop), boundary(1, 0.0)], config(), cx)), Err(GraphError::InvalidBranch { .. })));
    }
    assert!(LinearLossElement::from_darcy("sealed", barrier(0.0, 0.0, 1.0),
        gas(300.0), 0.02, 0.01, 0.0, source(), ToleranceBasis::Analytic).is_err());
    let mut loss = material_loss(1e-9, 300.0, 5.0);
    loss.resistance = LinearLossResistance::new(2.0 * loss.resistance.value());
    assert!(LossGraph::new_mixed(2, vec![PassiveBranch::Linear { from: 0, to: 1, loss }]).is_err());
}

#[test]
fn material_pressure_reference_admits_series_and_refuses_large_gauge_offsets() {
    let loss = material_loss(1e-9, 300.0, 5.0);
    let mut second = loss.clone();
    second.name = "second filter".to_owned();
    let series = LossGraph::new_mixed(3, vec![
        PassiveBranch::Linear { from: 0, to: 1, loss: loss.clone() },
        PassiveBranch::Linear { from: 1, to: 2, loss: second },
    ]).unwrap();
    let solved = with_cx(|cx| series.solve(
        &[boundary(0, 200.0), boundary(2, 0.0)], config(), cx)).unwrap();
    close(solved.pressures[1].value(), 100.0);
    close(solved.branches[0].flow.value(), solved.branches[1].flow.value());

    let graph = LossGraph::new_mixed(2, vec![
        PassiveBranch::Linear { from: 0, to: 1, loss },
    ]).unwrap();
    // Every pair has the same supported 100 Pa drop and superficial speed.
    // The first violates absolute source support, the second violates only
    // the frozen-density allowance, and the third has negative absolute pressure.
    for (from, to) in [(200_000.0, 199_900.0), (2000.0, 1900.0), (-200_000.0, -200_100.0)] {
        assert!(matches!(with_cx(|cx| graph.solve(
            &[boundary(0, from), boundary(1, to)], config(), cx)),
            Err(GraphError::InvalidBranch { branch: 0, .. })));
    }
}

#[test]
fn material_edges_cannot_mix_different_frozen_gas_states() {
    let mut other_spec = GasSpec::dry_air_ussa1976();
    // Explicit synthetic helium-like source: the exact transport coefficient
    // is irrelevant to this refusal of a different coherent gas identity.
    other_spec.molar_mass = 0.004;
    other_spec.gamma = 5.0 / 3.0;
    other_spec.conductivity = fs_material::gas::ConductivityModel::Eucken;
    let other = GasState::try_new(&other_spec, 300.0, 100_000.0).unwrap();
    for other in [gas(600.0), other] {
        let second = LinearLossElement::from_darcy("second", barrier(1e-9, 0.0, 5.0),
            other, 0.02, 0.01, 0.0, source(), ToleranceBasis::Analytic).unwrap();
        assert!(matches!(LossGraph::new_mixed(3, vec![
            PassiveBranch::Linear { from: 0, to: 1, loss: material_loss(1e-9, 300.0, 5.0) },
            PassiveBranch::Linear { from: 1, to: 2, loss: second },
        ]), Err(GraphError::InvalidBranch { branch: 1, .. })));
    }
}

#[test]
fn darcy_and_mixed_fan_intersections_use_the_actual_pressure_law() {
    for (a, b, expected) in [(20.0, 0.0, 10.0 / 3.0),
        (2.0, 3.0, 200.0 / (12.0 + 1344.0_f64.sqrt()))]
    {
        let graph = LossGraph::new_mixed(2, vec![passive("filter", 0, 1, a, b)]).unwrap();
        let point = with_cx(|cx| graph.solve_with_fan(0, 1, &fan(), config(), fan_config(), cx)).unwrap();
        close(point.flow.value(), expected);
        close(point.network_pressure.value(), a * expected + b * expected * expected);
        assert!(point.pressure_residual.value().abs() <= fan_config().absolute_pressure_tolerance.value());
    }
}

#[test]
fn fan_trials_do_not_falsely_refuse_a_supported_material_operating_point() {
    let gas = gas(300.0);
    let permeability = gas.dynamic_viscosity * 0.01 / 20.0;
    let loss = LinearLossElement::from_darcy("filter", barrier(permeability, 0.0, 4.0),
        gas, 1.0, 0.01, 0.0, source(), ToleranceBasis::Analytic).unwrap();
    let graph = LossGraph::new_mixed(2, vec![PassiveBranch::Linear { from: 0, to: 1, loss }]).unwrap();
    // Peak-pressure trial has Q/A = 5 m/s, outside the 4 m/s source support;
    // the actual fan intersection is 10/3 m/s and must remain admissible.
    let point = with_cx(|cx| graph.solve_with_fan(0, 1, &fan(), config(), fan_config(), cx)).unwrap();
    close(point.flow.value(), 10.0 / 3.0);
}

#[test]
fn mixed_graph_tangent_and_adjoint_match_series_closed_form() {
    let graph = LossGraph::new_mixed(3, vec![
        passive("darcy", 0, 1, 2.0, 0.0), passive("orifice", 1, 2, 0.0, 3.0),
    ]).unwrap();
    with_cx(|cx| {
        let solved = graph.solve_linearized(&[boundary(0, 16.0), boundary(2, 0.0)],
            config(), sensitivity(), cx).unwrap();
        let tangent = solved.tangent(&[1.0, 0.0], &[Pressure::new(0.0); 2], cx).unwrap();
        for flow in &tangent.flows { close(flow.value(), -2.0 / 7.0); }
        close(tangent.pressures[1].value(), -24.0 / 7.0);
        let gradient = solved.adjoint(&[PressureWeight::new(0.0); 3],
            &[FlowWeight::new(1.0), FlowWeight::new(0.0)], cx).unwrap();
        close(gradient.log_resistances[0], -2.0 / 7.0);
        close(gradient.log_resistances[1], -6.0 / 7.0);
        close(gradient.boundary_pressures[0].value(), 1.0 / 14.0);
        close(gradient.boundary_pressures[1].value(), -1.0 / 14.0);
    });
}

#[test]
fn stagnant_linear_edges_have_their_physical_finite_derivative() {
    let graph = LossGraph::new_mixed(2, vec![passive("darcy", 0, 1, 2.0, 0.0)]).unwrap();
    with_cx(|cx| {
        let solved = graph.solve_linearized(&[boundary(0, 0.0), boundary(1, 0.0)],
            config(), sensitivity(), cx).unwrap();
        let tangent = solved.tangent(&[1.0], &[Pressure::new(2.0), Pressure::new(0.0)], cx).unwrap();
        close(tangent.flows[0].value(), 1.0);
        let gradient = solved.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx).unwrap();
        assert_eq!(gradient.log_resistances[0], 0.0);
        close(gradient.boundary_pressures[0].value(), 0.5);
    });
}

#[test]
fn combined_edge_resistance_derivative_scales_both_physical_terms() {
    let graph = LossGraph::new_mixed(2, vec![passive("combined", 0, 1, 2.0, 3.0)]).unwrap();
    with_cx(|cx| {
        let boundaries = [boundary(0, 16.0), boundary(1, 0.0)];
        let solved = graph.solve_linearized(&boundaries, config(), sensitivity(), cx).unwrap();
        let tangent = solved.tangent(&[1.0], &[Pressure::new(0.0); 2], cx).unwrap();
        close(tangent.flows[0].value(), -8.0 / 7.0);
        let gradient = solved.adjoint(&[PressureWeight::new(0.0); 2], &[FlowWeight::new(1.0)], cx).unwrap();
        close(gradient.log_resistances[0], -8.0 / 7.0);
        close(gradient.boundary_pressures[0].value(), 1.0 / 14.0);
        let h = 1e-5;
        let mut flows = [0.0; 2];
        for (index, scale) in [1.0 - h, 1.0 + h].into_iter().enumerate() {
            let perturbed = LossGraph::new_mixed(2,
                vec![passive("combined", 0, 1, 2.0 * scale, 3.0 * scale)]).unwrap();
            flows[index] = perturbed.solve(&boundaries, config(), cx).unwrap().branches[0].flow.value();
        }
        close(tangent.flows[0].value(), (flows[1] - flows[0]) / (2.0 * h));
    });
}

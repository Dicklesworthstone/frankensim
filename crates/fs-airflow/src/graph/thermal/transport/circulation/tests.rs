//! G1/G3: actual buoyant hydraulics followed by simultaneous closed-loop heat
//! transport, against independent exponential-law and derivative references.

use super::*;
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use crate::graph::{GraphSolution, LossGraph};
use fs_material::fluid_transport::FluidTransportState;

fn air() -> TransportAir {
    TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 1000.0 }
}

fn budgets() -> TransportConfig {
    TransportConfig {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-10),
        relative_flow_tolerance: 1e-10,
        absolute_heat_tolerance_w: 1e-7,
        relative_heat_tolerance: 1e-9,
    }
}

fn flow(cx: &Cx<'_>, reverse_cold_declaration: bool) -> GraphSolution {
    labelled_flow(cx, reverse_cold_declaration, [0, 1])
}

fn labelled_flow(cx: &Cx<'_>, reverse_cold_declaration: bool, nodes: [usize; 2]) -> GraphSolution {
    let state = |density| FluidTransportState::try_new(density, 1000.0, 0.1, 0.001).unwrap();
    let [bottom, top] = nodes;
    let cold = if reverse_cold_declaration { edge("cold", bottom, top, 1.0) }
        else { edge("cold", top, bottom, 1.0) };
    let mut positions = [[0.0; 3]; 2];
    positions[top][2] = 1.0;
    let mut controls = hydraulic_config();
    controls.absolute_flow_tolerance = VolumetricFlowRate::new(1e-13);
    controls.relative_flow_tolerance = 1e-13;
    LossGraph::new(2, vec![edge("hot", bottom, top, 1.0), cold]).unwrap()
        .with_boussinesq_gravity(&state(1.0), &[state(0.875), state(1.0)],
            &positions, [0.0, 0.0, -8.0], 0.2, cx).unwrap()
        .solve(&[boundary(bottom, 0.0)], controls, cx).unwrap()
}

fn models(capacity: f64, ntu: [f64; 2]) -> Vec<BranchThermalModel> {
    ["hot wall", "cold wall"].into_iter().zip(ntu).map(|(name, x)| {
        BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, capacity*x).unwrap()])
    }).collect()
}

fn close(actual: f64, expected: f64) {
    assert!((actual-expected).abs() < 2e-7,
        "actual {actual:.16e}, expected {expected:.16e}");
}

fn oracle(ntu: [f64; 2], walls: [f64; 2]) -> [f64; 2] {
    let a = (-ntu[0]).exp();
    let b = (-ntu[1]).exp();
    let cold = ((1.0-b)*walls[1] + b*(1.0-a)*walls[0]) / (1.0-a*b);
    [cold, walls[0] + (cold-walls[0])*a]
}

#[test]
fn g1_closed_buoyant_loop_transfers_equal_opposite_wall_heat_without_fresh_air() {
    with_cx(|cx| {
        let flow = flow(cx, false);
        let capacity = flow.branches[0].flow.value()*1000.0;
        close(flow.branches[0].flow.value(), 0.5_f64.sqrt());
        let net = TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [0.6, 1.1]),
            Temperature::new(325.0), budgets()).unwrap();
        let result = net.march(cx, &[360.0, 300.0]).unwrap();
        let expected = oracle([0.6, 1.1], [360.0, 300.0]);
        for (actual, expected) in result.node_temperatures_k.iter().zip(expected) {
            close(actual.unwrap(), expected);
        }
        let hot = result.branches[0].march.as_ref().unwrap().total_heat_rate_w;
        let cold = result.branches[1].march.as_ref().unwrap().total_heat_rate_w;
        assert!(hot > 1000.0 && cold < -1000.0);
        close(hot, capacity*(expected[1]-expected[0]));
        close(hot, -cold);
        close(result.wall_heat_rate_w, 0.0);
        assert_eq!(result.external_heat_gain_w, 0.0);
        assert!(result.recirculation.is_none());
        assert_eq!(net.initial_references(cx).unwrap(), vec![325.0; 2]);
        assert_eq!(result, net.march(cx, &[360.0, 300.0]).unwrap());

        // Numerical coordinates/initial guesses supply no physical heat bath.
        let other = TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [0.6, 1.1]),
            Temperature::new(450.0), budgets()).unwrap().march(cx, &[360.0, 300.0]).unwrap();
        assert_eq!(result.node_temperatures_k, other.node_temperatures_k);
        assert_eq!(result.branches, other.branches);
        assert!(TransportNetwork::new(cx, &flow, air(), models(capacity, [0.6, 1.1]),
            &[], budgets()).is_err(), "ordinary DAG admission remains explicit");
    });
}

#[test]
fn g3_reversed_edge_and_scaled_conductance_preserve_the_physical_loop() {
    with_cx(|cx| {
        for (reverse, nodes) in [(false, [0, 1]), (true, [0, 1]), (false, [1, 0]), (true, [1, 0])] {
            let flow = labelled_flow(cx, reverse, nodes);
            let capacity = flow.branches[0].flow.value()*1000.0;
            assert_eq!(flow.branches[1].flow.value().is_sign_negative(), reverse);
            let net = TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [0.6, 1.1]),
                Temperature::new(325.0), budgets()).unwrap();
            let scaled = net.scaled_conductance(cx, 0.5).unwrap();
            let result = scaled.march(cx, &[360.0, 300.0]).unwrap();
            for (node, expected) in nodes.into_iter().zip(oracle([0.3, 0.55], [360.0, 300.0])) {
                close(result.node_temperatures_k[node].unwrap(), expected);
            }
            assert_eq!(scaled.initial_references(cx).unwrap(), vec![325.0; 2]);
        }
        let flow = flow(cx, false);
        let capacity = flow.branches[0].flow.value()*1000.0;
        let weak = TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [1e-10; 2]),
            Temperature::new(280.0), budgets()).unwrap();
        let equilibrium = weak.march(cx, &[340.0; 2]).unwrap();
        assert_eq!(equilibrium.node_temperatures_k, vec![Some(340.0); 2]);
        assert_eq!(equilibrium.wall_heat_rate_w, 0.0);
    });
}

#[test]
fn g1_two_buoyant_paths_mix_by_capacity_inside_the_same_circulation() {
    with_cx(|cx| {
        let state = |density| FluidTransportState::try_new(density, 1000.0, 0.1, 0.001).unwrap();
        let mut controls = hydraulic_config();
        controls.absolute_flow_tolerance = VolumetricFlowRate::new(1e-14);
        controls.relative_flow_tolerance = 1e-14;
        // Deliberately not ordered by upstream node or by heat-flow path.
        let flow = LossGraph::new(3, vec![
            edge("return", 2, 0, 1.0), edge("right rise", 0, 2, 1.0),
            edge("left rise", 0, 1, 1.0), edge("top", 1, 2, 1.0),
        ]).unwrap().with_boussinesq_gravity(&state(1.0),
            &[state(1.0), state(0.9375), state(0.875), state(1.0)],
            &[[0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]],
            [0.0, 0.0, -8.0], 0.2, cx).unwrap()
            .solve(&[boundary(0, 0.0)], controls, cx).unwrap();
        assert!(flow.branches.iter().all(|branch| branch.flow.value() > 0.0));
        let c: Vec<_> = flow.branches.iter().map(|branch| 1000.0*branch.flow.value()).collect();
        let exchange = |name, conductance| BranchThermalModel::Exchange(
            vec![AirSegment::new(name, 1.0, conductance).unwrap()]);
        let net = TransportNetwork::new_circulating(cx, &flow, air(), vec![
            exchange("cold", 1.1*c[0]), exchange("right", 0.8*c[1]),
            exchange("left", 0.6*c[2]), BranchThermalModel::Adiabatic,
        ], Temperature::new(325.0), budgets()).unwrap();
        let result = net.march(cx, &[300.0, 340.0, 360.0]).unwrap();
        let weight = c[2]/(c[2]+c[1]);
        let a = weight*(-0.6_f64).exp() + (1.0-weight)*(-0.8_f64).exp();
        let b = weight*(1.0-(-0.6_f64).exp())*360.0
            + (1.0-weight)*(1.0-(-0.8_f64).exp())*340.0;
        let return_attenuation = (-1.1_f64).exp();
        let bottom = ((1.0-return_attenuation)*300.0 + return_attenuation*b)
            /(1.0-return_attenuation*a);
        close(result.node_temperatures_k[0].unwrap(), bottom);
        close(result.node_temperatures_k[1].unwrap(), 360.0+(bottom-360.0)*(-0.6_f64).exp());
        close(result.node_temperatures_k[2].unwrap(), a*bottom+b);
        close(result.wall_heat_rate_w, 0.0);
        assert!((weight-0.5).abs() > 0.1, "arithmetic-mean mixing must be falsifiable");
    });
}

#[test]
fn g1_circulation_tangents_match_finite_differences_and_the_full_adjoint_identity() {
    with_cx(|cx| {
        let flow = flow(cx, false);
        let capacity = flow.branches[0].flow.value()*1000.0;
        let net = TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [0.6, 1.1]),
            Temperature::new(325.0), budgets()).unwrap();
        let linear = net.linearize(cx, &[360.0, 300.0]).unwrap();
        let mut direction = linear.zero_direction();
        direction.walls_k = vec![0.4, -0.2];
        direction.log_conductances = vec![0.1, -0.15];
        let tangent = linear.apply(cx, &direction).unwrap();
        let epsilon = 1e-4;
        let perturbed = |sign: f64| {
            let ntu = [0.6*(sign*epsilon*0.1).exp(), 1.1*(sign*epsilon*(-0.15)).exp()];
            TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, ntu),
                Temperature::new(325.0), budgets()).unwrap()
                .march(cx, &[360.0+sign*epsilon*0.4, 300.0-sign*epsilon*0.2]).unwrap()
        };
        let plus = perturbed(1.0);
        let minus = perturbed(-1.0);
        for i in 0..2 {
            close(tangent.node_temperatures_k[i].unwrap(),
                (plus.node_temperatures_k[i].unwrap()-minus.node_temperatures_k[i].unwrap())/(2.0*epsilon));
            close(tangent.branch_outlets_k[i].unwrap(),
                (plus.branches[i].outlet_temperature_k.unwrap()-minus.branches[i].outlet_temperature_k.unwrap())/(2.0*epsilon));
            close(tangent.reference_temperatures_k[i],
                (plus.reference_temperatures_k[i]-minus.reference_temperatures_k[i])/(2.0*epsilon));
        }
        let mut objective = linear.zero_objective();
        objective.node_temperatures = vec![0.3, -0.4];
        objective.branch_outlets = vec![0.2, 0.5];
        objective.references = vec![-0.1, 0.7];
        objective.wall_heat_rate = 0.002;
        objective.heat_imbalance = 0.001;
        let gradient = linear.pullback(cx, &objective).unwrap();
        let output_dot: f64 = (0..2).map(|i| {
            objective.node_temperatures[i]*tangent.node_temperatures_k[i].unwrap()
                + objective.branch_outlets[i]*tangent.branch_outlets_k[i].unwrap()
                + objective.references[i]*tangent.reference_temperatures_k[i]
        }).sum::<f64>() + objective.wall_heat_rate*tangent.wall_heat_rate_w
            + objective.heat_imbalance*tangent.heat_imbalance_w;
        let input_dot: f64 = (0..2).map(|i| gradient.walls[i]*direction.walls_k[i]
            + gradient.log_conductances[i]*direction.log_conductances[i]).sum();
        assert!((output_dot-input_dot).abs() < 1e-10);
        assert!(linear.pullback_branch_flows(cx, &objective).is_err());
    });
}

#[test]
fn g0_closed_circulation_refuses_missing_anchors_mass_supplies_and_invalid_seeds() {
    with_cx(|cx| {
        let flow = flow(cx, false);
        let capacity = flow.branches[0].flow.value()*1000.0;
        assert!(TransportNetwork::new_circulating(cx, &flow, air(),
            vec![BranchThermalModel::Adiabatic; 2], Temperature::new(325.0), budgets()).is_err());
        // A heat anchor in one component cannot determine another loop's
        // temperature. Both sets of flows come from the original solver.
        let state = |density| FluidTransportState::try_new(density, 1000.0, 0.1, 0.001).unwrap();
        let disconnected = LossGraph::new(4, vec![
            edge("hot", 0, 1, 1.0), edge("cold", 1, 0, 1.0),
            edge("unanchored rise", 2, 3, 1.0), edge("unanchored return", 3, 2, 1.0),
        ]).unwrap().with_boussinesq_gravity(&state(1.0),
            &[state(0.875), state(1.0), state(0.875), state(1.0)],
            &[[0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
            [0.0, 0.0, -8.0], 0.2, cx).unwrap()
            .solve(&[boundary(0, 0.0), boundary(2, 0.0)], hydraulic_config(), cx).unwrap();
        let mut disconnected_models = models(capacity, [0.6, 1.1]);
        disconnected_models.extend([BranchThermalModel::Adiabatic, BranchThermalModel::Adiabatic]);
        assert!(matches!(TransportNetwork::new_circulating(cx, &disconnected, air(),
            disconnected_models, Temperature::new(325.0), budgets()),
            Err(TransportError::Node { node: 2 | 3, .. })));
        let mut leaking = flow.clone();
        leaking.branches[0].flow = VolumetricFlowRate::new(flow.branches[0].flow.value()*1.01);
        assert!(matches!(TransportNetwork::new_circulating(cx, &leaking, air(),
            models(capacity, [0.6, 1.1]), Temperature::new(325.0), budgets()),
            Err(TransportError::FlowImbalance { .. })));
        for seed in [0.0, -1.0, f64::NAN] {
            assert!(TransportNetwork::new_circulating(cx, &flow, air(), models(capacity, [0.6, 1.1]),
                Temperature::new(seed), budgets()).is_err());
        }
    });
}

#[test]
fn g4_circulation_and_its_derivatives_keep_bounded_cancellation() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 1, kernel_id: 1712, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let flow = flow(&cx, false);
        let net = TransportNetwork::new_circulating(&cx, &flow, air(),
            models(flow.branches[0].flow.value()*1000.0, [0.6, 1.1]),
            Temperature::new(325.0), budgets()).unwrap();
        let linear = net.linearize(&cx, &[360.0, 300.0]).unwrap();
        gate.request();
        assert!(matches!(net.march(&cx, &[360.0, 300.0]), Err(TransportError::Interrupted)));
        assert!(matches!(linear.apply(&cx, &linear.zero_direction()), Err(TransportError::Interrupted)));
        assert!(matches!(linear.pullback(&cx, &linear.zero_objective()), Err(TransportError::Interrupted)));
    });
}

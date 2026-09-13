//! Closed-form fixed points over real hydraulic topologies, including a
//! non-series-parallel bridge and downstream air warmed by upstream solids.

use super::*;
use super::super::transport::{BranchThermalModel, TransportAir, TransportConfig, TransportInlet};
use crate::conjugate::AirSegment;
use crate::graph::tests::{boundary, bridge, config as hydraulic_config, edge, with_cx};
use crate::graph::{GraphSolution, LossGraph};
use fs_qty::{Density, Temperature, VolumetricFlowRate};

fn flow() -> GraphSolution {
    with_cx(|cx| LossGraph::new(3, vec![edge("left", 0, 1, 1.0 / 9.0),
        edge("right", 0, 1, 1.0), edge("tail", 1, 2, 1.0 / 16.0)]).unwrap()
        .solve(&[boundary(0, 2.0), boundary(2, 0.0)], hydraulic_config(), cx).unwrap())
}
fn network<'a>(cx: &Cx<'_>, flow: &'a GraphSolution, conductances: &[f64]) -> TransportNetwork<'a> {
    TransportNetwork::new(cx, flow,
        TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 1.0 },
        conductances.iter().enumerate().map(|(i, &g)| BranchThermalModel::Exchange(
            vec![AirSegment::new(&format!("wall-{i}"), 1.0, g).unwrap()])).collect(),
        &[TransportInlet { node: 0, temperature: Temperature::new(300.0) }],
        TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(1e-9), relative_flow_tolerance: 1e-9,
            absolute_heat_tolerance_w: 1e-7, relative_heat_tolerance: 1e-8 }).unwrap()
}
fn states(references: &[f64], conductances: &[f64], walls: &[f64]) -> Vec<SolidRegionState> {
    references.iter().zip(conductances).zip(walls).enumerate().map(|(i, ((&reference, &g), &wall))| SolidRegionState {
        region: format!("wall-{i}"), area_m2: 1.0, mean_wall_temperature_k: wall,
        heat_rate_w: g * (wall - reference), mean_reference_temperature_k: Some(reference),
    }).collect()
}
fn common_solid(references: &[f64], conductances: &[f64], power: f64) -> Vec<SolidRegionState> {
    let wall = (power + references.iter().zip(conductances).map(|(r, g)| r * g).sum::<f64>())
        / conductances.iter().sum::<f64>();
    states(references, conductances, &vec![wall; references.len()])
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-7, "actual {actual:.16e}, expected {expected:.16e}");
}

#[test]
fn common_solid_split_merge_fixed_point_matches_total_effectiveness() {
    let flow = flow();
    with_cx(|cx| {
        let g = [3.0, 2.0, 4.0];
        let net = network(cx, &flow, &g);
        let mut calls = 0;
        let result = solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, references| {
            calls += 1;
            Ok(common_solid(references, &g, 40.0))
        }).unwrap();
        let effective = 1.0 - (-1.0_f64).exp() * (3.0 * (-1.0_f64).exp() + (-2.0_f64).exp()) / 4.0;
        let wall = 300.0 + 40.0 / (4.0 * effective);
        close(result.solid[0].mean_wall_temperature_k, wall);
        close(result.transport.node_temperatures_k[2].unwrap(), 310.0);
        close(result.transport.external_heat_gain_w, 40.0);
        assert_eq!(result.iterations, calls);
        for (solid, &used) in result.solid.iter().zip(&result.reference_temperatures_k) {
            assert_eq!(solid.mean_reference_temperature_k.unwrap().to_bits(), used.to_bits());
        }
        for (&used, &next) in result.reference_temperatures_k.iter().zip(&result.transport.reference_temperatures_k) {
            assert!((used - next).abs() <= ConjugateConfig::default().temperature_tolerance_k);
        }
    });
}

#[test]
fn upstream_heat_changes_the_downstream_solid_solution() {
    let flow = flow();
    with_cx(|cx| {
        let g = [3.0, 2.0, 4.0];
        let powers = [10.0, 20.0, 30.0];
        let net = network(cx, &flow, &g);
        let result = solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, references| {
            let walls: Vec<f64> = references.iter().zip(g).zip(powers).map(|((&r, g), p)| r + p / g).collect();
            Ok(states(references, &g, &walls))
        }).unwrap();
        let mixed = 300.0 + (powers[0] + powers[1]) / 4.0;
        let tail_wall = mixed + powers[2] / (4.0 * (1.0 - (-1.0_f64).exp()));
        close(result.transport.node_temperatures_k[1].unwrap(), mixed);
        close(result.transport.node_temperatures_k[2].unwrap(), 315.0);
        close(result.solid[2].mean_wall_temperature_k, tail_wall);
        close(result.transport.wall_heat_rate_w, 60.0);
        assert!(result.solid[2].mean_wall_temperature_k > 319.0);
    });
}

#[test]
fn non_series_parallel_bridge_couples_the_reverse_cross_branch() {
    with_cx(|cx| {
        let flow = bridge().solve(&[boundary(0, 25.0), boundary(3, 0.0)], hydraulic_config(), cx).unwrap();
        let g = [1.0; 5];
        let net = network(cx, &flow, &g);
        let result = solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, references| Ok(common_solid(references, &g, 20.0))).unwrap();
        let node1 = 1.0 - (-1.0_f64 / 3.0).exp();
        let cross_out = 1.0 - (1.0 - node1) * (-1.0_f64).exp();
        let direct_out = 1.0 - (-1.0_f64).exp();
        let node2 = (cross_out + direct_out) / 2.0;
        let left_out = 1.0 - (1.0 - node1) * (-0.5_f64).exp();
        let right_out = 1.0 - (1.0 - node2) * (-0.5_f64).exp();
        let effective = (left_out + right_out) / 2.0;
        close(result.solid[0].mean_wall_temperature_k, 300.0 + 20.0 / (4.0 * effective));
        close(result.transport.node_temperatures_k[3].unwrap(), 305.0);
        assert!(flow.branches[2].flow.value() < 0.0);
        close(result.transport.branches[2].inlet_temperature_k.unwrap(), result.transport.node_temperatures_k[1].unwrap());
    });
}

#[test]
fn fixed_relaxation_cancellation_resumes_the_exact_numerical_tail() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let flow = flow();
    let g = [3.0, 2.0, 4.0];
    let config = ConjugateConfig::default();
    let full = with_cx(|cx| solve_coupled_transport(cx, &network(cx, &flow, &g), &config,
        |_, references| Ok(common_solid(references, &g, 40.0))).unwrap());
    let gate = CancelGate::new_clock_free();
    let cancelled = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 1, kernel_id: 713, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        let net = network(&cx, &flow, &g);
        let mut calls = 0;
        solve_coupled_transport(&cx, &net, &config, |_, references| {
            calls += 1;
            if calls == 6 { gate.request(); }
            Ok(common_solid(references, &g, 40.0))
        }).unwrap_err()
    });
    let TransportError::Airflow(AirflowError::Cancelled { iteration, references_k }) = cancelled else { panic!("complete cancellation state required") };
    assert_eq!(iteration, 5);
    assert_eq!(references_k.len(), 3);
    let resumed = with_cx(|cx| solve_coupled_transport_from(cx, &network(cx, &flow, &g), &config, &references_k,
        |_, references| Ok(common_solid(references, &g, 40.0))).unwrap());
    assert_eq!(full.reference_temperatures_k, resumed.reference_temperatures_k);
    assert_eq!(full.solid, resumed.solid);
    assert_eq!(full.transport, resumed.transport);
    assert_eq!(full.iterations, iteration + resumed.iterations);
}

#[test]
fn budget_exhaustion_does_not_publish_a_partially_coupled_field() {
    let flow = flow();
    with_cx(|cx| {
        let config = ConjugateConfig { max_iterations: 1, ..ConjugateConfig::default() };
        let result = solve_coupled_transport(cx, &network(cx, &flow, &[3.0, 2.0, 4.0]), &config,
            |_, references| Ok(common_solid(references, &[3.0, 2.0, 4.0], 40.0)));
        assert!(matches!(result, Err(TransportError::Airflow(AirflowError::ConjugateNotConverged { iterations: 1, .. }))));
    });
}

#[test]
fn malformed_resume_and_configuration_refuse_before_solid_work() {
    let flow = flow();
    with_cx(|cx| {
        let net = network(cx, &flow, &[3.0, 2.0, 4.0]);
        for references in [vec![], vec![300.0, 300.0], vec![300.0, f64::NAN, 300.0], vec![0.0; 3]] {
            assert!(solve_coupled_transport_from(cx, &net, &ConjugateConfig::default(), &references,
                |_, _| panic!("invalid resume must not execute solid work")).is_err());
        }
        let bad = ConjugateConfig { max_iterations: 0, ..ConjugateConfig::default() };
        assert!(solve_coupled_transport(cx, &net, &bad, |_, _| panic!("invalid budget must not execute solid work")).is_err());
        assert!(matches!(solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, _| Ok(vec![])),
            Err(TransportError::Airflow(AirflowError::SolidResponseArity { .. }))));
    });
}

#[test]
fn a_wrong_solid_heat_rate_cannot_pass_a_temperature_only_convergence_test() {
    let flow = flow();
    with_cx(|cx| {
        let g = [3.0, 2.0, 4.0];
        let net = network(cx, &flow, &g);
        let result = solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, references| {
            let mut response = states(references, &g, &[350.0; 3]);
            response[0].heat_rate_w += 1.0;
            Ok(response)
        });
        assert!(matches!(result, Err(TransportError::Airflow(AirflowError::ConjugateBalanceUnclosed { .. }))));
    });
}

#[test]
fn applied_reference_and_region_ownership_are_checked_by_the_existing_driver() {
    let flow = flow();
    with_cx(|cx| {
        let g = [3.0, 2.0, 4.0];
        let net = network(cx, &flow, &g);
        for wrong_name in [false, true] {
            let result = solve_coupled_transport(cx, &net, &ConjugateConfig::default(), |_, references| {
                let mut response = common_solid(references, &g, 40.0);
                if wrong_name { response[0].region = "foreign".into(); }
                else { response[0].mean_reference_temperature_k = Some(references[0] + 1.0); }
                Ok(response)
            });
            match result {
                Err(TransportError::Airflow(AirflowError::SegmentRegionMismatch { .. })) if wrong_name => {}
                Err(TransportError::Airflow(AirflowError::ReferenceTemperatureMismatch { .. })) if !wrong_name => {}
                other => panic!("wrong solid binding must refuse: {other:?}"),
            }
        }
    });
}

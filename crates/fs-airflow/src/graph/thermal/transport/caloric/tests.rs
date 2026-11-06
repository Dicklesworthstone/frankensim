//! Source enthalpy through real hydraulic splits, merges and solid coupling.

use super::*;
use crate::graph::LossGraph;
use crate::graph::tests::{boundary, config as hydraulic_config, edge, with_cx};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityQuery;
use fs_qty::{Pressure, QuantitySpec};

fn source(intercept: f64, slope: f64, datum: f64) -> Arc<IsobaricHeatCapacityLaw> {
    let key = PropertyKey::new("Cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [250.0, 400.0, 700.0]
                    .into_iter()
                    .map(|t| (t, intercept + slope * t))
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    700.0,
                )
                .with_quantity(
                    "p",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "analytic network Cp fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured single coolant".into(),
            phase: "single-phase".into(),
            process: "fixed pressure".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
            .unwrap()
    };
    Arc::new(
        IsobaricHeatCapacityLaw::try_new(
            &card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &key,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &point(250.0),
                upper: &point(700.0),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            datum,
        )
        .unwrap(),
    )
}

fn budgets() -> TransportConfig {
    TransportConfig {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-9),
        relative_flow_tolerance: 1e-9,
        absolute_heat_tolerance_w: 1e-6,
        relative_heat_tolerance: 1e-10,
    }
}
fn inlet(node: usize, temperature: f64) -> TransportInlet {
    TransportInlet {
        node,
        temperature: Temperature::new(temperature),
    }
}
fn heated(name: &str, conductance: f64) -> BranchThermalModel {
    BranchThermalModel::Exchange(vec![AirSegment::new(name, 1.0, conductance).unwrap()])
}
fn network<'a>(
    cx: &Cx<'_>,
    flow: &'a GraphSolution,
    law: Arc<IsobaricHeatCapacityLaw>,
    models: Vec<BranchThermalModel>,
    inlets: &[TransportInlet],
) -> CaloricTransportNetwork<'a> {
    CaloricTransportNetwork::new(
        cx,
        flow,
        Density::new(1.0),
        Temperature::new(300.0),
        law,
        models,
        inlets,
        budgets(),
    )
    .unwrap()
}
fn merge_flow(cx: &Cx<'_>) -> GraphSolution {
    LossGraph::new(
        3,
        vec![
            edge("left", 0, 1, 1.0 / 9.0),
            edge("right", 0, 1, 1.0),
            edge("tail", 1, 2, 1.0 / 16.0),
        ],
    )
    .unwrap()
    .solve(
        &[boundary(0, 2.0), boundary(2, 0.0)],
        hydraulic_config(),
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}

// Independent calorimetry for Cp(T)=800+2T; all oracles use this primitive,
// not the source evaluator or the native network's enthalpy coordinates.
fn delta_h(initial: f64, final_t: f64) -> f64 {
    let delta = final_t - initial;
    (800.0 + 2.0 * initial) * delta + delta * delta
}
fn after_heat(initial: f64, heat: f64) -> f64 {
    let cp = 800.0 + 2.0 * initial;
    initial + 2.0 * heat / (cp + (cp * cp + 4.0 * heat).sqrt())
}

#[test]
fn supplies_mix_original_enthalpy_instead_of_temperature_or_endpoint_cp() {
    with_cx(|cx| {
        let flow = LossGraph::new(3, vec![edge("cold", 0, 2, 1.0), edge("hot", 1, 2, 1.0)])
            .unwrap()
            .solve(
                &[boundary(0, 9.0), boundary(1, 4.0), boundary(2, 0.0)],
                hydraulic_config(),
                cx,
            )
            .unwrap();
        let law = source(800.0, 2.0, 300.0);
        let net = network(
            cx,
            &flow,
            law.clone(),
            vec![BranchThermalModel::Adiabatic; 2],
            &[inlet(0, 300.0), inlet(1, 600.0)],
        );
        let result = net.march(cx, &[]).unwrap();
        let cold_mass = flow.branches[0].flow.value();
        let hot_mass = flow.branches[1].flow.value();
        let expected = after_heat(
            300.0,
            hot_mass / (cold_mass + hot_mass) * delta_h(300.0, 600.0),
        );
        close(result.node_temperatures_k[2].unwrap(), expected, 2e-10);
        close(result.external_heat_gain_w, 0.0, 2e-8);
        assert!(
            (expected - (cold_mass * 300.0 + hot_mass * 600.0) / (cold_mass + hot_mass)).abs()
                > 10.0
        );
        assert_eq!(net.caloric_law().identity(), law.identity());
        assert_eq!(net.reference_temperature(), Temperature::new(300.0));
        assert_eq!(net.reference_density(), Density::new(1.0));
    });
}

#[test]
fn split_merge_and_downstream_exchange_share_the_original_channel_and_watt_ledger() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let law = source(800.0, 2.0, 275.0);
        let net = network(
            cx,
            &flow,
            law.clone(),
            vec![
                heated("left-wall", 3000.0),
                heated("right-wall", 1000.0),
                heated("tail-wall", 2500.0),
            ],
            &[inlet(0, 300.0)],
        );
        let walls = [500.0, 650.0, 320.0];
        let result = net.march(cx, &walls).unwrap();
        let masses: Vec<_> = flow
            .branches
            .iter()
            .map(|branch| branch.flow.value())
            .collect();
        let left = result.branches[0].outlet_temperature_k.unwrap();
        let right = result.branches[1].outlet_temperature_k.unwrap();
        let mixed = after_heat(
            300.0,
            (masses[0] * delta_h(300.0, left) + masses[1] * delta_h(300.0, right))
                / (masses[0] + masses[1]),
        );
        close(result.node_temperatures_k[1].unwrap(), mixed, 2e-10);
        close(
            result.branches[2].inlet_temperature_k.unwrap(),
            mixed,
            2e-10,
        );
        let tail = CaloricAirPath::new(
            mixed,
            masses[2],
            law,
            vec![AirSegment::new("tail-wall", 1.0, 2500.0).unwrap()],
        )
        .unwrap()
        .march(&walls[2..])
        .unwrap();
        close(
            result.node_temperatures_k[2].unwrap(),
            tail.outlet_temperature_k,
            2e-10,
        );
        close(
            result.wall_heat_rate_w,
            masses[2] * delta_h(300.0, result.node_temperatures_k[2].unwrap()),
            2e-7,
        );
        close(result.external_heat_gain_w, result.wall_heat_rate_w, 2e-7);
        assert_eq!(net.initial_references(cx).unwrap(), vec![300.0; 3]);
        assert_eq!(net.regions(), vec!["left-wall", "right-wall", "tail-wall"]);
        assert_eq!(result, net.march(cx, &walls).unwrap());
    });
}

#[test]
fn native_constant_cp_results_and_reversed_source_rows_remain_consistent() {
    with_cx(|cx| {
        let flow = LossGraph::new(2, vec![edge("reverse", 1, 0, 1.0)])
            .unwrap()
            .solve(
                &[boundary(0, 1.0), boundary(1, 0.0)],
                hydraulic_config(),
                cx,
            )
            .unwrap();
        let models = vec![BranchThermalModel::Exchange(vec![
            AirSegment::new("near-sink", 1.0, 700.0).unwrap(),
            AirSegment::new("near-source", 1.0, 300.0).unwrap(),
        ])];
        let native = TransportNetwork::new(
            cx,
            &flow,
            TransportAir {
                density: Density::new(1.0),
                specific_heat_j_kg_k: 1000.0,
            },
            models.clone(),
            &[inlet(0, 300.0)],
            budgets(),
        )
        .unwrap();
        let expected = native.march(cx, &[600.0, 320.0]).unwrap();
        for datum in [275.0, 650.0] {
            let net = network(
                cx,
                &flow,
                source(1000.0, 0.0, datum),
                models.clone(),
                &[inlet(0, 300.0)],
            );
            let result = net.march(cx, &[600.0, 320.0]).unwrap();
            assert_eq!(net.regions(), vec!["near-source", "near-sink"]);
            close(
                result.node_temperatures_k[1].unwrap(),
                expected.node_temperatures_k[1].unwrap(),
                2e-10,
            );
            close(result.wall_heat_rate_w, expected.wall_heat_rate_w, 2e-7);
            for (actual, reference) in result
                .reference_temperatures_k
                .iter()
                .zip(&expected.reference_temperatures_k)
            {
                close(*actual, *reference, 2e-10);
            }
        }
        let a = network(
            cx,
            &flow,
            source(800.0, 2.0, 275.0),
            models.clone(),
            &[inlet(0, 300.0)],
        )
        .march(cx, &[600.0, 320.0])
        .unwrap();
        let b = network(
            cx,
            &flow,
            source(800.0, 2.0, 650.0),
            models,
            &[inlet(0, 300.0)],
        )
        .march(cx, &[600.0, 320.0])
        .unwrap();
        close(
            a.node_temperatures_k[1].unwrap(),
            b.node_temperatures_k[1].unwrap(),
            2e-10,
        );
        close(a.wall_heat_rate_w, b.wall_heat_rate_w, 2e-7);
        assert_eq!(a.enthalpy_reference_k, 300.0);
        assert_eq!(b.enthalpy_reference_k, 300.0);
    });
}

#[test]
fn retained_weak_exchange_heat_reaches_the_external_network_ledger() {
    with_cx(|cx| {
        let flow = LossGraph::new(2, vec![edge("weak", 0, 1, 1.0)])
            .unwrap()
            .solve(
                &[boundary(0, 1.0), boundary(1, 0.0)],
                hydraulic_config(),
                cx,
            )
            .unwrap();
        let net = CaloricTransportNetwork::new(
            cx,
            &flow,
            Density::new(1.0),
            Temperature::new(300.0),
            source(1000.0, 0.0, 650.0),
            vec![heated("wall", 1e-6)],
            &[inlet(0, 300.0)],
            TransportConfig {
                absolute_heat_tolerance_w: 1e-14,
                relative_heat_tolerance: 1e-12,
                ..budgets()
            },
        )
        .unwrap();
        let result = net.march(cx, &[400.0]).unwrap();
        let expected = 1e-4 * (1.0 - 0.5e-9);
        close(result.wall_heat_rate_w, expected, 2e-16);
        close(result.external_heat_gain_w, expected, 2e-16);
        close(result.heat_imbalance_w, 0.0, 1e-14);
    });
}

#[test]
fn source_support_and_actual_hydraulic_residuals_refuse_without_repair() {
    with_cx(|cx| {
        let mut flow = merge_flow(cx);
        let law = source(800.0, 2.0, 300.0);
        assert!(matches!(
            CaloricTransportNetwork::new(
                cx,
                &flow,
                Density::new(1.0),
                Temperature::new(300.0),
                law.clone(),
                vec![BranchThermalModel::Adiabatic; 3],
                &[inlet(0, 750.0)],
                budgets()
            ),
            Err(TransportError::Airflow(AirflowError::CaloricSource(_)))
        ));
        let net = network(
            cx,
            &flow,
            law.clone(),
            vec![
                heated("a", 100_000.0),
                BranchThermalModel::Adiabatic,
                BranchThermalModel::Adiabatic,
            ],
            &[inlet(0, 300.0)],
        );
        assert!(matches!(
            net.march(cx, &[1000.0]),
            Err(TransportError::Airflow(AirflowError::CaloricSource(_)))
        ));
        assert!(net.march(cx, &[]).is_err());
        assert!(net.march(cx, &[f64::NAN]).is_err());

        flow.branches[2].flow = VolumetricFlowRate::new(flow.branches[2].flow.value() * 1.01);
        let loose = CaloricTransportNetwork::new(
            cx,
            &flow,
            Density::new(1.0),
            Temperature::new(300.0),
            law,
            vec![
                heated("a", 3000.0),
                heated("b", 1000.0),
                BranchThermalModel::Adiabatic,
            ],
            &[inlet(0, 300.0)],
            TransportConfig {
                absolute_flow_tolerance: VolumetricFlowRate::new(0.1),
                ..budgets()
            },
        )
        .unwrap();
        match loose.march(cx, &[500.0, 500.0]) {
            Err(TransportError::HeatImbalance {
                residual_w,
                hydraulic_defect_w,
                ..
            }) => {
                assert!(residual_w > 1.0);
                close(residual_w, hydraulic_defect_w, 2e-7);
            }
            other => panic!("the source heat defect must remain visible: {other:?}"),
        }
    });
}

fn solid(references: &[f64], conductances: &[f64], power: f64) -> Vec<SolidRegionState> {
    let wall = (power
        + references
            .iter()
            .zip(conductances)
            .map(|(r, g)| r * g)
            .sum::<f64>())
        / conductances.iter().sum::<f64>();
    references
        .iter()
        .zip(conductances)
        .enumerate()
        .map(|(i, (&reference, &g))| SolidRegionState {
            region: format!("wall-{i}"),
            area_m2: 1.0,
            mean_wall_temperature_k: wall,
            heat_rate_w: g * (wall - reference),
            mean_reference_temperature_k: Some(reference),
        })
        .collect()
}
fn coupled_models(conductances: &[f64]) -> Vec<BranchThermalModel> {
    conductances
        .iter()
        .enumerate()
        .map(|(i, &g)| heated(&format!("wall-{i}"), g))
        .collect()
}
fn coupling() -> ConjugateConfig {
    ConjugateConfig {
        max_iterations: 200,
        temperature_tolerance_k: 1e-10,
        balance_tolerance_w: 1e-6,
        balance_relative_tolerance: 1e-12,
        ..ConjugateConfig::default()
    }
}

#[test]
fn shared_solid_and_iqn_close_the_source_caloric_split_merge_network() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let conductances = [3000.0, 1000.0, 4000.0];
        let net = network(
            cx,
            &flow,
            source(800.0, 2.0, 650.0),
            coupled_models(&conductances),
            &[inlet(0, 300.0)],
        );
        let power = 60_000.0;
        let expected = after_heat(300.0, power / flow.branches[2].flow.value());
        let mut calls = 0;
        let solved = net
            .solve(cx, &coupling(), |_, references| {
                calls += 1;
                Ok(solid(references, &conductances, power))
            })
            .unwrap();
        assert_eq!(calls, solved.iterations);
        close(
            solved.transport.node_temperatures_k[2].unwrap(),
            expected,
            2e-8,
        );
        close(solved.transport.wall_heat_rate_w, power, 3e-6);
        for (state, &used) in solved.solid.iter().zip(&solved.reference_temperatures_k) {
            assert_eq!(state.mean_reference_temperature_k, Some(used));
        }
        for balance in &solved.region_balances {
            assert!(balance.imbalance_w.abs() <= coupling().balance_tolerance_w);
        }
        let accelerated = net
            .solve_iqn(cx, &coupling(), IqnIlsConfig::default(), |_, references| {
                Ok(solid(references, &conductances, power))
            })
            .unwrap();
        close(
            accelerated.transport.node_temperatures_k[2].unwrap(),
            expected,
            2e-8,
        );
        close(accelerated.transport.wall_heat_rate_w, power, 3e-6);
    });
}

#[test]
fn source_coupling_cancellation_resumes_the_native_fixed_point() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let flow = with_cx(merge_flow);
    let conductances = [3000.0, 1000.0, 4000.0];
    let law = source(800.0, 2.0, 300.0);
    let full = with_cx(|cx| {
        network(
            cx,
            &flow,
            law.clone(),
            coupled_models(&conductances),
            &[inlet(0, 300.0)],
        )
        .solve(cx, &coupling(), |_, references| {
            Ok(solid(references, &conductances, 60_000.0))
        })
        .unwrap()
    });
    let gate = CancelGate::new_clock_free();
    let refused = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 713,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let net = network(
            &cx,
            &flow,
            law.clone(),
            coupled_models(&conductances),
            &[inlet(0, 300.0)],
        );
        let mut calls = 0;
        net.solve(&cx, &coupling(), |_, references| {
            calls += 1;
            if calls == 4 {
                gate.request();
            }
            Ok(solid(references, &conductances, 60_000.0))
        })
        .unwrap_err()
    });
    let TransportError::Airflow(AirflowError::Cancelled {
        iteration,
        references_k,
    }) = refused
    else {
        panic!("source coupling must retain its complete reference vector");
    };
    assert_eq!(iteration, 3);
    let resumed = with_cx(|cx| {
        network(
            cx,
            &flow,
            law,
            coupled_models(&conductances),
            &[inlet(0, 300.0)],
        )
        .solve_from(cx, &coupling(), &references_k, |_, references| {
            Ok(solid(references, &conductances, 60_000.0))
        })
        .unwrap()
    });
    assert_eq!(
        full.reference_temperatures_k,
        resumed.reference_temperatures_k
    );
    assert_eq!(full.transport, resumed.transport);
    assert_eq!(full.iterations, iteration + resumed.iterations);
}

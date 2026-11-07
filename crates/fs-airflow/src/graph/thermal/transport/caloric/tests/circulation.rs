//! Original variable-Cp heat through actual buoyant closed hydraulic graphs.

use super::*;
use fs_material::fluid_transport::FluidTransportState;

fn controls() -> CaloricCirculationConfig {
    CaloricCirculationConfig {
        temperature_tolerance_k: 1e-9,
        max_iterations: 32,
    }
}

fn loop_flow(cx: &Cx<'_>, reference: f64, reverse: bool) -> GraphSolution {
    let cp = 800.0 + 2.0 * reference;
    let state = |rho| FluidTransportState::try_new(rho, cp, 0.1, 0.001).unwrap();
    let cold = if reverse {
        edge("cold", 0, 1, 1.0)
    } else {
        edge("cold", 1, 0, 1.0)
    };
    let mut config = hydraulic_config();
    config.absolute_flow_tolerance = VolumetricFlowRate::new(1e-13);
    config.relative_flow_tolerance = 1e-13;
    LossGraph::new(2, vec![edge("hot", 0, 1, 1.0), cold])
        .unwrap()
        .with_boussinesq_gravity(
            &state(1.0),
            &[state(0.875), state(1.0)],
            &[[0.0; 3], [0.0, 0.0, 1.0]],
            [0.0, 0.0, -8.0],
            0.2,
            cx,
        )
        .unwrap()
        .solve(&[boundary(0, 0.0)], config, cx)
        .unwrap()
}

fn loop_models(mass: f64, exposures: [f64; 2]) -> Vec<BranchThermalModel> {
    vec![
        heated("hot wall", mass * exposures[0]),
        heated("cold wall", mass * exposures[1]),
    ]
}

fn closed<'a>(
    cx: &Cx<'_>,
    flow: &'a GraphSolution,
    datum: f64,
    reference: f64,
    models: Vec<BranchThermalModel>,
) -> CaloricTransportNetwork<'a> {
    CaloricTransportNetwork::new_circulating(
        cx,
        flow,
        Density::new(1.0),
        Temperature::new(reference),
        source(800.0, 2.0, datum),
        models,
        budgets(),
        controls(),
    )
    .unwrap()
}

// Independent analytic integral for Cp(T)=800+2T. Bisection solves the
// original wall-exchange equation without querying the material or channels.
fn exchange(inlet: f64, wall: f64, exposure: f64) -> f64 {
    if inlet == wall {
        return inlet;
    }
    let (mut lo, mut hi) = (inlet.min(wall), inlet.max(wall));
    for _ in 0..100 {
        let t = lo + 0.5 * (hi - lo);
        let integral =
            -(800.0 + 2.0 * wall) * ((wall - t) / (wall - inlet)).ln() - 2.0 * (t - inlet);
        if (integral < exposure) == (wall > inlet) {
            lo = t;
        } else {
            hi = t;
        }
    }
    lo + 0.5 * (hi - lo)
}

fn fixed_point(lower: f64, upper: f64, mut returned: impl FnMut(f64) -> f64) -> f64 {
    let (mut lo, mut hi) = (lower, upper);
    for _ in 0..100 {
        let t = lo + 0.5 * (hi - lo);
        if returned(t) > t {
            lo = t;
        } else {
            hi = t;
        }
    }
    lo + 0.5 * (hi - lo)
}

#[test]
fn variable_cp_closed_loop_matches_independent_calorimetry_and_reversed_flow() {
    with_cx(|cx| {
        let exposure = [850.0, 1200.0];
        let walls = [500.0, 280.0];
        let bottom = fixed_point(walls[1], walls[0], |t| {
            exchange(exchange(t, walls[0], exposure[0]), walls[1], exposure[1])
        });
        let top = exchange(bottom, walls[0], exposure[0]);
        for reverse in [false, true] {
            let flow = loop_flow(cx, 300.0, reverse);
            let mass = flow.branches[0].flow.value();
            let net = closed(cx, &flow, 275.0, 300.0, loop_models(mass, exposure));
            let result = net.march(cx, &walls).unwrap();
            close(result.node_temperatures_k[0].unwrap(), bottom, 2e-8);
            close(result.node_temperatures_k[1].unwrap(), top, 2e-8);
            let hot = result.branches[0].march.as_ref().unwrap().total_heat_rate_w;
            let cold = result.branches[1].march.as_ref().unwrap().total_heat_rate_w;
            close(hot, mass * delta_h(bottom, top), 3e-6);
            close(hot, -cold, 2e-6);
            close(result.wall_heat_rate_w, 0.0, 2e-6);
            assert_eq!(result.external_heat_gain_w, 0.0);
            assert!(result.recirculation.is_none());
            assert_eq!(net.initial_references(cx).unwrap(), vec![300.0; 2]);
            assert_eq!(result, net.march(cx, &walls).unwrap());
            let frozen = TransportNetwork::new_circulating(
                cx,
                &flow,
                TransportAir {
                    density: Density::new(1.0),
                    specific_heat_j_kg_k: 1400.0,
                },
                loop_models(mass, exposure),
                Temperature::new(300.0),
                budgets(),
            )
            .unwrap()
            .march(cx, &walls)
            .unwrap();
            assert!((frozen.node_temperatures_k[0].unwrap() - bottom).abs() > 0.1);
        }
    });
}

#[test]
fn unequal_closed_paths_mix_mass_enthalpy_before_the_return_channel() {
    with_cx(|cx| {
        let state = |rho| FluidTransportState::try_new(rho, 1400.0, 0.1, 0.001).unwrap();
        let mut config = hydraulic_config();
        config.absolute_flow_tolerance = VolumetricFlowRate::new(1e-14);
        config.relative_flow_tolerance = 1e-14;
        let flow = LossGraph::new(
            3,
            vec![
                edge("return", 2, 0, 1.0),
                edge("right rise", 0, 2, 1.0),
                edge("left rise", 0, 1, 1.0),
                edge("top", 1, 2, 1.0),
            ],
        )
        .unwrap()
        .with_boussinesq_gravity(
            &state(1.0),
            &[state(1.0), state(0.9375), state(0.875), state(1.0)],
            &[[0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]],
            [0.0, 0.0, -8.0],
            0.2,
            cx,
        )
        .unwrap()
        .solve(&[boundary(0, 0.0)], config, cx)
        .unwrap();
        let mass: Vec<_> = flow.branches.iter().map(|b| b.flow.value()).collect();
        let weight = mass[2] / (mass[1] + mass[2]);
        assert!((weight - 0.5).abs() > 0.1);
        let mixed = |bottom| {
            let left = exchange(bottom, 620.0, 600.0);
            let right = exchange(bottom, 500.0, 800.0);
            after_heat(
                280.0,
                weight * delta_h(280.0, left) + (1.0 - weight) * delta_h(280.0, right),
            )
        };
        let bottom = fixed_point(280.0, 620.0, |t| exchange(mixed(t), 280.0, 1000.0));
        let expected_mix = mixed(bottom);
        let net = closed(
            cx,
            &flow,
            650.0,
            300.0,
            vec![
                heated("cold", mass[0] * 1000.0),
                heated("right", mass[1] * 800.0),
                heated("left", mass[2] * 600.0),
                BranchThermalModel::Adiabatic,
            ],
        );
        let result = net.march(cx, &[280.0, 500.0, 620.0]).unwrap();
        close(result.node_temperatures_k[0].unwrap(), bottom, 3e-8);
        close(
            result.node_temperatures_k[1].unwrap(),
            exchange(bottom, 620.0, 600.0),
            3e-8,
        );
        close(result.node_temperatures_k[2].unwrap(), expected_mix, 3e-8);
        close(result.wall_heat_rate_w, 0.0, 2e-6);
        let mean_temperature = weight * exchange(bottom, 620.0, 600.0)
            + (1.0 - weight) * exchange(bottom, 500.0, 800.0);
        assert!((expected_mix - mean_temperature).abs() > 0.05);
    });
}

#[test]
fn caloric_datum_and_physical_snapshot_reference_do_not_create_a_heat_reservoir() {
    with_cx(|cx| {
        let mut previous = None;
        for (datum, reference) in [(275.0, 300.0), (650.0, 300.0), (275.0, 350.0)] {
            let flow = loop_flow(cx, reference, false);
            let mass = flow.branches[0].flow.value();
            let net = closed(
                cx,
                &flow,
                datum,
                reference,
                loop_models(mass, [850.0, 1200.0]),
            );
            let result = net.march(cx, &[500.0, 280.0]).unwrap();
            assert_eq!(result.enthalpy_reference_k, reference);
            assert_eq!(net.reference_temperature().value(), reference);
            assert_eq!(net.caloric_law().datum_temperature_k(), datum);
            if let Some((bottom, top, heat)) = previous {
                close(result.node_temperatures_k[0].unwrap(), bottom, 3e-8);
                close(result.node_temperatures_k[1].unwrap(), top, 3e-8);
                close(
                    result.branches[0].march.as_ref().unwrap().total_heat_rate_w,
                    heat,
                    4e-6,
                );
            }
            previous = Some((
                result.node_temperatures_k[0].unwrap(),
                result.node_temperatures_k[1].unwrap(),
                result.branches[0].march.as_ref().unwrap().total_heat_rate_w,
            ));
        }
    });
}

#[test]
fn weak_wall_anchors_require_a_small_newton_correction_not_only_small_mixer_watts() {
    with_cx(|cx| {
        let flow = loop_flow(cx, 300.0, false);
        let mass = flow.branches[0].flow.value();
        let mut heat_budget = budgets();
        heat_budget.absolute_heat_tolerance_w = 1.0;
        let net = CaloricTransportNetwork::new_circulating(
            cx,
            &flow,
            Density::new(1.0),
            Temperature::new(300.0),
            source(800.0, 2.0, 650.0),
            loop_models(mass, [1e-5, 2e-5]),
            heat_budget,
            CaloricCirculationConfig {
                temperature_tolerance_k: 1e-5,
                max_iterations: 32,
            },
        )
        .unwrap();
        // The initial 450 K seed already satisfies a 1 W mixer gate and a
        // 1e-5 K raw mixer-temperature gate, but is 50 K from equilibrium.
        // In the weak-coupling limit UA_hot(T_hot-T)=UA_cold(T-T_cold).
        let result = net.march(cx, &[600.0, 300.0]).unwrap();
        for value in result.node_temperatures_k {
            close(value.unwrap(), 400.0, 2e-5);
        }
    });
}

#[test]
fn closed_source_admission_preserves_anchor_mass_support_and_work_refusals() {
    with_cx(|cx| {
        let flow = loop_flow(cx, 300.0, false);
        let mass = flow.branches[0].flow.value();
        let construct = |models, controls| {
            CaloricTransportNetwork::new_circulating(
                cx,
                &flow,
                Density::new(1.0),
                Temperature::new(300.0),
                source(800.0, 2.0, 275.0),
                models,
                budgets(),
                controls,
            )
        };
        assert!(construct(vec![BranchThermalModel::Adiabatic; 2], controls()).is_err());
        for control in [
            CaloricCirculationConfig {
                temperature_tolerance_k: 0.0,
                max_iterations: 32,
            },
            CaloricCirculationConfig {
                temperature_tolerance_k: 1e-9,
                max_iterations: 0,
            },
        ] {
            assert!(construct(loop_models(mass, [850.0, 1200.0]), control).is_err());
        }
        let net = construct(loop_models(mass, [850.0, 1200.0]), controls()).unwrap();
        assert!(net.march(cx, &[900.0; 2]).is_err());
        assert!(net.march(cx, &[500.0]).is_err());
        let limited = construct(
            loop_models(mass, [850.0, 1200.0]),
            CaloricCirculationConfig {
                max_iterations: 1,
                ..controls()
            },
        )
        .unwrap();
        assert!(matches!(
            limited.march(cx, &[500.0, 280.0]),
            Err(TransportError::CaloricCirculationNotConverged { iterations: 1, .. })
        ));
        assert!(
            net.with_recirculation(
                cx,
                vec![],
                CaloricRecirculationConfig {
                    temperature_tolerance_k: 1e-9,
                    max_iterations: 32
                }
            )
            .is_err()
        );
        let open = LossGraph::new(2, vec![edge("open", 0, 1, 1.0)])
            .unwrap()
            .solve(
                &[boundary(0, 1.0), boundary(1, 0.0)],
                hydraulic_config(),
                cx,
            )
            .unwrap();
        assert!(
            CaloricTransportNetwork::new_circulating(
                cx,
                &open,
                Density::new(1.0),
                Temperature::new(300.0),
                source(800.0, 2.0, 275.0),
                vec![heated("wall", 1000.0)],
                budgets(),
                controls()
            )
            .is_err()
        );
    });
}

#[test]
fn shared_solid_driver_closes_hot_input_and_cold_removal_on_the_same_closed_loop() {
    with_cx(|cx| {
        let flow = loop_flow(cx, 300.0, false);
        let mass = flow.branches[0].flow.value();
        let conductances = [mass * 850.0, mass * 1200.0];
        let net = closed(cx, &flow, 275.0, 300.0, loop_models(mass, [850.0, 1200.0]));
        let power = 10_000.0;
        let mut coupling = coupling();
        coupling.temperature_tolerance_k = 1e-7;
        coupling.balance_tolerance_w = 2e-4;
        let result = net
            .solve(cx, &coupling, |_, references| {
                Ok(vec![
                    SolidRegionState {
                        region: "hot wall".into(),
                        area_m2: 1.0,
                        mean_wall_temperature_k: references[0] + power / conductances[0],
                        heat_rate_w: power,
                        mean_reference_temperature_k: Some(references[0]),
                    },
                    SolidRegionState {
                        region: "cold wall".into(),
                        area_m2: 1.0,
                        mean_wall_temperature_k: 300.0,
                        heat_rate_w: conductances[1] * (300.0 - references[1]),
                        mean_reference_temperature_k: Some(references[1]),
                    },
                ])
            })
            .unwrap();
        close(
            result.transport.branches[0]
                .march
                .as_ref()
                .unwrap()
                .total_heat_rate_w,
            power,
            2e-4,
        );
        close(
            result.transport.branches[1]
                .march
                .as_ref()
                .unwrap()
                .total_heat_rate_w,
            -power,
            2e-4,
        );
        close(result.transport.wall_heat_rate_w, 0.0, 2e-6);
    });
}

#[test]
fn cancellation_publishes_no_closed_loop_and_preserves_replay() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    with_cx(|cx| {
        let flow = loop_flow(cx, 300.0, false);
        let mass = flow.branches[0].flow.value();
        let net = closed(cx, &flow, 275.0, 300.0, loop_models(mass, [850.0, 1200.0]));
        let expected = net.march(cx, &[500.0, 280.0]).unwrap();
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 719,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            gate.request();
            assert!(matches!(
                net.march(&cancelled, &[500.0, 280.0]),
                Err(TransportError::Interrupted)
            ));
        });
        assert_eq!(expected, net.march(cx, &[500.0, 280.0]).unwrap());
    });
}

#[test]
fn outside_source_baths_use_a_feasible_declared_reference_seed() {
    with_cx(|cx| {
        let flow = loop_flow(cx, 300.0, false);
        let mass = flow.branches[0].flow.value();
        // The source supports 250..700 K. The preferred wall midpoint clips
        // to 700 K and immediately overheats that trial, but the weak heater's
        // real circulating fluid remains near the 300 K cooler.
        let net = closed(cx, &flow, 650.0, 300.0, loop_models(mass, [1.0, 1000.0]));
        let bottom = fixed_point(300.0, 350.0, |t| {
            exchange(exchange(t, 2000.0, 1.0), 300.0, 1000.0)
        });
        let top = exchange(bottom, 2000.0, 1.0);
        let result = net.march(cx, &[2000.0, 300.0]).unwrap();
        close(result.node_temperatures_k[0].unwrap(), bottom, 3e-8);
        close(result.node_temperatures_k[1].unwrap(), top, 3e-8);
        assert!(bottom > 300.0 && top < 305.0);
        close(result.wall_heat_rate_w, 0.0, 2e-6);
    });
}

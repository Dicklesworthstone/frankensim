//! Native finite exhaust recovery: physical references, mixing, solids and refusals.

use super::*;
use crate::conjugate::caloric::exchanger::FlowArrangement;
use crate::graph::thermal::transport::recirculation::RecirculationLink;

fn recovery_config() -> CaloricHeatRecoveryConfig {
    CaloricHeatRecoveryConfig {
        temperature_tolerance_k: 1e-9,
        ..CaloricHeatRecoveryConfig::default()
    }
}

fn link(supply_node: usize, exhaust_node: usize, conductance: f64) -> CaloricHeatRecoveryLink {
    CaloricHeatRecoveryLink {
        supply_node,
        exhaust_node,
        arrangement: FlowArrangement::Counterflow,
        conductance_w_per_k: conductance,
    }
}

fn single_flow(cx: &Cx<'_>, reversed: bool) -> GraphSolution {
    let (from, to) = if reversed { (1, 0) } else { (0, 1) };
    LossGraph::new(2, vec![edge("channel", from, to, 1.0)])
        .unwrap()
        .solve(
            &[boundary(0, 1.0), boundary(1, 0.0)],
            hydraulic_config(),
            cx,
        )
        .unwrap()
}

fn one<'a>(
    cx: &Cx<'_>,
    flow: &'a GraphSolution,
    fresh: f64,
    wall_conductance: f64,
    law: Arc<IsobaricHeatCapacityLaw>,
) -> CaloricTransportNetwork<'a> {
    network(
        cx,
        flow,
        law,
        vec![heated("wall-0", wall_conductance)],
        &[inlet(0, fresh)],
    )
    .with_heat_recovery(cx, vec![link(0, 1, 1000.0)], recovery_config())
    .unwrap()
}

fn check_report(march: &TransportMarch, tolerance: f64, config: CaloricHeatRecoveryConfig) {
    let recovery = march
        .heat_recovery
        .as_ref()
        .expect("native heat recovery report");
    assert!(march.recirculation.is_none());
    assert!(march.heat_imbalance_w.abs() <= tolerance);
    assert!(recovery.heat_imbalance_w.abs() <= tolerance);
    assert!(recovery.supply_heat_residual_w <= tolerance);
    assert!(recovery.temperature_correction_bound_k <= config.temperature_tolerance_k);
    assert!(
        recovery.exchanger_entropy_generation_w_per_k
            >= -recovery.exchanger_entropy_roundoff_allowance_w_per_k
    );
    for stream in &recovery.streams {
        assert!(
            stream.exhaust_inlet_energy_residual_w.abs()
                <= tolerance + stream.exhaust_inlet_energy_roundoff_allowance_w
        );
        assert!(
            stream.exchanger.first_endpoint_energy_residual_w.abs()
                <= tolerance + stream.exchanger.first_endpoint_energy_roundoff_allowance_w
        );
        assert!(
            stream.exchanger.second_endpoint_energy_residual_w.abs()
                <= tolerance + stream.exchanger.second_endpoint_energy_roundoff_allowance_w
        );
    }
}

#[test]
fn finite_recovery_matches_analytical_heating_cooling_and_reversed_hydraulics() {
    with_cx(|cx| {
        for reversed in [false, true] {
            let flow = single_flow(cx, reversed);
            assert_eq!(flow.branches[0].flow.value() < 0.0, reversed);
            for (fresh, wall) in [(300.0, 600.0), (600.0, 300.0)] {
                let net = one(cx, &flow, fresh, 1000.0, source(1000.0, 0.0, 650.0));
                assert_eq!(net.heat_recovery_links(), &[link(0, 1, 1000.0)]);
                assert!(net.recirculation_links().is_empty());
                let march = net.march(cx, &[wall]).unwrap();
                let report = march.heat_recovery.as_ref().unwrap();
                let recovered = &report.streams[0];
                let attenuation = (-1.0_f64).exp();
                let effectiveness = 0.5;
                let expected_supply = ((1.0 - effectiveness) * fresh
                    + effectiveness * (1.0 - attenuation) * wall)
                    / (1.0 - effectiveness * attenuation);
                let expected_exhaust = wall + (expected_supply - wall) * attenuation;
                let expected_heat = 1000.0 * (expected_supply - fresh);
                let wall_heat = 1000.0 * (expected_exhaust - expected_supply);
                close(recovered.supply_temperature_k, expected_supply, 3e-8);
                close(march.node_temperatures_k[0].unwrap(), expected_supply, 3e-8);
                close(
                    march.node_temperatures_k[1].unwrap(),
                    expected_exhaust,
                    3e-8,
                );
                close(
                    recovered.exchanger.first_outlet_temperature_k,
                    expected_exhaust - expected_heat / 1000.0,
                    4e-8,
                );
                close(
                    recovered.exchanger.heat_rate_first_to_second_w,
                    expected_heat,
                    5e-5,
                );
                close(march.wall_heat_rate_w, wall_heat, 5e-5);
                close(report.external_heat_gain_w, wall_heat, 5e-5);
                close(recovered.supply_mass_flow_kg_per_s, 1.0, 1e-14);
                close(recovered.exhaust_mass_flow_kg_per_s, 1.0, 1e-14);
                let entropy = 1000.0
                    * ((expected_exhaust - expected_heat / 1000.0) / expected_exhaust).ln()
                    + 1000.0 * (expected_supply / fresh).ln();
                close(report.exchanger_entropy_generation_w_per_k, entropy, 2e-8);
                check_report(
                    &march,
                    budgets().absolute_heat_tolerance_w,
                    recovery_config(),
                );
                assert!(report.iterations >= 2);
            }
        }
    });
}

// Independent original Cp(T)=800+2T quadratic inverse and Simpson integral.
// This manufactures a 300 -> 400 K fresh path and a 500 K exhaust at Q=150 kW.
fn recovery_ua() -> f64 {
    let heat = 150_000.0;
    let integrand = |progress: f64| {
        let hot = after_heat(500.0, -progress);
        let cold = after_heat(300.0, heat - progress);
        1.0 / (hot - cold)
    };
    let step = heat / 16_384.0;
    let mut total = integrand(0.0) + integrand(heat);
    for i in 1..16_384 {
        total += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(step * i as f64);
    }
    total * step / 3.0
}

#[test]
fn variable_cp_recovery_preserves_original_calorimetry_and_datum_independence() {
    with_cx(|cx| {
        let flow = single_flow(cx, false);
        let branch_ua = 2000.0 * 2.0_f64.ln() - 200.0;
        let recovered_ua = recovery_ua();
        let mut previous: Option<TransportMarch> = None;
        for datum in [250.0, 650.0] {
            // The declared 0.1 mW physical gate is kept separate from source
            // endpoint rounding and each exchanger's constitutive allowance.
            let config = TransportConfig {
                absolute_heat_tolerance_w: 1e-4,
                ..budgets()
            };
            let recovery = CaloricHeatRecoveryConfig {
                temperature_tolerance_k: 1e-7,
                ..recovery_config()
            };
            let net = CaloricTransportNetwork::new(
                cx,
                &flow,
                Density::new(1.0),
                Temperature::new(300.0),
                source(800.0, 2.0, datum),
                vec![heated("wall-0", branch_ua)],
                &[inlet(0, 300.0)],
                config,
            )
            .unwrap()
            .with_heat_recovery(cx, vec![link(0, 1, recovered_ua)], recovery)
            .unwrap();
            let march = net.march(cx, &[600.0]).unwrap();
            let report = march.heat_recovery.as_ref().unwrap();
            close(march.node_temperatures_k[0].unwrap(), 400.0, 2e-7);
            close(march.node_temperatures_k[1].unwrap(), 500.0, 2e-7);
            close(
                report.streams[0].exchanger.heat_rate_first_to_second_w,
                150_000.0,
                5e-4,
            );
            close(march.wall_heat_rate_w, 170_000.0, 5e-4);
            close(report.external_heat_gain_w, 170_000.0, 5e-4);
            close(
                report.streams[0].exchanger.first_outlet_temperature_k,
                after_heat(500.0, -150_000.0),
                3e-7,
            );
            check_report(&march, config.absolute_heat_tolerance_w, recovery);
            if let Some(before) = previous {
                close(
                    march.node_temperatures_k[0].unwrap(),
                    before.node_temperatures_k[0].unwrap(),
                    2e-7,
                );
                close(
                    report.external_heat_gain_w,
                    before.heat_recovery.unwrap().external_heat_gain_w,
                    5e-4,
                );
            }
            previous = Some(march);
        }
    });
}

#[test]
fn native_split_merge_and_unequal_ports_keep_their_actual_mass_and_remaining_exhaust() {
    with_cx(|cx| {
        let flow = merge_flow(cx);
        let net = network(
            cx,
            &flow,
            source(1000.0, 0.0, 650.0),
            coupled_models(&[3000.0, 1000.0, 4000.0]),
            &[inlet(0, 300.0)],
        )
        .with_heat_recovery(cx, vec![link(0, 2, 4000.0)], recovery_config())
        .unwrap();
        let march = net.march(cx, &[600.0; 3]).unwrap();
        let a = (-2.0_f64).exp();
        let expected = (150.0 + 0.5 * (1.0 - a) * 600.0) / (1.0 - 0.5 * a);
        close(march.node_temperatures_k[0].unwrap(), expected, 3e-8);
        close(
            march.branches[0].outlet_temperature_k.unwrap(),
            march.branches[1].outlet_temperature_k.unwrap(),
            2e-11,
        );
        assert_eq!(net.regions(), vec!["wall-0", "wall-1", "wall-2"]);
        close(
            march.heat_recovery.as_ref().unwrap().streams[0].supply_mass_flow_kg_per_s,
            4.0,
            1e-13,
        );
        check_report(
            &march,
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        );

        // Three kg/s exhaust is recovered into the whole four kg/s supply.
        // The other one kg/s leaves through an independent unheated exhaust.
        let flow = LossGraph::new(
            3,
            vec![edge("heated", 0, 1, 1.0 / 9.0), edge("bypass", 0, 2, 1.0)],
        )
        .unwrap()
        .solve(
            &[boundary(0, 1.0), boundary(1, 0.0), boundary(2, 0.0)],
            hydraulic_config(),
            cx,
        )
        .unwrap();
        let net = network(
            cx,
            &flow,
            source(1000.0, 0.0, 650.0),
            vec![heated("wall-0", 3000.0), BranchThermalModel::Adiabatic],
            &[inlet(0, 300.0)],
        )
        .with_heat_recovery(cx, vec![link(0, 1, 1500.0)], recovery_config())
        .unwrap();
        let march = net.march(cx, &[600.0]).unwrap();
        let decay = (-0.125_f64).exp();
        let epsilon = (1.0 - decay) / (1.0 - 0.75 * decay);
        let receiving_fraction = 0.75 * epsilon;
        let attenuation = (-1.0_f64).exp();
        let expected_supply = ((1.0 - receiving_fraction) * 300.0
            + receiving_fraction * (1.0 - attenuation) * 600.0)
            / (1.0 - receiving_fraction * attenuation);
        let before = 600.0 + (expected_supply - 600.0) * attenuation;
        let heat = 4000.0 * (expected_supply - 300.0);
        let post = before - heat / 3000.0;
        let report = march.heat_recovery.as_ref().unwrap();
        close(report.streams[0].supply_mass_flow_kg_per_s, 4.0, 1e-13);
        close(report.streams[0].exhaust_mass_flow_kg_per_s, 3.0, 1e-13);
        close(march.node_temperatures_k[0].unwrap(), expected_supply, 3e-8);
        close(march.node_temperatures_k[2].unwrap(), expected_supply, 3e-8);
        close(
            report.external_heat_gain_w,
            3000.0 * (post - 300.0) + 1000.0 * (expected_supply - 300.0),
            5e-5,
        );
        check_report(
            &march,
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        );
    });
}

#[test]
fn the_original_shared_solid_and_iqn_drivers_close_with_finite_recovery() {
    with_cx(|cx| {
        let flow = single_flow(cx, false);
        let net = one(cx, &flow, 300.0, 1000.0, source(1000.0, 0.0, 650.0));
        let power = 60_000.0;
        let mut calls = 0;
        let solved = net
            .solve(cx, &coupling(), |_, references| {
                calls += 1;
                Ok(solid(references, &[1000.0], power))
            })
            .unwrap();
        assert_eq!(calls, solved.iterations);
        close(
            solved.transport.node_temperatures_k[0].unwrap(),
            360.0,
            2e-7,
        );
        close(
            solved.transport.node_temperatures_k[1].unwrap(),
            420.0,
            2e-7,
        );
        close(solved.transport.wall_heat_rate_w, power, 3e-6);
        close(
            solved
                .transport
                .heat_recovery
                .as_ref()
                .unwrap()
                .external_heat_gain_w,
            power,
            3e-6,
        );
        assert_eq!(
            solved.solid[0].mean_reference_temperature_k,
            Some(solved.reference_temperatures_k[0])
        );
        assert!(solved.region_balances[0].imbalance_w.abs() <= coupling().balance_tolerance_w);
        check_report(
            &solved.transport,
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        );
        let accelerated = net
            .solve_iqn(cx, &coupling(), IqnIlsConfig::default(), |_, references| {
                Ok(solid(references, &[1000.0], power))
            })
            .unwrap();
        close(
            accelerated.transport.node_temperatures_k[0].unwrap(),
            360.0,
            2e-7,
        );
        close(accelerated.transport.wall_heat_rate_w, power, 3e-6);
        check_report(
            &accelerated.transport,
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        );
    });
}

#[test]
fn retained_exchanger_gain_resolves_an_adiabatic_dag_with_a_very_weak_fresh_anchor() {
    with_cx(|cx| {
        let flow = single_flow(cx, false);
        let net = network(
            cx,
            &flow,
            source(1000.0, 0.0, 650.0),
            vec![BranchThermalModel::Adiabatic],
            &[inlet(0, 300.0)],
        )
        .with_heat_recovery(cx, vec![link(0, 1, 1e22)], recovery_config())
        .unwrap();
        let march = net.march(cx, &[]).unwrap();
        let report = march.heat_recovery.as_ref().unwrap();
        assert_eq!(report.streams[0].exchanger.heat_rate_first_to_second_w, 0.0);
        assert_eq!(
            report.streams[0]
                .response
                .heat_from_first_inlet_enthalpy_kg_per_s,
            1.0
        );
        assert!(report.streams[0].response.first_outlet_enthalpy_gain > 0.0);
        assert!(report.streams[0].response.first_outlet_enthalpy_gain < f64::EPSILON);
        assert_eq!(march.node_temperatures_k, vec![Some(300.0), Some(300.0)]);
        assert_eq!(report.iterations, 1);
        check_report(
            &march,
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        );
    });
}

#[test]
fn full_port_ownership_modes_source_limits_and_nonlinear_budget_refuse_transactionally() {
    with_cx(|cx| {
        let flow = single_flow(cx, false);
        let make = || {
            network(
                cx,
                &flow,
                source(1000.0, 0.0, 650.0),
                vec![heated("wall-0", 1000.0)],
                &[inlet(0, 300.0)],
            )
        };
        for links in [
            vec![link(1, 0, 1000.0)],
            vec![link(0, 2, 1000.0)],
            vec![link(0, 1, 1000.0), link(0, 1, 500.0)],
            vec![link(0, 1, -1.0)],
        ] {
            assert!(
                make()
                    .with_heat_recovery(cx, links, recovery_config())
                    .is_err()
            );
        }
        let bounded = make()
            .with_heat_recovery(
                cx,
                vec![link(0, 1, 1000.0)],
                CaloricHeatRecoveryConfig {
                    max_iterations: 1,
                    ..recovery_config()
                },
            )
            .unwrap();
        assert!(matches!(
            bounded.march(cx, &[600.0]),
            Err(TransportError::CaloricHeatRecoveryNotConverged { iterations: 1, .. })
        ));
        let returned = || {
            vec![RecirculationLink {
                supply_node: 0,
                return_node: 1,
                fraction: 0.2,
            }]
        };
        let return_config = CaloricRecirculationConfig {
            temperature_tolerance_k: 1e-9,
            max_iterations: 64,
        };
        assert!(
            make()
                .with_heat_recovery(cx, vec![link(0, 1, 1000.0)], recovery_config())
                .unwrap()
                .with_recirculation(cx, returned(), return_config)
                .is_err()
        );
        assert!(
            make()
                .with_recirculation(cx, returned(), return_config)
                .unwrap()
                .with_heat_recovery(cx, vec![link(0, 1, 1000.0)], recovery_config())
                .is_err()
        );
        let limited = one(cx, &flow, 650.0, 1000.0, source(1000.0, 0.0, 300.0));
        assert!(matches!(
            limited.march(cx, &[1000.0]),
            Err(TransportError::Airflow(AirflowError::CaloricSource(_)))
        ));
    });
}

#[test]
fn recovery_cancellation_preserves_the_native_error_and_a_later_march() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let flow = with_cx(|cx| single_flow(cx, false));
    let net = with_cx(|cx| one(cx, &flow, 300.0, 1000.0, source(1000.0, 0.0, 650.0)));
    let gate = CancelGate::new_clock_free();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 20,
                kernel_id: 720,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(matches!(
            net.march(&cx, &[600.0]),
            Err(TransportError::Interrupted)
        ));
    });
    with_cx(|cx| {
        check_report(
            &net.march(cx, &[600.0]).unwrap(),
            budgets().absolute_heat_tolerance_w,
            recovery_config(),
        )
    });
}

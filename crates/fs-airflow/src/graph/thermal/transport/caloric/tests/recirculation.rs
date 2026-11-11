//! G1 analytical return-loop balances and G3 source/order invariance.

use super::*;
use crate::graph::thermal::transport::recirculation::RecirculationLink;

fn feedback_config() -> CaloricRecirculationConfig {
    CaloricRecirculationConfig {
        temperature_tolerance_k: 1e-10,
        max_iterations: 64,
    }
}

fn link(supply_node: usize, return_node: usize, fraction: f64) -> RecirculationLink {
    RecirculationLink {
        supply_node,
        return_node,
        fraction,
    }
}

fn single_flow(cx: &Cx<'_>) -> GraphSolution {
    LossGraph::new(2, vec![edge("channel", 0, 1, 1.0)])
        .unwrap()
        .solve(
            &[boundary(0, 1.0), boundary(1, 0.0)],
            hydraulic_config(),
            cx,
        )
        .unwrap()
}

// Independent calorimetry for Cp(T) = a + bT. Neither expected temperature
// nor expected watts calls the material source or production channel solver.
fn affine_heat(coefficients: (f64, f64), initial: f64, final_t: f64) -> f64 {
    let (a, b) = coefficients;
    let delta = final_t - initial;
    (a + b * initial) * delta + 0.5 * b * delta * delta
}

fn return_oracle(
    coefficients: (f64, f64),
    fresh: f64,
    wall: f64,
    fraction: f64,
    mass: f64,
    conductance: f64,
) -> (f64, f64, f64) {
    let (a, b) = coefficients;
    let cp_fresh = a + b * fresh;
    let mixture = |outlet| {
        let heat = fraction * affine_heat(coefficients, fresh, outlet);
        fresh + 2.0 * heat / (cp_fresh + (cp_fresh * cp_fresh + 2.0 * b * heat).sqrt())
    };
    let mut lower = 0.0;
    let mut upper = 1.0;
    for _ in 0..90 {
        let fraction_to_wall = 0.5 * (lower + upper);
        let outlet = fresh + fraction_to_wall * (wall - fresh);
        let inlet = mixture(outlet);
        let change = outlet - inlet;
        // Integral from inlet to outlet of Cp(T)/(wall - T) dT = hA/m.
        // Its logarithmic primitive applies to heating and cooling alike.
        let integral = -b * change - (a + b * wall) * (-change / (wall - inlet)).ln_1p();
        if integral < conductance / mass {
            lower = fraction_to_wall;
        } else {
            upper = fraction_to_wall;
        }
    }
    let outlet = fresh + 0.5 * (lower + upper) * (wall - fresh);
    (
        mixture(outlet),
        outlet,
        (1.0 - fraction) * mass * affine_heat(coefficients, fresh, outlet),
    )
}

#[test]
fn constant_cp_returns_match_the_native_solver_and_analytical_feedback() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let mass = flow.branches[0].flow.value();
        let cp = 1000.0;
        let conductance = 1500.0;
        for fraction in [0.0, 0.25, 0.8, 0.9999] {
            let models = vec![heated("wall", conductance)];
            let native = TransportNetwork::new(
                cx,
                &flow,
                TransportAir {
                    density: Density::new(1.0),
                    specific_heat_j_kg_k: cp,
                },
                models.clone(),
                &[inlet(0, 300.0)],
                budgets(),
            )
            .unwrap()
            .with_recirculation(cx, vec![link(0, 1, fraction)], 1e-10)
            .unwrap()
            .march(cx, &[650.0])
            .unwrap();
            let net = network(
                cx,
                &flow,
                source(cp, 0.0, 650.0),
                models,
                &[inlet(0, 300.0)],
            )
            .with_recirculation(cx, vec![link(0, 1, fraction)], feedback_config())
            .unwrap();
            let result = net.march(cx, &[650.0]).unwrap();
            let attenuation = (-conductance / (mass * cp)).exp();
            let rise = (1.0 - attenuation) * 350.0 / (1.0 - fraction * attenuation);
            close(
                result.node_temperatures_k[0].unwrap(),
                300.0 + fraction * rise,
                2e-8,
            );
            close(result.node_temperatures_k[1].unwrap(), 300.0 + rise, 2e-8);
            close(
                result.wall_heat_rate_w,
                (1.0 - fraction) * mass * cp * rise,
                2e-5,
            );
            close(result.wall_heat_rate_w, native.wall_heat_rate_w, 2e-5);
            for (actual, expected) in result
                .node_temperatures_k
                .iter()
                .zip(&native.node_temperatures_k)
            {
                close(actual.unwrap(), expected.unwrap(), 2e-8);
            }
            close(net.initial_references(cx).unwrap()[0], 300.0, 1e-10);
            if fraction == 0.0 {
                assert!(result.recirculation.is_none());
            } else {
                let report = result.recirculation.as_ref().unwrap();
                assert_eq!(report.caloric_reference_temperature_k, Some(300.0));
                assert_eq!(
                    native
                        .recirculation
                        .as_ref()
                        .unwrap()
                        .caloric_reference_temperature_k,
                    None
                );
                close(report.external_heat_gain_w, result.wall_heat_rate_w, 2e-5);
                close(
                    report.supplies[0].fresh_mass_flow_kg_per_s,
                    (1.0 - fraction) * mass,
                    1e-14,
                );
                close(
                    report.supplies[0].return_mass_flow_kg_per_s,
                    fraction * mass,
                    1e-14,
                );
                close(report.streams[0].mass_flow_kg_per_s, fraction * mass, 1e-14);
                assert!(report.max_mixing_residual_k <= feedback_config().temperature_tolerance_k);
            }
        }
    });
}

#[test]
fn complete_exhaust_recovery_uses_mass_independently_of_reference_cp() {
    with_cx(|cx| {
        let flow = LossGraph::new(
            3,
            vec![edge("small", 0, 1, 25.0), edge("large", 0, 2, 100.0 / 49.0)],
        )
        .unwrap()
        .solve(
            &[boundary(0, 1.0), boundary(1, 0.0), boundary(2, 0.0)],
            hydraulic_config(),
            cx,
        )
        .unwrap();
        let small = flow.branches[0].flow.value();
        let large = flow.branches[1].flow.value();
        assert_eq!((small, large), (0.2, 0.7));
        let total = small + large;
        let fraction = small / total;
        assert_eq!(fraction * total, small);
        assert!(fraction.next_up() * total > small);
        // The same Cp(T)=600+T source has Cp=900 or 901 at the two
        // physical references. At Cp=900 the rounded capacity draw is
        // 180.00000000000003 W/K, above the available 180 W/K, despite
        // exact mass availability. A diagnostic Cp must not reject this loop.
        let law = source(600.0, 1.0, 300.0);
        let build = |reference| {
            CaloricTransportNetwork::new(
                cx,
                &flow,
                Density::new(1.0),
                Temperature::new(reference),
                law.clone(),
                vec![BranchThermalModel::Adiabatic; 2],
                &[inlet(0, 300.0)],
                budgets(),
            )
            .unwrap()
        };
        for (reference, cp) in [(300.0, 900.0), (301.0, 901.0)] {
            let result = build(reference)
                .with_recirculation(cx, vec![link(0, 1, fraction)], feedback_config())
                .unwrap()
                .march(cx, &[])
                .unwrap();
            assert_eq!(result.node_temperatures_k, vec![Some(300.0); 3]);
            assert_eq!(result.wall_heat_rate_w, 0.0);
            assert_eq!(result.external_heat_gain_w, 0.0);
            let report = result.recirculation.unwrap();
            assert_eq!(report.caloric_reference_temperature_k, Some(reference));
            close(report.supplies[0].fresh_mass_flow_kg_per_s, large, 1e-14);
            assert_eq!(report.supplies[0].return_mass_flow_kg_per_s, small);
            assert_eq!(report.streams[0].mass_flow_kg_per_s, small);
            // The stream capacity and small * cp are one product reached through different
            // roundings; 1 ulp at 180 W/K is 2.8e-14, so 1e-14 absolute could only be met bit for
            // bit. Two eps relative names the allowance.
            close(
                report.streams[0].capacity_w_per_k,
                small * cp,
                2.0 * f64::EPSILON * (small * cp).abs(),
            );
            assert_eq!(report.external_heat_gain_w, 0.0);
            assert_eq!(report.heat_imbalance_w, 0.0);
            // Admit the complete draw exactly, without an overdraw tolerance
            // or clipping that would also admit the next larger fraction.
            assert!(matches!(
                build(reference).with_recirculation(
                    cx,
                    vec![link(0, 1, fraction.next_up())],
                    feedback_config(),
                ),
                Err(TransportError::Node { node: 1, .. })
            ));
        }
    });
}

#[test]
fn weak_return_exchange_retains_heat_below_absolute_temperature_resolution() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let mass = flow.branches[0].flow.value();
        let cp = 1000.0;
        let conductance = 1e-6;
        let fraction = 0.5;
        let result = CaloricTransportNetwork::new(
            cx,
            &flow,
            Density::new(1.0),
            Temperature::new(300.0),
            source(cp, 0.0, 650.0),
            vec![heated("wall", conductance)],
            &[inlet(0, 300.0)],
            TransportConfig {
                absolute_heat_tolerance_w: 1e-14,
                relative_heat_tolerance: 1e-12,
                ..budgets()
            },
        )
        .unwrap()
        .with_recirculation(cx, vec![link(0, 1, fraction)], feedback_config())
        .unwrap()
        .march(cx, &[400.0])
        .unwrap();
        // Avoid subtracting exp(-NTU) from one in this independent oracle.
        let effectiveness = -(-conductance / (mass * cp)).exp_m1();
        let expected = mass * cp * (1.0 - fraction) * effectiveness * 100.0
            / ((1.0 - fraction) + fraction * effectiveness);
        close(result.wall_heat_rate_w, expected, 2e-16);
        close(result.external_heat_gain_w, expected, 2e-16);
        let report = result.recirculation.as_ref().unwrap();
        close(report.external_heat_gain_w, expected, 2e-16);
        close(report.heat_imbalance_w, 0.0, 1e-14);
        // Rounded K endpoints cannot reconstruct the retained source watts.
        let rounded_endpoint_heat =
            (1.0 - fraction) * mass * cp * (result.node_temperatures_k[1].unwrap() - 300.0);
        assert!((rounded_endpoint_heat - expected).abs() > 1e-12);
    });
}

#[test]
fn affine_cp_heating_and_cooling_close_original_enthalpy_and_wall_integrals() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let mass = flow.branches[0].flow.value();
        for (coefficients, fresh, wall, fraction) in [
            ((-900.0, 6.0), 300.0, 650.0, 0.75),
            ((-900.0, 6.0), 650.0, 300.0, 0.75),
            ((4500.0, -6.0), 300.0, 650.0, 0.85),
        ] {
            let result = network(
                cx,
                &flow,
                source(coefficients.0, coefficients.1, 650.0),
                vec![heated("wall", 1500.0)],
                &[inlet(0, fresh)],
            )
            .with_recirculation(cx, vec![link(0, 1, fraction)], feedback_config())
            .unwrap()
            .march(cx, &[wall])
            .unwrap();
            let (expected_inlet, expected_outlet, expected_heat) =
                return_oracle(coefficients, fresh, wall, fraction, mass, 1500.0);
            close(result.node_temperatures_k[0].unwrap(), expected_inlet, 2e-8);
            close(
                result.node_temperatures_k[1].unwrap(),
                expected_outlet,
                2e-8,
            );
            close(result.wall_heat_rate_w, expected_heat, 2e-5);
            let report = result.recirculation.as_ref().unwrap();
            close(report.external_heat_gain_w, expected_heat, 2e-5);
            close(report.heat_imbalance_w, 0.0, 2e-5);
            assert!(report.max_mixing_residual_k <= feedback_config().temperature_tolerance_k);
            let actual_inlet = result.node_temperatures_k[0].unwrap();
            let actual_outlet = result.node_temperatures_k[1].unwrap();
            close(
                mass * (affine_heat(coefficients, fresh, actual_inlet)
                    - fraction * affine_heat(coefficients, fresh, actual_outlet)),
                0.0,
                2e-5,
            );
        }
    });
}

#[test]
fn adiabatic_returns_mix_fresh_mass_enthalpy_independently_of_link_order() {
    with_cx(|cx| {
        let flow = LossGraph::new(3, vec![edge("cold", 0, 2, 1.0), edge("hot", 1, 2, 1.0)])
            .unwrap()
            .solve(
                &[boundary(0, 9.0), boundary(1, 4.0), boundary(2, 0.0)],
                hydraulic_config(),
                cx,
            )
            .unwrap();
        let run = |links| {
            network(
                cx,
                &flow,
                source(800.0, 2.0, 650.0),
                vec![BranchThermalModel::Adiabatic; 2],
                &[inlet(0, 300.0), inlet(1, 600.0)],
            )
            .with_recirculation(cx, links, feedback_config())
            .unwrap()
            .march(cx, &[])
            .unwrap()
        };
        let result = run(vec![link(0, 2, 0.5), link(1, 2, 0.25)]);
        assert_eq!(result, run(vec![link(1, 2, 0.25), link(0, 2, 0.5)]));
        let cold_mass = flow.branches[0].flow.value();
        let hot_mass = flow.branches[1].flow.value();
        let fresh_cold = 0.5 * cold_mass;
        let fresh_hot = 0.75 * hot_mass;
        let return_enthalpy = fresh_hot / (fresh_cold + fresh_hot) * delta_h(300.0, 600.0);
        let expected_return = after_heat(300.0, return_enthalpy);
        close(
            result.node_temperatures_k[2].unwrap(),
            expected_return,
            2e-8,
        );
        close(
            result.node_temperatures_k[0].unwrap(),
            after_heat(300.0, 0.5 * return_enthalpy),
            2e-8,
        );
        close(
            result.node_temperatures_k[1].unwrap(),
            after_heat(300.0, 0.75 * delta_h(300.0, 600.0) + 0.25 * return_enthalpy),
            2e-8,
        );
        assert!((expected_return - 450.0).abs() > 10.0);
        let report = result.recirculation.as_ref().unwrap();
        close(result.wall_heat_rate_w, 0.0, 1e-12);
        close(report.external_heat_gain_w, 0.0, 2e-6);
        close(report.heat_imbalance_w, 0.0, 2e-6);
        assert_eq!(report.caloric_reference_temperature_k, Some(300.0));
        for ((supply, stream), (mass, fraction)) in report
            .supplies
            .iter()
            .zip(&report.streams)
            .zip([(cold_mass, 0.5), (hot_mass, 0.25)])
        {
            close(
                supply.fresh_mass_flow_kg_per_s,
                (1.0 - fraction) * mass,
                1e-14,
            );
            close(supply.return_mass_flow_kg_per_s, fraction * mass, 1e-14);
            close(stream.mass_flow_kg_per_s, fraction * mass, 1e-14);
            // These W/K fields are explicitly reference-Cp diagnostics.
            close(
                supply.fresh_capacity_w_per_k,
                1400.0 * supply.fresh_mass_flow_kg_per_s,
                1e-10,
            );
            close(
                stream.capacity_w_per_k,
                1400.0 * stream.mass_flow_kg_per_s,
                1e-10,
            );
        }
    });
}

#[test]
fn changing_the_caloric_datum_preserves_the_complete_return_solution() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let run = |datum| {
            network(
                cx,
                &flow,
                source(-900.0, 6.0, datum),
                vec![heated("wall", 1500.0)],
                &[inlet(0, 300.0)],
            )
            .with_recirculation(cx, vec![link(0, 1, 0.75)], feedback_config())
            .unwrap()
            .march(cx, &[650.0])
            .unwrap()
        };
        let a = run(275.0);
        let b = run(650.0);
        for (left, right) in a.node_temperatures_k.iter().zip(&b.node_temperatures_k) {
            close(left.unwrap(), right.unwrap(), 2e-10);
        }
        for (&left, &right) in a
            .reference_temperatures_k
            .iter()
            .zip(&b.reference_temperatures_k)
        {
            close(left, right, 2e-10);
        }
        close(a.wall_heat_rate_w, b.wall_heat_rate_w, 2e-7);
        close(a.external_heat_gain_w, b.external_heat_gain_w, 2e-7);
        close(
            a.recirculation.as_ref().unwrap().external_heat_gain_w,
            b.recirculation.as_ref().unwrap().external_heat_gain_w,
            2e-7,
        );
        assert_eq!(a.enthalpy_reference_k, b.enthalpy_reference_k);
    });
}

#[test]
fn changing_only_cp_away_from_the_reference_changes_return_cooling_capacity() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let run = |intercept, slope| {
            network(
                cx,
                &flow,
                source(intercept, slope, 300.0),
                vec![heated("wall", 1500.0)],
                &[inlet(0, 300.0)],
            )
            .with_recirculation(cx, vec![link(0, 1, 0.75)], feedback_config())
            .unwrap()
            .march(cx, &[650.0])
            .unwrap()
        };
        // Both laws have Cp(300 K) = 900 J/(kg K): identical hydraulic,
        // density and reference-capacity admission, different actual storage.
        let variable = run(-900.0, 6.0);
        let frozen = run(900.0, 0.0);
        assert!(
            frozen.node_temperatures_k[1].unwrap() - variable.node_temperatures_k[1].unwrap()
                > 30.0
        );
        assert!(variable.wall_heat_rate_w - frozen.wall_heat_rate_w > 40_000.0);
        // Both laws reach the same reference Cp (900) through different arithmetic, so the two
        // reference capacities agree to a rounding, not bit for bit (225 and 225.00000000000003):
        // 2 eps relative.
        let capacities = [
            variable.recirculation.as_ref().unwrap().supplies[0].fresh_capacity_w_per_k,
            frozen.recirculation.as_ref().unwrap().supplies[0].fresh_capacity_w_per_k,
        ];
        assert!(
            (capacities[0] - capacities[1]).abs()
                <= 2.0 * f64::EPSILON * capacities[0].abs().max(capacities[1].abs())
        );
    });
}

#[test]
fn unsupported_return_states_and_exhausted_iteration_budgets_refuse() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let build = |conductance, config| {
            network(
                cx,
                &flow,
                source(-900.0, 6.0, 300.0),
                vec![heated("wall", conductance)],
                &[inlet(0, 300.0)],
            )
            .with_recirculation(cx, vec![link(0, 1, 0.75)], config)
            .unwrap()
        };
        for wall in [100.0, 1000.0] {
            assert!(matches!(
                build(100_000.0, feedback_config()).march(cx, &[wall]),
                Err(TransportError::Airflow(AirflowError::CaloricSource(_)))
            ));
        }
        let short = CaloricRecirculationConfig {
            max_iterations: 1,
            ..feedback_config()
        };
        match build(1500.0, short).march(cx, &[650.0]) {
            Err(TransportError::CaloricRecirculationNotConverged {
                iterations,
                max_mixing_residual_k,
                mixing_heat_residual_w,
            }) => {
                assert_eq!(iterations, 1);
                assert!(max_mixing_residual_k > short.temperature_tolerance_k);
                assert!(mixing_heat_residual_w > budgets().absolute_heat_tolerance_w);
            }
            other => panic!("a supported loop must exhaust its nonlinear budget: {other:?}"),
        }
    });
}

#[test]
fn the_shared_solid_fixed_point_closes_fresh_exhaust_enthalpy() {
    with_cx(|cx| {
        let flow = single_flow(cx);
        let conductances = [1500.0];
        let fraction = 0.65;
        let power = 60_000.0;
        let mass = flow.branches[0].flow.value();
        let net = network(
            cx,
            &flow,
            source(800.0, 2.0, 650.0),
            coupled_models(&conductances),
            &[inlet(0, 300.0)],
        )
        .with_recirculation(cx, vec![link(0, 1, fraction)], feedback_config())
        .unwrap();
        let mut calls = 0;
        let solved = net
            .solve(cx, &coupling(), |_, references| {
                calls += 1;
                Ok(solid(references, &conductances, power))
            })
            .unwrap();
        let exhaust_enthalpy = power / (mass * (1.0 - fraction));
        close(
            solved.transport.node_temperatures_k[1].unwrap(),
            after_heat(300.0, exhaust_enthalpy),
            2e-8,
        );
        close(
            solved.transport.node_temperatures_k[0].unwrap(),
            after_heat(300.0, fraction * exhaust_enthalpy),
            2e-8,
        );
        close(solved.transport.wall_heat_rate_w, power, 3e-6);
        close(
            solved
                .transport
                .recirculation
                .as_ref()
                .unwrap()
                .external_heat_gain_w,
            power,
            3e-6,
        );
        assert_eq!(calls, solved.iterations);
        assert_eq!(
            solved.solid[0].mean_reference_temperature_k,
            Some(solved.reference_temperatures_k[0])
        );
        assert!(solved.region_balances[0].imbalance_w.abs() <= coupling().balance_tolerance_w);
    });
}

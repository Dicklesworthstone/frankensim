//! Surface radiation through the same heat ports and step driver as other thermal devices.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortTangent, Receipt};
use fs_system::elements::{
    ConductanceMatrix, ElementError, HeatCapacity, Heater, SourceAccount, Waveform,
};
use fs_system::{
    AnalysisError, BlockSolverPolicy, MonolithicConfig, NodeId, StateMap, StationaryConfig, System,
    SystemBuilder,
};

const SIGMA: f64 = 5.670_374_419e-8;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 52,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn close(actual: f64, expected: f64) {
    let tolerance = 128.0 * f64::EPSILON * actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} against {expected:e}, allowance {tolerance:e}"
    );
}

fn step(element: &ConductanceMatrix, kelvin: [f64; 2], dt_s: f64) -> ([f64; 2], Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s,
        inputs: &kelvin,
    };
    with_cx(|cx| match check_component(element, &case, cx).unwrap() {
        Observed::Stepped {
            receipt, outputs, ..
        } => ([outputs[0], outputs[1]], receipt),
        other => panic!("radiation did not take the step: {other:?}"),
    })
}

/// Independent radiosity reduction and direct fourth powers; no production secant.
fn coefficient(area_m2: f64, emissivities: [f64; 2]) -> f64 {
    SIGMA * area_m2 / (1.0 / emissivities[0] + 1.0 / emissivities[1] - 1.0)
}

#[test]
fn surfaces_match_stefan_boltzmann_reversal_and_equilibrium() {
    for emissivities in [[1.0, 1.0], [0.8, 0.6], [1.0e-200, 0.5]] {
        let area = 0.125;
        let element = ConductanceMatrix::gray_radiation(area, emissivities).unwrap();
        assert_eq!(element.widths(), [1, 1]);
        assert_eq!(element.class(), Class::Smooth);
        for temperatures in [[900.0_f64, 300.0], [300.0, 900.0], [450.0, 450.0]] {
            let dt = 0.25;
            let (loads, receipt) = step(&element, temperatures, dt);
            let expected = coefficient(area, emissivities)
                * (temperatures[0].powi(4) - temperatures[1].powi(4));
            close(loads[0], expected);
            close(loads[1], -expected);
            close(loads.iter().sum(), 0.0);
            close(
                receipt.entropy_generated_j_per_k,
                dt * expected * (temperatures[1].recip() - temperatures[0].recip()),
            );
            assert!(receipt.entropy_generated_j_per_k >= 0.0);
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
            assert_eq!(receipt.numerical_loss_j, 0.0);
        }
    }
}

#[test]
fn tangents_and_transposes_follow_the_fourth_power_law() {
    let area = 0.4;
    let emissivities = [0.75, 0.625];
    let element = ConductanceMatrix::gray_radiation(area, emissivities).unwrap();
    with_cx(|cx| {
        for temperatures in [[300.0_f64, 900.0], [900.0, 300.0], [450.0, 450.0]] {
            let case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.1,
                inputs: &temperatures,
            };
            check_port_tangent(&element, &case, 1.0e-8, 1.0e-8, cx).unwrap();
            let slopes = [
                4.0 * coefficient(area, emissivities) * temperatures[0].powi(3),
                -4.0 * coefficient(area, emissivities) * temperatures[1].powi(3),
            ];
            let direction = [-0.4, 0.7];
            let mut output = [0.0; 2];
            element
                .jvp(&(), 0.1, &temperatures, &direction, &mut output, cx)
                .unwrap();
            let expected = slopes[0] * direction[0] + slopes[1] * direction[1];
            close(output[0], expected);
            close(output[1], -expected);
            let covector = [1.3, -0.6];
            element
                .vjp(&(), 0.1, &temperatures, &covector, &mut output, cx)
                .unwrap();
            close(output[0], slopes[0] * (covector[0] - covector[1]));
            close(output[1], slopes[1] * (covector[0] - covector[1]));
        }
        // Observe each derivative separately: a hot-dominated directional sum
        // would hide cancellation of the colder surface's small derivative.
        for temperatures in [[300.0_f64, 0.01], [1.0e8, 1.0]] {
            let slopes = [
                4.0 * coefficient(area, emissivities) * temperatures[0].powi(3),
                -4.0 * coefficient(area, emissivities) * temperatures[1].powi(3),
            ];
            let mut output = [0.0; 2];
            for column in 0..2 {
                let mut direction = [0.0; 2];
                direction[column] = 1.0;
                element
                    .jvp(&(), 0.1, &temperatures, &direction, &mut output, cx)
                    .unwrap();
                close(output[0], slopes[column]);
                close(output[1], -slopes[column]);
            }
            element
                .vjp(&(), 0.1, &temperatures, &[1.0, 0.0], &mut output, cx)
                .unwrap();
            close(output[0], slopes[0]);
            close(output[1], slopes[1]);
        }
    });
    // At the reference temperature the dimensionless factor is exactly one.
    let mut slope = [0.0; 2];
    with_cx(|cx| {
        element
            .jvp(&(), 0.1, &[1.0, 1.0], &[1.0, 0.0], &mut slope, cx)
            .unwrap();
    });
    assert_eq!(
        slope[0].to_bits(),
        element.points()[0].conductance_w_k.to_bits()
    );
}

#[test]
fn one_ulp_temperature_difference_still_transfers_heat() {
    let area = 0.5;
    let emissivities = [0.7, 0.9];
    let element = ConductanceMatrix::gray_radiation(area, emissivities).unwrap();
    let cold = 500.0_f64;
    let hot = f64::from_bits(cold.to_bits() + 1);
    let delta = hot - cold;
    // Independent binomial expansion about the colder surface, avoiding the
    // cancellation the direct fourth-power oracle would suffer here.
    let expected = coefficient(area, emissivities)
        * delta
        * (4.0 * cold.powi(3) + delta * (6.0 * cold * cold + delta * (4.0 * cold + delta)));
    let (forward, _) = step(&element, [hot, cold], 0.1);
    let (reverse, _) = step(&element, [cold, hot], 0.1);
    assert!(forward[0] > 0.0);
    close(forward[0], expected);
    close(reverse[0], -expected);
}

#[test]
fn parameters_and_unrepresentable_live_rates_refuse() {
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ConductanceMatrix::gray_radiation(area, [0.8, 0.6]),
            Err(ElementError::Parameter { .. })
        ));
    }
    for value in [0.0, -0.1, 1.001, f64::NAN, f64::INFINITY] {
        for emissivities in [[value, 0.8], [0.8, value]] {
            assert!(matches!(
                ConductanceMatrix::gray_radiation(1.0, emissivities),
                Err(ElementError::Parameter { .. })
            ));
        }
    }
    assert!(ConductanceMatrix::gray_radiation(f64::MIN_POSITIVE, [1.0e-200; 2]).is_err());
    let element = ConductanceMatrix::gray_radiation(1.0, [1.0; 2]).unwrap();
    with_cx(|cx| {
        let mut out = [0.0; 2];
        assert!(matches!(
            element.propose(&(), 0.1, &[f64::MAX, 300.0], &mut (), &mut out, cx),
            Err(ElementError::Law { .. })
        ));
        assert!(matches!(
            element.propose(&(), 0.1, &[0.0, 300.0], &mut (), &mut out, cx),
            Err(ElementError::Temperature { .. })
        ));
        assert!(matches!(
            element.propose(&(), 0.1, &[f64::MIN_POSITIVE; 2], &mut (), &mut out, cx),
            Err(ElementError::Law { .. })
        ));
    });
}

#[test]
fn identity_keeps_authored_area_and_ordered_emissivities() {
    let words = |area, emissivities| {
        let element = ConductanceMatrix::gray_radiation(area, emissivities).unwrap();
        let mut words = Vec::new();
        element.identity_words(&mut words);
        words
    };
    let identity = words(0.5, [0.8, 0.6]);
    assert_eq!(identity, words(0.5, [0.8, 0.6]));
    assert_ne!(identity, words(0.6, [0.8, 0.6]));
    assert_ne!(identity, words(0.5, [0.6, 0.8]));
    let a = ConductanceMatrix::gray_radiation(0.5, [0.8, 0.6]).unwrap();
    let b = ConductanceMatrix::gray_radiation(0.5, [0.6, 0.8]).unwrap();
    assert_eq!(step(&a, [600.0, 300.0], 0.1), step(&b, [600.0, 300.0], 0.1));
}

#[test]
fn finite_bodies_radiate_on_the_generic_graph_against_an_independent_midpoint_solve() {
    for policy in [
        BlockSolverPolicy::InterfaceNewton,
        BlockSolverPolicy::MonolithicNewtonKrylov(MonolithicConfig::default()),
    ] {
        check_radiating_pair(policy);
    }
}

fn check_radiating_pair(policy: BlockSolverPolicy) {
    let capacity = [40.0, 80.0];
    let mut reference = [600.0_f64, 300.0];
    let initial_energy = capacity[0] * reference[0] + capacity[1] * reference[1];
    let area = 0.2;
    let emissivities = [0.8, 0.6];
    let mut builder = SystemBuilder::new();
    let a = builder
        .add_with_tangent_and_stage_residual(
            "a",
            HeatCapacity::new(capacity[0]).unwrap(),
            reference[0],
        )
        .unwrap();
    let b = builder
        .add_with_tangent_and_stage_residual(
            "b",
            HeatCapacity::new(capacity[1]).unwrap(),
            reference[1],
        )
        .unwrap();
    let radiation = builder
        .add_with_tangent_and_stage_residual(
            "radiation",
            ConductanceMatrix::gray_radiation(area, emissivities).unwrap(),
            (),
        )
        .unwrap();
    let radiation_a = builder.port(radiation, "a").unwrap();
    let radiation_b = builder.port(radiation, "b").unwrap();
    builder
        .bond(builder.port(a, "heat").unwrap(), radiation_a)
        .unwrap();
    builder
        .bond(builder.port(b, "heat").unwrap(), radiation_b)
        .unwrap();
    let mut system = builder.seal().unwrap();
    system.set_solver_policy(policy);
    assert!(system.guess(radiation_a, &[reference[0]]));
    assert!(system.guess(radiation_b, &[reference[1]]));
    let dt = 0.1;
    with_cx(|cx| {
        for _ in 0..64 {
            // Solve q = sigma A epsilon_eff (mean(Ta)^4 - mean(Tb)^4)
            // by bisection, using energy conservation for the two new temperatures.
            let coefficient = coefficient(area, emissivities);
            let mut low = 0.0;
            let mut high = coefficient * (reference[0].powi(4) - reference[1].powi(4));
            for _ in 0..80 {
                let q = 0.5 * (low + high);
                let means = [
                    reference[0] - 0.5 * dt * q / capacity[0],
                    reference[1] + 0.5 * dt * q / capacity[1],
                ];
                let residual = q - coefficient * (means[0].powi(4) - means[1].powi(4));
                if residual > 0.0 {
                    high = q;
                } else {
                    low = q;
                }
            }
            let exchanged = dt * 0.5 * (low + high);
            reference[0] -= exchanged / capacity[0];
            reference[1] += exchanged / capacity[1];
            system.step(dt, cx).unwrap();
            for (node, expected) in [(a, reference[0]), (b, reference[1])] {
                let actual = *system.state::<f64>(node).unwrap();
                assert!(
                    (actual - expected).abs() <= 1.0e-8,
                    "{actual} against {expected} K"
                );
            }
        }
    });
    let final_a = *system.state::<f64>(a).unwrap();
    let final_b = *system.state::<f64>(b).unwrap();
    assert!(400.0 < final_a && final_a < 600.0);
    assert!(300.0 < final_b && final_b < 400.0);
    // Forward roundoff allowance for 64 pairs of state updates and their sums.
    let allowance = 64.0 * 32.0 * f64::EPSILON * initial_energy;
    assert!((capacity[0] * final_a + capacity[1] * final_b - initial_energy).abs() <= allowance);
    assert!(system.ledger().energy_residual_j().abs() <= allowance);
    assert!(system.ledger().entropy_generated_j_per_k() > 0.0);
}

struct TemperatureState(NodeId);

impl StateMap for TemperatureState {
    fn dim(&self) -> usize {
        1
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system.state::<f64>(self.0).ok_or(AnalysisError::Invalid {
            what: "thermal body state",
        })?;
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        system.set_accepted_state(self.0, state[0])?;
        Ok(())
    }
}

#[test]
fn heated_body_reaches_the_stefan_boltzmann_stationary_temperature() {
    let (ambient, power_w, area, emissivities) = (300.0_f64, 75.0, 0.2, [0.8, 0.6]);
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent_and_stage_residual("body", HeatCapacity::new(50.0).unwrap(), 350.0)
        .unwrap();
    let radiation = builder
        .add_with_tangent_and_stage_residual(
            "radiation",
            ConductanceMatrix::gray_radiation(area, emissivities).unwrap(),
            (),
        )
        .unwrap();
    let heater = builder
        .add_with_tangent_and_stage_residual(
            "heater",
            Heater::new(Waveform::constant(power_w).unwrap()).unwrap(),
            SourceAccount::default(),
        )
        .unwrap();
    let (body_heat, radiation_a, radiation_b, heater_heat) = (
        builder.port(body, "heat").unwrap(),
        builder.port(radiation, "a").unwrap(),
        builder.port(radiation, "b").unwrap(),
        builder.port(heater, "heat").unwrap(),
    );
    builder
        .common_effort(body_heat, &[radiation_a, heater_heat])
        .unwrap();
    builder
        .boundary(radiation_b, Box::new(move |_, _, out| out[0] = ambient))
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(radiation_a, &[350.0]) && system.guess(heater_heat, &[350.0]));
    let expected = (ambient.powi(4) + power_w / coefficient(area, emissivities))
        .sqrt()
        .sqrt();
    with_cx(|cx| {
        let solution = system
            .solve_stationary(&TemperatureState(body), StationaryConfig::default(), cx)
            .unwrap();
        assert!(solution.residual_norm <= 1.0e-10);
        assert!(
            (solution.state[0] - expected).abs() < 1.0e-7,
            "stationary temperature {} K against {expected} K",
            solution.state[0]
        );
    });
    // A noncommitting analysis leaves the thermal state, source work and clock untouched.
    assert_eq!(*system.state::<f64>(body).unwrap(), 350.0);
    assert_eq!(
        *system.state::<SourceAccount>(heater).unwrap(),
        SourceAccount::default()
    );
    assert_eq!(system.time_s(), 0.0);
}

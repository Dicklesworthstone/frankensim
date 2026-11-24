//! Live passive dissipation on the shared component and graph paths.
#![allow(clippy::float_cmp)]

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortTangent, Proposal, Receipt};
use fs_system::SystemBuilder;
use fs_system::elements::{Damper, DissipationLaw, ElementError, HeatCapacity, Mass, Motion};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 11,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn close(actual: f64, expected: f64) {
    let tolerance = 32.0 * f64::EPSILON * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}

/// An authored analytic fixture, not a fitted material dataset.
#[derive(Debug)]
struct Authored {
    thermal_slope: f64,
    flow_slope: f64,
}

impl DissipationLaw for Authored {
    fn factor(&self, kelvin: f64, flow: f64) -> Result<(f64, [f64; 2]), String> {
        if !(295.0..=350.0).contains(&kelvin) {
            return Err("authored support is 295..=350 K".into());
        }
        let thermal = 1.0 + self.thermal_slope * (kelvin - 300.0);
        let flow_factor = 1.0 + self.flow_slope * flow.abs();
        let sign = if flow == 0.0 { 0.0 } else { flow.signum() };
        Ok((
            thermal * flow_factor,
            [
                self.thermal_slope * flow_factor,
                thermal * self.flow_slope * sign,
            ],
        ))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            u64::from_le_bytes(*b"testlaw1"),
            self.thermal_slope.to_bits(),
            self.flow_slope.to_bits(),
            295.0_f64.to_bits(),
            350.0_f64.to_bits(),
        ]);
    }
}

fn authored() -> Arc<dyn DissipationLaw> {
    Arc::new(Authored {
        thermal_slope: 0.2,
        flow_slope: 0.5,
    })
}

#[derive(Debug)]
struct Fixed {
    factor: f64,
    partials: [f64; 2],
}

impl DissipationLaw for Fixed {
    fn factor(&self, _: f64, _: f64) -> Result<(f64, [f64; 2]), String> {
        Ok((self.factor, self.partials))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            u64::from_le_bytes(*b"fixedlaw"),
            self.factor.to_bits(),
            self.partials[0].to_bits(),
            self.partials[1].to_bits(),
        ]);
    }
}

fn receipt(damper: &Damper, inputs: &[f64], outputs: &mut [f64], cx: &Cx<'_>) -> Receipt {
    match damper
        .propose(&(), 0.125, inputs, &mut (), outputs, cx)
        .unwrap()
    {
        Proposal::Step(receipt) => receipt,
        Proposal::Shorten { .. } => panic!("a stateless damper does not shorten a step"),
    }
}

#[test]
fn live_temperature_and_signed_flow_preserve_power_and_have_exact_tangents() {
    let coefficient = 2.5;
    let damper = Damper::new(coefficient).unwrap().following(authored());
    with_cx(|cx| {
        for inputs in [
            [3.0, -1.0, 300.0],
            [-2.0, 1.0, 315.0],
            [1.5, -0.25, 340.0],
            [1.25, 1.25, 300.0],
        ] {
            let case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.125,
                inputs: &inputs,
            };
            let Observed::Stepped {
                outputs, receipt, ..
            } = check_component(&damper, &case, cx).unwrap()
            else {
                panic!("a stateless damper does not shorten a step");
            };
            let rate = inputs[0] - inputs[1];
            let thermal = 1.0 + 0.2 * (inputs[2] - 300.0);
            let flow = 1.0 + 0.5 * rate.abs();
            // Expand the physical effort and heating independently of the
            // product-rule implementation in elements::dissipation.
            let effort = coefficient * thermal * rate * flow;
            let power = coefficient * thermal * rate * rate * flow;
            close(outputs[0], effort);
            close(outputs[1], -effort);
            close(outputs[2], -power);
            close(receipt.dissipated_j, 0.125 * power);
            close(receipt.entropy_generated_j_per_k, 0.125 * power / inputs[2]);

            let effort_by_flow = coefficient * thermal * (1.0 + rate.abs());
            let effort_by_temperature = coefficient * 0.2 * rate * flow;
            let power_by_flow = coefficient * thermal * (2.0 * rate + 1.5 * rate * rate.abs());
            let power_by_temperature = coefficient * 0.2 * rate * rate * flow;
            let direction = [0.4, -0.7, 0.25];
            let mut jvp = [0.0; 3];
            damper
                .jvp(&(), case.dt_s, &inputs, &direction, &mut jvp, cx)
                .unwrap();
            let d_effort = effort_by_flow * 1.1 + effort_by_temperature * 0.25;
            let d_power = power_by_flow * 1.1 + power_by_temperature * 0.25;
            close(jvp[0], d_effort);
            close(jvp[1], -d_effort);
            close(jvp[2], -d_power);
            let mut vjp = [0.0; 3];
            damper
                .vjp(&(), case.dt_s, &inputs, &[1.5, -0.25, -0.5], &mut vjp, cx)
                .unwrap();
            let along = 1.75 * effort_by_flow + 0.5 * power_by_flow;
            close(vjp[0], along);
            close(vjp[1], -along);
            close(
                vjp[2],
                1.75 * effort_by_temperature + 0.5 * power_by_temperature,
            );

            // At zero flow, F is C1 but has a curvature corner. Its exact
            // derivative is checked above; the kit's fixed finite difference
            // step only has first-order truncation at that corner.
            if rate != 0.0 {
                check_port_tangent(&damper, &case, 2.0e-6, 1.0e-8, cx).unwrap();
            }
        }
    });
}

#[test]
fn an_absent_law_keeps_constant_bits_and_law_identity_is_durable() {
    let plain = Damper::new(2.5).unwrap();
    let unit = plain.clone().following(Arc::new(Fixed {
        factor: 1.0,
        partials: [0.0; 2],
    }));
    with_cx(|cx| {
        for inputs in [[3.0, -1.0, 300.0], [-2.0, 1.0, 315.0], [1.25, 1.25, 300.0]] {
            let mut expected = [0.0; 3];
            let mut actual = [0.0; 3];
            let a = receipt(&plain, &inputs, &mut expected, cx);
            let b = receipt(&unit, &inputs, &mut actual, cx);
            assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
            assert_eq!(a.dissipated_j.to_bits(), b.dissipated_j.to_bits());
            assert_eq!(
                a.entropy_generated_j_per_k.to_bits(),
                b.entropy_generated_j_per_k.to_bits()
            );
            assert_eq!(
                a.energy_tolerance_j.to_bits(),
                b.energy_tolerance_j.to_bits()
            );
        }
    });
    let identity = |element: &Damper| {
        let mut words = Vec::new();
        element.identity_words(&mut words);
        words
    };
    assert_eq!(identity(&plain), [2.5_f64.to_bits()]);
    assert_ne!(identity(&plain), identity(&unit));
    let other = plain.clone().following(Arc::new(Authored {
        thermal_slope: 0.2,
        flow_slope: 0.25,
    }));
    assert_ne!(identity(&plain.following(authored())), identity(&other));
}

#[derive(Debug)]
struct NegativeDifferential;

impl DissipationLaw for NegativeDifferential {
    fn factor(&self, _: f64, flow: f64) -> Result<(f64, [f64; 2]), String> {
        let denominator = 1.0 + flow * flow;
        Ok((
            1.0 / denominator,
            [0.0, -2.0 * flow / (denominator * denominator)],
        ))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(u64::from_le_bytes(*b"passive1"));
    }
}

#[test]
fn passivity_allows_negative_differential_resistance_and_explicit_zero_loss() {
    let damper = Damper::new(1.0)
        .unwrap()
        .following(Arc::new(NegativeDifferential));
    assert_eq!(damper.class(), Class::Smooth);
    with_cx(|cx| {
        let mut outputs = [0.0; 3];
        let observed = receipt(&damper, &[2.0, 0.0, 300.0], &mut outputs, cx);
        close(outputs[0], 0.4);
        close(outputs[2], -0.8);
        assert!(observed.entropy_generated_j_per_k > 0.0);
        damper
            .jvp(
                &(),
                0.125,
                &[2.0, 0.0, 300.0],
                &[1.0, 0.0, 0.0],
                &mut outputs,
                cx,
            )
            .unwrap();
        close(outputs[0], -0.12);
        close(outputs[1], 0.12);
        close(outputs[2], -0.16);

        let zero = Damper::new(2.0).unwrap().following(Arc::new(Fixed {
            factor: 0.0,
            partials: [0.0; 2],
        }));
        let observed = receipt(&zero, &[3.0, -1.0, 300.0], &mut outputs, cx);
        assert!(outputs.iter().all(|value| *value == 0.0));
        assert_eq!(observed.dissipated_j, 0.0);
        zero.jvp(
            &(),
            0.125,
            &[3.0, -1.0, 300.0],
            &[1.0, -1.0, 1.0],
            &mut outputs,
            cx,
        )
        .unwrap();
        assert!(outputs.iter().all(|value| *value == 0.0));
    });
}

#[test]
fn unsupported_or_unrepresentable_laws_refuse_before_writing_outputs() {
    with_cx(|cx| {
        for (coefficient, factor, partials) in [
            (1.0, -1.0, [0.0, 0.0]),
            (1.0, f64::NAN, [0.0, 0.0]),
            (1.0, f64::INFINITY, [0.0, 0.0]),
            (1.0, 1.0, [f64::NAN, 0.0]),
            (1.0, 1.0, [0.0, f64::INFINITY]),
            (1.0e-300, 1.0e-300, [0.0, 0.0]),
            (f64::MAX, 1.0, [0.0, 0.0]),
        ] {
            let damper = Damper::new(coefficient)
                .unwrap()
                .following(Arc::new(Fixed { factor, partials }));
            let inputs = [2.0, 0.0, 300.0];
            let mut outputs = [11.0, 12.0, 13.0];
            assert!(matches!(
                damper.propose(&(), 0.125, &inputs, &mut (), &mut outputs, cx),
                Err(ElementError::Law { .. })
            ));
            assert_eq!(outputs, [11.0, 12.0, 13.0]);
            assert!(matches!(
                damper.jvp(&(), 0.125, &inputs, &[1.0; 3], &mut outputs, cx),
                Err(ElementError::Law { .. })
            ));
            assert_eq!(outputs, [11.0, 12.0, 13.0]);
            assert!(matches!(
                damper.vjp(&(), 0.125, &inputs, &[1.0; 3], &mut outputs, cx),
                Err(ElementError::Law { .. })
            ));
            assert_eq!(outputs, [11.0, 12.0, 13.0]);
        }
        let damper = Damper::new(1.0).unwrap().following(authored());
        let mut outputs = [11.0, 12.0, 13.0];
        let error = damper
            .propose(&(), 0.125, &[2.0, 0.0, 351.0], &mut (), &mut outputs, cx)
            .unwrap_err();
        assert!(
            matches!(error, ElementError::Law { reason } if reason == "authored support is 295..=350 K")
        );
        assert_eq!(outputs, [11.0, 12.0, 13.0]);
    });
}

#[test]
fn self_heating_changes_the_live_resistance_on_the_shared_midpoint_graph() {
    let (mass_kg, capacity, coefficient, dt_s) = (2.0, 2.0, 2.0, 0.05);
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 4.0,
            },
        )
        .unwrap();
    let damper = builder
        .add_with_tangent(
            "damper",
            Damper::new(coefficient).unwrap().following(authored()),
            (),
        )
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(capacity).unwrap(), 300.0)
        .unwrap();
    let (p, a) = (
        builder.port(mass, "p").unwrap(),
        builder.port(damper, "a").unwrap(),
    );
    builder.common_flow(p, &[a]).unwrap();
    let b = builder.port(damper, "b").unwrap();
    builder
        .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    let body_heat = builder.port(body, "heat").unwrap();
    let damper_heat = builder.port(damper, "heat").unwrap();
    builder.common_effort(body_heat, &[damper_heat]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(damper_heat, &[300.0]));

    let (mut velocity, mut kelvin) = (4.0, 300.0);
    with_cx(|cx| {
        for _ in 0..32 {
            // Independently eliminate T_new using exact kinetic/thermal
            // energy balance, then bracket the scalar midpoint momentum
            // equation. No graph, Damper, or callback code enters this oracle.
            let mut lower = 0.0;
            let mut upper = velocity;
            for _ in 0..80 {
                let next = 0.5 * (lower + upper);
                let mean = 0.5 * (velocity + next);
                let mean_temperature =
                    kelvin + 0.25 * mass_kg * (velocity * velocity - next * next) / capacity;
                let resistance = coefficient
                    * (1.0 + 0.2 * (mean_temperature - 300.0))
                    * (1.0 + 0.5 * mean.abs());
                let residual = next - velocity + dt_s * resistance * mean / mass_kg;
                if residual > 0.0 {
                    upper = next;
                } else {
                    lower = next;
                }
            }
            let next = 0.5 * (lower + upper);
            kelvin += 0.5 * mass_kg * (velocity * velocity - next * next) / capacity;
            velocity = next;
            system.step(dt_s, cx).unwrap();
            let actual_velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let actual_kelvin = *system.state::<f64>(body).unwrap();
            assert!(
                (actual_velocity - velocity).abs() < 1.0e-8,
                "{actual_velocity} vs {velocity}"
            );
            assert!(
                (actual_kelvin - kelvin).abs() < 1.0e-8,
                "{actual_kelvin} vs {kelvin}"
            );
            let energy = 0.5 * mass_kg * actual_velocity * actual_velocity
                + capacity * (actual_kelvin - 300.0);
            assert!((energy - 16.0).abs() < 1.0e-8, "energy = {energy}");
        }
    });
    assert!(velocity < 0.1);
    assert!(kelvin > 307.9);
    assert!(system.ledger().energy_residual_j().abs() < 1.0e-8);
    assert!(system.ledger().entropy_generated_j_per_k() > 0.0);
}

//! G0/G1/G3: original layer laws, analytic tangents and coupled thermal storage.

#![cfg(feature = "component")]
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::thermal::conduction::component::LayerConductor;
use fs_material::thermal::conduction::{ConductivityCurve, LayerConductionConfig, LayerGeometry};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortTangent, Receipt};
use fs_system::SystemBuilder;
use fs_thermal::enthalpy::DepositionError;

fn in_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 80,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 2e-13,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn response(
    conductor: &LayerConductor,
    temperatures: &[f64; 2],
    cx: &Cx<'_>,
) -> (Vec<f64>, Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: 0.25,
        inputs: temperatures,
    };
    let Observed::Stepped {
        outputs, receipt, ..
    } = check_component(conductor, &case, cx).unwrap()
    else {
        panic!("a quasisteady conductor has no clock to shorten");
    };
    (outputs, receipt)
}

#[test]
fn constant_layers_reproduce_series_resistance_reversal_and_equilibrium() {
    let curves = [
        ConductivityCurve::new(vec![[250.0, 5.0], [700.0, 5.0]]).unwrap(),
        ConductivityCurve::new(vec![[250.0, 10.0], [700.0, 10.0]]).unwrap(),
    ];
    let geometry = LayerGeometry::planar([0.03, 0.04], 0.2).unwrap();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for contact in [0.0, 0.08] {
            let conductor =
                LayerConductor::new(curves.clone(), geometry, contact, config()).unwrap();
            assert_eq!(conductor.class(), Class::Smooth);
            assert_eq!(conductor.config(), config());
            assert_eq!(conductor.geometry(), geometry);
            assert_eq!(conductor.conductivity_curves(), &curves);
            let conductance = 0.2 / (0.03 / 5.0 + contact + 0.04 / 10.0);
            for temperatures in [[600.0, 300.0], [300.0, 600.0], [400.0, 400.0]] {
                let (output, receipt) = response(&conductor, &temperatures, cx);
                let heat_w = conductance * (temperatures[0] - temperatures[1]);
                close(output[0], heat_w, 5e-13, 1e-11);
                assert_eq!(output[1], -output[0]);
                assert_eq!(receipt.stored_change_j, 0.0);
                assert_eq!(receipt.dissipated_j, 0.0);
                assert_eq!(receipt.numerical_loss_j, 0.0);
                close(
                    receipt.entropy_generated_j_per_k,
                    0.25 * heat_w * (1.0 / temperatures[1] - 1.0 / temperatures[0]),
                    5e-13,
                    1e-13,
                );
                check_port_tangent(
                    &conductor,
                    &Case {
                        state: &(),
                        stale: &(),
                        dt_s: 0.25,
                        inputs: &temperatures,
                    },
                    1e-7,
                    1e-9,
                    cx,
                )
                .unwrap();
            }
        }
    });
}

#[test]
fn nonlinear_flux_and_nonsymmetric_tangent_match_independent_kirchhoff_integral() {
    // k(T)=0.4+0.002 T. Both layers have the same chart, so
    // Q=A/L * integral(Tb..Ta, k(T) dT) needs no interface solve.
    let curve = ConductivityCurve::new(vec![[250.0, 0.9], [700.0, 1.8]]).unwrap();
    let conductor = LayerConductor::new(
        [curve.clone(), curve],
        LayerGeometry::planar([0.2, 0.3], 0.05).unwrap(),
        0.0,
        config(),
    )
    .unwrap();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for temperatures in [
            [600.0, 300.0],
            [300.0, 600.0],
            [450.0, 450.0],
            [430.0, 490.0],
        ] {
            let [a, b] = temperatures;
            let expected = 0.1 * (a - b) * (0.4 + 0.001 * (a + b));
            let (output, receipt) = response(&conductor, &temperatures, cx);
            close(output[0], expected, 5e-13, 2e-11);
            assert!(receipt.residual <= receipt.tolerance);
            let jacobian = [0.1 * (0.4 + 0.002 * a), -0.1 * (0.4 + 0.002 * b)];
            let mut forward = [f64::NAN; 2];
            conductor
                .jvp(&(), 0.25, &temperatures, &[2.0, -3.0], &mut forward, cx)
                .unwrap();
            close(
                forward[0],
                2.0 * jacobian[0] - 3.0 * jacobian[1],
                5e-13,
                1e-12,
            );
            assert_eq!(forward[1], -forward[0]);
            let mut transpose = [f64::NAN; 2];
            conductor
                .vjp(&(), 0.25, &temperatures, &[2.0, -3.0], &mut transpose, cx)
                .unwrap();
            close(transpose[0], 5.0 * jacobian[0], 5e-13, 1e-12);
            close(transpose[1], 5.0 * jacobian[1], 5e-13, 1e-12);
            check_port_tangent(
                &conductor,
                &Case {
                    state: &(),
                    stale: &(),
                    dt_s: 0.25,
                    inputs: &temperatures,
                },
                1e-7,
                1e-9,
                cx,
            )
            .unwrap();
        }
    });
}

#[test]
fn each_material_keeps_its_own_support_across_a_contact_jump() {
    let curves = [
        ConductivityCurve::new(vec![[500.0, 5.0], [700.0, 5.0]]).unwrap(),
        ConductivityCurve::new(vec![[250.0, 10.0], [400.0, 10.0]]).unwrap(),
    ];
    let geometry = LayerGeometry::planar([0.03, 0.04], 0.2).unwrap();
    let conductor = LayerConductor::new(curves.clone(), geometry, 0.08, config()).unwrap();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let (output, _) = response(&conductor, &[600.0, 300.0], cx);
        close(output[0], 0.2 * 300.0 / 0.09, 5e-13, 1e-11);
        // Bulk temperatures are supported, but perfect contact cannot
        // choose a common face temperature in these disjoint source spans.
        let unsupported = LayerConductor::new(curves, geometry, 0.0, config()).unwrap();
        let mut output = [7.0, 8.0];
        assert_eq!(
            unsupported.propose(&(), 0.25, &[600.0, 300.0], &mut (), &mut output, cx),
            Err(DepositionError::MaterialDomain)
        );
        assert_eq!(output, [7.0, 8.0]);
    });
}

#[test]
fn invalid_inputs_support_limits_and_cancellation_leave_outputs_untouched() {
    let curve = ConductivityCurve::new(vec![[250.0, 0.9], [700.0, 1.8]]).unwrap();
    let curves = [curve.clone(), curve];
    let geometry = LayerGeometry::planar([0.2, 0.3], 0.05).unwrap();
    assert!(
        LayerConductor::new(
            curves.clone(),
            geometry,
            0.0,
            LayerConductionConfig {
                max_iterations: 0,
                ..config()
            }
        )
        .is_err()
    );
    assert!(
        LayerConductor::new(
            curves.clone(),
            geometry,
            0.0,
            LayerConductionConfig {
                relative_heat_tolerance: f64::NAN,
                ..config()
            }
        )
        .is_err()
    );
    assert!(LayerConductor::new(curves.clone(), geometry, -1.0, config()).is_err());
    let conductor = LayerConductor::new(curves, geometry, 0.0, config()).unwrap();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let attempt = |dt, inputs: &[f64]| {
            let mut output = [7.0, 8.0];
            assert!(
                conductor
                    .propose(&(), dt, inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0]);
            assert!(
                conductor
                    .jvp(&(), dt, inputs, &[1.0, 1.0], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0]);
            assert!(
                conductor
                    .vjp(&(), dt, inputs, &[1.0, 1.0], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0]);
        };
        for inputs in [
            vec![],
            vec![300.0],
            vec![300.0; 3],
            vec![f64::NAN, 300.0],
            vec![249.0, 300.0],
            vec![300.0, 701.0],
        ] {
            attempt(0.25, &inputs);
        }
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            attempt(dt, &[600.0, 300.0]);
        }
        let mut short = [7.0];
        assert!(
            conductor
                .propose(&(), 0.25, &[600.0, 300.0], &mut (), &mut short, cx)
                .is_err()
        );
        assert_eq!(short, [7.0]);
        let mut output = [7.0, 8.0];
        assert!(
            conductor
                .jvp(&(), 0.25, &[600.0, 300.0], &[1.0], &mut output, cx)
                .is_err()
        );
        assert!(
            conductor
                .vjp(
                    &(),
                    0.25,
                    &[600.0, 300.0],
                    &[f64::NAN, 1.0],
                    &mut output,
                    cx
                )
                .is_err()
        );
        assert_eq!(output, [7.0, 8.0]);
        cx.cancel_gate().request();
        assert_eq!(
            conductor.propose(&(), 0.25, &[600.0, 300.0], &mut (), &mut output, cx),
            Err(DepositionError::Cancelled)
        );
        attempt(0.25, &[600.0, 300.0]);
    });
}

fn body(capacity: f64, initial: f64) -> (UniformEnthalpyBody, EquilibriumPhaseState) {
    let knot = |temperature| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: capacity * (temperature - 250.0),
        temperature_k: temperature,
        liquid_mass_fraction: 0.0,
        bulk_density_kg_m3: 1000.0,
    };
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x39; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(250.0), knot(700.0)],
    )
    .unwrap();
    let state = chart
        .state_at_specific_enthalpy(capacity * (initial - 250.0))
        .unwrap();
    (UniformEnthalpyBody::new(chart, 1.0).unwrap(), state)
}

fn mean_temperature(old: f64, next: f64) -> f64 {
    if old == next {
        old
    } else {
        (next - old) / ((next - old) / old).ln_1p()
    }
}

// Independent elimination of two caloric stores, logarithmic-mean
// temperatures and the affine-k integral; no production layer solve.
fn coupled_oracle(old: [f64; 2], dt: f64) -> [f64; 2] {
    let mut lower = 0.0;
    let mut upper = (old[0] - old[1]) / (1.0 / 12.0 + 1.0 / 20.0);
    for _ in 0..64 {
        let heat = 0.5 * lower + 0.5 * upper;
        let a = mean_temperature(old[0], old[0] - heat / 12.0);
        let b = mean_temperature(old[1], old[1] + heat / 20.0);
        let predicted = dt * 0.1 * (a - b) * (0.4 + 0.001 * (a + b));
        if heat < predicted {
            lower = heat;
        } else {
            upper = heat;
        }
    }
    let heat = 0.5 * lower + 0.5 * upper;
    [old[0] - heat / 12.0, old[1] + heat / 20.0]
}

#[test]
fn nonlinear_conductor_and_enthalpy_bodies_share_one_conservative_graph() {
    let mut builder = SystemBuilder::new();
    let (owner, state) = body(12.0, 600.0);
    let a = builder.add_with_tangent("a", owner, state).unwrap();
    let (owner, state) = body(20.0, 300.0);
    let b = builder.add_with_tangent("b", owner, state).unwrap();
    // The local chart is moved into the component. Its scope ends before
    // sealing/stepping: production use needs no leaked or static curves.
    let wall = {
        let curve = ConductivityCurve::new(vec![[250.0, 0.9], [700.0, 1.8]]).unwrap();
        builder
            .add_with_tangent(
                "wall",
                LayerConductor::new(
                    [curve.clone(), curve],
                    LayerGeometry::planar([0.2, 0.3], 0.05).unwrap(),
                    0.0,
                    config(),
                )
                .unwrap(),
                (),
            )
            .unwrap()
    };
    let wall_a = builder.port(wall, "a").unwrap();
    let wall_b = builder.port(wall, "b").unwrap();
    builder
        .bond(builder.port(a, "heat").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(builder.port(b, "heat").unwrap(), wall_b)
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &[600.0]));
    assert!(system.guess(wall_b, &[300.0]));
    let mut expected = [600.0, 300.0];
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..60 {
            expected = coupled_oracle(expected, 0.4);
            system.step(0.4, cx).unwrap();
            for (node, temperature) in [a, b].into_iter().zip(expected) {
                close(
                    system
                        .state::<EquilibriumPhaseState>(node)
                        .unwrap()
                        .temperature_k(),
                    temperature,
                    0.0,
                    2e-9,
                );
            }
        }
    });
    let ta = system
        .state::<EquilibriumPhaseState>(a)
        .unwrap()
        .temperature_k();
    let tb = system
        .state::<EquilibriumPhaseState>(b)
        .unwrap()
        .temperature_k();
    close(
        12.0 * ta + 20.0 * tb,
        12.0 * 600.0 + 20.0 * 300.0,
        0.0,
        1e-8,
    );
    let ledger = system.ledger();
    close(ledger.energy_residual_j(), 0.0, 0.0, 1e-9);
    close(
        ledger.nodes[wall.index()].entropy_generated_j_per_k.value(),
        12.0 * (ta / 600.0).ln() + 20.0 * (tb / 300.0).ln(),
        0.0,
        1e-10,
    );
    assert_eq!(ledger.nodes[wall.index()].numerical_loss_j.value(), 0.0);
    let mean_a = system.last_input(wall_a)[0];
    let mean_b = system.last_input(wall_b)[0];
    let frozen_k_heat = 0.13 * (mean_a - mean_b);
    assert!(
        (system.last_output(wall_a)[0] - frozen_k_heat).abs() > 0.1,
        "the actual transfer must follow the changing material conductivity"
    );
}

#[test]
fn material_domain_refusal_rolls_back_the_coupled_state_clock_and_ledger() {
    let curve = ConductivityCurve::new(vec![[250.0, 1.0], [400.0, 2.0]]).unwrap();
    let mut builder = SystemBuilder::new();
    let (owner, state) = body(12.0, 350.0);
    let body = builder.add_with_tangent("body", owner, state).unwrap();
    let wall = builder
        .add_with_tangent(
            "wall",
            LayerConductor::new(
                [curve.clone(), curve],
                LayerGeometry::planar([0.2, 0.3], 0.05).unwrap(),
                0.0,
                config(),
            )
            .unwrap(),
            (),
        )
        .unwrap();
    let port = builder.port(wall, "a").unwrap();
    builder
        .bond(builder.port(body, "heat").unwrap(), port)
        .unwrap();
    builder
        .boundary(
            builder.port(wall, "b").unwrap(),
            Box::new(|t, _, values| {
                values[0] = if t < 0.25 { 300.0 } else { 500.0 };
            }),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(port, &[350.0]));
    in_cx(&CancelGate::new_clock_free(), |cx| {
        system.step(0.25, cx).unwrap();
        let state = *system.state::<EquilibriumPhaseState>(body).unwrap();
        let ledger = system.ledger().clone();
        let output = system.last_output(port).to_vec();
        assert!(system.step(0.25, cx).is_err());
        assert_eq!(system.state::<EquilibriumPhaseState>(body), Some(&state));
        assert_eq!(system.time_s(), 0.25);
        assert_eq!(system.ledger(), &ledger);
        assert_eq!(system.last_output(port), output);
    });
}

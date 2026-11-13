//! Native thermoviscous pipes, their complete tangent and an R-C-R graph.

#![cfg(feature = "component")]
#![allow(clippy::float_cmp)]

use core::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::pipe::component::CarreauPipeComponent;
use fs_material::fluid_transport::rheology::pipe::{CarreauPipeControls, LaminarCarreauPipe};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, EndOfStep, PortTangent, Proposal};
use fs_system::SystemBuilder;
use fs_system::elements::{Force, SourceAccount, Storage, Waveform};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as GAS_R;

fn in_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 57,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + absolute,
        "actual {actual:e}, expected {expected:e}, relative {relative:e}, absolute {absolute:e}"
    );
}

fn law(mode: ViscosityShiftMode, newtonian: bool, support: f64) -> ThermalCarreauYasuda {
    let (zero, infinity, n, a) = if newtonian {
        (3.0, 3.0, 0.45, 2.0)
    } else {
        (12.0, 0.0, 0.5, 1.0)
    };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(zero, infinity, 0.6, n, a, [0.0, support]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}

fn pipe(mode: ViscosityShiftMode, newtonian: bool) -> LaminarCarreauPipe {
    LaminarCarreauPipe::new(law(mode, newtonian, 1e8), 300.0, 0.01, 0.5, 1e6).unwrap()
}

fn controls() -> CarreauPipeControls {
    CarreauPipeControls {
        relative_constitutive_tolerance: 1e-13,
        relative_quadrature_tolerance: 1e-10,
        ..CarreauPipeControls::default()
    }
}

fn shift(temperature: f64) -> f64 {
    ((12_000.0 / GAS_R) * (1.0 / temperature - 1.0 / 300.0)).exp()
}

fn observed(
    component: &CarreauPipeComponent,
    inputs: &[f64],
    cx: &Cx<'_>,
) -> (Vec<f64>, fs_ports::component::Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: 0.125,
        inputs,
    };
    let Observed::Stepped {
        outputs, receipt, ..
    } = check_component(component, &case, cx).unwrap()
    else {
        panic!("a quasisteady pipe does not shorten time");
    };
    (outputs, receipt)
}

#[test]
fn poiseuille_limit_reversal_and_pressure_reference_preserve_all_port_accounts() {
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for mode in [
            ViscosityShiftMode::Vertical,
            ViscosityShiftMode::VerticalAndRate,
        ] {
            let native = pipe(mode, true);
            let component = CarreauPipeComponent::new(native, controls()).unwrap();
            assert_eq!(component.class(), Class::Smooth);
            assert_eq!(component.end_of_step(0), EndOfStep::Untold);
            assert_eq!(
                component
                    .ports()
                    .iter()
                    .map(|p| p.id().as_str())
                    .collect::<Vec<_>>(),
                ["a", "b", "heat"]
            );
            for temperature in [290.0, 300.0, 340.0] {
                let viscosity = 3.0 * shift(temperature);
                let conductance = PI * 0.01_f64.powi(4) / (8.0 * viscosity * 0.5);
                for pressure in [-6000.0, 0.0, 6000.0] {
                    let (outputs, receipt) =
                        observed(&component, &[pressure, 0.0, temperature], cx);
                    close(outputs[0], conductance * pressure, 5e-14, 0.0);
                    assert_eq!(outputs[1], -outputs[0]);
                    close(outputs[2], -conductance * pressure * pressure, 5e-14, 0.0);
                    close(receipt.dissipated_j, -0.125 * outputs[2], 0.0, 0.0);
                    close(
                        receipt.entropy_generated_j_per_k,
                        receipt.dissipated_j / temperature,
                        0.0,
                        0.0,
                    );
                    assert_eq!(receipt.stored_change_j, 0.0);
                    assert_eq!(receipt.numerical_loss_j, 0.0);
                    let (offset, _) = observed(&component, &[pressure + 1e6, 1e6, temperature], cx);
                    assert_eq!(offset, outputs);
                    let mut tangent = [f64::NAN; 3];
                    component
                        .jvp(
                            &(),
                            0.125,
                            &[pressure, 0.0, temperature],
                            &[1.0, 1.0, 0.0],
                            &mut tangent,
                            cx,
                        )
                        .unwrap();
                    assert_eq!(tangent, [0.0; 3]);
                }
            }
        }
    });
}

// Independent stress-coordinate oracle for mu(g)=12/sqrt(1+0.6g).
// The quadratic inverse g=(s²+s sqrt(s²+4))/(2*0.6), s=tau*0.6/12,
// replaces the production stress tangent, root solve and shear-coordinate
// quadrature. Positive Simpson quadrature here integrates pi*r²*g(r).
fn oracle_flow(mode: ViscosityShiftMode, pressure: f64, temperature: f64) -> f64 {
    const PANELS: usize = 4096;
    let factor = shift(temperature);
    let rate = |x: f64| {
        let stress = pressure.abs() * 0.01 / (2.0 * 0.5) * x;
        let stress = if mode == ViscosityShiftMode::Vertical {
            stress / factor
        } else {
            stress
        };
        let s = stress * 0.6 / 12.0;
        let reference = (s * s + s * (s * s + 4.0).sqrt()) / (2.0 * 0.6);
        if mode == ViscosityShiftMode::VerticalAndRate {
            reference / factor
        } else {
            reference
        }
    };
    let mut integral = rate(1.0);
    for i in 1..PANELS {
        let x = i as f64 / PANELS as f64;
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * x * x * rate(x);
    }
    (PI * 0.01_f64.powi(3) * integral / (3.0 * PANELS as f64)).copysign(pressure)
}

#[test]
fn nonlinear_original_profiles_match_independent_oracle_and_complete_thermal_tangent() {
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for mode in [
            ViscosityShiftMode::Vertical,
            ViscosityShiftMode::VerticalAndRate,
        ] {
            let native = pipe(mode, false);
            let component = CarreauPipeComponent::new(native, controls()).unwrap();
            for (pressure, temperature) in [(-8000.0, 290.0), (0.0, 310.0), (12000.0, 330.0)] {
                let inputs = [pressure, 0.0, temperature];
                let (outputs, _) = observed(&component, &inputs, cx);
                let profile = native
                    .at_temperature(temperature)
                    .unwrap()
                    .profile_at_pressure(pressure, controls())
                    .unwrap();
                assert_eq!(
                    outputs[0].to_bits(),
                    profile.flow().volumetric_flow_m3_s.to_bits()
                );
                assert_eq!(
                    outputs[2].to_bits(),
                    (-profile.flow().dissipated_power_w).to_bits()
                );
                close(
                    outputs[0],
                    oracle_flow(mode, pressure, temperature),
                    3e-10,
                    1e-20,
                );
                let case = Case {
                    state: &(),
                    stale: &(),
                    dt_s: 0.125,
                    inputs: &inputs,
                };
                check_port_tangent(&component, &case, 2e-5, 1e-10, cx).unwrap();
                let mut pushed = [0.0; 3];
                component
                    .jvp(&(), 0.125, &inputs, &[0.0, 0.0, 1.0], &mut pushed, cx)
                    .unwrap();
                let h = 1e-3;
                let derivative = (oracle_flow(mode, pressure, temperature + h)
                    - oracle_flow(mode, pressure, temperature - h))
                    / (2.0 * h);
                close(pushed[0], derivative, 2e-8, 1e-18);
                close(pushed[2], -pressure * derivative, 2e-8, 1e-14);
                assert_eq!(pushed[1], -pushed[0]);
                assert_eq!(
                    pushed[0].to_bits(),
                    profile
                        .temperature_flow_derivative_m3_s_k()
                        .unwrap()
                        .to_bits()
                );
            }
        }
    });
}

#[test]
fn original_constitutive_residual_remains_visible_without_replacing_physical_heat() {
    let controls = CarreauPipeControls {
        relative_constitutive_tolerance: 1e-4,
        ..controls()
    };
    let native = pipe(ViscosityShiftMode::Vertical, false);
    let component = CarreauPipeComponent::new(native, controls).unwrap();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let mut nonzero = false;
        for pressure in [50.0, 500.0, 5000.0] {
            let profile = native.profile_at_pressure(pressure, controls).unwrap();
            let flow = profile.flow();
            let (outputs, receipt) = observed(&component, &[pressure, 0.0, 300.0], cx);
            assert_eq!(outputs[2], -flow.dissipated_power_w);
            assert_eq!(receipt.residual, flow.pressure_residual_pa.abs());
            assert_eq!(receipt.numerical_loss_j, 0.0);
            let defect = (0.125 * (pressure * outputs[0] + outputs[2])).abs();
            assert!(defect <= receipt.energy_tolerance_j);
            nonzero |= flow.pressure_residual_pa != 0.0;
        }
        assert!(
            nonzero,
            "the fixture must exercise the admitted source residual"
        );
    });
}

#[test]
fn source_support_work_budget_bad_ports_and_cancellation_never_publish_partial_results() {
    let native = pipe(ViscosityShiftMode::Vertical, false);
    let component = CarreauPipeComponent::new(native, controls()).unwrap();
    assert!(
        CarreauPipeComponent::new(
            native,
            CarreauPipeControls {
                maximum_constitutive_iterations: 0,
                ..controls()
            }
        )
        .is_err()
    );
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let attempt = |owner: &CarreauPipeComponent, dt: f64, inputs: &[f64]| {
            let mut output = [7.0, 8.0, 9.0];
            assert!(
                owner
                    .propose(&(), dt, inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0, 9.0]);
            assert!(
                owner
                    .jvp(&(), dt, inputs, &[1.0; 3], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0, 9.0]);
            assert!(
                owner
                    .vjp(&(), dt, inputs, &[1.0; 3], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [7.0, 8.0, 9.0]);
        };
        for inputs in [
            vec![],
            vec![1.0, 0.0],
            vec![1.0, 0.0, 300.0, 2.0],
            vec![f64::NAN, 0.0, 300.0],
            vec![1.0, 0.0, f64::INFINITY],
            vec![1.0, 0.0, 279.0],
            vec![1.0, 0.0, 361.0],
            vec![1e99, 0.0, 300.0],
            vec![f64::MAX, -f64::MAX, 300.0],
        ] {
            attempt(&component, 0.1, &inputs);
        }
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            attempt(&component, dt, &[1.0, 0.0, 300.0]);
        }
        let exhausted = CarreauPipeComponent::new(
            native,
            CarreauPipeControls {
                maximum_quadrature_refinements: 1,
                ..controls()
            },
        )
        .unwrap();
        attempt(&exhausted, 0.1, &[5000.0, 0.0, 300.0]);
        // A colder vertical-and-rate preparation cannot silently clip its
        // authored physical ceiling to fit the original reduced-rate support.
        let bounded = LaminarCarreauPipe::new(
            law(ViscosityShiftMode::VerticalAndRate, false, 1e6),
            300.0,
            0.01,
            0.5,
            1e6,
        )
        .unwrap();
        let bounded = CarreauPipeComponent::new(bounded, controls()).unwrap();
        attempt(&bounded, 0.1, &[0.0, 0.0, 280.0]);
        let mut short = [7.0, 8.0];
        assert!(
            component
                .propose(&(), 0.1, &[1.0, 0.0, 300.0], &mut (), &mut short, cx)
                .is_err()
        );
        assert_eq!(short, [7.0, 8.0]);
        let mut output = [7.0, 8.0, 9.0];
        assert!(
            component
                .jvp(&(), 0.1, &[1.0, 0.0, 300.0], &[1.0, 0.0], &mut output, cx)
                .is_err()
        );
        assert_eq!(output, [7.0, 8.0, 9.0]);
        cx.cancel_gate().request();
        attempt(&component, 0.1, &[1.0, 0.0, 300.0]);
    });
}

// Independent elimination of the midpoint hydraulic compliance and the
// sensible body's entropy-conjugate (logarithmic mean) temperature. This
// oracle uses only Poiseuille's formula and the authored Arrhenius constants.
fn rcr_oracle(volume: f64, temperature: f64, dt: f64) -> (f64, f64, f64, f64) {
    let response = |next_temperature: f64| {
        let increment = next_temperature - temperature;
        let mean = if increment == 0.0 {
            temperature
        } else {
            increment / (increment / temperature).ln_1p()
        };
        let g1 = PI * 0.01_f64.powi(4) / (8.0 * (3.0 * shift(mean)) * 0.5);
        let g2 = 0.5 * g1;
        let midpoint_pressure =
            (volume / 1e-8 + 0.5 * dt * g1 * 10_000.0 / 1e-8) / (1.0 + 0.5 * dt * (g1 + g2) / 1e-8);
        let first = g1 * (10_000.0 - midpoint_pressure);
        let second = g2 * midpoint_pressure;
        let heat = (10_000.0 - midpoint_pressure) * first + midpoint_pressure * second;
        (volume + dt * (first - second), heat, dt * 10_000.0 * first)
    };
    let (mut lower, mut upper) = (temperature, temperature + 1.0);
    assert!(upper - temperature > dt * response(upper).1);
    for _ in 0..64 {
        let trial = 0.5 * lower + 0.5 * upper;
        if trial - temperature < dt * response(trial).1 {
            lower = trial;
        } else {
            upper = trial;
        }
    }
    let next_temperature = 0.5 * lower + 0.5 * upper;
    let (next_volume, heat, work) = response(next_temperature);
    (next_volume, next_temperature, dt * heat, work)
}

#[test]
fn two_native_pipes_compliance_and_enthalpy_body_share_the_kernel_and_thermal_feedback() {
    let native = pipe(ViscosityShiftMode::VerticalAndRate, true);
    let knot = |temperature| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: 10.0 * (temperature - 300.0),
        temperature_k: temperature,
        liquid_mass_fraction: 1.0,
        bulk_density_kg_m3: 1000.0,
    };
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x87; 32]),
        SolidLiquidPhase::Liquid,
        vec![knot(280.0), knot(360.0)],
    )
    .unwrap();
    let initial = chart.state_at_specific_enthalpy(0.0).unwrap();
    let mut builder = SystemBuilder::new();
    let source = builder
        .add_with_tangent(
            "pressure",
            Force::in_domain(
                PortKind::FluidPressureFlux,
                Waveform::constant(10_000.0).unwrap(),
            ),
            SourceAccount::default(),
        )
        .unwrap();
    let first = builder
        .add_with_tangent(
            "first",
            CarreauPipeComponent::new(native, controls()).unwrap(),
            (),
        )
        .unwrap();
    let second = builder
        .add_with_tangent(
            "second",
            CarreauPipeComponent::new(native.with_length(1.0).unwrap(), controls()).unwrap(),
            (),
        )
        .unwrap();
    let compliance = builder
        .add_with_tangent(
            "compliance",
            Storage::in_domain(PortKind::FluidPressureFlux, 1e-8).unwrap(),
            0.0,
        )
        .unwrap();
    // Cp=10 J/(kg K), mass=0.1 kg: exactly 1 J/K of sensible storage.
    let body = builder
        .add_with_tangent(
            "thermal",
            UniformEnthalpyBody::new(chart, 0.1).unwrap(),
            initial,
        )
        .unwrap();
    builder
        .bond(
            builder.port(source, "p").unwrap(),
            builder.port(first, "a").unwrap(),
        )
        .unwrap();
    builder
        .common_effort(
            builder.port(compliance, "p").unwrap(),
            &[
                builder.port(first, "b").unwrap(),
                builder.port(second, "a").unwrap(),
            ],
        )
        .unwrap();
    builder
        .boundary(
            builder.port(second, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    let first_heat = builder.port(first, "heat").unwrap();
    let second_heat = builder.port(second, "heat").unwrap();
    builder
        .common_effort(
            builder.port(body, "heat").unwrap(),
            &[first_heat, second_heat],
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(first_heat, &[300.0]));
    assert!(system.guess(second_heat, &[300.0]));
    let (mut volume, mut temperature, mut heat, mut work) = (0.0, 300.0, 0.0, 0.0);
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..60 {
            let next = rcr_oracle(volume, temperature, 0.2);
            volume = next.0;
            temperature = next.1;
            heat += next.2;
            work += next.3;
            system.step(0.2, cx).unwrap();
            close(
                *system.state::<f64>(compliance).unwrap(),
                volume,
                2e-10,
                1e-15,
            );
            close(
                system
                    .state::<EquilibriumPhaseState>(body)
                    .unwrap()
                    .temperature_k(),
                temperature,
                0.0,
                2e-9,
            );
        }
    });
    assert!(
        temperature > 301.0,
        "the fixture must heat enough to change the viscosity"
    );
    let ledger = system.ledger();
    let dissipated = ledger.nodes[first.index()].dissipated_j.value()
        + ledger.nodes[second.index()].dissipated_j.value();
    close(dissipated, heat, 2e-10, 1e-10);
    close(
        ledger.nodes[body.index()].stored_change_j.value(),
        heat,
        2e-10,
        1e-10,
    );
    close(
        -system.state::<SourceAccount>(source).unwrap().energy_in_j,
        work,
        2e-10,
        1e-10,
    );
    close(work, heat + 0.5 * volume * volume / 1e-8, 2e-10, 1e-10);
    close(ledger.energy_residual_j(), 0.0, 0.0, 1e-9);
    for connection in &ledger.connections {
        close(connection.energy_defect_j.value(), 0.0, 0.0, 1e-9);
    }
    let endpoint = system.port(second, "a").unwrap();
    let hot_flow = system.last_output(endpoint)[0];
    let pressure = system.last_input(endpoint)[0];
    let cold_flow = PI * 0.01_f64.powi(4) * pressure / (8.0 * 3.0 * 1.0);
    assert!(
        hot_flow > 1.01 * cold_flow,
        "source temperature must affect the actual network flow"
    );
}

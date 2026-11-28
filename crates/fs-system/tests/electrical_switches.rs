//! G0/G2/G5: ideal switching laws and circuits with independent electrical balances.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, check_component, check_load_line};
use fs_ports::component::{Component, LoadLine, PortLayout, Proposal, audit_step};
use fs_system::elements::{
    Admittance, ElementError, Force, IdealDiode, IdealSwitch, Mass, Motion, SignalClock,
    SignalSource, SourceAccount, Spring, Storage, Waveform,
};
use fs_system::{NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, heat_sink};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 27,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

/// Unlike the kit's zero-compliance comparison, this also admits operating
/// points where removing the circuit impedance would prescribe an impossible
/// current into an open circuit. It uses the shared first-law audit.
fn audited_voltage<C>(component: &C, inputs: &[f64], compliance: f64, cx: &Cx<'_>) -> f64
where
    C: LoadLine<State = (), Error = ElementError>,
{
    let dt_s = 0.125;
    let mut output = [f64::NAN];
    let first = component
        .propose_on_line(&(), dt_s, inputs, &[compliance], &mut (), &mut output, cx)
        .unwrap();
    let mut again = [123.0];
    let second = component
        .propose_on_line(&(), dt_s, inputs, &[compliance], &mut (), &mut again, cx)
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(output[0].to_bits(), again[0].to_bits());
    let Proposal::Step(receipt) = first else {
        panic!("ideal switching does not shorten a step");
    };
    assert_eq!(receipt.dissipated_j, 0.0);
    assert_eq!(receipt.numerical_loss_j, 0.0);
    let mut actual = inputs.to_vec();
    actual[0] -= compliance * output[0];
    let layout = PortLayout::new(component.ports()).unwrap();
    audit_step(
        component.ports(),
        &layout,
        dt_s,
        &actual,
        &output,
        &receipt,
        0.0,
        0.0,
    )
    .unwrap();
    output[0]
}

#[test]
fn diode_complementarity_and_switch_open_short_laws_close_the_first_law() {
    let diode = IdealDiode::new();
    let switch = IdealSwitch::new();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for current in [-13.7, -1.0e-100, 0.0, 1.0e-100, 13.7] {
            for compliance in [1.0e-5, 0.3, 2.0e90] {
                let diode_v = audited_voltage(&diode, &[current], compliance, cx);
                let expected = if current < 0.0 {
                    current / compliance
                } else {
                    0.0
                };
                assert_eq!(diode_v, expected);
                assert!(diode_v <= 0.0);
                let conducting = audited_voltage(&switch, &[current, 1.0], compliance, cx);
                assert_eq!(conducting, 0.0);
                let blocked = audited_voltage(&switch, &[current, 0.0], compliance, cx);
                assert_eq!(blocked, current / compliance);
            }
        }
        for current in [0.0, 0.25, 10.0] {
            let diode_inputs = [current];
            let diode_case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.25,
                inputs: &diode_inputs,
            };
            check_component(&diode, &diode_case, cx).unwrap();
            check_load_line(&diode, &diode_case, &[0.75], cx).unwrap();
        }
        let latch = IdealSwitch::in_domain(PortKind::MechanicalForceVelocity).unwrap();
        for inputs in [[-3.0, 1.0], [0.0, 1.0], [3.0, 1.0], [0.0, 0.0]] {
            let case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.25,
                inputs: &inputs,
            };
            check_component(&switch, &case, cx).unwrap();
            check_load_line(&switch, &case, &[0.75], cx).unwrap();
            check_component(&latch, &case, cx).unwrap();
            check_load_line(&latch, &case, &[0.75], cx).unwrap();
        }
    });
}

/// A Thevenin source: source.p -> resistor.a, resistor.b -> device.p.
/// The source fixes its voltage, the device fixes the load voltage, and
/// Kirchhoff's equations give i = (source voltage - device voltage) / R.
fn circuit<C>(element: C, voltage_v: f64, resistance_ohm: f64) -> (SystemBuilder, NodeId, NodeId)
where
    C: LoadLine<State = (), Error = ElementError> + 'static,
{
    let electrical = PortKind::ElectricalVoltageCurrent;
    let mut builder = SystemBuilder::new();
    let source = builder
        .add_with_tangent(
            "source",
            Force::in_domain(electrical, Waveform::constant(voltage_v).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let resistor = builder
        .add_with_tangent(
            "resistor",
            Admittance::in_domain(electrical, 1.0 / resistance_ohm).unwrap(),
            (),
        )
        .unwrap();
    let device = builder.add_on_line("device", element, ()).unwrap();
    builder
        .bond(
            builder.port(source, "p").unwrap(),
            builder.port(resistor, "a").unwrap(),
        )
        .unwrap();
    builder
        .bond(
            builder.port(device, "p").unwrap(),
            builder.port(resistor, "b").unwrap(),
        )
        .unwrap();
    builder
        .heat_route(
            builder.port(resistor, "heat").unwrap(),
            &[],
            Some(heat_sink("ambient", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    (builder, source, device)
}

fn seal(builder: SystemBuilder) -> System {
    let mut system = builder.seal().unwrap();
    // A blocked current is exactly zero physically; its loop still has
    // subtraction roundoff from finite terminal voltages.
    system.set_config(SolveConfig {
        absolute_tolerance: 1.0e-12,
        ..SolveConfig::default()
    });
    system
}

#[test]
fn diode_rectifies_a_real_resistive_circuit_without_spurious_heat() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for voltage in [-12.3_f64, 0.0, 12.3] {
            let (builder, source, diode) = circuit(IdealDiode::new(), voltage, 4.1);
            let diode_port = builder.port(diode, "p").unwrap();
            let mut system = seal(builder);
            for _ in 0..8 {
                system.step(0.05, cx).unwrap();
                near(
                    system.last_input(diode_port)[0],
                    voltage.max(0.0) / 4.1,
                    2.0e-12,
                );
                near(system.last_output(diode_port)[0], voltage.min(0.0), 2.0e-12);
            }
            let expected_j = voltage.max(0.0).powi(2) / 4.1 * 0.4;
            near(
                system.state::<SourceAccount>(source).unwrap().energy_in_j,
                -expected_j,
                2.0e-11,
            );
            near(system.ledger().sinks[0].heat_j.value(), expected_j, 2.0e-11);
            let account = system
                .ledger()
                .nodes
                .iter()
                .find(|node| node.name == "device")
                .unwrap();
            assert_eq!(account.dissipated_j.value(), 0.0);
            assert_eq!(account.numerical_loss_j.value(), 0.0);
            near(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
        }
    });
}

#[test]
fn scheduled_switch_changes_at_exact_boundaries_and_replays_its_work() {
    let (mut builder, source, switch) = circuit(IdealSwitch::new(), 9.0, 3.0);
    let p = builder.port(switch, "p").unwrap();
    builder
        .scheduled_boundary(
            builder.port(switch, "closed").unwrap(),
            Box::new(|time, _, out| {
                out[0] = if (0.37..0.91).contains(&time) {
                    1.0
                } else {
                    0.0
                }
            }),
            vec![0.37, 0.91],
        )
        .unwrap();
    let mut system = seal(builder);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let checkpoint = system.checkpoint().unwrap();
        let run = |system: &mut System| {
            let mut reports = Vec::new();
            for (end, current, voltage) in [(0.37, 0.0, 9.0), (0.91, 3.0, 0.0), (1.2, 0.0, 9.0)] {
                reports.push(system.step(1.2 - system.time_s(), cx).unwrap());
                near(system.time_s(), end, 4.0e-16);
                near(system.last_input(p)[0], current, 2.0e-12);
                near(system.last_output(p)[0], voltage, 2.0e-12);
            }
            reports
        };
        let reports = run(&mut system);
        let ledger = system.ledger().clone();
        let account = *system.state::<SourceAccount>(source).unwrap();
        let expected_j = 9.0_f64.powi(2) / 3.0 * (0.91 - 0.37);
        near(account.energy_in_j, -expected_j, 2.0e-11);
        near(ledger.sinks[0].heat_j.value(), expected_j, 2.0e-11);
        near(ledger.energy_residual_j(), 0.0, 2.0e-11);
        system.restore(&checkpoint).unwrap();
        assert_eq!(run(&mut system), reports);
        assert_eq!(system.ledger(), &ledger);
        assert_eq!(system.state::<SourceAccount>(source).unwrap(), &account);
    });
}

#[test]
fn a_fractional_command_refuses_without_publishing_and_retry_matches_a_fresh_step() {
    let command = Rc::new(Cell::new(1.0));
    let (mut builder, source, switch) = circuit(IdealSwitch::new(), -8.0, 2.0);
    let command_in = Rc::clone(&command);
    builder
        .boundary(
            builder.port(switch, "closed").unwrap(),
            Box::new(move |_, _, out| out[0] = command_in.get()),
        )
        .unwrap();
    let mut system = seal(builder);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        system.step(0.1, cx).unwrap();
        let checkpoint = system.checkpoint().unwrap();
        let before = system.ledger().clone();
        let account = *system.state::<SourceAccount>(source).unwrap();
        command.set(0.5);
        let refusal = system.step(0.1, cx).unwrap_err().to_string();
        assert!(refusal.contains("closed command"), "{refusal}");
        assert_eq!(system.time_s(), 0.1);
        assert_eq!(system.ledger(), &before);
        assert_eq!(system.state::<SourceAccount>(source).unwrap(), &account);
        command.set(0.0);
        let report = system.step(0.1, cx).unwrap();
        let after = system.ledger().clone();
        system.restore(&checkpoint).unwrap();
        assert_eq!(system.step(0.1, cx).unwrap(), report);
        assert_eq!(system.ledger(), &after);
    });
}

#[test]
fn impossible_currents_bad_shapes_and_cancellation_refuse_without_outputs() {
    let diode = IdealDiode::new();
    let switch = IdealSwitch::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for (inputs, compliance) in [
            (vec![-1.0], vec![0.0]),
            (vec![1.0], vec![-1.0]),
            (vec![1.0], vec![f64::NAN]),
            (vec![f64::INFINITY], vec![1.0]),
            (vec![], vec![1.0]),
            (vec![1.0, 2.0], vec![1.0]),
            (vec![1.0], vec![]),
            (vec![-1.0], vec![f64::from_bits(1)]),
        ] {
            let mut output = [123.0];
            assert!(
                diode
                    .propose_on_line(&(), 0.1, &inputs, &compliance, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [123.0]);
        }
        for inputs in [[1.0, 0.0], [-1.0, 0.0], [1.0, 0.5], [1.0, 2.0]] {
            let mut output = [123.0];
            assert!(
                switch
                    .propose(&(), 0.1, &inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [123.0]);
        }
        let mut no_outputs = [];
        assert!(
            diode
                .propose(&(), 0.1, &[1.0], &mut (), &mut no_outputs, cx)
                .is_err()
        );
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                diode
                    .propose(&(), dt, &[1.0], &mut (), &mut [0.0], cx)
                    .is_err()
            );
        }
        gate.request();
        let mut output = [123.0];
        assert_eq!(
            switch.propose_on_line(&(), 0.1, &[1.0, 1.0], &[1.0], &mut (), &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(output, [123.0]);
    });
}

#[test]
fn native_codecs_are_stateless_and_the_two_laws_have_distinct_identities() {
    fn codec<C: Component<State = ()>>(component: C) {
        let mut words = vec![17];
        assert!(component.encode_state(&(), &mut words));
        assert_eq!(words, [17]);
        assert_eq!(component.decode_state(&[]), Some(()));
        assert_eq!(component.decode_state(&[0]), None);
        assert!(C::STATE_CHARTER.is_some());
    }
    codec(IdealDiode::default());
    codec(IdealSwitch::default());
    assert_ne!(IdealDiode::STATE_CHARTER, IdealSwitch::STATE_CHARTER);
    // The domain is part of the identity: electrical appends none, as before.
    let identity = |switch: IdealSwitch| {
        let mut words = Vec::new();
        switch.identity_words(&mut words);
        words
    };
    assert!(identity(IdealSwitch::new()).is_empty());
    let mechanical = IdealSwitch::in_domain(PortKind::MechanicalForceVelocity);
    let rotational = IdealSwitch::in_domain(PortKind::RotationalTorqueAngularVelocity);
    let words = [mechanical, rotational].map(|switch| identity(switch.unwrap()));
    assert_eq!(words, [[1], [2]]);
    for kind in [
        PortKind::ThermalTemperatureEntropy,
        PortKind::FluidPressureFlux,
        PortKind::MagneticMmfFluxRate,
        PortKind::ChemicalPotentialAmountFlow,
    ] {
        let reason = IdealSwitch::in_domain(kind).unwrap_err().to_string();
        assert!(reason.contains(&format!("{kind:?} domain")), "{reason}");
    }
    // Open, given a flow with no compliance, or with one so small that the effort overflows,
    // a switch refuses in its domain's words.
    let kinds = [
        PortKind::ElectricalVoltageCurrent,
        PortKind::MechanicalForceVelocity,
        PortKind::RotationalTorqueAngularVelocity,
    ];
    let words = [
        ["open circuit", "current", "voltage [V]"],
        ["open latch", "velocity", "force [N]"],
        ["open brake", "angular velocity", "torque [N m]"],
    ];
    for (kind, [open, flow, effort]) in kinds.into_iter().zip(words) {
        let switch = IdealSwitch::in_domain(kind).unwrap();
        let refuse = |compliance: f64| {
            let (gate, inputs, mut out) = (CancelGate::new_clock_free(), [1.0, 0.0], [0.0]);
            let step = with_cx(&gate, |cx| {
                switch.propose_on_line(&(), 0.1, &inputs, &[compliance], &mut (), &mut out, cx)
            });
            step.unwrap_err().to_string()
        };
        let reasons = [refuse(0.0), refuse(f64::from_bits(1))];
        let flowed = format!(
            "{open} cannot accept a prescribed nonzero {flow}; \
             supply an impedance or a physical switching model"
        );
        let named = [reasons[0].contains(&flowed), reasons[1].contains(effort)];
        assert_eq!(named, [true; 2], "{reasons:?}");
    }
}

/// A mass held still against a stretched spring by an open mechanical switch
/// is freed when the command closes at a scheduled event (binary times, so
/// steps end on it); then the midpoint rule turns `(omega x, v)` by
/// `phi = 2 atan(omega dt / 2)` a step: `x0 cos(n phi)`, `omega x0 sin(n phi)`.
#[test]
fn a_mechanical_switch_holds_a_mass_until_its_command_frees_it() {
    let (mass_kg, stiffness_n_m, stretch_m) = (0.5, 200.0, 0.01);
    let (dt_s, held, steps) = (1.0 / 1024.0, 64_u32, 464_u32);
    let release_s = f64::from(held) * dt_s;
    let mechanical = PortKind::MechanicalForceVelocity;
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let mut b = SystemBuilder::new();
    let body = Mass::new(mass_kg).unwrap();
    let mass = b.add_with_tangent("mass", body, rest).unwrap();
    let spring = Spring::in_domain(mechanical, stiffness_n_m, 0.0).unwrap();
    let spring = b.add_with_tangent("spring", spring, stretch_m).unwrap();
    let latch = IdealSwitch::in_domain(mechanical).unwrap();
    let latch = b.add_on_line("latch", latch, ()).unwrap();
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let members = [p(&b, spring, "a"), p(&b, latch, "p")];
    b.common_flow(p(&b, mass, "p"), &members).unwrap();
    let command = move |time: f64, _: f64, out: &mut [f64]| out[0] = f64::from(time >= release_s);
    b.scheduled_boundary(p(&b, latch, "closed"), Box::new(command), vec![release_s])
        .unwrap();
    b.boundary(p(&b, spring, "b"), Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    let mut system = seal(b);
    let omega = f64::sqrt(stiffness_n_m / mass_kg);
    let (phi, speed) = (2.0 * (0.5 * omega * dt_s).atan(), omega * stretch_m);
    let rounding = |size: f64, k: u32| 16.0 * f64::EPSILON * size * f64::from(k);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for k in 1..=steps {
            system.step(dt_s, cx).unwrap();
            let x = *system.state::<f64>(spring).unwrap();
            let v = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let turn = f64::from(k.saturating_sub(held)) * phi;
            near(x, stretch_m * turn.cos(), rounding(stretch_m, k));
            near(v, speed * turn.sin(), rounding(speed, k));
        }
        let residual_j = system.ledger().energy_residual_j();
        near(residual_j, 0.0, rounding(mass_kg * speed * speed, steps));
    });
}

fn rc_branch<C>(
    builder: &mut SystemBuilder,
    name: &str,
    element: C,
    supply: Waveform,
) -> (NodeId, NodeId, NodeId, fs_system::PortRef)
where
    C: LoadLine<State = (), Error = ElementError> + 'static,
{
    let electrical = PortKind::ElectricalVoltageCurrent;
    let source = builder
        .add_with_tangent(
            &format!("{name}_source"),
            Force::in_domain(electrical, supply),
            SourceAccount::default(),
        )
        .unwrap();
    let resistor = builder
        .add_with_tangent(
            &format!("{name}_resistor"),
            Admittance::in_domain(electrical, 0.1).unwrap(),
            (),
        )
        .unwrap();
    let capacitor = builder
        .add_with_tangent(
            &format!("{name}_capacitor"),
            Storage::in_domain(electrical, 0.01).unwrap(),
            0.0,
        )
        .unwrap();
    let device = builder.add_on_line(name, element, ()).unwrap();
    builder
        .common_flow(
            builder.port(resistor, "a").unwrap(),
            &[
                builder.port(source, "p").unwrap(),
                builder.port(capacitor, "p").unwrap(),
                builder.port(device, "p").unwrap(),
            ],
        )
        .unwrap();
    builder
        .boundary(
            builder.port(resistor, "b").unwrap(),
            Box::new(|_, _, out| out[0] = 0.0),
        )
        .unwrap();
    (
        source,
        capacitor,
        device,
        builder.port(resistor, "heat").unwrap(),
    )
}

#[test]
fn rc_branches_follow_midpoint_and_hold_charge_when_open_or_reverse_biased() {
    const R: f64 = 10.0;
    const C: f64 = 0.01;
    let mut builder = SystemBuilder::new();
    let (diode_source, diode_cap, diode, diode_heat) = rc_branch(
        &mut builder,
        "diode",
        IdealDiode::new(),
        Waveform::held_table(vec![(0.0, -12.0), (0.7, 6.0), (1.0, -18.0)]).unwrap(),
    );
    let (switch_source, switch_cap, switch, switch_heat) = rc_branch(
        &mut builder,
        "switch",
        IdealSwitch::new(),
        Waveform::constant(-12.0).unwrap(),
    );
    let thermal = builder
        .add_with_tangent(
            "thermal",
            fs_system::elements::HeatCapacity::new(10.0).unwrap(),
            300.0,
        )
        .unwrap();
    builder
        .common_effort(
            builder.port(thermal, "heat").unwrap(),
            &[diode_heat, switch_heat],
        )
        .unwrap();
    let command = builder
        .add_with_tangent(
            "command",
            SignalSource::new(vec![
                Waveform::held_table(vec![(0.0, 1.0), (0.2, 0.0), (0.35, 1.0)]).unwrap(),
            ])
            .unwrap(),
            SignalClock::default(),
        )
        .unwrap();
    let diode_p = builder.port(diode, "p").unwrap();
    let switch_p = builder.port(switch, "p").unwrap();
    builder
        .signal(
            builder.port(command, "output").unwrap(),
            builder.port(switch, "closed").unwrap(),
        )
        .unwrap();
    let mut system = seal(builder);
    for heat_port in [diode_heat, switch_heat] {
        assert!(system.guess(heat_port, &[300.0]));
    }
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut charges = [0.0; 2];
        let mut source_work = [0.0; 2];
        let mut resistor_heat = 0.0;
        let mut interrupted = 0;
        let mut reversed = 0;
        while system.time_s() < 1.4 {
            let start = system.time_s();
            let report = system.step(0.007_f64.min(1.4 - start), cx).unwrap();
            let dt = report.dt_s;
            let supply = if start < 0.7 {
                12.0
            } else if start < 1.0 {
                -6.0
            } else {
                18.0
            };
            let open = (0.2..0.35).contains(&start);
            let blocked = supply < charges[0] / C;
            for (index, capacitor, port, voltage, stopped) in [
                (0, diode_cap, diode_p, supply, blocked),
                (1, switch_cap, switch_p, 12.0, open),
            ] {
                let old_charge = charges[index];
                if !stopped {
                    // Independent implicit-midpoint solution of R I + Q/C = Vs.
                    let ratio = dt / (2.0 * R * C);
                    charges[index] =
                        ((1.0 - ratio) * old_charge + dt * voltage / R) / (1.0 + ratio);
                }
                let current = (charges[index] - old_charge) / dt;
                source_work[index] += dt * voltage * current;
                resistor_heat += dt * R * current * current;
                near(
                    *system.state::<f64>(capacitor).unwrap(),
                    charges[index],
                    5.0e-12,
                );
                near(system.last_input(port)[0], current, 2.0e-11);
                near(
                    system.last_output(port)[0],
                    if stopped {
                        voltage - charges[index] / C
                    } else {
                        0.0
                    },
                    2.0e-11,
                );
            }
            if open {
                interrupted += 1;
            }
            if blocked {
                reversed += 1;
            }
        }
        assert!(interrupted > 10 && reversed > 20);
        for (index, source) in [diode_source, switch_source].into_iter().enumerate() {
            near(
                system.state::<SourceAccount>(source).unwrap().energy_in_j,
                -source_work[index],
                5.0e-11,
            );
        }
        near(
            10.0 * (system.state::<f64>(thermal).unwrap() - 300.0),
            resistor_heat,
            5.0e-10,
        );
        near(
            source_work[0] + source_work[1] - resistor_heat,
            (charges[0] * charges[0] + charges[1] * charges[1]) / (2.0 * C),
            5.0e-11,
        );
        near(system.ledger().energy_residual_j(), 0.0, 5.0e-11);
    });
}

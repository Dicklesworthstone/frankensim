//! Signal-driven sources on real kernel storage, with independent recurrences.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};
use fs_system::SystemBuilder;
use fs_system::elements::{ControlledSource, ElementError, Mass, Motion, SourceAccount, Storage};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let key = StreamKey {
            seed: 17,
            kernel_id: 9,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            gate,
            arena,
            key,
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

#[test]
fn commanded_effort_and_flow_conform_in_each_domain_with_their_full_tangent() {
    let gate = CancelGate::new_clock_free();
    let mut state = SourceAccount::starting_at(0.25);
    state.energy_in_j = 12.0;
    let stale = SourceAccount::starting_at(-1.0);
    for kind in [
        PortKind::MechanicalForceVelocity,
        PortKind::RotationalTorqueAngularVelocity,
        PortKind::ElectricalVoltageCurrent,
        PortKind::FluidPressureFlux,
    ] {
        for source in [ControlledSource::effort(kind), ControlledSource::flow(kind)] {
            for conjugate in [-2.0, 0.0, 5.0] {
                for command in [-4.0, 0.0, 3.0] {
                    let inputs = [conjugate, command];
                    let case = Case {
                        state: &state,
                        stale: &stale,
                        dt_s: 0.125,
                        inputs: &inputs,
                    };
                    with_cx(&gate, |cx| {
                        let Observed::Stepped {
                            receipt, outputs, ..
                        } = check_component(&source, &case, cx).unwrap()
                        else {
                            panic!("controlled source does not shorten");
                        };
                        assert_eq!(outputs, [command]);
                        assert_eq!(receipt.stored_change_j, 0.125 * conjugate * command);
                        assert_eq!(receipt.dissipated_j, 0.0);
                        check_port_tangent(&source, &case, 1.0e-7, 1.0e-8, cx).unwrap();
                    });
                }
            }
        }
    }
}

#[test]
fn reciprocal_effort_and_flow_scales_retain_representable_work() {
    let gate = CancelGate::new_clock_free();
    let tiny = 2.0_f64.powi(-600);
    let large = 2.0_f64.powi(600);
    let state = SourceAccount::default();
    let stale = SourceAccount::starting_at(1.0);
    for source in [
        ControlledSource::effort(PortKind::ElectricalVoltageCurrent),
        ControlledSource::flow(PortKind::ElectricalVoltageCurrent),
    ] {
        for (dt_s, conjugate, command, expected_work) in [
            (tiny, large, tiny, tiny),
            (large, tiny, large, large),
            (tiny, large, -tiny, -tiny),
            (large, tiny, -large, -large),
        ] {
            // Power is exactly +/-1 W. Forming dt*command first would
            // instead underflow to zero or overflow before the conjugate
            // factor can restore the represented physical scale.
            let inputs = [conjugate, command];
            let case = Case {
                state: &state,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            with_cx(&gate, |cx| {
                let Observed::Stepped { receipt, audit, .. } =
                    check_component(&source, &case, cx).unwrap()
                else {
                    panic!("controlled source does not shorten");
                };
                assert_eq!(receipt.stored_change_j, expected_work);
                assert_eq!(audit.port_energy_in_j, expected_work);
            });
        }
    }
}

#[test]
fn current_commands_charge_and_recover_capacitor_energy_with_the_source_sign() {
    let gate = CancelGate::new_clock_free();
    let kind = PortKind::ElectricalVoltageCurrent;
    let mut builder = SystemBuilder::new();
    let capacitor = builder
        .add_with_tangent("capacitor", Storage::in_domain(kind, 0.25).unwrap(), 0.5)
        .unwrap();
    let source = builder
        .add_with_tangent(
            "source",
            ControlledSource::flow(kind),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .bond(
            builder.port(capacitor, "p").unwrap(),
            builder.port(source, "p").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(source, "command").unwrap(),
            Box::new(|t, _, out| {
                out[0] = if t < 0.1 { 2.0 } else { -3.0 };
            }),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let (mut charge, mut received) = (0.5, 0.0);
    with_cx(&gate, |cx| {
        for command in [2.0, -3.0] {
            let next = charge - 0.1 * command;
            received += 0.1 * command * (charge + next) / (2.0 * 0.25);
            system.step(0.1, cx).unwrap();
            charge = next;
            near(*system.state::<f64>(capacitor).unwrap(), charge, 1.0e-13);
            let account = system.state::<SourceAccount>(source).unwrap();
            near(account.energy_in_j, received, 1.0e-13);
            assert_eq!(account.time_s(), system.time_s());
            near(system.ledger().energy_residual_j(), 0.0, 1.0e-13);
        }
    });
    near(
        received + charge * charge / (2.0 * 0.25),
        0.5 * 0.5 / (2.0 * 0.25),
        1.0e-13,
    );
}

fn signal(name: &str, output: bool) -> PortSpec {
    let id = StableId::new(name).unwrap();
    if output {
        PortSpec::signal_out(id, NonZeroUsize::MIN)
    } else {
        PortSpec::signal_in(id, NonZeroUsize::MIN)
    }
}

// This observation adds no integration or physics: the existing Mass owns
// the motion and work, while the signal reports its accepted speed.
struct ObservedMass {
    mass: Mass,
    ports: Vec<PortSpec>,
}

impl Component for ObservedMass {
    type State = Motion;
    type Error = ElementError;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        self.mass.class()
    }
    fn stored_energy_j(&self, state: &Motion) -> Result<f64, ElementError> {
        self.mass.stored_energy_j(state)
    }
    fn propose(
        &self,
        old: &Motion,
        dt: f64,
        inputs: &[f64],
        new: &mut Motion,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let proposal = self
            .mass
            .propose(old, dt, inputs, new, &mut outputs[..1], cx)?;
        outputs[1] = old.velocity_m_s;
        Ok(proposal)
    }
}

// Authored proportional control, entirely on zero-power signal ports. It
// reads accepted speed, so its recurrence is explicit left-end feedback.
struct Feedback {
    ports: [PortSpec; 3],
    gain_n_s_per_m: f64,
}

impl Component for Feedback {
    type State = ();
    type Error = ElementError;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }
    fn propose(
        &self,
        _: &(),
        _: f64,
        inputs: &[f64],
        _: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        outputs[0] = self.gain_n_s_per_m * (inputs[0] - inputs[1]);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn a_signal_controller_drives_and_brakes_real_mass_through_the_generic_kernel() {
    let gate = CancelGate::new_clock_free();
    let (mass_kg, gain, dt) = (2.0, 4.0, 0.01);
    let mass = Mass::new(mass_kg).unwrap();
    let mut ports = mass.ports().to_vec();
    ports.push(signal("speed", true));
    let mut builder = SystemBuilder::new();
    let body = builder
        .add(
            "body",
            ObservedMass { mass, ports },
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let source = builder
        .add_with_tangent(
            "source",
            ControlledSource::effort(PortKind::MechanicalForceVelocity),
            SourceAccount::default(),
        )
        .unwrap();
    let feedback = builder
        .add(
            "feedback",
            Feedback {
                ports: [
                    signal("target", false),
                    signal("measured", false),
                    signal("force", true),
                ],
                gain_n_s_per_m: gain,
            },
            (),
        )
        .unwrap();
    builder
        .bond(
            builder.port(body, "p").unwrap(),
            builder.port(source, "p").unwrap(),
        )
        .unwrap();
    builder
        .signal(
            builder.port(body, "speed").unwrap(),
            builder.port(feedback, "measured").unwrap(),
        )
        .unwrap();
    builder
        .signal(
            builder.port(feedback, "force").unwrap(),
            builder.port(source, "command").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(feedback, "target").unwrap(),
            Box::new(|t, _, out| {
                out[0] = if t < 0.5 { 2.0 } else { -1.0 };
            }),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let (mut velocity, mut position, mut supplied) = (0.0, 0.0, 0.0);
    let (mut drove, mut recovered) = (false, false);
    with_cx(&gate, |cx| {
        for _ in 0..120 {
            let target = if system.time_s() < 0.5 { 2.0 } else { -1.0 };
            let force = gain * (target - velocity);
            let next = velocity + dt * force / mass_kg;
            let mean = 0.5 * (velocity + next);
            position += dt * mean;
            let work = dt * force * mean;
            supplied += work;
            drove |= work > 0.0;
            recovered |= work < 0.0;
            system.step(dt, cx).unwrap();
            velocity = next;
            let actual = system.state::<Motion>(body).unwrap();
            near(actual.velocity_m_s, velocity, 1.0e-10);
            near(actual.position_m, position, 1.0e-10);
        }
    });
    assert!(
        drove && recovered,
        "the same source delivers and recovers physical work"
    );
    near(supplied, 0.5 * mass_kg * velocity * velocity, 1.0e-12);
    near(
        system.state::<SourceAccount>(source).unwrap().energy_in_j,
        -supplied,
        1.0e-10,
    );
    near(system.ledger().energy_residual_j(), 0.0, 1.0e-10);
    assert_eq!(
        system.ledger().nodes[feedback.index()]
            .stored_change_j
            .value(),
        0.0
    );
}

#[test]
fn source_refusal_and_cancellation_publish_neither_state_nor_ports() {
    let gate = CancelGate::new_clock_free();
    let source = ControlledSource::effort(PortKind::MechanicalForceVelocity);
    let old = SourceAccount::default();
    let stale = SourceAccount::starting_at(4.0);
    with_cx(&gate, |cx| {
        for (dt, inputs) in [
            (0.0, &[1.0, 2.0][..]),
            (-1.0, &[1.0, 2.0][..]),
            (f64::NAN, &[1.0, 2.0][..]),
            (1.0, &[][..]),
            (1.0, &[1.0][..]),
            (1.0, &[1.0, 2.0, 3.0][..]),
            (1.0, &[1.0, f64::INFINITY][..]),
            (1.0, &[f64::MAX, f64::MAX][..]),
        ] {
            let mut new = stale.clone();
            let mut out = [91.0];
            assert!(
                source
                    .propose(&old, dt, inputs, &mut new, &mut out, cx)
                    .is_err()
            );
            assert_eq!(new, stale);
            assert_eq!(out, [91.0]);
        }
        let mut tangent = [91.0];
        assert!(
            source
                .jvp(&old, 1.0, &[1.0, 2.0], &[3.0], &mut tangent, cx)
                .is_err()
        );
        assert_eq!(tangent, [91.0]);
        let mut covector = [91.0, 92.0];
        assert!(
            source
                .vjp(&old, 1.0, &[1.0, 2.0], &[f64::NAN], &mut covector, cx)
                .is_err()
        );
        assert_eq!(covector, [91.0, 92.0]);
        cx.cancel_gate().request();
        let mut new = stale.clone();
        let mut out = [91.0];
        assert_eq!(
            source.propose(&old, 0.1, &[1.0, 2.0], &mut new, &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            source.jvp(&old, 0.1, &[1.0, 2.0], &[3.0, 4.0], &mut tangent, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            source.vjp(&old, 0.1, &[1.0, 2.0], &[3.0], &mut covector, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(new, stale);
        assert_eq!(out, [91.0]);
        assert_eq!(tangent, [91.0]);
        assert_eq!(covector, [91.0, 92.0]);
    });
}

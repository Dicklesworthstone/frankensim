//! G0/G2/G3: signal-modulated power conversion, independently solved graphs,
//! complete port derivatives and refusal without partial publication.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, EndOfStep, PortSpec, PortTangent, Proposal};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{ElementError, Mass, ModulatedTransformer, Motion, Storage};
use fs_system::{NodeId, System, SystemBuilder};

const A: PortKind = PortKind::MechanicalForceVelocity;
const B: PortKind = PortKind::ElectricalVoltageCurrent;
const MASS: f64 = 2.0;
const CAPACITY: f64 = 0.125;
const INITIAL: Motion = Motion {
    position_m: 0.13,
    velocity_m_s: -0.4,
};
const CHARGE: f64 = 0.2;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 13,
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
        "{actual} != {expected} (tolerance {tolerance})"
    );
}

#[test]
fn signed_and_zero_modulation_conform_with_the_complete_ratio_derivative() {
    let gate = CancelGate::new_clock_free();
    for kind_b in [
        A,
        B,
        PortKind::RotationalTorqueAngularVelocity,
        PortKind::FluidPressureFlux,
    ] {
        let transformer = ModulatedTransformer::new(A, kind_b);
        assert_eq!(transformer.class(), Class::Smooth);
        assert_eq!(transformer.end_of_step(1), EndOfStep::Untold);
        with_cx(&gate, |cx| {
            for ratio in [-2.75, 0.0, 0.31, 2.0] {
                for (flow, effort) in [(-1.5, 3.0), (2.0, -4.0), (0.0, 5.0)] {
                    let inputs = [flow, effort, ratio];
                    let case = Case {
                        state: &(),
                        stale: &(),
                        dt_s: 0.125,
                        inputs: &inputs,
                    };
                    let Observed::Stepped {
                        outputs, receipt, ..
                    } = check_component(&transformer, &case, cx).unwrap()
                    else {
                        panic!("a transformer never shortens the step");
                    };
                    assert_eq!(outputs, [ratio * effort, -ratio * flow]);
                    assert_eq!(receipt.stored_change_j, 0.0);
                    assert_eq!(receipt.dissipated_j, 0.0);
                    assert_eq!(receipt.numerical_loss_j, 0.0);
                    check_port_tangent(&transformer, &case, 1.0e-7, 1.0e-8, cx).unwrap();

                    // Vary only the ratio. It must influence BOTH ports,
                    // including at a zero ratio where physical columns vanish.
                    let mut direction = [f64::NAN; 2];
                    transformer
                        .jvp(&(), 0.125, &inputs, &[0.0, 0.0, 1.0], &mut direction, cx)
                        .unwrap();
                    assert_eq!(direction, [effort, -flow]);
                }
            }
        });
    }
}

fn graph(ratio: f64, modulated: bool) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Mass::new(MASS).unwrap(), INITIAL)
        .unwrap();
    let capacitor = builder
        .add_with_tangent(
            "capacitor",
            Storage::in_domain(B, CAPACITY).unwrap(),
            CHARGE,
        )
        .unwrap();
    if modulated {
        let transformer = builder
            .add_with_tangent("transformer", ModulatedTransformer::new(A, B), ())
            .unwrap();
        wire(&mut builder, mass, capacitor, transformer);
        builder
            .boundary(
                builder.port(transformer, "ratio").unwrap(),
                Box::new(move |_, _, values| values[0] = ratio),
            )
            .unwrap();
    } else {
        builder
            .transformer(
                builder.port(mass, "p").unwrap(),
                builder.port(capacitor, "p").unwrap(),
                ratio,
            )
            .unwrap();
    }
    (builder.seal().unwrap(), mass, capacitor)
}

fn wire(builder: &mut SystemBuilder, mass: NodeId, capacitor: NodeId, transformer: NodeId) {
    builder
        .bond(
            builder.port(mass, "p").unwrap(),
            builder.port(transformer, "a").unwrap(),
        )
        .unwrap();
    builder
        .bond(
            builder.port(transformer, "b").unwrap(),
            builder.port(capacitor, "p").unwrap(),
        )
        .unwrap();
}

fn energy(velocity: f64, charge: f64) -> f64 {
    0.5 * MASS * velocity * velocity + charge * charge / (2.0 * CAPACITY)
}

#[test]
fn a_constant_signal_matches_the_existing_fixed_transformer() {
    let gate = CancelGate::new_clock_free();
    for ratio in [-2.75, 1.8] {
        let (mut fixed, fixed_mass, fixed_capacitor) = graph(ratio, false);
        let (mut modulated, mass, capacitor) = graph(ratio, true);
        with_cx(&gate, |cx| {
            for _ in 0..500 {
                fixed.step(0.01, cx).unwrap();
                modulated.step(0.01, cx).unwrap();
                let a = fixed.state::<Motion>(fixed_mass).unwrap();
                let b = modulated.state::<Motion>(mass).unwrap();
                near(a.position_m, b.position_m, 1.0e-11);
                near(a.velocity_m_s, b.velocity_m_s, 1.0e-11);
                near(
                    *fixed.state::<f64>(fixed_capacitor).unwrap(),
                    *modulated.state::<f64>(capacitor).unwrap(),
                    1.0e-11,
                );
            }
        });
        near(modulated.ledger().energy_residual_j(), 0.0, 1.0e-11);
        assert_eq!(fixed.time_s(), modulated.time_s());
    }
}

#[test]
fn changing_sign_and_zero_ratio_follow_the_independent_midpoint_rotation() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Mass::new(MASS).unwrap(), INITIAL)
        .unwrap();
    let capacitor = builder
        .add_with_tangent(
            "capacitor",
            Storage::in_domain(B, CAPACITY).unwrap(),
            CHARGE,
        )
        .unwrap();
    let transformer = builder
        .add_with_tangent("transformer", ModulatedTransformer::new(A, B), ())
        .unwrap();
    wire(&mut builder, mass, capacitor, transformer);
    let ratio_at = |t: f64| {
        if t < 1.0 {
            1.0 + t
        } else if t < 2.0 {
            0.0
        } else {
            -0.5 * t
        }
    };
    builder
        .boundary(
            builder.port(transformer, "ratio").unwrap(),
            Box::new(move |t, _, values| values[0] = ratio_at(t)),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let (dt, mut angle) = (0.01, 0.0_f64);
    let omega = 1.0 / (MASS * CAPACITY).sqrt();
    let (initial_v, initial_q) = (MASS.sqrt() * INITIAL.velocity_m_s, CHARGE / CAPACITY.sqrt());
    with_cx(&gate, |cx| {
        for _ in 0..350 {
            // In energy-normalized coordinates each held ratio rotates by
            // 2 atan(r*omega*dt/2), including zero and negative rotations.
            let ratio = ratio_at(system.time_s());
            angle += 2.0 * (0.5 * ratio * omega * dt).atan();
            system.step(dt, cx).unwrap();
            let (sin, cos) = angle.sin_cos();
            let expected_v = (initial_v * cos + initial_q * sin) / MASS.sqrt();
            let expected_q = (initial_q * cos - initial_v * sin) * CAPACITY.sqrt();
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let charge = *system.state::<f64>(capacitor).unwrap();
            near(velocity, expected_v, 2.0e-11);
            near(charge, expected_q, 2.0e-11);
            near(
                energy(velocity, charge),
                energy(INITIAL.velocity_m_s, CHARGE),
                2.0e-11,
            );
        }
    });
    near(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
}

// A sensor on the existing storage owner. Its physical step and derivative
// are delegated unchanged; only an authored ratio from accepted charge is
// published as a zero-power signal, held over the step.
struct SensedStorage {
    owner: Storage,
    ports: [PortSpec; 2],
}

impl SensedStorage {
    fn new() -> Self {
        let owner = Storage::in_domain(B, CAPACITY).unwrap();
        let ports = [
            owner.ports()[0].clone(),
            PortSpec::signal_out(StableId::new("ratio").unwrap(), NonZeroUsize::MIN),
        ];
        Self { owner, ports }
    }
}

impl Component for SensedStorage {
    type State = f64;
    type Error = ElementError;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        self.owner.class()
    }
    fn stored_energy_j(&self, state: &f64) -> Result<f64, ElementError> {
        self.owner.stored_energy_j(state)
    }
    fn propose(
        &self,
        old: &f64,
        dt: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let proposal = self
            .owner
            .propose(old, dt, inputs, new, &mut outputs[..1], cx)?;
        outputs[1] = 0.75 + 0.5 * old;
        Ok(proposal)
    }
}

impl PortTangent for SensedStorage {
    fn jvp(
        &self,
        old: &f64,
        dt: f64,
        inputs: &[f64],
        direction: &[f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.owner
            .jvp(old, dt, inputs, direction, &mut outputs[..1], cx)?;
        outputs[1] = 0.0;
        Ok(())
    }
    fn vjp(
        &self,
        old: &f64,
        dt: f64,
        inputs: &[f64],
        covector: &[f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.owner.vjp(old, dt, inputs, &covector[..1], outputs, cx)
    }
}

#[test]
fn accepted_state_feedback_modulates_the_graph_and_preserves_its_energy() {
    let gate = CancelGate::new_clock_free();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent("mass", Mass::new(MASS).unwrap(), INITIAL)
        .unwrap();
    let capacitor = builder
        .add_with_tangent("capacitor", SensedStorage::new(), CHARGE)
        .unwrap();
    let transformer = builder
        .add_with_tangent("transformer", ModulatedTransformer::new(A, B), ())
        .unwrap();
    wire(&mut builder, mass, capacitor, transformer);
    builder
        .signal(
            builder.port(capacitor, "ratio").unwrap(),
            builder.port(transformer, "ratio").unwrap(),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let (dt, mut velocity, mut charge, mut position) =
        (0.01, INITIAL.velocity_m_s, CHARGE, INITIAL.position_m);
    let (mut minimum_ratio, mut maximum_ratio) = (f64::INFINITY, f64::NEG_INFINITY);
    with_cx(&gate, |cx| {
        for _ in 0..500 {
            let ratio = 0.75 + 0.5 * charge;
            minimum_ratio = minimum_ratio.min(ratio);
            maximum_ratio = maximum_ratio.max(ratio);
            // Eliminate charge from the two independent midpoint equations:
            // m*(v1-v0) = dt*r*(q0+q1)/(2C), q1-q0 = -dt*r*(v0+v1)/2.
            let mean_velocity = (velocity + dt * ratio * charge / (2.0 * MASS * CAPACITY))
                / (1.0 + dt * dt * ratio * ratio / (4.0 * MASS * CAPACITY));
            velocity = 2.0 * mean_velocity - velocity;
            charge -= dt * ratio * mean_velocity;
            position += dt * mean_velocity;
            system.step(dt, cx).unwrap();
            let motion = system.state::<Motion>(mass).unwrap();
            near(motion.velocity_m_s, velocity, 2.0e-11);
            near(motion.position_m, position, 2.0e-11);
            near(*system.state::<f64>(capacitor).unwrap(), charge, 2.0e-11);
            near(
                energy(velocity, charge),
                energy(INITIAL.velocity_m_s, CHARGE),
                2.0e-11,
            );
        }
    });
    assert!(
        maximum_ratio - minimum_ratio > 0.1,
        "feedback must actually change the ratio"
    );
    near(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
}

#[test]
fn malformed_unrepresentable_and_cancelled_calls_preserve_buffers() {
    let gate = CancelGate::new_clock_free();
    let transformer = ModulatedTransformer::new(A, B);
    with_cx(&gate, |cx| {
        for (dt, inputs) in [
            (0.0, vec![1.0, 2.0, 3.0]),
            (f64::NAN, vec![1.0, 2.0, 3.0]),
            (0.1, vec![1.0, 2.0]),
            (0.1, vec![1.0, 2.0, 3.0, 4.0]),
            (0.1, vec![1.0, 2.0, f64::INFINITY]),
            (0.1, vec![f64::NAN, 2.0, 1.0]),
            (1.0, vec![1.0, f64::MAX, 2.0]),
            (2.0, vec![1.0, f64::MAX, 1.0]),
        ] {
            let mut output = [17.0; 2];
            assert!(
                transformer
                    .propose(&(), dt, &inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [17.0; 2]);
            assert!(
                transformer
                    .jvp(&(), dt, &inputs, &[1.0; 3], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [17.0; 2]);
            let mut pull = [19.0; 3];
            assert!(
                transformer
                    .vjp(&(), dt, &inputs, &[1.0; 2], &mut pull, cx)
                    .is_err()
            );
            assert_eq!(pull, [19.0; 3]);
        }
        let inputs = [1.0, 2.0, 3.0];
        let mut short = [23.0];
        assert!(
            transformer
                .propose(&(), 0.1, &inputs, &mut (), &mut short, cx)
                .is_err()
        );
        assert_eq!(short, [23.0]);
        let mut output = [29.0; 2];
        for direction in [
            &[1.0, 2.0][..],
            &[1.0, 2.0, f64::NAN][..],
            &[f64::MAX, 0.0, 0.0][..],
        ] {
            assert!(
                transformer
                    .jvp(&(), 0.1, &inputs, direction, &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [29.0; 2]);
        }
        let mut pull = [31.0; 3];
        for covector in [&[1.0][..], &[f64::INFINITY, 0.0][..], &[0.0, f64::MAX][..]] {
            assert!(
                transformer
                    .vjp(&(), 0.1, &inputs, covector, &mut pull, cx)
                    .is_err()
            );
            assert_eq!(pull, [31.0; 3]);
        }
    });
    gate.request();
    with_cx(&gate, |cx| {
        let mut output = [37.0; 2];
        assert_eq!(
            transformer.propose(&(), 0.1, &[1.0; 3], &mut (), &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(output, [37.0; 2]);
        assert_eq!(
            transformer.jvp(&(), 0.1, &[1.0; 3], &[1.0; 3], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(output, [37.0; 2]);
        let mut pull = [41.0; 3];
        assert_eq!(
            transformer.vjp(&(), 0.1, &[1.0; 3], &[1.0; 2], &mut pull, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(pull, [41.0; 3]);
    });
}

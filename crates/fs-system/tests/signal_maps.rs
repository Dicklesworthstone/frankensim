//! G0/G2/G3: vector signal routing, loop equations and sampled physical control.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    AffineSignal, ControlledSource, ControllerState, ElementError, Mass, Motion, SampledController,
    SourceAccount,
};
use fs_system::{AdmissionError, NodeId, System, SystemBuilder};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let gate = CancelGate::new_clock_free();
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2.0e-11,
        "{actual} != {expected}"
    );
}

/// A zero-power observation of a vector input; no independent numerical law.
struct Capture([PortSpec; 1]);

impl Capture {
    fn new(width: usize) -> Self {
        Self([PortSpec::signal_in(
            StableId::new("input").unwrap(),
            NonZeroUsize::new(width).unwrap(),
        )])
    }
}

impl Component for Capture {
    type State = Vec<f64>;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.0
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &Vec<f64>) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _: &Vec<f64>,
        _: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        _: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        new.copy_from_slice(inputs);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn rectangular_rank_deficient_maps_route_vectors_and_fan_out_without_power() {
    let mut builder = SystemBuilder::new();
    let source = builder
        .add_with_tangent(
            "source",
            AffineSignal::new(1, 4, vec![0.0; 4], vec![-2.0, 3.0, 4.0, 6.0]).unwrap(),
            (),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(source, "input0").unwrap(),
            Box::new(|_, _, out| out[0] = 0.0),
        )
        .unwrap();
    let captured = builder
        .add("capture", Capture::new(3), vec![0.0; 3])
        .unwrap();
    let selected = builder.add("selected", Capture::new(1), vec![0.0]).unwrap();
    let copied = builder
        .add("copied", Capture::new(4), vec![0.0; 4])
        .unwrap();
    let output = builder.port(source, "output").unwrap();
    builder
        .signal_map(
            output,
            builder.port(captured, "input").unwrap(),
            vec![
                2.0, -1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, -4.0, 2.0, 0.0, -1.0,
            ],
        )
        .unwrap();
    builder
        .signal_map(
            output,
            builder.port(selected, "input").unwrap(),
            vec![0.0, 0.0, 1.0, 0.0],
        )
        .unwrap();
    builder
        .signal(output, builder.port(copied, "input").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.125, cx).unwrap();
    });
    assert_eq!(
        system.state::<Vec<f64>>(captured).unwrap(),
        &[-4.0, 0.0, 8.0]
    );
    assert_eq!(system.state::<Vec<f64>>(selected).unwrap(), &[4.0]);
    assert_eq!(
        system.state::<Vec<f64>>(copied).unwrap(),
        &[-2.0, 3.0, 4.0, 6.0]
    );
    assert_eq!(system.plan_shape().1, 0);
    assert_eq!(system.ledger().energy_residual_j(), 0.0);
    for node in &system.ledger().nodes {
        assert_eq!(node.stored_change_j.value(), 0.0);
        assert_eq!(node.dissipated_j.value(), 0.0);
    }
}

fn loop_graph(map: Vec<f64>, tangent: bool) -> (System, NodeId) {
    let mut builder = SystemBuilder::new();
    let law = AffineSignal::new(1, 2, vec![0.25, -0.5], vec![1.0, 3.0]).unwrap();
    let node = if tangent {
        builder.add_with_tangent("law", law, ()).unwrap()
    } else {
        builder.add("law", law, ()).unwrap()
    };
    builder
        .signal_map(
            builder.port(node, "output").unwrap(),
            builder.port(node, "input0").unwrap(),
            map,
        )
        .unwrap();
    (builder.seal().unwrap(), node)
}

#[test]
fn mapped_loop_uses_the_full_row_and_zero_gains_create_no_dependency() {
    // u = 0.5 (1 + 0.25 u) - 0.25 (3 - 0.5 u), so u = -1/3.
    with_cx(|cx| {
        for tangent in [true, false] {
            let (mut system, node) = loop_graph(vec![0.5, -0.25], tangent);
            assert_eq!(system.plan_shape().1, 1);
            for _ in 0..3 {
                system.step(0.125, cx).unwrap();
                near(
                    system.last_input(system.port(node, "input0").unwrap())[0],
                    -1.0 / 3.0,
                );
                let output = system.last_output(system.port(node, "output").unwrap());
                near(output[0], 11.0 / 12.0);
                near(output[1], 19.0 / 6.0);
            }
        }
        let (mut zero, node) = loop_graph(vec![0.0, -0.0], true);
        assert_eq!(zero.plan_shape().1, 0);
        zero.step(0.125, cx).unwrap();
        assert_eq!(zero.last_input(zero.port(node, "input0").unwrap()), &[0.0]);
        assert_eq!(
            zero.last_output(zero.port(node, "output").unwrap()),
            &[1.0, 3.0]
        );
    });
}

#[test]
fn mapped_signal_admission_is_atomic_and_coefficients_bind_problem_identity() {
    let mut builder = SystemBuilder::new();
    let law = builder
        .add_with_tangent(
            "law",
            AffineSignal::new(1, 2, vec![0.25, -0.5], vec![1.0, 3.0]).unwrap(),
            (),
        )
        .unwrap();
    let input = builder.port(law, "input0").unwrap();
    let output = builder.port(law, "output").unwrap();
    assert!(matches!(
        builder.signal_map(output, input, vec![1.0]),
        Err(AdmissionError::MapShape {
            rows: 1,
            columns: 2,
            entries: 1
        })
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            builder.signal_map(output, input, vec![0.0, value]),
            Err(AdmissionError::MapCoefficient { index: 1, .. })
        ));
    }
    assert!(matches!(
        builder.signal_map(input, output, vec![1.0, 0.0]),
        Err(AdmissionError::Role { .. })
    ));
    let body = builder
        .add_with_tangent(
            "body",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let physical = builder.port(body, "p").unwrap();
    assert!(matches!(
        builder.signal_map(output, physical, vec![1.0, 0.0]),
        Err(AdmissionError::Role { .. })
    ));
    builder
        .boundary(physical, Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    // All refusals leave the destination free, and only a successful map claims it.
    builder.signal_map(output, input, vec![0.5, -0.25]).unwrap();
    assert!(matches!(
        builder.signal(output, input),
        Err(AdmissionError::Mismatch { .. })
    ));
    assert!(matches!(
        builder.signal_map(output, input, vec![1.0, 0.0]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    builder.seal().unwrap();

    let (first, _) = loop_graph(vec![0.5, -0.25], true);
    let (same, _) = loop_graph(vec![0.5, -0.25], true);
    let (changed, _) = loop_graph(vec![0.5, -0.5], true);
    assert_eq!(first.problem_id(), same.problem_id());
    assert_ne!(first.problem_id(), changed.problem_id());
}

const PERIOD: f64 = 0.125;
const MASSES: [f64; 2] = [2.0, 3.0];
const INITIAL: [Motion; 2] = [
    Motion {
        position_m: 0.25,
        velocity_m_s: 1.0,
    },
    Motion {
        position_m: -0.5,
        velocity_m_s: -0.5,
    },
];

fn two_axis() -> (System, [NodeId; 2], [NodeId; 2], NodeId) {
    let mut builder = SystemBuilder::new();
    let controller = SampledController::new(
        PERIOD,
        0,
        2,
        2,
        vec![],
        vec![],
        vec![],
        vec![-2.0, 0.5, 1.0, -3.0],
    )
    .unwrap();
    let state = ControllerState::new(&controller, vec![]).unwrap();
    let control = builder
        .add_with_tangent("controller", controller, state)
        .unwrap();
    let mut bodies = Vec::new();
    let mut sources = Vec::new();
    for axis in 0..2 {
        let body = builder
            .add_with_tangent(
                &format!("body{axis}"),
                Mass::new(MASSES[axis]).unwrap().publishing(),
                INITIAL[axis],
            )
            .unwrap();
        let source = builder
            .add_with_tangent(
                &format!("source{axis}"),
                ControlledSource::effort(PortKind::MechanicalForceVelocity),
                SourceAccount::default(),
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
                builder.port(body, "velocity").unwrap(),
                builder.port(control, &format!("input{axis}")).unwrap(),
            )
            .unwrap();
        let mut projection = vec![0.0; 2];
        projection[axis] = 1.0;
        builder
            .signal_map(
                builder.port(control, "output").unwrap(),
                builder.port(source, "command").unwrap(),
                projection,
            )
            .unwrap();
        bodies.push(body);
        sources.push(source);
    }
    (
        builder.seal().unwrap(),
        bodies.try_into().unwrap(),
        sources.try_into().unwrap(),
        control,
    )
}

#[test]
fn vector_controller_drives_two_actuators_with_independent_work_and_replay() {
    with_cx(|cx| {
        for subdivisions in [1, 4] {
            let (mut system, bodies, sources, control) = two_axis();
            let checkpoint = system.checkpoint().unwrap();
            let mut expected = INITIAL;
            for sample in 1..=16 {
                let forces = [
                    -2.0 * expected[0].velocity_m_s + 0.5 * expected[1].velocity_m_s,
                    expected[0].velocity_m_s - 3.0 * expected[1].velocity_m_s,
                ];
                for axis in 0..2 {
                    let acceleration = forces[axis] / MASSES[axis];
                    expected[axis].position_m +=
                        PERIOD * expected[axis].velocity_m_s + 0.5 * PERIOD * PERIOD * acceleration;
                    expected[axis].velocity_m_s += PERIOD * acceleration;
                }
                for _ in 0..subdivisions {
                    system.step(PERIOD / f64::from(subdivisions), cx).unwrap();
                }
                assert_eq!(
                    system.state::<ControllerState>(control).unwrap().samples(),
                    sample
                );
                for axis in 0..2 {
                    let actual = system.state::<Motion>(bodies[axis]).unwrap();
                    near(actual.position_m, expected[axis].position_m);
                    near(actual.velocity_m_s, expected[axis].velocity_m_s);
                    near(
                        system
                            .state::<SourceAccount>(sources[axis])
                            .unwrap()
                            .energy_in_j,
                        0.5 * MASSES[axis]
                            * (INITIAL[axis].velocity_m_s.powi(2)
                                - expected[axis].velocity_m_s.powi(2)),
                    );
                    near(
                        system.last_input(system.port(sources[axis], "command").unwrap())[0],
                        forces[axis],
                    );
                }
                near(system.ledger().energy_residual_j(), 0.0);
            }
            let final_states = bodies.map(|node| *system.state::<Motion>(node).unwrap());
            let final_sources = sources.map(|node| *system.state::<SourceAccount>(node).unwrap());
            let final_ledger = system.ledger().clone();
            let final_control = system.state::<ControllerState>(control).unwrap().clone();
            let steps = system.steps();
            system.restore(&checkpoint).unwrap();
            for _ in 0..steps {
                system.step(PERIOD / f64::from(subdivisions), cx).unwrap();
            }
            assert_eq!(
                bodies.map(|node| *system.state::<Motion>(node).unwrap()),
                final_states
            );
            assert_eq!(
                sources.map(|node| *system.state::<SourceAccount>(node).unwrap()),
                final_sources
            );
            assert_eq!(
                system.state::<ControllerState>(control).unwrap(),
                &final_control
            );
            assert_eq!(system.ledger(), &final_ledger);
        }
    });
}

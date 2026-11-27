//! G0/G2/G3: reusable signal feedback driving real storage through the kernel.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{
    Case, Observed, check_component, check_exact_linear, check_port_tangent, check_stage_residual,
};
use fs_ports::component::{Class, Component, PortRole, PortTangent, StageResidual};
use fs_system::elements::{
    AffineSignal, ControlledSource, ElementError, Mass, Motion, SourceAccount,
};
use fs_system::{
    AnalysisError, BlockSolveMode, BlockSolverPolicy, MonolithicConfig, NodeId, StateMap,
    StationaryConfig, System, SystemBuilder,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 17,
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
        "{actual} != {expected}; tolerance {tolerance}"
    );
}

#[test]
fn rectangular_signal_map_conforms_and_uses_the_full_transpose() {
    let map =
        AffineSignal::new(3, 2, vec![2.0, -3.0, 0.0, -0.5, 4.0, 1.5], vec![0.25, -2.0]).unwrap();
    assert_eq!(map.class(), Class::ExactLinear);
    assert_eq!(map.ports().len(), 4);
    for (index, port) in map.ports()[..3].iter().enumerate() {
        assert_eq!(port.id().as_str(), format!("input{index}"));
        assert_eq!(port.role(), PortRole::SignalIn);
        assert_eq!(port.width(), 1);
    }
    assert_eq!(map.ports()[3].id().as_str(), "output");
    assert_eq!(map.ports()[3].role(), PortRole::SignalOut);
    assert_eq!(map.ports()[3].width(), 2);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let other = Case {
            state: &(),
            stale: &(),
            dt_s: 0.125,
            inputs: &[1.0, -2.0, 0.5],
        };
        for inputs in [[-1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.5, -0.25, 4.0]] {
            let case = Case {
                inputs: &inputs,
                ..other
            };
            let Observed::Stepped {
                receipt, outputs, ..
            } = check_component(&map, &case, cx).unwrap()
            else {
                panic!("a static signal map does not shorten");
            };
            assert_eq!(
                outputs,
                [
                    0.25 + 2.0 * inputs[0] - 3.0 * inputs[1],
                    -2.0 - 0.5 * inputs[0] + 4.0 * inputs[1] + 1.5 * inputs[2],
                ]
            );
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
            assert_eq!(receipt.energy_tolerance_j, 0.0);
            check_port_tangent(&map, &case, 1.0e-6, 1.0e-8, cx).unwrap();
            check_exact_linear(&map, &case, &other, cx).unwrap();
            check_stage_residual(&map, &case, &[], 1.0e-8, 1.0e-10, cx).unwrap();
        }
        let mut transpose = [f64::NAN; 3];
        map.vjp(&(), 0.125, other.inputs, &[3.0, -2.0], &mut transpose, cx)
            .unwrap();
        assert_eq!(transpose, [7.0, -17.0, -3.0]);
    });
}

#[test]
fn mass_observations_are_accepted_values_and_do_not_change_the_physical_step() {
    let old = Motion {
        position_m: 0.25,
        velocity_m_s: -0.5,
    };
    let stale = Motion {
        position_m: -2.0,
        velocity_m_s: 3.0,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for domain in [
            PortKind::MechanicalForceVelocity,
            PortKind::RotationalTorqueAngularVelocity,
            PortKind::ElectricalVoltageCurrent,
        ] {
            let plain = Mass::in_domain(domain, 2.0).unwrap();
            let published = plain.clone().publishing().publishing();
            assert_eq!(plain.ports().len(), 1);
            assert_eq!(published.ports().len(), 3);
            assert_eq!(published.ports()[0], plain.ports()[0]);
            assert_eq!(published.ports()[1].id().as_str(), "position");
            assert_eq!(published.ports()[2].id().as_str(), "velocity");
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s: 0.125,
                inputs: &[4.0],
            };
            check_component(&published, &case, cx).unwrap();
            check_port_tangent(&published, &case, 1.0e-6, 1.0e-8, cx).unwrap();
            let other = Case {
                state: &stale,
                stale: &old,
                inputs: &[-5.0],
                ..case
            };
            check_exact_linear(&published, &case, &other, cx).unwrap();
            let (mut a, mut b) = (stale, stale);
            let (mut plain_out, mut observed_out) = ([f64::NAN], [f64::NAN; 3]);
            let first = plain
                .propose(&old, case.dt_s, case.inputs, &mut a, &mut plain_out, cx)
                .unwrap();
            let second = published
                .propose(&old, case.dt_s, case.inputs, &mut b, &mut observed_out, cx)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(first, second);
            assert_eq!(
                observed_out,
                [plain_out[0], old.position_m, old.velocity_m_s]
            );
            assert_ne!(observed_out[2], b.velocity_m_s);
            let mut push = [f64::NAN; 3];
            published
                .jvp(&old, case.dt_s, case.inputs, &[1.0], &mut push, cx)
                .unwrap();
            assert_eq!(push, [0.5 * case.dt_s / 2.0, 0.0, 0.0]);
            let mut pull = [f64::NAN];
            published
                .vjp(
                    &old,
                    case.dt_s,
                    case.inputs,
                    &[0.0, 7.0, -9.0],
                    &mut pull,
                    cx,
                )
                .unwrap();
            assert_eq!(pull, [0.0]);
        }
    });
}

const MASS: f64 = 2.0;
const KP: f64 = 4.0;
const KV: f64 = 3.0;
const DT: f64 = 0.01;
const INITIAL: Motion = Motion {
    position_m: 0.25,
    velocity_m_s: -0.5,
};

fn target(time: f64) -> f64 {
    if time < 0.5 { 1.0 } else { -0.75 }
}

fn feedback_graph() -> (System, NodeId, NodeId, NodeId, Rc<Cell<bool>>) {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent_and_stage_residual("body", Mass::new(MASS).unwrap().publishing(), INITIAL)
        .unwrap();
    let source = builder
        .add_with_tangent_and_stage_residual(
            "source",
            ControlledSource::effort(PortKind::MechanicalForceVelocity),
            SourceAccount::default(),
        )
        .unwrap();
    let control = builder
        .add_with_tangent_and_stage_residual(
            "control",
            AffineSignal::new(3, 1, vec![KP, -KP, -KV], vec![0.0]).unwrap(),
            (),
        )
        .unwrap();
    builder
        .bond(
            builder.port(body, "p").unwrap(),
            builder.port(source, "p").unwrap(),
        )
        .unwrap();
    for (from, to) in [("position", "input1"), ("velocity", "input2")] {
        builder
            .signal(
                builder.port(body, from).unwrap(),
                builder.port(control, to).unwrap(),
            )
            .unwrap();
    }
    builder
        .signal(
            builder.port(control, "output").unwrap(),
            builder.port(source, "command").unwrap(),
        )
        .unwrap();
    let fail = Rc::new(Cell::new(false));
    let fail_input = Rc::clone(&fail);
    builder
        .boundary(
            builder.port(control, "input0").unwrap(),
            Box::new(move |t, _, out| {
                out[0] = if fail_input.get() {
                    // The boundary is finite; the affine law itself refuses
                    // the unrepresentable command after the graph is planned.
                    f64::MAX
                } else {
                    target(t)
                };
            }),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    system.set_solver_policy(BlockSolverPolicy::InterfaceNewton);
    (system, body, source, control, fail)
}

#[test]
fn monolithic_feedback_matches_interface_solve_and_independent_work() {
    let (mut interface, body, source, _, _) = feedback_graph();
    let (mut monolithic, _, _, _, _) = feedback_graph();
    monolithic.set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(
        MonolithicConfig::default(),
    ));
    assert!(
        monolithic
            .execution_plan_for_dt(DT)
            .blocks
            .iter()
            .any(|block| { block.solve_mode == BlockSolveMode::MonolithicNewtonKrylov })
    );
    let mut expected = INITIAL;
    let mut supplied_j = 0.0;
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..120 {
            // Left-end observed feedback and midpoint physical work, derived
            // independently of either kernel solve mode.
            let force = KP * (target(monolithic.time_s()) - expected.position_m)
                - KV * expected.velocity_m_s;
            let next_velocity = expected.velocity_m_s + DT * force / MASS;
            let displacement = DT * 0.5 * (expected.velocity_m_s + next_velocity);
            expected.position_m += displacement;
            expected.velocity_m_s = next_velocity;
            supplied_j += force * displacement;
            interface.step(DT, cx).unwrap();
            monolithic.step(DT, cx).unwrap();
            let actual = monolithic.state::<Motion>(body).unwrap();
            let reference = interface.state::<Motion>(body).unwrap();
            near(actual.position_m, expected.position_m, 2.0e-10);
            near(actual.velocity_m_s, expected.velocity_m_s, 2.0e-10);
            near(actual.position_m, reference.position_m, 2.0e-10);
            near(actual.velocity_m_s, reference.velocity_m_s, 2.0e-10);
            near(
                monolithic
                    .state::<SourceAccount>(source)
                    .unwrap()
                    .energy_in_j,
                -supplied_j,
                2.0e-10,
            );
            near(monolithic.ledger().energy_residual_j(), 0.0, 2.0e-10);
        }
    });
}

struct MotionCoordinates(NodeId);

impl StateMap for MotionCoordinates {
    fn dim(&self) -> usize {
        2
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        let state = system
            .state::<Motion>(self.0)
            .ok_or(AnalysisError::Invalid { what: "body state" })?;
        out.copy_from_slice(&[state.position_m, state.velocity_m_s]);
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        system.set_accepted_state(
            self.0,
            Motion {
                position_m: state[0],
                velocity_m_s: state[1],
            },
        )?;
        Ok(())
    }
}

#[test]
fn controlled_mass_stationary_analysis_recovers_the_authored_setpoint() {
    let (mut system, body, _, _, _) = feedback_graph();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let solution = system
            .solve_stationary(&MotionCoordinates(body), StationaryConfig::default(), cx)
            .unwrap();
        near(solution.state[0], target(0.0), 1.0e-9);
        near(solution.state[1], 0.0, 1.0e-9);
        assert!(solution.residual_norm < 1.0e-10);
    });
}

#[test]
fn authored_position_and_velocity_feedback_matches_independent_dynamics_and_work() {
    let (mut system, body, source, control, _) = feedback_graph();
    let (mut position, mut velocity) = (INITIAL.position_m, INITIAL.velocity_m_s);
    let mut work = 0.0;
    let (mut drove, mut braked) = (false, false);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..200 {
            let force = KP * (target(system.time_s()) - position) - KV * velocity;
            let next = velocity + DT * force / MASS;
            let displacement = DT * 0.5 * (velocity + next);
            let exchanged = force * displacement;
            drove |= exchanged > 0.0;
            braked |= exchanged < 0.0;
            system.step(DT, cx).unwrap();
            near(
                system.last_output(system.port(body, "position").unwrap())[0],
                position,
                2.0e-12,
            );
            near(
                system.last_output(system.port(body, "velocity").unwrap())[0],
                velocity,
                2.0e-12,
            );
            position += displacement;
            velocity = next;
            work += exchanged;
            let actual = system.state::<Motion>(body).unwrap();
            near(actual.position_m, position, 2.0e-12);
            near(actual.velocity_m_s, velocity, 2.0e-12);
            near(
                system.state::<SourceAccount>(source).unwrap().energy_in_j,
                -work,
                5.0e-12,
            );
            near(system.ledger().energy_residual_j(), 0.0, 5.0e-12);
        }
    });
    assert!(drove && braked);
    near(
        work,
        0.5 * MASS * (velocity * velocity - INITIAL.velocity_m_s * INITIAL.velocity_m_s),
        5.0e-12,
    );
    assert_eq!(
        system.ledger().nodes[control.index()]
            .stored_change_j
            .value(),
        0.0
    );
    assert_eq!(
        system.ledger().nodes[control.index()].dissipated_j.value(),
        0.0
    );
}

#[test]
fn rejected_and_cancelled_feedback_steps_replay_from_the_same_checkpoint() {
    let (mut system, body, source, _, fail) = feedback_graph();
    system.set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(
        MonolithicConfig::default(),
    ));
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..12 {
            system.step(DT, cx).unwrap();
        }
        let checkpoint = system.checkpoint().unwrap();
        let old_motion = *system.state::<Motion>(body).unwrap();
        let old_source = *system.state::<SourceAccount>(source).unwrap();
        let old_ledger = system.ledger().clone();
        fail.set(true);
        assert!(system.step(DT, cx).is_err());
        assert_eq!(system.state::<Motion>(body), Some(&old_motion));
        assert_eq!(system.state::<SourceAccount>(source), Some(&old_source));
        assert_eq!(system.ledger(), &old_ledger);
        assert_eq!(system.steps(), checkpoint.steps());
        assert_eq!(system.time_s(), checkpoint.time_s());
        fail.set(false);
        for _ in 0..25 {
            system.step(DT, cx).unwrap();
        }
        let expected_motion = *system.state::<Motion>(body).unwrap();
        let expected_source = *system.state::<SourceAccount>(source).unwrap();
        let expected_ledger = system.ledger().clone();
        system.restore(&checkpoint).unwrap();
        for _ in 0..25 {
            system.step(DT, cx).unwrap();
        }
        assert_eq!(system.state::<Motion>(body), Some(&expected_motion));
        assert_eq!(
            system.state::<SourceAccount>(source),
            Some(&expected_source)
        );
        assert_eq!(system.ledger(), &expected_ledger);
        cx.cancel_gate().request();
        assert!(system.step(DT, cx).is_err());
        assert_eq!(system.state::<Motion>(body), Some(&expected_motion));
        assert_eq!(
            system.state::<SourceAccount>(source),
            Some(&expected_source)
        );
        assert_eq!(system.ledger(), &expected_ledger);
    });
}

#[test]
fn invalid_shapes_nonfinite_arithmetic_and_cancelled_calls_refuse() {
    for (n, m, matrix, bias) in [
        (0, 1, vec![], vec![0.0]),
        (1, 0, vec![], vec![]),
        (usize::MAX, 2, vec![], vec![]),
        (2, 2, vec![1.0; 3], vec![0.0; 2]),
        (1, 2, vec![1.0; 2], vec![0.0]),
        (1, 1, vec![f64::NAN], vec![0.0]),
        (1, 1, vec![1.0], vec![f64::INFINITY]),
    ] {
        assert!(AffineSignal::new(n, m, matrix, bias).is_err());
    }
    let map = AffineSignal::new(1, 1, vec![2.0], vec![0.0]).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (dt, input) in [
            (0.0, &[1.0][..]),
            (-1.0, &[1.0][..]),
            (f64::NAN, &[1.0][..]),
            (1.0, &[][..]),
            (1.0, &[1.0, 2.0][..]),
            (1.0, &[f64::NAN][..]),
            (1.0, &[f64::MAX][..]),
        ] {
            let mut output = [91.0];
            assert!(
                map.propose(&(), dt, input, &mut (), &mut output, cx)
                    .is_err()
            );
            assert!(
                map.residual(&(), dt, input, &[], &mut [], &mut output, cx)
                    .is_err()
            );
            assert!(map.accept(&(), dt, input, &[], &mut (), cx).is_err());
        }
        assert!(map.initial_guess(&(), 1.0, &[1.0], &mut [0.0]).is_err());
        assert!(map.accept(&(), 1.0, &[1.0], &[0.0], &mut (), cx).is_err());
        assert!(map.propose(&(), 1.0, &[1.0], &mut (), &mut [], cx).is_err());
        assert!(
            map.jvp(&(), 1.0, &[f64::MAX], &[0.0], &mut [0.0], cx)
                .is_err()
        );
        assert!(
            map.vjp(&(), 1.0, &[f64::MAX], &[0.0], &mut [0.0], cx)
                .is_err()
        );
        assert!(
            map.jvp(&(), 1.0, &[1.0], &[f64::MAX], &mut [0.0], cx)
                .is_err()
        );
        assert!(
            map.vjp(&(), 1.0, &[1.0], &[f64::MAX], &mut [0.0], cx)
                .is_err()
        );
        cx.cancel_gate().request();
        let mut output = [91.0];
        assert_eq!(
            map.propose(&(), 1.0, &[1.0], &mut (), &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            map.jvp(&(), 1.0, &[1.0], &[1.0], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            map.vjp(&(), 1.0, &[1.0], &[1.0], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            map.residual(&(), 1.0, &[1.0], &[], &mut [], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            map.accept(&(), 1.0, &[1.0], &[], &mut (), cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(output, [91.0]);
    });
}

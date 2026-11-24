//! G0: a smooth feedback solve stays inside its owner's source-law domain.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};
use fs_system::{NodeId, Phase, PortRef, SolveConfig, System, SystemBuilder, SystemError};

#[derive(Debug, PartialEq)]
struct OutsideSupport(f64);

/// x -> ln(x) + 2, on its real-valued domain. The signal state records
/// accepted steps and their input, so rejected proposals must not leak.
struct Logarithm {
    ports: [PortSpec; 2],
    cancel_on_recovery: bool,
    shorten_above_s: Option<f64>,
}

impl Logarithm {
    fn new() -> Self {
        Self {
            ports: [
                PortSpec::signal_in(StableId::new("input").unwrap(), NonZeroUsize::MIN),
                PortSpec::signal_out(StableId::new("output").unwrap(), NonZeroUsize::MIN),
            ],
            cancel_on_recovery: false,
            shorten_above_s: None,
        }
    }
}

impl Component for Logarithm {
    type State = [f64; 2];
    type Error = OutsideSupport;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _: &Self::State) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        let input = inputs[0];
        // A refused owner may have written candidate scratch before it
        // discovered the refusal. Every later proposal must replace it.
        *new = [f64::NAN; 2];
        outputs[0] = f64::NAN;
        if !(input.is_finite() && input > 0.0) {
            if let Some(limit) = self.shorten_above_s
                && dt_s > limit
            {
                return Ok(Proposal::Shorten { dt_s: limit });
            }
            return Err(OutsideSupport(input));
        }
        if self.cancel_on_recovery && input < 0.1 {
            cx.cancel_gate().request();
        }
        cx.checkpoint().map_err(|_| OutsideSupport(input))?;
        *new = [old[0] + 1.0, input];
        outputs[0] = input.ln() + 2.0;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

impl PortTangent for Logarithm {
    fn jvp(
        &self,
        _: &Self::State,
        _: f64,
        inputs: &[f64],
        direction: &[f64],
        outputs: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        outputs[0] = direction[0] / inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _: &Self::State,
        _: f64,
        inputs: &[f64],
        output_bar: &[f64],
        input_bar: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        input_bar[0] = output_bar[0] / inputs[0];
        Ok(())
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let gate = CancelGate::new_clock_free();
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn feedback(law: Logarithm, tangent: bool, start: f64) -> (System, NodeId, PortRef, PortRef) {
    let mut builder = SystemBuilder::new();
    let node = if tangent {
        builder.add_with_tangent("law", law, [0.0, start])
    } else {
        builder.add("law", law, [0.0, start])
    }
    .unwrap();
    let input = builder.port(node, "input").unwrap();
    let output = builder.port(node, "output").unwrap();
    builder.signal(output, input).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(input, &[start]));
    (system, node, input, output)
}

fn observed(
    system: &System,
    node: NodeId,
    input: PortRef,
    output: PortRef,
) -> ([u64; 2], u64, u64) {
    (
        (*system.state::<[f64; 2]>(node).unwrap()).map(f64::to_bits),
        system.last_input(input)[0].to_bits(),
        system.last_output(output)[0].to_bits(),
    )
}

fn reference() -> f64 {
    // Independent bracketing, on the interval where x - ln(x) - 2 is
    // strictly decreasing. No Newton iterate is used as the oracle.
    let (mut a, mut b) = (0.1_f64, 0.5_f64);
    for _ in 0..60 {
        let x = a + 0.5 * (b - a);
        if x - x.ln() - 2.0 > 0.0 {
            a = x;
        } else {
            b = x;
        }
    }
    a + 0.5 * (b - a)
}

#[test]
fn a_smooth_newton_update_outside_support_is_backtracked_to_the_valid_root() {
    // At x=0.5 the full correction is -0.80685..., outside the domain.
    let start = 0.5_f64;
    let correction = -(start - start.ln() - 2.0) / (1.0 - 1.0 / start);
    assert!(start + correction < 0.0);
    assert!(start + 0.5 * correction > 0.0);
    for tangent in [true, false] {
        let (mut system, node, input, output) = feedback(Logarithm::new(), tangent, start);
        let report = with_cx(|cx| system.step(0.125, cx)).expect("valid supported root");
        let state = system.state::<[f64; 2]>(node).unwrap();
        assert_eq!(state[0], 1.0);
        assert!((state[1] - reference()).abs() < 2.0e-13, "{state:?}");
        assert_eq!(system.last_input(input)[0], state[1]);
        assert!((system.last_output(output)[0] - state[1]).abs() < 2.0e-13);
        assert_eq!(
            (system.steps(), system.time_s(), report.shortenings),
            (1, 0.125, 0)
        );
        assert_eq!(system.ledger().energy_residual_j(), 0.0);
    }
}

#[test]
fn exhausted_backtracking_keeps_the_original_refusal_and_all_accepted_state() {
    let (mut system, node, input, output) = feedback(Logarithm::new(), true, 0.5);
    system.set_config(SolveConfig {
        max_iterations: 1,
        ..SolveConfig::default()
    });
    let ledger = system.ledger().clone();
    let before = observed(&system, node, input, output);
    let error = with_cx(|cx| system.step(0.125, cx)).expect_err("no budget for a damped update");
    let SystemError::Component {
        node: name,
        phase,
        error,
    } = error
    else {
        panic!("the first owner's typed refusal must be preserved");
    };
    assert_eq!((name.as_str(), phase), ("law", Phase::Propose));
    let refused = error.downcast_ref::<OutsideSupport>().unwrap();
    assert!((refused.0 - (-0.306_852_819_440_054_6)).abs() < f64::EPSILON);
    assert_eq!(observed(&system, node, input, output), before);
    assert_eq!(system.ledger(), &ledger);
    assert_eq!((system.steps(), system.time_s()), (0, 0.0));

    system.set_config(SolveConfig::default());
    let checkpoint = system.checkpoint().unwrap();
    let first = with_cx(|cx| system.step(0.125, cx)).unwrap();
    let after = observed(&system, node, input, output);
    let ledger = system.ledger().clone();
    system.restore(&checkpoint).unwrap();
    let replay = with_cx(|cx| system.step(0.125, cx)).unwrap();
    assert_eq!(first.iterations, replay.iterations);
    assert_eq!(observed(&system, node, input, output), after);
    assert_eq!(system.ledger(), &ledger);
}

#[test]
fn invalid_initial_inputs_remain_owner_refusals() {
    let (mut system, node, _, _) = feedback(Logarithm::new(), true, -0.5);
    let error =
        with_cx(|cx| system.step(0.125, cx)).expect_err("initial input is not in the domain");
    let SystemError::Component { error, .. } = error else {
        panic!("expected the owner's refusal");
    };
    assert_eq!(
        error.downcast_ref::<OutsideSupport>(),
        Some(&OutsideSupport(-0.5))
    );
    assert_eq!(*system.state::<[f64; 2]>(node).unwrap(), [0.0, -0.5]);
    assert_eq!((system.steps(), system.time_s()), (0, 0.0));
}

#[test]
fn cancellation_during_a_backtracked_proposal_is_not_recovered() {
    let mut law = Logarithm::new();
    law.cancel_on_recovery = true;
    let (mut system, node, _, _) = feedback(law, true, 0.5);
    let ledger = system.ledger().clone();
    let error = with_cx(|cx| system.step(0.125, cx)).expect_err("cancelled during recovery");
    assert!(matches!(error, SystemError::Cancelled), "{error}");
    assert_eq!(*system.state::<[f64; 2]>(node).unwrap(), [0.0, 0.5]);
    assert_eq!(system.ledger(), &ledger);
    assert_eq!((system.steps(), system.time_s()), (0, 0.0));
}

#[test]
fn shortening_at_a_newton_iterate_still_restarts_the_physical_step() {
    let mut law = Logarithm::new();
    law.shorten_above_s = Some(0.25);
    let (mut system, node, _, _) = feedback(law, true, 0.5);
    let report = with_cx(|cx| system.step(0.5, cx)).unwrap();
    assert_eq!((report.dt_s, report.shortenings), (0.25, 1));
    assert_eq!((system.steps(), system.time_s()), (1, 0.25));
    let state = system.state::<[f64; 2]>(node).unwrap();
    assert_eq!(state[0], 1.0);
    assert!((state[1] - reference()).abs() < 2.0e-13);
}

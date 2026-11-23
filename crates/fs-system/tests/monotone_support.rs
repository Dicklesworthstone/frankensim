//! G0: monotone feedback keeps its audit inside a narrow source-card support.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_system::{PortRef, SolveConfig, System, SystemBuilder, SystemError};

struct NarrowLaw {
    ports: [PortSpec; 2],
    center: f64,
    half_width: f64,
    cancel_on_shrunk_probe: bool,
}

impl NarrowLaw {
    fn new(cancel_on_shrunk_probe: bool) -> Self {
        let width = NonZeroUsize::new(1).unwrap();
        Self {
            ports: [
                PortSpec::signal_in(StableId::new("input").unwrap(), width),
                PortSpec::signal_out(StableId::new("output").unwrap(), width),
            ],
            center: 1.0,
            half_width: 4.0 * f64::EPSILON,
            cancel_on_shrunk_probe,
        }
    }
}

impl Component for NarrowLaw {
    type State = ();
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Monotone1d
    }

    fn stored_energy_j(&self, _: &()) -> Result<f64, Self::Error> {
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
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| "cancelled")?;
        let input = inputs[0];
        let support = (self.center - self.half_width)..=(self.center + self.half_width);
        if !support.contains(&input) {
            // The first ordinary probes are farther away. Cancellation
            // here exercises a later retry, while the real iterate is valid.
            if self.cancel_on_shrunk_probe && (input - self.center).abs() < 1.0e-10 {
                cx.cancel_gate().request();
            }
            return Err("outside source-card support");
        }
        outputs[0] = 2.0 * self.center - input;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

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

fn feedback(cancel_on_shrunk_probe: bool) -> (System, PortRef, PortRef) {
    feedback_from(NarrowLaw::new(cancel_on_shrunk_probe), 1.0 + f64::EPSILON)
}

fn feedback_from(law: NarrowLaw, guess: f64) -> (System, PortRef, PortRef) {
    let mut builder = SystemBuilder::new();
    let law = builder.add("bounded monotone law", law, ()).unwrap();
    let input = builder.port(law, "input").unwrap();
    let output = builder.port(law, "output").unwrap();
    builder.signal(output, input).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(input, &[guess]));
    (system, input, output)
}

#[test]
fn near_root_inside_narrow_support_survives_out_of_domain_tangent_probes() {
    let (mut system, input, output) = feedback(false);
    with_cx(|cx| system.step(0.125, cx).unwrap());
    assert_eq!(system.steps(), 1);
    assert_eq!(system.time_s(), 0.125);
    // x = 2 - x has the independent root x = 1. The starting residual is
    // small but nonzero, so the ordinary propagated-residual audit runs.
    assert!((system.last_input(input)[0] - 1.0).abs() <= 2.0 * f64::EPSILON);
    assert!((system.last_output(output)[0] - 1.0).abs() <= 2.0 * f64::EPSILON);
}

#[test]
fn small_supported_feedback_brackets_in_its_own_units_from_zero() {
    let root = 1.0e-15;
    let mut law = NarrowLaw::new(false);
    law.center = root;
    law.half_width = 1.0e-12;
    let (mut system, input, output) = feedback_from(law, 0.0);
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-14,
        absolute_tolerance: 0.0,
        ..SolveConfig::default()
    });
    // The initial residual is 2e-15, versus a relative allowance of 2e-29.
    // A sqrt(eps) unit floor wastes the entire default search budget on
    // out-of-support points; the physical correction brackets this root.
    with_cx(|cx| system.step(0.125, cx).unwrap());
    assert!((system.last_input(input)[0] - root).abs() <= 1.0e-14 * root);
    assert!((system.last_output(output)[0] - root).abs() <= 1.0e-14 * root);
}

#[test]
fn cancellation_during_smaller_probes_keeps_the_accepted_clock_and_guess() {
    let (mut system, input, _) = feedback(true);
    let guess = system.last_input(input)[0].to_bits();
    with_cx(|cx| {
        assert!(matches!(
            system.step(0.125, cx),
            Err(SystemError::Cancelled)
        ));
    });
    assert_eq!(system.steps(), 0);
    assert_eq!(system.time_s(), 0.0);
    assert_eq!(system.last_input(input)[0].to_bits(), guess);
}

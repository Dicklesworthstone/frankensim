//! G0: the load-line tier of the Component contract on reference
//! implementations.
//!
//! The reference is dry friction with one limit force: a pad on one
//! mechanical port, and a link between two. Neither has a state. On the line
//! `flow = released - compliance * effort` the pad holds the force that
//! brings its flow to zero if that force is inside the limit, and otherwise
//! slides at the limit. Variants broken on purpose each fail the law they
//! break.

// Exact float comparisons here assert what the relation fixes exactly: a
// held pad has no flow and no loss, a sliding one is at its limit.
#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Failure, Observed, check_component, check_load_line};
use fs_ports::component::{
    Class, Component, ContractViolation, LoadLine, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Which law, if any, the component breaks on purpose.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fault {
    None,
    /// Accounts for the loss at the released flow, not at the flow on the
    /// line.
    LossAtReleased,
    /// Lists its heat port as a port of the line.
    LineOnHeat,
}

fn mechanical(name: &str) -> PortSpec {
    PortSpec::power(
        StableId::new(name).expect("id"),
        PortKind::MechanicalForceVelocity,
        PortRole::FlowIn,
        NonZeroUsize::new(1).expect("one"),
    )
    .expect("port")
}

fn heat() -> PortSpec {
    PortSpec::heat(StableId::new("heat").expect("id"), PortRole::EffortIn).expect("heat")
}

/// The friction force [N] and the flow [m/s] on the line `released - c e`
/// for a limit force: held at zero flow inside the limit, sliding at it
/// outside. With `c` zero the flow is the released one.
fn on_line(limit_n: f64, released: f64, c: f64) -> (f64, f64) {
    if c > 0.0 && released.abs() <= c * limit_n {
        (released / c, 0.0)
    } else {
        let force = if released == 0.0 { 0.0 } else { limit_n * released.signum() };
        (force, released - c * force)
    }
}

fn friction_receipt(dt_s: f64, force_n: f64, flow: f64, kelvin: f64) -> Receipt {
    let dissipated_j = dt_s * force_n * flow;
    Receipt {
        stored_change_j: 0.0,
        dissipated_j,
        numerical_loss_j: 0.0,
        entropy_generated_j_per_k: dissipated_j / kelvin,
        energy_tolerance_j: 64.0 * f64::EPSILON * dissipated_j.abs() + f64::MIN_POSITIVE,
        residual: 0.0,
        tolerance: 0.0,
    }
}

/// Dry friction between a body and the ground: port `p` receives the
/// body's velocity and returns the force that resists it, port `heat`
/// receives a temperature and gives the loss up.
#[derive(Debug)]
struct DryPad {
    limit_n: f64,
    fault: Fault,
    ports: [PortSpec; 2],
    line: [usize; 1],
}

impl DryPad {
    fn new(limit_n: f64, fault: Fault) -> Self {
        Self {
            limit_n,
            fault,
            ports: [mechanical("p"), heat()],
            line: [usize::from(fault == Fault::LineOnHeat)],
        }
    }
}

impl Component for DryPad {
    type State = ();
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &(),
        dt_s: f64,
        inputs: &[f64],
        new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        self.propose_on_line(old, dt_s, inputs, &[0.0], new, outputs, cx)
    }
}

impl LoadLine for DryPad {
    fn line_ports(&self) -> &[usize] {
        &self.line
    }

    fn propose_on_line(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| "cancelled")?;
        let (released, kelvin, c) = (inputs[0], inputs[1], compliance[0]);
        let (force, mut flow) = on_line(self.limit_n, released, c);
        if self.fault == Fault::LossAtReleased {
            flow = released;
        }
        outputs[0] = force;
        outputs[1] = -force * flow;
        Ok(Proposal::Step(friction_receipt(dt_s, force, flow, kelvin)))
    }
}

/// [`DryPad`] whose mandatory tier does not agree with its line tier at a
/// zero compliance.
#[derive(Debug)]
struct Disagreeing(DryPad);

impl Component for Disagreeing {
    type State = ();
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        self.0.ports()
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        // Half the limit, where the line tier slides at all of it.
        let (force, flow) = on_line(0.5 * self.0.limit_n, inputs[0], 0.0);
        outputs[0] = force;
        outputs[1] = -force * flow;
        Ok(Proposal::Step(friction_receipt(dt_s, force, flow, inputs[1])))
    }
}

impl LoadLine for Disagreeing {
    fn line_ports(&self) -> &[usize] {
        self.0.line_ports()
    }

    fn propose_on_line(
        &self,
        old: &(),
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        self.0.propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

/// Dry friction between two bodies: ports `a` and `b` receive their
/// velocities and return the force and its opposite, the slip is their
/// difference. The line spans both ports.
#[derive(Debug)]
struct DryLink {
    limit_n: f64,
    ports: [PortSpec; 3],
}

impl DryLink {
    fn new(limit_n: f64) -> Self {
        Self {
            limit_n,
            ports: [mechanical("a"), mechanical("b"), heat()],
        }
    }
}

impl Component for DryLink {
    type State = ();
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &(),
        dt_s: f64,
        inputs: &[f64],
        new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 4], new, outputs, cx)
    }
}

impl LoadLine for DryLink {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn propose_on_line(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        // The force acts as (f, -f) on the two ends, so the slip it takes
        // away per unit is (1, -1) C (1, -1)^T.
        let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
        let (force, slip) = on_line(self.limit_n, inputs[0] - inputs[1], c);
        outputs[0] = force;
        outputs[1] = -force;
        outputs[2] = -force * slip;
        // The two ends' powers are each as large as the common velocity
        // makes them, and their sum is rounded at that size, not at the
        // size of the loss.
        let at_a = inputs[0] - (compliance[0] - compliance[1]) * force;
        let at_b = inputs[1] - (compliance[2] - compliance[3]) * force;
        let mut receipt = friction_receipt(dt_s, force, slip, inputs[2]);
        receipt.energy_tolerance_j +=
            64.0 * f64::EPSILON * dt_s * force.abs() * (at_a.abs() + at_b.abs());
        Ok(Proposal::Step(receipt))
    }
}

fn case(dt_s: f64, inputs: &[f64]) -> Case<'_, ()> {
    Case {
        state: &(),
        stale: &(),
        dt_s,
        inputs,
    }
}

#[test]
fn a_dry_pad_conforms_on_its_line_held_and_sliding() {
    let pad = DryPad::new(12.0, Fault::None);
    // Compliance of a 2 kg body over a 1 ms step, dt / (2 m).
    let c = 2.5e-4;
    let (mut held, mut sliding) = (0_u32, 0_u32);
    for step in -40..=40 {
        // Released flows from -8 to 8 mm/s; the pad holds up to 3 mm/s.
        let released = 2.0e-4 * f64::from(step);
        let inputs = [released, 300.0];
        let seen = with_cx(|cx| check_load_line(&pad, &case(1.0e-3, &inputs), &[c], cx)).unwrap();
        let Observed::Stepped { outputs, receipt, .. } = seen.observed else {
            panic!("the pad takes the step");
        };
        let (force, flow) = (outputs[0], seen.flows[0]);
        assert!(force.abs() <= 12.0, "{force}");
        if released.abs() <= c * 12.0 {
            // Held: no flow, no loss, and the force is what the line asks.
            assert_eq!(flow, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
            assert!((force - released / c).abs() <= 1.0e-12 * 12.0, "{force}");
            held += 1;
        } else {
            // Sliding at the limit, the way it was released.
            assert_eq!(force, 12.0 * released.signum());
            assert!(flow * released > 0.0 && flow.abs() < released.abs(), "{flow}");
            assert!(receipt.dissipated_j > 0.0);
            assert_eq!(outputs[1], -force * flow);
            sliding += 1;
        }
        // The mandatory tier at the same inputs is the line with no
        // compliance, and conforms too.
        with_cx(|cx| check_component(&pad, &case(1.0e-3, &inputs), cx)).unwrap();
    }
    assert!(held >= 20 && sliding >= 20, "{held} {sliding}");
}

#[test]
fn a_dry_link_conforms_on_a_line_over_two_ports() {
    let link = DryLink::new(5.0);
    // Two bodies of 2 kg and 0.5 kg over 1 ms, and a coupling between them
    // that is not symmetric: the kit does not ask for symmetry.
    let compliance = [2.5e-4, 3.0e-5, -1.0e-5, 1.0e-3];
    let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
    let (mut held, mut sliding) = (0_u32, 0_u32);
    for step in -30..=30 {
        let slip = 5.0e-4 * f64::from(step);
        // The same released slip about several common velocities.
        for common in [0.0, 0.3, -2.0] {
            let inputs = [common + 0.5 * slip, common - 0.5 * slip, 320.0];
            let released = inputs[0] - inputs[1];
            let seen =
                with_cx(|cx| check_load_line(&link, &case(1.0e-3, &inputs), &compliance, cx))
                    .unwrap();
            let Observed::Stepped { outputs, audit, .. } = seen.observed else {
                panic!("the link takes the step");
            };
            assert_eq!(outputs[1], -outputs[0]);
            let on_the_line = seen.flows[0] - seen.flows[1];
            if released.abs() <= c * 5.0 {
                assert!(on_the_line.abs() <= 1.0e-15, "{on_the_line}");
                held += 1;
            } else {
                assert_eq!(outputs[0], 5.0 * released.signum());
                assert!(on_the_line * released > 0.0);
                sliding += 1;
            }
            // What the two ends give up is what leaves as heat, to the
            // rounding of the ends' own powers.
            let ends_j = 1.0e-3 * 5.0 * (seen.flows[0].abs() + seen.flows[1].abs());
            assert!(
                audit.port_energy_in_j.abs() <= 64.0 * f64::EPSILON * ends_j + 1.0e-300,
                "{}",
                audit.port_energy_in_j
            );
        }
    }
    assert!(held >= 30 && sliding >= 30, "{held} {sliding}");
}

#[test]
fn the_kit_refuses_what_breaks_the_line_tier() {
    let inputs = [8.0e-3, 300.0];
    let point = case(1.0e-3, &inputs);

    // The loss accounted at the released flow: the ports, audited at the
    // flow on the line, carry less than the receipt says left as heat.
    let wrong_loss = DryPad::new(12.0, Fault::LossAtReleased);
    let refusal = with_cx(|cx| check_load_line(&wrong_loss, &point, &[2.5e-4], cx)).unwrap_err();
    assert!(
        matches!(refusal, Failure::Contract(ContractViolation::EnergyBalance { .. })),
        "{refusal}"
    );
    // The same component conforms where line and released flow coincide.
    with_cx(|cx| check_load_line(&wrong_loss, &point, &[0.0], cx)).unwrap();

    // A mandatory tier that is not the line tier at a zero compliance.
    let disagreeing = Disagreeing(DryPad::new(12.0, Fault::None));
    let refusal = with_cx(|cx| check_load_line(&disagreeing, &point, &[2.5e-4], cx)).unwrap_err();
    assert!(
        matches!(&refusal, Failure::Law { law, .. } if law.contains("zero compliance")),
        "{refusal}"
    );

    // A line on a port that carries no power.
    let on_heat = DryPad::new(12.0, Fault::LineOnHeat);
    let refusal = with_cx(|cx| check_load_line(&on_heat, &point, &[2.5e-4], cx)).unwrap_err();
    assert!(
        matches!(&refusal, Failure::Law { law, .. } if law.contains("power ports that receive a flow")),
        "{refusal}"
    );

    // A compliance that is not the line's square matrix.
    let pad = DryPad::new(12.0, Fault::None);
    let refusal = with_cx(|cx| check_load_line(&pad, &point, &[2.5e-4, 0.0], cx)).unwrap_err();
    assert!(
        matches!(
            refusal,
            Failure::Contract(ContractViolation::VectorLength {
                vector: "compliance",
                expected: 1,
                actual: 2,
            })
        ),
        "{refusal}"
    );
    // Inputs of the wrong length are refused before anything is called.
    let refusal =
        with_cx(|cx| check_load_line(&pad, &case(1.0e-3, &[8.0e-3]), &[2.5e-4], cx)).unwrap_err();
    assert!(
        matches!(
            refusal,
            Failure::Contract(ContractViolation::VectorLength { vector: "inputs", .. })
        ),
        "{refusal}"
    );
}

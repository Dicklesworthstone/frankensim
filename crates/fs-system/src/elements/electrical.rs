//! Ideal switching relations on the kernel's condensed load line.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::Carrier;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use super::{ElementError, finite, inputs_finite, nonnegative, port, positive, roundoff};

/// A lossless ideal diode, positive current flowing from anode to cathode.
///
/// Electrical port `p` receives current [A] and returns anode-minus-cathode
/// voltage [V]. Its relation is `current >= 0`, `voltage <= 0`, and
/// `current * voltage = 0`: forward conduction has zero voltage drop and
/// reverse bias has zero current. Add it with [`crate::SystemBuilder::add_on_line`]
/// so its voltage is determined by the circuit, without a large resistance
/// approximation or a fabricated source of heat.
///
/// This is a law of the representative values of one step. It does not locate
/// a zero crossing inside that step, model reverse recovery, breakdown or a
/// forward voltage drop, or resolve impulses from ideal storage constraints.
#[derive(Debug, Clone)]
pub struct IdealDiode {
    ports: [PortSpec; 1],
}

impl IdealDiode {
    /// Construct an ideal diode with no parameters or stored energy.
    #[must_use]
    pub fn new() -> Self {
        Self {
            ports: [port(
                "p",
                Some(PortKind::ElectricalVoltageCurrent),
                PortRole::FlowIn,
            )],
        }
    }
}

impl Default for IdealDiode {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for IdealDiode {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::IdealDiode",
        state_family: "stateless",
        schema_grammar: "(): ideal diode stores no energy or charge",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, ElementError> {
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
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0], new, outputs, cx)
    }
}

impl LoadLine for IdealDiode {
    fn line_ports(&self) -> &[usize] {
        &[0]
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
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        validate(dt_s, inputs, 1, compliance, outputs)?;
        let (voltage, receipt) = relation(dt_s, inputs[0], compliance[0], inputs[0] >= 0.0)?;
        cx.checkpoint()?;
        outputs[0] = voltage;
        Ok(Proposal::Step(receipt))
    }
}

/// A lossless ideal switch whose `closed` signal is exactly zero or one.
///
/// Port `p` receives a flow and returns an effort. `closed = 1` sets zero
/// effort, either flow direction allowed; `closed = 0` sets zero flow, the
/// connected system determining the effort. Electrically ([`Self::new`])
/// closed is a short; on a mechanical or rotational port ([`Self::in_domain`])
/// closed is free (a cut: no force or torque) and open holds the port still
/// (a latch or a brake). Add with [`crate::SystemBuilder::add_on_line`].
/// Inputs are ordered `[flow, closed]`; the only output is the effort.
///
/// The command is held over a step. Supply prescribed changes through
/// [`crate::SystemBuilder::scheduled_boundary`] or a source that announces
/// its events; an interval-mean fractional command is refused. This element
/// adds no contact resistance, arcing, switching loss, event locator or
/// impulse model. Instantaneous feedback that puts the binary command inside
/// the circuit's nonlinear solve is not a continuous control law.
#[derive(Debug, Clone)]
pub struct IdealSwitch {
    ports: [PortSpec; 2],
}

impl IdealSwitch {
    /// Construct an ideal switch with electrical `p` and signal `closed`.
    #[must_use]
    pub fn new() -> Self {
        Self::on(PortKind::ElectricalVoltageCurrent)
    }

    /// The same law on a port of the power domain `kind`, read as above.
    ///
    /// # Errors
    /// A thermal, fluid, magnetic or chemical `kind`, whose zero effort is
    /// not a short or a free joint (zero kelvin, zero absolute pressure).
    pub fn in_domain(kind: PortKind) -> Result<Self, ElementError> {
        match kind {
            PortKind::ElectricalVoltageCurrent
            | PortKind::MechanicalForceVelocity
            | PortKind::RotationalTorqueAngularVelocity => Ok(Self::on(kind)),
            _ => Err(ElementError::Law {
                reason: format!("an ideal switch is not admitted in the {kind:?} domain"),
            }),
        }
    }

    fn on(kind: PortKind) -> Self {
        Self {
            ports: [
                port("p", Some(kind), PortRole::FlowIn),
                PortSpec::signal_in(
                    StableId::new("closed").expect("static switch signal name"),
                    NonZeroUsize::MIN,
                ),
            ],
        }
    }
}

impl Default for IdealSwitch {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for IdealSwitch {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::IdealSwitch",
        state_family: "stateless",
        schema_grammar: "(): commanded ideal switch stores no energy or charge",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn identity_words(&self, words: &mut Vec<u64>) {
        // An electrical switch appends nothing, as before it took a domain.
        match self.ports[0].carrier() {
            Carrier::Power(PortKind::MechanicalForceVelocity) => words.push(1),
            Carrier::Power(PortKind::RotationalTorqueAngularVelocity) => words.push(2),
            _ => {}
        }
    }

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, ElementError> {
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
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0], new, outputs, cx)
    }
}

impl LoadLine for IdealSwitch {
    fn line_ports(&self) -> &[usize] {
        &[0]
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
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        validate(dt_s, inputs, 2, compliance, outputs)?;
        if inputs[1] != 0.0 && inputs[1] != 1.0 {
            return Err(ElementError::Parameter {
                name: "ideal switch closed command (0 or 1)",
                value: inputs[1],
            });
        }
        let (voltage, receipt) = relation(dt_s, inputs[0], compliance[0], inputs[1] == 1.0)?;
        cx.checkpoint()?;
        outputs[0] = voltage;
        Ok(Proposal::Step(receipt))
    }
}

fn validate(
    dt_s: f64,
    inputs: &[f64],
    input_width: usize,
    compliance: &[f64],
    outputs: &[f64],
) -> Result<(), ElementError> {
    for (name, actual, expected) in [
        ("ideal electrical input width", inputs.len(), input_width),
        ("ideal electrical output width", outputs.len(), 1),
        ("ideal electrical load-line width", compliance.len(), 1),
    ] {
        if actual != expected {
            return Err(ElementError::Parameter {
                name,
                value: actual as f64,
            });
        }
    }
    positive("ideal electrical timestep [s]", dt_s)?;
    inputs_finite(inputs)?;
    nonnegative("ideal electrical compliance [A/V]", compliance[0])?;
    Ok(())
}

/// Solve i = released - compliance * v, either v = 0 or i = 0.
fn relation(
    dt_s: f64,
    released_a: f64,
    compliance_s: f64,
    closed: bool,
) -> Result<(f64, Receipt), ElementError> {
    if closed || released_a == 0.0 {
        // At C = i = 0 every open-circuit voltage is possible. Zero is the
        // deterministic representative; no claim of a uniquely determined
        // voltage is made for an isolated ideal open circuit.
        return Ok((0.0, Receipt::lossless(0.0, 0.0)));
    }
    if compliance_s == 0.0 {
        return Err(ElementError::Law {
            reason: "an ideal open circuit cannot accept a prescribed nonzero current; \
                     supply circuit impedance or a physical switching model"
                .to_owned(),
        });
    }
    let voltage_v = finite("ideal open-circuit voltage [V]", released_a / compliance_s)?;
    if voltage_v == 0.0 {
        return Err(ElementError::Law {
            reason: "the nonzero open-circuit voltage is below representable resolution".to_owned(),
        });
    }
    // Only division and subtraction roundoff is permitted: the ideal law
    // dissipates nothing. Include one subnormal current quantum because
    // scaling a subnormal released current by eps can itself underflow.
    let current_roundoff_a = roundoff(released_a.abs()) + f64::from_bits(1);
    let allowance_j = finite(
        "ideal electrical energy allowance [J]",
        (current_roundoff_a * voltage_v.abs()) * dt_s,
    )?;
    Ok((voltage_v, Receipt::lossless(0.0, allowance_j)))
}

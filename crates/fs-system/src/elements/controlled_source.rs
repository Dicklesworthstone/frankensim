//! Ideal power sources commanded through the kernel's signal connections.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal};
use fs_ports::{PortKind, StableId};

use super::{ElementError, SourceAccount, finite, inputs_finite, port, positive, source_step};

/// An ideal effort or flow source driven by a scalar signal.
///
/// Port `p` exchanges physical power; signal input `command` selects its
/// returned effort or flow in the selected domain's SI unit. The signal
/// itself carries no energy. Both positive and negative commands are valid.
/// The [`SourceAccount`] records the signed work received at `p`, including
/// negative work when the source supplies energy to the connected device.
///
/// The command is held over the current step, just like every Component
/// input. Sampling, limiting, delay and controller state belong to the
/// component producing the signal. This source adds no control law, current
/// limit, internal impedance or power-supply model. Callers must resolve any
/// prescribed command events with their existing step/event policy.
#[derive(Debug, Clone)]
pub struct ControlledSource {
    ports: [PortSpec; 2],
}

impl ControlledSource {
    /// Prescribe effort at `p`, receiving its conjugate flow.
    ///
    /// Examples of `command` units are force [N], torque [N m], voltage [V]
    /// and pressure [Pa]. On a common-flow junction with a mass, a positive
    /// force command accelerates the mass in its positive direction; the
    /// source receives the negative of the mass's velocity.
    #[must_use]
    pub fn effort(kind: PortKind) -> Self {
        Self::new(kind, PortRole::FlowIn)
    }

    /// Prescribe flow at `p`, receiving its conjugate effort.
    ///
    /// Examples of `command` units are velocity [m/s], angular velocity
    /// [rad/s], current [A] and volume rate [m³/s]. Flow is positive INTO
    /// this source: a negative electrical command supplies positive current
    /// to a connected capacitor, matching [`super::FlowSource`].
    #[must_use]
    pub fn flow(kind: PortKind) -> Self {
        Self::new(kind, PortRole::EffortIn)
    }

    fn new(kind: PortKind, role: PortRole) -> Self {
        Self {
            ports: [
                port("p", Some(kind), role),
                PortSpec::signal_in(
                    StableId::new("command").expect("static source signal name"),
                    NonZeroUsize::MIN,
                ),
            ],
        }
    }

    fn validate(old: &SourceAccount, dt_s: f64, inputs: &[f64]) -> Result<(), ElementError> {
        width("controlled source input width", inputs.len(), 2)?;
        inputs_finite(inputs)?;
        positive("controlled source timestep [s]", dt_s)?;
        finite("controlled source clock [s]", old.time_s())?;
        finite("controlled source energy [J]", old.energy_in_j)?;
        Ok(())
    }
}

fn width(name: &'static str, actual: usize, expected: usize) -> Result<(), ElementError> {
    if actual == expected {
        Ok(())
    } else {
        Err(ElementError::Parameter {
            name,
            value: actual as f64,
        })
    }
}

impl Component for ControlledSource {
    type State = SourceAccount;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &SourceAccount) -> Result<f64, ElementError> {
        finite("controlled source energy [J]", state.energy_in_j)
    }

    fn propose(
        &self,
        old: &SourceAccount,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SourceAccount,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        Self::validate(old, dt_s, inputs)?;
        width("controlled source output width", outputs.len(), 1)?;
        let command = inputs[1];
        // Match the port's power pairing before applying duration. Multiplying
        // dt by the command first can erase or overflow otherwise representable
        // work when effort and flow have reciprocal physical scales.
        let work_j = finite(
            "controlled source work [J]",
            self.ports[0].energy_in_j(dt_s, &inputs[..1], &[command]),
        )?;
        let (next, receipt) = source_step(old, dt_s, work_j);
        finite("controlled source energy [J]", next.energy_in_j)?;
        finite("controlled source clock [s]", next.time_s())?;
        finite(
            "controlled source energy allowance [J]",
            receipt.energy_tolerance_j,
        )?;
        cx.checkpoint()?;
        *new = next;
        outputs[0] = command;
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for ControlledSource {
    fn jvp(
        &self,
        old: &SourceAccount,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(old, dt_s, inputs)?;
        width("controlled source input direction width", d_inputs.len(), 2)?;
        width(
            "controlled source output direction width",
            d_outputs.len(),
            1,
        )?;
        inputs_finite(d_inputs)?;
        cx.checkpoint()?;
        d_outputs[0] = d_inputs[1];
        Ok(())
    }

    fn vjp(
        &self,
        old: &SourceAccount,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(old, dt_s, inputs)?;
        width(
            "controlled source output covector width",
            bar_outputs.len(),
            1,
        )?;
        width(
            "controlled source input covector width",
            bar_inputs.len(),
            2,
        )?;
        inputs_finite(bar_outputs)?;
        cx.checkpoint()?;
        bar_inputs.copy_from_slice(&[0.0, bar_outputs[0]]);
        Ok(())
    }
}

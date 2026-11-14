//! Lossless power conversion whose ratio is an ordinary signal input.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{
    Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use super::{ElementError, finite, inputs_finite, port, positive, roundoff};

/// A lossless scalar transformer controlled by the signal `ratio`.
///
/// Port `a` receives flow and returns effort; port `b` receives effort and
/// returns flow. With all flows positive into this component, its law is
/// `e_a = ratio * e_b`, `f_b = -ratio * f_a`. Both products use the same
/// admitted ratio, so their powers cancel, including when the ratio changes
/// sign or passes through zero. The signal carries no power.
///
/// Join `a` to an inertia and `b` to a compliance, resistance or effort
/// source. A constant ratio reproduces [`crate::SystemBuilder::transformer`]
/// between those owners. The power domains may differ; the ratio then has
/// the required conversion units between their SI effort/flow coordinates.
///
/// The ratio is held over the step. Its controller, geometry or material
/// owner must supply the appropriate value and resolve changes in it. This
/// interface law does not derive crank/cam geometry, enforce a displacement
/// constraint, model a ratio-changing actuator, or connect two integral-
/// causality inertias. It returns only a mean flow, not its endpoint; a load
/// line that needs the endpoint cannot be connected at `b`.
#[derive(Debug, Clone)]
pub struct ModulatedTransformer {
    ports: [PortSpec; 3],
}

impl ModulatedTransformer {
    /// Build a transformer between two scalar power domains.
    ///
    /// Its input order is `[flow_a, effort_b, ratio]`, and its output order
    /// is `[effort_a, flow_b]`. Any finite signed ratio is admitted, including
    /// zero: this causal map uses no inverse or division by the ratio.
    #[must_use]
    pub fn new(kind_a: PortKind, kind_b: PortKind) -> Self {
        Self {
            ports: [
                port("a", Some(kind_a), PortRole::FlowIn),
                port("b", Some(kind_b), PortRole::EffortIn),
                PortSpec::signal_in(
                    StableId::new("ratio").expect("static transformer signal name"),
                    NonZeroUsize::MIN,
                ),
            ],
        }
    }

    fn validate(dt_s: f64, inputs: &[f64]) -> Result<(), ElementError> {
        extent("transformer input width", inputs.len(), 3)?;
        positive("transformer timestep [s]", dt_s)?;
        inputs_finite(inputs)
    }

    fn response(inputs: &[f64]) -> Result<[f64; 2], ElementError> {
        let [flow_a, effort_b, ratio] = [inputs[0], inputs[1], inputs[2]];
        Ok([
            finite("transformer effort at a", ratio * effort_b)?,
            finite("transformer flow at b", -ratio * flow_a)?,
        ])
    }

    fn allowance(
        &self,
        dt_s: f64,
        inputs: &[f64],
        response: &[f64; 2],
    ) -> Result<f64, ElementError> {
        let work_a = finite(
            "transformer work at a [J]",
            self.ports[0].energy_in_j(dt_s, &inputs[..1], &response[..1]),
        )?;
        let work_b = finite(
            "transformer work at b [J]",
            self.ports[1].energy_in_j(dt_s, &inputs[1..2], &response[1..]),
        )?;
        // Scale each work before adding, so two large opposing exchanges
        // need not form an overflowing magnitude sum to bound roundoff.
        finite(
            "transformer energy allowance [J]",
            roundoff(work_a.abs()) + roundoff(work_b.abs()),
        )
    }
}

fn extent(name: &'static str, actual: usize, expected: usize) -> Result<(), ElementError> {
    if actual == expected {
        Ok(())
    } else {
        Err(ElementError::Parameter {
            name,
            value: actual as f64,
        })
    }
}

impl Component for ModulatedTransformer {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::ModulatedTransformer",
        state_family: "no-state",
        schema_grammar: "(): an ideal signal-modulated transformer stores nothing",
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
        // Each power output is bilinear in its received value and ratio.
        // Constant-ratio graphs are affine, but a signal may change it.
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        EndOfStep::Untold
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, inputs)?;
        extent("transformer output width", outputs.len(), 2)?;
        let response = Self::response(inputs)?;
        let allowance = self.allowance(dt_s, inputs, &response)?;
        cx.checkpoint()?;
        outputs.copy_from_slice(&response);
        Ok(Proposal::Step(Receipt::lossless(0.0, allowance)))
    }
}

impl PortTangent for ModulatedTransformer {
    fn jvp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, inputs)?;
        extent("transformer input direction width", d_inputs.len(), 3)?;
        extent("transformer output direction width", d_outputs.len(), 2)?;
        inputs_finite(d_inputs)?;
        let response = Self::response(inputs)?;
        self.allowance(dt_s, inputs, &response)?;
        let push = [
            finite(
                "transformer effort direction",
                inputs[2] * d_inputs[1] + inputs[1] * d_inputs[2],
            )?,
            finite(
                "transformer flow direction",
                -inputs[2] * d_inputs[0] - inputs[0] * d_inputs[2],
            )?,
        ];
        cx.checkpoint()?;
        d_outputs.copy_from_slice(&push);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        Self::validate(dt_s, inputs)?;
        extent("transformer output covector width", bar_outputs.len(), 2)?;
        extent("transformer input covector width", bar_inputs.len(), 3)?;
        inputs_finite(bar_outputs)?;
        let response = Self::response(inputs)?;
        self.allowance(dt_s, inputs, &response)?;
        let pull = [
            finite("transformer flow covector", -inputs[2] * bar_outputs[1])?,
            finite("transformer effort covector", inputs[2] * bar_outputs[0])?,
            finite(
                "transformer ratio covector",
                inputs[1] * bar_outputs[0] - inputs[0] * bar_outputs[1],
            )?,
        ];
        cx.checkpoint()?;
        bar_inputs.copy_from_slice(&pull);
        Ok(())
    }
}

//! Replayable point-sampled references on zero-power signal ports.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};

use super::{ElementError, Waveform, finite, positive};
use crate::ledger::Sum;

/// Accepted clock of a [`SignalSource`], including its compensation.
///
/// The default starts at zero. Like a physical source's clock, this state
/// belongs to the System and rewinds with rejected steps and checkpoints.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SignalClock(Sum);

impl SignalClock {
    /// The accepted source time [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.0.value()
    }

    fn validate(&self) -> Result<(), ElementError> {
        for word in self.0.to_words() {
            finite("signal clock word [s]", f64::from_bits(word))?;
        }
        finite("signal clock [s]", self.time_s())?;
        if self.time_s() < 0.0 {
            return Err(ElementError::Parameter {
                name: "signal clock [s]",
                value: self.time_s(),
            });
        }
        Ok(())
    }
}

/// One or more authored waveforms, evaluated at the beginning of each step.
///
/// There are no inputs. Signal port `output` has one coefficient per waveform.
/// Values are point samples `waveform.value_at(old.time_s())`, never averages
/// over a trial interval: a [`super::SampledController`] therefore receives the
/// same reference at a sample boundary regardless of physical substep size.
/// The signal carries no work. A [`super::ControlledSource`] accounts for any
/// physical energy an actuator commanded by it delivers or recovers.
///
/// Table knots are declared events, so the kernel does not step across them.
/// Sine waves do not prescribe a sampling rate. Resolution, signal units and
/// anti-alias filtering remain authored choices. The waveforms and their order
/// are bound into the problem identity, and the clock has a native state codec.
#[derive(Debug, Clone)]
pub struct SignalSource {
    waveforms: Vec<Waveform>,
    ports: [PortSpec; 1],
}

enum PreparedStep {
    Advance(SignalClock),
    Shorten(f64),
}

impl SignalSource {
    /// Build a nonempty vector reference using already admitted waveforms.
    pub fn new(waveforms: Vec<Waveform>) -> Result<Self, ElementError> {
        let width = NonZeroUsize::new(waveforms.len()).ok_or(ElementError::Parameter {
            name: "signal source width",
            value: 0.0,
        })?;
        Ok(Self {
            waveforms,
            ports: [PortSpec::signal_out(
                StableId::new("output").expect("static signal name"),
                width,
            )],
        })
    }

    fn validate(&self, old: &SignalClock, dt_s: f64, inputs: &[f64]) -> Result<(), ElementError> {
        old.validate()?;
        positive("signal source timestep [s]", dt_s)?;
        extent("signal source input width", inputs.len(), 0)
    }

    fn value(&self, index: usize, old: &SignalClock) -> Result<f64, ElementError> {
        finite(
            "signal source output",
            self.waveforms[index].value_at(old.time_s()),
        )
    }

    // A tangent describes this same admitted interval and point sample. A
    // zero-input map still has no derivative where its primal cannot run.
    fn prepare_step(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<PreparedStep, ElementError> {
        cx.checkpoint()?;
        self.validate(old, dt_s, inputs)?;
        if let Some(event_s) = self.next_event_s(old).filter(|event_s| dt_s > *event_s) {
            return Ok(PreparedStep::Shorten(event_s));
        }
        let mut next = *old;
        next.0.add(dt_s);
        next.validate()?;
        if next.time_s() <= old.time_s() {
            return Err(ElementError::Parameter {
                name: "signal source clock advance [s]",
                value: dt_s,
            });
        }
        for index in 0..self.waveforms.len() {
            if index % 256 == 0 {
                cx.checkpoint()?;
            }
            self.value(index, old)?;
        }
        cx.checkpoint()?;
        Ok(PreparedStep::Advance(next))
    }

    fn admit_stage(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SignalClock, ElementError> {
        match self.prepare_step(old, dt_s, inputs, cx)? {
            PreparedStep::Advance(next) => Ok(next),
            PreparedStep::Shorten(_) => Err(ElementError::Parameter {
                name: "signal source step crosses next event [s]",
                value: dt_s,
            }),
        }
    }

    fn publish_outputs(&self, old: &SignalClock, outputs: &mut [f64]) {
        for (index, output) in outputs.iter_mut().enumerate() {
            *output = self.waveforms[index].value_at(old.time_s());
        }
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

impl Component for SignalSource {
    type State = SignalClock;
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::SignalSource",
        state_family: "signal-clock",
        schema_grammar: "clock_s: finite nonnegative Sum running value and compensation",
        codec_grammar: "two u64 words: Sum::to_words of clock_s",
        codec_version: 1,
    });

    fn encode_state(&self, state: &SignalClock, words: &mut Vec<u64>) -> bool {
        if state.validate().is_err() {
            return false;
        }
        words.extend(state.0.to_words());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<SignalClock> {
        let words: [u64; 2] = words.try_into().ok()?;
        let state = SignalClock(Sum::from_words(words));
        state.validate().ok()?;
        Some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.waveforms.len() as u64);
        for waveform in &self.waveforms {
            waveform.identity_words(words);
        }
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // The map has no inputs, so its tangent is always the same empty map.
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &SignalClock) -> Result<f64, ElementError> {
        state.validate()?;
        Ok(0.0)
    }

    fn limit_dt_s(&self, state: &SignalClock, proposed_s: f64) -> Result<f64, ElementError> {
        state.validate()?;
        let limit = self.waveforms.iter().fold(proposed_s, |limit, waveform| {
            waveform.limit_dt_s(state.0, limit)
        });
        Ok(limit)
    }

    fn next_event_s(&self, state: &SignalClock) -> Option<f64> {
        if state.validate().is_err() {
            return Some(f64::NAN);
        }
        self.waveforms
            .iter()
            .filter_map(|waveform| waveform.next_event_s(state.0))
            .reduce(f64::min)
    }

    fn propose(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SignalClock,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        extent(
            "signal source output width",
            outputs.len(),
            self.waveforms.len(),
        )?;
        let next = match self.prepare_step(old, dt_s, inputs, cx)? {
            PreparedStep::Advance(next) => next,
            PreparedStep::Shorten(dt_s) => return Ok(Proposal::Shorten { dt_s }),
        };
        // All arithmetic has been admitted; publication evaluates those same
        // deterministic values, with no allocation or second state.
        self.publish_outputs(old, outputs);
        *new = next;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

// The clock is accepted bookkeeping, not an internal solve unknown. Reference
// values remain frozen at the accepted start time during every stage trial.
impl StageResidual for SignalSource {
    fn unknowns(&self) -> usize {
        0
    }

    fn initial_guess(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), ElementError> {
        extent("signal source stage width", z.len(), 0)?;
        self.validate(old, dt_s, inputs)
    }

    fn residual(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("signal source stage width", z.len(), 0)?;
        extent("signal source residual width", r.len(), 0)?;
        extent(
            "signal source output width",
            outputs.len(),
            self.waveforms.len(),
        )?;
        self.admit_stage(old, dt_s, inputs, cx)?;
        self.publish_outputs(old, outputs);
        Ok(())
    }

    fn jvp_z(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("signal source stage width", z.len(), 0)?;
        extent("signal source stage direction width", dz.len(), 0)?;
        extent("signal source residual direction width", dr.len(), 0)?;
        extent(
            "signal source output direction width",
            d_outputs.len(),
            self.waveforms.len(),
        )?;
        self.admit_stage(old, dt_s, inputs, cx)?;
        d_outputs.fill(0.0);
        Ok(())
    }

    fn vjp_z(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("signal source stage width", z.len(), 0)?;
        extent("signal source residual covector width", bar_r.len(), 0)?;
        extent("signal source stage covector width", bar_z.len(), 0)?;
        extent(
            "signal source output covector width",
            bar_outputs.len(),
            self.waveforms.len(),
        )?;
        super::inputs_finite(bar_outputs)?;
        self.admit_stage(old, dt_s, inputs, cx).map(|_| ())
    }

    fn accept(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut SignalClock,
        cx: &Cx<'_>,
    ) -> Result<Receipt, ElementError> {
        cx.checkpoint()?;
        extent("signal source stage width", z.len(), 0)?;
        let next = self.admit_stage(old, dt_s, inputs, cx)?;
        *new = next;
        Ok(Receipt::lossless(0.0, 0.0))
    }
}

impl PortTangent for SignalSource {
    fn jvp(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("signal source input direction width", d_inputs.len(), 0)?;
        extent(
            "signal source output direction width",
            d_outputs.len(),
            self.waveforms.len(),
        )?;
        self.admit_stage(old, dt_s, inputs, cx)?;
        d_outputs.fill(0.0);
        Ok(())
    }

    fn vjp(
        &self,
        old: &SignalClock,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent(
            "signal source output covector width",
            bar_outputs.len(),
            self.waveforms.len(),
        )?;
        extent("signal source input covector width", bar_inputs.len(), 0)?;
        super::inputs_finite(bar_outputs)?;
        self.admit_stage(old, dt_s, inputs, cx).map(|_| ())
    }
}

//! Explicit bounds on authored control signals.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};

use super::{ElementError, finite, positive};

/// The stateless signal relation `output[i] = clamp(input[i], lower[i], upper[i])`.
///
/// Both `input` and `output` are vector signal ports, with one pair of authored
/// bounds per coordinate. Equal bounds hold a coordinate constant. This is an
/// explicit control law: it neither changes a physical state nor relaxes an
/// actuator's admission rules. No signal carries energy, mass or dissipation.
///
/// Add with [`crate::SystemBuilder::add_with_tangent`] when the signal is part
/// of a feedback loop. The tangent is diagonal, one strictly inside the bounds
/// and zero outside. At either bound it selects zero from the generalized
/// Jacobian; a classical derivative does not exist there unless both bounds
/// coincide. The same selection is used by the transpose. This piecewise law
/// supplies no smoothing, event location, anti-windup or stability guarantee.
#[derive(Debug, Clone)]
pub struct SignalLimiter {
    lower: Vec<f64>,
    upper: Vec<f64>,
    ports: [PortSpec; 2],
}

impl SignalLimiter {
    /// Construct a limiter with nonempty, equal-width, finite bound vectors.
    ///
    /// # Errors
    /// Refuses an empty or mismatched extent, nonfinite bounds, or a lower
    /// bound greater than its corresponding upper bound.
    pub fn new(lower: Vec<f64>, upper: Vec<f64>) -> Result<Self, ElementError> {
        let width = NonZeroUsize::new(lower.len()).ok_or(ElementError::Parameter {
            name: "signal limiter width",
            value: 0.0,
        })?;
        extent("signal limiter upper width", upper.len(), lower.len())?;
        for (&low, &high) in lower.iter().zip(&upper) {
            finite("signal limiter lower bound", low)?;
            finite("signal limiter upper bound", high)?;
            if low > high {
                return Err(ElementError::Parameter {
                    name: "signal limiter lower bound exceeds upper bound",
                    value: low,
                });
            }
        }
        Ok(Self {
            lower,
            upper,
            ports: [
                PortSpec::signal_in(StableId::new("input").expect("static signal name"), width),
                PortSpec::signal_out(StableId::new("output").expect("static signal name"), width),
            ],
        })
    }

    fn validate(&self, dt_s: f64, inputs: &[f64], cx: &Cx<'_>) -> Result<(), ElementError> {
        cx.checkpoint()?;
        positive("signal limiter timestep [s]", dt_s)?;
        extent("signal limiter input width", inputs.len(), self.lower.len())?;
        finite_vector(inputs, cx)
    }

    fn tangent(
        &self,
        dt_s: f64,
        inputs: &[f64],
        seed: &[f64],
        result: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.validate(dt_s, inputs, cx)?;
        extent(
            "signal limiter tangent seed width",
            seed.len(),
            self.lower.len(),
        )?;
        extent(
            "signal limiter tangent result width",
            result.len(),
            self.lower.len(),
        )?;
        finite_vector(seed, cx)?;
        for (index, output) in result.iter_mut().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint()?;
            }
            *output = if inputs[index] > self.lower[index] && inputs[index] < self.upper[index] {
                seed[index]
            } else {
                0.0
            };
        }
        Ok(())
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

fn finite_vector(values: &[f64], cx: &Cx<'_>) -> Result<(), ElementError> {
    for (index, &value) in values.iter().enumerate() {
        if index % 256 == 0 {
            cx.checkpoint()?;
        }
        finite("signal limiter input or tangent seed", value)?;
    }
    Ok(())
}

impl Component for SignalLimiter {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::SignalLimiter",
        state_family: "no-state",
        schema_grammar: "(): an authored signal limiter stores nothing",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.lower.len() as u64);
        words.extend(self.lower.iter().map(|value| value.to_bits()));
        words.extend(self.upper.iter().map(|value| value.to_bits()));
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
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.validate(dt_s, inputs, cx)?;
        extent(
            "signal limiter output width",
            outputs.len(),
            self.lower.len(),
        )?;
        for (index, output) in outputs.iter_mut().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint()?;
            }
            *output = inputs[index].clamp(self.lower[index], self.upper[index]);
        }
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

impl PortTangent for SignalLimiter {
    fn jvp(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.tangent(dt_s, inputs, d_inputs, d_outputs, cx)
    }

    fn vjp(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.tangent(dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

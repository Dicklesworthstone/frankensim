//! Authored static signal maps for composing feedback without a new engine.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};

use super::{ElementError, finite, inputs_finite, positive};

/// The stateless signal relation `output = bias + matrix * input`.
///
/// Each input is its own scalar signal port, `input0`, `input1`, and so on,
/// so independent references and observations can be wired by the ordinary
/// signal connection. Port `output` is a vector of `output_width` values.
/// Matrix entries are row-major, and the contiguous input vector follows
/// input-port order. Negative and zero gains and rectangular matrices are
/// admitted. No signal carries physical power, storage or dissipation.
///
/// The caller authors the gains and bias in the input and output coordinates'
/// required units. This law provides no unit inference, saturation, delay,
/// dynamic controller state or stability claim. Feedback from
/// [`super::Mass::publishing`] uses the accepted state at the start of each
/// step; connect its output to a [`super::ControlledSource`] to exchange
/// physical work with a device.
#[derive(Debug, Clone)]
pub struct AffineSignal {
    input_width: usize,
    output_width: usize,
    matrix: Vec<f64>,
    bias: Vec<f64>,
    ports: Vec<PortSpec>,
}

impl AffineSignal {
    /// Build a finite dense map from `input_width` scalar signals to one
    /// `output_width`-wide signal, with exactly one bias per output.
    ///
    /// # Errors
    /// Refuses zero dimensions, overflowing or mismatched matrix extents,
    /// a mismatched bias extent, or any nonfinite weight or bias.
    pub fn new(
        input_width: usize,
        output_width: usize,
        matrix: Vec<f64>,
        bias: Vec<f64>,
    ) -> Result<Self, ElementError> {
        let output = NonZeroUsize::new(output_width).ok_or(ElementError::Parameter {
            name: "affine signal output width",
            value: output_width as f64,
        })?;
        if input_width == 0 {
            return Err(ElementError::Parameter {
                name: "affine signal input width",
                value: 0.0,
            });
        }
        let entries = input_width
            .checked_mul(output_width)
            .ok_or(ElementError::Parameter {
                name: "affine signal matrix extent",
                value: matrix.len() as f64,
            })?;
        extent("affine signal matrix entries", matrix.len(), entries)?;
        extent("affine signal bias width", bias.len(), output_width)?;
        for &value in &matrix {
            finite("affine signal matrix entry", value)?;
        }
        for &value in &bias {
            finite("affine signal bias", value)?;
        }
        let mut ports = Vec::with_capacity(input_width + 1);
        for index in 0..input_width {
            ports.push(PortSpec::signal_in(
                StableId::new(format!("input{index}")).expect("generated signal name"),
                NonZeroUsize::MIN,
            ));
        }
        ports.push(PortSpec::signal_out(
            StableId::new("output").expect("static signal name"),
            output,
        ));
        Ok(Self {
            input_width,
            output_width,
            matrix,
            bias,
            ports,
        })
    }

    fn validate(&self, dt_s: f64, inputs: &[f64]) -> Result<(), ElementError> {
        positive("affine signal timestep [s]", dt_s)?;
        extent("affine signal input width", inputs.len(), self.input_width)?;
        inputs_finite(inputs)
    }

    fn row(
        &self,
        row: usize,
        inputs: &[f64],
        biased: bool,
        cx: &Cx<'_>,
    ) -> Result<f64, ElementError> {
        let mut value = if biased { self.bias[row] } else { 0.0 };
        for (column, (&weight, &input)) in self.matrix
            [row * self.input_width..(row + 1) * self.input_width]
            .iter()
            .zip(inputs)
            .enumerate()
        {
            if column % 256 == 0 {
                cx.checkpoint()?;
            }
            value += weight * input;
        }
        finite("affine signal output", value)
    }

    fn check_response(&self, inputs: &[f64], cx: &Cx<'_>) -> Result<(), ElementError> {
        for row in 0..self.output_width {
            self.row(row, inputs, true, cx)?;
        }
        Ok(())
    }

    fn evaluate(
        &self,
        dt_s: f64,
        inputs: &[f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        self.validate(dt_s, inputs)?;
        extent(
            "affine signal output width",
            outputs.len(),
            self.output_width,
        )?;
        for (row, output) in outputs.iter_mut().enumerate() {
            *output = self.row(row, inputs, true, cx)?;
        }
        cx.checkpoint()?;
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

impl Component for AffineSignal {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::AffineSignal",
        state_family: "no-state",
        schema_grammar: "(): a static affine signal map stores nothing",
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
        words.extend([self.input_width as u64, self.output_width as u64]);
        words.extend(self.matrix.iter().map(|value| value.to_bits()));
        words.extend(self.bias.iter().map(|value| value.to_bits()));
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
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
        self.evaluate(dt_s, inputs, outputs, cx)?;
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

// The affine map has no internal state unknown. Its input/output equations
// still participate in a coupled algebraic loop through the retained inputs.
impl StageResidual for AffineSignal {
    fn unknowns(&self) -> usize {
        0
    }

    fn initial_guess(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), ElementError> {
        extent("affine signal stage width", z.len(), 0)?;
        self.validate(dt_s, inputs)
    }

    fn residual(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("affine signal stage width", z.len(), 0)?;
        extent("affine signal residual width", r.len(), 0)?;
        self.evaluate(dt_s, inputs, outputs, cx)
    }

    fn jvp_z(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("affine signal stage width", z.len(), 0)?;
        extent("affine signal stage direction width", dz.len(), 0)?;
        extent("affine signal residual direction width", dr.len(), 0)?;
        extent(
            "affine signal output direction width",
            d_outputs.len(),
            self.output_width,
        )?;
        self.validate(dt_s, inputs)?;
        self.check_response(inputs, cx)?;
        cx.checkpoint()?;
        d_outputs.fill(0.0);
        Ok(())
    }

    fn vjp_z(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        extent("affine signal stage width", z.len(), 0)?;
        extent("affine signal residual covector width", bar_r.len(), 0)?;
        extent(
            "affine signal output covector width",
            bar_outputs.len(),
            self.output_width,
        )?;
        extent("affine signal stage covector width", bar_z.len(), 0)?;
        self.validate(dt_s, inputs)?;
        inputs_finite(bar_outputs)?;
        self.check_response(inputs, cx)?;
        cx.checkpoint()?;
        Ok(())
    }

    fn accept(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        _new: &mut (),
        cx: &Cx<'_>,
    ) -> Result<Receipt, ElementError> {
        cx.checkpoint()?;
        extent("affine signal stage width", z.len(), 0)?;
        self.validate(dt_s, inputs)?;
        self.check_response(inputs, cx)?;
        cx.checkpoint()?;
        Ok(Receipt::lossless(0.0, 0.0))
    }
}

impl PortTangent for AffineSignal {
    fn jvp(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        self.validate(dt_s, inputs)?;
        extent(
            "affine signal input direction width",
            d_inputs.len(),
            self.input_width,
        )?;
        extent(
            "affine signal output direction width",
            d_outputs.len(),
            self.output_width,
        )?;
        inputs_finite(d_inputs)?;
        self.check_response(inputs, cx)?;
        for (row, output) in d_outputs.iter_mut().enumerate() {
            *output = self.row(row, d_inputs, false, cx)?;
        }
        cx.checkpoint()?;
        Ok(())
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
        cx.checkpoint()?;
        self.validate(dt_s, inputs)?;
        extent(
            "affine signal output covector width",
            bar_outputs.len(),
            self.output_width,
        )?;
        extent(
            "affine signal input covector width",
            bar_inputs.len(),
            self.input_width,
        )?;
        inputs_finite(bar_outputs)?;
        self.check_response(inputs, cx)?;
        for (column, output) in bar_inputs.iter_mut().enumerate() {
            let mut value = 0.0;
            for (row, &input) in bar_outputs.iter().enumerate() {
                if row % 256 == 0 {
                    cx.checkpoint()?;
                }
                value += self.matrix[row * self.input_width + column] * input;
            }
            *output = finite("affine signal input covector", value)?;
        }
        cx.checkpoint()?;
        Ok(())
    }
}

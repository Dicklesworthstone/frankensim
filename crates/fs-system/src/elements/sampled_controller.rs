//! Sampled, held discrete state-space control on zero-power signals.
//!
//! The authored matrices are already a discrete-time law, not a continuous
//! model to be integrated here. At sample `k`, from the state `x[k]`,
//! `y[k] = C*x[k] + D*u[k]` and `x[k+1] = A*x[k] + B*u[k]`. The first
//! sample is at zero; subsequent samples are at the representable times
//! `k * sample_period_s`. The output is held between those events even if
//! the physics takes smaller steps. A refused coupled step accepts neither
//! the recurrence nor the held command.
//!
//! Inputs `input0`, `input1`, ... are scalar signals; `output` is one vector
//! signal. Their sources must supply point samples at the start of a sample
//! interval, not an interval mean. The controller cannot reconstruct a point
//! sample from an average. Gains, units, sample period and initial memory are
//! authored; there is no gain synthesis, saturation or stability claim.

use core::cell::RefCell;
use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};

use crate::ledger::Sum;

use super::{ElementError, finite, inputs_finite, positive};

/// Memory, held output and clock of a [`SampledController`].
///
/// After the first accepted part of sample interval `k`, `coordinates()`
/// holds `x[k+1]` and `held_output()` holds `y[k]`. The recurrence does not
/// run again until the following sample boundary. None of this memory is
/// physical stored energy.
#[derive(Debug, Clone, PartialEq)]
pub struct ControllerState {
    coordinates: Vec<f64>,
    held_output: Vec<f64>,
    clock_s: Sum,
    samples: u64,
}

impl ControllerState {
    /// An unsampled controller at time zero with the supplied `x[0]`.
    /// The held output is zero until the first accepted proposal samples it.
    ///
    /// # Errors
    /// The wrong state width or a nonfinite initial coordinate.
    pub fn new(controller: &SampledController, initial: Vec<f64>) -> Result<Self, ElementError> {
        extent(
            "controller initial state width",
            initial.len(),
            controller.state_dim,
        )?;
        inputs_finite(&initial)?;
        Ok(Self {
            coordinates: initial,
            held_output: vec![0.0; controller.output_dim],
            clock_s: Sum::default(),
            samples: 0,
        })
    }

    /// Accepted discrete memory, in the authored state-coordinate order.
    #[must_use]
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }

    /// The most recently sampled output, held for its entire period.
    #[must_use]
    pub fn held_output(&self) -> &[f64] {
        &self.held_output
    }

    /// Accepted physical time [s], on the same compensated clock as the system.
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.clock_s.value()
    }

    /// Number of accepted samples, independent of physical substep count.
    #[must_use]
    pub fn samples(&self) -> u64 {
        self.samples
    }
}

/// An authored discrete state-space controller with an explicit sample clock.
///
/// Construction prepares all temporary storage. Proposals and tangent calls
/// allocate nothing when the candidate state has the admitted widths (as the
/// system's candidate slot does). Scratch is not accepted state, is not exposed
/// by `Debug`, and never contributes to the recurrence or a snapshot.
#[derive(Clone)]
pub struct SampledController {
    sample_period_s: f64,
    state_dim: usize,
    input_dim: usize,
    output_dim: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    ports: Vec<PortSpec>,
    scratch: RefCell<Vec<f64>>,
}

impl fmt::Debug for SampledController {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SampledController")
            .field("sample_period_s", &self.sample_period_s)
            .field(
                "dimensions",
                &(self.state_dim, self.input_dim, self.output_dim),
            )
            .field("a", &self.a)
            .field("b", &self.b)
            .field("c", &self.c)
            .field("d", &self.d)
            .finish()
    }
}

#[derive(Clone, Copy)]
struct Phase {
    sample: bool,
    samples: u64,
    boundary_s: f64,
    until_s: f64,
}

impl SampledController {
    /// Build `x_next = A*x + B*u`, `y = C*x + D*u`, with row-major matrices.
    ///
    /// For dimensions `(n, m, p)`, the matrix extents are `n*n`, `n*m`,
    /// `p*n`, and `p*m`, respectively. `n` may be zero for a sampled gain;
    /// input and output widths must be positive. Scalar input ports are
    /// named `input0` through `input{m-1}`; the vector output is `output`.
    ///
    /// # Errors
    /// A nonpositive/nonfinite period, zero input/output width, an overflowing
    /// extent, a wrong matrix shape or any nonfinite matrix coefficient.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sample_period_s: f64,
        state_dim: usize,
        input_dim: usize,
        output_dim: usize,
        a: Vec<f64>,
        b: Vec<f64>,
        c: Vec<f64>,
        d: Vec<f64>,
    ) -> Result<Self, ElementError> {
        positive("controller sample period [s]", sample_period_s)?;
        let output_width = NonZeroUsize::new(output_dim)
            .ok_or_else(|| parameter("controller output width", 0.0))?;
        if input_dim == 0 {
            return Err(parameter("controller input width", 0.0));
        }
        for (name, matrix, rows, columns) in [
            ("controller A shape", &a, state_dim, state_dim),
            ("controller B shape", &b, state_dim, input_dim),
            ("controller C shape", &c, output_dim, state_dim),
            ("controller D shape", &d, output_dim, input_dim),
        ] {
            let expected = rows
                .checked_mul(columns)
                .ok_or_else(|| parameter(name, f64::INFINITY))?;
            extent(name, matrix.len(), expected)?;
            inputs_finite(matrix)?;
        }
        let scratch_len = state_dim
            .checked_add(output_dim)
            .and_then(|n| n.checked_add(3))
            .ok_or_else(|| parameter("controller state extent", f64::INFINITY))?;
        let port_count = input_dim
            .checked_add(1)
            .ok_or_else(|| parameter("controller port count", f64::INFINITY))?;
        let mut ports = Vec::with_capacity(port_count);
        for input in 0..input_dim {
            let id = StableId::new(format!("input{input}"))
                .map_err(|_| parameter("controller input identifier", input as f64))?;
            ports.push(PortSpec::signal_in(id, NonZeroUsize::MIN));
        }
        let id =
            StableId::new("output").map_err(|_| parameter("controller output identifier", 0.0))?;
        ports.push(PortSpec::signal_out(id, output_width));
        Ok(Self {
            sample_period_s,
            state_dim,
            input_dim,
            output_dim,
            a,
            b,
            c,
            d,
            ports,
            scratch: RefCell::new(vec![0.0; scratch_len.max(input_dim)]),
        })
    }

    /// The fixed interval between samples [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 {
        self.sample_period_s
    }

    /// Number of discrete memory coordinates.
    #[must_use]
    pub fn state_dimension(&self) -> usize {
        self.state_dim
    }

    /// Number of scalar input ports.
    #[must_use]
    pub fn input_dimension(&self) -> usize {
        self.input_dim
    }

    /// Width of the output signal.
    #[must_use]
    pub fn output_dimension(&self) -> usize {
        self.output_dim
    }

    // Integer multiples avoid accumulating period-rounding drift. Refuse
    // before distinct samples collapse to one representable event time.
    fn sample_time(&self, index: u64) -> Result<f64, ElementError> {
        if index > (1_u64 << 53) {
            return Err(parameter(
                "controller sample counter resolution",
                index as f64,
            ));
        }
        let time_s = finite(
            "controller sample time [s]",
            index as f64 * self.sample_period_s,
        )?;
        if index > 0 && time_s <= (index - 1) as f64 * self.sample_period_s {
            return Err(parameter("controller sample clock resolution [s]", time_s));
        }
        Ok(time_s)
    }

    fn validate_state(&self, state: &ControllerState) -> Result<(), ElementError> {
        extent(
            "controller state width",
            state.coordinates.len(),
            self.state_dim,
        )?;
        extent(
            "controller held output width",
            state.held_output.len(),
            self.output_dim,
        )?;
        inputs_finite(&state.coordinates)?;
        inputs_finite(&state.held_output)?;
        for word in state.clock_s.to_words() {
            finite("controller clock word [s]", f64::from_bits(word))?;
        }
        let now = finite("controller clock [s]", state.time_s())?;
        if state.samples == 0 {
            if state.clock_s.to_words() != [0, 0]
                || state.held_output.iter().any(|value| value.to_bits() != 0)
            {
                return Err(parameter("controller unsampled state", now));
            }
        } else {
            let start = self.sample_time(state.samples - 1)?;
            let end = self.sample_time(state.samples)?;
            if !(now > start && now <= end) {
                return Err(parameter(
                    "controller clock outside sample interval [s]",
                    now,
                ));
            }
        }
        Ok(())
    }

    fn phase(&self, state: &ControllerState) -> Result<Phase, ElementError> {
        self.validate_state(state)?;
        let sample = state.samples == 0 || state.time_s() == self.sample_time(state.samples)?;
        let samples = state.samples + u64::from(sample);
        let boundary_s = self.sample_time(samples)?;
        let until_s = positive(
            "controller next sample interval [s]",
            state.clock_s.interval_to(boundary_s),
        )?;
        Ok(Phase {
            sample,
            samples,
            boundary_s,
            until_s,
        })
    }

    fn validate_call(
        &self,
        old: &ControllerState,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<Phase, ElementError> {
        positive("controller timestep [s]", dt_s)?;
        extent("controller input width", inputs.len(), self.input_dim)?;
        inputs_finite(inputs)?;
        self.phase(old)
    }

    fn advance_clock(old: &ControllerState, dt_s: f64, phase: Phase) -> Result<Sum, ElementError> {
        let mut clock_s = old.clock_s;
        clock_s.add(dt_s);
        let now = clock_s.value();
        if !(now.is_finite() && now > old.time_s() && now <= phase.boundary_s) {
            return Err(parameter(
                "controller clock cannot advance inside sample interval [s]",
                now,
            ));
        }
        Ok(clock_s)
    }

    fn sample(
        &self,
        old: &ControllerState,
        inputs: &[f64],
        scratch: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        for row in 0..self.state_dim {
            let ax = dot(
                &self.a[row * self.state_dim..][..self.state_dim],
                &old.coordinates,
                0.0,
                cx,
            )?;
            scratch[row] = dot(
                &self.b[row * self.input_dim..][..self.input_dim],
                inputs,
                ax,
                cx,
            )?;
        }
        for row in 0..self.output_dim {
            let cx_value = dot(
                &self.c[row * self.state_dim..][..self.state_dim],
                &old.coordinates,
                0.0,
                cx,
            )?;
            scratch[self.state_dim + row] = dot(
                &self.d[row * self.input_dim..][..self.input_dim],
                inputs,
                cx_value,
                cx,
            )?;
        }
        Ok(())
    }
}

fn parameter(name: &'static str, value: f64) -> ElementError {
    ElementError::Parameter { name, value }
}

fn extent(name: &'static str, actual: usize, expected: usize) -> Result<(), ElementError> {
    if actual == expected {
        Ok(())
    } else {
        Err(parameter(name, actual as f64))
    }
}

fn dot(row: &[f64], values: &[f64], mut sum: f64, cx: &Cx<'_>) -> Result<f64, ElementError> {
    cx.checkpoint()?;
    for (index, (&coefficient, &value)) in row.iter().zip(values).enumerate() {
        if index % 64 == 0 {
            cx.checkpoint()?;
        }
        let term = finite("controller matrix product", coefficient * value)?;
        sum = finite("controller matrix sum", sum + term)?;
    }
    Ok(sum)
}

impl Component for SampledController {
    type State = ControllerState;
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::SampledController",
        state_family: "sampled-state-space-memory",
        schema_grammar: "sample count u64; compensated clock; finite discrete memory and held output with admitted dimensions and clock inside its sample interval",
        codec_grammar: "sample count u64; clock Sum::to_words two u64; state_dim f64 bits; output_dim f64 bits",
        codec_version: 1,
    });

    fn encode_state(&self, state: &ControllerState, words: &mut Vec<u64>) -> bool {
        if self.validate_state(state).is_err() {
            return false;
        }
        words.push(state.samples);
        words.extend(state.clock_s.to_words());
        words.extend(
            state
                .coordinates
                .iter()
                .chain(&state.held_output)
                .map(|x| x.to_bits()),
        );
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<ControllerState> {
        if words.len() != 3 + self.state_dim + self.output_dim {
            return None;
        }
        let state = ControllerState {
            samples: words[0],
            clock_s: Sum::from_words([words[1], words[2]]),
            coordinates: words[3..3 + self.state_dim]
                .iter()
                .map(|&x| f64::from_bits(x))
                .collect(),
            held_output: words[3 + self.state_dim..]
                .iter()
                .map(|&x| f64::from_bits(x))
                .collect(),
        };
        self.validate_state(&state).ok().map(|()| state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.sample_period_s.to_bits(),
            self.state_dim as u64,
            self.input_dim as u64,
            self.output_dim as u64,
        ]);
        words.extend(
            self.a
                .iter()
                .chain(&self.b)
                .chain(&self.c)
                .chain(&self.d)
                .map(|x| x.to_bits()),
        );
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // D at a sample, zero during its hold: this slope depends on the
        // accepted clock even when two consecutive physical steps are equal.
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &ControllerState) -> Result<f64, ElementError> {
        self.validate_state(state)?;
        Ok(0.0)
    }

    fn next_event_s(&self, state: &ControllerState) -> Option<f64> {
        // The fallible negotiation and proposal paths reject invalid state;
        // this infallible optional query does not encode an error as NaN.
        self.phase(state).ok().map(|phase| phase.until_s)
    }

    fn limit_dt_s(&self, state: &ControllerState, proposed_s: f64) -> Result<f64, ElementError> {
        positive("controller proposed timestep [s]", proposed_s)?;
        let phase = self.phase(state)?;
        Ok(state
            .clock_s
            .limit_before_event(proposed_s, phase.until_s)
            .min(phase.until_s))
    }

    fn propose(
        &self,
        old: &ControllerState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut ControllerState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        let phase = self.validate_call(old, dt_s, inputs)?;
        extent("controller output width", outputs.len(), self.output_dim)?;
        if dt_s > phase.until_s {
            return Ok(Proposal::Shorten {
                dt_s: phase.until_s,
            });
        }
        let clock_s = Self::advance_clock(old, dt_s, phase)?;
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| parameter("controller workspace already borrowed", 1.0))?;
        if phase.sample {
            self.sample(old, inputs, &mut scratch, cx)?;
        } else {
            scratch[..self.state_dim].copy_from_slice(&old.coordinates);
            scratch[self.state_dim..self.state_dim + self.output_dim]
                .copy_from_slice(&old.held_output);
        }
        cx.checkpoint()?;
        // All fallible work precedes publication, including the next memory
        // coordinate and every output; a late overflow publishes no prefix.
        new.coordinates.resize(self.state_dim, 0.0);
        new.coordinates.copy_from_slice(&scratch[..self.state_dim]);
        new.held_output.resize(self.output_dim, 0.0);
        new.held_output
            .copy_from_slice(&scratch[self.state_dim..self.state_dim + self.output_dim]);
        new.clock_s = clock_s;
        new.samples = phase.samples;
        outputs.copy_from_slice(&new.held_output);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

impl PortTangent for SampledController {
    fn jvp(
        &self,
        old: &ControllerState,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        let phase = self.validate_call(old, dt_s, inputs)?;
        extent(
            "controller input direction width",
            d_inputs.len(),
            self.input_dim,
        )?;
        extent(
            "controller output direction width",
            d_outputs.len(),
            self.output_dim,
        )?;
        inputs_finite(d_inputs)?;
        if dt_s > phase.until_s {
            return Err(parameter(
                "controller tangent step crosses sample boundary [s]",
                dt_s,
            ));
        }
        Self::advance_clock(old, dt_s, phase)?;
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| parameter("controller workspace already borrowed", 1.0))?;
        if phase.sample {
            self.sample(old, inputs, &mut scratch, cx)?;
            for row in 0..self.output_dim {
                scratch[row] = dot(
                    &self.d[row * self.input_dim..][..self.input_dim],
                    d_inputs,
                    0.0,
                    cx,
                )?;
            }
        } else {
            scratch[..self.output_dim].fill(0.0);
        }
        cx.checkpoint()?;
        d_outputs.copy_from_slice(&scratch[..self.output_dim]);
        Ok(())
    }

    fn vjp(
        &self,
        old: &ControllerState,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        cx.checkpoint()?;
        let phase = self.validate_call(old, dt_s, inputs)?;
        extent(
            "controller output covector width",
            bar_outputs.len(),
            self.output_dim,
        )?;
        extent(
            "controller input covector width",
            bar_inputs.len(),
            self.input_dim,
        )?;
        inputs_finite(bar_outputs)?;
        if dt_s > phase.until_s {
            return Err(parameter(
                "controller tangent step crosses sample boundary [s]",
                dt_s,
            ));
        }
        Self::advance_clock(old, dt_s, phase)?;
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| parameter("controller workspace already borrowed", 1.0))?;
        if phase.sample {
            self.sample(old, inputs, &mut scratch, cx)?;
            for column in 0..self.input_dim {
                let mut sum = 0.0;
                for (row, &bar) in bar_outputs.iter().enumerate() {
                    if row % 64 == 0 {
                        cx.checkpoint()?;
                    }
                    let term = finite(
                        "controller transpose product",
                        self.d[row * self.input_dim + column] * bar,
                    )?;
                    sum = finite("controller transpose sum", sum + term)?;
                }
                scratch[column] = sum;
            }
        } else {
            scratch[..self.input_dim].fill(0.0);
        }
        cx.checkpoint()?;
        bar_inputs.copy_from_slice(&scratch[..self.input_dim]);
        Ok(())
    }
}

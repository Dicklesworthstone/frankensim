//! Accepted-state observations on any physical Component, without another state owner.

use core::fmt;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;

use super::{
    Class, Component, ContractViolation, EndOfStep, LoadLine, PortLayout, PortRole, PortSpec,
    PortTangent, Proposal, Receipt, StageResidual,
};

/// A physical owner's refusal or a failure of its added state observations.
#[derive(Debug)]
pub enum ObservationError<E> {
    /// Original owner or observation callback error, retained without conversion.
    Owner(E),
    /// A caller supplied a vector with the wrong length.
    Contract(ContractViolation),
    /// An observation was not written or was not finite.
    InvalidSignal {
        /// Coefficient index among the appended observations.
        index: usize,
        /// Rejected value; an unwritten coefficient is NaN.
        value: f64,
    },
    /// Cancellation during observation evaluation or checking.
    Cancelled,
}

impl<E: fmt::Debug> fmt::Display for ObservationError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "observed component owner refused: {error:?}"),
            Self::Contract(error) => fmt::Display::fmt(error, f),
            Self::InvalidSignal { index, value } => {
                write!(f, "state observation {index} must be finite, got {value}")
            }
            Self::Cancelled => write!(f, "state observation cancelled"),
        }
    }
}

impl<E: fmt::Debug + 'static> std::error::Error for ObservationError<E> {}

/// Append pure, held start-of-step signal observations to an existing Component.
///
/// The callback reads the accepted `old` state, never the candidate or step-mean
/// ports, and writes every appended coefficient in port order. It must obey
/// [`Component::propose`]'s purity and bounded-cancellation rules. Observations
/// may be nonlinear functions of the state; their input derivative at fixed
/// `old` is zero. The owner's physics, state, receipts and optional solver tiers
/// are unchanged. This supplies explicit sampled feedback, not implicit
/// within-step material feedback or a higher-order feedback discretization.
///
/// Durable codecs are forwarded only after [`Self::with_observation_identity`]
/// identifies the complete observation law. No identity is inferred from a
/// callback address, type name or sample values.
pub struct StateSignals<C: Component, F> {
    inner: C,
    observe: F,
    ports: Vec<PortSpec>,
    inputs: usize,
    native_outputs: usize,
    outputs: usize,
    observation_identity: Option<[u8; 32]>,
}

impl<C, F> StateSignals<C, F>
where
    C: Component,
    F: Fn(&C, &C::State, &mut [f64], &Cx<'_>) -> Result<(), C::Error>,
{
    /// Append the supplied signal-output ports in their declared order.
    /// Original port indices and input layout are preserved. Duplicate IDs or
    /// any appended port other than `SignalOut` refuse admission.
    pub fn new(inner: C, signals: Vec<PortSpec>, observe: F) -> Result<Self, ContractViolation> {
        for signal in &signals {
            if signal.role() != PortRole::SignalOut {
                return Err(ContractViolation::PortRole {
                    port: signal.id().as_str().to_owned(),
                });
            }
        }
        let native = PortLayout::new(inner.ports())?;
        let mut ports = inner.ports().to_vec();
        ports.extend(signals);
        let layout = PortLayout::new(&ports)?;
        Ok(Self {
            inner,
            observe,
            ports,
            inputs: native.input_len(),
            native_outputs: native.output_len(),
            outputs: layout.output_len(),
            observation_identity: None,
        })
    }

    /// Declare the stable identity of the complete original observation law,
    /// including its parameters, units, support and refusal behavior. This is
    /// an owner declaration, not authentication of an arbitrary callback.
    #[must_use]
    pub fn with_observation_identity(mut self, identity: [u8; 32]) -> Self {
        self.observation_identity = Some(identity);
        self
    }

    /// The unchanged physical owner.
    #[must_use]
    pub fn inner(&self) -> &C {
        &self.inner
    }

    fn length(
        &self,
        vector: &'static str,
        expected: usize,
        actual: usize,
    ) -> Result<(), ObservationError<C::Error>> {
        if expected == actual {
            Ok(())
        } else {
            Err(ObservationError::Contract(
                ContractViolation::VectorLength {
                    vector,
                    expected,
                    actual,
                },
            ))
        }
    }

    fn fill(
        &self,
        outputs: &mut [f64],
        value: f64,
        cx: &Cx<'_>,
    ) -> Result<(), ObservationError<C::Error>> {
        cx.checkpoint().map_err(|_| ObservationError::Cancelled)?;
        for chunk in outputs.chunks_mut(256) {
            chunk.fill(value);
            cx.checkpoint().map_err(|_| ObservationError::Cancelled)?;
        }
        Ok(())
    }

    fn sample(
        &self,
        old: &C::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ObservationError<C::Error>> {
        self.fill(outputs, f64::NAN, cx)?;
        (self.observe)(&self.inner, old, outputs, cx).map_err(ObservationError::Owner)?;
        for (index, &value) in outputs.iter().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint().map_err(|_| ObservationError::Cancelled)?;
            }
            if !value.is_finite() {
                return Err(ObservationError::InvalidSignal { index, value });
            }
        }
        cx.checkpoint().map_err(|_| ObservationError::Cancelled)?;
        Ok(())
    }

    fn observations(
        &self,
        proposal: Proposal,
        old: &C::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ObservationError<C::Error>> {
        if matches!(proposal, Proposal::Step(_)) {
            self.sample(old, outputs, cx)?;
        }
        Ok(proposal)
    }
}

impl<C, F> Component for StateSignals<C, F>
where
    C: Component,
    F: Fn(&C, &C::State, &mut [f64], &Cx<'_>) -> Result<(), C::Error>,
{
    type State = C::State;
    type Error = ObservationError<C::Error>;
    const STATE_CHARTER: Option<StateIdentityCharterV2> = C::STATE_CHARTER;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        self.inner.class()
    }
    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        self.inner
            .stored_energy_j(state)
            .map_err(ObservationError::Owner)
    }
    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        self.inner
            .stored_mass_kg(state)
            .map_err(ObservationError::Owner)
    }
    fn stored_entropy_j_per_k(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        self.inner
            .stored_entropy_j_per_k(state)
            .map_err(ObservationError::Owner)
    }
    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner
            .stored_carried_quantity(state, carried)
            .map_err(ObservationError::Owner)
    }
    fn limit_dt_s(&self, state: &Self::State, proposed_s: f64) -> Result<f64, Self::Error> {
        self.inner
            .limit_dt_s(state, proposed_s)
            .map_err(ObservationError::Owner)
    }
    fn next_event_s(&self, state: &Self::State) -> Option<f64> {
        self.inner.next_event_s(state)
    }
    fn end_of_step(&self, port: usize) -> EndOfStep {
        if port < self.inner.ports().len() {
            self.inner.end_of_step(port)
        } else {
            EndOfStep::Untold
        }
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
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("outputs", self.outputs, outputs.len())?;
        let (native, observed) = outputs.split_at_mut(self.native_outputs);
        let proposal = self
            .inner
            .propose(old, dt_s, inputs, new, native, cx)
            .map_err(ObservationError::Owner)?;
        self.observations(proposal, old, observed, cx)
    }
    fn admit(
        &self,
        old: &Self::State,
        new: &Self::State,
        receipt: &Receipt,
    ) -> Result<(), Self::Error> {
        self.inner
            .admit(old, new, receipt)
            .map_err(ObservationError::Owner)
    }
    fn encode_state(&self, state: &Self::State, words: &mut Vec<u64>) -> bool {
        self.observation_identity.is_some() && self.inner.encode_state(state, words)
    }
    fn decode_state(&self, words: &[u64]) -> Option<Self::State> {
        self.observation_identity
            .and_then(|_| self.inner.decode_state(words))
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        // Length-prefix the original owner's grammar before appending this
        // adapter's version and optional declared source digest.
        let start = words.len();
        words.push(0);
        self.inner.identity_words(words);
        words[start] = (words.len() - start - 1) as u64;
        words.push(0x7374_6174_6573_6967); // "statesig"
        words.push(1);
        words.push(u64::from(self.observation_identity.is_some()));
        if let Some(identity) = self.observation_identity {
            for bytes in identity.chunks_exact(8) {
                let mut word = [0; 8];
                word.copy_from_slice(bytes);
                words.push(u64::from_le_bytes(word));
            }
        }
    }
}

impl<C, F> PortTangent for StateSignals<C, F>
where
    C: PortTangent,
    F: Fn(&C, &C::State, &mut [f64], &Cx<'_>) -> Result<(), C::Error>,
{
    fn jvp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("d_inputs", self.inputs, d_inputs.len())?;
        self.length("d_outputs", self.outputs, d_outputs.len())?;
        let (native, observed) = d_outputs.split_at_mut(self.native_outputs);
        self.inner
            .jvp(old, dt_s, inputs, d_inputs, native, cx)
            .map_err(ObservationError::Owner)?;
        self.fill(observed, 0.0, cx)
    }
    fn vjp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("bar_inputs", self.inputs, bar_inputs.len())?;
        self.length("bar_outputs", self.outputs, bar_outputs.len())?;
        self.inner
            .vjp(
                old,
                dt_s,
                inputs,
                &bar_outputs[..self.native_outputs],
                bar_inputs,
                cx,
            )
            .map_err(ObservationError::Owner)
    }
}

impl<C, F> LoadLine for StateSignals<C, F>
where
    C: LoadLine,
    F: Fn(&C, &C::State, &mut [f64], &Cx<'_>) -> Result<(), C::Error>,
{
    fn line_ports(&self) -> &[usize] {
        self.inner.line_ports()
    }
    fn line_ends(&self) -> &[Option<usize>] {
        self.inner.line_ends()
    }
    fn reads_mean_flows_only(&self) -> bool {
        self.inner.reads_mean_flows_only()
    }
    fn propose_on_line(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("outputs", self.outputs, outputs.len())?;
        let (native, observed) = outputs.split_at_mut(self.native_outputs);
        let proposal = self
            .inner
            .propose_on_line(old, dt_s, inputs, compliance, new, native, cx)
            .map_err(ObservationError::Owner)?;
        self.observations(proposal, old, observed, cx)
    }
}

impl<C, F> StageResidual for StateSignals<C, F>
where
    C: StageResidual,
    F: Fn(&C, &C::State, &mut [f64], &Cx<'_>) -> Result<(), C::Error>,
{
    fn unknowns(&self) -> usize {
        self.inner.unknowns()
    }
    fn initial_guess(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.inner
            .initial_guess(old, dt_s, inputs, z)
            .map_err(ObservationError::Owner)
    }
    fn residual(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("outputs", self.outputs, outputs.len())?;
        let (native, observed) = outputs.split_at_mut(self.native_outputs);
        self.inner
            .residual(old, dt_s, inputs, z, r, native, cx)
            .map_err(ObservationError::Owner)?;
        self.sample(old, observed, cx)
    }
    fn jvp_z(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("d_outputs", self.outputs, d_outputs.len())?;
        let (native, observed) = d_outputs.split_at_mut(self.native_outputs);
        self.inner
            .jvp_z(old, dt_s, inputs, z, dz, dr, native, cx)
            .map_err(ObservationError::Owner)?;
        self.fill(observed, 0.0, cx)
    }
    fn vjp_z(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.length("bar_outputs", self.outputs, bar_outputs.len())?;
        self.inner
            .vjp_z(
                old,
                dt_s,
                inputs,
                z,
                bar_r,
                &bar_outputs[..self.native_outputs],
                bar_z,
                cx,
            )
            .map_err(ObservationError::Owner)
    }
    fn accept(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut Self::State,
        cx: &Cx<'_>,
    ) -> Result<Receipt, Self::Error> {
        self.length("inputs", self.inputs, inputs.len())?;
        self.inner
            .accept(old, dt_s, inputs, z, new, cx)
            .map_err(ObservationError::Owner)
    }
}

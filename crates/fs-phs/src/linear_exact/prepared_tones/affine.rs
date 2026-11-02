//! Affine envelopes on the same complete physical source and power owners.
use super::{
    DenseExponentialBudget, LinearSplitBudget, PhsError, PortHamiltonian, PreparedExactTones,
    PreparedStepError, PreparedStepRecord, ResistanceWitness, Sum, System, finite, refuse,
    uncancelled,
};

impl PreparedExactTones<'static> {
    /// Own a complete source with affine cosine/sine envelopes. End coefficients
    /// use the interval-start carrier phase; no physical state is reset.
    pub fn new_owned_affine(
        system: PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
    ) -> Result<Self, PhsError> {
        Self::new_owned_affine_with_resistance_channels(
            system,
            frequencies,
            duration,
            split_budget,
            dense_budget,
            &[],
        )
    }

    /// Affine envelopes with independently integrated original resistance
    /// channels. The four-coordinate source extent is admitted before allocation.
    pub fn new_owned_affine_with_resistance_channels(
        system: PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
    ) -> Result<Self, PhsError> {
        Self::new_owned_affine_with_resistance_channels_controlled(
            system,
            frequencies,
            duration,
            split_budget,
            dense_budget,
            channels,
            || false,
        )
        .map_err(uncancelled)
    }

    /// Same original coupled preparation with bounded cooperative cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_owned_affine_with_resistance_channels_controlled(
        system: PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            System::Owned(Box::new(system)),
            frequencies,
            duration,
            split_budget,
            dense_budget,
            channels,
            &mut cancelled,
            true,
        )
    }
}

impl<'a> PreparedExactTones<'a> {
    /// Borrow the unchanged system for a complete affine-envelope source.
    pub fn new_affine(
        system: &'a PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
    ) -> Result<Self, PhsError> {
        Self::prepare(
            System::Borrowed(system),
            frequencies,
            duration,
            split_budget,
            dense_budget,
            &[],
            &mut || false,
            true,
        )
        .map_err(uncancelled)
    }

    /// Whether this owner admits affine endpoints (four quadratures per tone).
    /// Constant preparation retains its original two-coordinate fast path.
    pub fn is_affine(&self) -> bool {
        self.affine
    }

    /// Exact motion/work/heat for
    /// `u(t)=sum[(c0+t/h*(c1-c0))*cos(w*t)+(s0+t/h*(s1-s0))*sin(w*t)]`.
    /// All four port-major arrays share the interval-start carrier phase.
    /// Rotate the END coefficients with `advance_force_phase` for continuation.
    /// This requires affine preparation; a constant plan refuses. All caller
    /// outputs remain unchanged on failure and warm evaluation allocates nothing.
    #[allow(clippy::too_many_arguments)]
    pub fn step_affine_into(
        &mut self,
        initial: &[f64],
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.check_outputs(next, average_output)?;
        let record = self.evaluate(
            initial,
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
        )?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok(record)
    }

    /// Same affine step plus every prepared original resistance-channel heat.
    /// The original signed work, positive losses and physical balance gate
    /// include all envelope/carrier/interference terms on the same trajectory.
    #[allow(clippy::too_many_arguments)]
    pub fn step_affine_observed_into(
        &mut self,
        initial: &[f64],
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
        channel_heat: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.check_outputs(next, average_output)?;
        if channel_heat.len() != self.resistance_channel_count() {
            return Err(refuse("continuous-tone channel output extent"));
        }
        let record = self.evaluate(
            initial,
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
        )?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        channel_heat.copy_from_slice(&self.scratch.heat[1..]);
        Ok(record)
    }
}

// Coalesce complete endpoint values BEFORE subtraction. Forming per-source
// differences first loses small envelope changes among large cancelling tones.
#[allow(clippy::too_many_arguments)]
pub(super) fn combine(
    mapping: &[usize],
    unique: usize,
    c0: &[f64],
    s0: &[f64],
    c1: &[f64],
    s1: &[f64],
    sums: &mut [Sum],
    input: &mut [f64],
) -> Result<(), PhsError> {
    sums.fill(Sum::default());
    let k = mapping.len();
    for p in 0..c0.len() / k {
        for (source, &target) in mapping.iter().enumerate() {
            let at = p * k + source;
            let index = 4 * (p * unique + target);
            for (offset, value) in [c0[at], s0[at], c1[at], s1[at]].into_iter().enumerate() {
                sums[index + offset].add(value);
            }
        }
    }
    for (out, &sum) in input.iter_mut().zip(sums.iter()) {
        *out = sum.value();
    }
    for source in input.chunks_exact_mut(4) {
        source[2] -= source[0];
        source[3] -= source[1];
    }
    if !finite(input) {
        return Err(refuse("combined affine-tone endpoint/delta overflows"));
    }
    Ok(())
}

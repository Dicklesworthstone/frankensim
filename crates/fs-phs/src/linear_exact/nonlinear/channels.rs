//! Source-resolved positive heat on the accepted held nonlinear trajectory.
//!
//! These are read-only physical resistance observations. The held nonlinear
//! effort, complete generator and accepted endpoint are shared with the parent
//! owner. No extra integration or accepted physical state is introduced.
use super::*;

impl PreparedNonlinearZoh {
    /// Retain declared resistance channels on the SAME accepted trajectory.
    ///
    /// Each original witness integrates its positive square independently.
    /// Channels need not partition the system resistance; the caller owns
    /// their physical interpretation and must prevent duplicate thermalization.
    /// The system's original total heat still controls mechanical acceptance.
    /// Only the quadratic limit is a continuous exact physical trajectory;
    /// nonlinear effort and channel heat converge with time refinement.
    ///
    /// Preparation caches one existing positive integral per channel. Byte and
    /// multiply allowances cover the complete retained caches and conservative
    /// cold work. Warm stepping allocates no buffers in this owner.
    pub fn new_owned_with_resistance_channels(
        system: PortHamiltonian,
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
    ) -> Result<Self, PhsError> {
        Self::new_owned_with_resistance_channels_controlled(
            system,
            duration,
            split_budget,
            dense_budget,
            channels,
            || false,
        )
        .map_err(uncancelled)
    }

    /// Cooperatively cancelled preparation of complete source heat factors.
    /// Each bounded positive-integral preparation is one atomic cold operation.
    pub fn new_owned_with_resistance_channels_controlled(
        system: PortHamiltonian,
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        admit(&system, split_budget, dense_budget, channels)?;
        Self::prepare_owned(
            system,
            duration,
            split_budget,
            None,
            channels,
            dense_budget,
            &mut cancelled,
        )
    }

    /// Current numerical Newton-update allowance, independent of physical data.
    pub fn iteration_limit(&self) -> usize {
        self.max_iterations
    }

    /// Number of optional resistance observations in caller order.
    pub fn resistance_channel_count(&self) -> usize {
        self.model.channels.len()
    }

    /// Publish the original state, mean output and source heat on one success.
    pub fn step_observed_into(
        &mut self,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        mean_output: &mut [f64],
        channel_heat: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_observed_into_controlled(initial, input, next, mean_output, channel_heat, || {
            false
        })
        .map_err(uncancelled)
    }

    /// Every output remains untouched on a failed solve, channel, or cancellation.
    /// Heat is the integral of (Q*x+held_nonlinear_effort)^T R_channel
    /// (Q*x+held_nonlinear_effort), with the complete original dynamics.
    #[allow(clippy::too_many_arguments)]
    pub fn step_observed_into_controlled(
        &mut self,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        mean_output: &mut [f64],
        channel_heat: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        admit_outputs(&self.model, next, mean_output)?;
        if channel_heat.len() != self.resistance_channel_count() {
            return Err(refuse("held nonlinear heat needs every prepared channel").into());
        }
        let record = solve(
            &self.model,
            &mut self.scratch,
            self.max_iterations,
            initial,
            input,
            &mut None,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&self.scratch.x);
        mean_output.copy_from_slice(&self.scratch.output);
        channel_heat.copy_from_slice(&self.scratch.heat[1..]);
        Ok(record)
    }
}

fn admit(
    system: &PortHamiltonian,
    split: LinearSplitBudget,
    budget: DenseExponentialBudget,
    channels: &[ResistanceWitness],
) -> Result<(), PhsError> {
    let (n, m) = (system.n, system.m);
    if !PreparedNonlinearZoh::fits_reference_capacity(n, m)
        || n > split.max_modes
        || channels.len() > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
    {
        return Err(refuse(
            "held nonlinear state/port/channel dimensions exceed capacity",
        ));
    }
    let (f, d, c) = (
        (2 * n + m + 1) as u128,
        (2 * n) as u128,
        channels.len() as u128,
    );
    if f.max(d) > budget.max_dimension as u128 {
        return Err(refuse("held nonlinear complete dimensions exceed budget"));
    }
    // Retained flow/output maps, Newton scratch, total and per-channel loss
    // factors, and the largest cold preparation workspace. Use u128 so even an
    // extreme caller allowance cannot wrap the admission arithmetic.
    let bytes = 8
        * (f * f * (2 * m as u128 + 25) + (80 + c) * d * d + (8 + 3 * c) * (n * n) as u128 + c)
        + 128 * (n as u128 + m as u128 + c);
    let work = 224 * f * f * f * (m as u128).max(1) + 2048 * d * d * d * (c + 1);
    if bytes > budget.max_dense_bytes as u128 || work > budget.max_multiply_terms as u128 {
        return Err(refuse(
            "held nonlinear complete caches exceed explicit byte/work budget",
        ));
    }
    if channels.iter().any(|w| {
        w.matrix.len() != n * n || w.root.len() != n * n || !finite(&w.matrix) || !finite(&w.root)
    }) {
        return Err(refuse(
            "held nonlinear channel witness must match the original state",
        ));
    }
    Ok(())
}

pub(super) fn prepare(
    system: &PortHamiltonian,
    split: &super::super::LinearSplit,
    duration: f64,
    channels: &[ResistanceWitness],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<PositiveHeat>, PreparedStepError> {
    let mut prepared = Vec::with_capacity(channels.len());
    for channel in channels {
        poll(cancelled)?;
        prepared.push(PositiveHeat::from_factor(
            system,
            split,
            duration,
            &channel.matrix,
            channel.root.clone(),
        )?);
    }
    poll(cancelled)?;
    Ok(prepared)
}

//! Continuous sources through the original complete quadratic pHS owner.
#[cfg(test)]
mod accounting_tests;
mod affine;
mod observer;
pub(super) mod power;
pub use observer::PreparedToneObserver;

use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, balance_admitted,
    finite,
    prepared::PreparedQuadraticFlow,
    prepared_zoh::{System, energetic_scale_into},
    refuse, split_quadratic,
};
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, ResistanceWitness,
    storage_energy_change,
};
use fs_math::det;
use crate::interconnection::accumulation::Sum as PortIntegral;
use power::{Sum, TonePower};
use std::sync::Arc;

struct Scratch {
    input: Vec<f64>,
    sums: Vec<Sum>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
    port_integrals: Vec<PortIntegral>,
    output: Vec<f64>,
    effort: Vec<f64>,
    midpoint: Vec<f64>,
    derivative: Vec<f64>,
    coordinates: Vec<f64>,
    heat: Vec<f64>,
    work_input: Vec<f64>,
    port_work: Vec<f64>,
}

/// Complete original quadratic dynamics under simultaneous continuous tones.
///
/// Each physical port has caller-ordered cosine/sine coefficients:
/// `u[p](t)=sum_j(Fc[p*K+j]*cos(Omega[j]*t)+Fs[p*K+j]*sin(Omega[j]*t))`.
/// The convention is `exp(-i Omega t)`; DC is explicit. All original states,
/// coupling, damping and material memories remain. Source oscillators are
/// numerical forcing coordinates and acquire no fictitious stored energy.
/// The `new_affine` constructor family additionally interpolates each cosine
/// and sine envelope between supplied endpoints. Those plans require the
/// corresponding `step_affine_into` or `step_affine_observed_into` methods.
///
/// Motion and average outputs use the existing complete exponential/Gramian
/// owner. Original resistance squares and paired source/output factors give
/// independent heat and signed input work, including all cross-frequency
/// terms. The original storage and physical balance gate remain authoritative.
///
/// Preparation owns all numerical buffers. Warm calls allocate nothing and
/// evaluate no exponentials or transcendentals, assuming allocation-free
/// storage callbacks. Caller state and source phase are never mutated.
pub struct PreparedExactTones<'a> {
    identity: Arc<()>,
    affine: bool,
    system: System<'a>,
    duration: f64,
    frequencies: Vec<f64>,
    unique: Vec<f64>,
    mapping: Vec<usize>,
    phase: Vec<(f64, f64)>,
    flow: PreparedQuadraticFlow,
    power: TonePower,
    scratch: Scratch,
}

fn uncancelled(error: PreparedStepError) -> PhsError {
    match error {
        PreparedStepError::Solver(error) => error,
        PreparedStepError::Cancelled => {
            refuse("uncancelled continuous-tone preparation interrupted")
        }
    }
}

impl PreparedExactTones<'static> {
    /// Own the unchanged physical system and its complete fixed-clock caches.
    pub fn new_owned(
        system: PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
    ) -> Result<Self, PhsError> {
        Self::new_owned_with_resistance_channels(
            system,
            frequencies,
            duration,
            split_budget,
            dense_budget,
            &[],
        )
    }

    /// Additionally prepare original positive resistance channels. Each witness
    /// retains its physical matrix and source squares. Channels may overlap;
    /// the caller owns their interpretation and must not assume their sum is R.
    /// Total heat always uses the system's original complete R independently.
    pub fn new_owned_with_resistance_channels(
        system: PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
    ) -> Result<Self, PhsError> {
        Self::prepare(
            System::Owned(Box::new(system)),
            frequencies,
            duration,
            split_budget,
            dense_budget,
            channels,
            &mut || false,
            false,
        )
        .map_err(uncancelled)
    }

    /// Same cold preparation with polling between bounded shared numerical
    /// kernels and channel factors. No partly prepared owner is published.
    #[allow(clippy::too_many_arguments)]
    pub fn new_owned_with_resistance_channels_controlled(
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
            false,
        )
    }
}

impl<'a> PreparedExactTones<'a> {
    /// Borrow the immutable original quadratic system for one positive interval.
    /// Capacity counts complete source quadratures and every channel before
    /// allocation. No state/tone truncation or automatic clock change occurs.
    pub fn new(
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
            false,
        )
        .map_err(uncancelled)
    }

    /// Borrow the original system and retain its separately declared positive
    /// resistance channels. Preparation checks cancellation between bounded
    /// kernels; no reconstructed system or copied storage owner is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_resistance_channels_controlled(
        system: &'a PortHamiltonian,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        Self::prepare(
            System::Borrowed(system),
            frequencies,
            duration,
            split_budget,
            dense_budget,
            channels,
            &mut cancelled,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        owner: System<'a>,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
        cancelled: &mut impl FnMut() -> bool,
        affine: bool,
    ) -> Result<Self, PreparedStepError> {
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        let system = owner.get();
        let (n, m) = (system.n, system.m);
        admit(
            system,
            frequencies,
            duration,
            channels,
            split_budget,
            dense_budget,
            affine,
        )?;
        let mut unique = Vec::new();
        let mut mapping = Vec::with_capacity(frequencies.len());
        for &omega in frequencies {
            let index = if let Some(index) = unique.iter().position(|&old| old == omega) {
                index
            } else {
                unique.push(omega);
                unique.len() - 1
            };
            mapping.push(index);
        }
        let split = split_quadratic(system, split_budget)?;
        if !split.is_exact_quadratic() {
            return Err(
                refuse("continuous-tone evolution requires exactly quadratic storage").into(),
            );
        }
        let q = split.quadratic_matrix();
        let mut outputs = vec![vec![0.0; n]; m];
        for (p, row) in outputs.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = crate::port_output_component(&system.g, m, p, |i| q[i * n + j])
                    .ok_or_else(|| refuse("continuous-tone physical output map overflows"))?;
            }
        }
        let flow = PreparedQuadraticFlow::new_tones_controlled(
            system,
            &split,
            duration,
            &unique,
            affine,
            &outputs,
            dense_budget,
            &mut *cancelled,
        )?;
        let power = TonePower::new(
            system,
            &split,
            &unique,
            flow.source_scales(),
            duration,
            affine,
            channels,
            cancelled,
        )?;
        let sources = if affine { 4 } else { 2 } * m * unique.len();
        let scratch = Scratch {
            input: vec![0.0; sources],
            sums: vec![Sum::default(); sources],
            augmented: vec![0.0; n + sources + 1],
            balanced: vec![0.0; n],
            endpoint: vec![0.0; n],
            integrals: vec![0.0; m],
            port_integrals: vec![PortIntegral::default(); m],
            output: vec![0.0; m],
            effort: vec![0.0; n],
            midpoint: vec![0.0; n],
            derivative: vec![0.0; n],
            coordinates: vec![0.0; power.coordinate_len()],
            heat: vec![0.0; channels.len() + 1],
            work_input: vec![0.0; sources],
            port_work: vec![0.0; m],
        };
        let phase = unique
            .iter()
            .map(|omega| det::sin_cos(omega * duration))
            .collect();
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        Ok(Self {
            identity: Arc::new(()),
            affine,
            system: owner,
            duration,
            frequencies: frequencies.to_vec(),
            unique,
            mapping,
            phase,
            flow,
            power,
            scratch,
        })
    }

    /// Read-only original physical system, in its unchanged state/port order.
    pub fn system(&self) -> &PortHamiltonian {
        self.system.get()
    }
    /// Prepared interval duration [s].
    pub fn duration(&self) -> f64 {
        self.duration
    }
    /// Internal transfer view, consumed only after a successful native step.
    pub(crate) fn port_integrals(&self) -> &[PortIntegral] {
        &self.scratch.port_integrals
    }
    /// Number of original physical state coordinates.
    pub fn state_dim(&self) -> usize {
        self.system().n
    }
    /// Number of original physical source ports.
    pub fn port_dim(&self) -> usize {
        self.system().m
    }
    /// Number of coefficients required per port, in original caller order.
    pub fn tone_count(&self) -> usize {
        self.frequencies.len()
    }
    /// Original caller-ordered angular frequencies [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }
    /// Number of additional original resistance observations.
    pub fn resistance_channel_count(&self) -> usize {
        self.scratch.heat.len() - 1
    }

    /// Rotate one caller tone after an accepted physical step. Checkpoints must
    /// retain these external coefficients together with the physical state.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), PhsError> {
        let &index = self
            .mapping
            .get(tone)
            .ok_or_else(|| refuse("continuous-tone phase index"))?;
        let (sn, cs) = self.phase[index];
        let next = (cs * cosine + sn * sine, cs * sine - sn * cosine);
        if !finite(&[cosine, sine, next.0, next.1]) {
            return Err(refuse("continuous-tone phase rotation is unrepresentable"));
        }
        Ok(next)
    }

    /// Exact motion, actual average outputs and independent original work/heat.
    /// Output buffers remain unchanged on every refusal; scratch may change.
    pub fn step_into(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.check_outputs(next, average_output)?;
        let record = self.evaluate(initial, cosine, sine, None)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok(record)
    }

    /// Same step plus every prepared original resistance-channel integral.
    /// All caller outputs publish atomically after the unchanged physical gate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_observed_into(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
        channel_heat: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.check_outputs(next, average_output)?;
        if channel_heat.len() != self.resistance_channel_count() {
            return Err(refuse("continuous-tone channel output extent"));
        }
        let record = self.evaluate(initial, cosine, sine, None)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        channel_heat.copy_from_slice(&self.scratch.heat[1..]);
        Ok(record)
    }

    /// Advance continuous tones and independently integrate each physical
    /// port's signed work on the SAME trajectory, plus all resistance channels.
    /// Work includes startup and cross-frequency terms; it is not pressure
    /// average times flow average, nor a remainder of the energy balance.
    /// Caller buffers publish together only after every account is accepted.
    #[allow(clippy::too_many_arguments)]
    pub fn step_accounted_into(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
        channel_heat: &mut [f64],
        port_work: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_accounted_into_controlled(
            initial,
            cosine,
            sine,
            next,
            average_output,
            channel_heat,
            port_work,
            || false,
        )
        .map_err(uncancelled)
    }

    /// Same allocation-free step with checks before, between physical port
    /// accounts, and immediately before publication. Private scratch is reusable
    /// after refusal; all caller output buffers remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_accounted_into_controlled(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        next: &mut [f64],
        average_output: &mut [f64],
        channel_heat: &mut [f64],
        port_work: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        self.check_outputs(next, average_output)?;
        if channel_heat.len() != self.resistance_channel_count()
            || port_work.len() != self.port_dim()
        {
            return Err(refuse("continuous-tone heat/work output extent").into());
        }
        let record = self.evaluate(initial, cosine, sine, None)?;
        // evaluate() uses gradient scratch while checking storage. Restore the
        // actual initial effort required by the original source-work integral.
        self.system
            .get()
            .storage
            .gradient(initial, &mut self.scratch.effort);
        self.scratch.work_input.fill(0.0);
        let width = 2 * self.unique.len();
        let mut supplied = Sum::default();
        let mut magnitude = Sum::default();
        for port in 0..self.port_dim() {
            if cancelled() {
                return Err(PreparedStepError::Cancelled);
            }
            let range = port * width..(port + 1) * width;
            self.scratch.work_input[range.clone()]
                .copy_from_slice(&self.scratch.input[range.clone()]);
            let work = self.power.source_work_into(
                &self.scratch.effort,
                &self.scratch.derivative,
                &self.scratch.input,
                &self.scratch.work_input,
                &mut self.scratch.coordinates,
            )?;
            self.scratch.work_input[range].fill(0.0);
            self.scratch.port_work[port] = work;
            supplied.add(work);
            magnitude.add(work.abs());
        }
        let total = supplied.value();
        let scale = magnitude.value() + record.supplied.abs();
        if !finite(&self.scratch.port_work)
            || !finite(&[total, scale])
            || (total - record.supplied).abs() > 1e-12 * scale
        {
            return Err(refuse("continuous-tone original port work does not close supply").into());
        }
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        channel_heat.copy_from_slice(&self.scratch.heat[1..]);
        port_work.copy_from_slice(&self.scratch.port_work);
        Ok(record)
    }

    fn check_outputs(&self, next: &[f64], output: &[f64]) -> Result<(), PhsError> {
        if next.len() != self.state_dim() || output.len() != self.port_dim() {
            return Err(refuse(
                "continuous-tone output buffers must retain complete state and ports",
            ));
        }
        Ok(())
    }

    /// Internal composition seam: retain the original total-input trajectory and
    /// separately integrate one declared source's work on that trajectory. The
    /// caller supplies reusable packing scratch; no additional owner or heat law
    /// is introduced. Neither source work nor interface work is an energy remainder.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_source_work_into(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        end: Option<(&[f64], &[f64])>,
        source_cosine: &[f64],
        source_sine: &[f64],
        source_end: Option<(&[f64], &[f64])>,
        held_exchange: &[f64],
        work_input: &mut [f64],
        next: &mut [f64],
        average_output: &mut [f64],
    ) -> Result<(PreparedStepRecord, f64, f64), PhsError> {
        self.check_outputs(next, average_output)?;
        let count = self.port_dim() * self.tone_count();
        if source_cosine.len() != count
            || source_sine.len() != count
            || !finite(source_cosine)
            || !finite(source_sine)
            || self.affine != source_end.is_some()
            || source_end.is_some_and(|(c, s)| {
                c.len() != count || s.len() != count || !finite(c) || !finite(s)
            })
            || work_input.len() < self.scratch.input.len()
            || held_exchange.len() != self.port_dim()
            || !finite(held_exchange)
            || self.frequencies.first().copied() != Some(0.0)
        {
            return Err(refuse(
                "continuous-tone source-work input or scratch extent",
            ));
        }
        let record = self.evaluate(initial, cosine, sine, end)?;
        let work_input = &mut work_input[..self.scratch.input.len()];
        if let Some((end_cosine, end_sine)) = source_end {
            affine::combine(
                &self.mapping,
                self.unique.len(),
                source_cosine,
                source_sine,
                end_cosine,
                end_sine,
                &mut self.scratch.sums,
                work_input,
            )?;
        } else {
            combine_tones(
                &self.mapping,
                self.unique.len(),
                source_cosine,
                source_sine,
                &mut self.scratch.sums,
                work_input,
            )?;
        }
        let system = self.system.get();
        system.storage.gradient(initial, &mut self.scratch.effort);
        let source_work = self.power.source_work_into(
            &self.scratch.effort,
            &self.scratch.derivative,
            &self.scratch.input,
            work_input,
            &mut self.scratch.coordinates,
        )?;
        // The multirate adapter reserves its first tone for the held interface.
        // Use the same native power integral as the physical/source work; the
        // rounded mean outputs can have lost small opposing-port transfers.
        work_input.fill(0.0);
        let stride = if self.affine { 4 } else { 2 };
        for (port, &input) in held_exchange.iter().enumerate() {
            work_input[stride * (port * self.unique.len() + self.mapping[0])] = input;
        }
        let exchange_work = self.power.source_work_into(
            &self.scratch.effort,
            &self.scratch.derivative,
            &self.scratch.input,
            work_input,
            &mut self.scratch.coordinates,
        )?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        Ok((record, source_work, exchange_work))
    }

    fn evaluate(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
        end: Option<(&[f64], &[f64])>,
    ) -> Result<PreparedStepRecord, PhsError> {
        let system = self.system.get();
        let (n, m) = (system.n, system.m);
        let k = self.frequencies.len();
        if self.affine != end.is_some()
            || initial.len() != n
            || cosine.len() != m * k
            || sine.len() != m * k
            || !finite(initial)
            || !finite(cosine)
            || !finite(sine)
            || end.is_some_and(|(c, s)| {
                c.len() != m * k || s.len() != m * k || !finite(c) || !finite(s)
            })
        {
            return Err(refuse(
                "continuous-tone complete finite state and port-major inputs required",
            ));
        }
        let scratch = &mut self.scratch;
        if let Some((end_cosine, end_sine)) = end {
            affine::combine(
                &self.mapping,
                self.unique.len(),
                cosine,
                sine,
                end_cosine,
                end_sine,
                &mut scratch.sums,
                &mut scratch.input,
            )?;
        } else {
            combine_tones(
                &self.mapping,
                self.unique.len(),
                cosine,
                sine,
                &mut scratch.sums,
                &mut scratch.input,
            )?;
        }
        self.flow.evaluate_port_integrals_into(
            initial,
            &scratch.input,
            &mut scratch.augmented,
            &mut scratch.balanced,
            &mut scratch.endpoint,
            &mut scratch.integrals,
            &mut scratch.port_integrals,
        )?;
        for (out, &integral) in scratch.output.iter_mut().zip(&scratch.integrals) {
            *out = integral / self.duration;
        }
        system.storage.gradient(initial, &mut scratch.effort);
        let supplied = self.power.evaluate_into(
            system,
            &scratch.effort,
            &scratch.input,
            &mut scratch.derivative,
            &mut scratch.coordinates,
            &mut scratch.heat,
        )?;
        let dissipated = scratch.heat[0];
        let delta_h = storage_energy_change(
            system.storage.as_ref(),
            initial,
            &scratch.endpoint,
            &mut scratch.midpoint,
            &mut scratch.effort,
        );
        let scale = energetic_scale_into(system, initial, &mut scratch.effort)
            .max(energetic_scale_into(
                system,
                &scratch.endpoint,
                &mut scratch.effort,
            ))
            .max(delta_h.abs())
            .max(supplied.abs())
            .max(dissipated.abs());
        let residual = delta_h + dissipated - supplied;
        if !finite(&scratch.output)
            || !finite(&scratch.heat)
            || !system.hamiltonian(&scratch.endpoint).is_finite()
            || !finite(&[delta_h, dissipated, supplied, scale, residual])
            || scratch.heat.iter().any(|v| *v < 0.0)
            || !balance_admitted(delta_h, dissipated, supplied, scale, self.power.doublings())
        {
            return Err(refuse(
                "continuous-tone original physical energy balance refused",
            ));
        }
        Ok(PreparedStepRecord {
            delta_h,
            dissipated,
            supplied,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        })
    }
}

fn admit(
    system: &PortHamiltonian,
    frequencies: &[f64],
    dt: f64,
    channels: &[ResistanceWitness],
    split: LinearSplitBudget,
    budget: DenseExponentialBudget,
    affine: bool,
) -> Result<(), PhsError> {
    let (n, m) = (system.n, system.m);
    if n == 0
        || m == 0
        || frequencies.is_empty()
        || !dt.is_finite()
        || dt <= 0.0
        || (affine && !(1.0 / dt).is_finite())
        || frequencies
            .iter()
            .any(|w| !w.is_finite() || *w < 0.0 || !(w * dt).is_finite())
    {
        return Err(refuse(
            "continuous tones require finite frequencies, physical ports and positive clock",
        ));
    }
    let source = m
        .checked_mul(frequencies.len())
        .and_then(|v| v.checked_mul(if affine { 4 } else { 2 }))
        .ok_or_else(|| refuse("continuous-tone source dimension overflow"))?;
    let flow = n
        .checked_add(source)
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| refuse("continuous-tone flow dimension overflow"))?;
    let power = n
        .checked_add(n.max(m))
        .and_then(|v| v.checked_add(source))
        .ok_or_else(|| refuse("continuous-tone power dimension overflow"))?;
    let dimension = flow.max(power);
    if n > split.max_modes
        || dimension > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || dimension > budget.max_dimension
    {
        return Err(refuse(
            "complete continuous-tone state/source/power dimension exceeds budget",
        ));
    }
    // Conservative numerical payload/work, including every optional channel,
    // BEFORE any factor/candidate allocation. Existing bounded kernels retain
    // their own 32-term/48-doubling admission. Scores are not wall-clock claims.
    let (f, d, c) = (flow as u128, power as u128, channels.len() as u128);
    let bytes = 8 * (f * f * (2 * m as u128 + 25) + d * d * (c + 80) + (n * n) as u128 * c)
        + 128 * (n + m + source + frequencies.len()) as u128
        + 8 * (source + m) as u128;
    let work =
        224 * f * f * f * (m as u128).max(1) + 2048 * d * d * d * (c + 2) + 8 * m as u128 * d * d;
    if bytes > budget.max_dense_bytes as u128 || work > budget.max_multiply_terms as u128 {
        return Err(refuse(
            "complete continuous-tone channels exceed explicit byte/work budget",
        ));
    }
    if channels.iter().any(|w| {
        w.matrix.len() != n * n || w.root.len() != n * n || !finite(&w.matrix) || !finite(&w.root)
    }) {
        return Err(refuse(
            "continuous-tone channel witness must match original state layout",
        ));
    }
    Ok(())
}

// Shared original-order compensated source packing for motion and observers.
fn combine_tones(
    mapping: &[usize],
    unique: usize,
    cosine: &[f64],
    sine: &[f64],
    sums: &mut [Sum],
    input: &mut [f64],
) -> Result<(), PhsError> {
    sums.fill(Sum::default());
    let k = mapping.len();
    for p in 0..cosine.len() / k {
        for (source, &target) in mapping.iter().enumerate() {
            let index = 2 * (p * unique + target);
            sums[index].add(cosine[p * k + source]);
            sums[index + 1].add(sine[p * k + source]);
        }
    }
    for (out, &sum) in input.iter_mut().zip(sums.iter()) {
        *out = sum.value();
    }
    if !finite(input) {
        return Err(refuse("combined continuous-tone source overflows"));
    }
    Ok(())
}

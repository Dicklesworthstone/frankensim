//! Continuous sources through the original complete quadratic pHS owner.
mod power;

use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    prepared::PreparedQuadraticFlow,
    prepared_zoh::{System, energetic_scale_into},
    refuse, split_quadratic,
};
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, ResistanceWitness,
    storage_energy_change,
};
use fs_math::det;
use power::{Sum, TonePower};

struct Scratch {
    input: Vec<f64>,
    sums: Vec<Sum>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
    output: Vec<f64>,
    effort: Vec<f64>,
    midpoint: Vec<f64>,
    derivative: Vec<f64>,
    coordinates: Vec<f64>,
    heat: Vec<f64>,
}

/// Complete original quadratic dynamics under simultaneous continuous tones.
///
/// Each physical port has caller-ordered cosine/sine coefficients:
/// `u[p](t)=sum_j(Fc[p*K+j]*cos(Omega[j]*t)+Fs[p*K+j]*sin(Omega[j]*t))`.
/// The convention is `exp(-i Omega t)`; DC is explicit. All original states,
/// coupling, damping and material memories remain. Source oscillators are
/// numerical forcing coordinates and acquire no fictitious stored energy.
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
        )
        .map_err(uncancelled)
    }

    fn prepare(
        owner: System<'a>,
        frequencies: &[f64],
        duration: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        channels: &[ResistanceWitness],
        cancelled: &mut impl FnMut() -> bool,
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
        let outputs: Vec<Vec<f64>> = (0..m)
            .map(|p| {
                (0..n)
                    .map(|j| (0..n).map(|i| system.g[i * m + p] * q[i * n + j]).sum())
                    .collect()
            })
            .collect();
        let flow = PreparedQuadraticFlow::new_tones_controlled(
            system,
            &split,
            duration,
            &unique,
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
            channels,
            cancelled,
        )?;
        let sources = 2 * m * unique.len();
        let scratch = Scratch {
            input: vec![0.0; sources],
            sums: vec![Sum::default(); sources],
            augmented: vec![0.0; n + sources + 1],
            balanced: vec![0.0; n],
            endpoint: vec![0.0; n],
            integrals: vec![0.0; m],
            output: vec![0.0; m],
            effort: vec![0.0; n],
            midpoint: vec![0.0; n],
            derivative: vec![0.0; n],
            coordinates: vec![0.0; n + n.max(m) + sources],
            heat: vec![0.0; channels.len() + 1],
        };
        let phase = unique
            .iter()
            .map(|omega| det::sin_cos(omega * duration))
            .collect();
        if cancelled() {
            return Err(PreparedStepError::Cancelled);
        }
        Ok(Self {
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
        let record = self.evaluate(initial, cosine, sine)?;
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
        let record = self.evaluate(initial, cosine, sine)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_output.copy_from_slice(&self.scratch.output);
        channel_heat.copy_from_slice(&self.scratch.heat[1..]);
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

    fn evaluate(
        &mut self,
        initial: &[f64],
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        let system = self.system.get();
        let (n, m) = (system.n, system.m);
        let k = self.frequencies.len();
        if initial.len() != n
            || cosine.len() != m * k
            || sine.len() != m * k
            || !finite(initial)
            || !finite(cosine)
            || !finite(sine)
        {
            return Err(refuse(
                "continuous-tone complete finite state and port-major inputs required",
            ));
        }
        let scratch = &mut self.scratch;
        scratch.sums.fill(Sum::default());
        for p in 0..m {
            for (source, &target) in self.mapping.iter().enumerate() {
                let index = 2 * (p * self.unique.len() + target);
                scratch.sums[index].add(cosine[p * k + source]);
                scratch.sums[index + 1].add(sine[p * k + source]);
            }
        }
        for (out, &sum) in scratch.input.iter_mut().zip(&scratch.sums) {
            *out = sum.value();
        }
        self.flow.evaluate_into(
            initial,
            &scratch.input,
            &mut scratch.augmented,
            &mut scratch.balanced,
            &mut scratch.endpoint,
            &mut scratch.integrals,
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
            || residual.abs() > 1e-12 * scale
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
) -> Result<(), PhsError> {
    let (n, m) = (system.n, system.m);
    if n == 0
        || m == 0
        || frequencies.is_empty()
        || !dt.is_finite()
        || dt <= 0.0
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
        .and_then(|v| v.checked_mul(2))
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
    let bytes = 8 * (f * f * (2 * m as u128 + 24) + d * d * (c + 80) + (n * n) as u128 * c)
        + 128 * (n + m + source + frequencies.len()) as u128;
    let work = 224 * f * f * f * (m as u128).max(1) + 2048 * d * d * d * (c + 2);
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

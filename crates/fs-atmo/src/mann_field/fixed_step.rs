//! Mann amplitudes evolved on an explicitly declared fixed physical clock.
//!
//! Each step uses the exact OU transition `a' = exp(-dt/tau) a + sigma xi`,
//! with `sigma² = -expm1(-2dt/tau)` for unit stationary variance. See
//! Gillespie (1996), doi:10.1103/PhysRevE.54.2084. The supplied snapshot is
//! retained exactly at step zero; a distinct stream supplies future draws.
//! A realization is bound to its dt: changing dt does not preserve the same
//! path or establish common-random-number refinement equivalence.
//!
//! Spatial synthesis, analytic derivatives and observation filters retain the
//! original full Mann tensor. Authored OU times remain separate from Mann's
//! spatial-distortion lifetime. OU evolution is not temporally band-limited;
//! this supplies neither a temporal Nyquist guarantee nor resolved eddy turnover.

use super::*;
use fs_rand::{STREAM_SEMANTICS_VERSION, Stream, StreamCheckpoint};
use std::sync::Arc;

/// Separate innovation namespace; no existing 120 Hz stream is changed.
pub const FIXED_STEP_MANN_KERNEL: u32 = 0x4d46_5354;
/// Exact integer-to-f64 step range, also below the Philox counter limit.
pub const MAX_FIXED_STEP: u64 = (1u64 << 53) - 1;

#[derive(Clone, Debug, PartialEq)]
struct Definition {
    model: MannSpectralModel,
    times: Vec<f64>,
    transition: Vec<[f64; 2]>,
    dt_s: f64,
    origin_time_s: f64,
    origin_evolution_seed: u64,
    evolution_seed: u64,
}

/// A prepared temporal owner. Clones share immutable spectrum/clock data;
/// each owns accepted and trial amplitudes. Successful advances and samples
/// allocate no numerical storage. Refusals preserve the accepted state.
#[derive(Clone, Debug)]
pub struct FixedStepMannField {
    definition: Arc<Definition>,
    step: u64,
    amplitudes: [Vec<f64>; 4],
    trial: [Vec<f64>; 4],
}

impl PartialEq for FixedStepMannField {
    fn eq(&self, other: &Self) -> bool {
        self.definition == other.definition
            && self.step == other.step
            && self.amplitudes == other.amplitudes
    }
}

/// Complete accepted state. Scratch is deliberately absent. Version one
/// binds the fixed-step counter namespace and its four quadrature keys.
#[derive(Clone, Debug, PartialEq)]
pub struct FixedStepMannCheckpoint {
    pub version: u32,
    pub stream_semantics_version: u32,
    pub model: MannSpectralModel,
    pub correlation_times_s: Vec<f64>,
    pub dt_s: f64,
    pub origin_time_s: f64,
    pub origin_evolution_seed: u64,
    pub evolution_seed: u64,
    pub step: u64,
    pub amplitudes: [Vec<f64>; 4],
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use a matching finite fixed clock and explicit work/storage budgets",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-fixed-step-cancelled",
            "temporal advance cancelled",
        ))
    } else {
        Ok(())
    }
}

fn time_at(origin: f64, dt: f64, step: u64) -> Result<f64, Refusal> {
    let elapsed = step as f64 * dt;
    let next_elapsed = (step as f64 + 1.0) * dt;
    let time = origin + elapsed;
    let next_time = origin + next_elapsed;
    if step > MAX_FIXED_STEP
        || !time.is_finite()
        || !elapsed.is_finite()
        || !next_time.is_finite()
        || next_elapsed <= elapsed
        || next_time <= time
    {
        return Err(error(
            "mann-fixed-step-clock",
            "fixed-step time is not representable",
        ));
    }
    Ok(time)
}

fn zeros(n: usize) -> Result<Vec<f64>, Refusal> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| {
        error(
            "mann-fixed-step-allocation",
            "unable to reserve temporal workspace",
        )
    })?;
    values.resize(n, 0.0);
    Ok(values)
}

impl FixedStepMannField {
    /// Retain this snapshot's amplitudes and physical advection phase at step
    /// zero. The new seed addresses only future fixed-step innovations.
    pub fn from_snapshot(
        snapshot: HomogeneousMannField,
        evolution_seed: u64,
        dt_s: f64,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        let origin_time_s = snapshot.tick() as f64 / TICK_HZ;
        Self::resume(
            FixedStepMannCheckpoint {
                version: 1,
                stream_semantics_version: STREAM_SEMANTICS_VERSION,
                model: snapshot.model,
                correlation_times_s: snapshot.correlation_times_s,
                dt_s,
                origin_time_s,
                origin_evolution_seed: snapshot.evolution_seed,
                evolution_seed,
                step: 0,
                amplitudes: snapshot.paths.map(StationaryOuPath::into_amplitudes),
            },
            maximum_workspace_bytes,
        )
    }

    /// Nominal temporal numerical storage: 88 bytes per mode. Includes both
    /// amplitude banks, correlation times and transition pairs; excludes the
    /// immutable spectral geometry, fixed metadata and allocator metadata.
    pub fn required_workspace_bytes(modes: usize) -> Result<usize, Refusal> {
        modes
            .checked_mul(11 * core::mem::size_of::<f64>())
            .ok_or_else(|| {
                error(
                    "mann-fixed-step-budget",
                    "temporal workspace size overflows",
                )
            })
    }

    /// Actual capacities, counting shared immutable temporal arrays once for
    /// each owner. Summing this over clones is a conservative combined budget.
    pub fn workspace_bytes(&self) -> usize {
        (self.definition.times.capacity()
            + 2 * self.definition.transition.capacity()
            + self
                .amplitudes
                .iter()
                .chain(&self.trial)
                .map(Vec::capacity)
                .sum::<usize>())
            * core::mem::size_of::<f64>()
    }

    pub fn model(&self) -> &MannSpectralModel {
        &self.definition.model
    }
    pub fn dt_s(&self) -> f64 {
        self.definition.dt_s
    }
    pub const fn step(&self) -> u64 {
        self.step
    }
    pub fn elapsed_time_s(&self) -> f64 {
        self.step as f64 * self.dt_s()
    }
    /// Original snapshot's phase time; elapsed CFD time starts at zero.
    pub fn origin_time_s(&self) -> f64 {
        self.definition.origin_time_s
    }
    pub fn time_s(&self) -> f64 {
        self.origin_time_s() + self.elapsed_time_s()
    }
    pub fn advection_mps(&self) -> [f64; 3] {
        self.model().advection_mps
    }
    pub fn correlation_times_s(&self) -> &[f64] {
        &self.definition.times
    }

    /// Copy accepted amplitudes/step into existing storage only when spectrum,
    /// dt, phase origin, original realization and future stream all match.
    pub fn copy_temporal_state_from(&mut self, source: &Self) -> Result<(), Refusal> {
        if !Arc::ptr_eq(&self.definition, &source.definition)
            && self.definition != source.definition
        {
            return Err(error(
                "mann-fixed-step-source-mismatch",
                "temporal state belongs to a different model or clock",
            ));
        }
        for (target, source) in self.amplitudes.iter_mut().zip(&source.amplitudes) {
            target.copy_from_slice(source);
        }
        self.step = source.step;
        Ok(())
    }

    /// Atomically advance all four paths. The explicit work budget counts
    /// `4 * mode_count * (target_step - step)` normal draws. Cancellation is
    /// polled at most every 64 draws and once immediately before publication.
    pub fn advance_to(
        &mut self,
        target_step: u64,
        maximum_mode_transitions: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        poll(&mut cancelled)?;
        let span = target_step.checked_sub(self.step).ok_or_else(|| {
            error(
                "mann-fixed-step-backwards",
                "fixed-step clock cannot advance backwards",
            )
        })?;
        time_at(self.origin_time_s(), self.dt_s(), target_step)?;
        let work = usize::try_from(span)
            .ok()
            .and_then(|s| s.checked_mul(4 * self.model().mode_count()))
            .ok_or_else(|| error("mann-fixed-step-work", "transition work overflows"))?;
        if work > maximum_mode_transitions {
            return Err(error(
                "mann-fixed-step-work",
                "transition work exceeds the declared budget",
            ));
        }
        if span == 0 {
            return poll(&mut cancelled);
        }
        for (trial, accepted) in self.trial.iter_mut().zip(&self.amplitudes) {
            trial.copy_from_slice(accepted);
        }
        for step in self.step + 1..=target_step {
            for (path, amplitudes) in self.trial.iter_mut().enumerate() {
                for (mode, (a, &[rho, sigma])) in amplitudes
                    .iter_mut()
                    .zip(&self.definition.transition)
                    .enumerate()
                {
                    if mode % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let key = StreamKey {
                        seed: self.definition.evolution_seed ^ PATH_KEYS[path],
                        kernel: FIXED_STEP_MANN_KERNEL,
                        tile: mode as u32,
                    };
                    let mut rng = Stream::resume(StreamCheckpoint::current(key, 2 * step))
                        .expect("current-version checkpoint and admitted counter always resume");
                    *a = rho * *a + sigma * rng.next_normal();
                    if !a.is_finite() {
                        return Err(error(
                            "mann-fixed-step-state",
                            "evolved amplitude is not finite",
                        ));
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.amplitudes, &mut self.trial);
        self.step = target_step;
        Ok(())
    }

    /// Read the current amplitudes and physical advection phase without draws.
    pub fn sample(&self, point_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, None)
    }
    /// Read the accepted quadratures through the same harmonic slip-plane
    /// observation as the original clock owner; no additional state or draws.
    pub fn sample_wall_blocked(
        &self,
        point_m: [f64; 3],
        wall: MannPlaneWall,
    ) -> Result<FlowSample, Refusal> {
        wall::sample_modes(
            self.model(),
            self.time_s(),
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            wall,
        )
    }

    pub fn sample_filtered(
        &self,
        point_m: [f64; 3],
        filter: &MannBoxFilter,
    ) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, filter.gains(self.model())?)
    }
    fn sample_with_gains(
        &self,
        point_m: [f64; 3],
        gains: Option<&[f64]>,
    ) -> Result<FlowSample, Refusal> {
        sample_modes(
            self.model(),
            self.time_s(),
            core::array::from_fn(|i| self.amplitudes[i].as_slice()),
            point_m,
            gains,
        )
    }
    /// Underlying stationary space-time covariance; realized samples occur
    /// only at integer multiples of this owner's declared dt.
    pub fn covariance(&self, separation_m: [f64; 3], lag_s: f64) -> Result<[[f64; 3]; 3], Refusal> {
        covariance_modes(
            self.model(),
            self.correlation_times_s(),
            separation_m,
            lag_s,
            None,
        )
    }
    /// Stationary two-point wall-blocked covariance with this owner's
    /// original correlation times and phase advection. Both heights matter.
    pub fn covariance_wall_blocked(
        &self,
        first_m: [f64; 3],
        second_m: [f64; 3],
        lag_s: f64,
        wall: MannPlaneWall,
    ) -> Result<[[f64; 3]; 3], Refusal> {
        wall::covariance_modes(
            self.model(),
            Some(self.correlation_times_s()),
            first_m,
            second_m,
            lag_s,
            wall,
        )
    }

    pub fn covariance_filtered(
        &self,
        separation_m: [f64; 3],
        lag_s: f64,
        filter: &MannBoxFilter,
    ) -> Result<[[f64; 3]; 3], Refusal> {
        covariance_modes(
            self.model(),
            self.correlation_times_s(),
            separation_m,
            lag_s,
            filter.gains(self.model())?,
        )
    }

    pub fn checkpoint(&self) -> FixedStepMannCheckpoint {
        FixedStepMannCheckpoint {
            version: 1,
            stream_semantics_version: STREAM_SEMANTICS_VERSION,
            model: self.model().clone(),
            correlation_times_s: self.definition.times.clone(),
            dt_s: self.dt_s(),
            origin_time_s: self.origin_time_s(),
            origin_evolution_seed: self.definition.origin_evolution_seed,
            evolution_seed: self.definition.evolution_seed,
            step: self.step,
            amplitudes: self.amplitudes.clone(),
        }
    }

    /// Restore the exact accepted path and admit actual numerical capacities.
    pub fn resume(
        checkpoint: FixedStepMannCheckpoint,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        let n = checkpoint.model.mode_count();
        if checkpoint.version != 1
            || checkpoint.stream_semantics_version != STREAM_SEMANTICS_VERSION
            || !checkpoint.dt_s.is_finite()
            || checkpoint.dt_s <= 0.0
            || !checkpoint.origin_time_s.is_finite()
            || checkpoint.correlation_times_s.len() != n
            || checkpoint
                .correlation_times_s
                .iter()
                .any(|tau| !tau.is_finite() || *tau <= 0.0)
            || checkpoint
                .amplitudes
                .iter()
                .any(|a| a.len() != n || a.iter().any(|a| !a.is_finite()))
        {
            return Err(error(
                "mann-fixed-step-state",
                "invalid clock, model or checkpoint amplitudes",
            ));
        }
        time_at(checkpoint.origin_time_s, checkpoint.dt_s, checkpoint.step)?;
        if Self::required_workspace_bytes(n)? > maximum_workspace_bytes {
            return Err(error(
                "mann-fixed-step-budget",
                "temporal workspace exceeds its byte budget",
            ));
        }
        let mut transition = Vec::new();
        transition.try_reserve_exact(n).map_err(|_| {
            error(
                "mann-fixed-step-allocation",
                "unable to reserve transitions",
            )
        })?;
        for tau in &checkpoint.correlation_times_s {
            let ratio = checkpoint.dt_s / tau;
            let rho = det::exp(-ratio);
            let sigma = det::sqrt(-det::expm1(-2.0 * ratio));
            if ratio == 0.0 || !rho.is_finite() || rho == 1.0 || !sigma.is_finite() || sigma <= 0.0
            {
                return Err(error(
                    "mann-fixed-step-clock",
                    "mean-reverting OU transition is not representable",
                ));
            }
            transition.push([rho, sigma]);
        }
        let result = Self {
            definition: Arc::new(Definition {
                model: checkpoint.model,
                times: checkpoint.correlation_times_s,
                transition,
                dt_s: checkpoint.dt_s,
                origin_time_s: checkpoint.origin_time_s,
                origin_evolution_seed: checkpoint.origin_evolution_seed,
                evolution_seed: checkpoint.evolution_seed,
            }),
            step: checkpoint.step,
            amplitudes: checkpoint.amplitudes,
            trial: [zeros(n)?, zeros(n)?, zeros(n)?, zeros(n)?],
        };
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(error(
                "mann-fixed-step-budget",
                "allocated temporal capacities exceed the byte budget",
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;

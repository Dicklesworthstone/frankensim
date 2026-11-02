//! Bounded simultaneous continuous tones in the shared physical runtime.

mod power;

use super::super::{
    ExactZohKernel, ModalAcousticFrame, ModalAcousticState, ModalAcousticTimeError,
    ModalAcousticTimeModel, ModalAcousticWorkspace, ModalStepData, stage_modal_frame,
};
use super::{
    HarmonicForceTransition, advance_response, invalid, runtime::require_mechanical_observation,
};
use power::{TonePower, scaling};
use std::sync::Arc;

/// Hard bound on prescribed tones in one prepared modal bank. No truncation.
pub const MAX_HARMONIC_TONES: usize = 32;
const MAX_COORDINATES: usize = 2 + 2 * MAX_HARMONIC_TONES;

/// Explicit storage and work admission for simultaneous-tone preparation.
/// Work scores count fixed algebraic loop sizes, not elapsed time or CPU
/// instructions. They make expensive preparations/steps refuse before any
/// coefficient allocation; they are not a real-time scheduling certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModalToneBudget {
    /// Maximum input tones, additionally bounded by [`MAX_HARMONIC_TONES`].
    pub maximum_tones: usize,
    /// Upper allowance for owned coefficient payload plus peak preparation
    /// scratch [bytes]. Excludes allocator overhead, the already-owned model,
    /// and caller-owned forces/workspace. Step scratch is bounded on the stack.
    pub maximum_storage_bytes: usize,
    /// Maximum preparation work score: summed over modes, with `n=2+2K`
    /// for distinct frequencies and `d` resolved-interval doublings,
    /// `160*n^2 + d*(8*n^3+80*n^2) + K*(2000+80*d)`.
    pub maximum_preparation_work: usize,
    /// Maximum per-step work score: `modes*(16*n^2+128*input_tones)`.
    pub maximum_step_work: usize,
}

impl ModalToneBudget {
    /// Bounded reference ceilings; callers may lower them for their callback.
    #[must_use]
    pub const fn realtime_reference() -> Self {
        Self {
            maximum_tones: MAX_HARMONIC_TONES,
            maximum_storage_bytes: 64 * 1024 * 1024,
            maximum_preparation_work: 1_000_000_000,
            maximum_step_work: 10_000_000,
        }
    }
}

#[derive(Clone, Debug)]
struct PreparedToneMode {
    motion: Box<[HarmonicForceTransition]>,
    power: TonePower,
}

/// Prepared finite sum of continuous tones for every retained mode.
///
/// Input force rows use original caller order, in mass-normalized
/// `N/sqrt(kg)`: `F(t)=sum_j(Fc[j]*cos(Omega[j]*t)+Fs[j]*sin(Omega[j]*t))`.
/// Frequencies may be inharmonic, zero, repeated or arbitrarily close. Only
/// exactly equal frequencies are coalesced; close frequencies retain beats.
/// Each physical mode has one homogeneous state and complete cross-frequency
/// work/positive heat forms. Adding single-tone heats is never used.
///
/// Preparation allocates within explicit budgets. Motion, phase rotation and
/// physical steps allocate nothing and call no transcendentals. No tone or
/// mode is discarded. Numerical power allowances are estimates, not rigorous
/// interval certificates. Severe cancellation can lose relative accuracy;
/// neither force amplitudes nor physical damping are normalized to hide it.
///
/// Bound to the model's immutable kernel and its clones, just like
/// [`super::PreparedModalHarmonic`]. Legacy natural-frequency pressure
/// transfers must be zero; use a separate physical broadband observer.
#[derive(Clone, Debug)]
pub struct PreparedModalTones {
    owner: Arc<ExactZohKernel>,
    frequencies: Box<[f64]>,
    unique: Box<[f64]>,
    source_to_unique: Box<[usize]>,
    modes: Box<[PreparedToneMode]>,
    duration_s: f64,
}

impl PreparedModalTones {
    /// Prepare every requested tone and mode, including all power cross terms.
    /// All frequencies must lie within the model's nominal Nyquist guard.
    /// A longer integration window can contain many cycles; endpoints alone
    /// do not resolve observations within that window.
    ///
    /// # Errors
    /// Empty/oversized or invalid frequencies, invalid duration/budgets,
    /// storage/work admission, nonzero narrow-band transfers, or coefficients
    /// that cannot be represented in finite arithmetic. No model is mutated.
    pub fn try_new(
        model: &ModalAcousticTimeModel,
        angular_frequencies_rad_s: &[f64],
        duration_s: f64,
        budget: ModalToneBudget,
    ) -> Result<Self, ModalAcousticTimeError> {
        require_mechanical_observation(model)?;
        let count = angular_frequencies_rad_s.len();
        if count == 0 || count > MAX_HARMONIC_TONES || count > budget.maximum_tones {
            return Err(invalid(
                "tone count exceeds the explicit nonempty preparation limit",
            ));
        }
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid("tone duration must be positive and finite"));
        }
        let maximum =
            core::f64::consts::PI * model.budget.nyquist_guard_fraction / model.sample_period_s;
        let mut unique = [0.0; MAX_HARMONIC_TONES];
        let mut mapping = [0; MAX_HARMONIC_TONES];
        let mut distinct = 0;
        for (index, &frequency) in angular_frequencies_rad_s.iter().enumerate() {
            if !frequency.is_finite() || frequency < 0.0 || frequency > maximum {
                return Err(invalid(
                    "all tones must lie inside the model's nominal Nyquist guard",
                ));
            }
            mapping[index] =
                if let Some(found) = unique[..distinct].iter().position(|&x| x == frequency) {
                    found
                } else {
                    unique[distinct] = frequency;
                    distinct += 1;
                    distinct - 1
                };
        }
        admit_cost(model, &unique[..distinct], count, duration_s, budget)?;
        let mut modes = Vec::with_capacity(model.modes.len());
        for mode in &model.modes {
            let mut motion = Vec::with_capacity(distinct);
            for &frequency in &unique[..distinct] {
                motion.push(HarmonicForceTransition::try_new(
                    mode.angular_frequency_rad_s,
                    mode.damping_ratio,
                    frequency,
                    duration_s,
                )?);
            }
            let power = TonePower::try_new(
                mode.angular_frequency_rad_s,
                mode.damping_ratio,
                &unique[..distinct],
                duration_s,
            )?;
            modes.push(PreparedToneMode {
                motion: motion.into_boxed_slice(),
                power,
            });
        }
        Ok(Self {
            owner: Arc::clone(&model.transition_kernel),
            frequencies: angular_frequencies_rad_s.into(),
            unique: unique[..distinct].into(),
            source_to_unique: mapping[..count].into(),
            modes: modes.into_boxed_slice(),
            duration_s,
        })
    }

    /// Physical duration of each prepared interval [s].
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Original caller-ordered frequencies [rad/s], including duplicates/DC.
    #[must_use]
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }

    /// Number of force coefficients required in each modal tone row.
    #[must_use]
    pub fn tone_count(&self) -> usize {
        self.frequencies.len()
    }

    /// Number of retained physical modes.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }

    /// Rotate one original tone's interval-local coefficients after a successful
    /// step. Checkpoints include these caller-owned phases with modal state.
    ///
    /// # Errors
    /// Invalid tone index, nonfinite force, or unrepresentable rotated force.
    pub fn advance_force_phase(
        &self,
        tone_index: usize,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<(f64, f64), ModalAcousticTimeError> {
        let &unique = self
            .source_to_unique
            .get(tone_index)
            .ok_or_else(|| invalid("tone index must belong to the preparation"))?;
        self.modes[0].motion[unique].advance_force_phase(cosine_force, sine_force)
    }

    /// Propose one mode's exact endpoint for a continuous observer, without
    /// committing state or applying physical model budgets. Force rows contain
    /// exactly [`Self::tone_count`] values in original caller order.
    ///
    /// # Errors
    /// Invalid mode/row, nonfinite inputs, or unrepresentable motion.
    pub fn advance_mode(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        cosine_forces: &[f64],
        sine_forces: &[f64],
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        let mode = self
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("tone mode index must be admitted"))?;
        let (cosine, sine) = self.collect_forces(cosine_forces, sine_forces)?;
        mode.advance(
            state,
            &cosine[..self.unique.len()],
            &sine[..self.unique.len()],
        )
    }

    fn collect_forces(
        &self,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<([f64; MAX_HARMONIC_TONES], [f64; MAX_HARMONIC_TONES]), ModalAcousticTimeError>
    {
        let mut c = [Sum::default(); MAX_HARMONIC_TONES];
        let mut s = c;
        for forces in [cosine, sine] {
            if forces.len() != self.tone_count() {
                return Err(ModalAcousticTimeError::ForceCountMismatch {
                    expected: self.tone_count(),
                    found: forces.len(),
                });
            }
            if forces.iter().any(|x| !x.is_finite()) {
                return Err(invalid("tone force coefficients must be finite"));
            }
        }
        for (source, &target) in self.source_to_unique.iter().enumerate() {
            c[target].add(cosine[source]);
            s[target].add(sine[source]);
        }
        let c = c.map(Sum::value);
        let s = s.map(Sum::value);
        if c.iter().chain(&s).any(|x| !x.is_finite()) {
            return Err(invalid("coalesced tone forces must be representable"));
        }
        Ok((c, s))
    }
}

impl PreparedToneMode {
    fn advance(
        &self,
        state: ModalAcousticState,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<ModalAcousticState, ModalAcousticTimeError> {
        let mut held = Sum::default();
        let mut driven_q = Sum::default();
        let mut driven_v = Sum::default();
        let first = &self.motion[0];
        let mut has_held = false;
        for ((motion, &fc), &fs) in self.motion.iter().zip(cosine).zip(sine) {
            let fc = first.duration_s * fc;
            let fs = first.duration_s * fs;
            if !fc.is_finite() || !fs.is_finite() {
                return Err(invalid("scaled tone forces must be representable"));
            }
            let [qc, qs, vc, vs] = motion.force;
            driven_q.add(qc * fc);
            driven_q.add(qs * fs);
            driven_v.add(vc * fc);
            driven_v.add(vs * fs);
            if motion.cosine_remainder {
                has_held = true;
                held.add(fc);
            }
        }
        advance_response(
            state,
            first.duration_s,
            first.angular_frequency_rad_s,
            first.transition,
            has_held.then(|| held.value()),
            driven_q.value(),
            driven_v.value(),
        )
    }
}

impl ModalAcousticTimeModel {
    /// Advance all simultaneous tones with complete input-work/heat accounts.
    /// Force arrays are mode-major: `index = mode*tone_count + tone`, in the
    /// original preparation order. No internal phase advances; rotate external
    /// coefficients only after success. Existing model gates apply to every
    /// candidate, and a refusal leaves all accepted modal states unchanged.
    ///
    /// # Errors
    /// Foreign preparation, wrong/nonfinite arrays, nonzero narrow-band
    /// observer, incompatible workspace, state/energy budget or power balance.
    pub fn step_tones_into<'w>(
        &mut self,
        prepared: &PreparedModalTones,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        if !Arc::ptr_eq(&prepared.owner, &self.transition_kernel) {
            return Err(invalid(
                "tone preparation must belong to this model or its clone",
            ));
        }
        require_mechanical_observation(self)?;
        let expected = self.modes.len() * prepared.tone_count();
        for forces in [cosine_forces, sine_forces] {
            if forces.len() != expected {
                return Err(ModalAcousticTimeError::ForceCountMismatch {
                    expected,
                    found: forces.len(),
                });
            }
            if forces.iter().any(|x| !x.is_finite()) {
                return Err(invalid("tone force coefficients must be finite"));
            }
        }
        stage_modal_frame(
            &self.modes,
            &self.states,
            self.budget,
            workspace,
            "independent simultaneous-tone work, heat and modal energy disagree",
            "aggregate simultaneous-tone work, heat and modal energy disagree",
            |index, _mode, state| {
                let begin = index * prepared.tone_count();
                let end = begin + prepared.tone_count();
                let (c, s) = prepared
                    .collect_forces(&cosine_forces[begin..end], &sine_forces[begin..end])?;
                let c = &c[..prepared.unique.len()];
                let s = &s[..prepared.unique.len()];
                let mode = &prepared.modes[index];
                let next = mode.advance(state, c, s)?;
                let (work, work_allowance, loss, loss_allowance) =
                    mode.power.evaluate(state, c, s)?;
                Ok(ModalStepData {
                    next,
                    work,
                    work_allowance,
                    loss,
                    loss_allowance,
                })
            },
        )?;
        self.states.copy_from_slice(&workspace.candidate);
        Ok(&workspace.frame)
    }
}

fn admit_cost(
    model: &ModalAcousticTimeModel,
    frequencies: &[f64],
    input_tones: usize,
    h: f64,
    budget: ModalToneBudget,
) -> Result<(), ModalAcousticTimeError> {
    let k = frequencies.len();
    let n = 2 + 2 * k;
    let square = n * n;
    let count = model.modes.len();
    // The hard mode/tone caps keep these small integer products representable.
    // Include both retained forms, motion, shared maps, and preparation scratch.
    let storage = count
        * (core::mem::size_of::<PreparedToneMode>()
            + k * core::mem::size_of::<HarmonicForceTransition>()
            + 2 * square * 8)
        + 4 * square * 8
        + (10 * n + 20 * input_tones) * 8
        + core::mem::size_of::<PreparedModalTones>();
    let step_work = count * (16 * square + 128 * input_tones);
    let mut preparation = 0usize;
    for mode in &model.modes {
        let (_, _, _, _, doublings) = scaling(
            mode.angular_frequency_rad_s,
            mode.damping_ratio,
            frequencies,
            h,
        )?;
        let cost =
            160 * square + doublings * (8 * n * square + 80 * square) + k * (2000 + 80 * doublings);
        preparation = preparation
            .checked_add(cost)
            .ok_or_else(|| invalid("tone preparation work score overflow"))?;
    }
    for (what, value, limit) in [
        (
            "tone coefficient and preparation storage bytes",
            storage,
            budget.maximum_storage_bytes,
        ),
        (
            "tone preparation work score",
            preparation,
            budget.maximum_preparation_work,
        ),
        ("tone step work score", step_work, budget.maximum_step_work),
    ] {
        if value > limit {
            return Err(ModalAcousticTimeError::BudgetExceeded {
                what,
                value: value as f64,
                limit: limit as f64,
            });
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

//! Continuous forcing of the actual modal trajectory and retained receiver.
use super::{BroadbandModalError, BroadbandModalObserver, MAX_MODAL_AUGMENTED_MATRIX_ENTRIES};
use fs_modal::time::{
    MAX_HARMONIC_TONES, ModalAcousticState, ModalAcousticTimeError, ModalToneBudget,
    PreparedModalHarmonic, PreparedModalTones,
};
use std::sync::Arc;

#[derive(Clone, Debug)]
struct ForcingColumns {
    owner: Arc<()>,
    tone_count: usize,
    affine: bool,
    // Exact equal frequencies only, preserving original coefficient order.
    groups: Vec<Vec<usize>>,
    // Mode-major, then filter-state-major, then unique-frequency cosine/sine.
    columns: Vec<Vec<[f64; 2]>>,
    delta_columns: Vec<Vec<[f64; 2]>>,
}

/// One continuous force frequency on the original receiver's sample clock.
/// The retained acoustic state follows the full driven and transient motion.
/// Preparation allocates; stepping is allocation-free and source-bound.
#[derive(Clone, Debug)]
pub struct PreparedBroadbandHarmonic {
    mechanics: PreparedModalHarmonic,
    filters: ForcingColumns,
}

impl PreparedBroadbandHarmonic {
    /// Exact physical observation interval [s].
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.mechanics.duration_s()
    }
    /// Continuous drive angular frequency [rad/s], including DC.
    #[must_use]
    pub fn forcing_angular_frequency_rad_s(&self) -> f64 {
        self.mechanics.forcing_angular_frequency_rad_s()
    }
    /// Move interval-local force phase forward after a successful observation.
    /// A refused step must leave the caller's phase unchanged.
    pub fn advance_force_phase(
        &self,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), BroadbandModalError> {
        self.mechanics
            .advance_force_phase(cosine, sine)
            .map_err(BroadbandModalError::Modal)
    }
}

/// Simultaneous continuous tones through one original broadband receiver.
/// Forces are mode-major, with the original tone order inside each mode.
/// The mechanical owner retains cross-frequency work and viscous heat.
#[derive(Clone, Debug)]
pub struct PreparedBroadbandTones {
    mechanics: PreparedModalTones,
    filters: ForcingColumns,
}

impl PreparedBroadbandTones {
    /// Whether this source uses affine rotating-quadrature envelopes.
    #[must_use]
    pub fn is_affine(&self) -> bool {
        self.mechanics.is_affine()
    }
    /// Exact physical observation interval [s].
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.mechanics.duration_s()
    }
    /// Original continuous frequencies [rad/s]; equal frequencies remain visible.
    #[must_use]
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.mechanics.angular_frequencies_rad_s()
    }
    /// Original number of cosine/sine force coefficients for each mode.
    #[must_use]
    pub fn tone_count(&self) -> usize {
        self.filters.tone_count
    }
    /// Rotate one original tone after a successful complete observation step.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), BroadbandModalError> {
        self.mechanics
            .advance_force_phase(tone, cosine, sine)
            .map_err(BroadbandModalError::Modal)
    }
}

impl BroadbandModalObserver {
    /// Prepare continuous `Fc*cos(Omega*t) + Fs*sin(Omega*t)` modal loading.
    /// Each force pair starts at this interval's sample time, in exp(-i Omega t)
    /// convention. The interval is exactly the existing observer sample period.
    ///
    /// The original fitted continuous transfer, retained acoustic state, sampled
    /// propagation delay and modal budgets remain in force. Complex constant
    /// legacy quadrature is refused: it does not define a causal broadband
    /// response to a new forcing frequency and its startup transient. A real
    /// constant gain is a valid instantaneous transfer and remains available.
    pub fn prepare_harmonic(
        &self,
        omega: f64,
    ) -> Result<PreparedBroadbandHarmonic, BroadbandModalError> {
        let filters = self.admit_forcing_columns(&[omega], false)?;
        let mechanics =
            PreparedModalHarmonic::try_new(&self.model, omega, self.model.sample_period_s())
                .map_err(BroadbandModalError::Modal)?;
        let filters = self.prepare_forcing_columns(&[omega], filters)?;
        Ok(PreparedBroadbandHarmonic { mechanics, filters })
    }

    /// Prepare simultaneous continuous tones without enlarging a filter's
    /// dense exponential by the total number of tones. Each unique frequency
    /// supplies two exact forcing columns; all reuse the homogeneous block.
    /// Exact duplicates are summed before application, never near-frequency
    /// merged. The shared mechanical owner admits the explicit work/heat budget.
    pub fn prepare_tones(
        &self,
        omega: &[f64],
        budget: ModalToneBudget,
    ) -> Result<PreparedBroadbandTones, BroadbandModalError> {
        let filters = self.admit_forcing_columns(omega, false)?;
        let mechanics =
            PreparedModalTones::try_new(&self.model, omega, self.model.sample_period_s(), budget)
                .map_err(BroadbandModalError::Modal)?;
        let filters = self.prepare_forcing_columns(omega, filters)?;
        Ok(PreparedBroadbandTones { mechanics, filters })
    }

    /// Prepare continuous linear force envelopes on each carrier. Both force
    /// endpoints use the same interval-start phase basis. Source and retained
    /// acoustic filter follow the actual affine-times-sinusoid trajectory;
    /// this does not interpolate endpoint velocities or fade output pressure.
    pub fn prepare_affine_tones(
        &self,
        omega: &[f64],
        budget: ModalToneBudget,
    ) -> Result<PreparedBroadbandTones, BroadbandModalError> {
        let filters = self.admit_forcing_columns(omega, true)?;
        let mechanics = PreparedModalTones::try_new_affine(
            &self.model,
            omega,
            self.model.sample_period_s(),
            budget,
        )
        .map_err(BroadbandModalError::Modal)?;
        let filters = self.prepare_forcing_columns(omega, filters)?;
        Ok(PreparedBroadbandTones { mechanics, filters })
    }

    /// Advance the exact driven trajectory and its physical receiver together.
    /// Both force arrays contain one coefficient per retained mode. Mechanical
    /// states, published frame, filter history and delay commit only on success.
    /// The returned pressure is at the interval endpoint after sampled delay.
    pub fn step_harmonic(
        &mut self,
        prepared: &PreparedBroadbandHarmonic,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<f64, BroadbandModalError> {
        let (raw, pressure) = self.preview_forcing(
            &prepared.filters,
            cosine,
            sine,
            None,
            |i, state, c, s, _, _| prepared.mechanics.advance_mode(i, state, c[0], s[0]),
        )?;
        self.model
            .step_harmonic_into(
                &prepared.mechanics,
                cosine,
                sine,
                &mut self.candidate_workspace,
            )
            .map_err(BroadbandModalError::Modal)?;
        self.commit_observation(raw);
        Ok(pressure)
    }

    /// Advance all declared tones coherently on the original sample clock.
    /// Arrays are `[mode * tone_count + tone]`. Cross-tone interference is
    /// retained in pressure and in the shared mechanical work/heat accounting.
    /// A refusal advances no accepted state or caller-owned force phase.
    pub fn step_tones(
        &mut self,
        prepared: &PreparedBroadbandTones,
        cosine: &[f64],
        sine: &[f64],
    ) -> Result<f64, BroadbandModalError> {
        let (raw, pressure) = self.preview_forcing(
            &prepared.filters,
            cosine,
            sine,
            None,
            |i, state, c, s, _, _| prepared.mechanics.advance_mode(i, state, c, s),
        )?;
        self.model
            .step_tones_into(
                &prepared.mechanics,
                cosine,
                sine,
                &mut self.candidate_workspace,
            )
            .map_err(BroadbandModalError::Modal)?;
        self.commit_observation(raw);
        Ok(pressure)
    }

    /// Commit one affine force interval, including causal acoustic startup,
    /// sampled delay and original mechanical work/heat. Endpoint quadratures
    /// are mode-major and share the interval-start carrier basis. All accepted
    /// histories and caller phase remain unchanged on any refusal.
    pub fn step_affine_tones(
        &mut self,
        prepared: &PreparedBroadbandTones,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
    ) -> Result<f64, BroadbandModalError> {
        let (raw, pressure) = self.preview_forcing(
            &prepared.filters,
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
            |i, state, c0, s0, c1, s1| {
                prepared
                    .mechanics
                    .advance_mode_affine(i, state, c0, s0, c1, s1)
            },
        )?;
        self.model
            .step_tones_affine_into(
                &prepared.mechanics,
                cosine_start,
                sine_start,
                cosine_end,
                sine_end,
                &mut self.candidate_workspace,
            )
            .map_err(BroadbandModalError::Modal)?;
        self.commit_observation(raw);
        Ok(pressure)
    }

    fn admit_forcing_columns(
        &self,
        omega: &[f64],
        affine: bool,
    ) -> Result<ForcingColumns, BroadbandModalError> {
        let maximum = core::f64::consts::PI * self.model.budget().nyquist_guard_fraction
            / self.model.sample_period_s();
        if omega.is_empty()
            || omega
                .iter()
                .any(|w| !w.is_finite() || *w < 0.0 || *w > maximum)
        {
            return Err(BroadbandModalError::Invalid(
                "continuous forcing must lie inside the original modal Nyquist guard",
            ));
        }
        if omega.len() > MAX_HARMONIC_TONES {
            return Err(BroadbandModalError::Cap {
                what: "continuous observer tones",
                found: omega.len(),
                cap: MAX_HARMONIC_TONES,
            });
        }
        let count = self.model.modes().len().checked_mul(omega.len()).ok_or(
            BroadbandModalError::Invalid("continuous observer force shape overflow"),
        )?;
        if count > MAX_MODAL_AUGMENTED_MATRIX_ENTRIES {
            return Err(BroadbandModalError::Cap {
                what: "continuous observer force coefficients",
                found: count,
                cap: MAX_MODAL_AUGMENTED_MATRIX_ENTRIES,
            });
        }
        if self
            .constant
            .as_ref()
            .is_some_and(|h| h.iter().any(|z| z.im != 0.0))
        {
            return Err(BroadbandModalError::Invalid(
                "continuous forcing requires a real constant gain or a fitted causal broadband transfer",
            ));
        }
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for (index, &w) in omega.iter().enumerate() {
            if let Some(group) = groups.iter_mut().find(|g| omega[g[0]] == w) {
                group.push(index);
            } else {
                groups.push(vec![index]);
            }
        }
        let entries = self
            .filters
            .iter()
            .try_fold(0usize, |sum, filter| {
                let n = filter.n + if affine { 3 } else { 1 };
                let new = if filter.state.is_empty() {
                    0
                } else {
                    groups.len().checked_mul(n * n)?
                };
                sum.checked_add(filter.n * filter.n)?.checked_add(new)
            })
            .ok_or(BroadbandModalError::Invalid(
                "continuous augmented observer shape overflow",
            ))?;
        if entries > MAX_MODAL_AUGMENTED_MATRIX_ENTRIES {
            return Err(BroadbandModalError::Cap {
                what: "continuous augmented observer matrix entries",
                found: entries,
                cap: MAX_MODAL_AUGMENTED_MATRIX_ENTRIES,
            });
        }
        Ok(ForcingColumns {
            owner: Arc::clone(&self.forcing_owner),
            tone_count: omega.len(),
            affine,
            groups,
            columns: Vec::new(),
            delta_columns: Vec::new(),
        })
    }

    fn prepare_forcing_columns(
        &self,
        omega: &[f64],
        mut prepared: ForcingColumns,
    ) -> Result<ForcingColumns, BroadbandModalError> {
        let groups = &prepared.groups;
        let dt = self.model.sample_period_s();
        let mut columns = Vec::with_capacity(self.filters.len());
        let mut delta_columns = if prepared.affine {
            Vec::with_capacity(self.filters.len())
        } else {
            Vec::new()
        };
        for ((filter, mode), &damping) in self.filters.iter().zip(self.model.modes())
            .zip(self.model.damping_rates_per_s()) {
            let ss = filter.model.state_space();
            let mut response = vec![[0.0; 2]; ss.n * groups.len()];
            let mut delta_response = vec![
                [0.0; 2];
                if prepared.affine {
                    ss.n * groups.len()
                } else {
                    0
                }
            ];
            if ss.n != 0 {
                // [scaled displacement, velocity, filter state, Fc, Fs].
                // Fc' = Omega*Fs, Fs' = -Omega*Fc; force(t)=Fc(t).
                let n = ss.n + if prepared.affine { 6 } else { 4 };
                let source = ss.n + 2;
                let mut a = vec![0.0; n * n];
                let w = mode.angular_frequency_rad_s;
                a[1] = filter.q_scale * dt;
                a[n] = -w * w / filter.q_scale * dt;
                a[n + 1] = -damping * dt;
                a[n + source] = dt;
                for i in 0..ss.n {
                    a[(i + 2) * n + 1] = ss.b[i] * dt;
                    for j in 0..ss.n {
                        a[(i + 2) * n + j + 2] = ss.a[i * ss.n + j] * dt;
                    }
                }
                for (tone, group) in groups.iter().enumerate() {
                    a[source * n + source + 1] = omega[group[0]] * dt;
                    a[(source + 1) * n + source] = -omega[group[0]] * dt;
                    if prepared.affine {
                        // Normalized endpoint deltas, never force/time rates.
                        // C'=Omega*S+Dc/h, S'=-Omega*C+Ds/h.
                        a[source * n + source + 2] = 1.0;
                        a[(source + 1) * n + source + 3] = 1.0;
                        a[(source + 2) * n + source + 3] = omega[group[0]] * dt;
                        a[(source + 3) * n + source + 2] = -omega[group[0]] * dt;
                    }
                    let transition = fs_phs::linear_exact::exponential(&a, n)
                        .map_err(BroadbandModalError::Exact)?;
                    for row in 0..ss.n {
                        response[row * groups.len() + tone] = [
                            transition[(row + 2) * n + source],
                            transition[(row + 2) * n + source + 1],
                        ];
                        if prepared.affine {
                            delta_response[row * groups.len() + tone] = [
                                transition[(row + 2) * n + source + 2],
                                transition[(row + 2) * n + source + 3],
                            ];
                        }
                    }
                }
            }
            columns.push(response);
            if prepared.affine {
                delta_columns.push(delta_response);
            }
        }
        prepared.columns = columns;
        prepared.delta_columns = delta_columns;
        Ok(prepared)
    }

    fn preview_forcing(
        &mut self,
        prepared: &ForcingColumns,
        cosine: &[f64],
        sine: &[f64],
        ends: Option<(&[f64], &[f64])>,
        mut advance: impl FnMut(
            usize,
            ModalAcousticState,
            &[f64],
            &[f64],
            &[f64],
            &[f64],
        ) -> Result<ModalAcousticState, ModalAcousticTimeError>,
    ) -> Result<(f64, f64), BroadbandModalError> {
        if !Arc::ptr_eq(&prepared.owner, &self.forcing_owner) || prepared.affine != ends.is_some() {
            return Err(BroadbandModalError::Invalid(
                "forcing preparation must match this broadband receiver and envelope law",
            ));
        }
        let tones = prepared.tone_count;
        let count = self.model.modes().len() * tones;
        if cosine.len() != count
            || sine.len() != count
            || cosine.iter().chain(sine).any(|x| !x.is_finite())
        {
            return Err(BroadbandModalError::Invalid(
                "finite complete mode-major continuous force arrays required",
            ));
        }
        if ends.is_some_and(|(c, s)| {
            c.len() != count
                || s.len() != count
                || c.iter().chain(s).any(|value| !value.is_finite())
        }) {
            return Err(BroadbandModalError::Invalid(
                "finite complete affine force endpoint arrays required",
            ));
        }
        let mut raw = 0.0;
        for (index, &state) in self.model.states().iter().enumerate() {
            let c = &cosine[index * tones..(index + 1) * tones];
            let s = &sine[index * tones..(index + 1) * tones];
            let (ce, se) = ends.map_or((c, s), |(ce, se)| {
                (
                    &ce[index * tones..(index + 1) * tones],
                    &se[index * tones..(index + 1) * tones],
                )
            });
            let next = advance(index, state, c, s, ce, se).map_err(BroadbandModalError::Modal)?;
            if let Some(h) = &self.constant {
                raw += h[index].re * next.velocity_m_sqrt_kg_per_s;
                continue;
            }
            let filter = &mut self.filters[index];
            for row in 0..filter.state.len() {
                let r = (row + 2) * filter.n;
                let mut x = filter.transition[r] * filter.q_scale * state.displacement_m_sqrt_kg
                    + filter.transition[r + 1] * state.velocity_m_sqrt_kg_per_s;
                for (j, &value) in filter.state.iter().enumerate() {
                    x += filter.transition[r + j + 2] * value;
                }
                for (tone, group) in prepared.groups.iter().enumerate() {
                    let column = prepared.columns[index][row * prepared.groups.len() + tone];
                    x += column[0] * grouped_force(c, group) + column[1] * grouped_force(s, group);
                    if prepared.affine {
                        let delta =
                            prepared.delta_columns[index][row * prepared.groups.len() + tone];
                        // Match mechanics: independently coalesce both
                        // endpoints before forming a finite envelope delta.
                        x += delta[0] * (grouped_force(ce, group) - grouped_force(c, group))
                            + delta[1] * (grouped_force(se, group) - grouped_force(s, group));
                    }
                }
                if !x.is_finite() {
                    return Err(BroadbandModalError::NonFinite {
                        what: "continuously forced filter state",
                    });
                }
                filter.candidate[row] = x;
            }
            let mut y = filter.d * next.velocity_m_sqrt_kg_per_s;
            for (&coefficient, &value) in filter.c.iter().zip(&filter.candidate) {
                y += coefficient * value;
            }
            raw += y;
        }
        Ok((raw, self.preview_pressure(raw)?))
    }
}

// Match the mechanical owner's original-order compensated duplicate-tone sum.
fn grouped_force(forces: &[f64], group: &[usize]) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for &index in group {
        let value = forces[index];
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
#[path = "forcing_tests.rs"]
mod tests;

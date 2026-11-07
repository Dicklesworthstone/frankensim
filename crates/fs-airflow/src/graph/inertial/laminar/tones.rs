//! Continuous pressure waves on the original developed Stokes pipe.
//!
//! The native exact tone propagator retains every viscous-memory coordinate
//! and independently integrates pressure work and irreversible wall heat,
//! including cross-frequency terms. No sampled or endpoint pressure replaces
//! the forcing inside an interval. The finite wall-spectrum approximation and
//! the original incompressible, Newtonian, developed-flow assumptions remain.

use fs_phs::linear_exact::{DenseExponentialBudget, PreparedExactTones};

use super::*;

/// Pressure DIFFERENCE across the pipe at one prepared angular frequency.
/// At the start of an interval, `dp(t) = cosine*cos(omega*t) + sine*sin(omega*t)`.
/// DC is represented by `omega=0` and an exactly zero sine coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarPressureTone {
    pub cosine: Pressure,
    pub sine: Pressure,
}

impl<M: FluidTransport> LaminarPipe<M> {
    /// Prepare continuous pulsatile forcing on this pipe's unchanged history.
    ///
    /// Frequencies are nonnegative, strictly increasing rad/s. Every nonzero
    /// frequency must lie inside the original wall-response admission band;
    /// DC uses the exact Poiseuille limit. The caller supplies the complete
    /// native byte/work/dimension budget. No modes or tones are discarded.
    ///
    /// Preparation caches exponentials, work/heat integrals and phase rotations.
    /// Successful warm steps allocate nothing and evaluate no transcendentals.
    /// This removes pressure-hold error; it does not certify arbitrary startup
    /// against the radial continuum or infer transition from a Reynolds ceiling.
    pub fn prepare_tones(
        &self,
        cx: &Cx<'_>,
        angular_frequencies_rad_s: &[f64],
        duration_s: f64,
        budget: DenseExponentialBudget,
    ) -> Result<PreparedLaminarTones<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        if !positive(duration_s)
            || angular_frequencies_rad_s.is_empty()
            || angular_frequencies_rad_s
                .iter()
                .any(|w| !w.is_finite() || *w < 0.0)
            || angular_frequencies_rad_s.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(LaminarPipeError::Invalid(
                "positive interval and increasing nonnegative pressure frequencies required",
            ));
        }
        for &omega in angular_frequencies_rad_s {
            checkpoint(cx)?;
            if omega > 0.0 {
                self.impedance_at(omega)?;
            }
        }
        let native = PreparedExactTones::new(
            &self.system,
            angular_frequencies_rad_s,
            duration_s,
            LinearSplitBudget {
                max_modes: self.flow_scales.len(),
                max_omega_dt: 1.0,
            },
            budget,
        )?;
        let count = angular_frequencies_rad_s.len();
        let prepared = PreparedLaminarTones {
            pipe: self,
            native,
            cosine: vec![0.0; 2 * count],
            sine: vec![0.0; 2 * count],
            pending: vec![0.0; self.flow_scales.len()],
            pending_pressure: vec![
                LaminarPressureTone {
                    cosine: Pressure::new(0.0),
                    sine: Pressure::new(0.0),
                };
                count
            ],
        };
        checkpoint(cx)?;
        Ok(prepared)
    }
}

/// Cached continuous-pressure trajectory, work, heat and source-phase advance.
/// Accepted pipe history and pressure phases are published together.
pub struct PreparedLaminarTones<'a, M> {
    pipe: &'a LaminarPipe<M>,
    native: PreparedExactTones<'a>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    pending: Vec<f64>,
    pending_pressure: Vec<LaminarPressureTone>,
}

impl<M: FluidTransport> PreparedLaminarTones<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }

    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.native.angular_frequencies_rad_s()
    }

    /// Advance a complete continuous-pressure interval and its source phases.
    ///
    /// Pressure coefficients describe the start of THIS interval, at the
    /// accepted state's clock. On success they are rotated to its end; retain
    /// them with the physical state for restart. Reusing the same buffer on
    /// successive calls therefore produces one continuous forcing history.
    /// A source, regime, arithmetic or cancellation refusal changes neither.
    ///
    /// The sum of tone amplitudes bounds the pressure throughout the interval.
    /// The cooperative mean/profile-flow equations then bound every interior
    /// Reynolds number, including hidden peaks and reversals between endpoints.
    /// This conservative envelope may refuse a strongly cancelling waveform.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarPipeState,
        pressures: &mut [LaminarPressureTone],
    ) -> Result<LaminarPipeStep, LaminarPipeError> {
        checkpoint(cx)?;
        if !Arc::ptr_eq(&state.identity, &self.pipe.identity)
            || pressures.len() != self.native.tone_count()
        {
            return Err(LaminarPipeError::Invalid(
                "original pipe history and complete pressure tones required",
            ));
        }
        let duration_s = self.native.duration();
        let next_time = state.accepted_time_s + duration_s;
        if !next_time.is_finite() || next_time <= state.accepted_time_s {
            return Err(LaminarPipeError::Invalid("unrepresentable pipe clock"));
        }
        let mut maximum_pressure = 0.0;
        for (j, pressure) in pressures.iter().enumerate() {
            let c = pressure.cosine.value();
            let s = pressure.sine.value();
            if !c.is_finite()
                || !s.is_finite()
                || (self.native.angular_frequencies_rad_s()[j] == 0.0 && s != 0.0)
            {
                return Err(LaminarPipeError::Invalid(
                    "finite pressure coefficients and zero DC sine required",
                ));
            }
            maximum_pressure += det::hypot(c, s);
            // The second native pressure port stays zero. Common pressure is
            // an exact null direction and need not enter dense arithmetic.
            self.cosine[j] = c;
            self.sine[j] = s;
            let (next_c, next_s) = self.native.advance_force_phase(j, c, s)?;
            self.pending_pressure[j] = LaminarPressureTone {
                cosine: Pressure::new(next_c),
                sine: Pressure::new(next_s),
            };
        }
        if !maximum_pressure.is_finite() {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable pressure envelope",
            ));
        }
        let initial_max = state
            .coordinates
            .iter()
            .zip(&self.pipe.flow_scales)
            .map(|(x, scale)| (x / scale).abs())
            .fold(0.0_f64, f64::max);
        let maximum = initial_max
            .max(maximum_pressure / self.pipe.steady_resistance)
            .min(initial_max + (maximum_pressure / self.pipe.inertance) * duration_s);
        if !maximum.is_finite() || (initial_max == 0.0 && maximum_pressure > 0.0 && maximum == 0.0)
        {
            return Err(LaminarPipeError::Invalid("unrepresentable flow envelope"));
        }
        let reynolds_envelope = self.pipe.admit_flow(maximum)?;
        let initial_volume_flow_m3_s = state.volume_flow_m3_s();
        let initial_mean_kinetic_energy_j = state.mean_kinetic_energy_j();
        let initial_profile_memory_energy_j = state.profile_memory_energy_j();
        let mut output = [0.0; 2];
        let native = self.native.step_into(
            &state.coordinates,
            &self.cosine,
            &self.sine,
            &mut self.pending,
            &mut output,
        )?;
        let final_volume_flow_m3_s = self.pending[0] / self.pipe.flow_scales[0];
        let final_mean_kinetic_energy_j = 0.5 * self.pending[0].powi(2);
        let final_profile_memory_energy_j = energy(&self.pending[1..]);
        let displaced_volume_m3 = duration_s * output[0];
        if ![
            final_volume_flow_m3_s,
            final_mean_kinetic_energy_j,
            final_profile_memory_energy_j,
            displaced_volume_m3,
        ]
        .into_iter()
        .all(f64::is_finite)
            || (output[0] != 0.0 && displaced_volume_m3 == 0.0)
            || output[1] != -output[0]
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable flow, energy or port continuity",
            ));
        }
        checkpoint(cx)?;
        state.coordinates.copy_from_slice(&self.pending);
        state.accepted_time_s = next_time;
        pressures.copy_from_slice(&self.pending_pressure);
        Ok(LaminarPipeStep {
            duration_s,
            initial_volume_flow_m3_s,
            final_volume_flow_m3_s,
            average_volume_flow_m3_s: output[0],
            displaced_volume_m3,
            reynolds_envelope,
            initial_mean_kinetic_energy_j,
            final_mean_kinetic_energy_j,
            initial_profile_memory_energy_j,
            final_profile_memory_energy_j,
            native,
        })
    }
}

#[cfg(test)]
mod tests;

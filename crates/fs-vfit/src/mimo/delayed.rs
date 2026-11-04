//! Exact held-input runtime for entrywise retarded multiport impedances.
//!
//! For `tau_ab / dt = d + f`, a step sees the previous delayed input for
//! `f*dt`, then the next delayed input for `(1-f)*dt`. Both pole propagation
//! and output integrals are evaluated analytically on those two intervals.
//! No delay is rounded, and no interpolation of sampled pressures is used.
//! The resulting interval-mean effort exchanges work with the caller's held
//! flow. Restored continuous and actual sampled checks remain **GridOnly**;
//! this runtime does not turn a transcendental delay into a Hamiltonian proof.

use super::{DelayedMimoModel, MimoPassivityReport, check_delayed_mimo_passivity,
    hermitian_part, zoh_mean_phi, zoh_pole_denominator};
use super::load::{MimoLoadError, MimoLoadLimits, validate_source};
use crate::model::PoleTerm;
use crate::passivity::CertificateClass;
use fs_math::{c64::C64, det};

/// Separate limits for the source model, delayed pole coordinates and input
/// rings. `source.max_memory_bytes` covers their combined retained payload.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DelayedMimoLoadLimits {
    pub source: MimoLoadLimits,
    /// Real delayed state dimension, counting a conjugate pair twice per entry.
    /// At most 65536; this is separate from the original common-pole dimension.
    pub max_delayed_states: usize,
    /// Total retained held-flow values across input-port rings, at most 2^24.
    pub max_history_values: usize,
}

#[derive(Clone, Debug)]
struct DelayedPole {
    decay: C64,
    pair: bool,
    mean_residues: Vec<C64>,
    old_gains: Vec<C64>,
    new_gains: Vec<C64>,
}

/// A proper delayed impedance with exact fractional-arrival interval means.
///
/// Accepted poles and input histories advance together. Clone this owner for a
/// larger coupled transaction. Stepping allocates nothing; malformed or
/// nonfinite candidates leave accepted history and caller output unchanged.
#[derive(Clone, Debug)]
pub struct PreparedDelayedMimoLoad {
    source: DelayedMimoModel,
    limits: DelayedMimoLoadLimits,
    time_step_s: f64,
    integer_delays: Vec<usize>,
    terms: Vec<DelayedPole>,
    old_direct: Vec<f64>,
    new_direct: Vec<f64>,
    mean_direct: Vec<f64>,
    state: Vec<C64>,
    candidate_state: Vec<C64>,
    history: Vec<f64>,
    history_length: usize,
    cursor: usize,
    continuous: MimoPassivityReport,
    sampled: MimoPassivityReport,
}

fn exponential(pole: C64, duration: f64) -> (C64, C64, C64) {
    let argument = pole.scale(duration);
    let (s, c) = det::sin_cos(argument.im);
    let exponential = C64::new(c, s).scale(det::exp(argument.re));
    let (phi1, phi2) = zoh_mean_phi(argument, exponential);
    (exponential, phi1, phi2)
}
fn finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }
fn real_output(z: C64, pair: bool) -> f64 { if pair { 2.0 * z.re } else { z.re } }

impl PreparedDelayedMimoLoad {
    /// Admit the complete source and its restored delays before allocating
    /// runtime histories. The delay-free residual need not itself be passive.
    /// Both restored continuous and actual sampled mean-output diagnostics must
    /// pass; neither report is promoted beyond its original GridOnly class.
    pub fn new(
        source: &DelayedMimoModel, time_step_s: f64, limits: DelayedMimoLoadLimits,
    ) -> Result<Self, MimoLoadError> {
        let base_bytes = validate_source(&source.model, time_step_s, limits.source)?;
        let n = source.model.ports;
        let entries = n * n;
        if source.delays.len() != entries
            || source.delays.iter().any(|t| !t.is_finite() || *t < 0.0)
            || (source.model.reciprocal && (0..entries).any(|j| {
                source.delays[j] != source.delays[(j % n) * n + j / n]
            }))
            || limits.max_delayed_states > 65_536
            || !(1..=1 << 24).contains(&limits.max_history_values)
        {
            return Err(MimoLoadError("delayed MIMO requires finite causal delays and bounded exact shapes"));
        }
        let delayed_states = source.model.order().checked_mul(entries)
            .ok_or(MimoLoadError("delayed MIMO state size overflow"))?;
        if delayed_states > limits.max_delayed_states {
            return Err(MimoLoadError("delayed MIMO pole coordinates exceed their state allowance"));
        }
        let mut history_length = 1usize;
        for &delay in &source.delays {
            let samples = delay / time_step_s;
            if !samples.is_finite() || samples.floor() >= limits.max_history_values as f64 {
                return Err(MimoLoadError("delayed MIMO arrival exceeds its history allowance"));
            }
            history_length = history_length.max(samples.floor() as usize + 1);
        }
        let history_values = history_length.checked_mul(n)
            .ok_or(MimoLoadError("delayed MIMO history size overflow"))?;
        if history_values > limits.max_history_values {
            return Err(MimoLoadError("delayed MIMO input rings exceed their history allowance"));
        }
        let complex_states = source.model.poles.len() * entries;
        // Base source/report allowance plus three complex coefficient banks,
        // accepted/candidate coordinates, delay/direct arrays and input rings.
        // Includes conservative padding; cold eigensolver scratch is separate.
        let required = base_bytes.checked_add(complex_states * 80)
            .and_then(|b| b.checked_add(entries * 48 + history_values * 8))
            .ok_or(MimoLoadError("delayed MIMO retained payload overflow"))?;
        if required > limits.source.max_memory_bytes {
            return Err(MimoLoadError("delayed MIMO payload exceeds its memory allowance"));
        }
        let band = limits.source.band_rad_s;
        let continuous = check_delayed_mimo_passivity(source, (band[0], band[1]))
            .map_err(|_| MimoLoadError("restored delayed MIMO passivity admission failed"))?;
        if !continuous.passive || !continuous.worst.0.is_finite() {
            return Err(MimoLoadError("restored delayed MIMO impedance has a passivity violation"));
        }
        let integer_delays: Vec<usize> = source.delays.iter()
            .map(|t| (t / time_step_s).floor() as usize).collect();
        let fractions: Vec<f64> = source.delays.iter().zip(&integer_delays)
            .map(|(t, d)| t / time_step_s - *d as f64).collect();
        let mut old_direct: Vec<f64> = source.model.d.iter().zip(&fractions)
            .map(|(d, f)| d * f).collect();
        let mut new_direct: Vec<f64> = source.model.d.iter().zip(&fractions)
            .map(|(d, f)| d * (1.0 - f)).collect();
        let mut terms = Vec::with_capacity(source.model.poles.len());
        for (pole, residues) in source.model.poles.iter().zip(&source.model.residues) {
            let (p, pair) = match pole {
                PoleTerm::Real { pole, .. } => (C64::from_re(*pole), false),
                PoleTerm::Pair { pole, .. } => (*pole, true),
            };
            let (decay, mean_gain, _) = exponential(p, time_step_s);
            let mut term = DelayedPole { decay, pair,
                mean_residues: residues.iter().map(|r| *r * mean_gain).collect(),
                old_gains: Vec::with_capacity(entries), new_gains: Vec::with_capacity(entries),
            };
            for (j, &f) in fractions.iter().enumerate() {
                let old_duration = f * time_step_s;
                let new_duration = (1.0 - f) * time_step_s;
                let (_, old_phi1, old_phi2) = exponential(p, old_duration);
                let (new_decay, new_phi1, new_phi2) = exponential(p, new_duration);
                let old_gain = old_phi1.scale(old_duration);
                term.old_gains.push(new_decay * old_gain);
                term.new_gains.push(new_phi1.scale(new_duration));
                // Integral of the old piece includes its continued state
                // throughout the new piece. Scale before multiplying dt twice.
                let old_mean = old_phi2.scale(f * old_duration)
                    + (new_phi1 * old_phi1).scale(f * new_duration);
                let new_mean = new_phi2.scale((1.0 - f) * new_duration);
                old_direct[j] += real_output(residues[j] * old_mean, pair);
                new_direct[j] += real_output(residues[j] * new_mean, pair);
            }
            if !finite(decay) || term.mean_residues.iter().chain(&term.old_gains)
                .chain(&term.new_gains).any(|z| !finite(*z)) {
                return Err(MimoLoadError("delayed MIMO realization coefficients left the finite set"));
            }
            terms.push(term);
        }
        if old_direct.iter().chain(&new_direct).any(|v| !v.is_finite()) {
            return Err(MimoLoadError("delayed MIMO mean coefficients left the finite set"));
        }
        let mean_direct = new_direct.iter().zip(&integer_delays)
            .map(|(d, lag)| if *lag == 0 { *d } else { 0.0 }).collect();
        let mut prepared = Self {
            source: source.clone(), limits, time_step_s, integer_delays, terms,
            old_direct, new_direct, mean_direct,
            state: vec![C64::ZERO; complex_states],
            candidate_state: vec![C64::ZERO; complex_states],
            history: vec![0.0; history_values], history_length, cursor: 0,
            sampled: continuous.clone(), continuous,
        };
        prepared.sampled = prepared.sampled_passivity(limits.source.sampled_passivity_points);
        if !prepared.sampled.passive {
            return Err(MimoLoadError("actual sampled delayed MIMO impedance has a passivity violation"));
        }
        Ok(prepared)
    }

    #[must_use]
    pub fn source_model(&self) -> &DelayedMimoModel { &self.source }
    #[must_use]
    pub fn limits(&self) -> &DelayedMimoLoadLimits { &self.limits }
    #[must_use]
    pub fn ports(&self) -> usize { self.source.model.ports }
    #[must_use]
    pub fn time_step_s(&self) -> f64 { self.time_step_s }
    /// Restored continuous and actual sampled interval-mean reports, in order.
    #[must_use]
    pub fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        (&self.continuous, &self.sampled)
    }
    /// Current held-flow coefficient for the simultaneous interval solve.
    /// Only entries whose first arrival is within this interval contribute.
    #[must_use]
    pub fn mean_direct_term(&self) -> &[f64] { &self.mean_direct }
    /// Accepted delayed states, in pole-term then row-major matrix-entry order.
    #[must_use]
    pub fn pole_states(&self) -> &[C64] { &self.state }
    /// Accepted input-port-major ring storage; the cursor is the next write slot
    /// in each equal-length input ring. Initial prehistory is identically zero.
    #[must_use]
    pub fn input_history(&self) -> &[f64] { &self.history }
    #[must_use]
    pub fn history_cursor(&self) -> usize { self.cursor }

    fn past_input(&self, input: usize, lag: usize) -> f64 {
        debug_assert!(lag > 0 && lag <= self.history_length);
        let slot = (self.cursor + self.history_length - lag) % self.history_length;
        self.history[input * self.history_length + slot]
    }

    /// Accepted history's interval-mean effort, leaving state and output intact
    /// on refusal. Add mean_direct_term()*input for the complete coupled effort.
    pub fn mean_history_output_into(&self, output: &mut [f64]) -> Result<(), MimoLoadError> {
        let n = self.ports();
        if output.len() != n {
            return Err(MimoLoadError("delayed MIMO history requires one output per port"));
        }
        let mut candidate = [0.0; 128];
        for (a, effort) in candidate[..n].iter_mut().enumerate() {
            for b in 0..n {
                let j = a * n + b;
                let lag = self.integer_delays[j];
                *effort += self.old_direct[j] * self.past_input(b, lag + 1);
                if lag > 0 { *effort += self.new_direct[j] * self.past_input(b, lag); }
            }
        }
        let entries = n * n;
        for (k, term) in self.terms.iter().enumerate() {
            for (a, effort) in candidate[..n].iter_mut().enumerate() {
                let mut sum = C64::ZERO;
                for b in 0..n {
                    let j = a * n + b;
                    sum = sum + term.mean_residues[j] * self.state[k * entries + j];
                }
                *effort += real_output(sum, term.pair);
            }
        }
        if candidate[..n].iter().any(|v| !v.is_finite()) {
            return Err(MimoLoadError("delayed MIMO history output left the finite set"));
        }
        output.copy_from_slice(&candidate[..n]);
        Ok(())
    }

    /// Advance exact retarded pole states and held-flow rings, returning the
    /// interval-mean effort. Refusal changes no accepted history or output.
    pub fn step_mean_into(&mut self, input: &[f64], output: &mut [f64]) -> Result<(), MimoLoadError> {
        let n = self.ports();
        if input.len() != n || output.len() != n || input.iter().any(|v| !v.is_finite()) {
            return Err(MimoLoadError("delayed MIMO step requires one finite input and output per port"));
        }
        let mut candidate = [0.0; 128];
        self.mean_history_output_into(&mut candidate[..n])?;
        for (a, effort) in candidate[..n].iter_mut().enumerate() {
            *effort += self.mean_direct[a * n..(a + 1) * n].iter()
                .zip(input).map(|(d, u)| d * u).sum::<f64>();
        }
        if candidate[..n].iter().any(|v| !v.is_finite()) || !self.preview_state(input) {
            return Err(MimoLoadError("delayed MIMO candidate output or state left the finite set"));
        }
        self.accept_state(input);
        output.copy_from_slice(&candidate[..n]);
        Ok(())
    }

    pub(crate) fn preview_state(&mut self, input: &[f64]) -> bool {
        let n = self.ports();
        if input.len() != n || input.iter().any(|v| !v.is_finite()) { return false; }
        let entries = n * n;
        for (k, term) in self.terms.iter().enumerate() {
            for j in 0..entries {
                let b = j % n;
                let lag = self.integer_delays[j];
                let old = self.past_input(b, lag + 1);
                let new = if lag == 0 { input[b] } else { self.past_input(b, lag) };
                self.candidate_state[k * entries + j] = term.decay * self.state[k * entries + j]
                    + term.old_gains[j].scale(old) + term.new_gains[j].scale(new);
            }
        }
        self.candidate_state.iter().all(|z| finite(*z))
    }

    pub(crate) fn accept_state(&mut self, input: &[f64]) {
        self.state.copy_from_slice(&self.candidate_state);
        for (b, &u) in input.iter().enumerate() {
            self.history[b * self.history_length + self.cursor] = u;
        }
        self.cursor = (self.cursor + 1) % self.history_length;
    }

    /// Actual sampled transfer of the analytically integrated held-input map.
    /// Fractional delays generally differ from multiplying an ordinary sampled
    /// ZOH transfer by exp(-i*omega*tau), because the held signal has images.
    pub fn eval_mean_discrete_iw(&self, omega: f64) -> Result<Vec<C64>, MimoLoadError> {
        if !omega.is_finite() || omega < 0.0 || omega > core::f64::consts::PI / self.time_step_s {
            return Err(MimoLoadError("delayed MIMO frequency must lie between DC and Nyquist"));
        }
        let phase = omega * self.time_step_s;
        let (s, c) = det::sin_cos(phase);
        let z = C64::new(c, s);
        let mut response: Vec<C64> = self.new_direct.iter().zip(&self.old_direct)
            .map(|(new, old)| C64::from_re(*new) + z.conj().scale(*old)).collect();
        for (term, source_pole) in self.terms.iter().zip(&self.source.model.poles) {
            let pole = match source_pole {
                PoleTerm::Real { pole, .. } => C64::from_re(*pole),
                PoleTerm::Pair { pole, .. } => *pole,
            };
            let denominator = zoh_pole_denominator(pole, omega, self.time_step_s);
            let conjugate_denominator = if term.pair {
                zoh_pole_denominator(pole.conj(), omega, self.time_step_s)
            } else { C64::ONE };
            for (j, entry) in response.iter_mut().enumerate() {
                let first = term.mean_residues[j]
                    * (term.new_gains[j] + z.conj() * term.old_gains[j]) * z.conj()
                    / denominator;
                *entry = *entry + if term.pair {
                    first + term.mean_residues[j].conj()
                        * (term.new_gains[j].conj() + z.conj() * term.old_gains[j].conj())
                        * z.conj() / conjugate_denominator
                } else { first };
            }
        }
        for (entry, lag) in response.iter_mut().zip(&self.integer_delays) {
            let (s, c) = det::sin_cos(phase * *lag as f64);
            *entry = *entry * C64::new(c, -s);
        }
        if response.iter().any(|z| !finite(*z)) {
            return Err(MimoLoadError("delayed MIMO sampled response left the finite set"));
        }
        Ok(response)
    }

    fn sampled_passivity(&self, points: usize) -> MimoPassivityReport {
        let scale = self.old_direct.iter().chain(&self.new_direct)
            .chain(&self.source.model.d).fold(0.0f64, |m, v| m.max(v.abs()));
        let mut worst = (f64::INFINITY, 0.0);
        let mut violation_freqs = Vec::new();
        for i in 0..points {
            let w = core::f64::consts::PI / self.time_step_s * (i as f64 / (points - 1) as f64);
            let lambda = self.eval_mean_discrete_iw(w).map(|z| hermitian_part(&z, self.ports()).0[0]);
            if let Ok(lambda) = lambda {
                if lambda < worst.0 { worst = (lambda, w); }
                if !lambda.is_finite() || lambda < -1.0e-12 * scale { violation_freqs.push(w); }
            } else { violation_freqs.push(w); }
        }
        MimoPassivityReport { passive: violation_freqs.is_empty(), class: CertificateClass::GridOnly,
            worst, crossings: Vec::new(), violation_freqs }
    }
}

#[cfg(test)]
mod tests;

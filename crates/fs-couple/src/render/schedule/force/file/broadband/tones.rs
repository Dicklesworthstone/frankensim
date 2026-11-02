//! Prescribed continuous physical forces through the existing signed ports.
use super::super::Reader;
use super::super::tone_envelope::{GainEnvelope, ToneEnvelopePlan};
use super::super::tone_frequency::{FrequencySweep, ToneFrequencyPlan};
use super::super::tone_record::PhysicalTone;
use crate::modal_acoustic_time::broadband::{BroadbandModalObserver, PreparedBroadbandTones};
use fs_modal::time::{MAX_HARMONIC_TONES, ModalToneBudget};

pub(super) struct ToneProgram {
    pub sources: Vec<Option<ToneSource>>,
    pub count: usize,
    pub envelope_count: usize,
    pub frequency_sweep_count: usize,
    pub maximum_sweep_force_error_bound_n: f64,
    pub budget: ModalToneBudget,
}

pub(super) struct ToneSource {
    omega: Vec<f64>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    envelopes: Option<Vec<Option<GainEnvelope>>>,
    sweeps: Option<Vec<Option<FrequencySweep>>>,
}

/// One interval-local phase bank. The original performance starts at phase zero;
/// callback partitioning never resets it or samples a tone into held forces.
pub(super) struct ToneDriver {
    prepared: PreparedBroadbandTones,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    next_cosine: Vec<f64>,
    next_sine: Vec<f64>,
    envelope: Option<EnvelopeDriver>,
}

struct EnvelopeDriver {
    curves: Vec<Option<GainEnvelope>>,
    sweeps: Option<Vec<Option<FrequencySweep>>>,
    sample: u64,
    start_cosine: Vec<f64>,
    start_sine: Vec<f64>,
    end_cosine: Vec<f64>,
    end_sine: Vec<f64>,
}

/// Locate a complete optional source section before the receiver section.
pub(super) fn section(text: &str) -> Option<(usize, usize)> {
    let mut offset = 0;
    for (line, raw) in text.split_inclusive('\n').enumerate() {
        if raw.split_ascii_whitespace().next() == Some("tone_limits") {
            return Some((offset, line));
        }
        offset += raw.len();
    }
    None
}

impl ToneProgram {
    pub fn read(
        text: &str,
        first_line: usize,
        columns: &[&[Vec<f64>]],
        sample_rate_hz: u32,
        samples: u64,
        nyquist_guard_fraction: f64,
        projection_terms: &mut usize,
    ) -> Result<Self, String> {
        let mut reader = Reader {
            lines: text.lines(),
            line: first_line,
        };
        let mut row = reader.row("tone_limits").map_err(|e| e.to_string())?;
        let ceiling = ModalToneBudget::realtime_reference();
        let mut budget = ModalToneBudget {
            maximum_tones: row.count(MAX_HARMONIC_TONES).map_err(|e| e.to_string())?,
            maximum_storage_bytes: row
                .count(ceiling.maximum_storage_bytes)
                .map_err(|e| e.to_string())?,
            maximum_preparation_work: row
                .count(ceiling.maximum_preparation_work)
                .map_err(|e| e.to_string())?,
            maximum_step_work: row
                .count(ceiling.maximum_step_work)
                .map_err(|e| e.to_string())?,
        };
        row.finish().map_err(|e| e.to_string())?;
        if budget.maximum_tones < 2 || sample_rate_hz == 0 || columns.is_empty() {
            return Err(
                "continuous tones need a clock, voices and capacity for held DC plus a tone".into(),
            );
        }
        let count: usize = reader.one("tones").map_err(|e| e.to_string())?;
        if count == 0 || count > columns.len() * (budget.maximum_tones - 1) {
            return Err(
                "continuous tone count exceeds the explicit per-voice capacity including DC".into(),
            );
        }
        let parse_storage = count
            .checked_mul(2)
            .and_then(|n| n.checked_add(columns.len() * 4))
            .and_then(|n| n.checked_mul(core::mem::size_of::<(usize, usize, f64, f64, f64)>()))
            .ok_or("continuous tone record size overflow")?;
        if parse_storage > budget.maximum_storage_bytes {
            return Err("continuous tone records exceed tone_limits storage".into());
        }
        // Keep original physical rows intact for sweep-error admission; group
        // only indices so interleaved voices retain their authored tone IDs.
        let mut records = Vec::with_capacity(count);
        let mut grouped = vec![Vec::new(); columns.len()];
        for original_tone in 0..count {
            let tone = PhysicalTone::read(&mut reader, columns, sample_rate_hz, projection_terms)?;
            if grouped[tone.voice].len() + 2 > budget.maximum_tones {
                return Err(
                    "continuous tones plus held DC exceed the voice's explicit capacity".into(),
                );
            }
            grouped[tone.voice].push(original_tone);
            records.push(tone);
        }
        budget.maximum_storage_bytes -= parse_storage;
        let mut envelopes = ToneEnvelopePlan::read_optional(
            &mut reader,
            count,
            samples,
            &mut budget.maximum_storage_bytes,
        )?;
        let envelope_count = envelopes.as_ref().map_or(0, ToneEnvelopePlan::count);
        let mut sweeps = ToneFrequencyPlan::read_optional(
            &mut reader,
            &records,
            samples,
            sample_rate_hz,
            nyquist_guard_fraction,
            envelopes.as_ref(),
            &mut budget.maximum_storage_bytes,
        )?;
        let frequency_sweep_count = sweeps.as_ref().map_or(0, ToneFrequencyPlan::count);
        let maximum_sweep_force_error_bound_n = sweeps
            .as_ref()
            .map_or(0.0, ToneFrequencyPlan::maximum_force_error_bound_n);
        if reader.lines.next().is_some() {
            return Err("unexpected trailing continuous-tone record".into());
        }
        // These file ceilings cover all affected voices together. Reserve the
        // driver's four force arrays and frequency arrays before dividing the
        // remaining mechanics allowance equally between active voices.
        let active = grouped.iter().filter(|rows| !rows.is_empty()).count();
        let mut storage = 0usize;
        let mut envelope_step_work = 0usize;
        for (ports, rows) in columns
            .iter()
            .zip(&grouped)
            .filter(|(_, rows)| !rows.is_empty())
        {
            let tones = rows.len() + 1;
            let modes = ports.first().map_or(0, Vec::len);
            let mut bytes = modes
                .checked_mul(tones)
                .and_then(|n| n.checked_mul(32))
                .and_then(|n| n.checked_add(tones * 8))
                .ok_or("tone buffer size overflow")?;
            let swept = rows
                .iter()
                .filter(|&&tone| sweeps.as_ref().is_some_and(|plan| plan.contains(tone)))
                .count();
            if swept != 0
                || envelopes
                    .as_ref()
                    .is_some_and(|plan| rows.iter().any(|&tone| plan.contains(tone)))
            {
                // Four scaled endpoint arrays are additional to the unchanged
                // rotating carrier/next-carrier arrays. Curve metadata is moved
                // into voices while the original parse index still exists.
                let extra = modes
                    .checked_mul(tones)
                    .and_then(|n| n.checked_mul(32))
                    .and_then(|n| {
                        n.checked_add(tones * core::mem::size_of::<Option<GainEnvelope>>())
                    })
                    .and_then(|n| n.checked_add(core::mem::size_of::<EnvelopeDriver>()))
                    .ok_or("envelope endpoint buffer size overflow")?;
                bytes = bytes
                    .checked_add(extra)
                    .ok_or("envelope endpoint buffer size overflow")?;
                if swept != 0 {
                    bytes = bytes
                        .checked_add(tones * core::mem::size_of::<Option<FrequencySweep>>())
                        .ok_or("frequency sweep state size overflow")?;
                    // One modulation/phase plan per swept tone, then bounded
                    // matrix-vector arithmetic per modal projection. Trig is
                    // outside the mode loop and no runtime allocation occurs.
                    envelope_step_work = envelope_step_work
                        .checked_add(
                            modes
                                .checked_mul(16)
                                .and_then(|n| n.checked_add(128))
                                .and_then(|n| n.checked_mul(swept))
                                .ok_or("frequency sweep work overflow")?,
                        )
                        .ok_or("frequency sweep work overflow")?;
                }
                envelope_step_work = envelope_step_work
                    .checked_add(
                        modes
                            .checked_mul(tones)
                            .and_then(|n| n.checked_mul(4))
                            .and_then(|n| n.checked_add(12 * tones))
                            .ok_or("envelope work overflow")?,
                    )
                    .ok_or("envelope work overflow")?;
            }
            storage = storage
                .checked_add(bytes)
                .ok_or("tone buffer size overflow")?;
        }
        budget.maximum_storage_bytes = budget
            .maximum_storage_bytes
            .checked_sub(storage)
            .ok_or("continuous force buffers exceed tone_limits storage")?
            / active;
        budget.maximum_preparation_work /= active;
        budget.maximum_step_work = budget
            .maximum_step_work
            .checked_sub(envelope_step_work)
            .ok_or("physical tone modulation exceeds tone_limits step work")?
            / active;
        let mut sources = Vec::with_capacity(columns.len());
        for (ports, rows) in columns.iter().zip(grouped) {
            if rows.is_empty() {
                sources.push(None);
                continue;
            }
            let tones = rows.len() + 1;
            let modes = ports[0].len();
            let mut omega = Vec::with_capacity(tones);
            omega.push(0.0); // Existing held schedule is the work-conjugate DC input.
            let mut cosine = vec![0.0; modes * tones];
            let mut sine = cosine.clone();
            let swept = sweeps
                .as_ref()
                .is_some_and(|plan| rows.iter().any(|&tone| plan.contains(tone)));
            let affine = swept
                || envelopes
                    .as_ref()
                    .is_some_and(|plan| rows.iter().any(|&tone| plan.contains(tone)));
            let mut curves = affine.then(|| Vec::with_capacity(tones));
            let mut frequency_curves = swept.then(|| Vec::with_capacity(tones));
            if let Some(curves) = &mut curves {
                curves.push(None);
            } // Held DC is not enveloped.
            if let Some(curves) = &mut frequency_curves {
                curves.push(None);
            } // Held DC is not frequency-modulated.
            for (tone, original_tone) in rows.into_iter().enumerate() {
                let record = &records[original_tone];
                if let Some(curves) = &mut curves {
                    curves.push(envelopes.as_mut().and_then(|plan| plan.take(original_tone)));
                }
                if let Some(curves) = &mut frequency_curves {
                    curves.push(sweeps.as_mut().and_then(|plan| plan.take(original_tone)));
                }
                omega.push(record.omega);
                for (mode, &shape) in ports[record.port].iter().enumerate() {
                    let index = mode * tones + tone + 1;
                    cosine[index] = shape * record.cosine;
                    sine[index] = shape * record.sine;
                    if !cosine[index].is_finite() || !sine[index].is_finite() {
                        return Err("projected continuous force is not representable".into());
                    }
                }
            }
            sources.push(Some(ToneSource {
                omega,
                cosine,
                sine,
                envelopes: curves,
                sweeps: frequency_curves,
            }));
        }
        Ok(Self {
            sources,
            count,
            envelope_count,
            frequency_sweep_count,
            maximum_sweep_force_error_bound_n,
            budget,
        })
    }
}

impl ToneSource {
    pub fn prepare(
        self,
        observer: &BroadbandModalObserver,
        budget: ModalToneBudget,
    ) -> Result<ToneDriver, String> {
        let prepared = if self.envelopes.is_some() {
            observer.prepare_affine_tones(&self.omega, budget)
        } else {
            observer.prepare_tones(&self.omega, budget)
        }
        .map_err(|e| e.to_string())?;
        let envelope = self.envelopes.map(|curves| EnvelopeDriver {
            curves,
            sweeps: self.sweeps,
            sample: 0,
            start_cosine: vec![0.0; self.cosine.len()],
            start_sine: vec![0.0; self.cosine.len()],
            end_cosine: vec![0.0; self.cosine.len()],
            end_sine: vec![0.0; self.cosine.len()],
        });
        Ok(ToneDriver {
            prepared,
            next_cosine: vec![0.0; self.cosine.len()],
            next_sine: vec![0.0; self.sine.len()],
            cosine: self.cosine,
            sine: self.sine,
            envelope,
        })
    }
}

impl ToneDriver {
    pub fn step(
        &mut self,
        observer: &mut BroadbandModalObserver,
        held: &[f64],
    ) -> Result<f64, String> {
        let tones = self.prepared.tone_count();
        if held.len() * tones != self.cosine.len() {
            return Err(
                "held and continuous forces must use the original complete mode basis".into(),
            );
        }
        let pressure = if let Some(envelope) = &mut self.envelope {
            let next_sample = envelope
                .sample
                .checked_add(1)
                .ok_or("envelope sample clock overflow")?;
            for (tone, curve) in envelope.curves.iter().enumerate() {
                let gain = curve
                    .as_ref()
                    .map_or(Ok([1.0, 1.0]), |curve| curve.endpoints(envelope.sample))?;
                // The complete modulation plan is evaluated once for the
                // original physical tone, independently of modal dimension.
                let sweep = envelope
                    .sweeps
                    .as_ref()
                    .and_then(|curves| curves[tone].as_ref())
                    .map(|curve| curve.step(envelope.sample, gain))
                    .transpose()?;
                for (mode, &held_force) in held.iter().enumerate() {
                    let index = mode * tones + tone;
                    let (c, s) = if tone == 0 {
                        (held_force, 0.0)
                    } else {
                        (self.cosine[index], self.sine[index])
                    };
                    let (endpoints, next) = if let Some(step) = sweep {
                        (
                            step.endpoints(c, s).map_err(|e| e.to_string())?,
                            step.advance_force_phase(c, s).map_err(|e| e.to_string())?,
                        )
                    } else {
                        let endpoints = [gain[0] * c, gain[0] * s, gain[1] * c, gain[1] * s];
                        if endpoints.iter().any(|x| !x.is_finite()) {
                            return Err("enveloped physical force is not representable".into());
                        }
                        (
                            endpoints,
                            self.prepared
                                .advance_force_phase(tone, c, s)
                                .map_err(|e| e.to_string())?,
                        )
                    };
                    envelope.start_cosine[index] = endpoints[0];
                    envelope.start_sine[index] = endpoints[1];
                    envelope.end_cosine[index] = endpoints[2];
                    envelope.end_sine[index] = endpoints[3];
                    self.next_cosine[index] = next.0;
                    self.next_sine[index] = next.1;
                }
            }
            let pressure = observer
                .step_affine_tones(
                    &self.prepared,
                    &envelope.start_cosine,
                    &envelope.start_sine,
                    &envelope.end_cosine,
                    &envelope.end_sine,
                )
                .map_err(|e| e.to_string())?;
            for curve in envelope.curves.iter_mut().flatten() {
                curve.commit(next_sample);
            }
            if let Some(curves) = &mut envelope.sweeps {
                for curve in curves.iter_mut().flatten() {
                    curve.commit(next_sample);
                }
            }
            envelope.sample = next_sample;
            pressure
        } else {
            // Preserve the existing constant-tone numerical path exactly.
            for (mode, &force) in held.iter().enumerate() {
                self.cosine[mode * tones] = force;
                self.next_cosine[mode * tones] = force;
                for tone in 1..tones {
                    let index = mode * tones + tone;
                    let (c, s) = self
                        .prepared
                        .advance_force_phase(tone, self.cosine[index], self.sine[index])
                        .map_err(|e| e.to_string())?;
                    self.next_cosine[index] = c;
                    self.next_sine[index] = s;
                }
            }
            observer
                .step_tones(&self.prepared, &self.cosine, &self.sine)
                .map_err(|e| e.to_string())?
        };
        // Only an accepted physical and acoustic step advances source phase.
        core::mem::swap(&mut self.cosine, &mut self.next_cosine);
        core::mem::swap(&mut self.sine, &mut self.next_sine);
        Ok(pressure)
    }
}

#[cfg(test)]
#[path = "tones/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "tones/envelope_tests.rs"]
mod envelope_tests;

#[cfg(test)]
#[path = "tones/frequency_tests.rs"]
mod frequency_tests;

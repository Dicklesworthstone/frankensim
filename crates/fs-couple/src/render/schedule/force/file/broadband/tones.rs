//! Prescribed continuous physical forces through the existing signed ports.
use super::super::Reader;
use super::super::tone_envelope::{GainEnvelope, ToneEnvelopePlan};
use super::super::tone_record::PhysicalTone;
use crate::modal_acoustic_time::broadband::{BroadbandModalObserver, PreparedBroadbandTones};
use fs_modal::time::{MAX_HARMONIC_TONES, ModalToneBudget};

pub(super) struct ToneProgram {
    pub sources: Vec<Option<ToneSource>>,
    pub count: usize,
    pub envelope_count: usize,
    pub budget: ModalToneBudget,
}

pub(super) struct ToneSource {
    omega: Vec<f64>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    envelopes: Option<Vec<Option<GainEnvelope>>>,
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
        let mut records = vec![Vec::new(); columns.len()];
        for original_tone in 0..count {
            let PhysicalTone {
                voice,
                port,
                omega,
                cosine,
                sine,
            } = PhysicalTone::read(&mut reader, columns, sample_rate_hz, projection_terms)?;
            if records[voice].len() + 2 > budget.maximum_tones {
                return Err(
                    "continuous tones plus held DC exceed the voice's explicit capacity".into(),
                );
            }
            records[voice].push((original_tone, port, omega, cosine, sine));
        }
        budget.maximum_storage_bytes -= parse_storage;
        let mut envelopes = ToneEnvelopePlan::read_optional(
            &mut reader,
            count,
            samples,
            &mut budget.maximum_storage_bytes,
        )?;
        let envelope_count = envelopes.as_ref().map_or(0, ToneEnvelopePlan::count);
        if reader.lines.next().is_some() {
            return Err("unexpected trailing continuous-tone record".into());
        }
        // These file ceilings cover all affected voices together. Reserve the
        // driver's four force arrays and frequency arrays before dividing the
        // remaining mechanics allowance equally between active voices.
        let active = records.iter().filter(|rows| !rows.is_empty()).count();
        let mut storage = 0usize;
        let mut envelope_step_work = 0usize;
        for (ports, rows) in columns
            .iter()
            .zip(&records)
            .filter(|(_, rows)| !rows.is_empty())
        {
            let tones = rows.len() + 1;
            let modes = ports.first().map_or(0, Vec::len);
            let mut bytes = modes
                .checked_mul(tones)
                .and_then(|n| n.checked_mul(32))
                .and_then(|n| n.checked_add(tones * 8))
                .ok_or("tone buffer size overflow")?;
            if envelopes
                .as_ref()
                .is_some_and(|plan| rows.iter().any(|row| plan.contains(row.0)))
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
            .ok_or("physical gain envelopes exceed tone_limits step work")?
            / active;
        let mut sources = Vec::with_capacity(columns.len());
        for (ports, rows) in columns.iter().zip(records) {
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
            let affine = envelopes
                .as_ref()
                .is_some_and(|plan| rows.iter().any(|row| plan.contains(row.0)));
            let mut curves = affine.then(|| Vec::with_capacity(tones));
            if let Some(curves) = &mut curves {
                curves.push(None);
            } // Held DC is not enveloped.
            for (tone, (original_tone, port, frequency, c, s)) in rows.into_iter().enumerate() {
                if let Some(curves) = &mut curves {
                    curves.push(envelopes.as_mut().and_then(|plan| plan.take(original_tone)));
                }
                omega.push(frequency);
                for (mode, &shape) in ports[port].iter().enumerate() {
                    let index = mode * tones + tone + 1;
                    cosine[index] = shape * c;
                    sine[index] = shape * s;
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
            }));
        }
        Ok(Self {
            sources,
            count,
            envelope_count,
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
        let pressure = if let Some(envelope) = &mut self.envelope {
            let next_sample = envelope
                .sample
                .checked_add(1)
                .ok_or("envelope sample clock overflow")?;
            for (tone, curve) in envelope.curves.iter().enumerate() {
                let [start, end] = curve
                    .as_ref()
                    .map_or(Ok([1.0, 1.0]), |curve| curve.endpoints(envelope.sample))?;
                for mode in 0..held.len() {
                    let index = mode * tones + tone;
                    envelope.start_cosine[index] = start * self.cosine[index];
                    envelope.start_sine[index] = start * self.sine[index];
                    envelope.end_cosine[index] = end * self.cosine[index];
                    envelope.end_sine[index] = end * self.sine[index];
                    if [
                        envelope.start_cosine[index],
                        envelope.start_sine[index],
                        envelope.end_cosine[index],
                        envelope.end_sine[index],
                    ]
                    .iter()
                    .any(|x| !x.is_finite())
                    {
                        return Err("enveloped physical force is not representable".into());
                    }
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
            envelope.sample = next_sample;
            pressure
        } else {
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

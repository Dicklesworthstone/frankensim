//! Version-2 physical tone forcing of the complete coupled mechanical owner.
use super::{
    ForceInitialization, MAX_MODAL_PERFORMANCE_BYTES, MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN,
    MODAL_COUPLED_PERFORMANCE_SCHEMA, ModalPerformanceInfo, ParsedPerformance, Reader,
    tone_envelope::ToneEnvelopePlan, tone_record::PhysicalTone,
};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{
    RenderError,
    schedule::force::{
        coupled::{
            CoupledModalSystem,
            tones::{CoupledModalTones, CoupledToneBudget},
        },
        program::{PhysicalForceInput, ProjectedForceEvent, ProjectedForceProgram},
    },
};
use fs_exec::CancelGate;

#[path = "coupled_tones/broadband.rs"]
mod broadband;

const MAX_TONES: usize = 32;
const MAX_STORAGE_BYTES: usize = 128 * 1024 * 1024;
const MAX_PREPARATION_WORK: u64 = 100_000_000_000;

/// A v2 assembly with continuous physical port forces and its original held
/// controls. Every spring, dashpot, state and work/heat channel stays in the
/// complete coupled owner. Observation uses the supplied real instantaneous
/// pressure/velocity map or an explicitly sampled complete-network causal
/// receiver; neither implies instrument validity or an acoustic reaction load.
pub struct CoupledTonePerformance {
    info: ModalPerformanceInfo,
    source: broadband::Source,
    held: Vec<f64>,
    events: Vec<ProjectedForceEvent>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    next_cosine: Vec<f64>,
    next_sine: Vec<f64>,
    envelopes: Option<ToneEnvelopePlan>,
    force_cosine: Vec<f64>,
    force_sine: Vec<f64>,
    force_cosine_end: Vec<f64>,
    force_sine_end: Vec<f64>,
    next_event: usize,
    completed: u64,
    max_block: usize,
    poisoned: bool,
}

impl CoupledTonePerformance {
    /// Detect the explicit extension without changing the held-only file path.
    pub fn has_section(bytes: &[u8]) -> bool {
        std::str::from_utf8(bytes).is_ok_and(|text| section(text).is_some())
    }

    /// Parse the existing v2 mechanical file, followed by:
    ///
    /// ```text
    /// coupled_tone_limits MAX_TONES MAX_DIMENSION STORAGE_BYTES MULTIPLY_TERMS
    /// tones COUNT
    /// tone VOICE PORT OMEGA_RAD_S COSINE_N SINE_N
    /// ```
    ///
    /// Tones are physical forces continuous from t=0 over the declared window.
    /// Optional `envelopes`/`envelope`/`gain` records give piecewise-affine
    /// physical force gains without resetting carrier phase or resonant state.
    /// MAX_TONES includes the held DC column and is at most 32. The dense cold
    /// dimension is at most 256; total tone storage is at most 128 MiB and
    /// preparation work at most 100 billion multiply-add terms. Driver phase
    /// buffers, records and frequencies are reserved before passing the remaining
    /// byte allowance to the existing dense owner. Original file/model/force
    /// limits remain binding. The existing signed port projection is charged to
    /// compile_limits, including both quadratures.
    ///
    /// Static preload solves the complete original spring equilibrium under
    /// pre-window held loads. It does not assume tones ran before the window.
    /// An optional `coupled_broadband` suffix supplies one causal receiver over
    /// all flattened modal velocities; see `examples/COUPLED_CONTINUOUS_FORCES.md`.
    /// Without it, imaginary constant pressure maps refuse. Contact/friction
    /// schemas are never replaced by this admitted linear trajectory.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, String> {
        if bytes.len() > MAX_MODAL_PERFORMANCE_BYTES {
            return Err("coupled tone input exceeds 4 MiB".into());
        }
        let text = std::str::from_utf8(bytes).map_err(|_| "coupled tone input must be UTF-8")?;
        let (offset, line) = section(text).ok_or("explicit coupled_tone_limits required")?;
        let mut parsed = ParsedPerformance::from_bytes(&bytes[..offset], max_block)
            .map_err(|e| e.to_string())?;
        if parsed.info.schema != MODAL_COUPLED_PERFORMANCE_SCHEMA
            || parsed.contact.is_some()
            || parsed.multiple.is_some()
            || parsed.frictions.is_some()
        {
            return Err("continuous coupled tones require the original linear v2 mechanics".into());
        }
        let (connections, coupling) = parsed.coupled.take().ok_or("v2 connections required")?;
        let mut reader = Reader {
            lines: text[offset..].lines(),
            line,
        };
        let mut row = reader
            .row("coupled_tone_limits")
            .map_err(|e| e.to_string())?;
        let maximum_tones = row.count(MAX_TONES).map_err(|e| e.to_string())?;
        let mut budget = CoupledToneBudget {
            max_dimension: row
                .count(fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM)
                .map_err(|e| e.to_string())?,
            max_dense_bytes: row.count(MAX_STORAGE_BYTES).map_err(|e| e.to_string())?,
            max_multiply_terms: row.parse().map_err(|e| e.to_string())?,
        };
        row.finish().map_err(|e| e.to_string())?;
        let count: usize = reader.one("tones").map_err(|e| e.to_string())?;
        if maximum_tones < 2
            || count == 0
            || count >= maximum_tones
            || budget.max_dimension == 0
            || budget.max_multiply_terms == 0
            || budget.max_multiply_terms > MAX_PREPARATION_WORK
        {
            return Err(
                "coupled tones exceed their explicit count, dimension or preparation budget".into(),
            );
        }
        let tones = count + 1;
        let coefficients = parsed
            .info
            .modes
            .checked_mul(tones)
            .ok_or("tone extent overflow")?;
        // Conservative simultaneous driver/parse payload, before any source
        // array allocation. Original model and force-compiler storage keeps its
        // existing separate bounds; allocator metadata and caller bytes excluded.
        let driver_bytes = coefficients
            .checked_mul(4)
            .and_then(|n| n.checked_add(2 * parsed.info.modes + tones + 3 * parsed.info.voices))
            .and_then(|n| n.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(count * core::mem::size_of::<PhysicalTone>()))
            .ok_or("tone storage overflow")?;
        budget.max_dense_bytes = budget
            .max_dense_bytes
            .checked_sub(driver_bytes)
            .ok_or("coupled tone phase buffers exceed declared storage")?;
        let columns: Vec<_> = parsed
            .voices
            .iter()
            .map(|v| v.port_shapes.as_slice())
            .collect();
        let mut records = Vec::with_capacity(count);
        for _ in 0..count {
            records.push(PhysicalTone::read(
                &mut reader,
                &columns,
                parsed.info.sample_rate_hz,
                &mut parsed.force_config.max_projection_terms,
            )?);
        }
        let envelopes = ToneEnvelopePlan::read_optional(
            &mut reader,
            count,
            parsed.info.samples,
            &mut budget.max_dense_bytes,
        )?;
        let envelope_coefficients = if envelopes.is_some() { coefficients } else { 0 };
        let envelope_bytes = envelope_coefficients
            .checked_mul(4 * core::mem::size_of::<f64>())
            .ok_or("coupled envelope buffer size overflow")?;
        budget.max_dense_bytes = budget
            .max_dense_bytes
            .checked_sub(envelope_bytes)
            .ok_or("coupled envelope force buffers exceed source storage")?;
        let receiver = broadband::ReceiverSpec::read_optional(&mut reader, parsed.info.modes)?;
        let mut offsets = Vec::with_capacity(parsed.voices.len());
        let mut modes = 0;
        for voice in &parsed.voices {
            offsets.push(modes);
            modes += voice.model.modes().len();
        }
        let mut omega = Vec::with_capacity(tones);
        omega.push(0.0);
        let mut cosine = vec![0.0; coefficients];
        let mut sine = vec![0.0; coefficients];
        for (tone, record) in records.into_iter().enumerate() {
            omega.push(record.omega);
            for (mode, &shape) in columns[record.voice][record.port].iter().enumerate() {
                let index = (offsets[record.voice] + mode) * tones + tone + 1;
                cosine[index] = shape * record.cosine;
                sine[index] = shape * record.sine;
                if !cosine[index].is_finite() || !sine[index].is_finite() {
                    return Err("projected coupled tone is not representable".into());
                }
            }
        }
        let preload = parsed
            .voices
            .iter()
            .any(|v| v.initialization == ForceInitialization::StaticPreload);
        if preload
            && parsed
                .voices
                .iter()
                .any(|v| v.initialization != ForceInitialization::StaticPreload)
        {
            return Err("coupled preload must cover every component".into());
        }
        let mut models = Vec::with_capacity(parsed.voices.len());
        let mut inputs = Vec::with_capacity(parsed.voices.len());
        for voice in parsed.voices {
            inputs.push(
                PhysicalForceInput::new(voice.port_shapes, voice.initial_force_n)
                    .map_err(|e| e.to_string())?,
            );
            models.push(if receiver.is_some() {
                // The sampled receiver owns observation. Clearing only the
                // legacy map permits its complex narrow-band metadata without
                // altering any physical coordinate, preload or heat channel.
                voice.model.without_pressure_transfers()
            } else {
                voice.model
            });
        }
        let (initial, mut events) =
            ProjectedForceProgram::compile(inputs, parsed.events, parsed.force_config)
                .map_err(|e| e.to_string())?
                .into_parts();
        let held: Vec<_> = initial.into_iter().flatten().collect();
        for event in &mut events {
            event.mode += offsets[event.voice];
            event.voice = 0;
        }
        for (mode, &force) in held.iter().enumerate() {
            cosine[mode * tones] = force;
        }
        let gate = CancelGate::new();
        let mut network = CoupledModalSystem::new(models, connections, coupling, &gate)
            .map_err(|e| e.to_string())?;
        if preload {
            network
                .initialize_static_equilibrium(&held, &gate)
                .map_err(|e| e.to_string())?;
        }
        let source = if envelopes.is_some() {
            CoupledModalTones::new_affine(network, &omega, budget, &gate)
        } else {
            CoupledModalTones::new(network, &omega, budget, &gate)
        }
        .map_err(|e| e.to_string())?;
        let source = broadband::Source::new(source, receiver, parsed.info.sample_rate_hz, &gate)?;
        let mut info = parsed.info;
        info.input_hash = fs_blake3::hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, bytes);
        Ok(Self {
            info,
            source,
            held,
            events,
            next_cosine: vec![0.0; coefficients],
            next_sine: vec![0.0; coefficients],
            cosine,
            sine,
            envelopes,
            force_cosine: vec![0.0; envelope_coefficients],
            force_sine: vec![0.0; envelope_coefficients],
            force_cosine_end: vec![0.0; envelope_coefficients],
            force_sine_end: vec![0.0; envelope_coefficients],
            next_event: 0,
            completed: 0,
            max_block,
            poisoned: false,
        })
    }

    pub const fn info(&self) -> ModalPerformanceInfo {
        self.info
    }
    /// Original accepted mechanics, clock, and independent work/heat witness.
    pub fn network(&self) -> &CoupledModalSystem {
        self.source.network()
    }
    pub fn compiled_controls(&self) -> usize {
        self.events.len()
    }
    /// Authored physical tone records, excluding the held DC coordinate.
    pub fn continuous_force_tones(&self) -> usize {
        self.source.tone_count() - 1
    }
    /// Explicit physical force gain curves, indexed by original tone records.
    pub fn envelope_count(&self) -> usize {
        self.envelopes.as_ref().map_or(0, ToneEnvelopePlan::count)
    }
    /// The explicit sampled receiver's fit, delay and held-out evidence, when
    /// present. Producer samples and this interpolation remain EstimateOnly.
    pub fn broadband_report(
        &self,
    ) -> Option<&crate::modal_acoustic_time::broadband::BroadbandModalReport> {
        self.source.report()
    }

    fn step(&mut self) -> Result<f64, RenderError> {
        let tones = self.source.tone_count();
        for (mode, &held) in self.held.iter().enumerate() {
            self.next_cosine[mode * tones] = held;
            for tone in 1..tones {
                let index = mode * tones + tone;
                let (c, s) = self
                    .source
                    .advance_force_phase(tone, self.cosine[index], self.sine[index])
                    .map_err(RenderError::Coupled)?;
                self.next_cosine[index] = c;
                self.next_sine[index] = s;
            }
        }
        let mut next_event = self.next_event;
        while next_event < self.events.len() && self.events[next_event].sample == self.completed {
            let event = self.events[next_event];
            self.next_cosine[event.mode * tones] = event.force_n_per_sqrt_kg;
            next_event += 1;
        }
        let candidate = if self.envelopes.is_some() {
            // Scale independent carrier phases, including zero/re-attack.
            // Stage both endpoints before changing any accepted held input.
            self.stage_envelope_forces()?;
            self.source.step_affine(
                &self.force_cosine,
                &self.force_sine,
                &self.force_cosine_end,
                &self.force_sine_end,
            )
        } else {
            for mode in 0..self.held.len() {
                self.cosine[mode * tones] = self.next_cosine[mode * tones];
            }
            self.source.step(&self.cosine, &self.sine)
        };
        let pressure = match candidate {
            Ok(pressure) => pressure,
            Err(error) => {
                for (mode, &held) in self.held.iter().enumerate() {
                    self.cosine[mode * tones] = held;
                }
                return Err(error);
            }
        };
        for (mode, held) in self.held.iter_mut().enumerate() {
            *held = self.next_cosine[mode * tones];
        }
        core::mem::swap(&mut self.cosine, &mut self.next_cosine);
        core::mem::swap(&mut self.sine, &mut self.next_sine);
        self.next_event = next_event;
        self.completed += 1;
        if let Some(envelopes) = &mut self.envelopes {
            envelopes.commit(self.completed);
        }
        Ok(pressure)
    }

    fn stage_envelope_forces(&mut self) -> Result<(), RenderError> {
        let envelopes = self.envelopes.as_ref().ok_or_else(refused)?;
        let tones = self.source.tone_count();
        // Query each original tone only once, then apply its gain to every
        // signed modal projection. The held DC column keeps its event value.
        for mode in 0..self.held.len() {
            let dc = mode * tones;
            self.force_cosine[dc] = self.next_cosine[dc];
            self.force_cosine_end[dc] = self.next_cosine[dc];
            self.force_sine[dc] = 0.0;
            self.force_sine_end[dc] = 0.0;
        }
        for tone in 1..tones {
            let [start, end] = envelopes.endpoints(tone - 1, self.completed).map_err(|_| {
                RenderError::Control {
                    what: "physical tone gain cursor or interpolation refused",
                }
            })?;
            for mode in 0..self.held.len() {
                let i = mode * tones + tone;
                let values = [
                    start * self.cosine[i],
                    start * self.sine[i],
                    end * self.cosine[i],
                    end * self.sine[i],
                ];
                if values.iter().any(|x| !x.is_finite()) {
                    return Err(RenderError::Control {
                        what: "physical tone gain endpoints overflow",
                    });
                }
                [
                    self.force_cosine[i],
                    self.force_sine[i],
                    self.force_cosine_end[i],
                    self.force_sine_end[i],
                ] = values;
            }
        }
        Ok(())
    }
}

fn section(text: &str) -> Option<(usize, usize)> {
    let mut offset = 0;
    for (line, raw) in text.split_inclusive('\n').enumerate() {
        if raw.split_ascii_whitespace().next() == Some("coupled_tone_limits") {
            return Some((offset, line));
        }
        offset += raw.len();
    }
    None
}
fn refused() -> RenderError {
    RenderError::Control {
        what: "coupled tone performance exceeded its admitted clock, callback or window",
    }
}
impl PressureRenderer for CoupledTonePerformance {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }
    fn max_block_len(&self) -> usize {
        self.max_block
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned || rate != self.info.sample_rate_hz {
            return Err(refused());
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned || samples > self.info.samples - self.completed {
            return Err(refused());
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() || output.len() > self.max_block {
            return Err(refused());
        }
        for slot in output {
            match self.step() {
                Ok(pressure) => *slot = pressure,
                Err(error) => {
                    self.poisoned = true;
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "coupled_tones/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "coupled_tones/broadband_tests.rs"]
mod broadband_tests;

#[cfg(test)]
#[path = "coupled_tones/envelope_tests.rs"]
mod envelope_tests;

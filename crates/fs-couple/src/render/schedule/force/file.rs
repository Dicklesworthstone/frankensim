//! Authored reduced structural images -> physical force schedules -> live renderer.
//!
//! `frankensim-modal-performance-v1` is a line-oriented input, NOT a certificate
//! or a CAD/material-card importer. Modes and port columns must already share
//! one mass-normalized basis. Observer transfers are the existing narrow-band
//! complex transfers under exp(-i omega t). Nothing converts note numbers to
//! frequencies or invents acoustic gains. All numbers are in the SI units below.
//!
//! Record order (ASCII whitespace separates fields; no ignored fields/records):
//! ```text
//! frankensim-modal-performance-v1
//! sample_rate_hz RATE
//! samples COUNT
//! full_scale_pa PRESSURE
//! limits NYQUIST_FRACTION MAX_Q MAX_V MAX_ENERGY MAX_PRESSURE
//! compile_limits MAX_CONTROLS MAX_PROJECTION_TERMS
//! voices COUNT
//! voice retain-state|static-preload MODE_COUNT PORT_COUNT
//! mode OMEGA_RAD_S ZETA H_RE H_IM Q V
//! ... one mode row per mode ...
//! port INITIAL_FORCE_N SHAPE_0 ... SHAPE_N_MINUS_1
//! ... one port row per port, then remaining voices ...
//! events COUNT
//! force SAMPLE VOICE_INDEX PORT_INDEX FORCE_N
//! ... one force row per event ...
//! ```
//! Q is m sqrt(kg), V is m sqrt(kg)/s, H is Pa s/(m sqrt(kg)), and
//! port shapes are 1/sqrt(kg). `limits` applies independently to every voice.
//! `static-preload` requires zero Q/V fields: the actual initial state is the
//! declared force equilibrium, not a silently discarded supplied vibration.
//! Every event must occur inside [0, samples); no authored event is dropped.
//! Input byte identity, rather than formatting-insensitive semantic identity,
//! is retained so a sidecar can bind the exact supplied artifact.

use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_blake3::{ContentHash, hash_domain};
use fs_math::c64::C64;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::RenderError;
use super::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use super::super::ScheduledRenderer;

/// Schema token at the beginning of each file.
pub const MODAL_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v1";
/// Byte-read limit to apply BEFORE allocating or decoding an input file.
pub const MAX_MODAL_PERFORMANCE_BYTES: usize = 4 * 1024 * 1024;
/// A bounded offline performance; 600 seconds at 48 kHz.
pub const MAX_MODAL_PERFORMANCE_SAMPLES: u64 = 28_800_000;
/// Domain for the exact input bytes, not a physical-validity assertion.
pub const MODAL_PERFORMANCE_HASH_DOMAIN: &str = "org.frankensim.fs-couple.modal-performance-input.v1";
const MAX_VOICES: usize = 64;
const MAX_MODES: usize = 4096;
const MAX_PORT_WEIGHTS: usize = 65_536;
const MAX_EVENTS: usize = 65_536;
const MAX_CONTROLS: usize = 262_144;
const MAX_PROJECTION_TERMS: usize = 16_777_216;

/// Immutable description of the admitted input and its actual render clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalPerformanceInfo {
    /// Audio samples per second; not inferred from a note or output extension.
    pub sample_rate_hz: u32,
    /// Exact number of output samples, including a final short callback.
    pub samples: u64,
    /// Declared pascals mapped to positive full-scale PCM.
    pub full_scale_pa: f64,
    /// Input bytes under [`MODAL_PERFORMANCE_HASH_DOMAIN`].
    pub input_hash: ContentHash,
    /// Independently simulated voices summed in file order.
    pub voices: usize,
    /// Total number of retained modes across all voices.
    pub modes: usize,
    /// Number of authored physical-force assignments (not expanded modal controls).
    pub force_events: usize,
}

/// A completely decoded and admitted performance, ready for block rendering.
/// The runtime owns the data: no source file is reopened during the render.
pub struct ModalPerformance {
    info: ModalPerformanceInfo,
    renderer: ScheduledRenderer,
}

/// Syntax, resource admission or the existing physics owner's refusal.
#[derive(Debug)]
pub enum ModalPerformanceError {
    /// First offending one-based line; empty/truncated files name the expected line.
    Input {
        /// First offending or missing one-based line.
        line: usize,
        /// Violated syntax or resource rule.
        what: &'static str,
    },
    /// The same model, preload, projection or scheduling refusal as the native API.
    Render(RenderError),
}
impl core::fmt::Display for ModalPerformanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "modal performance line {line}: {what}"),
            Self::Render(error) => write!(f, "modal performance: {error}"),
        }
    }
}
impl std::error::Error for ModalPerformanceError {}
impl From<RenderError> for ModalPerformanceError {
    fn from(value: RenderError) -> Self { Self::Render(value) }
}

impl ModalPerformance {
    /// Decode, bind and compile the whole performance before returning a runtime.
    ///
    /// Host admission caps are 64 voices, 4096 TOTAL modes, 65536 TOTAL port
    /// coefficients/events, 262144 expanded controls and 16777216 projection
    /// terms. Input-authored compilation budgets may be smaller, never enlarged.
    /// Rates are in 1..=192000 Hz, samples in 1..=28800000, and max_block in
    /// 1..=65536. Counts are checked BEFORE count-derived allocations.
    ///
    /// No random source is used. No authority is inferred from imported numbers:
    /// this advances the supplied reduced model, not a validated instrument.
    ///
    /// # Errors
    /// Malformed, oversized or incomplete input; invalid physical data; and any
    /// existing modal-state, preload, projection or schedule admission refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, ModalPerformanceError> {
        if bytes.len() > MAX_MODAL_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(input(1, "input exceeds 4 MiB or callback size is outside 1..=65536"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| input(1, "input must be UTF-8"))?;
        let mut reader = Reader { lines: text.lines(), line: 0 };
        reader.row(MODAL_PERFORMANCE_SCHEMA)?.finish()?;
        let sample_rate_hz: u32 = reader.one("sample_rate_hz")?;
        if !(1..=192_000).contains(&sample_rate_hz) {
            return Err(input(reader.line, "sample_rate_hz must be in 1..=192000"));
        }
        let samples: u64 = reader.one("samples")?;
        if samples == 0 || samples > MAX_MODAL_PERFORMANCE_SAMPLES {
            return Err(input(reader.line, "samples must be in 1..=28800000"));
        }
        let full_scale_pa: f64 = reader.one("full_scale_pa")?;
        if !full_scale_pa.is_finite() || full_scale_pa <= 0.0 {
            return Err(input(reader.line, "full_scale_pa must be finite and positive"));
        }
        let mut row = reader.row("limits")?;
        let budget = ModalAcousticTimeBudget {
            nyquist_guard_fraction: row.scalar()?,
            maximum_abs_displacement_m_sqrt_kg: row.scalar()?,
            maximum_abs_velocity_m_sqrt_kg_per_s: row.scalar()?,
            maximum_total_energy_j: row.scalar()?,
            maximum_abs_pressure_pa: row.scalar()?,
        };
        row.finish()?;
        let mut row = reader.row("compile_limits")?;
        let max_controls = row.count(MAX_CONTROLS)?;
        let max_projection_terms = row.count(MAX_PROJECTION_TERMS)?;
        row.finish()?;
        let voice_count: usize = reader.one("voices")?;
        if voice_count == 0 || voice_count > MAX_VOICES {
            return Err(input(reader.line, "voices must be in 1..=64"));
        }
        let mut voices = Vec::with_capacity(voice_count);
        let mut total_modes = 0;
        let mut total_weights = 0;
        for _ in 0..voice_count {
            let mut row = reader.row("voice")?;
            let initialization = match row.word()? {
                "retain-state" => ForceInitialization::RetainState,
                "static-preload" => ForceInitialization::StaticPreload,
                _ => return Err(input(row.line, "expected retain-state or static-preload")),
            };
            let modes = row.count(MAX_MODES - total_modes)?;
            let ports = row.count(MAX_PORT_WEIGHTS)?;
            if modes == 0 || ports == 0 {
                return Err(input(row.line, "each voice needs modes and force ports"));
            }
            let weights = modes.checked_mul(ports)
                .filter(|n| *n <= MAX_PORT_WEIGHTS - total_weights)
                .ok_or_else(|| input(row.line, "total port-weight budget exceeded"))?;
            row.finish()?;
            total_modes += modes;
            total_weights += weights;
            let mut model_modes = Vec::with_capacity(modes);
            let mut states = Vec::with_capacity(modes);
            for _ in 0..modes {
                let mut row = reader.row("mode")?;
                model_modes.push(ModalAcousticMode {
                    angular_frequency_rad_s: row.scalar()?,
                    damping_ratio: row.scalar()?,
                    pressure_per_modal_velocity: C64::new(row.scalar()?, row.scalar()?),
                });
                let state = ModalAcousticState {
                    displacement_m_sqrt_kg: row.scalar()?,
                    velocity_m_sqrt_kg_per_s: row.scalar()?,
                };
                if initialization == ForceInitialization::StaticPreload
                    && (state.displacement_m_sqrt_kg != 0.0 || state.velocity_m_sqrt_kg_per_s != 0.0) {
                    return Err(input(row.line, "static preload cannot discard nonzero initial Q/V"));
                }
                states.push(state);
                row.finish()?;
            }
            let mut columns = Vec::with_capacity(ports);
            let mut initial_forces = Vec::with_capacity(ports);
            for _ in 0..ports {
                let mut row = reader.row("port")?;
                initial_forces.push(row.scalar()?);
                let mut column = Vec::with_capacity(modes);
                for _ in 0..modes { column.push(row.scalar()?); }
                row.finish()?;
                columns.push(column);
            }
            let mut model = ModalAcousticTimeModel::try_new(sample_rate_hz, model_modes, budget)
                .map_err(RenderError::Modal)?;
            model.restore_states(&states).map_err(RenderError::Modal)?;
            voices.push(ModalForceVoice::new(model, columns, initial_forces, initialization)?);
        }
        let event_count: usize = reader.one("events")?;
        if event_count > MAX_EVENTS {
            return Err(input(reader.line, "force event count exceeds 65536"));
        }
        let mut events = Vec::with_capacity(event_count);
        for _ in 0..event_count {
            let mut row = reader.row("force")?;
            let event = ModalForceEvent {
                sample: row.parse()?, voice: row.parse()?, port: row.parse()?, force_n: row.scalar()?,
            };
            if event.sample >= samples {
                return Err(input(row.line, "force event is outside the declared half-open render window"));
            }
            row.finish()?;
            events.push(event);
        }
        if reader.lines.next().is_some() {
            return Err(input(reader.line + 1, "unexpected trailing record"));
        }
        let renderer = ScheduledRenderer::from_modal_forces(voices, events, ForceRenderConfig {
            sample_rate_hz, max_block, max_events: event_count, max_controls, max_projection_terms,
        })?;
        Ok(Self {
            info: ModalPerformanceInfo { sample_rate_hz, samples, full_scale_pa,
                input_hash: hash_domain(MODAL_PERFORMANCE_HASH_DOMAIN, bytes),
                voices: voice_count, modes: total_modes, force_events: event_count },
            renderer,
        })
    }

    /// Read the frozen input description without advancing the renderer.
    #[must_use]
    pub const fn info(&self) -> ModalPerformanceInfo { self.info }

    /// Move the admitted runtime to the existing scheduler/audio export APIs.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }
}

struct Reader<'a> { lines: Lines<'a>, line: usize }
impl<'a> Reader<'a> {
    fn row(&mut self, key: &str) -> Result<Row<'a>, ModalPerformanceError> {
        self.line += 1;
        let line = self.lines.next().ok_or_else(|| input(self.line, "missing record"))?;
        let mut fields = line.split_ascii_whitespace();
        if fields.next() != Some(key) { return Err(input(self.line, "unexpected record kind/order")); }
        Ok(Row { fields, line: self.line })
    }
    fn one<T: FromStr>(&mut self, key: &str) -> Result<T, ModalPerformanceError> {
        let mut row = self.row(key)?;
        let value = row.parse()?;
        row.finish()?;
        Ok(value)
    }
}
struct Row<'a> { fields: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, ModalPerformanceError> {
        self.fields.next().ok_or_else(|| input(self.line, "missing field"))
    }
    fn parse<T: FromStr>(&mut self) -> Result<T, ModalPerformanceError> {
        self.word()?.parse().map_err(|_| input(self.line, "invalid number"))
    }
    fn scalar(&mut self) -> Result<f64, ModalPerformanceError> {
        let value: f64 = self.parse()?;
        if !value.is_finite() { return Err(input(self.line, "physical values must be finite")); }
        Ok(value)
    }
    fn count(&mut self, maximum: usize) -> Result<usize, ModalPerformanceError> {
        let count = self.parse()?;
        if count > maximum { return Err(input(self.line, "count exceeds the admitted budget")); }
        Ok(count)
    }
    fn finish(mut self) -> Result<(), ModalPerformanceError> {
        if self.fields.next().is_some() { return Err(input(self.line, "unexpected extra field")); }
        Ok(())
    }
}
fn input(line: usize, what: &'static str) -> ModalPerformanceError {
    ModalPerformanceError::Input { line, what }
}

#[cfg(test)]
mod tests;

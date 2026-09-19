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
//!
//! Version 2 keeps these records and inserts `coupling_limits`, `connections`,
//! then `connection K C REST` / `left COMPONENT SHAPES...` /
//! `right COMPONENT SHAPES...` before `events`. All limits are explicit; the
//! example and complete units are in examples/COUPLED_MODAL_PERFORMANCES.md.
//! Its preloads solve the complete network; mixed retain/preload modes refuse.
//!
//! An additive voice form declares one untethered translation explicitly:
//! `voice free-mass 1 PORT_COUNT`, followed by `mass KG X_M V_M_S`, then
//! the usual port records. Its one coordinate is q=sqrt(m)*x. Port and
//! attachment shapes retain their original 1/sqrt(kg) units; a unit physical
//! translation uses 1/sqrt(m). Free motion has no direct pressure transfer.
//! Existing `mode` rows still require strictly positive natural frequencies.
//! `voice free-mass-preload 1 PORT_COUNT` uses the same mass row with zero X/V,
//! then solves its declared supporting network before the window. Every other
//! component must also request static preload. An unsupported mass still refuses.

// Version 3 adds exactly one implicit compliant contact before the events.
// Its initial states are retained; nonlinear static preload is not inferred.
// Version 4 adds multi_contact_limits and contacts COUNT, then repeats the
// unchanged version-3 contact records for a simultaneously solved contact set.
mod contact;
// Version 5 adds an explicit friction declaration for every normal contact.
mod friction;

use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_blake3::{ContentHash, hash_domain};
use fs_math::c64::C64;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::RenderError;
use super::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use super::super::ScheduledRenderer;
use super::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_exec::CancelGate;

/// Schema token at the beginning of each file.
pub const MODAL_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v1";
/// Version 2 adds explicit bilateral spring/damper connections between components.
pub const MODAL_COUPLED_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v2";
/// Domain-separated exact byte identity for the coupled model schema.
pub const MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN: &str = "org.frankensim.fs-couple.modal-performance-input.v2";
/// Version 3 adds one unilateral contact to the bilateral network.
pub const MODAL_CONTACT_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v3";
/// Exact byte identity for the contact-enabled input.
pub const MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN: &str = "org.frankensim.fs-couple.modal-performance-input.v3";
/// Version 4 adds simultaneously resolved compliant normal contacts.
pub const MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v4";
/// Exact byte identity for the multi-contact input.
pub const MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN: &str = "org.frankensim.fs-couple.modal-performance-input.v4";
/// Version 5 adds explicitly authored 1-D regularized Coulomb friction.
pub const MODAL_FRICTION_PERFORMANCE_SCHEMA: &str = "frankensim-modal-performance-v5";
/// Exact byte identity including all friction coefficients, shapes and sources.
pub const MODAL_FRICTION_PERFORMANCE_HASH_DOMAIN: &str = "org.frankensim.fs-couple.modal-performance-input.v5";
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
    /// Exact admitted schema, distinguishing independent and coupled mechanics.
    pub schema: &'static str,
    /// Declared bilateral connections, separate from normal contacts.
    pub connections: usize,
    /// Declared normal contacts: zero in v1/v2, one in v3, a bounded set in v4/v5.
    pub contacts: usize,
    /// Authored friction laws, including zero-coefficient controls; zero in v1-v4.
    /// This is not the number of contacts currently touching or sliding.
    pub friction_contacts: usize,
    /// Audio samples per second; not inferred from a note or output extension.
    pub sample_rate_hz: u32,
    /// Exact number of output samples, including a final short callback.
    pub samples: u64,
    /// Declared pascals mapped to positive full-scale PCM.
    pub full_scale_pa: f64,
    /// Input bytes under the matching version-specific hash domain.
    pub input_hash: ContentHash,
    /// Source components, independent in v1 or connected mechanically thereafter.
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
    /// Decode, bind and compile thewhole performance before returning a runtime.
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
        let schema = match text.lines().next().and_then(|line| line.split_ascii_whitespace().next()) {
            Some(MODAL_PERFORMANCE_SCHEMA) => MODAL_PERFORMANCE_SCHEMA,
            Some(MODAL_COUPLED_PERFORMANCE_SCHEMA) => MODAL_COUPLED_PERFORMANCE_SCHEMA,
            Some(MODAL_CONTACT_PERFORMANCE_SCHEMA) => MODAL_CONTACT_PERFORMANCE_SCHEMA,
            Some(MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA) => MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA,
            Some(MODAL_FRICTION_PERFORMANCE_SCHEMA) => MODAL_FRICTION_PERFORMANCE_SCHEMA,
            _ => return Err(input(1, "unsupported modal performance schema")),
        };
        reader.row(schema)?.finish()?;
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
            let kind = row.word()?;
            let is_free_mass = matches!(kind, "free-mass" | "free-mass-preload");
            let initialization = match kind {
                "free-mass" | "retain-state" => ForceInitialization::RetainState,
                "static-preload" | "free-mass-preload" => ForceInitialization::StaticPreload,
                _ => return Err(input(row.line, "expected retain-state, static-preload, free-mass or free-mass-preload")),
            };
            let modes = row.count(MAX_MODES - total_modes)?;
            let ports = row.count(MAX_PORT_WEIGHTS)?;
            if modes == 0 || ports == 0 {
                return Err(input(row.line, "each voice needs modes and force ports"));
            }
            if is_free_mass && modes != 1 {
                return Err(input(row.line, "free-mass requires exactly one translational coordinate"));
            }
            let weights = modes.checked_mul(ports)
                .filter(|n| *n <= MAX_PORT_WEIGHTS - total_weights)
                .ok_or_else(|| input(row.line, "total port-weight budget exceeded"))?;
            row.finish()?;
            total_modes += modes;
            total_weights += weights;
            let model = if is_free_mass {
                let mut row = reader.row("mass")?;
                let mass_kg = row.scalar()?;
                let displacement_m = row.scalar()?;
                let velocity_m_s = row.scalar()?;
                if initialization == ForceInitialization::StaticPreload
                    && (displacement_m != 0.0 || velocity_m_s != 0.0) {
                    return Err(input(row.line, "free-mass preload cannot discard nonzero initial position/velocity"));
                }
                row.finish()?;
                ModalAcousticTimeModel::try_free_mass(sample_rate_hz, mass_kg,
                    displacement_m, velocity_m_s, budget).map_err(RenderError::Modal)?
            } else {
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
                let mut model = ModalAcousticTimeModel::try_new(sample_rate_hz, model_modes, budget)
                    .map_err(RenderError::Modal)?;
                model.restore_states(&states).map_err(RenderError::Modal)?;
                model
            };
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
            voices.push(ModalForceVoice::new(model, columns, initial_forces, initialization)?);
        }
        // V2/v3/v4/v5 have connection records. V1 keeps its exact compiler and does
        // not reinterpret any formerly accepted independent performance.
        let coupled = if schema != MODAL_PERFORMANCE_SCHEMA {
            let mut row = reader.row("coupling_limits")?;
            let max_connections = row.count(64)?;
            let max_setup_terms = row.count(MAX_PROJECTION_TERMS)?;
            let coupling = ModalCouplingConfig {
                max_modes: MAX_MODES, max_connections, max_setup_terms,
                nyquist_guard_fraction: row.scalar()?,
                maximum_total_energy_j: row.scalar()?,
                maximum_abs_pressure_pa: row.scalar()?,
                maximum_abs_connection_force_n: row.scalar()?,
                solve_relative_tolerance: row.scalar()?,
                energy_absolute_tolerance_j: row.scalar()?,
                energy_relative_tolerance: row.scalar()?,
            };
            row.finish()?;
            let count: usize = reader.one("connections")?;
            if count > max_connections { return Err(input(reader.line, "connection count exceeds coupling_limits")); }
            let mut connections = Vec::with_capacity(count);
            for _ in 0..count {
                let mut row = reader.row("connection")?;
                let stiffness_n_m = row.scalar()?;
                let damping_n_s_m = row.scalar()?;
                let rest_extension_m = row.scalar()?;
                row.finish()?;
                let left = read_attachment(&mut reader, "left", &voices, &mut total_weights)?;
                let right = read_attachment(&mut reader, "right", &voices, &mut total_weights)?;
                connections.push(ModalConnection { left, right, stiffness_n_m, damping_n_s_m, rest_extension_m });
            }
            Some((connections, coupling))
        } else { None };
        let contact = if schema == MODAL_CONTACT_PERFORMANCE_SCHEMA {
            Some(contact::read(&mut reader, &voices, &mut total_weights)?)
        } else { None };
        let multiple = if schema == MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA
            || schema == MODAL_FRICTION_PERFORMANCE_SCHEMA {
            Some(contact::read_set(&mut reader, &voices, &mut total_weights)?)
        } else { None };
        let contact_count = if contact.is_some() { 1 } else {
            multiple.as_ref().map_or(0, |(contacts, _)| contacts.len())
        };
        let frictions = if schema == MODAL_FRICTION_PERFORMANCE_SCHEMA {
            let (contacts, _) = multiple.as_ref()
                .ok_or_else(|| input(reader.line, "friction requires a complete normal-contact set"))?;
            Some(friction::read(&mut reader, &voices, contacts, &mut total_weights)?)
        } else { None };
        let friction_count = frictions.as_ref().map_or(0, |items| items.iter().flatten().count());
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
        let force_config = ForceRenderConfig {
            sample_rate_hz, max_block, max_events: event_count, max_controls, max_projection_terms,
        };
        let (renderer, connection_count, domain) = match coupled {
            Some((connections, coupling)) => {
                let count = connections.len();
                match multiple {
                    Some((contacts, contact_set)) => match frictions {
                        Some(frictions) => (ScheduledRenderer::from_frictional_modal_forces(
                            voices, events, force_config, connections, coupling, contacts, contact_set,
                            frictions, &CancelGate::new())?, count, MODAL_FRICTION_PERFORMANCE_HASH_DOMAIN),
                        None => (ScheduledRenderer::from_multi_contact_modal_forces(
                            voices, events, force_config, connections, coupling, contacts, contact_set,
                            &CancelGate::new())?, count, MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN),
                    },
                    None => match contact {
                        Some((contact, contact_config)) => (ScheduledRenderer::from_contact_modal_forces(
                            voices, events, force_config, connections, coupling, contact, contact_config,
                            &CancelGate::new())?, count, MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN),
                        None => (ScheduledRenderer::from_coupled_modal_forces(voices, events, force_config,
                            connections, coupling, &CancelGate::new())?, count, MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN),
                    },
                }
            }
            None => (ScheduledRenderer::from_modal_forces(voices, events, force_config)?, 0, MODAL_PERFORMANCE_HASH_DOMAIN),
        };
        Ok(Self {
            info: ModalPerformanceInfo { schema, connections: connection_count, contacts: contact_count, friction_contacts: friction_count,
                sample_rate_hz, samples, full_scale_pa, input_hash: hash_domain(domain, bytes),
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

fn read_attachment(reader: &mut Reader<'_>, key: &str, voices: &[ModalForceVoice], total_weights: &mut usize)
    -> Result<ModalAttachment, ModalPerformanceError>
{
    let mut row = reader.row(key)?;
    let component: usize = row.parse()?;
    let model = voices.get(component).ok_or_else(|| input(row.line, "connection attachment names an unknown voice"))?;
    let count = model.model.modes().len();
    if count > MAX_PORT_WEIGHTS - *total_weights {
        return Err(input(row.line, "combined actuator/connection shape budget exceeded"));
    }
    *total_weights += count;
    let mut shapes = Vec::with_capacity(count);
    for _ in 0..count { shapes.push(row.scalar()?); }
    row.finish()?;
    Ok(ModalAttachment { component, shapes })
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

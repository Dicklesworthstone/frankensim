//! File-driven jet, gas, cylindrical resonator and physical outlet pressure.
//!
//! The versioned source selects existing owners, including fractional travel,
//! optional gas-derived wide-tube loss, and the multirate Rayleigh receiver.
//! It supplies no instrument-name defaults. This is a rigid cylindrical pipe
//! with an estimated tonal jet and a prescribed terminal reflection, not a
//! measured Boehm flute, resolved turbulent jet or matched radiation load.

use super::{FlueOutletConfig, FlueOutletPerformance, FlueRadiationInput};
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::viscothermal::ViscothermalSectionSpec;
use crate::flue_loop::{FlueError, FlueGeometry, FlueVoice, network::FlueNetworkSpec};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::render::schedule::GestureCompileError;
use fs_aeroac::jetcard::{JetCard, mint_tonal_interim_card};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources, read_gas_record};
use fs_scenario::gesture::GestureSchedule;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};

/// Versioned SI source grammar; the suffix is the existing gesture format.
pub const FLUE_PERFORMANCE_SCHEMA: &str = "frankensim-flue-performance-v1";
/// Complete source including an optional inline card and gesture suffix.
pub const MAX_FLUE_PERFORMANCE_BYTES: usize = 1024 * 1024;

/// Input failures retain their source line or original physics owner.
#[derive(Debug)]
pub enum FlueInputError {
    /// Malformed, unsupported or over-budget declaration.
    Input { /// One-based line, or zero for whole-file admission.
        line: usize, /// Refused condition.
        what: &'static str },
    /// Original material parsing, composition or phase refusal.
    Gas { /// One-based declaration line.
        line: usize, /// Shared material-owner error.
        source: GasInputError },
    /// Original jet-card schema or authority refusal.
    Card { /// One-based card line.
        line: usize, /// Shared source-card error.
        source: fs_aeroac::AeroacError },
    /// Geometry, jet, propagation or numerical admission failed.
    Physics(FlueError),
    /// Gesture, receiver or finite performance admission failed.
    Performance(GestureCompileError),
}
impl core::fmt::Display for FlueInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "flue input line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "flue gas line {line}: {source}"),
            Self::Card { line, source } => write!(f, "flue card line {line}: {source}"),
            Self::Physics(source) => write!(f, "flue physics: {source}"),
            Self::Performance(source) => write!(f, "flue performance: {source}"),
        }
    }
}
impl std::error::Error for FlueInputError {}
fn bad(line: usize, what: &'static str) -> FlueInputError { FlueInputError::Input { line, what } }
fn poll(gate: &CancelGate) -> Result<(), FlueInputError> {
    if gate.is_requested() { Err(bad(0, "flue source construction cancelled")) } else { Ok(()) }
}

/// Exact source identity and clocks; none is a measured-instrument certificate.
#[derive(Clone, Copy, Debug)]
pub struct FluePerformanceInfo {
    /// Hash of all original source bytes, including inline card and gestures.
    pub input_hash: ContentHash,
    /// Jet and waveguide sample clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Physical microphone output clock [Hz].
    pub sample_rate_hz: u32,
    /// Finite OUTPUT samples, including the caller's declared tail window.
    pub samples: u64,
    /// Positive pascals mapped to PCM full scale, never peak normalization.
    pub full_scale_pa: f64,
    /// Existing card owner's content identity.
    pub jet_card_hash: u64,
    /// Exact explicit selector, or an inline canonical card.
    pub jet_card_source: &'static str,
    /// Blowing-pressure assignments on the mechanical clock.
    pub compiled_controls: usize,
}

/// A complete source-to-microphone runtime and its original material records.
pub struct FlueFilePerformance {
    info: FluePerformanceInfo,
    renderer: FlueOutletPerformance,
    geometry: FlueGeometry,
    card: JetCard,
    gas_records: Vec<GasRecord>,
}
impl FlueFilePerformance {
    /// Parse the fixed-order source shown in `examples/flue.performance`.
    /// `jet_card tonal-interim` explicitly selects the retained stage-I card;
    /// `jet_card json {CANONICAL_JSON}` supplies a complete one-line card.
    /// Gains, seed, dimensions, limits, clocks and the pressure phrase are data.
    /// `exterior same` shares the bore gas; another gas record changes only the
    /// receiver environment at the same static pressure.
    ///
    /// The original jet card keeps its own narrow lab scope. The complete
    /// authored composition remains EstimateOnly, including outside-card
    /// diagnostics; loading never promotes a flute accuracy or SPL claim.
    ///
    /// # Errors
    /// Malformed/over-budget input, missing card, original material, geometry,
    /// jet, fractional/loss/receiver domain, gesture or cancellation refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize, gate: &CancelGate)
        -> Result<Self, FlueInputError>
    {
        Self::from_bytes_with_gas_sources(bytes, max_block, gate, None)
    }

    /// The same path with explicitly loaded NASA-9/transport source packs.
    /// All gas forms use `fs_material::gas_species::input::read_gas_record`.
    pub fn from_bytes_with_gas_sources(bytes: &[u8], max_block: usize, gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>)
        -> Result<Self, FlueInputError>
    {
        poll(gate)?;
        if bytes.len() > MAX_FLUE_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "flue source exceeds its byte or callback budget"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "flue source must be UTF-8"))?;
        let (configuration, gestures) = text.split_once("\ngesture\n")
            .ok_or_else(|| bad(0, "missing LF-delimited gesture section"))?;
        let (header, body) = configuration.split_once('\n').ok_or_else(|| bad(1, "missing flue schema header"))?;
        if header != FLUE_PERFORMANCE_SCHEMA { return Err(bad(1, "unsupported flue source schema")); }
        let mut rows = Rows { lines: body.lines(), line: 1 };
        let mut row = rows.row("clock")?;
        let mechanics_rate_hz: u32 = row.number()?;
        let sample_rate_hz: u32 = row.number()?;
        row.end()?;
        if !(8_000..=192_000).contains(&sample_rate_hz) || mechanics_rate_hz == 0
            || mechanics_rate_hz % sample_rate_hz != 0
            || !(1..=16).contains(&(mechanics_rate_hz/sample_rate_hz)) {
            return Err(bad(rows.line, "output clock must be 8..192 kHz and mechanics an integer multiple in 1..=16"));
        }
        let mut row = rows.row("output")?;
        let samples = row.number()?;
        let full_scale_pa = row.scalar()?;
        let maximum_abs_pressure_pa = row.scalar()?;
        row.end()?;
        if full_scale_pa <= 0.0 { return Err(bad(rows.line, "full scale must be positive pascals")); }
        let mut gas_records = Vec::new();
        let mut row = rows.row("gas")?;
        let gas = row.gas(sources, &mut gas_records)?;
        row.end()?;
        let exterior_record = rows.record("exterior")?;
        let exterior = if exterior_record == "same" {
            gas
        } else {
            let mut row = Row { fields: exterior_record.split_ascii_whitespace(), line: rows.line };
            let state = row.gas(sources, &mut gas_records)?;
            row.end()?;
            state
        };
        let mut row = rows.row("tube")?;
        let length_m = row.scalar()?;
        let radius_m = row.scalar()?;
        let max_length_error_m = row.scalar()?;
        let reflection = row.scalar()?;
        row.end()?;
        let mut row = rows.row("flue")?;
        let geometry = FlueGeometry {
            flue_width_m: row.scalar()?, flue_height_m: row.scalar()?, cut_up_m: row.scalar()?,
            labium_offset_m: row.scalar()?, jet_angle_rad: row.scalar()?, mouth_loss_rel: row.scalar()?,
            bore_area_m2: core::f64::consts::PI * radius_m * radius_m,
        };
        row.end()?;
        let card_source = rows.record("jet_card")?;
        let (card, jet_card_source) = if card_source == "tonal-interim" {
            (mint_tonal_interim_card(), "tonal-interim")
        } else if let Some(json) = card_source.strip_prefix("json ") {
            (JetCard::from_json(json), "inline-json")
        } else { return Err(bad(rows.line, "jet_card requires tonal-interim or json followed by a complete canonical card")); };
        let card = card.map_err(|source| FlueInputError::Card { line: rows.line, source })?;
        let mut row = rows.row("jet")?;
        let displacement_gain = row.scalar()?;
        let seed_amplitude_rel = row.scalar()?;
        let voice_seed = row.number()?;
        let voice_index = row.number()?;
        row.end()?;
        let mut row = rows.row("jet_limits")?;
        let minimum_jet_speed_m_s = row.scalar()?;
        let maximum_jet_mach = row.scalar()?;
        let max_jet_history_samples = row.number()?;
        row.end()?;
        let mut row = rows.row("wave_memory_bytes")?;
        let max_wave_memory_bytes: usize = row.number()?;
        row.end()?;
        if max_wave_memory_bytes > 256*1024*1024 { return Err(bad(rows.line, "flue wave memory exceeds 256 MiB")); }
        let mut row = rows.row("fractional")?;
        let policy = FractionalTubePolicy { maximum_frequency_hz: row.scalar()?,
            maximum_summed_phase_error_rad: row.scalar()? };
        row.end()?;
        let loss_record = rows.record("viscothermal")?;
        let viscothermal = if loss_record == "none" {
            Vec::new()
        } else {
            let mut row = Row { fields: loss_record.split_ascii_whitespace(), line: rows.line };
            let loss = ViscothermalSectionSpec { section: 0, cells: row.number()?,
                minimum_frequency_hz: row.scalar()?, maximum_frequency_hz: row.scalar()? };
            row.end()?;
            vec![loss]
        };
        let mut row = rows.row("receiver")?;
        let receiver = CircularOutletReceiver { position_m: [row.scalar()?,row.scalar()?,row.scalar()?],
            radial_rings: row.number()?, angular_points: row.number()?, maximum_frequency_hz: row.scalar()? };
        row.end()?;
        if receiver.maximum_frequency_hz > policy.maximum_frequency_hz {
            return Err(bad(rows.line, "fractional phase band must cover the receiver's declared band"));
        }
        let mut row = rows.row("controls")?;
        let max_controls: usize = row.number()?;
        let max_compile_work = row.number()?;
        row.end()?;
        if max_controls > 262_144 { return Err(bad(rows.line, "flue control budget exceeds 262144 ticks")); }
        if rows.lines.next().is_some() { return Err(bad(rows.line+1,"unexpected records before gesture section")); }
        for line in gestures.lines() {
            if let Some(count) = line.strip_prefix("events\t") {
                if count.parse::<usize>().ok().is_none_or(|n| n>16384) {
                    return Err(bad(rows.line+1,"gesture event count exceeds 16384-event input budget"));
                }
            }
        }
        let schedule = GestureSchedule::from_canonical_bytes(gestures.as_bytes())
            .map_err(|_| bad(rows.line+1, "invalid canonical blowing-pressure gesture suffix"))?;
        if schedule.to_canonical_bytes() != gestures.as_bytes() {
            return Err(bad(rows.line+1,"gesture suffix must be canonical without ignored trailing data"));
        }
        poll(gate)?;
        let source = FlueNetworkSpec {
            graph: TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
                sections: vec![TubeSection { nodes: [0,1], length_m, radius_m, max_length_error_m }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes },
            gases: vec![gas], viscothermal, minimum_jet_speed_m_s, maximum_jet_mach, max_jet_history_samples,
        };
        let voice = FlueVoice::from_fractional_network(Some(&card), displacement_gain, seed_amplitude_rel,
            geometry, source, None, policy, mechanics_rate_hz, voice_seed, voice_index, gate)
            .map_err(FlueInputError::Physics)?;
        let renderer = FlueOutletPerformance::with_multirate_radiation(voice, schedule,
            FlueOutletConfig { sample_rate_hz, samples, max_block, max_controls, max_compile_work,
                outlet_node: 1, receiver, maximum_abs_pressure_pa },
            FlueRadiationInput { exterior, additional_outlets: Vec::new() }, mechanics_rate_hz)
            .map_err(FlueInputError::Performance)?;
        poll(gate)?;
        let info = FluePerformanceInfo { input_hash: hash_domain("org.frankensim.fs-couple.flue-performance.v1",bytes),
            mechanics_rate_hz, sample_rate_hz, samples, full_scale_pa,
            jet_card_hash: card.content_hash(), jet_card_source, compiled_controls: renderer.controls.len() };
        Ok(Self { info, renderer, geometry, card, gas_records })
    }
    /// Original identity and finite clocks, separate from the evolving solver.
    #[must_use] pub const fn info(&self) -> FluePerformanceInfo { self.info }
    /// Actual physical runtime; inspect its original and lowered gas/tube inputs.
    #[must_use] pub const fn renderer(&self) -> &FlueOutletPerformance { &self.renderer }
    /// Supplied mouth geometry, with bore area derived from the same tube radius.
    #[must_use] pub const fn geometry(&self) -> FlueGeometry { self.geometry }
    /// The complete selected card, including its narrower lab-validity scope.
    #[must_use] pub const fn card(&self) -> &JetCard { &self.card }
    /// Gas records and thermochemical receipts, in source order.
    #[must_use] pub fn gas_records(&self) -> &[GasRecord] { &self.gas_records }
    /// Move the original physical pressure producer into streaming/PCM output.
    #[must_use] pub fn into_renderer(self) -> FlueOutletPerformance { self.renderer }
}

struct Rows<'a> { lines: Lines<'a>, line: usize }
impl<'a> Rows<'a> {
    fn record(&mut self, expected: &'static str) -> Result<&'a str, FlueInputError> {
        self.line += 1;
        let text = self.lines.next().ok_or_else(|| bad(self.line,"missing source record"))?.trim_ascii();
        let key = text.split_ascii_whitespace().next().unwrap_or("");
        if key != expected { return Err(bad(self.line, expected)); }
        Ok(text[key.len()..].trim_ascii())
    }
    fn row(&mut self, expected: &'static str) -> Result<Row<'a>, FlueInputError> {
        let text = self.record(expected)?;
        Ok(Row { fields: text.split_ascii_whitespace(), line: self.line })
    }
}
struct Row<'a> { fields: SplitAsciiWhitespace<'a>, line: usize }
impl Row<'_> {
    fn number<T: FromStr>(&mut self) -> Result<T, FlueInputError> {
        self.fields.next().and_then(|s| s.parse().ok()).ok_or_else(|| bad(self.line,"missing or invalid numeric field"))
    }
    fn scalar(&mut self) -> Result<f64, FlueInputError> {
        let x: f64 = self.number()?;
        if !x.is_finite() { return Err(bad(self.line,"physical scalars must be finite")); }
        Ok(x)
    }
    fn end(mut self) -> Result<(), FlueInputError> {
        if self.fields.next().is_some() { Err(bad(self.line,"unexpected trailing fields")) } else { Ok(()) }
    }
    fn gas(&mut self, sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>)
        -> Result<GasState, FlueInputError>
    {
        let record = read_gas_record(&mut self.fields, None, sources)
            .map_err(|source| FlueInputError::Gas { line: self.line, source })?;
        let gas = *record.state();
        records.push(record);
        Ok(gas)
    }
}

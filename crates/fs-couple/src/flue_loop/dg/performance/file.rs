//! Explicit DG bore and original jet source through the shared wind renderer.
//!
//! This grammar uses the same gas records, DG source loader, flue control
//! compiler, source equations and PressureRenderer as the existing instrument
//! paths. Observation is declared internal mean-inlet pressure. The jet card
//! remains an exploration source; it is never a resolved-flow qualification.

use super::{FlueDgConfig, FlueDgPerformance};
use crate::dg_port::input::{DgSourceError, DgSourceInput};
use crate::flue_loop::dg::FlueDgSpec;
use crate::flue_loop::network::performance::file::{
    FlueInputError, FluePerformanceInfo, MAX_FLUE_PERFORMANCE_BYTES, Rows,
};
use crate::flue_loop::{FlueGeometry, FlueVoice};
use fs_aeroac::jetcard::{JetCard, mint_tonal_interim_card};
use fs_blake3::hash_domain;
use fs_exec::CancelGate;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources};
use fs_scenario::gesture::GestureSchedule;

pub const FLUE_DG_PERFORMANCE_SCHEMA: &str = "frankensim-flue-dg-performance-v1";

/// Original source/parser refusals stay attached to their physical owner.
#[derive(Debug)]
pub enum FlueDgInputError {
    Source(FlueInputError),
    Duct(DgSourceError),
}
impl core::fmt::Display for FlueDgInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "{error}"),
            Self::Duct(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for FlueDgInputError {}
impl From<FlueInputError> for FlueDgInputError {
    fn from(error: FlueInputError) -> Self { Self::Source(error) }
}
impl From<DgSourceError> for FlueDgInputError {
    fn from(error: DgSourceError) -> Self { Self::Duct(error) }
}
fn bad(line: usize, what: &'static str) -> FlueDgInputError {
    FlueInputError::Input { line, what }.into()
}
fn poll(gate: &CancelGate) -> Result<(), FlueDgInputError> {
    if gate.is_requested() { Err(bad(0, "DG flue source construction cancelled")) } else { Ok(()) }
}

pub struct FlueDgFilePerformance {
    info: FluePerformanceInfo,
    renderer: FlueDgPerformance,
    geometry: FlueGeometry,
    card: JetCard,
    gas_records: Vec<GasRecord>,
}

impl FlueDgFilePerformance {
    pub fn from_bytes(
        bytes: &[u8], max_block: usize, gate: &CancelGate,
    ) -> Result<Self, FlueDgInputError> {
        Self::from_bytes_with_gas_sources(bytes, max_block, gate, None)
    }

    /// All geometry, gas, solver limits, source controls and pressure scale are
    /// authored. dg_tube and dg_bore use the exact shared DG input grammar.
    /// The bore's actual inlet area supplies the mouth's characteristic basis.
    /// No exterior gas or receiver is inferred from internal pressure.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, FlueDgInputError> {
        poll(gate)?;
        if bytes.len() > MAX_FLUE_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "DG flue source exceeds its byte or callback budget"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "DG flue source must be UTF-8"))?;
        let (configuration, gestures) = text.split_once("\ngesture\n")
            .ok_or_else(|| bad(0, "missing LF-delimited gesture section"))?;
        let (header, body) = configuration.split_once('\n')
            .ok_or_else(|| bad(1, "missing DG flue schema header"))?;
        if header != FLUE_DG_PERFORMANCE_SCHEMA {
            return Err(bad(1, "unsupported DG flue source schema"));
        }
        let mut rows = Rows { lines: body.lines(), line: 1 };
        let mut row = rows.row("clock")?;
        let sample_rate_hz: u32 = row.number()?;
        row.end()?;
        if !(8_000..=3_072_000).contains(&sample_rate_hz) {
            return Err(bad(rows.line, "DG flue mechanical clock must be 8..3072 kHz"));
        }
        let mut row = rows.row("output")?;
        let samples: u64 = row.number()?;
        let full_scale_pa = row.scalar()?;
        let maximum_abs_pressure_pa = row.scalar()?;
        row.end()?;
        if samples == 0 || samples > 28_800_000 || full_scale_pa <= 0.0 {
            return Err(bad(rows.line, "DG flue needs 1..28800000 source samples and positive PCM pascals"));
        }
        let mut gas_records = Vec::new();
        let mut row = rows.row("gas")?;
        let gas = row.gas(sources, &mut gas_records)?;
        row.end()?;
        let kind = if rows.next_is("dg_tube") { "dg_tube" } else { "dg_bore" };
        let mut row = rows.row(kind)?;
        let source = DgSourceInput::read(
            kind, row.line, &mut row.fields,
            &mut || {
                rows.line += 1;
                let text = rows.lines.next()
                    .ok_or_else(|| bad(rows.line, "missing DG source geometry/gas record"))?;
                Ok::<_, FlueDgInputError>((rows.line, text.split_ascii_whitespace()))
            },
            &gas, sources, &mut gas_records, gate,
        )?;
        row.end()?;
        let mut row = rows.row("flue")?;
        let mut geometry = FlueGeometry {
            flue_width_m: row.scalar()?,
            flue_height_m: row.scalar()?,
            cut_up_m: row.scalar()?,
            labium_offset_m: row.scalar()?,
            jet_angle_rad: row.scalar()?,
            mouth_loss_rel: row.scalar()?,
            // Filled from the already admitted original DG inlet below.
            bore_area_m2: 0.0,
        };
        row.end()?;
        let mut row = rows.row("windway")?;
        let windway = match row.fields.next() {
            Some("none") => None,
            Some("slit") => Some((row.scalar()?, None)),
            Some("unsteady") => Some((row.scalar()?, Some((row.number()?, row.number()?)))),
            _ => return Err(bad(row.line, "windway requires none, slit LENGTH or unsteady LENGTH TRANSVERSE_CELLS CONVECTION_CELLS")),
        };
        row.end()?;
        let card_source = rows.record("jet_card")?;
        let (card, jet_card_source) = if card_source == "tonal-interim" {
            (mint_tonal_interim_card(), "tonal-interim")
        } else if let Some(json) = card_source.strip_prefix("json ") {
            (JetCard::from_json(json), "inline-json")
        } else {
            return Err(bad(rows.line, "jet_card requires tonal-interim or json CANONICAL_CARD"));
        };
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
        let mut row = rows.row("controls")?;
        let max_controls: usize = row.number()?;
        let max_compile_work = row.number()?;
        row.end()?;
        if max_controls > 262_144 {
            return Err(bad(rows.line, "DG flue control budget exceeds 262144 ticks"));
        }
        if rows.record("observation")? != "inlet-pressure" {
            return Err(bad(rows.line, "DG flue currently observes internal mean inlet pressure only"));
        }
        if rows.lines.next().is_some() {
            return Err(bad(rows.line + 1, "unexpected records before gesture section"));
        }
        for line in gestures.lines() {
            if let Some(count) = line.strip_prefix("events\t") {
                if count.parse::<usize>().ok().is_none_or(|n| n > 16_384) {
                    return Err(bad(rows.line + 1, "gesture event count exceeds 16384-event input budget"));
                }
            }
        }
        let schedule = GestureSchedule::from_canonical_bytes(gestures.as_bytes())
            .map_err(|_| bad(rows.line + 1, "invalid canonical pressure/angle gesture suffix"))?;
        if schedule.to_canonical_bytes() != gestures.as_bytes() {
            return Err(bad(rows.line + 1, "gesture suffix must be canonical without ignored trailing data"));
        }
        poll(gate)?;
        let load = source.prepare(&gas, 1.0 / f64::from(sample_rate_hz), gate)?;
        geometry.bore_area_m2 = load.duct().inlet_area_m2();
        let voice = FlueVoice::from_dg(
            Some(&card), displacement_gain, seed_amplitude_rel, geometry,
            FlueDgSpec { load, minimum_jet_speed_m_s, maximum_jet_mach, max_jet_history_samples },
            sample_rate_hz, voice_seed, voice_index, gate,
        ).map_err(FlueInputError::Physics)?;
        let voice = match windway {
            None => voice,
            Some((length, None)) => voice.with_viscous_windway(length)
                .map_err(FlueInputError::Physics)?,
            Some((length, Some((transverse, convection)))) => voice
                .with_unsteady_windway(length, transverse, convection)
                .map_err(FlueInputError::Physics)?,
        };
        let renderer = FlueDgPerformance::new(voice, schedule, FlueDgConfig {
            sample_rate_hz, samples, max_block, max_controls, max_compile_work,
            maximum_abs_pressure_pa,
        }).map_err(FlueInputError::Performance)?;
        poll(gate)?;
        let info = FluePerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.flue-dg-performance.v1", bytes),
            mechanics_rate_hz: sample_rate_hz, sample_rate_hz, samples, full_scale_pa,
            jet_card_hash: card.content_hash(), jet_card_source,
            compiled_controls: renderer.compiled_controls(),
        };
        Ok(Self { info, renderer, geometry, card, gas_records })
    }

    pub const fn info(&self) -> FluePerformanceInfo { self.info }
    pub const fn renderer(&self) -> &FlueDgPerformance { &self.renderer }
    pub const fn geometry(&self) -> FlueGeometry { self.geometry }
    pub const fn card(&self) -> &JetCard { &self.card }
    pub fn gas_records(&self) -> &[GasRecord] { &self.gas_records }
    pub fn into_renderer(self) -> FlueDgPerformance { self.renderer }
}

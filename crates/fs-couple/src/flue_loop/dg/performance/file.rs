//! Explicit DG bore and original jet source through the shared wind renderer.
//!
//! This grammar uses the same gas records, DG source loader, flue control
//! compiler, source equations and PressureRenderer as the existing instrument
//! paths. Observation is explicitly internal mean-inlet pressure or the
//! computed outlet's baffled microphone in its original exterior gas. The jet
//! card remains an exploration source, never a resolved-flow qualification.
//! After observation, optional `observer_gas stokes-kirchhoff|lossless
//! MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS` selects the shared physical
//! exterior propagation for every observed opening, including the optional
//! mouth. Omission retains the original lossless receiver.
//! After a baffled observation, optional `mouth_radiation X Y WIDTH_CELLS
//! CUT_UP_CELLS` adds the actual rectangular mouth to that same microphone.
//! Optional `mouth_reaction none|compact-baffled|mutual-baffled` selects no
//! replacement, the original compact self reaction, or one complete shared
//! rectangular-mouth/circular-terminal reaction. The latter includes all
//! reciprocal cross terms and requires source quadrature over mechanical Nyquist.
//! Windway, mouth-contraction and jet-source selections share the characteristic
//! network's parser and original physical owners. Thus `jet_transfer
//! nolle-price-2015 CHANNEL CHAMFER_DEPTH_M`, the forced Bickley response and
//! card centreline speed retain the same source-profile validity boundaries.

use super::{FlueDgConfig, FlueDgPerformance, FlueMouthRadiation};
use crate::dg_port::input::{DgSourceError, DgSourceInput};
use crate::dg_port::input::wall::DgMaterialWallRecord;
use crate::dg_port::input::wall::regions::{DgMaterialWallLayout, DgMaterialWallLayoutInput};
use crate::flue_loop::dg::FlueDgSpec;
use crate::flue_loop::network::performance::file::{
    FlueInputError, FluePerformanceInfo, MAX_FLUE_PERFORMANCE_BYTES, Rows,
    held_blowing_pressure,
};
use crate::flue_loop::network::performance::file::excitation::{FlueFlowInput, FlueJetInput};
use crate::flue_loop::{FlueGeometry, FlueVoice};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_aeroac::jetcard::{JetCard, mint_tonal_interim_card};
use fs_blake3::hash_domain;
use fs_exec::CancelGate;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources};
use fs_matdb::NormalizedMaterialCardPack;
use fs_scenario::gesture::GestureSchedule;
use std::sync::Arc;
use crate::pcm_wav::baffled::gas::GasPropagation;

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

fn read_observation(rows: &mut Rows<'_>) -> Result<Option<CircularOutletReceiver>, FlueDgInputError> {
    let mut row = rows.row("observation")?;
    let receiver = match row.fields.next() {
        Some("inlet-pressure") => None,
        Some("baffled") => Some(CircularOutletReceiver {
            position_m: [row.scalar()?, row.scalar()?, row.scalar()?],
            radial_rings: row.number()?,
            angular_points: row.number()?,
            maximum_frequency_hz: row.scalar()?,
        }),
        _ => return Err(bad(row.line,
            "observation requires inlet-pressure or baffled X Y Z RINGS ANGLES FMAX_HZ")),
    };
    row.end()?;
    Ok(receiver)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MouthReaction { None, CompactBaffled, MutualBaffled }

fn read_mouth(rows: &mut Rows<'_>, exterior: bool)
    -> Result<(Option<FlueMouthRadiation>, MouthReaction), FlueDgInputError>
{
    let mouth = if rows.next_is("mouth_radiation") {
        let mut row = rows.row("mouth_radiation")?;
        let input = FlueMouthRadiation {
            center_m: [row.scalar()?, row.scalar()?],
            aperture_size_m: None,
            cells: [row.number()?, row.number()?],
        };
        row.end()?;
        if !exterior {
            return Err(bad(rows.line, "DG mouth radiation requires observation baffled"));
        }
        Some(input)
    } else { None };
    let reaction = if rows.next_is("mouth_reaction") {
        let mut row = rows.row("mouth_reaction")?;
        let selected = match row.fields.next() {
            Some("none") => MouthReaction::None,
            Some("compact-baffled") => MouthReaction::CompactBaffled,
            Some("mutual-baffled") => MouthReaction::MutualBaffled,
            _ => return Err(bad(row.line, "mouth_reaction requires none, compact-baffled or mutual-baffled")),
        };
        row.end()?;
        if selected != MouthReaction::None && mouth.is_none() {
            return Err(bad(rows.line, "DG mouth reaction requires an explicit mouth_radiation rectangle"));
        }
        selected
    } else { MouthReaction::None };
    Ok((mouth, reaction))
}

pub struct FlueDgFilePerformance {
    info: FluePerformanceInfo,
    renderer: FlueDgPerformance,
    geometry: FlueGeometry,
    card: JetCard,
    gas_records: Vec<GasRecord>,
    wall: DgMaterialWallLayout,
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
        Self::from_bytes_with_sources(bytes, max_block, gate, sources,
            |_| Err("DG solid wall requires a normalized material source loader".into()))
    }

    /// Optional `dg_wall` after the complete DG bore block resolves a solid
    /// source through the same local wall law used by moving-valve sources.
    /// The original gas diffusion, source equations and observer are retained.
    pub fn from_bytes_with_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
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
        let source = if rows.next_is("dg_molecular") {
            let mut row = rows.row("dg_molecular")?;
            let law = fs_material::relaxation::input::read_relaxation(&mut row.fields)
                .map_err(|source| FlueInputError::Gas { line: row.line,
                    source: GasInputError::Material(source) })?;
            row.end()?;
            source.with_molecular_relaxation(law, &gas_records[0])?
        } else { source };
        let wall = if rows.next_is("dg_wall") {
            let mut row = rows.row("dg_wall")?;
            let declaration = DgMaterialWallLayoutInput::read(
                row.line, &mut row.fields, &mut || {
                    rows.line += 1;
                    let text = rows.lines.next()
                        .ok_or_else(|| bad(rows.line, "missing DG material wall region"))?;
                    Ok::<_, FlueDgInputError>((rows.line, text.split_ascii_whitespace()))
                },
            )?;
            row.end()?;
            declaration.resolve(gate, &mut |path, line| {
                let pack = load_material(path).map_err(|source|
                    FlueInputError::MaterialSource { line, source })?;
                Ok::<_, FlueDgInputError>(Arc::new(pack))
            })?
        } else { DgMaterialWallLayout::default() };
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
        let flow_input = FlueFlowInput::read(&mut rows)?;
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
        let jet_input = FlueJetInput::read(&mut rows)?;
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
        let observation = read_observation(&mut rows)?;
        let observer_gas = if rows.next_is("observer_gas") {
            let mut row = rows.row("observer_gas")?;
            let exterior = source.exterior_gas_record()
                .ok_or_else(|| bad(row.line, "observer_gas requires an original DG exterior gas record"))?;
            Some(GasPropagation::read_record_for_gas(&mut row.fields, exterior)
                .map_err(|source| FlueInputError::Gas { line: row.line,
                    source: GasInputError::Material(source) })?)
        } else { None };
        let (mouth_radiation, mouth_reaction) = read_mouth(&mut rows, observation.is_some())?;
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
        let load = source.prepare_with_material_wall_layout(&gas,
            1.0 / f64::from(sample_rate_hz), gate, &wall)?;
        geometry.bore_area_m2 = load.duct().inlet_area_m2();
        let voice = FlueVoice::from_dg(
            Some(&card), displacement_gain, seed_amplitude_rel, geometry,
            FlueDgSpec { load, minimum_jet_speed_m_s, maximum_jet_mach, max_jet_history_samples },
            sample_rate_hz, voice_seed, voice_index, gate,
        ).map_err(FlueInputError::Physics)?;
        let held_pressure_pa = held_blowing_pressure(&schedule, None);
        let voice = flow_input.attach(voice, held_pressure_pa)?;
        let voice = jet_input.attach(voice, held_pressure_pa, gate)?;
        let config = FlueDgConfig {
            sample_rate_hz, samples, max_block, max_controls, max_compile_work,
            maximum_abs_pressure_pa,
        };
        let mut renderer = match observation {
            Some(receiver) => FlueDgPerformance::with_baffled_outlet(
                voice, schedule, config, receiver,
            ),
            None => FlueDgPerformance::new(voice, schedule, config),
        }.map_err(FlueInputError::Performance)?;
        if let Some(mouth) = mouth_radiation {
            renderer = renderer.with_mouth_radiation(mouth)
                .map_err(FlueInputError::Performance)?;
        }
        renderer = match mouth_reaction {
            MouthReaction::None => renderer,
            MouthReaction::CompactBaffled => renderer.with_mouth_reaction(gate)
                .map_err(FlueInputError::Performance)?,
            MouthReaction::MutualBaffled => renderer.with_mutual_mouth_reaction(gate)
                .map_err(FlueInputError::Performance)?,
        };
        if let Some(propagation) = observer_gas {
            renderer = renderer.with_observer_gas(&propagation, gate)
                .map_err(FlueInputError::Performance)?;
        }
        poll(gate)?;
        let info = FluePerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.flue-dg-performance.v1", bytes),
            mechanics_rate_hz: sample_rate_hz, sample_rate_hz, samples, full_scale_pa,
            jet_card_hash: card.content_hash(), jet_card_source,
            compiled_controls: renderer.compiled_controls(),
        };
        Ok(Self { info, renderer, geometry, card, gas_records, wall })
    }

    pub const fn info(&self) -> FluePerformanceInfo { self.info }
    pub const fn renderer(&self) -> &FlueDgPerformance { &self.renderer }
    pub const fn geometry(&self) -> FlueGeometry { self.geometry }
    pub const fn card(&self) -> &JetCard { &self.card }
    pub fn gas_records(&self) -> &[GasRecord] { &self.gas_records }
    pub const fn wall_record(&self) -> Option<&DgMaterialWallRecord> { self.wall.uniform_record() }
    /// Original physical material partition, independent of DG refinement.
    pub const fn wall_layout(&self) -> &DgMaterialWallLayout { &self.wall }
    pub fn into_renderer(self) -> FlueDgPerformance { self.renderer }
}

#[cfg(test)]
mod tests;


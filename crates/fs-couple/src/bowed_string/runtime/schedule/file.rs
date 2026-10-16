//! Authored string mechanics and a canonical bow track, using the existing
//! persistent solver and gesture sampler. This is an input adapter, not another
//! friction law or interpolator. The one-way compact-body bridge approximation
//! remains explicit; supplied numbers are not a material-validation certificate.

use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_blake3::{ContentHash, hash_domain};
use fs_material::gas::GasState;
use fs_scenario::{RadiatingPlate, gesture::{GestureSchedule, GestureTarget, GestureValue}};
use crate::bowed_string::{BowGesture, BowedRunConfig, BowedRunError, BowedStringCard, FrictionIsland, Termination};
use crate::render::{RenderContext, RenderError, RenderVoice, schedule::ScheduledRenderer};
use crate::stribeck_friction::StribeckFriction;
use crate::thin_plate::CompactBody;
use super::{BowedScheduleError, ScheduledBowedRenderer};
use super::super::BowedStringState;

/// Versioned header; the trailing schedule uses the existing canonical grammar.
pub const BOWED_PERFORMANCE_SCHEMA: &str = "frankensim-bowed-performance-v1";
/// Bound bytes BEFORE reading or parsing a source file.
pub const MAX_BOWED_PERFORMANCE_BYTES: usize = 1024 * 1024;
const MAX_MODES: usize = 512;
const MAX_EVENTS: usize = 16_384;
const MAX_CONTROLS: usize = 262_144;
const MAX_WORK: u64 = 16_777_216;

/// Syntax and existing physical/scheduling owner refusals remain distinguishable.
#[derive(Debug)]
pub enum BowedPerformanceError {
    /// First offending header line, or zero for whole-file/schedule admission.
    Input { /// One-based header line, zero for non-header admission.
        line: usize, /// Violated input condition.
        what: &'static str },
    /// The existing string/bridge/friction solver rejected the supplied state.
    Model(BowedRunError),
    /// The existing gesture compiler rejected the source or its budgets.
    Schedule(BowedScheduleError),
    /// The existing common renderer rejected the voice.
    Render(RenderError),
}
impl core::fmt::Display for BowedPerformanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "bowed input line {line}: {what}"),
            Self::Model(e) => write!(f, "bowed model: {e:?}"),
            Self::Schedule(e) => write!(f, "bowed input: {e}"),
            Self::Render(e) => write!(f, "bowed render: {e}"),
        }
    }
}
impl std::error::Error for BowedPerformanceError {}
fn bad(line: usize, what: &'static str) -> BowedPerformanceError {
    BowedPerformanceError::Input { line, what }
}

/// Source facts, not measured material identity or a physical validation claim.
#[derive(Clone, Copy, Debug)]
pub struct BowedPerformanceInfo {
    /// Domain-separated exact bytes of the complete performance file.
    pub input_hash: ContentHash,
    /// Original mechanics/output clock of the bowed solver [Hz].
    pub sample_rate_hz: u32,
    /// Complete, half-open physical window [0, samples).
    pub samples: u64,
    /// Declared PCM scale [Pa], not an acoustic gain in an ensemble.
    pub full_scale_pa: f64,
    /// Retained transverse string modes; the observer body is separate.
    pub string_modes: usize,
    /// Friction-refresh substeps per mechanical sample, not oversampled PCM.
    pub subsamples: usize,
    /// Authored bow commands in the canonical track.
    pub gesture_events: usize,
    /// Actual stored bow controls, including the initial held value.
    pub compiled_controls: usize,
}

/// File-driven performance ready for the common renderer and pressure ensembles.
pub struct BowedPerformance {
    info: BowedPerformanceInfo,
    renderer: ScheduledRenderer,
}
impl BowedPerformance {
    /// Decode the complete input before constructing a physical runtime.
    /// Header records are fixed-order and reject extra fields. The schedule
    /// must be exact canonical bytes with one bow track targeting string zero.
    /// Count caps precede decoder reservations. All authored commands must start
    /// by the last observed control tick; a ramp may continue beyond the window.
    /// No sample advances at admission, and no oscillator, gain or interpolation
    /// is inserted. The legacy positive initial-load admission is retained;
    /// a step at time zero may explicitly release the bow before its first step.
    ///
    /// # Errors
    /// Malformed/oversized input, invalid physics, unsupported tracks, unobserved
    /// commands, exhausted compilation budgets or the original owner's refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, BowedPerformanceError> {
        if bytes.len() > MAX_BOWED_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "source exceeds 1 MiB or callback capacity is outside 1..=65536"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "source must be UTF-8"))?;
        let (header, schedule_text) = text.split_once("\nschedule\n")
            .ok_or_else(|| bad(0, "expected an LF-delimited schedule record and canonical gesture bytes"))?;
        let mut r = Reader { lines: header.lines(), line: 0 };
        r.row(BOWED_PERFORMANCE_SCHEMA)?.finish()?;
        let mut row = r.row("audio")?;
        let rate: u32 = row.parse()?;
        let samples: u64 = row.parse()?;
        let full_scale_pa = row.scalar()?;
        row.finish()?;
        if !(1..=192_000).contains(&rate) || samples == 0 || samples > 28_800_000
            || samples > u64::from(rate) * 600 || full_scale_pa <= 0.0
        { return Err(bad(r.line, "audio requires a rate in 1..=192000, positive scale and at most 600 seconds/28800000 samples")); }
        let mut row = r.row("string")?;
        let mut card = BowedStringCard {
            length_m: row.scalar()?, tension_n: row.scalar()?, linear_density_kg_m: row.scalar()?,
            bending_stiffness_n_m2: row.scalar()?, viscous_bending_n_m2_s: row.scalar()?,
            mode_count: row.count(MAX_MODES)?, zetas: Vec::new(), sample_rate_hz: rate,
        };
        row.finish()?;
        if card.mode_count == 0 { return Err(bad(r.line, "the string needs at least one retained mode")); }
        let mut row = r.row("damping")?;
        for _ in 0..card.mode_count { card.zetas.push(row.scalar()?); }
        row.finish()?;
        let mut row = r.row("stribeck")?;
        let friction = StribeckFriction::try_new(row.scalar()?, row.scalar()?, row.scalar()?)
            .map_err(|what| bad(r.line, what))?;
        row.finish()?;
        let mut row = r.row("subsamples")?;
        let subsamples = row.count(256)?; row.finish()?;
        if subsamples == 0 { return Err(bad(r.line, "subsamples must be explicit and positive")); }
        let mut row = r.row("body")?;
        let radiator = RadiatingPlate { area_m2: row.scalar()?, mass_kg: row.scalar()?,
            frequency_hz: row.scalar()?, damping_ratio: row.scalar()? };
        row.finish()?;
        let mut row = r.row("ambient")?;
        let ambient = GasState::try_new_moist_air(row.scalar()?, row.scalar()?, row.scalar()?)
            .map_err(|_| bad(r.line, "ambient state is outside the shared moist-air model"))?;
        row.finish()?;
        let mut row = r.row("listener_m")?;
        let listener_m = row.scalar()?; row.finish()?;
        let mut row = r.row("compile_limits")?;
        let max_work: u64 = row.parse()?;
        let max_controls = row.count(MAX_CONTROLS)?; row.finish()?;
        if max_work > MAX_WORK { return Err(bad(r.line, "bow compilation exceeds 16777216 source visits")); }
        if r.lines.next().is_some() { return Err(bad(r.line + 1, "unexpected trailing header record")); }

        // The shared canonical decoder reserves from declared item counts.
        // Check ALL such records first, including malformed or trailing ones.
        let lines = schedule_text.lines().count();
        for line in schedule_text.lines() {
            if let Some(value) = line.strip_prefix("tracks\t") {
                if value.parse::<usize>().ok() != Some(1) {
                    return Err(bad(0, "exactly one canonical bow track is required"));
                }
            }
            if let Some(value) = line.strip_prefix("events\t") {
                let count = value.parse::<usize>().map_err(|_| bad(0, "invalid gesture event count"))?;
                if count > MAX_EVENTS || count > lines { return Err(bad(0, "gesture count exceeds source/16384-event budget")); }
            }
        }
        let schedule = GestureSchedule::from_canonical_bytes(schedule_text.as_bytes())
            .map_err(|e| BowedPerformanceError::Schedule(BowedScheduleError::Gesture(e)))?;
        if schedule.to_canonical_bytes() != schedule_text.as_bytes() {
            return Err(bad(0, "gesture bytes must be canonical without ignored fields or trailing records"));
        }
        let [track] = schedule.tracks() else { return Err(bad(0, "exactly one bow track is required")); };
        if track.target != (GestureTarget::BowStroke { string: 0 }) || schedule.control_rate_hz > rate {
            return Err(bad(0, "the only track must target bow string zero on a control rate no greater than mechanics"));
        }
        // Use the same last tick and tick/rate comparison as the actual sampler.
        let last_tick = (samples - 1) * u64::from(schedule.control_rate_hz) / u64::from(rate);
        let last_time = last_tick as f64 / f64::from(schedule.control_rate_hz);
        if track.events.iter().any(|event| event.time_s > last_time) {
            return Err(bad(0, "an authored bow command starts after the last observed control tick"));
        }
        let GestureValue::Bow { velocity_m_per_s, normal_force_n, station } = track.initial else {
            return Err(bad(0, "the initial control must be a complete bow value"));
        };
        let gesture = BowGesture::admit(velocity_m_per_s, normal_force_n, station)
            .map_err(|e| BowedPerformanceError::Model(BowedRunError::Gesture(e)))?;
        let body = CompactBody::from_radiator(radiator)
            .map_err(|_| bad(0, "the existing compact body rejected its area/mass/frequency/damping"))?;
        let config = BowedRunConfig { card, island: FrictionIsland::Stribeck(friction), gesture,
            steps: samples as usize, subsamples,
            termination: Termination::PlateOnePort { body: Box::new(body), ambient }, listener_m };
        let state = BowedStringState::new(&config, max_block).map_err(BowedPerformanceError::Model)?;
        let bow = ScheduledBowedRenderer::new(state, &schedule, &track.id, samples, max_work, max_controls)
            .map_err(BowedPerformanceError::Schedule)?;
        let info = BowedPerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.bowed-performance.v1", bytes),
            sample_rate_hz: rate, samples, full_scale_pa, string_modes: config.card.mode_count,
            subsamples, gesture_events: track.events.len(), compiled_controls: bow.pending_controls().len(),
        };
        let context = RenderContext::new(vec![RenderVoice::BowedString(Box::new(bow))], max_block);
        let renderer = ScheduledRenderer::new(context, vec![], 0).map_err(BowedPerformanceError::Render)?;
        renderer.validate_sample_rate(rate).map_err(BowedPerformanceError::Render)?;
        renderer.validate_sample_count(samples).map_err(BowedPerformanceError::Render)?;
        Ok(Self { info, renderer })
    }

    /// Actual source clocks and counts, unchanged by output conversion.
    #[must_use]
    pub const fn info(&self) -> BowedPerformanceInfo { self.info }
    /// Transfer retained physics and controls to the existing rendering owners.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }
}

struct Reader<'a> { lines: Lines<'a>, line: usize }
struct Row<'a> { words: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Reader<'a> {
    fn row(&mut self, key: &str) -> Result<Row<'a>, BowedPerformanceError> {
        self.line += 1;
        let mut words = self.lines.next().ok_or_else(|| bad(self.line, "missing header record"))?.split_ascii_whitespace();
        if words.next() != Some(key) { return Err(bad(self.line, "unexpected header record")); }
        Ok(Row { words, line: self.line })
    }
}
impl Row<'_> {
    fn parse<T: FromStr>(&mut self) -> Result<T, BowedPerformanceError> {
        self.words.next().and_then(|s| s.parse().ok()).ok_or_else(|| bad(self.line, "missing or malformed field"))
    }
    fn scalar(&mut self) -> Result<f64, BowedPerformanceError> {
        let x: f64 = self.parse()?;
        if !x.is_finite() { return Err(bad(self.line, "physical scalars must be finite")); }
        Ok(x)
    }
    fn count(&mut self, cap: usize) -> Result<usize, BowedPerformanceError> {
        let n = self.parse()?;
        if n > cap { return Err(bad(self.line, "count exceeds the input budget")); }
        Ok(n)
    }
    fn finish(mut self) -> Result<(), BowedPerformanceError> {
        if self.words.next().is_some() { return Err(bad(self.line, "unexpected extra field")); }
        Ok(())
    }
}

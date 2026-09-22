//! Authored reed/duct performances over the existing physical and gesture owners.
//!
//! The file supplies primitive SI parameters, a static segmented duct, moist air,
//! and one canonical blowing-pressure track. It does not infer material identity
//! or turn a pitch into a medium property. Pressure commands replace the fixture
//! attack envelope; all retained reed and bore state survives release and resume.
//! Output retains `ReedBoreVoice`'s bore-pressure plus compact-jet observation,
//! NOT a calibrated exterior microphone or a two-way radiation-load model.

use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_blake3::{ContentHash, hash_domain};
use fs_duct::{Duct, HoleState, Segment, Termination};
use fs_material::gas::GasState;
use fs_scenario::{BeatingReed, gesture::{GestureSchedule, GestureTarget}};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{ReedBoreVoice, RenderContext, RenderError, RenderVoice};
use crate::thin_plate::PlateBank;
use super::{GestureCompileError, PressureGestureBinding, ScheduledRenderer, compile_pressure_gestures};

/// Versioned header; the schedule retains the shared canonical wire format.
pub const REED_PERFORMANCE_SCHEMA: &str = "frankensim-reed-performance-v1";
/// Bound input bytes before reading or invoking count-driven decoders.
pub const MAX_REED_PERFORMANCE_BYTES: usize = 1024 * 1024;
const MAX_SEGMENTS: usize = 64;
const MAX_EVENTS: usize = 16_384;
const MAX_CONTROLS: usize = 262_144;
const MAX_WORK: u64 = 16_777_216;

/// Syntax, model and scheduling refusals retain their owning diagnostics.
#[derive(Debug)]
pub enum ReedPerformanceError {
    /// Invalid source record, with a one-based header line (zero for the file).
    Input {
        /// Header line, or zero for whole-file/schedule admission.
        line: usize,
        /// Failed condition.
        what: &'static str,
    },
    /// The existing reed, TMM or gas realization refused.
    Model(AcousticRealizeError),
    /// The existing pressure compiler refused.
    Schedule(GestureCompileError),
    /// The existing common renderer refused.
    Render(RenderError),
}
impl core::fmt::Display for ReedPerformanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "reed input line {line}: {what}"),
            Self::Model(e) => write!(f, "reed model: {e}"),
            Self::Schedule(e) => write!(f, "reed input: {e}"),
            Self::Render(e) => write!(f, "reed render: {e}"),
        }
    }
}
impl std::error::Error for ReedPerformanceError {}
fn bad(line: usize, what: &'static str) -> ReedPerformanceError {
    ReedPerformanceError::Input { line, what }
}

/// Source facts; neither a sourced-material receipt nor a fidelity certificate.
#[derive(Clone, Copy, Debug)]
pub struct ReedPerformanceInfo {
    /// Domain-separated hash of the complete exact input bytes.
    pub input_hash: ContentHash,
    /// Actual solver clock [Hz], independent of the observer's output clock.
    pub sample_rate_hz: u32,
    /// Finite half-open physical window [0, samples).
    pub samples: u64,
    /// Declared PCM scale [Pa]; never a per-part ensemble gain.
    pub full_scale_pa: f64,
    /// Static axial and side-branch records in their authored order.
    pub segments: usize,
    /// Side branches included in the same TMM reflectance.
    pub tone_holes: usize,
    /// Whether the source has nonzero reed mass and invokes massive mechanics.
    pub massive_reed: bool,
    /// Authored pressure commands, excluding the initial held value.
    pub gesture_events: usize,
    /// Compiled sample-addressed controls, including the initial value.
    pub compiled_controls: usize,
}

/// Finite pressure producer retaining the physical solver and schedule together.
///
/// It implements the common pressure interface directly, so decimation and
/// heterogeneous ensembles retain its horizon without extracting an unbounded
/// raw renderer. Admission never advances a sample. A late physics failure
/// inherits the existing renderer's poison-on-failure contract.
pub struct ReedPerformance {
    info: ReedPerformanceInfo,
    renderer: ScheduledRenderer,
}
impl ReedPerformance {
    /// Decode and admit all records before building or advancing the runtime.
    ///
    /// The schedule owns the entire pressure history, so the fixture's extra
    /// attack envelope is disabled. The six reed fields are opening, width,
    /// closing pressure, mass, stiffness and damping ratio. Zero mass/stiffness
    /// retain the owner's quasistatic/derived-stiffness limits. Static tone-hole
    /// openings are explicit fractions; they are never silently clamped here.
    ///
    /// # Errors
    /// Malformed/oversized input, unsupported or unobserved gestures, exhausted
    /// compilation budgets, or the existing physical owner's domain refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, ReedPerformanceError> {
        if bytes.len() > MAX_REED_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "source exceeds 1 MiB or callback capacity is outside 1..=65536"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "source must be UTF-8"))?;
        let (header, schedule_text) = text.split_once("\nschedule\n")
            .ok_or_else(|| bad(0, "expected an LF-delimited schedule record and canonical gesture bytes"))?;
        let mut r = Reader { lines: header.lines(), line: 0 };
        r.row(REED_PERFORMANCE_SCHEMA)?.finish()?;
        let mut row = r.row("audio")?;
        let rate: u32 = row.parse()?;
        let samples: u64 = row.parse()?;
        let full_scale_pa = row.positive()?;
        row.finish()?;
        if !(1..=192_000).contains(&rate) || samples == 0 || samples > 28_800_000
            || samples > u64::from(rate) * 600
        { return Err(bad(r.line, "audio requires a rate in 1..=192000 and at most 600 seconds/28800000 positive samples")); }
        let mut row = r.row("ambient")?;
        let air = GasState::try_new_moist_air(row.scalar()?, row.scalar()?, row.scalar()?)
            .map_err(|_| bad(r.line, "ambient is outside the shared moist-air model"))?;
        row.finish()?;
        let mut row = r.row("reed")?;
        let reed = BeatingReed {
            rest_opening_m: row.positive()?, width_m: row.positive()?, closing_pressure_pa: row.positive()?,
            mass_kg: row.scalar()?, stiffness_n_m: row.scalar()?, damping_ratio: row.scalar()?,
            blowing_pressure_pa: 0.0, attack_s: 0.0,
        };
        row.finish()?;
        if !crate::reed_bore::reed_parameters_valid(reed) {
            return Err(bad(r.line, "reed mass, stiffness and damping must be finite and nonnegative"));
        }
        let mut row = r.row("listener_m")?;
        let listener_m = row.positive()?; row.finish()?;
        let mut row = r.row("termination")?;
        let termination = match row.word()? {
            "closed" => Termination::Closed,
            "ideal-open" => Termination::IdealOpen,
            "unflanged" => Termination::UnflangedOpen,
            "flanged" => Termination::FlangedOpen,
            _ => return Err(bad(r.line, "unsupported static duct termination")),
        };
        row.finish()?;
        let mut row = r.row("segments")?;
        let count = row.count(MAX_SEGMENTS)?; row.finish()?;
        if count == 0 { return Err(bad(r.line, "at least one axial segment is required")); }
        let mut segments = Vec::with_capacity(count);
        let mut holes = 0;
        for _ in 0..count {
            let mut row = r.next()?;
            let segment = match row.word()? {
                "cylinder" => Segment::Cylinder { radius: row.positive()?, length: row.positive()? },
                "cone" => Segment::Cone {
                    inlet_radius: row.positive()?, outlet_radius: row.positive()?, length: row.positive()?,
                },
                "hole" => {
                    let hole_radius = row.positive()?;
                    let chimney_height = row.positive()?;
                    let bore_radius = row.positive()?;
                    let sigma = row.scalar()?;
                    if hole_radius >= bore_radius || !(0.0..=1.0).contains(&sigma) {
                        return Err(bad(r.line, "hole radius must be below bore radius and opening must be in [0,1]"));
                    }
                    holes += 1;
                    let state = if sigma == 0.0 { HoleState::Closed }
                        else if sigma == 1.0 { HoleState::Open } else { HoleState::Vent(sigma) };
                    Segment::ToneHole { hole_radius, chimney_height, bore_radius, state }
                }
                _ => return Err(bad(r.line, "expected cylinder, cone or hole geometry")),
            };
            row.finish()?;
            segments.push(segment);
        }
        if matches!(segments.first(), Some(Segment::ToneHole { .. }))
            || matches!(segments.last(), Some(Segment::ToneHole { .. }))
        { return Err(bad(r.line, "the inlet and termination must have explicit axial segments")); }
        let mut row = r.row("compile_limits")?;
        let max_work: u64 = row.parse()?;
        let max_controls = row.count(MAX_CONTROLS)?; row.finish()?;
        if max_work > MAX_WORK { return Err(bad(r.line, "compilation exceeds 16777216 source visits")); }
        if r.lines.next().is_some() { return Err(bad(r.line + 1, "unexpected trailing header record")); }
        let schedule = decode_schedule(schedule_text)?;
        let [track] = schedule.tracks() else { return Err(bad(0, "exactly one pressure track is required")); };
        if track.target != GestureTarget::BlowingPressure || schedule.control_rate_hz == 0
            || schedule.control_rate_hz > rate {
            return Err(bad(0, "only blowing pressure is supported, at a control rate no greater than mechanics"));
        }
        let last_tick = (samples - 1) * u64::from(schedule.control_rate_hz) / u64::from(rate);
        // One track emits at most one control per tick. Bound the compiler
        // before it can allocate, including long constant-valued schedules.
        if last_tick >= MAX_CONTROLS as u64 {
            return Err(bad(0, "pressure compilation exceeds the 262144 control-tick budget"));
        }
        let last_time = last_tick as f64 / f64::from(schedule.control_rate_hz);
        if track.events.iter().any(|event| event.time_s > last_time) {
            return Err(bad(0, "an authored pressure command starts after the last observed control tick"));
        }
        let controls = compile_pressure_gestures(&schedule, &[PressureGestureBinding {
            track: track.id.clone(), voice: 0,
        }], rate, samples, max_work).map_err(ReedPerformanceError::Schedule)?;
        if controls.len() > max_controls { return Err(bad(0, "compiled pressure controls exceed the declared storage budget")); }
        let info = ReedPerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.reed-performance.v1", bytes),
            sample_rate_hz: rate, samples, full_scale_pa, segments: count, tone_holes: holes,
            massive_reed: reed.mass_kg > 0.0, gesture_events: track.events.len(), compiled_controls: controls.len(),
        };
        let voice = ReedBoreVoice::new(&Duct { segments }, &air, reed, termination,
            PlateBank::default(), listener_m, rate, samples as usize, None)
            .map_err(ReedPerformanceError::Model)?;
        let context = RenderContext::new(vec![RenderVoice::ReedBore(voice)], max_block);
        let renderer = ScheduledRenderer::new(context, controls, max_controls).map_err(ReedPerformanceError::Render)?;
        renderer.validate_sample_rate(rate).map_err(ReedPerformanceError::Render)?;
        renderer.validate_sample_count(samples).map_err(ReedPerformanceError::Render)?;
        Ok(Self { info, renderer })
    }

    /// Original source facts, unchanged by output-clock conversion.
    #[must_use]
    pub const fn info(&self) -> ReedPerformanceInfo { self.info }
    /// Read-only access to sample-addressed controls and retained state clocks.
    #[must_use]
    pub const fn renderer(&self) -> &ScheduledRenderer { &self.renderer }
    /// Remaining samples on the mechanical clock.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.info.samples - self.renderer.samples_rendered() }
}
impl PressureRenderer for ReedPerformance {
    fn samples_rendered(&self) -> u64 { self.renderer.samples_rendered() }
    fn max_block_len(&self) -> usize { self.renderer.context().max_block_len() }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> { self.renderer.validate_sample_rate(rate) }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.renderer.validate_sample_count(samples)?;
        if samples > self.remaining_samples() {
            return Err(RenderError::Sizing { what: "requested pressure window exceeds the reed performance horizon" });
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        self.renderer.block(output)
    }
}

fn decode_schedule(text: &str) -> Result<GestureSchedule, ReedPerformanceError> {
    // The shared decoder reserves from declared counts, including trailing ones.
    let lines = text.lines().count();
    for line in text.lines() {
        if let Some(value) = line.strip_prefix("tracks\t") {
            if value.parse::<usize>().ok() != Some(1) { return Err(bad(0, "exactly one canonical pressure track is required")); }
        }
        if let Some(value) = line.strip_prefix("events\t") {
            let count = value.parse::<usize>().map_err(|_| bad(0, "invalid gesture event count"))?;
            if count > MAX_EVENTS || count > lines { return Err(bad(0, "gesture count exceeds source/16384-event budget")); }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(text.as_bytes())
        .map_err(|e| ReedPerformanceError::Schedule(GestureCompileError::Gesture(e)))?;
    if schedule.to_canonical_bytes() != text.as_bytes() {
        return Err(bad(0, "gesture bytes must be canonical without ignored fields or trailing records"));
    }
    Ok(schedule)
}
struct Reader<'a> { lines: Lines<'a>, line: usize }
struct Row<'a> { words: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Reader<'a> {
    fn next(&mut self) -> Result<Row<'a>, ReedPerformanceError> {
        self.line += 1;
        let words = self.lines.next().ok_or_else(|| bad(self.line, "missing header record"))?.split_ascii_whitespace();
        Ok(Row { words, line: self.line })
    }
    fn row(&mut self, key: &str) -> Result<Row<'a>, ReedPerformanceError> {
        let mut row = self.next()?;
        if row.word()? != key { return Err(bad(self.line, "unexpected header record")); }
        Ok(row)
    }
}
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, ReedPerformanceError> {
        self.words.next().ok_or_else(|| bad(self.line, "missing field"))
    }
    fn parse<T: FromStr>(&mut self) -> Result<T, ReedPerformanceError> {
        self.word()?.parse().map_err(|_| bad(self.line, "malformed field"))
    }
    fn scalar(&mut self) -> Result<f64, ReedPerformanceError> {
        let x: f64 = self.parse()?;
        if !x.is_finite() { return Err(bad(self.line, "physical scalars must be finite")); }
        Ok(x)
    }
    fn positive(&mut self) -> Result<f64, ReedPerformanceError> {
        let x = self.scalar()?;
        if x <= 0.0 { return Err(bad(self.line, "physical dimension/scale must be positive")); }
        Ok(x)
    }
    fn count(&mut self, cap: usize) -> Result<usize, ReedPerformanceError> {
        let n = self.parse()?;
        if n > cap { return Err(bad(self.line, "count exceeds the input budget")); }
        Ok(n)
    }
    fn finish(mut self) -> Result<(), ReedPerformanceError> {
        if self.words.next().is_some() { return Err(bad(self.line, "unexpected extra field")); }
        Ok(())
    }
}

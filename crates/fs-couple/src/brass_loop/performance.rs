//! A lip-valve performance file: a measured lip, a held bore, a pressure
//! phrase and a physical outlet microphone.
//!
//! Nothing here names an instrument, a note or a pitch. The file supplies
//! SI primitives and their provenance; the note that sounds is the lock the
//! loop finds between the lip and the bore ([`super::BrassVoice`]). The
//! radiated sound is the retarded derivative of the bore's MOUTH VOLUME FLOW
//! at a baffled point, never the mouthpiece pressure:
//!
//! `p+` (wave sent into the bore) -> [`mouth_flow_line`] (the transfer-matrix
//! duct's `U_mouth / p+`) -> [`BaffledPressure`] (Rayleigh integral of the
//! actual outlet disk).
//!
//! ```text
//! frankensim-lip-performance-v1
//! # comment lines and blank lines are ignored
//! audio <rate_hz> <samples>                  output clock (also the mechanics clock of the lines)
//! full_scale_pa <pa>
//! ambient moist_air <temperature_k> <pressure_pa> <relative_humidity>
//! ambient mixture <temperature_k> <pressure_pa> <mole|mass> <count> <species fraction>...
//! ambient thermochemical <temperature_k> <pressure_pa> <mole|mass> <count> <id fraction>...
//! lip_response <f_hz> <q> <inverse_mu_m2_per_kg> <width_m> <rest_gap_m> <face_area_m2>
//! lip <mass_kg> <stiffness_n_m> <damping_n_s_m> <width_m> <rest_gap_m> <face_area_m2>
//! lip_source <free text to end of line>      required provenance sentence
//! cup_volume_m3 <volume>
//! termination <flanged|unflanged|ideal_open>
//! modes <n_modes> <extra_slices>
//! mouth_flow_tail_tol <relative_energy>      optional; default 1e-10 (the line ladder's own)
//! cylinder <radius_m> <length_m>             bore records, mouthpiece to bell
//! cone <inlet_radius_m> <outlet_radius_m> <length_m>
//! bessel <inlet_radius_m> <outlet_radius_m> <length_m> <flare> <pieces>
//! bore_source <free text to end of line>     required provenance sentence
//! bore nonlinear <max_cell_m> <mechanics_per_output> <gas_steps_per_mechanics> <wall_poles> <wall_band_lo_hz> <wall_band_hi_hz> <wall_tolerance>
//! microphone <x_m> <y_m> <z_m> <radial_rings> <angular_points> <max_frequency_hz>
//! schedule
//! <canonical GestureSchedule bytes: blowing-pressure, optional pre-stress>
//! ```
//!
//! Exactly one of `lip_response` and `lip` is given. `lip_response` carries
//! the quantities lip-reed measurements publish: the outward-striking
//! resonance `f`, its quality factor `Q`, and `1/mu`, the pressure-collecting
//! area per effective mass (`mu = m / A`). The island follows
//! `m = A / (1/mu)`, `k = m (2 pi f)^2`, `c = m (2 pi f) / Q`; the face area
//! `A` cancels in every motion and only fixes the force normalization.
//!
//! A `bessel` record is a Bessel horn `r(s) = b (u - s)^(-flare)` through the
//! stated end radii over the stated length, expanded into `pieces` cones whose
//! end radii lie exactly on that curve (the piecewise-conic image every
//! transfer-matrix brass model uses; `pieces` is the refinement axis).
//!
//! The optional `bore nonlinear` record replaces the lines with the shared
//! nonlinear gas bore ([`super::nonlinear`], the `nonlinear-brass` feature):
//! the ambient gas's own NASA-9 laws evolve in cells no longer than
//! `max_cell_m`, with Zwikker-Kosten wall memory admitted on the stated band
//! and tolerance and a baffled-piston bell. It needs an `ambient
//! thermochemical` record and its packs. The mechanics clock is
//! `mechanics_per_output` (1..=16) times the output clock, each mechanics
//! sample taking `gas_steps_per_mechanics` fixed gas steps; the microphone
//! receives every mechanics sample's bell volume flow.
//!
//! Limits, not hidden: one held bore (no valve or slide motion), a plane-mode
//! mouth-flow transfer, a baffled-disk radiation integral, and the fs-duct
//! flanged-radiation ceiling (`ka <= 20` at the Nyquist frequency of the
//! declared clock at the mouth radius), which sets the highest usable clock.

use fs_blake3::{ContentHash, hash_domain};
use fs_duct::{Duct, LossModel, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{
    GasRecord, ThermochemicalSources, read_frozen_thermochemical_model, read_gas_record,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_math::det;
use fs_scenario::gesture::{GestureSchedule, GestureTarget};
use fs_vfit::discretize::DelayedFilter;

use super::{BrassBlockDiag, BrassControl, BrassVoice, LipIsland};
use crate::acoustic_realize::AcousticRealizeError;
use crate::driving_point::LineLadder;
use crate::driving_point::outlet::mouth_flow_line;
use crate::mm_line::{MmComboReceipt, MmLineConfig, MmLoad};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;

/// First line of every lip performance file.
pub const LIP_PERFORMANCE_SCHEMA: &str = "frankensim-lip-performance-v1";
/// Largest accepted file.
pub const MAX_LIP_PERFORMANCE_BYTES: usize = 1 << 20;
const MAX_SEGMENTS: usize = 512;
const MAX_BESSEL_PIECES: usize = 256;
const MAX_MODES: u64 = 16;
const MAX_CONTROL_TICKS: u64 = 262_144;
const MAX_EVENTS: usize = 16_384;
const MAX_SECONDS: u64 = 600;
const MAX_BLOCK: usize = 65_536;
const MIN_RATE_HZ: u32 = 8_000;
const MAX_RATE_HZ: u32 = 192_000;

/// The lip's lumped numbers before provenance is attached.
#[derive(Debug, Clone, Copy)]
struct LipCore {
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    width_m: f64,
    rest_gap_m: f64,
    face_area_m2: f64,
}

/// A `bore nonlinear` record: how the shared nonlinear gas bore is cut,
/// clocked and walled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearBoreRecord {
    /// Largest axial gas cell [m].
    pub maximum_cell_m: f64,
    /// Mechanics (microphone) samples per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Fixed gas steps per mechanics sample.
    pub gas_steps_per_mechanics: u32,
    /// Zwikker-Kosten diffusive poles per wall field.
    pub wall_poles: usize,
    /// Band on which the wall memory is admitted [Hz].
    pub wall_band_hz: [f64; 2],
    /// Largest sampled relative error of either wall operator on that band.
    pub wall_tolerance: f64,
}

const NONLINEAR_FEATURE: &str = "bore nonlinear needs fs-couple's nonlinear-brass feature (build music_render with --features nonlinear-brass)";

/// A parsed, validated, not yet realized performance.
#[derive(Debug, Clone)]
pub struct LipDescription {
    input_hash: ContentHash,
    rate_hz: u32,
    samples: u64,
    full_scale_pa: f64,
    gas: GasRecord,
    lip: LipIsland,
    lip_source: String,
    bore_source: String,
    cup_volume_m3: f64,
    termination: Termination,
    modes: usize,
    extra_slices: usize,
    duct: Duct,
    microphone: CircularOutletReceiver,
    mouth_flow_tail_tol: f64,
    schedule: GestureSchedule,
    blow_track: String,
    tension_track: Option<String>,
    nonlinear: Option<(NonlinearBoreRecord, ThermochemicalGasModel)>,
}

struct Fields<'a> {
    words: core::str::SplitAsciiWhitespace<'a>,
    line: usize,
    key: &'a str,
}

impl<'a> Fields<'a> {
    fn err<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("line {}: {}: {what}", self.line, self.key))
    }

    fn word(&mut self) -> Result<&'a str, String> {
        match self.words.next() {
            Some(word) => Ok(word),
            None => self.err("missing field"),
        }
    }

    fn parse<T: core::str::FromStr>(&mut self) -> Result<T, String> {
        let word = self.word()?;
        match word.parse::<T>() {
            Ok(value) => Ok(value),
            Err(_) => self.err("malformed field"),
        }
    }

    fn number(&mut self) -> Result<f64, String> {
        let x: f64 = self.parse()?;
        if x.is_finite() {
            Ok(x)
        } else {
            self.err("physical scalars must be finite")
        }
    }

    fn positive(&mut self) -> Result<f64, String> {
        let x = self.number()?;
        if x > 0.0 {
            Ok(x)
        } else {
            self.err("expected a positive number")
        }
    }

    fn count(&mut self, max: u64) -> Result<u64, String> {
        let n: u64 = self.parse()?;
        if (1..=max).contains(&n) {
            Ok(n)
        } else {
            self.err("count outside its admitted range")
        }
    }

    fn finish(mut self) -> Result<(), String> {
        if self.words.next().is_some() {
            return self.err("unexpected extra field");
        }
        Ok(())
    }
}

fn once<T>(slot: &mut Option<T>, value: T, f: &Fields<'_>) -> Result<(), String> {
    if slot.is_some() {
        return f.err("may only appear once");
    }
    *slot = Some(value);
    Ok(())
}

fn text_of(f: &Fields<'_>, rest: &str) -> Result<String, String> {
    if rest.is_empty() {
        return f.err("needs a provenance sentence");
    }
    Ok(rest.to_string())
}

/// Expand a Bessel horn `r(s) = r_in (u / (u - s))^flare` through `r_out`
/// at `s = length` into `pieces` cones with end radii on that curve.
fn bessel_horn(
    r_in: f64,
    r_out: f64,
    length: f64,
    flare: f64,
    pieces: usize,
    out: &mut Vec<Segment>,
) -> Result<(), String> {
    if !(r_out > r_in) || !(flare > 0.0) || !(1..=MAX_BESSEL_PIECES).contains(&pieces) {
        return Err(
            "a Bessel horn needs a growing radius, a positive flare and 1..=256 pieces".to_string(),
        );
    }
    // r_out / r_in = (u / (u - length))^flare  =>  u / (u - length) = q.
    let q = det::pow(r_out / r_in, 1.0 / flare);
    if !(q.is_finite() && q > 1.0) {
        return Err("Bessel horn end radii and flare admit no horn".to_string());
    }
    let u = q * length / (q - 1.0);
    let step = length / pieces as f64;
    let mut previous = r_in;
    for k in 1..=pieces {
        let radius = if k == pieces {
            r_out
        } else {
            r_in * det::pow(u / (u - step * k as f64), flare)
        };
        if !(radius.is_finite() && radius > previous) {
            return Err("Bessel horn pieces must have strictly growing radius".to_string());
        }
        out.push(Segment::Cone {
            inlet_radius: previous,
            outlet_radius: radius,
            length: step,
        });
        previous = radius;
    }
    Ok(())
}

fn decode_schedule(text: &str) -> Result<GestureSchedule, String> {
    // The shared decoder reserves from declared counts, including trailing ones.
    let lines = text.lines().count();
    for line in text.lines() {
        if let Some(value) = line.strip_prefix("tracks\t") {
            if !matches!(value.parse::<usize>(), Ok(1 | 2)) {
                return Err("one blowing-pressure track and an optional pre-stress track".into());
            }
        }
        if let Some(value) = line.strip_prefix("events\t") {
            let count = value
                .parse::<usize>()
                .map_err(|_| "invalid gesture event count".to_string())?;
            if count > MAX_EVENTS || count > lines {
                return Err("gesture count exceeds the source/16384-event budget".into());
            }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(text.as_bytes())
        .map_err(|e| format!("gesture schedule refused: {e:?}"))?;
    if schedule.to_canonical_bytes() != text.as_bytes() {
        return Err(
            "gesture bytes must be canonical without ignored fields or trailing records".into(),
        );
    }
    Ok(schedule)
}

/// The control tick whose value holds at mechanical sample `sample`: tick `k`
/// applies at `ceil(k * rate / control_rate)`, the clock the pressure
/// scheduler of the reed files uses.
fn tick_at(sample: u64, control_rate_hz: u64, rate_hz: u64) -> u64 {
    sample * control_rate_hz / rate_hz
}

/// First mechanical sample of control tick `tick`.
fn tick_start(tick: u64, control_rate_hz: u64, rate_hz: u64) -> u64 {
    (tick * rate_hz).div_ceil(control_rate_hz)
}

impl LipDescription {
    /// Decode and admit every record before any physics is realized.
    ///
    /// # Errors
    /// A named refusal for malformed, repeated, missing, unphysical or
    /// out-of-budget records.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        Self::from_bytes_with_gas_sources(bytes, None)
    }

    /// Admit the same shared gas records as ducts, reeds and flue sources.
    /// NASA-9 records require explicitly associated species/transport packs;
    /// their receipts stay attached to the exact state used by every consumer.
    /// The gas is frozen throughout this performance, including the exterior.
    ///
    /// # Errors
    /// Ordinary file refusals plus missing, ambiguous or out-of-domain gas
    /// sources. No material properties are inferred from a species identifier.
    #[allow(clippy::too_many_lines)] // one record switch
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, String> {
        if bytes.len() > MAX_LIP_PERFORMANCE_BYTES {
            return Err("source exceeds 1 MiB".into());
        }
        let text = core::str::from_utf8(bytes).map_err(|_| "source must be UTF-8".to_string())?;
        let (header, schedule_text) = text.split_once("\nschedule\n").ok_or_else(|| {
            "expected an LF-delimited schedule record and canonical gesture bytes".to_string()
        })?;
        let mut lines = header.lines().enumerate();
        match lines.next() {
            Some((_, first)) if first.trim_end() == LIP_PERFORMANCE_SCHEMA => {}
            _ => return Err(format!("expected the {LIP_PERFORMANCE_SCHEMA} header")),
        }
        let mut audio: Option<(u32, u64)> = None;
        let mut full_scale: Option<f64> = None;
        let mut gas: Option<GasRecord> = None;
        let mut lip_core: Option<LipCore> = None;
        let mut lip_source: Option<String> = None;
        let mut bore_source: Option<String> = None;
        let mut cup: Option<f64> = None;
        let mut termination: Option<Termination> = None;
        let mut modes: Option<(usize, usize)> = None;
        let mut microphone: Option<CircularOutletReceiver> = None;
        let mut mouth_flow_tail_tol: Option<f64> = None;
        let mut nonlinear: Option<NonlinearBoreRecord> = None;
        let mut ambient_words = None;
        let mut segments: Vec<Segment> = Vec::new();
        for (index, raw) in lines {
            let line = index + 1;
            let trimmed = raw.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let (key, rest) = match trimmed.split_once(char::is_whitespace) {
                Some((key, rest)) => (key, rest.trim()),
                None => (trimmed, ""),
            };
            let mut f = Fields {
                words: rest.split_ascii_whitespace(),
                line,
                key,
            };
            match key {
                "audio" => {
                    let rate: u32 = f.parse()?;
                    let samples: u64 = f.parse()?;
                    if !(MIN_RATE_HZ..=MAX_RATE_HZ).contains(&rate)
                        || samples == 0
                        || samples > u64::from(rate) * MAX_SECONDS
                    {
                        return f.err("needs a rate in 8000..=192000 Hz and 1..=600 s of samples");
                    }
                    once(&mut audio, (rate, samples), &f)?;
                    f.finish()?;
                }
                "full_scale_pa" => {
                    let value = f.positive()?;
                    once(&mut full_scale, value, &f)?;
                    f.finish()?;
                }
                "ambient" => {
                    if f.words.clone().next() == Some("moist_air") {
                        f.words.next();
                    }
                    ambient_words = Some(f.words.clone());
                    let record = read_gas_record(&mut f.words, None, sources)
                        .map_err(|e| format!("line {line}: ambient gas refused: {e}"))?;
                    once(&mut gas, record, &f)?;
                    f.finish()?;
                }
                "lip_response" => {
                    let (freq, q, inverse_mu) = (f.positive()?, f.positive()?, f.positive()?);
                    let (width, gap, area) = (f.positive()?, f.positive()?, f.positive()?);
                    let mass = area / inverse_mu;
                    let omega = core::f64::consts::TAU * freq;
                    let core = LipCore {
                        mass_kg: mass,
                        stiffness_n_m: mass * omega * omega,
                        damping_n_s_m: mass * omega / q,
                        width_m: width,
                        rest_gap_m: gap,
                        face_area_m2: area,
                    };
                    once(&mut lip_core, core, &f)?;
                    f.finish()?;
                }
                "lip" => {
                    let (mass, stiffness) = (f.positive()?, f.positive()?);
                    let damping = f.number()?;
                    let (width, gap, area) = (f.positive()?, f.positive()?, f.positive()?);
                    if damping < 0.0 {
                        return f.err("damping must be nonnegative");
                    }
                    let core = LipCore {
                        mass_kg: mass,
                        stiffness_n_m: stiffness,
                        damping_n_s_m: damping,
                        width_m: width,
                        rest_gap_m: gap,
                        face_area_m2: area,
                    };
                    once(&mut lip_core, core, &f)?;
                    f.finish()?;
                }
                "lip_source" => {
                    let value = text_of(&f, rest)?;
                    once(&mut lip_source, value, &f)?;
                }
                "bore_source" => {
                    let value = text_of(&f, rest)?;
                    once(&mut bore_source, value, &f)?;
                }
                "cup_volume_m3" => {
                    let value = f.positive()?;
                    once(&mut cup, value, &f)?;
                    f.finish()?;
                }
                "termination" => {
                    let kind = match f.word()? {
                        "flanged" => Termination::FlangedOpen,
                        "unflanged" => Termination::UnflangedOpen,
                        "ideal_open" => Termination::IdealOpen,
                        _ => return f.err("expected flanged, unflanged or ideal_open"),
                    };
                    once(&mut termination, kind, &f)?;
                    f.finish()?;
                }
                "modes" => {
                    let n_modes = f.count(MAX_MODES)? as usize;
                    let extra = f.count(MAX_MODES)? as usize;
                    once(&mut modes, (n_modes, extra), &f)?;
                    f.finish()?;
                }
                "cylinder" => {
                    let (radius, length) = (f.positive()?, f.positive()?);
                    segments.push(Segment::Cylinder { radius, length });
                    f.finish()?;
                }
                "cone" => {
                    let (inlet, outlet, length) = (f.positive()?, f.positive()?, f.positive()?);
                    segments.push(if inlet == outlet {
                        Segment::Cylinder {
                            radius: inlet,
                            length,
                        }
                    } else {
                        Segment::Cone {
                            inlet_radius: inlet,
                            outlet_radius: outlet,
                            length,
                        }
                    });
                    f.finish()?;
                }
                "bessel" => {
                    let (inlet, outlet, length) = (f.positive()?, f.positive()?, f.positive()?);
                    let flare = f.positive()?;
                    let pieces = f.count(MAX_BESSEL_PIECES as u64)? as usize;
                    bessel_horn(inlet, outlet, length, flare, pieces, &mut segments)
                        .map_err(|e| format!("line {line}: bessel: {e}"))?;
                    f.finish()?;
                }
                "bore" => {
                    if f.word()? != "nonlinear" {
                        return f.err("expected nonlinear");
                    }
                    let maximum_cell_m = f.positive()?;
                    let mechanics_per_output = f.count(16)? as u32;
                    let gas_steps_per_mechanics = f.count(64)? as u32;
                    let wall_poles = f.count(64)? as usize;
                    let wall_band_hz = [f.positive()?, f.positive()?];
                    let wall_tolerance = f.positive()?;
                    if wall_band_hz[0] >= wall_band_hz[1] || wall_tolerance >= 1.0 {
                        return f
                            .err("needs a rising wall band and a relative wall tolerance below 1");
                    }
                    let record = NonlinearBoreRecord {
                        maximum_cell_m,
                        mechanics_per_output,
                        gas_steps_per_mechanics,
                        wall_poles,
                        wall_band_hz,
                        wall_tolerance,
                    };
                    once(&mut nonlinear, record, &f)?;
                    f.finish()?;
                }
                "mouth_flow_tail_tol" => {
                    let value = f.positive()?;
                    if value >= 1.0 {
                        return f.err("a relative tail energy must be below 1");
                    }
                    once(&mut mouth_flow_tail_tol, value, &f)?;
                    f.finish()?;
                }
                "microphone" => {
                    let position = [f.number()?, f.number()?, f.number()?];
                    let rings = f.count(128)? as usize;
                    let points = f.count(512)? as usize;
                    let band = f.positive()?;
                    if position[2] <= 0.0 {
                        return f.err("the microphone must stand in front of the baffle (z > 0)");
                    }
                    let receiver = CircularOutletReceiver {
                        position_m: position,
                        radial_rings: rings,
                        angular_points: points,
                        maximum_frequency_hz: band,
                    };
                    once(&mut microphone, receiver, &f)?;
                    f.finish()?;
                }
                _ => return Err(format!("line {line}: unknown key {key:?}")),
            }
            if segments.len() > MAX_SEGMENTS {
                return Err(format!("line {line}: bore exceeds {MAX_SEGMENTS} segments"));
            }
        }
        let missing = |key: &str| format!("missing required record {key}");
        let (rate_hz, samples) = audio.ok_or_else(|| missing("audio"))?;
        let core = lip_core.ok_or_else(|| missing("lip_response or lip"))?;
        let lip_source = lip_source.ok_or_else(|| missing("lip_source"))?;
        let bore_source = bore_source.ok_or_else(|| missing("bore_source"))?;
        if segments.is_empty() {
            return Err(missing("cylinder, cone or bessel"));
        }
        let lip = LipIsland {
            mass_kg: core.mass_kg,
            stiffness_n_m: core.stiffness_n_m,
            damping_n_s_m: core.damping_n_s_m,
            width_m: core.width_m,
            rest_gap_m: core.rest_gap_m,
            face_area_m2: core.face_area_m2,
            provenance: lip_source.clone(),
        };
        let (modes_n, extra_slices) = modes.ok_or_else(|| missing("modes"))?;
        let schedule = decode_schedule(schedule_text)?;
        let mut blow_track = None;
        let mut tension_track = None;
        for track in schedule.tracks() {
            let slot = match &track.target {
                GestureTarget::BlowingPressure => &mut blow_track,
                GestureTarget::PreStress => &mut tension_track,
                _ => {
                    return Err(
                        "a lip performance accepts only blowing-pressure and pre-stress tracks"
                            .into(),
                    );
                }
            };
            if slot.replace(track.id.clone()).is_some() {
                return Err("each control track may appear only once".into());
            }
        }
        let blow_track =
            blow_track.ok_or_else(|| "a blowing-pressure track is required".to_string())?;
        let nonlinear = match nonlinear {
            None => None,
            Some(_) if !cfg!(feature = "nonlinear-brass") => return Err(NONLINEAR_FEATURE.into()),
            // Its bell is a rigid piston in an infinite baffle, nothing else.
            Some(_) if termination != Some(Termination::FlangedOpen) => {
                return Err(
                    "bore nonlinear radiates as a baffled piston: it needs termination flanged"
                        .into(),
                );
            }
            Some(record) => {
                // The bore evolves the ambient gas's own laws, never a snapshot.
                let model = match (ambient_words, sources) {
                    (Some(mut words), Some(sources))
                        if words.clone().next() == Some("thermochemical") =>
                    {
                        read_frozen_thermochemical_model(&mut words, sources)
                            .map_err(|e| {
                                format!("bore nonlinear: ambient gas source refused: {e}")
                            })?
                            .0
                    }
                    _ => {
                        return Err(
                            "bore nonlinear evolves the ambient gas's own NASA-9 laws: it needs an ambient thermochemical record and its packs"
                                .into(),
                        );
                    }
                };
                Some((record, model))
            }
        };
        let control_rate = u64::from(schedule.control_rate_hz);
        if control_rate == 0 || control_rate > u64::from(rate_hz) {
            return Err(
                "the control rate must be positive and no greater than the mechanics clock".into(),
            );
        }
        let last_tick = tick_at(samples - 1, control_rate, u64::from(rate_hz));
        if last_tick >= MAX_CONTROL_TICKS {
            return Err("control compilation exceeds the 262144 control-tick budget".into());
        }
        let last_time = last_tick as f64 / control_rate as f64;
        if schedule
            .tracks()
            .iter()
            .any(|track| track.events.iter().any(|event| event.time_s > last_time))
        {
            return Err("an authored command starts after the last observed control tick".into());
        }
        Ok(Self {
            input_hash: hash_domain("org.frankensim.fs-couple.lip-performance.v1", bytes),
            rate_hz,
            samples,
            full_scale_pa: full_scale.ok_or_else(|| missing("full_scale_pa"))?,
            gas: gas.ok_or_else(|| missing("ambient"))?,
            lip,
            lip_source,
            bore_source,
            cup_volume_m3: cup.ok_or_else(|| missing("cup_volume_m3"))?,
            termination: termination.ok_or_else(|| missing("termination"))?,
            modes: modes_n,
            extra_slices,
            duct: Duct { segments },
            microphone: microphone.ok_or_else(|| missing("microphone"))?,
            mouth_flow_tail_tol: mouth_flow_tail_tol.unwrap_or(LineLadder::default().tail_tol),
            schedule,
            blow_track,
            tension_track,
            nonlinear,
        })
    }

    /// The `bore nonlinear` record and the ambient source model it evolves,
    /// when the file declares one.
    #[must_use]
    pub const fn nonlinear_bore(&self) -> Option<&(NonlinearBoreRecord, ThermochemicalGasModel)> {
        self.nonlinear.as_ref()
    }

    /// The exact shared state used by the lip, bore, mouth-flow and receiver.
    #[must_use]
    pub fn ambient_gas(&self) -> &GasState {
        self.gas.state()
    }

    /// Original gas source receipts, retained without constitutive conversion.
    #[must_use]
    pub const fn gas_record(&self) -> &GasRecord {
        &self.gas
    }

    /// Retained transverse modes in each section of the bore's acoustic image.
    #[must_use]
    pub const fn transverse_mode_count(&self) -> usize {
        self.modes
    }

    /// Retained source control items, including each track's initial value.
    /// Ramps remain owned and evaluated by this source on its control clock.
    #[must_use]
    pub fn control_event_count(&self) -> usize {
        self.schedule.tracks().iter().map(|track| 1 + track.events.len()).sum()
    }

    /// The passive bore (no lip, no cup shunt, no radiation observer) scanned at `frequencies_hz`: the multimodal and the
    /// plane-wave input impedance and the mouth-to-input flow ratio, from the same duct, gas and termination the render uses.
    ///
    /// # Errors
    /// The duct owner's refusal at any frequency (loss-model or radiation validity).
    pub fn bore_scan(&self, frequencies_hz: &[f64]) -> Result<BoreScan, String> {
        let gas = self.gas.state();
        let radius = self
            .duct
            .segments
            .first()
            .map(Segment::inlet_radius)
            .ok_or_else(|| "empty bore".to_string())?;
        let [x, y, z] = self.microphone.position_m;
        let mut rows = Vec::with_capacity(frequencies_hz.len());
        for &frequency_hz in frequencies_hz {
            let omega = core::f64::consts::TAU * frequency_hz;
            let multimodal = fs_duct::modal::mm_input_impedance(
                &self.duct,
                gas,
                omega,
                LossModel::Bessel,
                self.termination,
                self.modes,
                self.extra_slices,
            )
            .map_err(|e| format!("multimodal bore at {frequency_hz} Hz refused: {e}"))?;
            let scalar = fs_duct::input_impedance(
                &self.duct,
                gas,
                omega,
                LossModel::Bessel,
                self.termination,
            )
            .map_err(|e| format!("plane-wave bore at {frequency_hz} Hz refused: {e}"))?;
            rows.push(BoreResponse {
                frequency_hz,
                impedance_multimodal: [
                    multimodal.plane_impedance.re,
                    multimodal.plane_impedance.im,
                ],
                impedance_scalar: [scalar.impedance.re, scalar.impedance.im],
                mouth_over_input_flow: [scalar.u_mouth_over_u_in.re, scalar.u_mouth_over_u_in.im],
            });
        }
        Ok(BoreScan {
            input_characteristic_impedance: gas.density * gas.sound_speed
                / (core::f64::consts::PI * radius * radius),
            density_kg_m3: gas.density,
            microphone_distance_m: x.hypot(y).hypot(z),
            rows,
        })
    }

    /// Realize the bore (the lines, or the nonlinear gas bore when the file
    /// declares one) and the microphone.
    ///
    /// # Errors
    /// Any owner's refusal: the multimodal bank (passivity, ladder
    /// convergence, radiation `ka` ceiling), the mouth-flow ladder, the
    /// nonlinear bore (gas source, mesh, wall binding, clock) or the
    /// receiver's resolved-band admission.
    pub fn build(&self) -> Result<LipPerformance, String> {
        let gas = self.ambient_gas();
        let mouth_radius_m = self
            .duct
            .segments
            .last()
            .map(Segment::outlet_radius)
            .ok_or_else(|| "empty bore".to_string())?;
        let medium = RayleighMedium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        #[cfg(feature = "nonlinear-brass")]
        if let Some((record, model)) = &self.nonlinear {
            return self.build_nonlinear(*record, model, mouth_radius_m, medium);
        }
        let config = MmLineConfig {
            sample_rate_hz: self.rate_hz,
            n_modes: self.modes,
            extra_slices: self.extra_slices,
            loss: LossModel::Bessel,
        };
        let load = MmLoad::Analytic(self.termination);
        let voice = BrassVoice::new(
            core::slice::from_ref(&self.duct),
            &["held"],
            gas,
            &load,
            &config,
            self.lip.clone(),
            self.cup_volume_m3,
        )
        .map_err(|e| format!("brass voice refused: {e}"))?;
        let line = voice
            .receipts()
            .first()
            .cloned()
            .ok_or_else(|| "the bank realized no line".to_string())?;
        let (mouth, mouth_receipt) = mouth_flow_line(
            &self.duct,
            gas,
            self.termination,
            self.rate_hz,
            voice.port_impedance(),
            None,
            &LineLadder {
                tail_tol: self.mouth_flow_tail_tol,
                ..LineLadder::default()
            },
        )
        .map_err(|e| format!("mouth-flow line refused: {e}"))?;
        let receiver = BaffledPressure::circular_outlet_multirate_differentiating(
            mouth_radius_m,
            self.rate_hz,
            self.rate_hz,
            self.microphone,
            medium,
        )?;
        let bore = LipBoreInfo::Linear(LinearBoreInfo {
            line,
            mouth_flow_taps: mouth_receipt.n_fft,
            mouth_flow_rungs: mouth_receipt.rungs.len(),
            mouth_flow_tail_tol: self.mouth_flow_tail_tol,
            mouth_flow_ladder: mouth_receipt
                .rungs
                .iter()
                .map(|r| (r.n_fft, r.tail_energy, r.tap_change))
                .collect(),
        });
        Ok(self.performance(
            mouth_radius_m,
            bore,
            Engine::Linear { voice, mouth },
            receiver,
            1,
        ))
    }

    #[cfg(feature = "nonlinear-brass")]
    fn build_nonlinear(
        &self,
        record: NonlinearBoreRecord,
        model: &ThermochemicalGasModel,
        mouth_radius_m: f64,
        medium: RayleighMedium,
    ) -> Result<LipPerformance, String> {
        use super::nonlinear::{NonlinearClock, NonlinearLip};
        use fs_duct::nonlinear::bore::{WallGas, WallMemory};
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        let gas = self.ambient_gas();
        let [lowest, highest] = room_region(model.thermodynamics(), gas.temperature)?;
        let source = FrozenGas::new(
            model.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(lowest), Temperature::new(highest)),
            || false,
        )
        .map_err(|e| format!("nonlinear bore gas refused: {e}"))?;
        let clock = NonlinearClock {
            maximum_cell_m: record.maximum_cell_m,
            output_rate_hz: self.rate_hz,
            mechanics_per_output: record.mechanics_per_output,
            gas_steps_per_mechanics: record.gas_steps_per_mechanics,
            wall: WallMemory {
                poles: record.wall_poles,
                band_hz: record.wall_band_hz,
                maximum_relative_error: record.wall_tolerance,
            },
        };
        let wall_gas = WallGas {
            density_kg_m3: gas.density,
            dynamic_viscosity_pa_s: gas.dynamic_viscosity,
            gamma: gas.gamma,
            prandtl: gas.prandtl,
        };
        let lip = NonlinearLip::new(
            &self.duct,
            source,
            gas.temperature,
            gas.pressure,
            wall_gas,
            self.lip.clone(),
            self.cup_volume_m3,
            clock,
        )?;
        let receiver = BaffledPressure::circular_outlet_multirate_differentiating(
            mouth_radius_m,
            clock.mechanics_rate_hz(),
            self.rate_hz,
            self.microphone,
            medium,
        )?;
        let wall = lip
            .wall_receipt()
            .ok_or_else(|| "the nonlinear bore bound no wall memory".to_string())?;
        let bore = LipBoreInfo::Nonlinear(NonlinearBoreInfo {
            record,
            mechanics_rate_hz: clock.mechanics_rate_hz(),
            gas_step_s: clock.gas_step_s(),
            cells: lip.bore().centers_m().len(),
            radiation_taps: lip.bore().radiation_taps(),
            wall_distinct_radii: wall.distinct_radii,
            wall_maximum_arms: wall.maximum_arms,
            wall_maximum_complex_error: wall.maximum_complex_error,
            wall_maximum_loss_error: wall.maximum_loss_error,
        });
        Ok(self.performance(
            mouth_radius_m,
            bore,
            Engine::Nonlinear(Box::new(lip)),
            receiver,
            record.mechanics_per_output as usize,
        ))
    }

    fn performance(
        &self,
        mouth_radius_m: f64,
        bore: LipBoreInfo,
        engine: Engine,
        receiver: BaffledPressure,
        ratio: usize,
    ) -> LipPerformance {
        let gas = self.ambient_gas();
        let info = LipInfo {
            input_hash: self.input_hash,
            sample_rate_hz: self.rate_hz,
            samples: self.samples,
            full_scale_pa: self.full_scale_pa,
            segments: self.duct.segments.len(),
            bore_length_m: self.duct.segments.iter().map(Segment::axial_length).sum(),
            mouth_radius_m,
            lip: self.lip.clone(),
            lip_source: self.lip_source.clone(),
            bore_source: self.bore_source.clone(),
            cup_volume_m3: self.cup_volume_m3,
            gas_temperature_k: gas.temperature,
            sound_speed_m_s: gas.sound_speed,
            microphone: self.microphone,
            bore,
            control_rate_hz: self.schedule.control_rate_hz,
            has_lip_tension_track: self.tension_track.is_some(),
        };
        LipPerformance {
            info,
            engine,
            receiver,
            schedule: self.schedule.clone(),
            blow_track: self.blow_track.clone(),
            tension_track: self.tension_track.clone(),
            completed: 0,
            applied_pressure: f64::NAN,
            applied_tension: f64::NAN,
            poisoned: false,
            mouthpiece: vec![0.0; MAX_BLOCK],
            outgoing: vec![0.0; MAX_BLOCK],
            trace: vec![0.0; MAX_BLOCK * ratio],
            ratio,
            window_len: u64::from(self.rate_hz / 10).max(1),
            acc: Vec::new(),
            tap: None,
        }
    }
}

/// The temperatures the nonlinear field may reach: for every species, the
/// NASA-9 region that holds the room. The packs state no caloric continuity
/// across their region boundaries, so a state that would cross one refuses
/// instead. A bound shared with the next region belongs to that region.
#[cfg(feature = "nonlinear-brass")]
fn room_region(
    source: &fs_thermochem::mixture::FrozenIdealGasMixtureModelV1,
    room_k: f64,
) -> Result<[f64; 2], String> {
    let mut domain = [0.0_f64, f64::INFINITY];
    for component in source.components() {
        let regions = component.regions();
        let shared = |index: usize, maximum: f64| {
            regions
                .get(index + 1)
                .is_some_and(|next| next.temperature_min().value().to_bits() == maximum.to_bits())
        };
        let (index, region) = regions
            .iter()
            .enumerate()
            .find(|(index, region)| {
                let (minimum, maximum) = (
                    region.temperature_min().value(),
                    region.temperature_max().value(),
                );
                room_k >= minimum
                    && (room_k < maximum || (!shared(*index, maximum) && room_k <= maximum))
            })
            .ok_or_else(|| {
                format!(
                    "bore nonlinear: the room ({room_k} K) lies outside every NASA-9 region of {}",
                    component.species()
                )
            })?;
        let maximum = region.temperature_max().value();
        domain[0] = domain[0].max(region.temperature_min().value());
        domain[1] = domain[1].min(if shared(index, maximum) {
            maximum.next_down()
        } else {
            maximum
        });
    }
    Ok(domain)
}

/// What a realized performance reports about itself.
#[derive(Debug, Clone)]
pub struct LipInfo {
    /// Domain-separated hash of the exact input bytes.
    pub input_hash: ContentHash,
    /// Mechanics and output clock [Hz].
    pub sample_rate_hz: u32,
    /// Declared samples.
    pub samples: u64,
    /// Pascals at PCM full scale.
    pub full_scale_pa: f64,
    /// Bore segments after Bessel expansion.
    pub segments: usize,
    /// Axial length of the bore [m].
    pub bore_length_m: f64,
    /// Radius of the radiating outlet [m].
    pub mouth_radius_m: f64,
    /// The lip as the loop uses it, with its provenance sentence.
    pub lip: LipIsland,
    /// Provenance sentence of the lip numbers.
    pub lip_source: String,
    /// Provenance sentence of the bore geometry.
    pub bore_source: String,
    /// Lumped mouthpiece-cup volume [m^3].
    pub cup_volume_m3: f64,
    /// Bore gas temperature [K].
    pub gas_temperature_k: f64,
    /// Bore sound speed [m/s].
    pub sound_speed_m_s: f64,
    /// The declared microphone.
    pub microphone: CircularOutletReceiver,
    /// How the bore was realized, with its receipts.
    pub bore: LipBoreInfo,
    /// Control clock of the pressure phrase [Hz].
    pub control_rate_hz: u32,
    /// Whether a pre-stress (lip tension) track drives the lip.
    pub has_lip_tension_track: bool,
}

/// The realized bore of a performance.
#[derive(Debug, Clone)]
pub enum LipBoreInfo {
    /// The multimodal reflectance line and the plane-mode mouth-flow line.
    Linear(LinearBoreInfo),
    /// The shared nonlinear gas bore.
    Nonlinear(NonlinearBoreInfo),
}

/// Receipts of the linear lines.
#[derive(Debug, Clone)]
pub struct LinearBoreInfo {
    /// Realization receipt of the multimodal reflectance line.
    pub line: MmComboReceipt,
    /// Taps of the mouth-flow line.
    pub mouth_flow_taps: usize,
    /// Ladder rungs the mouth-flow line traversed.
    pub mouth_flow_rungs: usize,
    /// Relative tail energy the mouth-flow ladder was asked to reach (the
    /// owner's default unless the file declares `mouth_flow_tail_tol`).
    pub mouth_flow_tail_tol: f64,
    /// Every mouth-flow rung: `(taps, tail energy, tap change vs the previous rung)`.
    pub mouth_flow_ladder: Vec<(usize, f64, f64)>,
}

/// Receipts of the nonlinear gas bore.
#[derive(Debug, Clone, Copy)]
pub struct NonlinearBoreInfo {
    /// The declared record.
    pub record: NonlinearBoreRecord,
    /// Mechanics (microphone) clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Fixed gas step [s].
    pub gas_step_s: f64,
    /// Gas cells.
    pub cells: usize,
    /// Past bell flows the piston radiation convolves.
    pub radiation_taps: usize,
    /// Distinct wall radii bound.
    pub wall_distinct_radii: usize,
    /// Largest viscous or thermal arm count in any cell.
    pub wall_maximum_arms: usize,
    /// Worst sampled complex error of the wall operators on the band.
    pub wall_maximum_complex_error: f64,
    /// Worst sampled error of their dissipative parts.
    pub wall_maximum_loss_error: f64,
}

/// The passive bore's response at one frequency, at its input plane (the mouthpiece cup rim).
#[derive(Debug, Clone, Copy)]
pub struct BoreResponse {
    /// Frequency [Hz].
    pub frequency_hz: f64,
    /// Multimodal plane-mode input impedance `[re, im]` [Pa s/m^3]: what the playing loop's reflectance is built from.
    pub impedance_multimodal: [f64; 2],
    /// Plane-wave (scalar transfer-matrix) input impedance `[re, im]` [Pa s/m^3]: what the radiation observer's mouth flow uses.
    pub impedance_scalar: [f64; 2],
    /// Mouth volume flow over input volume flow `[re, im]` (scalar chain).
    pub mouth_over_input_flow: [f64; 2],
}

/// A scan of the passive bore with the quantities needed to read it.
#[derive(Debug, Clone)]
pub struct BoreScan {
    /// Plane-wave characteristic impedance of the input plane `rho c / (pi r^2)` [Pa s/m^3].
    pub input_characteristic_impedance: f64,
    /// Gas density [kg/m^3].
    pub density_kg_m3: f64,
    /// Distance from the outlet plane to the microphone [m].
    pub microphone_distance_m: f64,
    /// One row per requested frequency.
    pub rows: Vec<BoreResponse>,
}

/// Lip and bore statistics over one tenth of a second of a render: whether the lip closes and how large
/// the pressure at the mouthpiece is. A diagnostic of the loop, not a graded quantity.
#[derive(Debug, Clone, Copy)]
pub struct LipWindow {
    /// Start of the window [s].
    pub start_s: f64,
    /// Mean pressure at the mouthpiece (the mouth's DC operating point) [Pa].
    pub bore_mean_pa: f64,
    /// RMS of the alternating mouthpiece pressure [Pa].
    pub bore_ac_rms_pa: f64,
    /// Half the peak-to-peak mouthpiece pressure [Pa].
    pub bore_ac_peak_pa: f64,
    /// Smallest lip gap in the window [m] (0 = the lip closed).
    pub gap_min_m: f64,
    /// Mean lip gap [m] (mean of the per-chunk means).
    pub gap_mean_m: f64,
    /// Largest lip gap [m].
    pub gap_max_m: f64,
    /// Fraction of the control-tick chunks in which the lip closed at least once.
    pub closed_chunk_fraction: f64,
}

#[derive(Debug, Clone, Copy)]
struct WindowAcc {
    n: u64,
    sum: f64,
    sum_sq: f64,
    p_min: f64,
    p_max: f64,
    chunks: u64,
    closed_chunks: u64,
    gap_min: f64,
    gap_max: f64,
    gap_mean_sum: f64,
}

impl WindowAcc {
    const fn new() -> Self {
        Self {
            n: 0,
            sum: 0.0,
            sum_sq: 0.0,
            p_min: f64::INFINITY,
            p_max: f64::NEG_INFINITY,
            chunks: 0,
            closed_chunks: 0,
            gap_min: f64::INFINITY,
            gap_max: 0.0,
            gap_mean_sum: 0.0,
        }
    }
}

enum Engine {
    Linear {
        voice: BrassVoice,
        mouth: DelayedFilter,
    },
    #[cfg(feature = "nonlinear-brass")]
    Nonlinear(Box<super::nonlinear::NonlinearLip>),
}

impl Engine {
    fn apply(&mut self, control: BrassControl) -> Result<(), String> {
        match self {
            Self::Linear { voice, .. } => voice.apply(control).map_err(|e| e.to_string()),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.apply(control),
        }
    }

    fn diagnostics(&self) -> &[BrassBlockDiag] {
        match self {
            Self::Linear { voice, .. } => voice.diagnostics(),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.diagnostics(),
        }
    }
}

/// A realized performance: the lip, the bore and the microphone, advanced
/// together on one physical clock.
pub struct LipPerformance {
    info: LipInfo,
    engine: Engine,
    receiver: BaffledPressure,
    schedule: GestureSchedule,
    blow_track: String,
    tension_track: Option<String>,
    completed: u64,
    applied_pressure: f64,
    applied_tension: f64,
    poisoned: bool,
    mouthpiece: Vec<f64>,
    outgoing: Vec<f64>,
    trace: Vec<f64>,
    ratio: usize,
    window_len: u64,
    acc: Vec<WindowAcc>,
    tap: Option<Vec<f64>>,
}

fn refusal(what: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(what))
}

impl LipPerformance {
    /// The realized description.
    #[must_use]
    pub const fn info(&self) -> &LipInfo {
        &self.info
    }

    /// Existing microphone filter latency, in output samples. Physical
    /// source-to-listener travel remains in the pressure and is not included.
    #[must_use]
    pub fn filter_delay_samples(&self) -> f64 {
        self.receiver.filter_delay_samples()
    }

    /// The loop's per-chunk diagnostics (lock estimate, gap and power ledger).
    #[must_use]
    pub fn diagnostics(&self) -> &[BrassBlockDiag] {
        self.engine.diagnostics()
    }

    /// Start recording the pressure at the mouthpiece (the internal pressure the
    /// lip sees, DC included) for every sample rendered from now on.
    pub fn record_mouthpiece(&mut self) {
        self.tap.get_or_insert_with(Vec::new);
    }

    /// The recorded mouthpiece pressure [Pa] (empty unless [`Self::record_mouthpiece`] ran).
    #[must_use]
    pub fn mouthpiece_pressure(&self) -> &[f64] {
        self.tap.as_deref().unwrap_or(&[])
    }

    /// Lip and bore statistics per tenth of a second rendered so far.
    #[must_use]
    pub fn windows(&self) -> Vec<LipWindow> {
        let rate = f64::from(self.info.sample_rate_hz);
        self.acc
            .iter()
            .enumerate()
            .filter(|(_, a)| a.n > 0 && a.chunks > 0)
            .map(|(k, a)| {
                let n = a.n as f64;
                let mean = a.sum / n;
                LipWindow {
                    start_s: (k as u64 * self.window_len) as f64 / rate,
                    bore_mean_pa: mean,
                    bore_ac_rms_pa: (a.sum_sq / n - mean * mean).max(0.0).sqrt(),
                    bore_ac_peak_pa: 0.5 * (a.p_max - a.p_min),
                    gap_min_m: a.gap_min,
                    gap_mean_m: a.gap_mean_sum / a.chunks as f64,
                    gap_max_m: a.gap_max,
                    closed_chunk_fraction: a.closed_chunks as f64 / a.chunks as f64,
                }
            })
            .collect()
    }

    fn apply_controls(&mut self, tick: u64) -> Result<(), String> {
        let pressure = self
            .schedule
            .sample(&self.blow_track, tick)
            .map_err(|e| format!("pressure schedule refused: {e:?}"))?;
        if pressure != self.applied_pressure {
            self.engine
                .apply(BrassControl::SetBlowingPressure(pressure))
                .map_err(|e| format!("blowing pressure refused: {e}"))?;
            self.applied_pressure = pressure;
        }
        if let Some(track) = &self.tension_track {
            let tension = self
                .schedule
                .sample(track, tick)
                .map_err(|e| format!("lip tension schedule refused: {e:?}"))?;
            if tension != self.applied_tension {
                self.engine
                    .apply(BrassControl::SetLipTension(tension))
                    .map_err(|e| format!("lip tension refused: {e}"))?;
                self.applied_tension = tension;
            }
        }
        Ok(())
    }

    fn advance(&mut self, output: &mut [f64]) -> Result<(), String> {
        let rate = u64::from(self.info.sample_rate_hz);
        let control = u64::from(self.info.control_rate_hz);
        let mut done = 0usize;
        while done < output.len() {
            let n = self.completed + done as u64;
            let tick = tick_at(n, control, rate);
            self.apply_controls(tick)?;
            let until_next = tick_start(tick + 1, control, rate) - n;
            let chunk = usize::try_from(until_next)
                .unwrap_or(usize::MAX)
                .min(output.len() - done)
                .min(MAX_BLOCK)
                .max(1);
            match &mut self.engine {
                Engine::Linear { voice, .. } => voice
                    .step_block_outgoing(&mut self.mouthpiece[..chunk], &mut self.outgoing[..chunk])
                    .map_err(|e| format!("brass loop refused: {e}"))?,
                #[cfg(feature = "nonlinear-brass")]
                Engine::Nonlinear(lip) => lip.step_block(
                    &mut self.mouthpiece[..chunk],
                    &mut self.trace[..chunk * self.ratio],
                )?,
            }
            let diag = self.engine.diagnostics().last().cloned();
            let window = usize::try_from(n / self.window_len).unwrap_or(usize::MAX);
            while self.acc.len() <= window {
                self.acc.push(WindowAcc::new());
            }
            if let Some(tap) = self.tap.as_mut() {
                tap.extend_from_slice(&self.mouthpiece[..chunk]);
            }
            let acc = &mut self.acc[window];
            for &p in &self.mouthpiece[..chunk] {
                acc.n += 1;
                acc.sum += p;
                acc.sum_sq += p * p;
                acc.p_min = acc.p_min.min(p);
                acc.p_max = acc.p_max.max(p);
            }
            if let Some(d) = diag {
                acc.chunks += 1;
                acc.closed_chunks += u64::from(d.gap_min_m <= 0.0);
                acc.gap_min = acc.gap_min.min(d.gap_min_m);
                acc.gap_max = acc.gap_max.max(d.gap_max_m);
                acc.gap_mean_sum += d.gap_mean_m;
            }
            for i in 0..chunk {
                let frame = match &mut self.engine {
                    Engine::Linear { mouth, .. } => {
                        self.trace[0] = mouth
                            .push(self.outgoing[i])
                            .map_err(|e| format!("mouth-flow line refused: {e:?}"))?;
                        &self.trace[..self.ratio]
                    }
                    #[cfg(feature = "nonlinear-brass")]
                    Engine::Nonlinear(_) => &self.trace[i * self.ratio..(i + 1) * self.ratio],
                };
                output[done + i] = self.receiver.step_trace(frame)?;
            }
            done += chunk;
        }
        self.completed += output.len() as u64;
        Ok(())
    }
}

impl PressureRenderer for LipPerformance {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }

    fn max_block_len(&self) -> usize {
        MAX_BLOCK
    }

    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if sample_rate_hz != self.info.sample_rate_hz {
            return Err(RenderError::Sizing {
                what: "lip performance output clock differs from the requested clock",
            });
        }
        Ok(())
    }

    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if samples > self.info.samples - self.completed {
            return Err(RenderError::Sizing {
                what: "request runs past the samples the lip performance declares",
            });
        }
        Ok(())
    }

    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > MAX_BLOCK {
            return Err(RenderError::Sizing {
                what: "lip callback exceeds its admitted block length",
            });
        }
        self.validate_sample_count(output.len() as u64)?;
        match self.advance(output) {
            Ok(()) => Ok(()),
            Err(what) => {
                // A failed physical callback leaves a prefix advanced.
                self.poisoned = true;
                Err(refusal(what))
            }
        }
    }
}

#[cfg(test)]
#[path = "../render/schedule/reed/field_input/tests/support.rs"]
mod gas_test_support;

#[cfg(test)]
mod tests {
    use super::*;

    // The canonical pressure phrase of examples/reed-duct.performance.
    const SCHEDULE: &str = "frankensim-gesture-schedule-v1\n\
        control_rate_hz\t700\n\
        tracks\t1\n\
        track\tblow\tblowing-pressure\n\
        initial\tpressure-pa\t0e0\n\
        events\t5\n\
        event\t0e0\t1e-2\tpressure-pa\t2.8e3\n\
        event\t5e-3\t1e-2\tpressure-pa\t3.5e3\n\
        event\t2e-2\t0e0\tpressure-pa\t0e0\n\
        event\t4e-2\t5e-3\tpressure-pa\t2.8e3\n\
        event\t6e-2\t0e0\tpressure-pa\t0e0\n";

    // Test fixture, not data: a short brass-like bore and a lip-shaped
    // response so every owner in the chain is exercised.
    fn fixture(extra: &str) -> String {
        format!(
            "frankensim-lip-performance-v1\n\
             # test fixture, not data\n\
             audio 24000 2400\n\
             full_scale_pa 30000\n\
             ambient moist_air 293.15 101325 0\n\
             lip_response 200 7 0.1 0.012 0.0005 0.0001\n\
             lip_source test fixture, not a measurement\n\
             bore_source test fixture, not a measurement\n\
             cup_volume_m3 1.5e-6\n\
             termination flanged\n\
             modes 3 1\n\
             cylinder 0.006 0.30\n\
             cone 0.006 0.012 0.25\n\
             cone 0.012 0.035 0.12\n\
             microphone 0 0 1 32 128 3000\n\
             {extra}schedule\n{SCHEDULE}"
        )
    }

    fn with_gas(record: &str) -> String {
        fixture("").replace("moist_air 293.15 101325 0", record)
    }

    fn pressure(description: &LipDescription) -> Vec<f64> {
        let mut performance = description.build().expect("gas reaches all physics owners");
        assert_eq!(
            performance.info().sound_speed_m_s,
            description.ambient_gas().sound_speed
        );
        let mut result = vec![0.0; 1200];
        for block in result.chunks_mut(37) {
            performance.block(block).expect("physical pressure");
        }
        assert!(result.iter().all(|value| value.is_finite()));
        assert!(result.iter().any(|value| value.abs() > 1e-12));
        result
    }

    #[test]
    fn g3_shared_mixture_composition_and_pressure_change_the_radiated_lip_sound() {
        let legacy = LipDescription::from_bytes(fixture("").as_bytes()).unwrap();
        assert_eq!(
            *legacy.ambient_gas(),
            GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap()
        );
        let records = [
            "mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
            "mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2",
            "mixture 293.15 80000 mole 2 N2 0.8 O2 0.2",
        ];
        let mut rendered = Vec::new();
        for record in records {
            let description = LipDescription::from_bytes(with_gas(record).as_bytes()).unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, None).unwrap();
            assert_eq!(description.gas_record(), &expected);
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "composition must reach the audible loop"
        );
        assert_ne!(
            rendered[0], rendered[2],
            "pressure must reach flow and radiation"
        );
    }

    #[test]
    fn g3_nasa9_source_swap_changes_lip_pressure_and_retains_the_original_packs() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let record = "thermochemical 293.15 101325 mole 1 specimen-gas 1";
        let text = with_gas(record);
        assert!(
            LipDescription::from_bytes(text.as_bytes())
                .unwrap_err()
                .contains("explicit source registry")
        );
        let mut rendered = Vec::new();
        for cp_r in [3.5, 4.5] {
            let packs = gas_test_support::packs("specimen-gas", 0.028, cp_r, 200.0);
            let bindings = [packs.binding()];
            let sources = ThermochemicalSources {
                species: &bindings,
                water_species: None,
                phase: ThermochemicalPhase::FrozenGas,
            };
            let description =
                LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, Some(&sources))
                    .unwrap();
            assert_eq!(description.gas_record(), &expected);
            let GasRecord::Thermochemical(gas) = description.gas_record() else {
                panic!("retained source");
            };
            assert_eq!(gas.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "data-only heat-capacity swap must reach sound"
        );
    }

    #[test]
    fn gas_refusals_keep_source_admission_strict() {
        for record in [
            "mixture 293.15 101325 mole 2 N2 0.7 O2 0.2",
            "mixture 293.15 101325 mole 1 unavailable 1",
            "mixture 293.15 101325 mole 1 N2 1 extra",
            "moist_air 293.15 101325 0 extra",
            "thermochemical 293.15 101325 mole 1 specimen-gas 1",
        ] {
            assert!(
                LipDescription::from_bytes(with_gas(record).as_bytes()).is_err(),
                "{record}"
            );
        }
    }

    #[test]
    fn bessel_pieces_lie_on_the_horn_and_reach_both_ends() {
        let (r_in, r_out, length, flare) = (0.01105, 0.01507, 0.15042, 1.26);
        let mut segments = Vec::new();
        bessel_horn(r_in, r_out, length, flare, 16, &mut segments).expect("horn");
        assert_eq!(segments.len(), 16);
        let mut radii = vec![r_in];
        let mut total = 0.0;
        for segment in &segments {
            let Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } = *segment
            else {
                panic!("a Bessel piece is a cone");
            };
            assert_eq!(inlet_radius, *radii.last().expect("previous"));
            assert!(outlet_radius > inlet_radius);
            radii.push(outlet_radius);
            total += length;
        }
        assert!((total - length).abs() < 1e-15);
        assert_eq!(*radii.last().expect("end"), r_out);
        // r = b (u - s)^-flare, so r^(-1/flare) is affine in s: equal steps
        // in s give equal steps in r^(-1/flare).
        let affine: Vec<f64> = radii.iter().map(|r| det::pow(*r, -1.0 / flare)).collect();
        let first = affine[1] - affine[0];
        for pair in affine.windows(2) {
            assert!(((pair[1] - pair[0]) - first).abs() < 1e-9 * first.abs());
        }
        assert!(bessel_horn(r_out, r_in, length, flare, 4, &mut segments).is_err());
        assert!(bessel_horn(r_in, r_out, length, 0.0, 4, &mut segments).is_err());
    }

    #[test]
    fn lip_response_reproduces_the_published_quantities() {
        let (f, q, inverse_mu, area) = (204.5, 7.0, 0.063, 1.0e-4);
        let text = fixture("").replace(
            "lip_response 200 7 0.1 0.012 0.0005 0.0001",
            &format!("lip_response {f} {q} {inverse_mu} 0.012 0.00053 {area}"),
        );
        let description = LipDescription::from_bytes(text.as_bytes()).expect("parses");
        let lip = &description.lip;
        let resonance = (lip.stiffness_n_m / lip.mass_kg).sqrt() / core::f64::consts::TAU;
        let quality = lip.mass_kg * core::f64::consts::TAU * f / lip.damping_n_s_m;
        assert!((resonance - f).abs() < 1e-9 * f);
        assert!((quality - q).abs() < 1e-9 * q);
        assert!((lip.face_area_m2 / lip.mass_kg - inverse_mu).abs() < 1e-12);
        assert_eq!(lip.provenance, "test fixture, not a measurement");
    }

    #[test]
    fn the_mouth_flow_tail_tolerance_defaults_to_the_ladders_own_and_is_declarable() {
        let default = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        assert_eq!(default.mouth_flow_tail_tol, LineLadder::default().tail_tol);
        let declared = LipDescription::from_bytes(fixture("mouth_flow_tail_tol 1e-6\n").as_bytes())
            .expect("parses");
        assert_eq!(declared.mouth_flow_tail_tol, 1.0e-6);
    }

    #[test]
    fn control_ticks_apply_at_the_ceiling_of_their_start() {
        let (control, rate) = (700, 19_200);
        assert_eq!(tick_at(0, control, rate), 0);
        assert_eq!(tick_start(1, control, rate), 28);
        assert_eq!(tick_at(27, control, rate), 0);
        assert_eq!(tick_at(28, control, rate), 1);
        for tick in 0..2_000u64 {
            let start = tick_start(tick, control, rate);
            assert_eq!(tick_at(start, control, rate), tick);
            if start > 0 {
                assert_eq!(tick_at(start - 1, control, rate), tick - 1);
            }
        }
    }

    #[test]
    fn malformed_and_non_physical_files_are_refused_by_name() {
        for (extra, needle) in [
            ("bogus 1\n", "unknown key"),
            ("audio 24000 10\n", "only appear once"),
            ("cylinder 0.006\n", "missing field"),
            ("cone 0.006 0.012 0.1 0.1\n", "extra field"),
            ("microphone 0 0 -1 4 16 1000\n", "in front of the baffle"),
            ("lip 1e-4 100 0.1 0.01 0.0005 0.0001\n", "only appear once"),
            (
                "lip 1e-4 100 -1 0.01 0.0005 0.0001\n",
                "damping must be nonnegative",
            ),
            ("modes 0 1\n", "outside its admitted range"),
            ("mouth_flow_tail_tol 1\n", "below 1"),
            ("mouth_flow_tail_tol 0\n", "positive"),
            (
                "mouth_flow_tail_tol 1e-6\nmouth_flow_tail_tol 1e-6\n",
                "only appear once",
            ),
            ("bessel 0.01 0.02 0.1 1 0\n", "outside its admitted range"),
            ("bore linear 0.0024\n", "expected nonlinear"),
            (
                "bore nonlinear 0.0024 17 1 24 20 8000 1e-2\n",
                "outside its admitted range",
            ),
            (
                "bore nonlinear 0.0024 16 1 24 8000 20 1e-2\n",
                "rising wall band",
            ),
            ("bore nonlinear 0.0024 16 1 24 20 8000 1\n", "below 1"),
        ] {
            let error = LipDescription::from_bytes(fixture(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let no_bore_source =
            fixture("").replace("bore_source test fixture, not a measurement\n", "");
        let error = LipDescription::from_bytes(no_bore_source.as_bytes()).expect_err("source");
        assert!(error.contains("bore_source"), "{error}");
        let no_source_text =
            fixture("").replace("lip_source test fixture, not a measurement", "lip_source");
        let error = LipDescription::from_bytes(no_source_text.as_bytes()).expect_err("text");
        assert!(error.contains("provenance"), "{error}");
        let noncanonical = fixture("").replace("1e-2", "0.01");
        assert!(LipDescription::from_bytes(noncanonical.as_bytes()).is_err());
        assert!(LipDescription::from_bytes(b"not a lip file\n").is_err());
    }

    #[test]
    fn window_diagnostics_cover_the_render_and_describe_a_possible_lip() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert!(performance.windows().is_empty());
        let mut block = vec![0.0; 600];
        for _ in 0..4 {
            performance.block(&mut block).expect("block");
        }
        let windows = performance.windows();
        // 2400 samples at 24 kHz are exactly one tenth of a second.
        assert_eq!(windows.len(), 1);
        let w = windows[0];
        assert_eq!(w.start_s, 0.0);
        assert!(w.gap_min_m >= 0.0 && w.gap_min_m <= w.gap_mean_m && w.gap_mean_m <= w.gap_max_m);
        assert!((0.0..=1.0).contains(&w.closed_chunk_fraction));
        assert!(w.bore_ac_rms_pa.is_finite() && w.bore_ac_rms_pa <= w.bore_ac_peak_pa + 1e-9);
        assert!(w.bore_mean_pa.is_finite());
    }

    #[test]
    fn the_mouthpiece_tap_records_every_rendered_sample_only_when_asked() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert!(performance.mouthpiece_pressure().is_empty());
        let mut block = vec![0.0; 600];
        performance.block(&mut block).expect("block");
        assert!(performance.mouthpiece_pressure().is_empty());
        performance.record_mouthpiece();
        for _ in 0..3 {
            performance.block(&mut block).expect("block");
        }
        let tap = performance.mouthpiece_pressure();
        assert_eq!(tap.len(), 1800);
        assert!(tap.iter().all(|p| p.is_finite()));
    }

    #[test]
    fn the_bore_scan_reads_the_passive_bore_without_building_the_loop() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let scan = description
            .bore_scan(&[150.0, 300.0, 450.0])
            .expect("scans");
        assert_eq!(scan.rows.len(), 3);
        assert!(scan.input_characteristic_impedance > 0.0 && scan.microphone_distance_m > 0.0);
        for row in &scan.rows {
            let z = row.impedance_multimodal[0].hypot(row.impedance_multimodal[1]);
            assert!(z.is_finite() && z > 0.0);
            // a passive load: the real part of both input impedances is not negative
            assert!(row.impedance_multimodal[0] >= -1e-6 * z);
            assert!(row.impedance_scalar[0] >= -1e-6 * z);
            assert!(
                row.mouth_over_input_flow[0]
                    .hypot(row.mouth_over_input_flow[1])
                    .is_finite()
            );
        }
    }

    #[test]
    fn the_chain_builds_and_renders_finite_pressure_on_its_own_clock() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert_eq!(performance.info().samples, 2400);
        let LipBoreInfo::Linear(lines) = &performance.info().bore else {
            panic!("a file without a bore record realizes the lines");
        };
        assert!(lines.mouth_flow_taps >= 256);
        assert!(performance.validate_sample_rate(48_000).is_err());
        assert!(performance.validate_sample_rate(24_000).is_ok());
        assert!(performance.validate_sample_count(2401).is_err());
        let mut block = vec![0.0; 600];
        for _ in 0..4 {
            performance.block(&mut block).expect("block");
            assert!(block.iter().all(|p| p.is_finite()));
        }
        assert_eq!(performance.samples_rendered(), 2400);
        assert!(performance.block(&mut block).is_err());
    }

    const NONLINEAR_RECORD: &str = "bore nonlinear 0.0024 16 1 64 20 8000 1e-2\n";
    #[cfg(feature = "nonlinear-brass")]
    const SPECIMEN: &str = "thermochemical 293.15 101325 mole 1 specimen-gas 1";

    #[cfg(not(feature = "nonlinear-brass"))]
    #[test]
    fn a_nonlinear_bore_is_refused_by_name_without_its_feature() {
        let error = LipDescription::from_bytes(fixture(NONLINEAR_RECORD).as_bytes())
            .expect_err("no nonlinear bore in this build");
        assert!(error.contains("nonlinear-brass"), "{error}");
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn a_nonlinear_bore_evolves_only_a_source_backed_ambient_gas_through_a_flanged_bell() {
        let error = LipDescription::from_bytes(fixture(NONLINEAR_RECORD).as_bytes())
            .expect_err("a snapshot gas cannot evolve");
        assert!(error.contains("ambient thermochemical"), "{error}");
        let unflanged =
            fixture(NONLINEAR_RECORD).replace("termination flanged", "termination unflanged");
        let error = LipDescription::from_bytes(unflanged.as_bytes())
            .expect_err("the nonlinear bell is a baffled piston");
        assert!(error.contains("termination flanged"), "{error}");
    }

    /// One NASA-9 species: each region is (lowest K, highest K, a0..=a8).
    #[cfg(feature = "nonlinear-brass")]
    fn nasa9_species(
        name: &str,
        molar_mass: f64,
        reference_pressure: f64,
        elemental_reference: &str,
        regions: &[(f64, f64, [f64; 9])],
        source: &str,
    ) -> fs_thermochem::Nasa9StandardStateModelV1 {
        use fs_evidence::ValidityDomain;
        use fs_matdb::{
            ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance,
        };
        use fs_qty::{MolarMass, Pressure};
        use fs_thermochem::{
            ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
            NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
            ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
            required_nasa9_parameter_dims,
        };
        let regions = regions
            .iter()
            .map(|&(low, high, a)| {
                let mut parameters = std::collections::BTreeMap::new();
                for ((name, dims), value) in
                    required_nasa9_parameter_dims().into_iter().take(9).zip(a)
                {
                    parameters.insert(name.to_owned(), LawParameter { value, dims });
                }
                parameters.insert(
                    "reference_pressure".into(),
                    LawParameter {
                        value: reference_pressure,
                        dims: Pressure::DIMS,
                    },
                );
                Nasa9RegionV1::from_card(ConstitutiveModelCard {
                    law: LawId(NASA9_LAW_ID_V1.into()),
                    law_version: NASA9_LAW_VERSION_V1,
                    parameters,
                    state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                    initial_state: InitialStatePolicy::ZeroInternalState,
                    validity: ValidityDomain::unconstrained().with("T", low, high),
                    sources: Vec::new(),
                    provenance: Provenance {
                        source: source.into(),
                        license: "test fixture".into(),
                        artifact: None,
                    },
                })
                .unwrap()
            })
            .collect();
        Nasa9StandardStateModelV1::new(
            SpeciesId::new(name).unwrap(),
            MolarMass::new(molar_mass),
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(reference_pressure),
                ElementalReferenceIdV1::new(elemental_reference).unwrap(),
            )
            .unwrap(),
            regions,
        )
        .unwrap()
    }

    /// A data/matdb/seed-v1 species: molar mass, reference pressure and
    /// elemental reference from species.tsv; every region's bounds and
    /// coefficients from the McBride nasa9.tsv.
    #[cfg(feature = "nonlinear-brass")]
    fn seed_species(species: &str, nasa9: &str) -> fs_thermochem::Nasa9StandardStateModelV1 {
        let row: Vec<&str> = species
            .lines()
            .find(|line| line.starts_with("species\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(
            (row[3], row[4], row[5], row[7]),
            ("g/mol", "gas", "ideal-gas", "kPa"),
            "{row:?}"
        );
        let reference_pressure = row[6].parse::<f64>().unwrap() * 1e3;
        let mut regions: Vec<(String, f64, f64, [f64; 9])> = Vec::new();
        for line in nasa9.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.first() {
                Some(&"region") => {
                    assert_eq!((fields[5], fields[7]), ("K", "Pa"), "{line}");
                    assert_eq!(fields[6].parse::<f64>().unwrap(), reference_pressure);
                    regions.push((
                        fields[2].to_owned(),
                        fields[3].parse().unwrap(),
                        fields[4].parse().unwrap(),
                        [f64::NAN; 9],
                    ));
                }
                Some(&"coefficient") => {
                    let region = regions.iter_mut().find(|r| r.0 == fields[2]).unwrap();
                    let index: usize = fields[3].strip_prefix('a').unwrap().parse().unwrap();
                    region.3[index] = fields[4].parse().unwrap();
                }
                _ => {}
            }
        }
        let regions: Vec<_> = regions
            .into_iter()
            .map(|(_, low, high, a)| {
                assert!(a.iter().all(|c| c.is_finite()), "{a:?}");
                (low, high, a)
            })
            .collect();
        nasa9_species(
            row[1],
            row[2].parse::<f64>().unwrap() * 1e-3,
            reference_pressure,
            row[8],
            &regions,
            "McBride, Zehe and Gordon (2002), NASA/TP-2002-211556 (data/matdb/seed-v1)",
        )
    }

    // A species split at 1000 K, like every McBride (2002) air pack: the field
    // keeps to the room's region, which the source admits; the old fixed
    // 200-1000 K domain crossed the shared bound and refused.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_nonlinear_field_keeps_to_the_nasa9_region_that_holds_the_room() {
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
        use fs_thermochem::{Composition, CompositionBasis};
        let mut a = [0.0; 9];
        a[2] = 3.5;
        let source = FrozenIdealGasMixtureModelV1::new(
            vec![nasa9_species(
                "split-gas",
                0.029,
                1.0e5,
                "split-fixture",
                &[(200.0, 1000.0, a), (1000.0, 6000.0, a)],
                "two-region test fixture, not data",
            )],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap();
        let crossing = FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
            || false,
        )
        .expect_err("1000 K belongs to the upper region");
        assert!(
            format!("{crossing}").contains("caloric continuity"),
            "{crossing}"
        );
        let [low, high] = room_region(&source, 293.15).unwrap();
        assert_eq!(low, 200.0);
        assert!(high < 1000.0 && high == 1000.0_f64.next_down(), "{high}");
        FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(low), Temperature::new(high)),
            || false,
        )
        .unwrap();
        // A hot room is held in the upper region; a room off the source refuses.
        assert_eq!(room_region(&source, 1200.0).unwrap(), [1000.0, 6000.0]);
        let error = room_region(&source, 150.0).unwrap_err();
        assert!(error.contains("outside every NASA-9 region"), "{error}");
    }

    // The air the trumpet A/B renders breathe: the McBride (2002) N2, O2, Ar
    // and CO2 tables of data/matdb/seed-v1, each split at 1000 K with no
    // stated caloric continuity. The old fixed 200-1000 K domain refused this
    // air by name; the room's region admits it, a room state round-trips at
    // the ideal-gas sound speed of dry air, and a state above the region
    // refuses rather than extrapolating. The sound-speed bound (0.5 m/s about
    // Cramer's 343.2 m/s for dry air at 20 C, JASA 93 (1993) 2510) was fixed
    // before the first run.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_nonlinear_field_admits_the_mcbride_air_in_the_rooms_region() {
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
        use fs_thermochem::{Composition, CompositionBasis};
        macro_rules! seed {
            ($name:literal) => {
                seed_species(
                    include_str!(concat!(
                        env!("CARGO_MANIFEST_DIR"),
                        "/../../data/matdb/seed-v1/",
                        $name,
                        "/species.tsv"
                    )),
                    include_str!(concat!(
                        env!("CARGO_MANIFEST_DIR"),
                        "/../../data/matdb/seed-v1/",
                        $name,
                        "-nasa9-mcbride2002/nasa9.tsv"
                    )),
                )
            };
        }
        let air = vec![
            seed!("nitrogen"),
            seed!("oxygen"),
            seed!("argon"),
            seed!("carbon-dioxide"),
        ];
        for species in &air {
            let regions = species.regions();
            assert_eq!(regions.len(), 2, "{:?}", species.species());
            assert_eq!(regions[0].temperature_max().value(), 1000.0);
            assert_eq!(regions[1].temperature_min().value(), 1000.0);
            assert!(species.caloric_reference().is_none());
        }
        // Dry air by mole: N2 78.084, O2 20.946, Ar 0.934, CO2 0.036 %.
        let source = FrozenIdealGasMixtureModelV1::new(
            air,
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![0.78084, 0.20946, 0.00934, 0.00036],
            )
            .unwrap(),
        )
        .unwrap();
        let crossing = FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
            || false,
        )
        .expect_err("the McBride air states no continuity at 1000 K");
        assert!(
            format!("{crossing}").contains("caloric continuity"),
            "{crossing}"
        );
        let [low, high] = room_region(&source, 293.15).unwrap();
        assert_eq!([low, high], [200.0, 1000.0_f64.next_down()]);
        let gas = FrozenGas::new(
            source,
            CaloricDomain::new(Temperature::new(low), Temperature::new(high)),
            || false,
        )
        .unwrap();
        let r = gas.gas_constant();
        let room = gas
            .conserved(101_325.0 / (r * 293.15), 0.0, 101_325.0)
            .unwrap();
        let state = gas.primitive(room, || false).unwrap();
        println!(
            "mcbride air: R {r:.4} J/(kg K), T {:.9} K, c {:.4} m/s",
            state.temperature, state.sound_speed
        );
        assert!((state.temperature - 293.15).abs() <= 1e-9, "{state:?}");
        assert!((state.sound_speed - 343.2).abs() <= 0.5, "{state:?}");
        let hot = gas
            .conserved(101_325.0 / (r * 1200.0), 0.0, 101_325.0)
            .expect_err("1200 K lies above the room's region");
        assert!(
            format!("{hot}").contains("outside declared interval"),
            "{hot}"
        );
    }

    // The same lip, bore, gas and steady blow through both bore engines: the
    // nonlinear bore must speak in the regime the linear loop finds. The
    // tolerance (3 % on the playing frequency) was fixed before the first run.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_lip_speaks_on_the_nonlinear_bore_in_the_linear_loops_regime() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        const STEADY: &str = "frankensim-gesture-schedule-v1\n\
            control_rate_hz\t700\n\
            tracks\t1\n\
            track\tblow\tblowing-pressure\n\
            initial\tpressure-pa\t0e0\n\
            events\t1\n\
            event\t0e0\t1e-2\tpressure-pa\t3e3\n";
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let text = |extra: &str| {
            fixture(extra)
                .replace("moist_air 293.15 101325 0", SPECIMEN)
                .replace("audio 24000 2400", "audio 24000 4800")
                .replace(SCHEDULE, STEADY)
        };
        let render = |text: String| {
            let description =
                LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let mut performance = description.build().unwrap();
            performance.record_mouthpiece();
            let started = std::time::Instant::now();
            let mut radiated = vec![0.0; 4800];
            for block in radiated.chunks_mut(600) {
                performance.block(block).unwrap();
            }
            let seconds = started.elapsed().as_secs_f64();
            // The last 0.1 s: AC level and the playing frequency from the
            // interpolated first and last rising mean crossings.
            let tail = &performance.mouthpiece_pressure()[2400..];
            let mean = tail.iter().sum::<f64>() / tail.len() as f64;
            let rms = (tail.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>()
                / tail.len() as f64)
                .sqrt();
            let rising: Vec<f64> = tail
                .windows(2)
                .enumerate()
                .filter(|(_, w)| w[0] < mean && w[1] >= mean)
                .map(|(i, w)| i as f64 + (mean - w[0]) / (w[1] - w[0]))
                .collect();
            let frequency = match rising.as_slice() {
                [first, .., last] => (rising.len() - 1) as f64 * 24_000.0 / (last - first),
                _ => 0.0,
            };
            let level = radiated.iter().map(|p| p * p).sum::<f64>().sqrt();
            assert!(radiated.iter().all(|p| p.is_finite()) && level > 0.0);
            (performance.info().bore.clone(), rms, frequency, seconds)
        };
        let (_, rms_linear, f_linear, s_linear) = render(text(""));
        let (bore, rms_nonlinear, f_nonlinear, s_nonlinear) = render(text(NONLINEAR_RECORD));
        let LipBoreInfo::Nonlinear(receipt) = bore else {
            panic!("the record realizes the nonlinear bore");
        };
        println!(
            "lip regime: linear {f_linear:.2} Hz rms {rms_linear:.1} Pa ({s_linear:.1} s) | nonlinear {f_nonlinear:.2} Hz rms {rms_nonlinear:.1} Pa ({s_nonlinear:.1} s) | {receipt:?}"
        );
        assert!(
            rms_linear > 150.0,
            "the linear loop must speak: {rms_linear} Pa"
        );
        assert!(
            rms_nonlinear > 150.0,
            "the nonlinear loop must speak: {rms_nonlinear} Pa"
        );
        assert!((f_nonlinear / f_linear - 1.0).abs() <= 0.03);
    }
}
